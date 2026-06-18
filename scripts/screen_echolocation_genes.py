#!/usr/bin/env python3
"""
screen_echolocation_genes.py
----------------------------
Runs a genome-wide screen across all 10,511 genes in the MEME database to identify
individual genes showing a statistically significant enrichment of positive selection
substitutions on canonical echolocating lineages compared to canonical controls.

Outputs:
  - CSV report: scratch/echolocation_screen_results.csv
  - Top 20 table in stdout
"""

import os
import sys
import sqlite3
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

DB_PATH = "meme_results.db"
OUTPUT_CSV = "scratch/echolocation_screen_results.csv"

# Taxonomic Genera sets
MEGABAT_GENERA = {"pteropus", "rousettus", "eidolon", "cynopterus", "macroglossus"}
ODONTOCETE_GENERA = {
    "mesoplodon", "ziphius", "hyperoodon", "stenella", "delphinus", "tursiops", "sousa",
    "pseudorca", "peponocephala", "globicephala", "grampus", "steno", "sagmatias",
    "cephalorhynchus", "lagenorhynchus", "leucopleurus", "delphinapterus", "phocoena",
    "neophocaena", "lipotes", "platanista", "kogia", "physeter", "monodon", "inia", "pontoporia"
}
MYSTICETE_GENERA = {"megaptera", "balaenoptera", "eubalaena", "balaena", "eschrichtius", "caperea"}

def classify_species(label, node_id):
    label_lower = label.lower()
    match = re.search(r'\(([^)]+)\)', label_lower)
    if not match:
        return 0
        
    sci_name = match.group(1).strip()
    parts = sci_name.split()
    if not parts:
        return 0
    genus = parts[0]
    
    # 1. Cetaceans
    if genus in ODONTOCETE_GENERA:
        return 1
    if genus in MYSTICETE_GENERA:
        return 0
        
    # 2. Bats
    is_bat = (re.search(r'\bbats?\b', label_lower) is not None or 
              "flying fox" in label_lower or 
              "pipistrelle" in label_lower or 
              "rousette" in label_lower)
    if "bat-eared fox" in label_lower:
        is_bat = False
        
    if is_bat:
        if genus in MEGABAT_GENERA:
            return 0
        return 1
        
    return 0

def manual_fdr(pvals, alpha=0.05):
    """Failsafe manual Benjamini-Hochberg FDR correction using only numpy."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    by_descend = pvals.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = np.arange(n, 0, -1)
    qvals = np.minimum.accumulate(pvals[by_descend] * n / steps)
    qvals = qvals[by_orig]
    qvals[qvals > 1.0] = 1.0
    return qvals < alpha, qvals

def main():
    print(f"[*] Connecting to database at '{DB_PATH}'...")
    if not os.path.exists(DB_PATH):
        print(f"[!] Error: Database file '{DB_PATH}' not found.")
        sys.exit(1)
        
    conn = sqlite3.connect(DB_PATH)
    
    # 1. Get Leaf Nodes and Classify
    print("[*] Retrieving and classifying species...")
    nodes_df = pd.read_sql_query("SELECT node_id, node_label FROM master_nodes WHERE is_leaf = 1", conn)
    nodes_df["is_echolocator"] = nodes_df.apply(lambda r: classify_species(r["node_label"], r["node_id"]), axis=1)
    nodes_df = nodes_df.dropna(subset=["is_echolocator"])
    nodes_df["is_echolocator"] = nodes_df["is_echolocator"].astype(int)
    
    echolocators = set(nodes_df[nodes_df["is_echolocator"] == 1]["node_id"])
    non_echolocators = set(nodes_df[nodes_df["is_echolocator"] == 0]["node_id"])
    target_species = echolocators | non_echolocators
    
    print(f"    - Echolocating species: {len(echolocators)}")
    print(f"    - Background / Control species: {len(non_echolocators)}")
    
    # 2. Query positive selection substitutions at leaf nodes using fast gene-by-gene lookup
    print("[*] Querying selection events from database using gene-by-gene lookup...")
    import time
    t_start_query = time.time()
    
    # Get all genes with selected sites
    selected_genes = [r[0] for r in conn.execute("SELECT DISTINCT gene_name FROM site_results WHERE p_value < 0.05").fetchall()]
    all_genes = pd.read_sql_query("SELECT DISTINCT gene_name FROM gene_results", conn)["gene_name"].tolist()
    print(f"    - Found {len(selected_genes)} genes with selected sites.")
    
    records = []
    for idx, gene in enumerate(selected_genes):
        if idx % 2000 == 0 and idx > 0:
            print(f"      - Querying gene {idx}/{len(selected_genes)}...")
            
        sites = [r[0] for r in conn.execute("SELECT site_index FROM site_results WHERE gene_name = ? AND p_value < 0.05", (gene,)).fetchall()]
        if not sites:
            continue
            
        placeholders = ','.join(['?'] * len(sites))
        sql = f"""
        SELECT bm.master_node_id, COUNT(*)
        FROM site_substitutions ss
        JOIN branch_mappings bm ON ss.gene_name = bm.gene_name AND ss.branch_name = bm.branch_name
        JOIN master_nodes mn ON bm.master_node_id = mn.node_id
        WHERE ss.gene_name = ?
          AND ss.site_index IN ({placeholders})
          AND ss.is_synonymous = 0
          AND mn.is_leaf = 1
        GROUP BY bm.master_node_id
        """
        res = conn.execute(sql, [gene] + sites).fetchall()
        for node_id, count in res:
            records.append({
                "gene_name": gene,
                "master_node_id": node_id,
                "count": count
            })
            
    df = pd.DataFrame(records)
    conn.close()
    
    print(f"    - Retrieved {len(df)} total non-synonymous substitutions at selected sites in {time.time()-t_start_query:.2f}s.")
    
    # Filter for target species (all classified species)
    df = df[df["master_node_id"].isin(target_species)]
    df["is_echolocator"] = df["master_node_id"].apply(lambda x: 1 if x in echolocators else 0)
    print(f"    - Filtered to {len(df)} substitutions on classified lineages.")
    
    # 3. Aggregate by gene
    print("[*] Aggregating substitutions by gene...")
    # Group by gene and calculate echo/control counts
    gene_counts = []
    for gene_name, group in df.groupby("gene_name"):
        echo_count = int(group[group["is_echolocator"] == 1]["count"].sum())
        control_count = int(group[group["is_echolocator"] == 0]["count"].sum())
        gene_counts.append({
            "gene_name": gene_name,
            "echo_subs": echo_count,
            "control_subs": control_count
        })
        
    counts_df = pd.DataFrame(gene_counts)
    if counts_df.empty:
        print("[!] No selection events found on leaf branches.")
        sys.exit(0)
        
    # Background totals across all selected site substitutions
    total_echo_subs = int(counts_df["echo_subs"].sum())
    total_control_subs = int(counts_df["control_subs"].sum())
    
    print(f"    - Total target positive-selection substitutions: Echolocating={total_echo_subs}, Background={total_control_subs}")
    
    # 4. Perform Fisher's Exact Test per gene
    print("[*] Running Fisher's exact tests for all genes...")
    results = []
    for _, row in counts_df.iterrows():
        a = row["echo_subs"]
        c = row["control_subs"]
        
        # Contingency table:
        #            | This Gene | Other Genes (Background)
        # Echo       | a         | total_echo_subs - a
        # Control    | c         | total_control_subs - c
        b = total_echo_subs - a
        d = total_control_subs - c
        
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append({
            "gene_name": row["gene_name"],
            "echo_subs": a,
            "control_subs": c,
            "odds_ratio": round(odds_ratio, 4),
            "p_value": p_value
        })
        
    res_df = pd.DataFrame(results)
    
    # 5. Multiple testing correction
    print("[*] Adjusting P-values for multiple testing (Benjamini-Hochberg)...")
    _, q_values = manual_fdr(res_df["p_value"].values, alpha=0.05)
    res_df["fdr_q_value"] = q_values
    
    # Sort and save
    res_df = res_df.sort_values(by="p_value")
    
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    res_df.to_csv(OUTPUT_CSV, index=False)
    print(f"🎉 Complete screen results saved to '{OUTPUT_CSV}'")
    
    # 6. Display top results
    print("\n" + "="*80)
    print("🏆 Top 25 Genes with Positive Selection Enriched in Echolocators")
    print("="*80)
    print(f"{'Rank':<5} {'Gene Name':<15} {'Echo Subs':<10} {'Control Subs':<13} {'Odds Ratio':<12} {'P-value':<10} {'FDR Q-value':<12}")
    print("-"*80)
    
    top_df = res_df.head(25)
    for idx, row in enumerate(top_df.itertuples()):
        print(f"{idx+1:<5} {row.gene_name:<15} {row.echo_subs:<10} {row.control_subs:<13} {row.odds_ratio:<12.3f} {row.p_value:<10.2e} {row.fdr_q_value:<12.4f}")
    print("="*80)

if __name__ == "__main__":
    main()
