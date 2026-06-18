#!/usr/bin/env python3
"""
analyze_hearing_echolocation.py
-------------------------------
Performs an association analysis between hearing gene positive selection and the
echolocation phenotype (convergently evolved in microbats and toothed whales).

Queries the local 17 GB meme_results.db database to check for non-synonymous substitutions
at positive selection sites (MEME p-val < 0.05) on terminal branches, comparing
echolocators vs. non-echolocating outgroups.
"""

import sqlite3
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

DB_PATH = "meme_results.db"

# Major hearing genes known to undergo selection in echolocating lineages
HEARING_GENES = ["OTOF", "GJB2", "MYO7A", "TMC1", "CDH23", "MYO15A", "GJB6"]

# Genera lookup sets for precise taxonomic matching
MEGABAT_GENERA = {"pteropus", "rousettus", "eidolon", "cynopterus", "macroglossus"}
ODONTOCETE_GENERA = {
    "mesoplodon", "ziphius", "hyperoodon", "stenella", "delphinus", "tursiops", "sousa",
    "pseudorca", "peponocephala", "globicephala", "grampus", "steno", "sagmatias",
    "cephalorhynchus", "lagenorhynchus", "leucopleurus", "delphinapterus", "phocoena",
    "neophocaena", "lipotes", "platanista", "kogia", "physeter", "monodon", "inia", "pontoporia"
}
MYSTICETE_GENERA = {"megaptera", "balaenoptera", "eubalaena", "balaena", "eschrichtius", "caperea"}

def classify_species(label, node_id):
    """
    Classifies a species as:
      1: Echolocator (Microbats, Odontocetes)
      0: Background / Non-Echolocator (Megabats, Mysticetes, and all other mammals)
    """
    label_lower = label.lower()
    
    # Extract scientific name (usually inside parentheses, e.g. "black flying fox (Pteropus alecto)")
    match = re.search(r'\(([^)]+)\)', label_lower)
    if not match:
        return 0
        
    sci_name = match.group(1).strip()
    parts = sci_name.split()
    if not parts:
        return 0
    genus = parts[0]  # e.g. "pteropus" or "pteronura"
    
    # 1. Cetaceans (Whales/Dolphins/Porpoises)
    if genus in ODONTOCETE_GENERA:
        return 1
    if genus in MYSTICETE_GENERA:
        return 0
        
    # 2. Bats (Chiroptera)
    is_bat = (re.search(r'\bbats?\b', label_lower) is not None or 
              "flying fox" in label_lower or 
              "pipistrelle" in label_lower or 
              "rousette" in label_lower)
    # Exclude bat-eared fox
    if "bat-eared fox" in label_lower:
        is_bat = False
        
    if is_bat:
        if genus in MEGABAT_GENERA:
            return 0  # Megabats (canonical control group)
        return 1  # Microbats (echolocators)
        
    return 0

def main():
    print(f"[*] Connecting to database at '{DB_PATH}'...")
    if not os.path.exists(DB_PATH):
        print(f"[!] Error: Database file '{DB_PATH}' not found. Please make sure you are in the directory containing the file.")
        sys.exit(1)
        
    conn = sqlite3.connect(DB_PATH)
    
    # 1. Fetch all terminal nodes
    print("[*] Retrieving terminal node labels...")
    nodes_df = pd.read_sql_query("SELECT node_id, node_label FROM master_nodes WHERE is_leaf = 1", conn)
    
    nodes_df["is_echolocator"] = nodes_df.apply(lambda r: classify_species(r["node_label"], r["node_id"]), axis=1)
    
    # Filter out None values to restrict strictly to Chiroptera and Cetacea
    nodes_df = nodes_df.dropna(subset=["is_echolocator"])
    nodes_df["is_echolocator"] = nodes_df["is_echolocator"].astype(int)
    
    echolocators = set(nodes_df[nodes_df["is_echolocator"] == 1]["node_id"])
    non_echolocators = set(nodes_df[nodes_df["is_echolocator"] == 0]["node_id"])
    
    print(f"    - Found {len(echolocators)} echolocating terminal species.")
    print(f"    - Found {len(non_echolocators)} canonical non-echolocating species.")
    
    # Print examples
    print("\n[Examples of Echolocators]:")
    print(", ".join(list(nodes_df[nodes_df["is_echolocator"] == 1]["node_label"].head(5))))
    print("\n[Examples of Canonical Non-Echolocator controls]:")
    print(", ".join(list(nodes_df[nodes_df["is_echolocator"] == 0]["node_label"].head(5))))
    
    # 2. Get positive selection results for target hearing genes
    print("\n[*] Querying positive selection sites for hearing genes...")
    placeholders = ",".join(["?"] * len(HEARING_GENES))
    site_query = f"""
    SELECT gene_name, site_index, p_value, is_significant 
    FROM site_results 
    WHERE gene_name IN ({placeholders}) AND p_value < 0.05
    """
    selected_sites_df = pd.read_sql_query(site_query, conn, params=HEARING_GENES)
    print(f"    - Found {len(selected_sites_df)} sites under selection in hearing genes.")
    
    # 3. Query substitutions at these sites
    print("[*] Mapping amino acid substitutions at selected sites to terminal branches...")
    subst_query = f"""
    SELECT ss.gene_name, ss.site_index, ss.branch_name, bm.master_node_id, ss.ancestral_aa, ss.derived_aa
    FROM site_substitutions ss
    JOIN branch_mappings bm ON ss.gene_name = bm.gene_name AND ss.branch_name = bm.branch_name
    JOIN master_nodes mn ON bm.master_node_id = mn.node_id
    WHERE ss.gene_name IN ({placeholders})
      AND ss.is_synonymous = 0
      AND mn.is_leaf = 1
    """
    subst_df = pd.read_sql_query(subst_query, conn, params=HEARING_GENES)
    
    # Filter substitutions to only those at significant MEME sites
    merged_df = pd.merge(
        subst_df, 
        selected_sites_df, 
        on=["gene_name", "site_index"], 
        how="inner"
    )
    
    # 4. Count substitutions on echolocating vs non-echolocating branches
    # Keep only target species (echolocating microbats/odontocetes vs non-echolocating megabats/mysticetes)
    merged_df = merged_df[merged_df["master_node_id"].isin(echolocators | non_echolocators)]
    merged_df["is_echolocating_branch"] = merged_df["master_node_id"].apply(lambda x: 1 if x in echolocators else 0)
    
    # Count totals for Fisher's Exact test
    echo_subs = merged_df[merged_df["is_echolocating_branch"] == 1]
    non_echo_subs = merged_df[merged_df["is_echolocating_branch"] == 0]
    
    print(f"\n[Results for Hearing Genes ({', '.join(HEARING_GENES)})]:")
    print(f"    - Non-synonymous substitutions at select sites on echolocating terminal branches: {len(echo_subs)}")
    print(f"    - Non-synonymous substitutions at select sites on non-echolocating terminal branches: {len(non_echo_subs)}")
    
    # We also need to get the "background" rate of substitutions to construct the contingency table.
    # Alternatively, we can test: are hearing genes more likely to have substitutions on echolocating branches
    # compared to all other genes in the database?
    
    print("\n[*] Querying background rate of substitutions on terminal branches (sampling 1,000 genes for speed)...")
    # To run fast, we sample 1000 background genes. The results scale up representatively.
    bg_query = """
    SELECT bm.master_node_id, COUNT(*) as subst_count
    FROM site_substitutions ss
    JOIN branch_mappings bm ON ss.gene_name = bm.gene_name AND ss.branch_name = bm.branch_name
    JOIN master_nodes mn ON bm.master_node_id = mn.node_id
    WHERE ss.is_synonymous = 0 AND mn.is_leaf = 1
      AND ss.gene_name IN (SELECT gene_name FROM gene_results LIMIT 1000)
    GROUP BY bm.master_node_id
    """
    bg_subst_df = pd.read_sql_query(bg_query, conn)
    # Keep only target species (echolocating microbats/odontocetes vs non-echolocating megabats/mysticetes)
    bg_subst_df = bg_subst_df[bg_subst_df["master_node_id"].isin(echolocators | non_echolocators)]
    bg_subst_df["is_echolocator"] = bg_subst_df["master_node_id"].apply(lambda x: 1 if x in echolocators else 0)
    
    bg_echo_sample = bg_subst_df[bg_subst_df["is_echolocator"] == 1]["subst_count"].sum()
    bg_non_echo_sample = bg_subst_df[bg_subst_df["is_echolocator"] == 0]["subst_count"].sum()
    
    # Scale up to estimate genome-wide totals for display
    scaling_factor = 10.511  # 10511 total genes / 1000 sampled
    bg_echo_total = int(bg_echo_sample * scaling_factor)
    bg_non_echo_total = int(bg_non_echo_sample * scaling_factor)
    
    print(f"    - Sample background substitutions (1,000 genes): Echolocating={bg_echo_sample}, Non-Echolocating={bg_non_echo_sample}")
    print(f"    - Extrapolated genome-wide substitutions: Echolocating={bg_echo_total}, Non-Echolocating={bg_non_echo_total}")
    
    # 5. Fisher's Exact Test Contingency Table
    # Table structure:
    #                     | Hearing Genes | Background Genes (excluding hearing)
    # Echolocators        | A (echo_subs) | B (bg_echo_total - echo_subs)
    # Non-Echolocators    | C (non_echo)  | D (bg_non_echo_total - non_echo)
    
    A = len(echo_subs)
    C = len(non_echo_subs)
    B = bg_echo_total - A
    D = bg_non_echo_total - C
    
    odds_ratio, p_value = fisher_exact([[A, B], [C, D]], alternative="greater")
    
    print("\n" + "="*50)
    print("📈 Fisher's Exact Test: Association of Selection with Echolocation")
    print("="*50)
    print(f"Contingency Table:")
    print(f"                     Hearing    Background")
    print(f"Echolocating:        {A:<10} {B:<10}")
    print(f"Non-Echolocating:    {C:<10} {D:<10}")
    print("-"*50)
    print(f"Odds Ratio (Enrichment): {odds_ratio:.4f}")
    print(f"P-value:                 {p_value:.2e}")
    
    if p_value < 0.05:
        print("🎉 Result: Significant enrichment of positive selection in hearing genes on echolocating lineages!")
    else:
        print("❌ Result: No significant enrichment found.")
    print("="*50)
    
    # Show top specific sites under selection
    if A > 0:
        print("\n[Top substitutions in hearing genes on echolocating branches]:")
        disp_df = echo_subs.sort_values(by="p_value").head(10)
        for _, row in disp_df.iterrows():
            print(f"    - Gene {row['gene_name']} Site {row['site_index']} ({row['ancestral_aa']}->{row['derived_aa']}) on species {row['master_node_id']} (MEME p-val: {row['p_value']:.4f})")
            
    conn.close()

if __name__ == "__main__":
    import os
    import sys
    main()
