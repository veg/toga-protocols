#!/usr/bin/env python3
"""
run_continuous_phenotype_association.py
---------------------------------------
Runs a phylogenetic regression (PGLS) genome-wide selection screen. It associates 
the occurrence of leaf-level selection events (from meme_results.db) with a continuous 
phenotype (e.g. body size) across species, controlling for phylogenetic shared ancestry.

Methodology:
1. Parse the species tree (Newick format) and construct the phylogenetic covariance matrix V.
2. Load species phenotypes. If no phenotype file is specified, simulates body mass (log g)
   along the tree under a Brownian Motion process (clade-stratified for realistic sizes).
3. Retrieve non-synonymous selection events on terminal branches for each gene in the database.
4. Perform PGLS regression for each gene using the Cholesky-transformed variables:
     Y* = L^-1 * Y
     X* = L^-1 * X
     W* = L^-1 * 1 (intercept)
   Then perform OLS: Y* ~ X* + W* - 1 (passing through the origin).
5. Output results ranked by significance, corrected for multiple testing using FDR.
"""

import os
import sys
import re
import sqlite3
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.stats as stats

DB_PATH = "meme_results.db"
TREE_PATH = "toga_protocols/data/pruned_species_tree.nhx"
OUTPUT_CSV = "scratch/phenotype_association_results.csv"

# Custom Newick Parser & Node structure to compute covariance without heavy bio packages
class Node:
    def __init__(self, name=None, length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None

def parse_newick(newick_str):
    # Strip NHX annotations and whitespace
    newick_str = re.sub(r'\[[^\]]*\]', '', newick_str)
    newick_str = newick_str.strip().rstrip(';')
    
    stack = [Node()]
    current = stack[0]
    i = 0
    while i < len(newick_str):
        char = newick_str[i]
        if char == '(':
            new_node = Node()
            current.children.append(new_node)
            new_node.parent = current
            stack.append(new_node)
            current = new_node
            i += 1
        elif char == ')':
            current = stack.pop()
            # Parse label/length
            i += 1
            label = ""
            while i < len(newick_str) and newick_str[i] not in (',', ')', '('):
                label += newick_str[i]
                i += 1
            if ':' in label:
                name, dist = label.split(':', 1)
                current.name = name.strip() if name.strip() else None
                try:
                    current.length = float(dist.strip())
                except:
                    current.length = 0.0
            else:
                current.name = label.strip() if label.strip() else None
        elif char == ',':
            stack.pop()
            parent = stack[-1]
            new_node = Node()
            parent.children.append(new_node)
            new_node.parent = parent
            stack.append(new_node)
            current = new_node
            i += 1
        else:
            label = ""
            while i < len(newick_str) and newick_str[i] not in (',', ')', '('):
                label += newick_str[i]
                i += 1
            if ':' in label:
                name, dist = label.split(':', 1)
                current.name = name.strip() if name.strip() else None
                try:
                    current.length = float(dist.strip())
                except:
                    current.length = 0.0
            else:
                current.name = label.strip() if label.strip() else None
    return stack[0]

def get_paths_from_root(root):
    paths = {}
    def dfs(node, current_path):
        path_len = current_path[-1][1] + node.length if current_path else node.length
        new_path = current_path + [(node, path_len)]
        if not node.children:  # Leaf node
            if node.name:
                paths[node.name] = new_path
        for child in node.children:
            dfs(child, new_path)
    dfs(root, [])
    return paths

def compute_covariance_matrix(paths, leaf_names):
    n = len(leaf_names)
    V = np.zeros((n, n))
    for i in range(n):
        path_i = paths[leaf_names[i]]
        for j in range(i, n):
            path_j = paths[leaf_names[j]]
            # Find the last common node in path_i and path_j
            common_len = 0.0
            min_len = min(len(path_i), len(path_j))
            for k in range(min_len):
                if path_i[k][0] == path_j[k][0]:
                    common_len = path_i[k][1]
                else:
                    break
            V[i, j] = common_len
            V[j, i] = common_len
    return V

def simulate_phenotypes(paths, leaf_names, seed=42):
    """Simulates realistic log body mass (g) for species under Brownian Motion."""
    np.random.seed(seed)
    n = len(leaf_names)
    V = compute_covariance_matrix(paths, leaf_names)
    
    # Regularize V to make sure it is strictly positive definite
    V += np.eye(n) * 1e-6
    
    # Generate traits under BM: Y ~ N(mean, sigma^2 * V)
    # Average log body mass for mammals is ~5.0 (approx 150g), std is ~2.0
    mean_val = 5.0
    sigma = 2.0
    
    # Cholesky decomposition of V
    L = la.cholesky(V, lower=True)
    z = np.random.normal(0, 1, n)
    traits = mean_val + sigma * (L @ z)
    
    return pd.DataFrame({
        "species_id": leaf_names,
        "log_body_mass": traits
    })

def perform_pgls(y, x, L):
    """
    Fits PGLS: y = b0 + b1*x + e, e ~ N(0, V)
    Uses the transformation L^-1 * y = b0*(L^-1 * 1) + b1*(L^-1 * x) + n
    """
    n = len(y)
    
    # Solve lower triangular systems L @ y_trans = y  => y_trans = L^-1 @ y
    y_trans = la.solve_triangular(L, y, lower=True)
    x_trans = la.solve_triangular(L, x, lower=True)
    w_trans = la.solve_triangular(L, np.ones(n), lower=True)
    
    # Build design matrix (X_mat)
    X_mat = np.column_stack((w_trans, x_trans))
    
    # Fit OLS: y_trans = X_mat @ beta + error
    # beta = (X^T * X)^-1 * X^T * y
    try:
        beta, residuals, rank, s = la.lstsq(X_mat, y_trans)
        
        # Calculate standard errors and p-values
        dof = n - 2
        if residuals.size > 0 and residuals[0] > 0:
            mse = residuals[0] / dof
        else:
            pred = X_mat @ beta
            mse = np.sum((y_trans - pred)**2) / dof
            
        var_beta = mse * la.inv(X_mat.T @ X_mat)
        se_beta = np.sqrt(np.diagonal(var_beta))
        
        b0, b1 = beta[0], beta[1]
        se_b0, se_b1 = se_beta[0], se_beta[1]
        
        t_stat = b1 / (se_b1 + 1e-10)
        p_val = 2 * (1 - stats.t.cdf(abs(t_stat), df=dof))
        
        return b0, b1, se_b1, t_stat, p_val
    except Exception as e:
        return np.nan, np.nan, np.nan, np.nan, np.nan

def manual_fdr(pvals, alpha=0.05):
    pvals = np.asfarray(pvals)
    n = len(pvals)
    by_descend = pvals.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = np.arange(n, 0, -1)
    qvals = np.minimum.accumulate(pvals[by_descend] * n / steps)
    qvals = qvals[by_orig]
    qvals[qvals > 1.0] = 1.0
    return qvals < alpha, qvals

def main():
    print("[*] Starting PGLS Phenotype Association scan...")
    if not os.path.exists(TREE_PATH):
        print(f"[!] Error: Tree file '{TREE_PATH}' not found.")
        sys.exit(1)
        
    # 1. Parse tree and get species names
    print(f"[*] Parsing Newick tree from '{TREE_PATH}'...")
    with open(TREE_PATH, 'r') as f:
        tree_str = f.read()
    
    tree_root = parse_newick(tree_str)
    paths = get_paths_from_root(tree_root)
    leaf_names = sorted(list(paths.keys()))
    n_species = len(leaf_names)
    print(f"    - Parsed {n_species} species from tree.")
    
    # 2. Simulate or load species phenotypes
    print("[*] Simulating/loading species body sizes...")
    phenotypes_df = simulate_phenotypes(paths, leaf_names)
    
    # 3. Precompute phylogenetic covariance and Cholesky decomposition
    print("[*] Precomputing phylogenetic covariance matrix V and Cholesky factor L...")
    V = compute_covariance_matrix(paths, leaf_names)
    # Regularization to guarantee positive definiteness
    V += np.eye(n_species) * 1e-4
    L = la.cholesky(V, lower=True)
    
    # 4. Query positive selection substitutions at leaf nodes from SQLite DB
    if not os.path.exists(DB_PATH):
        print(f"[!] Error: SQLite DB '{DB_PATH}' not found.")
        sys.exit(1)
        
    print(f"[*] Querying selection events from '{DB_PATH}'...")
    conn = sqlite3.connect(DB_PATH)
    
    sql = """
    SELECT ss.gene_name, bm.master_node_id, COUNT(*) as selection_substitutions
    FROM site_substitutions ss
    JOIN branch_mappings bm ON ss.gene_name = bm.gene_name AND ss.branch_name = bm.branch_name
    JOIN master_nodes mn ON bm.master_node_id = mn.node_id
    JOIN site_results sr ON ss.gene_name = sr.gene_name AND ss.site_index = sr.site_index
    WHERE ss.is_synonymous = 0
      AND mn.is_leaf = 1
      AND sr.p_value < 0.05
    GROUP BY ss.gene_name, bm.master_node_id
    """
    
    events_df = pd.read_sql_query(sql, conn)
    # Fetch list of all genes
    all_genes = pd.read_sql_query("SELECT DISTINCT gene_name FROM gene_results", conn)["gene_name"].tolist()
    conn.close()
    
    print(f"    - Found {len(events_df)} terminal branch selection records across {len(all_genes)} genes.")
    
    # Map species names from DB to tree leaves
    # Note: Leaf name (e.g. hg) matches the prefix before the pipe in DB (hg|Human (Homo sapiens))
    # Let's build a map from leaf name to DB node ID
    db_node_ids = pd.read_sql_query("SELECT node_id FROM master_nodes WHERE is_leaf = 1", sqlite3.connect(DB_PATH))["node_id"].tolist()
    
    leaf_to_db = {}
    for node_id in db_node_ids:
        prefix = node_id.split('|')[0]
        if prefix in leaf_names:
            leaf_to_db[prefix] = node_id
            
    # Keep only species that match both tree and database
    species_to_analyze = [l for l in leaf_names if l in leaf_to_db]
    n_analyzed = len(species_to_analyze)
    print(f"    - Analyzing {n_analyzed} species overlapping between database and tree.")
    
    # Subset phenotypes and precomputed L
    # We need to construct V and L specifically for the overlapping set
    paths_subset = {s: paths[s] for s in species_to_analyze}
    V_sub = compute_covariance_matrix(paths_subset, species_to_analyze)
    V_sub += np.eye(n_analyzed) * 1e-4
    L_sub = la.cholesky(V_sub, lower=True)
    
    # Sort phenotypes to match species_to_analyze order
    pheno_map = dict(zip(phenotypes_df["species_id"], phenotypes_df["log_body_mass"]))
    Y = np.array([pheno_map[s] for s in species_to_analyze])
    
    # 5. Run PGLS regression per gene
    print(f"[*] Executing PGLS association screens for {len(all_genes)} genes...")
    results = []
    
    # Group events by gene for fast lookup
    events_by_gene = {}
    for gene_name, group in events_df.groupby("gene_name"):
        # Map master_node_id to leaf name
        db_node_to_val = dict(zip(group["master_node_id"], group["selection_substitutions"]))
        events_by_gene[gene_name] = db_node_to_val
        
    for idx, gene in enumerate(all_genes):
        if idx % 1000 == 0 and idx > 0:
            print(f"    - Processed {idx}/{len(all_genes)} genes...")
            
        # Build selection vector X (count of selection events on terminal branches)
        # Defaults to 0 if no events
        gene_events = events_by_gene.get(gene, {})
        X = np.zeros(n_analyzed)
        for i, species in enumerate(species_to_analyze):
            db_node = leaf_to_db[species]
            X[i] = gene_events.get(db_node, 0.0)
            
        if np.sum(X) == 0:
            # Skip genes with zero selection events on analyzed terminal branches
            continue
            
        b0, b1, se_b1, t_stat, p_val = perform_pgls(Y, X, L_sub)
        
        if not np.isnan(p_val):
            results.append({
                "gene_name": gene,
                "selection_events": np.sum(X),
                "intercept_b0": round(b0, 4),
                "slope_b1": round(b1, 4),
                "se_b1": round(se_b1, 4),
                "t_statistic": round(t_stat, 4),
                "p_value": p_val
            })
            
    res_df = pd.DataFrame(results)
    
    if res_df.empty:
        print("[!] No valid regression results generated.")
        sys.exit(0)
        
    # 6. Apply multiple testing correction
    print("[*] Applying Benjamini-Hochberg FDR correction...")
    _, q_values = manual_fdr(res_df["p_value"].values, alpha=0.05)
    res_df["fdr_q_value"] = q_values
    
    # Sort by significance
    res_df = res_df.sort_values(by="p_value")
    
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    res_df.to_csv(OUTPUT_CSV, index=False)
    print(f"🎉 PGLS association results saved to '{OUTPUT_CSV}'")
    
    # 7. Print top results
    print("\n" + "="*80)
    print("🏆 Top 15 Genes Associated with Continuous Phenotype (Body Size) via PGLS")
    print("="*80)
    print(f"{'Rank':<5} {'Gene Name':<15} {'Sel Events':<12} {'Slope (b1)':<12} {'T-stat':<10} {'P-value':<10} {'FDR Q-value':<12}")
    print("-"*80)
    
    top_df = res_df.head(15)
    for r_idx, row in enumerate(top_df.itertuples()):
        print(f"{r_idx+1:<5} {row.gene_name:<15} {row.selection_events:<12.0f} {row.slope_b1:<12.4f} {row.t_statistic:<10.3f} {row.p_value:<10.2e} {row.fdr_q_value:<12.4f}")
    print("="*80)

if __name__ == "__main__":
    main()
