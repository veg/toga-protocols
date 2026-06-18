#!/usr/bin/env python3
"""
run_joint_comparison.py
-----------------------
A script to perform a side-by-side comparison of HyPhy MEME (Maximum Likelihood)
and AxoMeme (PhyloAxialTransformer, deep learning) on the same simulated dataset.

Steps:
  1. Simulates evolution along a 64-taxa tree using pyvolve (90 purifying sites, 10 positive selection sites).
  2. Converts to NEXUS format.
  3. Runs predict_regression_nexus.py (AxoMeme inference).
  4. Runs hyphy meme (HyPhy MLE inference).
  5. Parses the MEME JSON results and matches them side-by-side with AxoMeme's CSV.
  6. Computes statistics (correlation, TPR, FPR) and prints a side-by-side validation table.
"""

import os
import sys
import json
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
import pyvolve

def make_balanced_tree(n_leaves, branch_length=0.05):
    leaves = [f"sp{i}" for i in range(1, n_leaves + 1)]
    while len(leaves) > 1:
        next_level = []
        for i in range(0, len(leaves), 2):
            next_level.append(f"({leaves[i]}:{branch_length},{leaves[i+1]}:{branch_length})")
        leaves = next_level
    return leaves[0] + ";"

def write_nexus(seq_dict, taxlabels, tree_str, output_path):
    ntax = len(taxlabels)
    nchar = len(next(iter(seq_dict.values())))
    with open(output_path, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN TAXA;\n")
        f.write(f"  DIMENSIONS NTAX={ntax};\n")
        f.write("  TAXLABELS\n")
        f.write("    " + " ".join(taxlabels) + "\n")
        f.write("  ;\n")
        f.write("END;\n\n")
        
        f.write("BEGIN CHARACTERS;\n")
        f.write(f"  DIMENSIONS NCHAR={nchar};\n")
        f.write("  FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
        f.write("  MATRIX\n")
        for tax in taxlabels:
            f.write(f"    {tax} {seq_dict[tax]}\n")
        f.write("  ;\n")
        f.write("END;\n\n")
        
        f.write("BEGIN TREES;\n")
        f.write(f"  TREE tree_1 = {tree_str}\n")
        f.write("END;\n")

def simulate(args, scratch_dir):
    tree_str = make_balanced_tree(args.num_taxa, args.branch_length)
    my_tree = pyvolve.read_tree(tree=tree_str)
    
    model_neg = pyvolve.Model("codon", {"omega": args.omega_neg, "kappa": 2.0})
    model_pos = pyvolve.Model("codon", {"omega": args.omega_pos, "kappa": 2.0})
    
    part_neg = pyvolve.Partition(models=model_neg, size=args.num_neg_sites)
    part_pos = pyvolve.Partition(models=model_pos, size=args.num_pos_sites)
    
    fasta_path = os.path.join(scratch_dir, "joint_sim.fasta")
    evolver = pyvolve.Evolver(tree=my_tree, partitions=[part_neg, part_pos])
    evolver(seqfile=fasta_path, write_joint_states=False, infofile=None)
    
    seq_dict = {}
    taxlabels = [f"sp{i}" for i in range(1, args.num_taxa + 1)]
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[record.id] = str(record.seq).upper()
        
    nexus_path = os.path.join(scratch_dir, "joint_sim.nex")
    write_nexus(seq_dict, taxlabels, tree_str, nexus_path)
    return nexus_path, tree_str

def run_axomeme(args, nexus_path, scratch_dir):
    pred_csv_path = os.path.join(scratch_dir, "joint_preds.csv")
    cmd = [
        "python3", "scripts/predict_regression_nexus.py",
        "--alignment", nexus_path,
        "--model", args.model,
        "--output", pred_csv_path,
        "--tier1_lrt_gate", str(args.tier1_lrt_gate),
        "--tier2_lrt_gate", str(args.tier2_lrt_gate)
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return pred_csv_path

def run_hyphy_meme(nexus_path, scratch_dir):
    json_path = os.path.join(scratch_dir, "joint_meme.json")
    cmd = [
        "hyphy", "meme",
        "--alignment", nexus_path,
        "--branches", "All",
        "--pvalue", "0.1",
        "--output", json_path
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json_path

def parse_meme_json(json_path, num_sites):
    with open(json_path, 'r') as f:
        data = json.load(f)
    mle_content = data.get('MLE', {}).get('content', {}).get('0', [])
    
    meme_results = []
    for s_idx in range(num_sites):
        if s_idx < len(mle_content):
            row = mle_content[s_idx]
            meme_results.append({
                "site": s_idx + 1,
                "meme_lrt": row[5],
                "meme_pval": row[6]
            })
        else:
            meme_results.append({
                "site": s_idx + 1,
                "meme_lrt": 0.0,
                "meme_pval": 1.0
            })
    return pd.DataFrame(meme_results)

def main():
    parser = argparse.ArgumentParser(description="Run side-by-side MEME and AxoMeme validation comparison.")
    parser.add_argument("--model", type=str, default="/Users/sergei/Projects/TOGA_MEME/MEME_transformer_joint.pt",
                        help="Path to trained model weights.")
    parser.add_argument("--num_taxa", type=int, default=64, help="Number of taxa in simulation tree.")
    parser.add_argument("--branch_length", type=float, default=0.10, help="Branch lengths.")
    parser.add_argument("--omega_neg", type=float, default=0.1, help="dN/dS for negative selection.")
    parser.add_argument("--omega_pos", type=float, default=4.0, help="dN/dS for positive selection.")
    parser.add_argument("--num_neg_sites", type=int, default=90, help="Number of purifying codons.")
    parser.add_argument("--num_pos_sites", type=int, default=10, help="Number of positive selection codons.")
    parser.add_argument("--tier1_lrt_gate", type=float, default=5.0, help="LRT gate for Tier 1 calls")
    parser.add_argument("--tier2_lrt_gate", type=float, default=3.0, help="LRT gate for Tier 2 calls")
    parser.add_argument("--scratch_dir", type=str, default="/Users/sergei/.gemini/antigravity-cli/brain/b2e83fd0-f7f4-4c8a-9fe4-ae9d2217a661/scratch",
                        help="Folder to store temporary files.")
    
    args = parser.parse_args()
    os.makedirs(args.scratch_dir, exist_ok=True)
    total_sites = args.num_neg_sites + args.num_pos_sites
    
    # 1. Run simulation
    nexus_path, _ = simulate(args, args.scratch_dir)
    print(f"🌲 Evolved dataset simulated: {args.num_taxa} species, {total_sites} codons.")
    
    # 2. Run AxoMeme (DL model)
    pred_csv = run_axomeme(args, nexus_path, args.scratch_dir)
    df_axo = pd.read_csv(pred_csv)
    
    # 3. Run HyPhy MEME (MLE)
    json_path = run_hyphy_meme(nexus_path, args.scratch_dir)
    df_meme = parse_meme_json(json_path, total_sites)
    
    # 4. Merge results
    df_axo = df_axo.rename(columns={"codon_site": "site"})
    df_comparison = pd.merge(df_axo, df_meme, on="site")
    
    # Add simulated truth labels
    df_comparison["true_omega"] = args.omega_neg
    df_comparison.loc[df_comparison["site"] > args.num_neg_sites, "true_omega"] = args.omega_pos
    
    # Calculate correlation
    pearson = df_comparison["predicted_lrt"].corr(df_comparison["meme_lrt"], method="pearson")
    spearman = df_comparison["predicted_lrt"].corr(df_comparison["meme_lrt"], method="spearman")
    
    # Print results
    print("\n" + "=" * 80)
    print("🥊 Side-by-Side Comparison: HyPhy MEME vs. AxoMeme")
    print("=" * 80)
    print(f"LRT Prediction Correlation (n={total_sites} sites):")
    print(f"  - Pearson Correlation (R):     {pearson:.4f}")
    print(f"  - Spearman Rank Correlation:   {spearman:.4f}")
    
    print("\n" + "-" * 80)
    print(f"Positive Selection Partition Comparison (Sites {args.num_neg_sites + 1}-{total_sites}, Evolved Omega = {args.omega_pos}):")
    print("-" * 80)
    
    df_pos = df_comparison.iloc[args.num_neg_sites:]
    print_cols = ["site", "ref_codon", "ref_aa", "meme_lrt", "meme_pval", "predicted_lrt", "selection_call"]
    print(df_pos[print_cols].to_string(index=False))
    
    print("\n" + "-" * 80)
    print(f"False Positives in Purifying Partition (Sites 1-{args.num_neg_sites}, Evolved Omega = {args.omega_neg}):")
    print("-" * 80)
    
    df_neg = df_comparison.iloc[:args.num_neg_sites]
    # FP is if MEME pval <= 0.1 OR AxoMeme selection_call != "Neutral"
    fp_meme = df_neg[df_neg["meme_pval"] <= 0.1]
    fp_axo = df_neg[df_neg["selection_call"] != "Neutral"]
    
    fp_all = df_neg[(df_neg["meme_pval"] <= 0.1) | (df_neg["selection_call"] != "Neutral")]
    if len(fp_all) > 0:
        print(fp_all[print_cols].to_string(index=False))
    else:
        print("No false positives called by either model (FPR = 0.0%).")
    
    print("=" * 80)
    print(f"Summary of positive selection calls at p <= 0.1 (MEME) / Tier 1-2 (AxoMeme):")
    print(f"  - HyPhy MEME true positives:   {len(df_pos[df_pos['meme_pval'] <= 0.1])} / {args.num_pos_sites}")
    print(f"  - AxoMeme true positives:      {len(df_pos[df_pos['selection_call'] != 'Neutral'])} / {args.num_pos_sites}")
    print(f"  - HyPhy MEME false positives:  {len(fp_meme)} / {args.num_neg_sites}")
    print(f"  - AxoMeme false positives:     {len(fp_axo)} / {args.num_neg_sites}")
    print("=" * 80)

if __name__ == "__main__":
    main()
