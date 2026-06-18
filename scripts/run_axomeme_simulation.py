#!/usr/bin/env python3
"""
run_axomeme_simulation.py
-------------------------
A script to test the AxoMeme transformer model (PhyloAxialTransformer) on 
controlled simulated sequence datasets generated using pyvolve.

Steps:
  1. Generate a balanced phylogenetic tree of specified leaf count and branch lengths.
  2. Define two codon partitions:
     - Partition 1: Purifying selection (dN/dS = omega_neg)
     - Partition 2: Positive selection (dN/dS = omega_pos)
  3. Simulate sequence evolution using pyvolve.
  4. Write the simulated alignment and tree as a single NEXUS file.
  5. Run predict_regression_nexus.py using the specified trained transformer model.
  6. Read and analyze the predicted Likelihood Ratio Test (LRT) statistics:
     - Compare predicted selection values at simulated positive selection sites vs purifying selection sites.
     - Calculate sensitivity (TPR) and false discovery rates.
"""

import os
import re
import sys
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
import pyvolve

def make_balanced_tree(n_leaves, branch_length=0.05):
    """Generates a balanced tree string with n_leaves (must be power of 2)."""
    leaves = [f"sp{i}" for i in range(1, n_leaves + 1)]
    while len(leaves) > 1:
        next_level = []
        for i in range(0, len(leaves), 2):
            next_level.append(f"({leaves[i]}:{branch_length},{leaves[i+1]}:{branch_length})")
        leaves = next_level
    return leaves[0] + ";"

def write_nexus(seq_dict, taxlabels, tree_str, output_path):
    """Writes a standard NEXUS alignment and embedded tree."""
    ntax = len(taxlabels)
    first_seq = next(iter(seq_dict.values()))
    nchar = len(first_seq)
    
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

def run_simulation(args, scratch_dir):
    """Runs pyvolve simulation and converts output to NEXUS."""
    tree_str = make_balanced_tree(args.num_taxa, args.branch_length)
    my_tree = pyvolve.read_tree(tree=tree_str)
    
    print(f"🌲 Generated balanced tree with {args.num_taxa} species (branch length {args.branch_length})")
    
    # 1. Define models
    print(f"🧬 Setting up codon models: purifying omega = {args.omega_neg}, positive omega = {args.omega_pos}")
    model_neg = pyvolve.Model("codon", {"omega": args.omega_neg, "kappa": 2.0})
    model_pos = pyvolve.Model("codon", {"omega": args.omega_pos, "kappa": 2.0})
    
    # 2. Define partitions
    part_neg = pyvolve.Partition(models=model_neg, size=args.num_neg_sites)
    part_pos = pyvolve.Partition(models=model_pos, size=args.num_pos_sites)
    
    # 3. Simulate
    print(f"🚦 Running pyvolve evolution (total codons: {args.num_neg_sites + args.num_pos_sites})...")
    fasta_path = os.path.join(scratch_dir, "simulated_raw.fasta")
    
    # Evolver parameters: output fasta
    evolver = pyvolve.Evolver(tree=my_tree, partitions=[part_neg, part_pos])
    evolver(seqfile=fasta_path, write_joint_states=False, infofile=None)
    
    # 4. Parse FASTA and build NEXUS
    seq_dict = {}
    taxlabels = [f"sp{i}" for i in range(1, args.num_taxa + 1)]
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[record.id] = str(record.seq).upper()
        
    nexus_path = os.path.join(scratch_dir, "simulated_alignment.nex")
    write_nexus(seq_dict, taxlabels, tree_str, nexus_path)
    print(f"💾 NEXUS alignment and tree saved to: {nexus_path}")
    return nexus_path

def run_inference(args, nexus_path, scratch_dir):
    """Executes the predict_regression_nexus.py driver script."""
    print("🧠 Running AxoMeme transformer inference...")
    pred_csv_path = os.path.join(scratch_dir, "predictions.csv")
    
    cmd = [
        "python3", "scripts/predict_regression_nexus.py",
        "--alignment", nexus_path,
        "--model", args.model,
        "--output", pred_csv_path,
        "--tier1_lrt_gate", str(args.tier1_lrt_gate),
        "--tier2_lrt_gate", str(args.tier2_lrt_gate)
    ]
    
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"❌ Inference failed with exit code {res.returncode}")
        print(f"Stderr:\n{res.stderr}")
        sys.exit(1)
        
    print("✅ Inference complete.")
    return pred_csv_path

def analyze_predictions(args, pred_csv):
    """Analyzes the predictions CSV file and computes performance metrics."""
    df = pd.read_csv(pred_csv)
    
    # Split predictions back into purifying vs positive selection partitions
    df_neg = df.iloc[:args.num_neg_sites]
    df_pos = df.iloc[args.num_neg_sites:]
    
    mean_neg_lrt = df_neg["predicted_lrt"].mean()
    max_neg_lrt = df_neg["predicted_lrt"].max()
    
    mean_pos_lrt = df_pos["predicted_lrt"].mean()
    min_pos_lrt = df_pos["predicted_lrt"].min()
    max_pos_lrt = df_pos["predicted_lrt"].max()
    
    # Count true positives (TP), false positives (FP), etc.
    # Positive sites are in df_pos, negative sites in df_neg
    tp = len(df_pos[df_pos["selection_call"] != "Neutral"])
    fn = len(df_pos[df_pos["selection_call"] == "Neutral"])
    
    fp = len(df_neg[df_neg["selection_call"] != "Neutral"])
    tn = len(df_neg[df_neg["selection_call"] == "Neutral"])
    
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    fdr = fp / (fp + tp) if (fp + tp) > 0 else 0.0
    
    print("\n" + "=" * 60)
    print("📊 AxoMeme Simulation Performance Analysis")
    print("=" * 60)
    print(f"Negative Partition (dN/dS = {args.omega_neg}, {args.num_neg_sites} sites):")
    print(f"  - Mean predicted LRT:       {mean_neg_lrt:.4f}")
    print(f"  - Max predicted LRT:        {max_neg_lrt:.4f}")
    print(f"  - Called False Positives:   {fp} / {args.num_neg_sites} (FPR = {fpr*100:.2f}%)")
    print("")
    print(f"Positive Partition (dN/dS = {args.omega_pos}, {args.num_pos_sites} sites):")
    print(f"  - Mean predicted LRT:       {mean_pos_lrt:.4f}")
    print(f"  - Min predicted LRT:        {min_pos_lrt:.4f}")
    print(f"  - Max predicted LRT:        {max_pos_lrt:.4f}")
    print(f"  - Called True Positives:    {tp} / {args.num_pos_sites} (Sensitivity/TPR = {tpr*100:.2f}%)")
    print("")
    print(f"Global Statistics:")
    print(f"  - True Positive Rate:       {tpr:.4f}")
    print(f"  - False Positive Rate:      {fpr:.4f}")
    print(f"  - False Discovery Rate:     {fdr:.4f}")
    print("=" * 60)
    
    # Print the positive selection sites' predictions
    print(f"\nDetailed predictions for simulated positive selection sites (Sites {args.num_neg_sites + 1}-{args.num_neg_sites + args.num_pos_sites}):")
    cols = ["codon_site", "ref_codon", "ref_aa", "predicted_lrt", "selection_call"]
    print(df_pos[cols].to_string(index=False))

def main():
    parser = argparse.ArgumentParser(description="Test AxoMeme on pyvolve simulated alignments.")
    parser.add_argument("--model", type=str, default="/Users/sergei/Projects/TOGA_MEME/MEME_transformer_joint.pt",
                        help="Path to trained transformer model weights.")
    parser.add_argument("--num_taxa", type=int, default=64, help="Number of species in tree (must be power of 2)")
    parser.add_argument("--branch_length", type=float, default=0.10, help="Branch length in simulated tree")
    parser.add_argument("--omega_neg", type=float, default=0.1, help="dN/dS for negative selection partition")
    parser.add_argument("--omega_pos", type=float, default=4.0, help="dN/dS for positive selection partition")
    parser.add_argument("--num_neg_sites", type=int, default=90, help="Number of negative selection codons")
    parser.add_argument("--num_pos_sites", type=int, default=10, help="Number of positive selection codons")
    parser.add_argument("--tier1_lrt_gate", type=float, default=5.0, help="LRT gate for Tier 1 calls")
    parser.add_argument("--tier2_lrt_gate", type=float, default=3.0, help="LRT gate for Tier 2 calls")
    parser.add_argument("--scratch_dir", type=str, default="/Users/sergei/.gemini/antigravity-cli/brain/b2e83fd0-f7f4-4c8a-9fe4-ae9d2217a661/scratch",
                        help="Folder to store temporary simulation files.")
    
    args = parser.parse_args()
    
    # Ensure scratch directory exists
    os.makedirs(args.scratch_dir, exist_ok=True)
    
    # 1. Run simulation
    nexus_path = run_simulation(args, args.scratch_dir)
    
    # 2. Run inference
    pred_csv = run_inference(args, nexus_path, args.scratch_dir)
    
    # 3. Analyze results
    analyze_predictions(args, pred_csv)

if __name__ == "__main__":
    main()
