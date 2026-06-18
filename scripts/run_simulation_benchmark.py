#!/usr/bin/env python3
"""
run_simulation_benchmark.py
---------------------------
A comprehensive benchmarking pipeline for evaluating the AxoMeme selection model
(PhyloAxialTransformer) across 102 selective scenarios (null, homogeneous, episodic,
and misspecified MutSel) using tree topologies from HyPhy test data.

Generates:
  - SQLite database: sim_benchmark_results.db
  - Dynamic HTML dashboard: toga_protocols/reports/simulation_dashboard.html (synced to Dropbox)
"""

import os
import re
import sys
import json
import sqlite3
import argparse
import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from Bio import SeqIO, Phylo
import pyvolve
import gzip

# Alphabetical order of Pyvolve amino acids
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def sanitize_name(name):
    """Removes underscores from names to prevent Pyvolve flag parsing bugs."""
    return name.replace("_", "")

def get_biophysical_fitness(property_type, select_coef=4.0):
    """Generates 20-element amino acid fitness profiles."""
    fitness = [0.0] * 20
    if property_type == 'hydrophobic':
        targets = set("AFGILMPVWY")
    elif property_type == 'charged':
        targets = set("DEHKR")
    elif property_type == 'polar':
        targets = set("CNQST")
    elif property_type == 'small':
        targets = set("AGPS")
    else:
        targets = set()
        
    for i, aa in enumerate(AA_LIST):
        if aa in targets:
            fitness[i] = select_coef
    return fitness

def to_newick(node, target_names_to_flag, flag_str="#m1"):
    """Custom recursive Newick tree formatter that applies branch-specific model flags."""
    name_str = node.name if node.name else ""
    bl_str = f":{node.branch_length:.6f}" if node.branch_length is not None else ""
    if node.name and node.name in target_names_to_flag:
        bl_str = f"{bl_str}{flag_str}"
        
    if node.is_terminal():
        return f"{name_str}{bl_str}"
    else:
        children_str = ",".join(to_newick(child, target_names_to_flag, flag_str) for child in node.clades)
        return f"({children_str}){name_str}{bl_str}"

def write_nexus(seq_dict, taxlabels, tree_str, output_path):
    """Writes a standard NEXUS alignment and embedded tree (supports gzip if output_path ends with .gz)."""
    ntax = len(taxlabels)
    first_seq = next(iter(seq_dict.values()))
    nchar = len(first_seq)
    
    open_func = gzip.open if output_path.endswith('.gz') else open
    mode = "wt" if output_path.endswith('.gz') else "w"
    
    with open_func(output_path, mode) as f:
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

def setup_db(db_path):
    """Initializes the SQLite database."""
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("""
        CREATE TABLE IF NOT EXISTS results (
            scenario_id INTEGER,
            replicate INTEGER,
            tree_name TEXT,
            regime TEXT,
            model_type TEXT,
            omega_neg REAL,
            omega_pos REAL,
            num_neg_sites INTEGER,
            num_pos_sites INTEGER,
            branch_setting TEXT,
            tpr REAL,
            fpr REAL,
            fdr REAL,
            mean_neg_lrt REAL,
            mean_pos_lrt REAL,
            PRIMARY KEY (scenario_id, replicate)
        )
    """)
    conn.commit()
    conn.close()

def save_result_to_db(db_path, result_dict):
    """Saves a single replicate run result to the database."""
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("""
        INSERT OR REPLACE INTO results 
        (scenario_id, replicate, tree_name, regime, model_type, omega_neg, omega_pos, 
         num_neg_sites, num_pos_sites, branch_setting, tpr, fpr, fdr, mean_neg_lrt, mean_pos_lrt)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        result_dict['scenario_id'], result_dict['replicate'], result_dict['tree_name'],
        result_dict['regime'], result_dict['model_type'], result_dict['omega_neg'],
        result_dict['omega_pos'], result_dict['num_neg_sites'], result_dict['num_pos_sites'],
        result_dict['branch_setting'], result_dict['tpr'], result_dict['fpr'],
        result_dict['fdr'], result_dict['mean_neg_lrt'], result_dict['mean_pos_lrt']
    ))
    conn.commit()
    conn.close()

def generate_dashboard(db_path, report_path, total_scenarios, current_run):
    """Generates the static HTML progress and results dashboard."""
    if not os.path.exists(db_path):
        return
        
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM results", conn)
    conn.close()
    
    # Calculate aggregate statistics
    total_reps = len(df)
    if total_reps > 0:
        overall_tpr = df["tpr"].mean()
        overall_fpr = df["fpr"].mean()
        overall_fdr = df["fdr"].mean()
    else:
        overall_tpr, overall_fpr, overall_fdr = 0.0, 0.0, 0.0
        
    # Group by scenario to build the table
    summary_list = []
    if total_reps > 0:
        groups = df.groupby(["scenario_id", "tree_name", "regime", "model_type", "omega_pos", "branch_setting"])
        for (scen_id, tree_name, regime, m_type, o_pos, b_set), g in groups:
            summary_list.append({
                "scenario_id": int(scen_id),
                "tree_name": tree_name,
                "regime": regime,
                "model_type": m_type,
                "omega_pos": float(o_pos),
                "branch_setting": b_set,
                "reps_completed": len(g),
                "mean_tpr": float(g["tpr"].mean()),
                "mean_fpr": float(g["fpr"].mean()),
                "mean_neg_lrt": float(g["mean_neg_lrt"].mean()),
                "mean_pos_lrt": float(g["mean_pos_lrt"].mean())
            })
    if total_reps > 0:
        summary_df = pd.DataFrame(summary_list).sort_values(by="scenario_id")
    else:
        summary_df = pd.DataFrame(columns=["scenario_id", "tree_name", "regime", "model_type", "omega_pos", "branch_setting", "reps_completed", "mean_tpr", "mean_fpr", "mean_neg_lrt", "mean_pos_lrt"])
    summary_json = summary_df.to_json(orient="records")
    
    # Progress Calculation
    progress_pct = (current_run / (total_scenarios * 100)) * 100
    
    # HTML Content
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>AxoMeme Simulation Benchmark Dashboard</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        :root {{
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --accent: #38bdf8;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --border: #334155;
            --success: #10b981;
            --warning: #f59e0b;
        }}
        body {{
            background-color: var(--bg-color);
            color: var(--text-primary);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--border);
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        h1 {{ margin: 0; font-size: 28px; color: var(--text-primary); }}
        .badge {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 14px;
            color: var(--accent);
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 20px;
            text-align: center;
        }}
        .card-title {{ font-size: 14px; color: var(--text-secondary); margin-bottom: 10px; }}
        .card-value {{ font-size: 24px; font-weight: bold; color: var(--accent); }}
        .progress-bar {{
            background-color: var(--border);
            height: 10px;
            border-radius: 5px;
            overflow: hidden;
            margin-top: 10px;
        }}
        .progress-fill {{
            background-color: var(--success);
            height: 100%;
            width: {progress_pct:.2f}%;
            transition: width 0.3s;
        }}
        .charts {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(450px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .chart-card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 20px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background-color: var(--card-bg);
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid var(--border);
            margin-top: 20px;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}
        th {{
            background-color: #1e293b;
            color: var(--text-secondary);
            font-size: 12px;
            text-transform: uppercase;
        }}
        tr:hover {{ background-color: #273549; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div>
                <h1>AxoMeme Selection Simulation Benchmark</h1>
                <p style="color: var(--text-secondary); margin: 5px 0 0 0;">Evaluating model sensitivity and false-positive rates on 10,200 datasets</p>
            </div>
            <div class="badge">Running scenarios ({current_run}/{total_scenarios * 100} reps)</div>
        </header>
        
        <div class="grid">
            <div class="card" style="grid-column: span 2;">
                <div class="card-title">Overall Progress ({progress_pct:.1f}%)</div>
                <div class="card-value">{current_run} / {total_scenarios * 100} Replicates</div>
                <div class="progress-bar"><div class="progress-fill"></div></div>
            </div>
            <div class="card">
                <div class="card-title">Average Sensitivity (TPR)</div>
                <div class="card-value" style="color: var(--success);">{overall_tpr*100:.2f}%</div>
            </div>
            <div class="card">
                <div class="card-title">False Positive Rate (FPR)</div>
                <div class="card-value" style="color: var(--warning);">{overall_fpr*100:.2f}%</div>
            </div>
            <div class="card">
                <div class="card-title">False Discovery Rate (FDR)</div>
                <div class="card-value">{overall_fdr*100:.2f}%</div>
            </div>
        </div>

        <div class="charts">
            <div class="chart-card">
                <h3 style="margin-top:0;">Sensitivity by Evolutionary Backdrop</h3>
                <canvas id="backdropChart"></canvas>
            </div>
            <div class="chart-card">
                <h3 style="margin-top:0;">FPR by Selective Regime</h3>
                <canvas id="regimeChart"></canvas>
            </div>
        </div>

        <h2>Scenario Summary Details</h2>
        <div style="overflow-x: auto;">
            <table id="summaryTable">
                <thead>
                    <tr>
                        <th>ID</th>
                        <th>Tree Backdrop</th>
                        <th>Regime</th>
                        <th>Model</th>
                        <th>Omega (Pos)</th>
                        <th>Branches</th>
                        <th>Reps</th>
                        <th>TPR</th>
                        <th>FPR</th>
                        <th>Mean Neg LRT</th>
                        <th>Mean Pos LRT</th>
                    </tr>
                </thead>
                <tbody>
                </tbody>
            </table>
        </div>
    </div>

    <script>
        const data = {summary_json};
        
        // Populate table
        const tbody = document.querySelector("#summaryTable tbody");
        data.forEach(row => {{
            const tr = document.createElement("tr");
            tr.innerHTML = `
                <td>${{row.scenario_id}}</td>
                <td>${{row.tree_name}}</td>
                <td>${{row.regime}}</td>
                <td>${{row.model_type}}</td>
                <td>${{row.omega_pos}}</td>
                <td>${{row.branch_setting}}</td>
                <td>${{row.reps_completed}}</td>
                <td style="color: var(--success); font-weight: bold;">${{(row.mean_tpr*100).toFixed(1)}}%</td>
                <td style="color: var(--warning);">${{(row.mean_fpr*100).toFixed(1)}}%</td>
                <td>${{row.mean_neg_lrt.toFixed(3)}}</td>
                <td>${{row.mean_pos_lrt.toFixed(3)}}</td>
            `;
            tbody.appendChild(tr);
        }});

        // Backdrop Chart
        const backdrops = [...new Set(data.map(d => d.tree_name))];
        const tprByBackdrop = backdrops.map(b => {{
            const subset = data.filter(d => d.tree_name === b && d.regime !== "Null");
            return subset.length ? (subset.reduce((acc, curr) => acc + curr.mean_tpr, 0) / subset.length) * 100 : 0;
        }});

        new Chart(document.getElementById('backdropChart'), {{
            type: 'bar',
            data: {{
                labels: backdrops,
                datasets: [{{
                    label: 'Mean Sensitivity (TPR %)',
                    data: tprByBackdrop,
                    backgroundColor: '#38bdf8',
                    borderWidth: 0,
                    borderRadius: 4
                }}]
            }},
            options: {{
                responsive: true,
                scales: {{ y: {{ min: 0, max: 100 }} }}
            }}
        }});

        // Regime Chart
        const regimes = [...new Set(data.map(d => d.regime))];
        const fprByRegime = regimes.map(r => {{
            const subset = data.filter(d => d.regime === r);
            return subset.length ? (subset.reduce((acc, curr) => acc + curr.mean_fpr, 0) / subset.length) * 100 : 0;
        }});

        new Chart(document.getElementById('regimeChart'), {{
            type: 'bar',
            data: {{
                labels: regimes,
                datasets: [{{
                    label: 'Mean False Positive Rate (FPR %)',
                    data: fprByRegime,
                    backgroundColor: '#f59e0b',
                    borderWidth: 0,
                    borderRadius: 4
                }}]
            }},
            options: {{
                responsive: true,
                scales: {{ y: {{ min: 0, max: 100 }} }}
            }}
        }});
    </script>
</body>
</html>
"""
    
    with open(report_path, "w") as f:
        f.write(html)
        
    # Copy to Dropbox
    dropbox_report_path = "/Users/sergei/Dropbox/TOGA2026/toga_protocols/reports/simulation_dashboard.html"
    os.makedirs(os.path.dirname(dropbox_report_path), exist_ok=True)
    subprocess.run(["cp", report_path, dropbox_report_path])

def setup_scenarios():
    """Generates the list of 102 benchmark scenarios (25 per tree backdrop + extras)."""
    backdrops = [
        {"name": "bglobin.nex", "num_neg": 270, "num_pos": 30},
        {"name": "lysin.nex", "num_neg": 270, "num_pos": 30},
        {"name": "camelid.nex", "num_neg": 270, "num_pos": 30},
        {"name": "pruned_species_tree.nhx", "num_neg": 450, "num_pos": 50}
    ]
    
    scenarios = []
    scen_id = 1
    
    for bd in backdrops:
        # 1. Regime I: Nulls (3 scenarios)
        scenarios.append({
            "scenario_id": scen_id, "tree_name": bd["name"], "regime": "Null", "model_type": "GY94",
            "omega_neg": 0.05, "omega_pos": 0.05, "num_neg_sites": bd["num_neg"] + bd["num_pos"], "num_pos_sites": 0, "branch_setting": "None"
        })
        scen_id += 1
        scenarios.append({
            "scenario_id": scen_id, "tree_name": bd["name"], "regime": "Null", "model_type": "GY94",
            "omega_neg": 0.20, "omega_pos": 0.20, "num_neg_sites": bd["num_neg"] + bd["num_pos"], "num_pos_sites": 0, "branch_setting": "None"
        })
        scen_id += 1
        scenarios.append({
            "scenario_id": scen_id, "tree_name": bd["name"], "regime": "Null", "model_type": "GY94",
            "omega_neg": 1.0, "omega_pos": 1.0, "num_neg_sites": bd["num_neg"] + bd["num_pos"], "num_pos_sites": 0, "branch_setting": "None"
        })
        scen_id += 1
        
        # 2. Regime II: Constant (6 scenarios)
        for w_pos in [2.5, 5.0]:
            for pos_frac in [0.05, 0.10, 0.20]:
                pos_sites = int((bd["num_neg"] + bd["num_pos"]) * pos_frac)
                neg_sites = (bd["num_neg"] + bd["num_pos"]) - pos_sites
                scenarios.append({
                    "scenario_id": scen_id, "tree_name": bd["name"], "regime": "Constant", "model_type": "GY94",
                    "omega_neg": 0.1, "omega_pos": w_pos, "num_neg_sites": neg_sites, "num_pos_sites": pos_sites, "branch_setting": "All"
                })
                scen_id += 1
                
        # 3. Regime III: Episodic (9 scenarios)
        for b_set in ["internal", "terminal", "clade"]:
            for w_pos in [3.0, 6.0, 12.0]:
                scenarios.append({
                    "scenario_id": scen_id, "tree_name": bd["name"], "regime": "Episodic", "model_type": "Branch-site",
                    "omega_neg": 0.1, "omega_pos": w_pos, "num_neg_sites": bd["num_neg"], "num_pos_sites": bd["num_pos"], "branch_setting": b_set
                })
                scen_id += 1
                
        # 4. Regime IV: MutSel (7 scenarios)
        # IV-A: Pure MutSel (4 scenarios)
        for prop in ["hydrophobic", "charged", "polar", "small"]:
            scenarios.append({
                "scenario_id": scen_id, "tree_name": bd["name"], "regime": "MutSel", "model_type": f"MutSel-{prop}",
                "omega_neg": 0.0, "omega_pos": 0.0, "num_neg_sites": bd["num_neg"] + bd["num_pos"], "num_pos_sites": 0, "branch_setting": "All"
            })
            scen_id += 1
        # IV-B: Hybrid MutSel + Episodic Clade (3 scenarios)
        for prop in ["hydrophobic", "charged", "polar"]:
            scenarios.append({
                "scenario_id": scen_id, "tree_name": bd["name"], "regime": "Hybrid-MutSel", "model_type": f"MutSel-{prop}+GY94",
                "omega_neg": 0.0, "omega_pos": 4.0, "num_neg_sites": bd["num_neg"], "num_pos_sites": bd["num_pos"], "branch_setting": "clade"
            })
            scen_id += 1
            
    return scenarios

def select_branches_for_flagging(tree, setting, tree_name):
    """Programmatically selects specific non-random branches to experience positive selection."""
    terminals = tree.get_terminals()
    internals = [n for n in tree.get_nonterminals() if n != tree.root]
    
    # Ensure they have names
    for i, leaf in enumerate(terminals):
        if not leaf.name:
            leaf.name = f"leaf{i}"
        leaf.name = sanitize_name(leaf.name)
    for i, node in enumerate(tree.get_nonterminals()):
        if not node.name:
            node.name = f"node{i}"
        node.name = sanitize_name(node.name)
        
    target_names = []
    
    if setting == "internal" and len(internals) > 0:
        # Take the top 10% longest internal branches
        internals.sort(key=lambda x: x.branch_length if x.branch_length is not None else 0.0, reverse=True)
        K = max(1, int(len(internals) * 0.10))
        target_names = [n.name for n in internals[:K]]
        
    elif setting == "terminal" and len(terminals) > 0:
        # Take the top 10% longest terminal branches
        terminals.sort(key=lambda x: x.branch_length if x.branch_length is not None else 0.0, reverse=True)
        K = max(1, int(len(terminals) * 0.10))
        target_names = [n.name for n in terminals[:K]]
        
    elif setting == "clade":
        # Select an internal node containing 10% to 30% of total leaves
        total_leaves = len(terminals)
        clade_node = None
        for node in internals:
            leaves_in_clade = len(node.get_terminals())
            if int(total_leaves * 0.10) <= leaves_in_clade <= int(total_leaves * 0.30):
                clade_node = node
                break
        if clade_node is None and len(internals) > 0:
            # Fallback to the first internal node with at least 2 leaves
            for node in internals:
                if len(node.get_terminals()) >= 2:
                    clade_node = node
                    break
        if clade_node:
            target_names = [n.name for n in clade_node.find_clades() if n.name]
            
    return target_names

def estimate_branch_lengths(alignment_path, tree_obj):
    """Runs HyPhy to estimate branch lengths of tree_obj using the original alignment."""
    # Write temporary alignment FASTA containing the species
    os.makedirs("scratch", exist_ok=True)
    temp_fasta = "scratch/temp_hyphy_align_est.fa"
    
    from predict_regression_nexus import parse_nexus_alignment_and_embedded_tree
    seq_dict, taxlabels, _ = parse_nexus_alignment_and_embedded_tree(alignment_path)
    
    with open(temp_fasta, "w") as f:
        for spec in taxlabels:
            if spec in seq_dict:
                f.write(f">{spec}\n{seq_dict[spec]}\n")
                
    # Format tree topology Newick string (removing any branch lengths)
    out_stream = StringIO()
    Phylo.write(tree_obj, out_stream, 'newick')
    raw_tree_str = out_stream.getvalue().strip()
    pruned_tree_str = re.sub(r':[0-9.eE-]+', '', raw_tree_str).strip().rstrip(';') + ';'
    
    temp_bf = "scratch/temp_hyphy_est_sim.bf"
    bf_content = f"""
DataSet ds = ReadDataFile("{temp_fasta}");
DataSetFilter df = CreateFilter(ds, 1);
HarvestFrequencies(freqs, df, 1, 1, 1);
global kappa = 1.0;
HKY85RateMatrix = [
    [*, kappa*t, t, kappa*t]
    [kappa*t, *, kappa*t, t]
    [t, kappa*t, *, kappa*t]
    [kappa*t, t, kappa*t, *]
];
Model HKY85Model = (HKY85RateMatrix, freqs);
UseModel(HKY85Model);
Tree T = "{pruned_tree_str}";
LikelihoodFunction lf = (df, T);
Optimize(res, lf);
fprintf(stdout, Format(T, 1, 1));
"""
    with open(temp_bf, "w") as f:
        f.write(bf_content)
        
    res = subprocess.run(["hyphy", temp_bf], capture_output=True, text=True)
    
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)
    if os.path.exists(temp_bf):
        os.remove(temp_bf)
        
    if res.returncode == 0 and res.stdout.strip():
        estimated_tree_str = res.stdout.strip()
        if not estimated_tree_str.endswith(';'):
            estimated_tree_str += ';'
        new_tree = Phylo.read(StringIO(estimated_tree_str), 'newick')
        return new_tree
    else:
        print(f"[!] Warning: HyPhy branch length estimation failed: {res.stderr}")
        return None

def simulate_only(sc, rep_idx, scratch_dir, replicates_dir):
    """Simulates a single replicate dataset and saves the gzipped NEXUS file."""
    # 1. Load Tree and Alignment Metadata
    if sc['tree_name'].endswith('.nhx'):
        # It is a tree file. Load tree directly and get leaf names.
        raw_path = f"/Users/sergei/Projects/TOGA_MEME/toga_protocols/data/{sc['tree_name']}"
        tree = Phylo.read(raw_path, 'newick')
        original_taxlabels = [leaf.name for leaf in tree.get_terminals() if leaf.name]
    else:
        # Load from HyPhy test data NEXUS file
        from predict_regression_nexus import parse_nexus_alignment_and_embedded_tree
        raw_path = f"/Users/sergei/Development/hyphy/tests/data/{sc['tree_name']}"
        _, original_taxlabels, original_tree_str = parse_nexus_alignment_and_embedded_tree(raw_path)
        
        # Process Newick tree
        clean_tree_str = re.sub(r'\{[^}]*\}', '', original_tree_str)
        clean_tree_str = re.sub(r'\[.*?\]', '', clean_tree_str).strip().rstrip(';') + ';'
        tree = Phylo.read(StringIO(clean_tree_str), 'newick')
        
        # Estimate branch lengths from original alignment if missing or all zero
        has_branch_lengths = any(
            clade.branch_length is not None and clade.branch_length > 0.0 
            for clade in tree.find_clades() 
            if clade != tree.root
        )
        if not has_branch_lengths:
            est_tree = estimate_branch_lengths(raw_path, tree)
            if est_tree:
                tree = est_tree
    
    # Fill in any missing branch lengths to prevent pyvolve validation error
    for node in tree.find_clades():
        if node.branch_length is None:
            node.branch_length = 0.001
            
    # Sanitize tree branch names to prevent pyvolve parsing crash
    for node in tree.find_clades():
        if node.name:
            node.name = sanitize_name(node.name)
            
    # Compile taxlabels and sanitize
    taxlabels = [sanitize_name(name) for name in original_taxlabels]
    
    # 2. Select branches for episodic selection
    target_branches = []
    if sc['regime'] == 'Episodic' or sc['regime'] == 'Hybrid-MutSel':
        target_branches = select_branches_for_flagging(tree, sc['branch_setting'], sc['tree_name'])
        
    # Build flagged Newick tree string
    # We flag target branches with #m1
    flagged_tree_str = to_newick(tree.clade, target_branches, "#m1") + ";"
    my_pyvolve_tree = pyvolve.read_tree(tree=flagged_tree_str)
    
    # 3. Setup Models & Partitions
    partitions_list = []
    
    if sc['regime'] == 'Null' or sc['regime'] == 'Constant':
        model_neg = pyvolve.Model("codon", {"omega": sc['omega_neg'], "kappa": 2.0})
        part_neg = pyvolve.Partition(models=model_neg, size=sc['num_neg_sites'])
        partitions_list.append(part_neg)
        if sc['num_pos_sites'] > 0:
            model_pos = pyvolve.Model("codon", {"omega": sc['omega_pos'], "kappa": 2.0})
            part_pos = pyvolve.Partition(models=model_pos, size=sc['num_pos_sites'])
            partitions_list.append(part_pos)
            
    elif sc['regime'] == 'Episodic':
        # Branch heterogeneity: m1 is the foreground, rootmodel is background
        m1 = pyvolve.Model("codon", {"omega": sc['omega_pos'], "kappa": 2.0}, name="m1")
        rootmodel = pyvolve.Model("codon", {"omega": sc['omega_neg'], "kappa": 2.0}, name="rootmodel")
        
        # Partition 1: Purifying/neutral background throughout
        part_neg = pyvolve.Partition(models=rootmodel, size=sc['num_neg_sites'])
        partitions_list.append(part_neg)
        
        # Partition 2: Evolving under episodic selection (m1 on foreground, rootmodel elsewhere)
        part_pos = pyvolve.Partition(models=[m1, rootmodel], size=sc['num_pos_sites'], root_model_name="rootmodel")
        partitions_list.append(part_pos)
        
    elif sc['regime'] == 'MutSel':
        # Pure MutSel domain (e.g. Hydrophobic)
        prop = sc['model_type'].split('-')[1]
        fit = get_biophysical_fitness(prop, select_coef=4.0)
        model_mutsel = pyvolve.Model("mutsel", {"fitness": fit})
        part_neg = pyvolve.Partition(models=model_mutsel, size=sc['num_neg_sites'])
        partitions_list.append(part_neg)
        
    elif sc['regime'] == 'Hybrid-MutSel':
        # Hybrid MutSel background + episodic selection
        prop = sc['model_type'].split('-')[1].split('+')[0]
        fit = get_biophysical_fitness(prop, select_coef=4.0)
        model_mutsel = pyvolve.Model("mutsel", {"fitness": fit}, name="rootmodel")
        
        # Episodic positive selection model m1
        m1 = pyvolve.Model("codon", {"omega": sc['omega_pos'], "kappa": 2.0}, name="m1")
        
        # Partition 1: MutSel background (conserved structural domain)
        part_neg = pyvolve.Partition(models=model_mutsel, size=sc['num_neg_sites'])
        partitions_list.append(part_neg)
        
        # Partition 2: Evolving under episodic selection (m1 on foreground, MutSel on background)
        part_pos = pyvolve.Partition(models=[m1, model_mutsel], size=sc['num_pos_sites'], root_model_name="rootmodel")
        partitions_list.append(part_pos)
        
    # 4. Simulate Evolution
    fasta_path = os.path.join(scratch_dir, f"sim_temp_{sc['scenario_id']}_{rep_idx}.fasta")
    evolver = pyvolve.Evolver(tree=my_pyvolve_tree, partitions=partitions_list)
    evolver(seqfile=fasta_path, write_joint_states=False, infofile=None)
    
    # 5. Build NEXUS with embedded tree (unflagged tree to prevent model parsing issues)
    unflagged_tree_str = to_newick(tree.clade, [], "") + ";"
    
    seq_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[record.id] = str(record.seq).upper()
        
    # Save gzipped alignment in the replicates folder
    nexus_gz_path = os.path.join(replicates_dir, f"scen_{sc['scenario_id']}_rep_{rep_idx}.nex.gz")
    write_nexus(seq_dict, taxlabels, unflagged_tree_str, nexus_gz_path)
    
    # Cleanup temp fasta
    if os.path.exists(fasta_path):
        os.remove(fasta_path)
        
    return nexus_gz_path, taxlabels

def analyze_prediction(pred_csv_path, sc):
    """Analyzes prediction output file and returns result dictionary."""
    df = pd.read_csv(pred_csv_path)
    
    df_neg = df.iloc[:sc['num_neg_sites']]
    df_pos = df.iloc[sc['num_neg_sites']:]
    
    mean_neg_lrt = df_neg["predicted_lrt"].mean()
    mean_pos_lrt = df_pos["predicted_lrt"].mean() if len(df_pos) > 0 else 0.0
    
    tp = len(df_pos[df_pos["selection_call"] != "Neutral"]) if len(df_pos) > 0 else 0
    fn = len(df_pos[df_pos["selection_call"] == "Neutral"]) if len(df_pos) > 0 else 0
    
    fp = len(df_neg[df_neg["selection_call"] != "Neutral"])
    tn = len(df_neg[df_neg["selection_call"] == "Neutral"])
    
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    fdr = fp / (fp + tp) if (fp + tp) > 0 else 0.0
    
    return {
        "tpr": tpr,
        "fpr": fpr,
        "fdr": fdr,
        "mean_neg_lrt": mean_neg_lrt,
        "mean_pos_lrt": mean_pos_lrt
    }

def main():
    parser = argparse.ArgumentParser(description="AxoMeme Comprehensive Simulation Benchmark Runner.")
    parser.add_argument("--model", type=str, default="/Users/sergei/Projects/TOGA_MEME/MEME_transformer_joint.pt",
                        help="Path to trained model weights.")
    parser.add_argument("--num_replicates", type=int, default=100, help="Number of replicates per scenario.")
    parser.add_argument("--tier1_lrt_gate", type=float, default=5.0, help="LRT gate for Tier 1 calls.")
    parser.add_argument("--tier2_lrt_gate", type=float, default=3.0, help="LRT gate for Tier 2 calls.")
    parser.add_argument("--scratch_dir", type=str, default="/Users/sergei/Projects/TOGA_MEME/scratch",
                        help="Folder to store temporary files.")
    
    args = parser.parse_args()
    
    db_path = os.path.join(args.scratch_dir, "sim_benchmark_results.db")
    report_path = "/Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/simulation_dashboard.html"
    os.makedirs(args.scratch_dir, exist_ok=True)
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    
    replicates_dir = os.path.join(args.scratch_dir, "simulated_replicates")
    os.makedirs(replicates_dir, exist_ok=True)
    
    # 1. Setup DB
    print("[*] Initializing SQLite benchmark database...")
    setup_db(db_path)
    
    # 2. Setup Scenarios
    scenarios = setup_scenarios()
    total_scen = len(scenarios)
    print(f"[*] Configured {total_scen} simulation scenarios across 4 phylogenetic backdrops.")
    
    # Check current progress
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("SELECT COUNT(*) FROM results")
    completed_reps = c.fetchone()[0]
    conn.close()
    print(f"[*] Starting benchmark. Completed reps in DB: {completed_reps} / {total_scen * args.num_replicates}")
    
    # 3. Main benchmark loop
    current_run_idx = completed_reps
    
    # Append sys path to import parser from predict_regression_nexus
    sys.path.append("scripts")
    
    for sc in scenarios:
        scen_id = sc['scenario_id']
        print(f"[*] Processing Scenario {scen_id}/{total_scen} ({sc['tree_name']}, {sc['regime']})...")
        
        # Determine missing replicates
        missing_reps = []
        for rep_idx in range(1, args.num_replicates + 1):
            conn = sqlite3.connect(db_path)
            c = conn.cursor()
            c.execute("SELECT 1 FROM results WHERE scenario_id = ? AND replicate = ?", (scen_id, rep_idx))
            exists = c.fetchone()
            conn.close()
            if not exists:
                missing_reps.append(rep_idx)
                
        if not missing_reps:
            print(f"[*] Scenario {scen_id} already fully completed.")
            continue
            
        print(f"[*] Simulating {len(missing_reps)} missing replicates...")
        align_list = []
        output_list = []
        reps_to_predict = []
        max_sp = 256
        
        for rep_idx in missing_reps:
            try:
                # Simulate the replicate (step 1-5)
                nexus_gz_path, taxlabels = simulate_only(sc, rep_idx, args.scratch_dir, replicates_dir)
                max_sp = min(256, len(taxlabels))
                pred_csv_path = os.path.join(args.scratch_dir, f"pred_{scen_id}_{rep_idx}.csv")
                
                align_list.append(nexus_gz_path)
                output_list.append(pred_csv_path)
                reps_to_predict.append((rep_idx, pred_csv_path))
            except Exception as e:
                print(f"[!] Warning: Simulation failed for replicate {rep_idx} on Scenario {scen_id}: {e}")
                
        if not align_list:
            print(f"[!] No replicates simulated successfully for Scenario {scen_id}.")
            continue
            
        # Run batched prediction
        print(f"[*] Running batch predictions for {len(align_list)} replicates...")
        cmd = [
            "python3", "scripts/predict_regression_nexus.py",
            "--alignment"
        ] + align_list + [
            "--model", args.model,
            "--output"
        ] + output_list + [
            "--tier1_lrt_gate", str(args.tier1_lrt_gate),
            "--tier2_lrt_gate", str(args.tier2_lrt_gate),
            "--max_species", str(max_sp)
        ]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Analyze predictions and save to DB
            scen_reps_completed = args.num_replicates - len(missing_reps)
            for rep_idx, pred_csv_path in reps_to_predict:
                if not os.path.exists(pred_csv_path):
                    print(f"[!] Warning: Prediction output not found for replicate {rep_idx} on Scenario {scen_id}")
                    continue
                try:
                    metrics = analyze_prediction(pred_csv_path, sc)
                    res = {
                        "scenario_id": sc['scenario_id'],
                        "replicate": rep_idx,
                        "tree_name": sc['tree_name'],
                        "regime": sc['regime'],
                        "model_type": sc['model_type'],
                        "omega_neg": sc['omega_neg'],
                        "omega_pos": sc['omega_pos'],
                        "num_neg_sites": sc['num_neg_sites'],
                        "num_pos_sites": sc['num_pos_sites'],
                        "branch_setting": sc['branch_setting'],
                        **metrics
                    }
                    save_result_to_db(db_path, res)
                    current_run_idx += 1
                    scen_reps_completed += 1
                except Exception as e:
                    print(f"[!] Warning: Failed parsing replicate {rep_idx} on Scenario {scen_id}: {e}")
                finally:
                    if os.path.exists(pred_csv_path):
                        os.remove(pred_csv_path)
            
            generate_dashboard(db_path, report_path, total_scen, current_run_idx)
            print(f"[*] Completed Scenario {scen_id}: {scen_reps_completed}/{args.num_replicates} replicates done. Dashboard updated. (Total progress: {current_run_idx} / {total_scen * args.num_replicates})")
        except Exception as e:
            print(f"[!] Warning: Batched prediction call failed for Scenario {scen_id}: {e}")
            
    # Final dashboard update
    generate_dashboard(db_path, report_path, total_scen, current_run_idx)
    print(f"🎉 Benchmark complete! Final results compiled and saved to '{report_path}'")

if __name__ == "__main__":
    main()
