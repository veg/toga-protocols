# Standard Operating Procedure (SOP): Selection Screen and Model Predictions

This protocol details how to sync MEME results, update the central database, launch dashboards, and run transformer predictions.

---

## 1. Updating the Central Database
To download new MEME results from the clusters (`magilla`, `silverback`, `m3.local`), compile them, and build the local SQLite database:

1. Open a terminal and run the main synchronization script:
   ```bash
   python3 scripts/pull_and_update_all.py
   ```
2. This script runs a complete pipeline:
   * Syncs raw `.json.gz` MEME files from remote hosts.
   * Compiles them incrementally into `meme_results.db` on `m3.local`.
   * Runs curation filters to mark sites as `GOLD`, `SILVER`, or `LIKELY_ERROR`.
   * Downloads the compiled database locally to the repository root.
   * Enforces the default significance threshold ($q \le 0.10$).
   * Precomputes branch mappings and selection classifications.

---

## 2. Launching the Local Dashboards
We run three separate FastAPI dashboards to analyze the datasets:

* **MEME Results Dashboard** (Port 8082):
  ```bash
  python3 scripts/serve_meme_dashboard_fastapi.py 8082
  ```
* **aBSREL Results Dashboard** (Port 8083):
  ```bash
  python3 scripts/serve_absrel_dashboard.py 8083
  ```
* **Alignment Curation Dashboard**:
  ```bash
  python3 scripts/serve_alignments_dashboard.py
  ```

Once launched, access the dashboards in your web browser at `http://localhost:[port]`.

---

## 3. Running Transformer Selection Predictions
To predict codon-level selection strengths on new alignments using the trained model:

1. Use the inference script [predict_regression_nexus.py](../../scripts/predict_regression_nexus.py):
   ```bash
   python3 scripts/predict_regression_nexus.py \
       --alignment msa/your_gene.nex \
       --model MEME_transformer_joint.pt \
       --output your_gene_predictions.csv \
       --device cpu
   ```
2. **Device Enforcement**: Always enforce **`--device cpu`** when running on macOS/local environments. PyTorch MPS has an active kernel bug that causes mathematically incorrect eigenvector scaling, yielding false predictions.
3. **Thresholds & Gates**: The script combines local relative percentiles (to scale for baseline divergence differences) with absolute predicted LRT gates (to handle dense selection in highly diversifying genes like camelid/vif):
   * **Tier 1 (High Confidence)**: Top 2% of variable sites (`--tier1_percentile 98.0`) **or** absolute predicted LRT $\ge 5.0$ (`--tier1_lrt_gate 5.0`).
   * **Tier 2 (Medium Confidence)**: Top 3% of variable sites (`--tier2_percentile 97.0`) **or** absolute predicted LRT $\ge 3.0$ (`--tier2_lrt_gate 3.0`).
4. **Estimated Trees**: If the input tree lacks branch lengths, the script automatically estimates HKY85 branch lengths using HyPhy in the background and saves the estimated tree as `[output_prefix]_estimated_tree.nwk`.
