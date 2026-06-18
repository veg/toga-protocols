# LLM Integration & Automation Protocol: TOGA Selection Analyses

This protocol is optimized for LLM coding agents assisting with TOGA evolutionary genomics selection screens. It documents the codebase structure, critical numerical edge cases, database schemas, and integration steps.

---

## 1. Codebase & Data Layout
* **Primary Database**: `meme_results.db` (SQLite, stored in root).
* **Model Checkpoint**: `MEME_transformer_joint.pt` (PyTorch, stored in root or home).
* **Alignments**: NEXUS format (gzipped or plain text) stored under `msa/` or `tests/`.
* **Dashboard Servers**: fastapi/dashboards in `scripts/serve_meme_dashboard_fastapi.py`.

---

## 2. Critical Engineering & Mathematical Caveats

### A. PyTorch MPS Kernel Bug (macOS)
> [!WARNING]
> Running the continuous selection PhyloAxialTransformer model on `mps` (Metal Performance Shaders) yields mathematically incorrect predictions due to a silent kernel bug in PyTorch's MPS backend.
* **Instruction**: Always force CPU execution (`--device cpu`) on macOS. Do not auto-detect `mps` even if available.

### B. Eigenvector Sign Convention (MDS Coordinates)
* **Context**: Tree topologies are converted into Patristic distance matrices, which are projected into 4D coordinate space using Multidimensional Scaling (MDS).
* **Numerical Problem**: Eigenvector computations via `numpy.linalg.eigh` are unique only up to sign. Noise or library shifts can cause signs of coordinate columns to flip arbitrarily, causing prediction shifts in the transformer projections.
* **Resolution**: Enforce sign convention normalization:
  ```python
  for col in range(evecs.shape[1]):
      max_abs_idx = np.argmax(np.abs(evecs[:, col]))
      sign = np.sign(evecs[max_abs_idx, col])
      if sign < 0:
          evecs[:, col] *= -1.0
  ```

### C. Flat Tree Branch Length Fallback
* **Context**: NEXUS alignments sometimes embed trees with missing or zero branch lengths, preventing patristic distance calculations.
* **Resolution**: Write a temporary FASTA and a HyPhy batch file (`.bf`) to run HKY85 maximum-likelihood branch-length estimation via HyPhy, outputting estimated branch lengths to prune and calculate patristic distances.

### D. Joint Calling Strategy (Percentiles + Gates)
* **Goal**: Maximize precision and F1-score across both sparse-selection genes (like SLCs) and dense-selection genes (like camelid).
* **Logic**:
  * **Tier 1 (High)**: `local_percentile >= 98.0` OR `predicted_lrt >= 5.0`
  * **Tier 2 (Medium)**: (`local_percentile >= 97.0` OR `predicted_lrt >= 3.0`) AND NOT Tier 1

---

## 3. SQLite Database Schema & Curation Logic

### `site_results` Table Schema
* `gene_name` (TEXT)
* `site_index` (INTEGER, 1-based)
* `is_significant` (INTEGER, 0 or 1, threshold: Benjamini-Hochberg $q \le 0.10$)
* `q_value` (REAL)
* `classification` (TEXT: `'GOLD'`, `'SILVER'`, `'LIKELY_ERROR'`, or `'NONE'`)
* `filter_reason` (TEXT)

### Curation Pipeline Rules
1. **Spatial Run Override**: Identify consecutive significant sites with gaps $\le 3$ codons. If a run contains $\ge 5$ significant sites, override their classification to `LIKELY_ERROR` (preserving individually validated high-quality `GOLD`/`SILVER` sites).
2. **Multi-Nucleotide Substitution (MNS)**: Codon mutations with Hamming distance $\ge 2$. Silver status is granted only if the carrying sequence has zero gaps or Ns in the 5 flanking codons (15 nucleotides) upstream and downstream.
3. **Flanking Gap Density**: Exclude sites wherecarrying sequence gap/N density exceeds $30\%$ in a 21-codon window.
