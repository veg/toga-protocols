# Standard Operating Procedure (SOP): MEME SQLite Database Management

This standard operating procedure describes how to maintain, curate, and optimize the central SQLite database (`meme_results.db`) containing evolutionary selection stats.

---

## 1. Schema Design and Key Tables
The SQLite database integrates site-level selection scores, gene metadata, and substitution maps. It consists of three primary tables:
* **`site_results`**: Stores Likelihood Ratio Test (LRT) statistics, p-values, FDR-corrected q-values, and quality control classifications for each codon site.
* **`gene_results`**: Contains gene-level attributes and counts of positively selected sites.
* **`site_substitutions`**: Records parsed ancestral reconstruction events and amino acid change coordinates.

---

## 2. Database Update Pipeline
To update selection stats when new alignments or cluster runs complete, follow this pipeline:

### Step 1: Synchronize Results
Run the synchronization script (`pull_and_update_all.py`) to pull raw gzipped MEME json output files from processing nodes (`m2.local`, `m3.local`, `magilla`, etc.) and load them incrementally into the database tables.

### Step 2: Recalculate FDR Correction
Raw p-values must be corrected for multiple testing using the Benjamini-Hochberg False Discovery Rate (FDR) procedure. Run the FDR calculator (`populate_meme_qvals.py`) to update the `q_value` column. Codon sites with $q \le 0.10$ are marked as statistically significant (`is_significant = 1`).

### Step 3: Run False Positive Curation Filters
To weed out artifactual selection signals, run the curation filters script (`update_meme_db_filter.py`). It applies the following quality checks:
1. **Spatial Clustering Filter**: Flags contiguous runs of significant sites (e.g. $\ge 5$ sites in a close window) as `LIKELY_ERROR` (often caused by local sequence misalignment).
2. **Flanking Gap Density**: Flags sites residing in gap-heavy regions ($>30\%$ gap or masking characters within a 21-codon window).
3. **Multi-Nucleotide Substitutions (MNS)**: Reviews amino acid transitions requiring multiple base changes to ensure alignment support.

---

## 3. Database Indexes and Performance Optimization
To prevent queries from scanning tables linearly during training or dashboard queries, maintain the following composite indexes:
* `idx_site_gene_idx` on `site_results(gene_name, site_index)`
* `idx_subst_branch_gene` on `site_substitutions(branch_name, gene_name)`
