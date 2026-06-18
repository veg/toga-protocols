# LLM Protocol: MEME Database Curation, Schemas, & Optimization

This protocol defines the schemas, Benjamini-Hochberg calculations, and curation filtering parameters for the central selection database.

---

## 1. Schema Definitions

### A. `site_results` Table
```sql
CREATE TABLE IF NOT EXISTS site_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_name TEXT NOT NULL,
    site_index INTEGER NOT NULL,
    p_value REAL NOT NULL,
    lrt REAL NOT NULL,
    q_value REAL,
    is_significant INTEGER DEFAULT 0,
    classification TEXT DEFAULT 'UNCLASSIFIED',
    filter_reason TEXT,
    selection_type TEXT,
    UNIQUE(gene_name, site_index)
);
```

### B. `gene_results` Table
```sql
CREATE TABLE IF NOT EXISTS gene_results (
    gene_name TEXT PRIMARY KEY,
    num_selected_sites INTEGER DEFAULT 0
);
```

---

## 2. Indexes and Performance
High-speed lookup is critical. Ensure these indexes are applied:
```sql
CREATE INDEX IF NOT EXISTS idx_site_gene_idx ON site_results(gene_name, site_index);
CREATE INDEX IF NOT EXISTS idx_subst_branch_gene ON site_substitutions(branch_name, gene_name);
```

---

## 3. Benjamini-Hochberg FDR Algorithm
To compute or update False Discovery Rate (FDR) q-values:
1. Collect all variable-codon sites $M$ across all genes.
2. Sort site records in ascending order of their raw p-values: $P_{(1)} \le P_{(2)} \le \dots \le P_{(M)}$.
3. For each rank $i$, calculate:
   $$q_i = P_{(i)} \cdot \frac{M}{i}$$
4. Enforce monotonicity by walking backwards from $M$ to $1$:
   $$Q_i = \min_{j \ge i} q_j$$
5. Set `is_significant = 1` for any site where $Q_i \le 0.10$.

Implement in [populate_meme_qvals.py](file:///Users/sergei/Projects/TOGA_MEME/scratch/populate_meme_qvals.py).

---

## 4. Curation & Spatial Run Filtering Math
To weed out false positive selection hotspots caused by indels, misalignment, or frameshifts:

1. **Spatial Run Filter**:
   * For a given gene, retrieve all significant sites ($q \le 0.10$).
   * Run a sliding window of width $W = 10$ codons.
   * If the number of significant sites in the window $N_{sig} \ge 5$ and the gap spacing between consecutive significant sites is $\le 3$ codons:
     * Flag all significant sites in this window as `LIKELY_ERROR` with `filter_reason = "spatial_clustering_cluster_run"`.
2. **Flanking Gap Density**:
   * For each significant site, extract the alignment columns within a flanking window of $\pm 10$ codons (window size = 21 codons).
   * Count the fraction of gap/masking characters (`-`, `?`, `N`) across all species in this window.
   * If the gap fraction $>0.30$, flag the site as `LIKELY_ERROR` with `filter_reason = "high_flanking_gap_density"`.
3. **MNS Filter**:
   * For sites with $q \le 0.10$, inspect codon substitutions.
   * If the substitution requires $\ge 2$ nucleotide changes (Hamming distance $\ge 2$) in a single codon transition:
     * Set `classification = 'SILVER'` unless the surrounding window ($5$ codons flanking) contains zero gap characters across all species.
