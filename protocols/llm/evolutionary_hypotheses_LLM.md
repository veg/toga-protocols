# LLM Protocol: Algorithmic Testing of Evolutionary Hypotheses

This protocol specifies the mathematical controls, data thresholds, and database queries for testing the three primary evolutionary hypotheses.

---

## 1. Hypothesis 1: $N_e$ and Selection Efficiency (Branch Length Controls)
To test the impact of effective population size $N_e$ on purifying and positive selection:

1. **Clade Classification**:
   * **Large $N_e$ Clades**: Rodentia, Chiroptera.
   * **Small $N_e$ Clades**: Primates, Cetacea.
2. **Branch Length Diagnostics**:
   * Calculate median terminal branch lengths for each clade. Systematic differences in sampling density inflate selection signals in short-branched clades (primarily Primates and Cetacea).
3. **Purifying Constraint Control**:
   * Let $L_{branch}$ be the branch length.
   * Calculate background $\omega$ for each terminal branch. Contrast the median background $\omega$ between Large $N_e$ and Small $N_e$ clades using a Mann-Whitney U test. Enforce $L_{branch}$ stratification to compare branches of equivalent length.
4. **aBSREL Positive Selection Check**:
   * Calculate detection rate:
     $$\text{Rate} = \frac{N_{sig}}{N_{total}}$$
   * Run a Chi-squared test of independence on the significant vs. non-significant branch counts between population classes.

---

## 2. Hypothesis 2: Clade-Specific Shifts vs. Ancestral Constancy
To isolate clade-specific adaptive shifts:

1. **Site Classification Tiers**:
   * **Euarchontoglires-only**: Significant positive selection (MEME EBF > 20) on one or more branches in Euarchontoglires, and zero branches in Laurasiatheria.
   * **Laurasiatheria-only**: Significant positive selection on one or more branches in Laurasiatheria, and zero branches in Euarchontoglires.
   * **Shared/Parallel**: Positive selection detected in at least one branch of both clades.
2. **Stratified Null Permutation Test**:
   * Shuffling selection scores: For each gene, extract the set of branches under selection.
   * Stratify branches into "leaf" and "internal" pools.
   * Randomly permute the selection flags within each pool across the tree topology $N = 10,000$ times.
   * Contrast the observed clade-specific counts against the permuted null distributions to compute empirical p-values.

---

## 3. Hypothesis 3: Taxonomic Sampling Density and Power Saturation
To evaluate how sequence count and tree length impact positive selection detection:

1. **Housekeeping Gene Control**:
   * Densely sampled genes ($\ge 600$ species) are biologically distinct from sparse alignments. Filter genes by functional class before fitting regression models.
2. **Tree Length vs. Selection Density Correlation**:
   * Let $D_{selection}$ be the selection density (fraction of codon sites with $q \le 0.10$):
     $$D_{selection} = \frac{N_{sig\_sites}}{N_{total\_sites}}$$
   * Compute the Spearman rank correlation $r_s$ between $D_{selection}$ and overall Tree Length (sum of branch lengths, substitutions/site).
   * Fit a logarithmic curve:
     $$D_{selection} \sim \beta_0 + \beta_1 \ln(\text{Tree Length})$$

---

## 4. Pure Negative Selection Curation Filters
To extract genes evolving under strict purifying selection (for functional enrichment):

### A. Homogeneous Purifying Selection Genes
* **Taxonomic Depth**: Species count $\ge 600$.
* **Mutational Depth**: Tree length $\ge 10.0$ substitutions/site.
* **Positive Selection Constraint**: Number of significant sites ($q \le 0.10$) = 0.
* **Purifying Range**: Average background $0.10 \le dN/dS \le 0.50$.

### B. Ultraconserved Mammalian Genes
* **Taxonomic Depth**: Species count $\ge 600$.
* **Mutational Depth**: Tree length $\ge 10.0$ substitutions/site.
* **Positive Selection Constraint**: Number of significant sites ($q \le 0.10$) = 0.
* **Extreme Purifying Range**: Average background $dN/dS < 0.10$.
