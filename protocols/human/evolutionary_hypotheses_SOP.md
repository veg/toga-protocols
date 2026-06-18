# Standard Operating Procedure (SOP): Testing Evolutionary Hypotheses

This standard operating procedure describes how to perform, control, and interpret tests for the three primary evolutionary hypotheses using the comparative database.

---

## Hypothesis 1: Effective Population Size ($N_e$) and Selection Efficiency
This hypothesis tests if clades with large effective population sizes (such as Rodentia and Chiroptera) exhibit more efficient selection (higher purifying constraint and positive selection rates) than clades with small effective population sizes (Primates and Cetacea).

### Core Controls and Interpretation:
* **Branch Length Controls**: Denser taxonomic sampling in specific clades (such as primates and cetaceans) leads to systematically shorter branch lengths in the database. Shorter branches reduce the time for mutations to occur, causing mathematical artifacts.
* **Purifying Selection**: Background constraint ($\omega$) must be adjusted for branch lengths when comparing clades. In terminal branches, large $N_e$ clades show systematically lower median $\omega$ values (stronger constraint) compared to small $N_e$ clades, matching population genetics expectations.
* **Positive Selection (aBSREL)**: Clades with small $N_e$ often exhibit higher positive selection detection rates because they are heavily sampled, generating shorter branch length artifacts that inflate statistical significance. Control for branch length is mandatory to confirm true biological selection.

---

## Hypothesis 2: Clade-Specific Adaptive Shifts vs. Ancestral Constancy
This hypothesis tests whether positive selection is localized to specific mammalian clades (Euarchontoglires vs. Laurasiatheria) or represents ancestral constancy.

### Control Null Model:
* **Stratified Permutation**: Shuffle selected branches within leaf and internal branch pools across the tree topology.
* **Interpretation**: A high proportion of clade-specific selection sites (occurring exclusively in Euarchontoglires or Laurasiatheria) relative to the null model indicates true lineage-specific adaptive shifts rather than uniform selection.

---

## Hypothesis 3: Taxonomic Sampling Density and Saturation of Selection
This hypothesis tests whether adding species to alignments yields diminishing returns (saturation) in positive selection site discovery.

### Biphasic Curation Findings:
1. **Housekeeping Constraint**: Densely sampled genes (e.g. 700+ species) have systematically lower selection densities. This is because genes conserved across all 700 genomes represent core housekeeping genes under intense purifying selection, not selection saturation.
2. **Mutational Depth (Tree Length)**: Selection detection density is strongly correlated with total Tree Length (substitutions/site) rather than species count. Mutational depth, not sequence count, drives selection detection power. Logarithmic models of tree length provide the best predictive fit.

---

## Curation of Extreme Constraint: Homogeneous Purifying and Ultraconservation
To identify genes under pure negative selection, filter the database for genes with zero selected sites ($q \le 0.10$) across high mutational depth (Tree Length $\ge 10.0$):
* **Homogeneous Purifying Selection ($0.10 \le dN/dS \le 0.50$)**: Enriched for translation, gene expression, and mitochondrial pathways.
* **Ultraconservation ($dN/dS < 0.10$)**: Enriched for developmental cell junction assembly, cell-cell adhesion, and action potentials.
