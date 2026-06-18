# Standard Operating Procedure (SOP): Legacy vs. TOGA Alignment Comparison

This standard operating procedure describes how to perform and interpret comparative selection analyses between **legacy** mammalian alignments (~45 species) and modern, dense **VGP/TOGA-level** mammalian alignments (~600–700 species).

---

## 1. Overview of Comparative Sensitivity
Dense taxon representation dramatically increases statistical power to resolve episodic selection pressure. Historical comparative results show:
* **Sensitivity Drop**: Legacy datasets recover only **~20%** of the selection hotspots detected at $p \le 0.01$ in dense TOGA datasets.
* **Signal Erasure**: More than **55%** of strongly selected TOGA sites show no signal whatsoever ($p > 0.05$) in legacy datasets, indicating that historical studies misclassify these positions as entirely neutral.
* **Support Loss**: Likelihood Ratio Test (LRT) statistics drop by an average of **~3.9-fold** in legacy alignments, reflecting a lack of evolutionary replication across branches.

---

## 2. Classification of Discrepancies (Legacy-Only Detections)
When a site is highly significant in a legacy alignment ($p \le 0.01$) but is neutral in the full TOGA alignment ($p > 0.05$), the discrepancy must be classified into one of two categories:

### Category A: Small-Tree Ancestral Reconstruction Errors (LRT Spikes)
In small alignments, the model lacks closely-related sister lineages. During ancestral state reconstruction, this under-sampling leads to reconstruction errors over wide branches. 
* **Mechanism**: A simple synonymous mutation is misidentified as a non-synonymous change due to a misreconstructed ancestral node. On a small tree, this artifact creates a massive, false selection score.
* **Resolution**: Populating the tree with sister taxa (e.g., in TOGA) resolves the correct ancestral state, correcting the false positive (e.g. AOX1 Site 697/744).

### Category B: Neutral Null Model Stabilization and Localized Signal Dilution
On small trees, background variation is poorly estimated, causing minor fluctuations in non-synonymous rates on internal branches to appear as positive selection. 
* **Mechanism**: Dense alignments stabilize the neutral null model by adding dozens of synonymous mutations across other lineages. This establishes the true background mutation rate and reveals that the site's overall mutation pattern is consistent with neutral or weak purifying selection.
* **Dilution**: Additionally, if selection is real but localized to a specific branch/clade, scaling to a large tree "dilutes" the global score unless branch-specific tests (like aBSREL) are run (e.g. AOX1 Site 667/714).

---

## 3. Scientific Guidelines
When reviewing legacy comparisons:
1. **Prioritize dense alignments**: Do not treat legacy-only significant sites as biological truth without checking if they survive in the larger TOGA dataset.
2. **Review ancestral state reconstructions**: Inspect the phylogenetic placement of substitutions to confirm they are not reconstruction errors.
