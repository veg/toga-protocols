# TOGA Selection Analyses: Overview

This repository contains protocols, automated curation pipelines, and machine-learning models developed to screen and validate positive selection across mammalian genomes aligned via **TOGA (Tool for Genome Alignment)**.

---

## 1. Analysis Architecture

Our evolutionary selection pipeline consists of three core computational components:

### A. Codon-Level Selection Screening (HyPhy MEME)
* **Goal**: Identify individual codon sites subject to episodic diversifying selection.
* **Method**: Run **HyPhy MEME (Mixed Effects Model of Evolution)** across TOGA-reconciled alignments.
* **Significance Threshold**: Default selection calls are made at **FDR $q \le 0.10$** (computed via Benjamini-Hochberg correction on asymptotic LRT $p$-values) to control false discoveries across thousands of alignment scans.

### B. Branch-Level Selection Screening (HyPhy aBSREL)
* **Goal**: Detect positive selection occurring along specific lineages (internal nodes or terminal branches) rather than individual sites.
* **Method**: Run **HyPhy aBSREL (adaptive Branch-Site Random Effects Likelihood)** to test which branches in the phylogeny show evidence of $\omega > 1$.

### C. Deep Learning Validation (PhyloAxialTransformer)
* **Goal**: Validate and regularize MLE selection statistics using a trained regression transformer.
* **Method**: The model processes sliding window alignment tokens, patristic distances, and multidimensional scaling (MDS) tree coordinates to output regularized predicted Likelihood Ratio Test (LRT) statistics.
* **Aesthetic & Alignment**: It acts as a regularized predictor, smoothing out MLE asymptotic noise (such as sequencing errors, high-frequency gap regions, or branch-length noise) to predict robust selective landscapes.

---

## 2. Quality Control & Curation Filters

To filter out false-positive selection calls caused by alignment artifacts or sequencing errors, all significant sites ($q \le 0.10$) are classified into three quality tiers:

1. **GOLD (High-Quality Single Substitution)**: Verified single nucleotide substitution with clean flanking homology.
2. **SILVER (High-Quality Multi-Nucleotide Substitution)**: Multi-nucleotide change (Hamming distance $\ge$ 2) occurring in a region of high local flanking sequence homology (no gaps or Ns within 5 flanking codons).
3. **LIKELY_ERROR (Filtered Out)**: Sites that fail due to:
   * **Low descendant representation**: $< 50\%$ valid sequence data in carrying lineages.
   * **High gap density**: $> 30\%$ gap/N density in the flanking window.
   * **Isolated short fragments**: Substitution carried on a sequence fragment shorter than 30 amino acids.
   * **Spatial run contamination**: Contiguous runs of $\ge 5$ significant sites within close proximity (spacing $\le 3$ codons), representing alignment shifts or sequencing anomalies.

---

## 3. Directory Structure

* `protocols/`:
  * `human/`: Standard Operating Procedures (SOPs) for human scientists.
    * [genome_selection_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/genome_selection_SOP.md): Selecting the best assembly representation per species.
    * [tree_trimming_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/tree_trimming_SOP.md): Pruning the species tree and handling name collisions.
    * [alignment_curation_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/alignment_curation_SOP.md): Filtering MSAs and sequence records.
    * [meme_analysis_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/meme_analysis_SOP.md): Running selection screens.
    * [model_training_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/model_training_SOP.md): Training and checkpointing the AxoMeme model.
    * [database_management_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/database_management_SOP.md): Curation and optimization workflows for the SQLite database.
    * [dashboard_deployment_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/dashboard_deployment_SOP.md): Deploying, caching, and exporting the local and static dashboards.
    * [legacy_comparison_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/legacy_comparison_SOP.md): Comparing TOGA to historical alignments and classifying discrepancies.
    * [recent_genomes_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/recent_genomes_SOP.md): Evaluating the contribution of newly sequenced mammalian genomes to selection statistics.
    * [evolutionary_hypotheses_SOP.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/human/evolutionary_hypotheses_SOP.md): Standard procedures for testing and controlling evolutionary hypotheses.
  * `llm/`: Detailed, structured, and contextual protocols optimized for LLM agents.
    * [genome_selection_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/genome_selection_LLM.md)
    * [tree_trimming_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/tree_trimming_LLM.md)
    * [alignment_curation_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/alignment_curation_LLM.md)
    * [meme_analysis_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/meme_analysis_LLM.md)
    * [model_training_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/model_training_LLM.md)
    * [database_management_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/database_management_LLM.md)
    * [dashboard_deployment_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/dashboard_deployment_LLM.md)
    * [legacy_comparison_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/legacy_comparison_LLM.md)
    * [recent_genomes_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/recent_genomes_LLM.md)
    * [evolutionary_hypotheses_LLM.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/protocols/llm/evolutionary_hypotheses_LLM.md)
* `reports/`: Detailed evolutionary hypothesis testing reports, functional enrichment analyses, and diagnostic visualizations:
  * [hypothesis_1_detailed_report.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/hypothesis_1_detailed_report.md): Selection efficiency ($N_e$) vs branch length noise.
  * [hypothesis_2_detailed_report.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/hypothesis_2_detailed_report.md): Clade-specific selection shifts.
  * [hypothesis_3_detailed_report.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/hypothesis_3_detailed_report.md): Tree length and taxon density.
  * [ultraconserved_enrichment_report.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/ultraconserved_enrichment_report.md): Functional enrichment of ultraconserved genes.
  * [homogeneous_purifying_selection_enrichment_report.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/homogeneous_purifying_selection_enrichment_report.md): Functional enrichment of genes under global purifying selection.
  * [grant_proposal_innovation_section.md](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/reports/grant_proposal_innovation_section.md): Draft text block for grant proposals focusing on epistasis and phenotypic adaptation mapping.
  * `images/`: Diagnostic validation plots and heatmaps.
* `data/`:
  * `species_to_assembly.csv`: Mapping of species to their chosen representative genome assembly.
  * `assemblies_and_species.csv`: Original metadata for all 1,971 assemblies examined.
  * `problematic_genomes.csv`: Curation report of flagged problematic assemblies.
  * `pruned_species_tree.nhx`: Cleaned species tree topology pruned to 742 representative species.
  * `species_to_taxon_map.csv`: Mapping of assemblies to their 6-character taxonomic codes.
* `scripts/`: Production-ready automation scripts for compiling, dashboard building, training, and model testing:
  * [generate_dashboards.py](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/scripts/generate_dashboards.py): Curation and species tree dashboard builder.
  * [run_training_pipeline.sh](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/scripts/run_training_pipeline.sh): Automation script for training and validating models.
  * [run_axomeme_simulation.py](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/scripts/run_axomeme_simulation.py): Evolutionary simulation and performance testing pipeline using `pyvolve`.
  * [run_joint_comparison.py](file:///Users/sergei/Projects/TOGA_MEME/toga_protocols/scripts/run_joint_comparison.py): Joint comparison pipeline running both HyPhy MEME and AxoMeme side-by-side on simulated data.
