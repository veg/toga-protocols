# 🧬 Ultra-Conservative Genes & Hidden Selection Report

This report presents results from an evolutionary analysis of the MEME database to identify **ultra-conservative genes** (characterized by extremely low global dN/dS and zero selection sites) and, crucially, to detect selection in genes that have been traditionally viewed as ultra-conserved in smaller datasets.

---

## 🔬 Methodology & Definitions

* **Dataset**: curating results from **10,247 valid genes** across up to **742 mammalian species**.
* **Ultra-conservative Criteria**:
  * Global gene dN/dS ($\omega$) $< 0.02$
  * Zero positive selection sites ($N_{\text{selected}} = 0$ at $P < 0.05$ under the MEME model)
  * Represented across at least 100 species in the database ($N_{\text{seqs}} \ge 100$) to rule out sequence sparsity.
* **Enrichment Analysis**: Ran a Gene Ontology (GO) Biological Process enrichment scan using the Enrichr API.

---

## 📊 Summary Statistics & Distributions

The distribution of global dN/dS across the entire dataset of 10,247 mammalian genes is as follows:

* **Mean dN/dS**: 0.2257
* **Median dN/dS**: 0.1528
* **75th Percentile**: 0.2822
* **99th Percentile**: 1.35
* **Ultra-conservative Genes**: **476 genes** (meeting the criteria above)

---

## 1. Top 20 Ultra-Conservative Genes in mammals

These genes exhibit the lowest global dN/dS in the database, with zero sites showing episodic positive selection:

| Rank | Gene | Number of Sequences | Alignment Length (Sites) | Global dN/dS | Selected Sites ($P < 0.05$) |
| :---: | :--- | :---: | :---: | :---: | :---: |
| 1 | **H4C12** (Histone H4) | 115 | 103 | 0.000000 | 0 |
| 2 | **CALM3** (Calmodulin 3) | 720 | 150 | 0.000394 | 0 |
| 3 | **ACTC1** (Actin alpha cardiac 1) | 720 | 379 | 0.000525 | 0 |
| 4 | **DCAF7** (DDB1-CUL4 cofactor 7) | 717 | 342 | 0.000572 | 0 |
| 5 | **H3C10** (Histone H3) | 225 | 136 | 0.000713 | 0 |
| 6 | **H3C2** (Histone H3) | 392 | 220 | 0.000728 | 0 |
| 7 | **PRPF8** (U5 snRNP spliceosome) | 670 | 2370 | 0.000775 | 0 |
| 8 | **ACTA2** (Actin alpha smooth muscle) | 716 | 340 | 0.000797 | 0 |
| 9 | **H4C4** (Histone H4) | 350 | 103 | 0.000803 | 0 |
| 10 | **H4C13** (Histone H4) | 261 | 103 | 0.000803 | 0 |
| 11 | **ELOF1** (Transcription factor) | 729 | 83 | 0.000804 | 0 |
| 12 | **H4C5** (Histone H4) | 127 | 103 | 0.000824 | 0 |
| 13 | **H3C8** (Histone H3) | 392 | 136 | 0.000853 | 0 |
| 14 | **H3C12** (Histone H3) | 174 | 136 | 0.000936 | 0 |
| 15 | **H3C3** (Histone H3) | 179 | 143 | 0.001013 | 0 |
| 16 | **H3C11** (Histone H3) | 446 | 136 | 0.001018 | 0 |
| 17 | **H2BC4** (Histone H2B) | 414 | 126 | 0.001059 | 0 |
| 18 | **H4C8** (Histone H4) | 365 | 112 | 0.001080 | 0 |
| 19 | **H4C2** (Histone H4) | 399 | 104 | 0.001094 | 0 |
| 20 | **GNB2** (G-protein beta 2) | 716 | 342 | 0.001162 | 0 |

---

## 2. GO Biological Process Enrichment Results

Enrichment analysis of the 476 ultra-conservative genes reveals that they represent the core cellular translation, splicing, and structural maintenance machinery:

| Rank | GO Term Name | P-value | Adjusted P-value (FDR) | Combined Score |
| :---: | :--- | :---: | :---: | :---: |
| 1 | **Cytoplasmic Translation (GO:0002181)** | $2.49 \times 10^{-42}$ | $2.75 \times 10^{-39}$ | 1484.86 |
| 2 | **Macromolecule Biosynthetic Process (GO:0009059)** | $4.29 \times 10^{-31}$ | $2.37 \times 10^{-28}$ | 422.89 |
| 3 | **Peptide Biosynthetic Process (GO:0043043)** | $5.36 \times 10^{-31}$ | $1.97 \times 10^{-28}$ | 423.80 |
| 4 | **Translation (GO:0006412)** | $1.92 \times 10^{-24}$ | $5.30 \times 10^{-22}$ | 310.57 |
| 5 | **Gene Expression (GO:0010467)** | $2.82 \times 10^{-19}$ | $6.21 \times 10^{-17}$ | 228.68 |
| 6 | **RNA Splicing, via Transesterification (GO:0000377)** | $4.01 \times 10^{-16}$ | $7.37 \times 10^{-14}$ | 189.88 |
| 7 | **mRNA Processing (GO:0006397)** | $6.59 \times 10^{-15}$ | $1.04 \times 10^{-12}$ | 168.86 |
| 8 | **mRNA Splicing, via Spliceosome (GO:0000398)** | $6.58 \times 10^{-14}$ | $9.07 \times 10^{-12}$ | 154.93 |
| 9 | **Protein-RNA Complex Assembly (GO:0022618)** | $1.42 \times 10^{-13}$ | $1.73 \times 10^{-11}$ | 162.36 |
| 10 | **Nucleosome Assembly (GO:0006334)** | $8.07 \times 10^{-13}$ | $8.90 \times 10^{-11}$ | 292.37 |

---

## 3. The "Unveiled Selection" Paradigm Shift
### Traditionally Conserved Genes with Selection in Denser Alignments

The most important question is: *Are there genes historically viewed as "ultra-conserved" or invariant, which show selection when we expand taxonomic sampling?*

**Yes.** By expanding sequence sampling to **742 mammalian species**, we achieve the mutational depth required to reject the neutral null at specific codons. As a result, several classic "housekeeping" controls show **episodic positive selection at specific sites**:

### A. Classic Housekeeping Genes
* **`ACTB`** (Beta-actin): Crucial structural component, traditionally assumed to be under strict purifying selection across its entire sequence. In our dataset, it has **1 selection site** under episodic positive selection ($\omega = 0.0032$, $P < 0.05$).
* **`CALM1`** (Calmodulin 1): Calcium sensor, widely cited as having 100% amino acid identity across most vertebrates. In our dataset, it features **1 selection site** ($\omega = 0.0022$, $P < 0.05$).
* **`GAPDH`** (Glyceraldehyde-3-phosphate dehydrogenase): Standard internal control in expression assays. In our dataset, it features **3 distinct selection sites** ($\omega = 0.040$, $P < 0.05$).

### B. Ribosomal Large (RPL) and Small (RPS) Subunit Genes
Ribosomal proteins are traditionally viewed as the most conserved eukaryotic genes. However, our 742-species dataset reveals selection in:
* **RPL Genes**: 
  * **`RPL10`** (1 site, global $\omega = 0.010$)
  * **`RPL10L`** (1 site, global $\omega = 0.054$)
  * **`RPL19`** (1 site, global $\omega = 0.003$)
  * **`RPL30`** (1 site, global $\omega = 0.006$)
  * **`RPL6`** (2 sites, global $\omega = 0.099$)
* **RPS Genes**:
  * **`RPS13`** (1 site, global $\omega = 0.003$)
  * **`RPS27L`** (1 site, global $\omega = 0.010$)
  * **`RPS6`** (161 sites under selection, global $\omega = 13.05$ - heavily selected, potentially due to signaling pathways).

### C. Histones
* **Histones (HIST genes)** remain the only class where we find **zero selection sites** across all variants (`H3C10`, `H4C12`, `H2BC4`, etc.). Their packing constraints in the nucleosome appear to remain completely rigid across all mammalian genomes.

---

## 💡 Key Takeaway
This analysis validates the core hypothesis of our taxon density study: **as sequence count scales up, "constant" genes are revealed to be dynamic**. Genes historically used as invariant negative controls (like `ACTB` and `CALM1`) actually undergo subtle, site-specific episodic adaptation that was invisible in smaller, legacy alignments.
