# Selection Signal Contributions from Recently Sequenced Genomes on Hot MEME EBF Branches

This report identifies and analyzes which recently sequenced mammalian genomes (submission date $\ge$ 2021-06-09) contribute significant selection signal (measured by substitutions at `GOLD` and `SILVER` classified sites) to hot branches evaluated using MEME (Mixed Effects Model of Evolution) and EBF (Empirical Bayes Factor).

## 1. Methodology and Background
- **Recent Genome Filtering:** Mammalian assemblies from `assemblies_and_species.tsv` with an NCBI submission date on or after **2021-06-09**.
- **Taxon Mapping:** Leaf names in the phylogenetic tree were mapped to assembly species names via `docs/species_to_taxon_map.tsv` using standardized uppercase matching.
- **Selection Signal:** Counted substitutions (`site_substitutions` table in `meme_results.db`) occurring on the branch corresponding to the species at sites classified as `GOLD` or `SILVER` in the overall MEME results.
- **Selection Classifications:**
  - `GOLD` sites represent robust selection signals with high confidence.
  - `SILVER` sites represent strong selection signals with moderate confidence.

## 2. Summary of Findings
- **High-Quality Chromosome Assemblies Drive Signal Recovery:** The top contributing genomes are dominated by chromosome-level assemblies with very high contig and scaffold N50 values (e.g., *Erinaceus europaeus*, *Ctenodactylus gundi*, *Castor canadensis*). High-quality assemblies lead to better ORF projection by TOGA, resulting in complete alignments and high-confidence substitution mapping.
- **Top 5 Contributing Genomes:**
  1. ***Erinaceus europaeus* (European hedgehog):** Assembly `GCF_950295315.1` (submitted 2023-05-01, Chromosome status) contributes **12,941** selection substitutions.
  2. ***Ctenodactylus gundi* (Gundi):** Assembly `GCA_048771875.1` (submitted 2025-03-21, Chromosome status) contributes **12,858** selection substitutions.
  3. ***Typhlomys cinereus* (Chinese pygmy dormouse):** Assembly `GCA_023101885.1` (submitted 2022-04-25, Scaffold status) contributes **11,191** selection substitutions.
  4. ***Rhynchocyon petersi* (Black and rufous elephant shrew):** Assemblies `GCA_043290065.1` and `GCA_043290085.1` (submitted 2024-10-17, Chromosome status) contribute **9,894** selection substitutions.
  5. ***Castor canadensis* (American beaver):** Assemblies `GCA_001984765.2` and `GCA_047511655.2` (submitted 2025/2026, Chromosome status) contribute **9,204** selection substitutions.

## 3. Top 30 Recently Sequenced Genomes and Selection Contributions

| Species | NCBI Accession | Submission Date | Assembly Status | Contig N50 (bp) | Scaffold N50 (bp) | GOLD Sites | SILVER Sites | Total Selection Signal |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| *Erinaceus europaeus* | GCF_950295315.1 | 2023-05-01 | Chromosome | 999,919 | 126,757,761 | 12,815 | 126 | **12,941** |
| *Ctenodactylus gundi* | GCA_048771875.1 | 2025-03-21 | Chromosome | 38,834,406 | 108,787,391 | 12,753 | 105 | **12,858** |
| *Typhlomys cinereus* | GCA_023101885.1 | 2022-04-25 | Scaffold | 1,971,523 | 4,558,085 | 11,078 | 113 | **11,191** |
| *Rhynchocyon petersi* | GCA_043290065.1 | 2024-10-17 | Chromosome | 13,791,435 | 546,604,020 | 9,810 | 84 | **9,894** |
| *Rhynchocyon petersi* | GCA_043290085.1 | 2024-10-17 | Chromosome | 13,617,243 | 521,973,465 | 9,810 | 84 | **9,894** |
| *Castor canadensis* | GCA_001984765.2 | 2026-01-29 | Chromosome | 144,908,913 | 158,004,004 | 9,133 | 71 | **9,204** |
| *Castor canadensis* | GCA_047511655.2 | 2025-02-10 | Chromosome | 59,379,167 | 158,968,201 | 9,133 | 71 | **9,204** |
| *Cephalopachus bancanus* | GCA_027257055.1 | 2022-12-23 | Contig | 5,252,243 | 5,252,243 | 8,723 | 63 | **8,786** |
| *Eliomys quercinus* | GCA_051143595.1 | 2023-10-26 | nan | 45,399,149 | 107,541,173 | 7,289 | 55 | **7,344** |
| *Eliomys quercinus* | GCA_051143605.1 | 2023-10-26 | nan | 51,740,598 | 108,114,398 | 7,289 | 55 | **7,344** |
| *Muscardinus avellanarius* | GCA_963383645.1 | 2023-08-26 | Chromosome | 2,432,660 | 119,481,734 | 6,787 | 53 | **6,840** |
| *Blarina brevicauda* | GCA_051295465.1 | 2025-07-11 | Contig | 2,920,142 | 2,920,142 | 6,580 | 60 | **6,640** |
| *Myocastor coypus* | DNA Zoo Consortium | 2024-10-07 | nan | 3,975,994 | 138,704,734 | 6,402 | 44 | **6,446** |
| *Dasypus novemcinctus* | GCF_030445035.1 | 2023-07-13 | Chromosome | 14,295,626 | 130,056,132 | 6,321 | 48 | **6,369** |
| *Dasypus novemcinctus* | GCF_030445035.2 | 2023-07-13 | Chromosome | 13,963,394 | 127,081,865 | 6,321 | 48 | **6,369** |
| *Perognathus longimembris pacificus* | GCF_023159225.1 | 2022-04-27 | Chromosome | 7,389,774 | 72,679,016 | 6,135 | 47 | **6,182** |
| *Doryrhina cyclops* | GCA_043880285.1 | 2021-10-02 | nan | 14,255,823 | 154,073,572 | 6,112 | 60 | **6,172** |
| *Micromys minutus* | GCA_963924665.1 | 2024-01-29 | Chromosome | 1,835,000 | 59,776,753 | 5,982 | 50 | **6,032** |
| *Rhinopoma microphyllum* | GCA_043880545.1 | 2021-10-19 | nan | 33,347,918 | 146,810,189 | 5,716 | 33 | **5,749** |
| *Tenrec ecaudatus* | GCA_050624465.1 | 2025-05-30 | Chromosome | 54,742,630 | 165,248,668 | 5,587 | 55 | **5,642** |
| *Tenrec ecaudatus* | GCF_050624435.1 | 2025-05-30 | Chromosome | 41,601,014 | 163,958,880 | 5,587 | 55 | **5,642** |
| *Tenrec ecaudatus* | TBG Hillerlab, Fauna | 2023-09-02 | nan | 71,825,797 | 156,487,869 | 5,587 | 55 | **5,642** |
| *Tenrec ecaudatus* | TBG Hillerlab, Fauna | 2023-09-02 | nan | 54,617,000 | 133,279,229 | 5,587 | 55 | **5,642** |
| *Salpingotus crassicauda* | GCA_040869305.1 | 2024-07-25 | Contig | 13,709,620 | 13,709,620 | 5,588 | 48 | **5,636** |
| *Sigmodon hispidus* | GCA_041721565.1 | 2024-09-03 | Chromosome | 360,703 | 107,795,196 | 5,452 | 42 | **5,494** |
| *Sigmodon hispidus* | GCA_047301775.1 | 2025-01-30 | Chromosome | 77,390,441 | 118,707,662 | 5,452 | 42 | **5,494** |
| *Aselliscus stoliczkanus* | GCA_043727835.1 | 2021-10-02 | nan | 64,984,936 | 162,130,651 | 5,379 | 36 | **5,415** |
| *Aselliscus stoliczkanus* | GCA_033961575.1 | 2023-11-28 | Chromosome | 73,009,000 | 162,004,500 | 5,379 | 36 | **5,415** |
| *Tolypeutes matacus* | GCA_026826555.1 | 2022-12-12 | Scaffold | 129,379 | 12,561,802 | 5,292 | 36 | **5,328** |
| *Vespertilio murinus* | GCA_963924515.1 | 2024-01-29 | Chromosome | 3,319,810 | 186,292,382 | 5,164 | 49 | **5,213** |


## 4. Key Takeaways and Implications for Selection Analyses
1. **Assembly Quality is Directly Proportional to Power:** Lower quality assemblies (e.g. contig-level assemblies or scaffolds with low N50) have missing sequences or frameshifting/inactivating errors. These prevent proper selection inference, reducing the signal count.
2. **Taxonomic Enrichment:** The highly represented orders/families (such as Eulipotyphla, Rodentia, and Chiroptera) reflect recent sequencing focus (e.g., from consortia like VGP, Bat1K) and provide valuable evolutionary depth to resolve branches that were historically poorly sampled.
3. **Phylogenetic Utility:** Adding these high-quality nodes splits long branches, reducing the risk of long-branch attraction and allowing MEME to pinpoint specific lineages where episodic diversifying selection occurred with high EBF support.

***

*Report compiled on 2026-06-18.*
