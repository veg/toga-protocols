# Standard Operating Procedure (SOP): Contribution of Recently Sequenced Genomes

This standard operating procedure describes how to evaluate and monitor the impact of recently sequenced mammalian genomes on our positive selection detection pipeline.

---

## 1. Context and Objective
Large-scale vertebrate sequencing projects (such as the Vertebrate Genomes Project and Bat1K) have dramatically expanded the taxonomic density of mammalian alignments. By introducing high-quality genome assemblies for previously under-sampled or highly divergent branches of the mammalian tree, these assemblies provide new mutational depth that allows the HyPhy MEME pipeline to detect historically hidden selection signals.

We define "recently sequenced genomes" as those with NCBI assembly submission dates within the last 5 years (specifically, on or after **June 9, 2021**).

---

## 2. Measuring Selection Signal Contribution
A genome's contribution to positive selection detection is measured by aggregating the non-synonymous substitutions occurring on its terminal branch that meet strict quality curation standards:

1. **Terminal Branch Location**: Substitutions must map specifically to the leaf branch leading to that species.
2. **Curation Tier Validation**: The substitutions must occur at codon sites classified as **`GOLD`** (high-confidence single nucleotide changes) or **`SILVER`** (supported multi-nucleotide changes).

---

## 3. High-Contributing Lineages
Analysis of the dataset shows that the top-contributing recent genomes share two primary characteristics:
* **High Evolutionary Divergence**: Species situated on long, isolated terminal branches (such as Gundi, Hedgehog, and Elephant Shrew) accumulate more mutations over evolutionary time, yielding more selection statistics.
* **High Assembly Quality**: Chromosome-level assemblies allow for cleaner alignments with fewer gaps or sequencing artifacts. This prevents the curation pipeline from discarding these positions and results in higher numbers of valid selection calls.

Key groups that dominate selection contributions include:
* **Rodents** (Gundis, Dormice, Beavers, Jerboas) due to high basal mutation rates.
* **Bats and Insectivores** (Hedgehogs, Shrews, microbats) reflecting the dense evolutionary sampling of these clades.
