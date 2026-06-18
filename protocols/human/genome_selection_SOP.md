# Standard Operating Procedure (SOP): Genome Assembly Selection & Quality Control

This protocol outlines how to select the best representative genome assembly per species and identify problematic assemblies that should be excluded from downstream selection analyses.

---

## 1. Assembly Selection Hierarchy
When multiple genome assemblies are available for a given species, the single "best" representative assembly is selected using a strict hierarchical scoring model:

1. **Assembly Level**: Select based on scaffolding quality:
   * `Chromosome` or `Complete Genome` (highest priority)
   * `Scaffold`
   * `Contig`
   * `Unknown` (lowest priority)
2. **Scaffold N50**: If levels are equal, prefer the assembly with the larger Scaffold N50 (in base pairs).
3. **Contig N50**: If Scaffold N50 is equal, prefer the assembly with the larger Contig N50 (in base pairs).
4. **TOGA Projections**: If N50 values are equal, select the assembly with the higher count of intact ORFs (Open Reading Frames) projected.
5. **Submission Date**: As a final tiebreaker, select the newest assembly by NCBI submission date.

---

## 2. Problematic Genome Curation (Flags)
To prevent assembly artifacts from contaminating selection analyses, genomes meeting any of the following criteria are flagged as **potentially problematic**:

* **Lack of Scaffolding**: Explicit assembly status is `Contig`.
* **Low Scaffold N50**: Scaffold N50 is under **100 KB** ($100,000$ bp).
* **Low Contig N50**: Contig N50 is under **10 KB** ($10,000$ bp).
* **High ORF Incompleteness**: Fewer than **8,000** intact ORFs projected by TOGA.
* **High Sequence Fragmentation**: Missing sequence count exceeds **2,000**.
* **High Inactivating Mutation Frequency**: More than **2,000** inactivating coding mutations (indels/frameshifts) detected.

---

## 3. Mapping Outputs
The selection script outputs three curated TSV files (located in the repository root and copied to `docs/`):
* `species_to_assembly.tsv`: One-to-one mapping from unique species to their chosen best representative assembly.
* `assembly_to_species.tsv`: Reference lookup matching assembly IDs back to their species names.
* `problematic_genomes.tsv`: List of all flagged assemblies along with their failed metrics and flag reasons.
