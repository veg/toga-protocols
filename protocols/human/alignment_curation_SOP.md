# Standard Operating Procedure (SOP): Alignment and Sequence Curation

This protocol outlines how to filter multiple sequence alignments (MSAs) and exclude low-quality or redundant sequences before running HyPhy selection analyses (MEME/aBSREL).

---

## 1. Gene-Level Quality Filters
Before examining individual sequences, exclude entire genes if they meet any of the following criteria:
* **Outlier Genes**: Genes with anomalous codon lengths or sequence counts.
* **Pseudogenes**: Genes explicitly annotated as pseudogenes or having inactive/aberrant characteristics.
* **Insufficient Taxonomic Depth**: Genes containing fewer than **3 valid species sequences** after applying sequence-level curation.

---

## 2. Sequence-Level Selection Filters
For the remaining genes, filter out low-quality or redundant sequences according to the following strict criteria:

1. **Species-Level Taxonomic Restriction**: Retain only sequences belonging to the **742 species** present in the pruned species tree (`species_to_taxon_map.csv`).
2. **Redundant Assembly Deduplication**: If multiple genome assemblies exist for a single species, select the representative sequence from the "best" assembly (according to the assembly selection SOP).
   * *Exception*: If a non-representative assembly has a **$\ge 10\%$ higher sequence coverage** (fewer gaps) than the representative one, select it instead.
3. **Coverage Threshold**: Discard any sequence where resolved nucleotides (A, C, G, T) make up **$< 50\%$** of the alignment length.
4. **Premature Stop Codons**: Discard any sequence containing an in-frame stop codon occurring before the final codon.

---

## 3. Post-Curation Alignment Cleaning
After filtering out disallowed sequences:
* **Column Degapping**: Check all alignment columns. If a column contains only gap characters (`-`) across all remaining sequences, remove the entire column to maintain a compact, cleaned alignment.
* **Min-Species Re-Evaluation**: Re-verify that the final alignment still contains **$\ge 3$ sequences**. If not, drop the entire gene from molecular evolutionary analysis.
