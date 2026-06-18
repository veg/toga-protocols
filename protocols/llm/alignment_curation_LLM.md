# LLM Protocol: Alignment Curation & Sequence Filtration

This protocol documents the algorithmic rules and thresholds for filtering multiple sequence alignments (MSAs) and sequence entries.

---

## 1. Sequence Exclusion Metrics & Rules
To filter sequence records within an MSA database or fasta/NEXUS file:

1. **Taxonomic Check**:
   * Let $S$ be the set of 742 allowed species names (from `species_to_taxon_map.csv`).
   * Discard any sequence where `species_name` is not in $S$.
2. **Coverage Calculation**:
   * Let $L$ be the alignment length (in nucleotides).
   * Let $N_{resolved}$ be the count of standard bases (`A`, `C`, `G`, `T`, case-insensitive).
   * Discard sequence if:
     $$\frac{N_{resolved}}{L} < 0.50$$
3. **Premature Stop Codon Check**:
   * Translate the in-frame codon sequence (ignoring gaps).
   * If any codon translates to a stop codon (`*`) before the final codon, discard the sequence.
4. **Deduplication Logic**:
   * If multiple sequences exist for the same species $s$:
     * Let $C_{best}$ be the coverage of the sequence from the chosen best assembly.
     * Let $C_{alt}$ be the coverage of any alternative assembly sequence.
     * If $C_{alt} \ge C_{best} + 0.10$, select the alternative sequence; otherwise, keep the chosen best assembly sequence.

---

## 2. Alignment Degapping & Minimum Verification
Once sequence exclusion is complete:

1. **Degapping Column Check**:
   * Iterate over all columns in the alignment matrix.
   * If a column $i$ contains only gap characters (`-`) across all remaining sequence rows, delete column $i$.
2. **Min-Depth Threshold**:
   * If the number of remaining sequence rows is $< 3$, discard the entire alignment file and record it as `dropped_few_seqs`.
