# LLM Protocol: Legacy vs. TOGA Alignment Mapping & Comparative Selection

This protocol specifies the algorithmic mapping, filters, and criteria for comparing selection detection between legacy alignments (~45 species) and dense TOGA alignments.

---

## 1. Algorithmic Coordinate Mapping
To map codon sites between a legacy alignment and a modern TOGA alignment without running pairwise sequence alignments:

1. **Species Set Definition**:
   * Define the set of 48 legacy species leaf names:
     $$S_{legacy} = \{\text{'hg'}, \text{'panTro'}, \text{'gorGor'}, \dots, \text{'tupChi'}\}$$
2. **Sub-alignment Extraction**:
   * Parse the TOGA NEXUS file (format: `{gene_name}.gz`) to retrieve sequence data for all leaves.
   * Retain only sequences matching species in $S_{legacy}$.
3. **All-Gap Column Pruning**:
   * Let $C$ be the sequence matrix of the legacy subset.
   * For each codon index $c$ (columns of length 3):
     * If the codon columns consist entirely of gap or missing characters (`-`, `?`) across all species in $S_{legacy}$:
       * Discard the column.
     * Else:
       * Append $c$ to the list of mapped indices `keep_codon_indices`.
4. **Site Index Coordinate Conversion**:
   * The 0-based legacy alignment index $i$ maps to the 1-based TOGA alignment index $t$:
     $$t = \text{keep\_codon\_indices}[i] + 1$$

Implement in [analyze_legacy_only_sites.py](file:///Users/sergei/Dropbox/TOGA2026/artifacts/scratch/analyze_legacy_only_sites.py).

---

## 2. Selection Discrepancy Classification

### A. Identification of Legacy-Only Sites
To compile selection sites that are significant only in the legacy background:
* **Legacy Filter**: Legacy p-value $\le 0.01$.
* **TOGA Filter**: Mapped TOGA p-value $> 0.05$ (or missing/pruned in the TOGA alignment).

### B. Classification Logic
Analyze the ancestral state reconstructions at the target site:

1. **Category A (Reconstruction Error / LRT Spike)**:
   * **Rule**: Mapped TOGA LRT $\approx 0.0$ AND Legacy LRT $> 10.0$ with a high inferred positive selection rate $\beta_+ \ge 50.0$.
   * **Validation**: Inspect the ancestral state at the parent nodes of mutated terminal branches in the legacy tree. Verify if the inclusion of sister taxa in TOGA corrects the ancestral state reconstruction from a non-synonymous mutation to a synonymous change (as seen in AOX1 Site 744 vs Site 697).
   
2. **Category B (Stabilization & Localized Dilution)**:
   * **Rule**: Mapped TOGA LRT $> 1.0$ with $0.05 < \text{TOGA p-value} \le 0.15$.
   * **Validation**: Inspect if the non-synonymous to synonymous substitution ratio at the site is stabilized to the background neutral expectation ($dN/dS \approx 1.0$) under the dense TOGA alignment (as seen in AOX1 Site 714 vs Site 667).
