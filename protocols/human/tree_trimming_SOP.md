# Standard Operating Procedure (SOP): Species Tree Pruning & Taxonomic Mapping

This protocol outlines how to prune the mammalian species tree to retain only the best representative assemblies, handle non-monophyletic species clades, and format taxon leaf names.

---

## 1. Pruning Goals
* **Target**: Reduce the tree from the initial set of **882 genome assemblies** to the core set of **742 representative species** assemblies.
* **Method**: Keep only the leaf nodes that match the chosen assemblies listed in `species_to_assembly.csv`. All other leaves are systematically pruned using a Newick-parsing library (e.g. `Bio.Phylo` in BioPython or `ETE3`).

---

## 2. Resolving Non-Monophyletic Species
During pruning, verify that assemblies representing the same species form monophyletic groups in the unpruned tree. We flagged **6 species** as non-monophyletic:
* `Capra hircus` (Domestic Goat) nested with `Capra aegagrus` (Wild Goat).
* `Ovis aries` (Domestic Sheep) nested with `Ovis orientalis` (Mouflon).
* `Cervus albirostris` (White-lipped Deer) nested with several other `Cervus` species.
* `Giraffa tippelskirchi` (Masai Giraffe) nested with `Giraffa camelopardalis` (Northern Giraffe).
* `Balaenoptera ricei` (Rices Whale) nested with `Balaenoptera edeni` (Brydes Whale).
* `Artibeus jamaicensis` nested with `Artibeus lituratus`.

**Action**: For these species, select the representative assembly carefully according to the assembly selection protocol (preferring the scaffold/chromosome level and scaffold N50) and prune the nested/domestic contaminant branches out.

---

## 3. Formatting Leaf Names & Handling Collisions
Leaf names in the final tree must be formatted using their standardized 6-character taxonomic codes (e.g., `hg38` for human, `odoHem` for mule deer). 

Because multiple subspecies or highly related assemblies can map to the same 6-character code, name collisions can occur during renaming. To prevent duplicates in the tree topology:
* Detect name collisions (e.g., `odoHem`, `equQua`, `gorGor`).
* Append unique sequential suffixes to the colliding leaf names in the tree (e.g., `odoHem_2`, `equQua_2`, `gorGor_2`).
* Maintain a lookup map linking these suffixed tree labels back to their original NCBI accessions in `species_to_taxon_map.csv`.
