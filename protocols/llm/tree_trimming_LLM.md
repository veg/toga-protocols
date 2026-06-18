# LLM Protocol: Species Tree Pruning & Newick Manipulation

This document details the algorithmic logic and tree pruning rules for processing phylogenetic species trees.

---

## 1. Algorithmic Tree Pruning Logic
To prune a Newick species tree to a subset of representative leaf names:

1. **Leaf Selection**: Compile the list of allowed leaf IDs (NCBI assembly accessions or assembly directory names) from the representative mapping `species_to_assembly.csv`.
2. **Pruning Algorithm**:
   * Parse the Newick string into a tree data structure (using `ete3` or BioPython's `Phylo`).
   * Perform recursive pruning: for each internal node, if all its descendant leaves are not in the allowed set, prune the entire node.
   * If some descendant leaves are allowed and others are not, remove the disallowed leaf nodes and collapse any single-child internal nodes to maintain correct branch lengths (by summing the branch lengths of the collapsed parent and child).
3. **Format**: Save the final pruned tree in Newick/NHX format to `data/pruned_species_tree.nhx`.

---

## 2. Resolving Paraphyly & Monophyly Checks
For programmatic verification of monophyly:
* For each species in the metadata that has multiple assemblies in the original tree, find the Most Recent Common Ancestor (MRCA) node of those assembly leaf nodes.
* Check if any other species' leaves are descendants of this MRCA.
* If other species exist inside the MRCA clade, flag the target species as **non-monophyletic** and output the nested contaminants (as documented in `tree_pruning_report.md` for Capra, Ovis, Cervus, Giraffa, Balaenoptera, and Artibeus).

---

## 3. Rename & Suffix Allocation Rules
To map leaves to standardized 6-character taxonomic codes:

1. **Mapping**: Read the mapping of assemblies to 6-character codes from `species_to_taxon_map.csv`.
2. **Collision Detection & Renaming**:
   * Track assigned names in a set.
   * For each leaf node to be renamed:
     * If the target name (e.g. `odoHem`) is not in the set, rename the node and add the name to the set.
     * If the target name is already in the set, append a sequential suffix (e.g. `_2`, `_3`) until a unique name is formed (e.g. `odoHem_2`), then rename the node and add it to the set.
   * The final tree must contain only unique leaf labels to prevent parser errors in downstream tools (like HyPhy or model loaders).
