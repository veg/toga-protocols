# LLM Protocol: Evaluating Selection Contribution of Recent Genomes

This protocol specifies the queries, dates, and database criteria to measure the positive selection signal contributed by recently sequenced genomes.

---

## 1. Classification Criteria for Recent Genomes
An assembly is classified as "recently sequenced" if it meets the following criteria:
* **Submission Date**: NCBI Assembly submission date $\ge \text{2021-06-09}$ (within a 5-year window of the comparative study).
* **Taxonomic Mapping**: The assembly must map to one of the 742 representative species in `species_to_taxon_map.csv` (using the `Trimmed_Taxon_Name` key).

---

## 2. Selection Signal Calculation Algorithm
To measure the total selection signal $S_g$ contributed by a recently sequenced genome $g$ (where $g$ represents a species/leaf taxon):

1. **Substitution Extraction**:
   * Query the substitutions database to identify non-synonymous substitutions occurring on the terminal branch representing the species $g$ (i.e. where `branch_name = Trimmed_Taxon_Name` from `species_to_taxon_map.csv`).
2. **Quality Curation Filter**:
   * Intersect the substitution records with the site results database (`site_results` table).
   * Restrict selection counts to codon sites where the site-level classification is explicitly `'GOLD'` or `'SILVER'`.
3. **Aggregated Metric**:
   * Calculate $S_g = N_{GOLD} + N_{SILVER}$, where $N_{GOLD}$ and $N_{SILVER}$ are the total counts of substitutions on the terminal branch of species $g$ across all genes.

### Database Query Implementation:
To compute these counts across all processed genes, execute the following SQL:
```sql
SELECT 
    sub.branch_name AS taxon_code,
    sum(case when sit.classification = 'GOLD' then 1 else 0 end) AS gold_count,
    sum(case when sit.classification = 'SILVER' then 1 else 0 end) AS silver_count,
    count(*) AS total_selection_signal
FROM site_substitutions sub
JOIN site_results sit 
    ON sub.gene_name = sit.gene_name 
   AND sub.site_index = sit.site_index
WHERE sub.branch_name IN (SELECT Trimmed_Taxon_Name FROM species_to_taxon_map)
  AND sit.classification IN ('GOLD', 'SILVER')
GROUP BY sub.branch_name
ORDER BY total_selection_signal DESC;
```

---

## 3. High-Signal Lineage Attributes
For optimization and target sequencing selections, prioritize:
1. **Long Branch Lengths**: High cumulative terminal branch lengths in the species tree (which represent higher evolutionary distance and accumulate more mutations).
2. **High Scaffold/Contig N50**: Ensure assemblies are chromosome-level or high-scaffold N50 to minimize alignment errors and prevent the curation pipeline from filtering out the region.
