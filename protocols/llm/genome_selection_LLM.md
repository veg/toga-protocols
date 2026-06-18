# LLM Protocol: Genome Assembly Selection & Quality Control

This document provides system context and sorting logic rules for programmatically executing genome assembly selection and quality filtering.

---

## 1. Database & Metadata Input Format
The assembly metadata is compiled from NCBI Assembly databases and TOGA gene projection summaries.
Expected attributes per assembly record:
* `species_name` (TEXT)
* `assembly_accession` (TEXT, e.g. GCA_...)
* `assembly_level` (TEXT: `'Chromosome'`, `'Complete Genome'`, `'Scaffold'`, `'Contig'`, or `'Unknown'`)
* `scaffold_n50` (INTEGER)
* `contig_n50` (INTEGER)
* `intact_orfs` (INTEGER/FLOAT)
* `missing_seqs` (INTEGER)
* `inactivating_muts` (INTEGER)
* `submission_date` (TEXT, ISO format or YYYY-MM-DD)

---

## 2. Selection Hierarchy Sorting Rules (Python/Pandas Logic)
To programmatically select the "best" representative assembly per species:

1. **Mapping Priority**: Map `assembly_level` strings to numeric scores:
   ```python
   level_priority = {
       'Chromosome': 4,
       'Complete Genome': 4,
       'Scaffold': 3,
       'Contig': 2,
       'Unknown': 1
   }
   ```
2. **Multi-Key Sorting**: Group records by `species_name` and sort each group using pandas `sort_values` with ascending priorities (ascending=False for highest values):
   ```python
   df_sorted = df.sort_values(
       by=[
           'level_score', 
           'scaffold_n50', 
           'contig_n50', 
           'intact_orfs', 
           'submission_date'
       ], 
       ascending=[False, False, False, False, False]
   )
   ```
3. **Deduplication**: Keep the first record per species:
   ```python
   best_representation_df = df_sorted.drop_duplicates(subset=['species_name'], keep='first')
   ```

---

## 3. Curation Flag Filter Thresholds
Genomes must be marked as `problematic` if they satisfy any of the following logical predicates:

```python
is_problematic = (
    (df['assembly_level'] == 'Contig') |
    (df['scaffold_n50'] < 100000) |
    (df['contig_n50'] < 10000) |
    (df['intact_orfs'] < 8000) |
    (df['missing_seqs'] > 2000) |
    (df['inactivating_muts'] > 2000)
)
```
Flagged genomes should be output to `problematic_genomes.tsv` for manual vetting or systematic exclusion from downstream multi-species alignments.
