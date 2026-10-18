# Projection matrix scope

| Artifact | Families | Genomes | Status | Permitted use |
|---|---:|---:|---|---|
| `release/data/presence_absence_matrix.csv` | 250 | 3,192 | Legacy, accession-resolved | Independent projection resampling and quality sensitivity |
| `release/data/copynumber_matrix.csv` | 250 | 3,192 | Legacy, accession-resolved | Legacy copy-number summaries |
| `release/tables/corrected_prevalence.csv` | 250 | 3,526 denominator | Final aggregate summary | Final group-level prevalence reporting only |
| `release/tables/genome_accounting.csv` | — | 3,526 | Final accounting | Fate and denominator accounting |

The two scopes must not be merged as though they were the same matrix. The 3,192-genome analyses are explicitly sensitivity analyses of the accession-resolved legacy projection. The 3,526-genome values are final aggregate summaries whose underlying corrected matrix has not been deposited.
