# Data gaps and non-reproducible stages

This document prevents unavailable evidence from being mistaken for an omission in the workflow instructions.

## Blocking raw inputs

The complete 1,069-genome candidate pool and raw candidate/canonical HMM hit tables were not present in the supplied archive. Consequently, the study cannot yet repeat domain discovery, reciprocal exclusion, and de novo clustering on alternative panels. The independent projection analysis validates transfer of fixed family HMMs, not rediscovery.

The corrected 3,526-genome presence/absence and copy-number matrices were also absent. The deposited matrices contain 3,192 genomes and are explicitly labeled legacy projection matrices. Corrected aggregate prevalence and genome-accounting tables remain available. No missing 334 genome columns have been invented.

## Missing analysis source tables

The following manuscript-era files were named but not deposited: `domain_validation.csv`, `enrichment_block_permutation.csv`, `enrichment_dependence_corrected.csv`, `motif_reanalysis_deep.csv`, `novel_family_distribution.csv`, `pfam_independent_validation.csv`, `regulon_motifs.csv`, `self_consistency.csv`, and `structural_all_hits_expanded.csv`.

Figures or summaries that depend on these files are treated as archived results, not regenerated evidence. `ADDITIONAL_DATA_SCHEMA.md` defines what must be deposited for a fully reproducible release.

## Analyses that require new work outside this repository

- Biochemical DNA-binding measurements, perturbation experiments, and transcriptomic confirmation require laboratory experiments.
- Family-level phylogenetically blocked HMM validation requires the complete protein sequence set for all 250 families; the current release does not contain it.
- A DOI requires an actual archival deposit. `CITATION.cff` and `ARCHIVE_DEPOSITION.md` prepare that action but do not claim a DOI.
- Public release of the repository requires an explicit visibility decision by the owners.
