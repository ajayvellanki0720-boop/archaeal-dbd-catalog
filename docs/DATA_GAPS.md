# Data gaps and non-reproducible stages

This document prevents unavailable evidence from being mistaken for an omission in the workflow instructions.

## Blocking raw inputs

The complete 1,069-genome candidate pool, original nineteen-family/twenty-profile reference inventory, 747 candidate-profile accessions, exact Pfam release, and raw candidate/reference HMM hit tables were not present in the supplied archive. Consequently, the study cannot repeat domain discovery, reciprocal exclusion, or de novo clustering on alternative panels. The independent projection analysis validates transfer of fixed family HMMs, not rediscovery.

The corrected 3,526-genome presence/absence and copy-number matrices were also absent. The deposited matrices contain 3,192 genomes and are explicitly labeled legacy projection matrices. Corrected aggregate prevalence and genome-accounting tables remain available. No missing 334 genome columns have been invented.

## Missing analysis source tables

The following manuscript-era files were named but not deposited: `domain_validation.csv`, `enrichment_block_permutation.csv`, `enrichment_dependence_corrected.csv`, `motif_reanalysis_deep.csv`, `novel_family_distribution.csv`, `pfam_independent_validation.csv`, `regulon_motifs.csv`, `self_consistency.csv`, and `structural_all_hits_expanded.csv`.

Figures or summaries that depend on these files are treated as archived results, not regenerated evidence. `ADDITIONAL_DATA_SCHEMA.md` defines what must be deposited for a fully reproducible release.

## Analyses that require new work outside this repository

- Biochemical DNA-binding measurements, perturbation experiments, and transcriptomic confirmation require laboratory experiments.
- Complete-genus composition validation is now available for all deposited sequences. A taxonomically blocked HMM calibration with explicit negatives and score thresholds remains future work.
- The repository and version-independent Zenodo concept DOI are public. Version 0.5.0 is archived at DOI 10.5281/zenodo.23496681; later corrections are versioned separately.

`scripts/preflight_extended_validation.py` converts these gaps into machine-readable readiness checks. Its strict mode is expected to fail until the external inputs listed in the validation protocols are deposited.
