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
- The complete protein sequence set is now deposited, but a family-level phylogenetically blocked HMM validation remains future work because it requires new train/test design, negative sets, profile fitting, and calibration rather than only sequence recovery.
- The repository and version-independent Zenodo concept DOI are public. Version 0.5.0 must still be released and archived after the tested files and checksum manifest are finalized.

`scripts/preflight_extended_validation.py` converts these gaps into machine-readable readiness checks. Its strict mode is expected to fail until the external inputs listed in the validation protocols are deposited.
