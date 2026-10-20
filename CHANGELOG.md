# Changelog

## 0.6.0 - 2026-10-20

- Added catalog-wide complete-genus separability analysis for all 13,237 deposited sequences; 13,166 sequences from 245 families were evaluable.
- Reframed the original random sequence holdout as an optimistic check and reported taxonomically blocked macro-, micro-, family-, and sequence-level results.
- Corrected the `got_hit` field for 379 HMM-rescued genomes to agree with final accounting fate.
- Added explicit search-universe provenance and mutually exclusive taxonomic-group mapping tables.
- Restricted accession-level prevalence inference to the deposited 3,192-genome matrices and retained the 3,526-genome values as archived aggregate results.
- Added quantitative phylogenetic Results, clarified structural denominators, and removed the unsupported transcriptomic narrative.
- Replaced confirmatory terminology with secondary validation and updated the line-numbered Scientific Reports package.

## 0.5.0 - 2026-10-19

- Added a 500-replicate lineage-stratified catalog-accumulation analysis across 10–100% of the discovery panel.
- Distinguished observation of fixed released families from retention under the five-genome eligibility rule.
- Integrated the analysis, limitations, supplementary figure, and reproducible outputs into the Scientific Reports submission package.
- Recovered and deposited all 13,237 catalog-member protein sequences with accession-, URL-, checksum-, date-, and fallback-level provenance.
- Added automated sequence-completeness and sequence-length verification.
- Calibrated genomic-neighborhood claims to the deposited exploratory evidence because the dependence-corrected source table was not preserved.
- Updated the line-numbered manuscript, Supplementary Information, cover letter, and point-by-point reviewer response for Scientific Reports.

## 0.4.0 - 2026-10-18

- Added genus-clustered prevalence estimates and bootstrap intervals.
- Added leave-one-major-lineage-out motif validation with same-lineage controls.
- Added leave-one-sequence-out family separability for the released ten-family subset.
- Added preregistered protocols and machine-readable readiness gates for analyses requiring missing data, new database snapshots, independent analysts, or laboratory experiments.

## 0.3.0 - 2026-10-18

- Added genus-blocked motif validation with zero train/test genus overlap.
- Added source, assembly-level, and reported-quality sensitivity analyses.
- Added claim-calibrated main figures and a workflow figure.
- Added a minimal worked example, container recipe, locked Python requirements, provenance documentation, and SHA-256 manifest generation.
- Added a supplementary methods document, reporting checklist, experimental-validation plan, cover letter, and response-to-reviewers draft.
- Replaced the manuscript's most interpretive figures and incorporated the stricter validation results and submission declarations.

## 0.2.0

- Added 500 lineage-matched resamples from independent projection genomes, with explicit exclusion of discovery-panel accessions.
- Added genome-level 70/30 motif holdout validation against dinucleotide-shuffled and GC-matched other-family controls, with descriptive resampling summaries that avoid treating overlapping splits as independent observations.
- Added versioned replicate-level results, summaries, and supplementary figures for both analyses.
- Extended the command-line interface, tests, documentation, and manuscript to report the new evidence and its limits.

## 0.1.0

- Added a portable Python package and command-line interface.
- Added automated release-integrity tests and GitHub Actions configuration.
- Added deterministic within-panel sensitivity analysis with 500 replicates.
- Added an executable definition of reciprocal domain exclusion and its threshold grid.
- Added a complete reproducibility-status inventory.
- Preserved the original scripts under `legacy/code` without presenting them as portable.
- Added the revised manuscript with bounded claims and corrected method descriptions.
