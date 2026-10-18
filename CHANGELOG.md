# Changelog

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
