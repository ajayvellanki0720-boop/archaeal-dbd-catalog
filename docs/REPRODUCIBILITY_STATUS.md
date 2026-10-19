# Reproducibility status

## Reproducible from this repository

- Verification of the 250-family and 13,237-member catalog summaries
- Verification of the 38,153 targeted-domain protein count
- Verification of 3,526-genome accounting outcomes
- Verification of motif, structure, and uncharacterized-family summary counts
- Inspection of the deposited HMM, IQ-TREE, and Foldseek archives
- The new 500-replicate within-panel resampling analysis
- The new lineage-stratified catalog accumulation analysis
- The new 500-replicate lineage-matched independent-genome projection analysis
- The new 30-replicate genome-level motif holdout analysis
- A 30-replicate genus-blocked motif transfer analysis with zero train/test genus overlap
- Source, assembly-level, and reported-quality sensitivity on the 3,192 accession-resolved projection genomes
- Claim-calibrated main figures and a workflow figure
- A minimal worked example, locked dependencies, container recipe, and checksum manifest
- Genus-clustered prevalence intervals, leave-one-major-lineage-out motif transfer, and ten-family sequence separability
- Frozen protocols and input gates for independent discovery, HMM calibration, external annotation, a discovery multiverse, and blind reproduction
- Automated regression tests for these quantities

## Present as released results but not independently regenerated here

- Final corrected prevalence summaries across 3,526 genomes
- Dependence-corrected genomic-neighborhood figure
- Deep motif reanalysis figure
- Structural reference comparisons
- IQ-TREE phylogenies

## Missing inputs or intermediates

- The complete 1,069-genome discovery candidate pool
- Raw candidate-profile HMM hits
- Raw canonical-profile HMM hits
- Corrected 3,526-genome presence-absence matrix
- Corrected 3,526-genome copy-number matrix
- `domain_validation.csv`
- `enrichment_block_permutation.csv`
- `enrichment_dependence_corrected.csv`
- `motif_reanalysis_deep.csv`
- `novel_family_distribution.csv`
- `pfam_independent_validation.csv`
- `regulon_motifs.csv`
- `self_consistency.csv`
- `structural_all_hits_expanded.csv`

## Legacy inconsistencies retained for transparency

- The raw projection metadata and matrices contain 3,192 genomes, whereas final accounting contains 3,526.
- The `got_hit` field in `genome_accounting.csv` was not updated for the 379 HMM-rescued genomes; the `fate` column contains the corrected classification.
- The original projection script does not implement the final rescue stage.
- The original evolutionary script uses an adjacent-tip heuristic, while the archived final analysis uses IQ-TREE and Fitch parsimony.
- The original synthesis text described the transcriptomic analysis as edgeR-style; the released analysis script instead performs simple tests on log-transformed CPM values.

## Interpretation

This repository improves transparency and makes the available checks repeatable. The independent projection analysis strengthens evidence that the family HMMs transfer to new accessions, and the motif holdout analyses test out-of-sample sequence prediction. Neither replaces de novo discovery on independent panels or biochemical validation. The release must not be described as a complete raw-to-figure reproduction until the missing inputs and corrected matrices are deposited and the external-tool workflow is rerun in a clean environment. It must not claim a public archive DOI until one has actually been minted.
