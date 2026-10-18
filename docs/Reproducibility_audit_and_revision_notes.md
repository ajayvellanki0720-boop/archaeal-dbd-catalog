# Reproducibility audit and revision notes

## What was inspected

The linked Google Drive package was downloaded and audited in full. It contains 59 files: 14 Python scripts, 10 data archives or matrices, 19 result tables, and 16 figures. All scripts were read, all CSV and JSON outputs were checked for dimensions and internal consistency, both compressed analysis archives were inventoried, and all figures were visually reviewed.

## Changes made to the manuscript

- Reframed the paper as a bounded catalog of proteins within Pfam clans CL0123 and CL0057, rather than a comprehensive map of all archaeal DNA binding proteins or 250 confirmed transcription factor families.
- Rebuilt the abstract and narrative around one central result: the two target clans contain coherent archaeal family diversity beyond a nineteen family reference set, but most of that diversity resolves to known sequence or structural classes.
- Defined the 331 genome panel algorithm exactly. All five underrepresented groups were retained, Euryarchaeota was capped at 140, Crenarchaeota at 70, and proportional within group sampling used seed 42. The manuscript now states that these caps are pragmatic rather than biologically optimal.
- Added a new 500 replicate within panel resampling analysis. Moderate 80 percent lineage stratified resampling detected all 250 families in every replicate, but smaller panels reduced the number satisfying the five genome eligibility threshold.
- Defined the reciprocal exclusion rule precisely. The “regulatory region” in the original text was actually the aligned candidate domain interval. A canonical hit excludes a candidate only when it overlaps more than half of that interval and scores within 5 bits of the candidate. Multiple candidates are resolved by the highest scoring retained candidate.
- Identified the 5 bit value as a heuristic score tolerance and stated that a full sensitivity analysis is impossible from the released files because the raw candidate and canonical hit tables are absent.
- Replaced “predicted regulons” with genomic neighborhood hypotheses except when discussing legacy file names. The direct regulation claim was removed because flagship motifs generally did not enrich upstream of predicted targets.
- Separated evidence for a DNA binding domain protein from evidence for a transcription factor.
- Reduced the main visual narrative to six figures, with the new panel sensitivity plot included.
- Corrected the Haloferax volcanii method description. The released script uses simple tests on transformed counts with three replicates per group; it does not implement edgeR style dispersion estimation. The manuscript now treats this analysis as exploratory.
- Added a limitations section that records the archive's reproducibility gaps and legacy projection matrices.

## New panel sensitivity result

| Resampling scheme | Genomes | Families detected at least once, median and 95 percent interval | Families present in at least five genomes, median and 95 percent interval |
|---|---:|---:|---:|
| 80 percent within every lineage | 264 | 250, 250 to 250 | 212, 204 to 220 |
| All rare lineages and 80 percent of major caps | 289 | 250, 250 to 250 | 224, 218 to 231 |
| Maximum 30 genomes per lineage | 157 | 242, 235 to 247 | 137, 127 to 148 |

This is a retention test using fixed released family memberships. It is not a full rediscovery analysis from alternative draws of the original candidate genome universe.

## Package findings that require repair before a reproducibility claim

### Portability

- Four scripts depend on an internal `host.artifacts` runtime.
- Eight scripts contain absolute paths under `/Users/ajayvellanki/.claude-science/`.
- No script exposes a command line interface or a single documented workflow entry point.
- The package does not pin tool and library versions in an environment file.

### Files named in the manuscript but absent from the release

- `domain_validation.csv`
- `enrichment_block_permutation.csv`
- `enrichment_dependence_corrected.csv`
- `motif_reanalysis_deep.csv`
- `novel_family_distribution.csv`
- `pfam_independent_validation.csv`
- `regulon_motifs.csv`
- `self_consistency.csv`
- `structural_all_hits_expanded.csv`

### Projection inconsistencies

- `projection_genome_metadata.csv`, `presence_absence_matrix.csv`, and `copynumber_matrix.csv` cover 3,192 genomes, while the final accounting and corrected composition tables cover 3,526 genomes.
- `genome_accounting.csv` lists 3,144 DIAMOND hit genomes, 379 HMM rescued genomes, and 3 zero hit genomes, but its legacy `got_hit` field remains false for all 382 initially negative genomes.
- `corrected_prevalence.csv` differs from `presence_prevalence.csv`; the corrected table should be the sole source for final prevalence claims.

### Method and text mismatches

- The released projection script uses DIAMOND and does not itself reproduce the final HMM rescued 3,526 genome projection.
- The released evolutionary script uses an adjacent tip heuristic, while the archive also contains the later IQ TREE and Fitch parsimony outputs used in the revision.
- The released transcriptomic script performs simple two sample tests after log transformed CPM normalization. The synthesis script incorrectly describes this as edgeR style dispersion estimation.
- The released genomic neighborhood scripts do not include all dependence corrected source tables, although the corresponding summary figure is present.

## Recommended release repairs

1. Deposit the complete 1,069 genome candidate pool and the raw canonical and candidate domain hit tables.
2. Add every missing table listed above or remove claims that cannot be regenerated from deposited data.
3. Regenerate both presence absence and copy number matrices for all 3,526 genomes and remove or clearly label the 3,192 genome legacy matrices.
4. Replace absolute paths and internal runtime calls with relative paths and ordinary file input and output arguments.
5. Add a locked environment file, exact external database versions, checksums, and a one command workflow.
6. Add automated checks that reconcile genome counts, family counts, denominators, and figure source tables.
7. Rerun the full discovery workflow on multiple independently sampled panels from the complete candidate pool. The new resampling analysis is useful but cannot answer this stronger robustness question.

## Bottom line

The archive is sufficient to inspect the main family catalog and verify many reported summaries, but it is not yet a complete, portable end to end reproducibility package. The revised manuscript now states that boundary explicitly and avoids claims that rely on unavailable intermediates.
