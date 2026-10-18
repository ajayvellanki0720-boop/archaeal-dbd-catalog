# Noncanonical archaeal DNA binding domain catalog

This repository accompanies the manuscript **A Catalog of Noncanonical HTH and RHH DNA Binding Domain Families in Archaea**. It organizes the released family catalog, records its evidentiary limits, and provides portable checks, panel-sensitivity analyses, independent-genome projection validation, and held-out motif validation.

## Scope

The study is a bounded survey of proteins detected through Pfam clans CL0123 and CL0057 after removing a nineteen-family reference set. It is **not** a comprehensive census of all archaeal DNA-binding proteins, and a family assignment is **not** proof that a protein is a transcription factor.

## Key findings represented in this release

- 38,153 proteins in the 331-genome panel carried a targeted clan domain.
- 16,535 were assigned to the reference set, 14,451 to another characterized family, and 7,167 were uncharacterized.
- Clustering and filtering defined 250 candidate families containing 13,237 proteins.
- Held-out family assignment recovered 97.3 percent of withheld members.
- In 500 lineage-matched samples of 331 projection genomes with no discovery-panel accessions, a median of 245 families were detected at least once and 206 occurred in at least five genomes.
- Of 79 Pfam-uncharacterized families, 56 were resolved by CDD and 14 more by structural comparison; nine remained unresolved by both resources tested.
- Eleven of twelve tested upstream motifs passed all three sequence null models.
- In genome-level 70/30 holdout tests, 10 of 12 motifs discriminated family upstreams from GC-matched upstreams of other families in at least 26 of 30 resamples.
- Genomic-neighborhood results are functional hypotheses, not confirmed regulons.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
make test
make all
```

Generated audit and validation outputs are written to `results/generated/`. The exact outputs used in the manuscript are versioned in `results/`.

## Repository layout

```text
src/archaeal_dbd/     Portable validation and sensitivity code
tests/                 Automated checks for headline invariants
release/data/          Released matrices and analysis archives
release/tables/        Released result tables
release/figures/       Released figures
legacy/code/           Original scripts preserved without claiming portability
results/               Versioned sensitivity and validation results
manuscript/            Revised manuscript
docs/                  Study design and reproducibility limitations
```

## Reproducibility status

`make all` reproduces the release audit, within-panel sensitivity analysis, independent-genome projection analysis, and genome-level motif holdout analysis. It does not rerun raw domain discovery because the released archive lacks the complete 1,069-genome candidate pool and raw canonical and candidate HMM hit tables. The 3,526-genome corrected prevalence summaries are present, but the raw presence and copy-number matrices cover the earlier 3,192-genome pass.

See [docs/REPRODUCIBILITY_STATUS.md](docs/REPRODUCIBILITY_STATUS.md) for a precise inventory. No unavailable step is silently replaced with synthetic data.

## Within-panel sensitivity analysis

The analysis performs 500 seeded replicates under three sampling schemes. All 250 families were detected in every stratified 80 percent replicate, while a median of 212 retained membership in at least five sampled genomes. This measures retention of already defined families. It is not equivalent to repeating discovery on independent panels.

## Independent-genome projection validation

The legacy projection matrix contains 3,192 genomes outside the 331-genome discovery panel. In 500 seeded, lineage-matched samples of 331 projection genomes, a median of 245 families were detected at least once (95 percent interval 242 to 248), and 206 occurred in at least five genomes (198 to 215). The accession overlap with the discovery panel was zero. This validates transfer of the released family HMMs to independent genomes; it does not repeat de novo family discovery.

## Held-out motif validation

For each of 12 tested families, genomes were split 70/30 for 30 seeded replicates. A position weight matrix learned only from training genomes was scored on held-out upstreams, dinucleotide-shuffled versions of those upstreams, and GC-matched upstreams from other families. All 12 families had a median AUC above 0.5 against dinucleotide shuffles. Ten of 12 exceeded an AUC of 0.5 against other-family upstreams in at least 26 of 30 resamples; nArcR-20 and nArcR-24 did not. The resamples overlap, so their distribution is reported descriptively rather than treated as 30 independent statistical observations. These results support predictive sequence signal for most tested families but do not establish protein binding or regulatory activity.

## Reciprocal exclusion rule

The portable implementation in `archaeal_dbd.thresholds` defines the reported rule explicitly. For a protein's best candidate interval, a canonical hit excludes the candidate when:

1. the aligned canonical interval overlaps more than the chosen fraction of the candidate interval; and
2. the canonical domain score is at least the candidate score minus the selected bit-score tolerance.

The manuscript configuration used an overlap fraction greater than 0.5 and a tolerance of 5 bits. The raw hit tables required for a full grid analysis are not in the current release.

## Tests and continuous integration

The test suite checks catalog size, membership counts, projection accounting, motif and structure summaries, legacy matrix dimensions, deterministic resampling, accession independence, motif helper behavior, and the reciprocal-exclusion definition. GitHub Actions runs these checks on every push and pull request.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. Database records and third-party derived data remain subject to their original terms.
