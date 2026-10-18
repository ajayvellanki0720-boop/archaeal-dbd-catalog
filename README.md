# Noncanonical archaeal DNA binding domain catalog

This repository accompanies the manuscript **A Catalog of Noncanonical HTH and RHH DNA Binding Domain Families in Archaea**. It organizes the released family catalog, records its evidentiary limits, and provides portable checks and a reproducible within-panel robustness analysis.

## Scope

The study is a bounded survey of proteins detected through Pfam clans CL0123 and CL0057 after removing a nineteen-family reference set. It is **not** a comprehensive census of all archaeal DNA-binding proteins, and a family assignment is **not** proof that a protein is a transcription factor.

## Key findings represented in this release

- 38,153 proteins in the 331-genome panel carried a targeted clan domain.
- 16,535 were assigned to the reference set, 14,451 to another characterized family, and 7,167 were uncharacterized.
- Clustering and filtering defined 250 candidate families containing 13,237 proteins.
- Held-out family assignment recovered 97.3 percent of withheld members.
- Of 79 Pfam-uncharacterized families, 56 were resolved by CDD and 14 more by structural comparison; nine remained unresolved by both resources tested.
- Eleven of twelve tested upstream motifs passed all three sequence null models.
- Genomic-neighborhood results are functional hypotheses, not confirmed regulons.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
make test
make all
```

Generated audit and resampling outputs are written to `results/generated/`.

## Repository layout

```text
src/archaeal_dbd/     Portable validation and sensitivity code
tests/                 Automated checks for headline invariants
release/data/          Released matrices and analysis archives
release/tables/        Released result tables
release/figures/       Released figures
legacy/code/           Original scripts preserved without claiming portability
results/               New within-panel resampling result
manuscript/            Revised manuscript
docs/                  Study design and reproducibility limitations
```

## Reproducibility status

`make all` reproduces the release audit and the new within-panel sensitivity analysis. It does not rerun raw domain discovery because the released archive lacks the complete 1,069-genome candidate pool and raw canonical and candidate HMM hit tables. The 3,526-genome corrected prevalence summaries are present, but the raw presence and copy-number matrices cover the earlier 3,192-genome pass.

See [docs/REPRODUCIBILITY_STATUS.md](docs/REPRODUCIBILITY_STATUS.md) for a precise inventory. No unavailable step is silently replaced with synthetic data.

## Within-panel sensitivity analysis

The analysis performs 500 seeded replicates under three sampling schemes. All 250 families were detected in every stratified 80 percent replicate, while a median of 212 retained membership in at least five sampled genomes. This measures retention of already defined families. It is not equivalent to repeating discovery on independent panels.

## Reciprocal exclusion rule

The portable implementation in `archaeal_dbd.thresholds` defines the reported rule explicitly. For a protein's best candidate interval, a canonical hit excludes the candidate when:

1. the aligned canonical interval overlaps more than the chosen fraction of the candidate interval; and
2. the canonical domain score is at least the candidate score minus the selected bit-score tolerance.

The manuscript configuration used an overlap fraction greater than 0.5 and a tolerance of 5 bits. The raw hit tables required for a full grid analysis are not in the current release.

## Tests and continuous integration

The test suite checks catalog size, membership counts, projection accounting, motif and structure summaries, legacy matrix dimensions, deterministic resampling, and the reciprocal-exclusion definition. GitHub Actions runs these checks on every push and pull request.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. Database records and third-party derived data remain subject to their original terms.

