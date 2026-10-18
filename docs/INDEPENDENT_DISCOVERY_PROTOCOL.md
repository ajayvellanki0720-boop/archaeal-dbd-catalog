# Preregistered independent discovery benchmark

## Purpose

This benchmark tests whether the frozen discovery procedure finds coherent noncanonical DBD families in genuinely unseen archaeal genomes. Projection with the existing 250 family HMMs is not counted as independent discovery.

## Frozen inputs and cutoff

- Freeze repository version, HMM library, canonical reference set, overlap rule, score tolerance, clustering settings, and software containers before downloading the benchmark genomes.
- Record the exact NCBI/GTDB release dates and every assembly accession with version.
- Eligible genomes must have first public release dates after the discovery-panel cutoff and must not share an accession with either the 331-genome discovery panel or the 3,192-genome legacy projection.
- Deduplicate at the species level before analysis; retain one assembly according to a prespecified quality ranking.

## Primary outcomes

1. Number and fraction of genomes yielding at least one candidate after canonical-family exclusion.
2. Recovery of frozen catalog families at their frozen family-specific thresholds.
3. Number of candidate clusters not assigned to the frozen catalog.
4. Fraction of new clusters supported by at least five genomes and at least two genera.
5. Empirical false-positive rate against the decoys specified in `FDR_CALIBRATION_PROTOCOL.md`.

Report effect sizes and 95% cluster-bootstrap intervals. Do not tune thresholds on the benchmark set. Any post hoc threshold is a separately labelled exploratory analysis.

## Required deposited artifacts

Genome manifest; proteome files or stable URLs and checksums; raw canonical and candidate HMM hits; retained intervals; family assignments; cluster memberships; command log; software versions; random seeds; complete failure log; and all tables used in figures.

## Status

Not executed. The release does not contain the full original raw discovery inputs, and downloading a new time-split genome panel would create a new external dataset that must be frozen and reviewed before analysis.
