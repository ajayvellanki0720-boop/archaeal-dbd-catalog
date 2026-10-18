# Preregistered discovery-parameter multiverse

The primary catalog must remain frozen. The multiverse evaluates whether its conclusions depend on defensible alternative settings.

| Decision | Frozen value | Sensitivity grid |
|---|---:|---|
| Candidate/canonical overlap fraction | >0.50 | >0.30, >0.40, >0.50, >0.60, >0.70 |
| Canonical score tolerance | 5 bits | 0, 2.5, 5, 10 bits |
| Clustering identity | manuscript setting | 20%, 25%, 30%, 35%, 40% |
| Clustering coverage | manuscript setting | 50%, 60%, 70%, 80% |
| Minimum family support | 5 genomes | 3, 5, 10 genomes |
| Upstream window | 250 bp | 100, 150, 200, 250, 300 bp |

Primary robustness outcomes are the adjusted Rand index of family assignments, fraction of the 250 families recovered, changes in the unresolved-family count, changes in lineage breadth, and changes in motif conclusions. Report the entire grid, not only settings that agree with the manuscript.

The released panel resampling, independent projection, genome-quality stratification, genus-blocked motif validation, and leave-one-lineage-out motif validation cover sampling and validation choices. They do not replace the raw discovery multiverse above.

## Status

Partially executed. Sampling and validation dimensions are versioned. Discovery-threshold and clustering dimensions are blocked by the absent candidate pool, raw HMM hits, and complete member sequences.
