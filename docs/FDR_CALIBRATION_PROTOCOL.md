# Preregistered HMM error-rate calibration

Calibrate each family separately; a universal HMM cutoff is not assumed to have a universal error rate.

## Positive sets

- Family members withheld by genus, with no sequence from the held-out genus used to construct its profile.
- A second evaluation withholding entire major lineages when sample size permits.

## Negative sets

- Length- and composition-matched archaeal proteins outside the two target Pfam clans.
- Canonical archaeal DBD families excluded during discovery.
- Within-protein sequence reversals.
- Dipeptide-preserving shuffled sequences.
- Proteins drawn from held-out genera so taxonomic overlap cannot inflate specificity.

## Frozen analysis

For every family and threshold, report sensitivity, specificity, precision, empirical false-discovery proportion, AUROC, precision–recall AUC, and bootstrap intervals. Select a threshold only on training folds, targeting an empirical false-discovery proportion of at most 5%; evaluate it once on the untouched taxonomic holdout. Apply Benjamini–Hochberg correction to any family-level hypothesis tests.

Families without enough positives or negatives receive an `insufficient_data` result rather than a pooled threshold.

## Status

Not executed. Complete protein sequences, raw HMM hits, and the canonical negative collection are absent from the release. The supplied HMM archive and ten-family sequence subset are insufficient for catalog-wide calibration.
