# Study design and claim hierarchy

## Confirmatory claims

1. Within Pfam clans CL0123 and CL0057, archaeal proteins outside a nineteen-family reference set form reproducible sequence families.
2. Most Pfam-uncharacterized families can be connected to an existing sequence or structural class using broader resources.
3. A small unresolved set remains after the tested sequence and structural comparisons.

## Exploratory claims

- Genomic neighborhoods nominate possible biological contexts.
- Upstream motifs nominate candidate binding specificities.
- Phylum clustering is compatible with vertical inheritance or ancient origin.
- Haloferax volcanii expression shows that mapped loci are transcribed.

These analyses do not establish direct transcriptional regulation, evolutionary mechanism, or condition-specific family function.

## Discovery panel rule

All genomes from Thaumarchaeota, DPANN, Asgardarchaeota, other deep TACK lineages, and Korarchaeota were retained. Euryarchaeota and Crenarchaeota were capped at 140 and 70. Within each capped group, subgroup allocations were proportional to subgroup abundance and sampling used seed 42. The resulting total of 331 was the output of this rule, not a prespecified target.

The two caps are pragmatic controls for database imbalance. They are not claimed to be biologically optimal.

## Required independent-panel validation

The present release adds an independent-genome projection test. Five hundred lineage-matched samples of 331 genomes were drawn from the legacy 3,192-genome projection matrix after excluding discovery-panel accessions. A median of 245 of 250 families were detected at least once (95 percent interval 242 to 248), and 206 occurred in at least five genomes (198 to 215). This tests whether the released family profiles transfer beyond the discovery accessions.

This analysis does not repeat candidate scanning, reciprocal canonical exclusion, clustering, or family definition. When the full 1,069-genome candidate pool becomes available, the complete discovery pipeline should still be repeated on multiple independent balanced panels. Families should be matched across runs by reciprocal profile searches and member overlap, with family recovery, splits, merges, boundary variation, and sensitivity to lineage caps and the minimum-genome threshold reported.

## Held-out motif validation

The revised analysis uses genome-level 70/30 train and test splits for each of 12 motif-bearing families. Position weight matrices are learned only from training genomes and evaluated on held-out family upstreams. Dinucleotide-shuffled versions of each held-out sequence control for composition and local sequence dependence, while GC-matched upstreams from other families test family specificity. Across 30 seeded splits, 10 of 12 motifs produced an AUC above 0.5 against the other-family control in at least 26 resamples. nArcR-20 and nArcR-24 did not. Because repeated splits overlap, the resampling distributions are reported descriptively and are not treated as independent replicates for hypothesis testing.

The additional genus-blocked analysis repeats this logic while holding out entire genera, with zero train/test genus overlap by construction. Seven of 12 families exceeded AUC 0.5 against other-family controls in at least 26 of 30 resamples. This asks whether an upstream signal transfers beyond the genera used to fit it. It remains a predictive sequence analysis and does not establish protein-DNA binding.

## Genome-quality sensitivity

The legacy accession-resolved projection was stratified by genome origin, assembly level, and reported completeness/contamination. Family prevalence and detection counts are reported for every stratum. The comparison is descriptive because these factors are confounded with lineage and with each other; it chiefly tests whether the catalog disappears when lower-quality assemblies are excluded.

This is predictive sequence validation, not a biochemical binding assay. The nArcR-6 model showed weak held-out discrimination from other families (median AUC 0.580) despite failing the earlier shuffled-family assignment null, so it remains a low-confidence motif.

## Required experimental validation

A focused experimental program should select a broadly distributed family, a family with strong dependence-corrected context, and an unresolved family. Direct binding, sequence specificity, and regulatory activity require assays such as EMSA or fluorescence anisotropy, footprinting or SELEX-seq, genetic perturbation, RNA-seq, and ideally ChIP-based localization.
