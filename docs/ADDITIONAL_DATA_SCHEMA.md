# Additional data required for full rigor

## Candidate genome pool

A tab-separated or CSV manifest with one row per candidate assembly and at least:

- assembly accession
- species taxid
- major group and subgroup
- assembly level
- RefSeq or GenBank source
- representative ranking fields
- selection status and exclusion reason

## Raw HMM hit tables

Both canonical and candidate-profile tables must include:

- target protein identifier
- genome accession
- query profile name and accession
- full-sequence score and E value
- domain score and conditional E value
- aligned coordinates `ali_from` and `ali_to`
- HMM coordinates
- profile gathering threshold
- HMMER and Pfam versions

These fields are sufficient for the implemented reciprocal-exclusion threshold grid.

## Corrected projection matrices

Presence-absence and copy-number matrices should contain one family row and one column for every genome in `genome_accounting.csv`, in the same order and with stable identifiers. Automated checks should require exactly 3,526 unique genome columns and 250 family rows.

## Independent-panel outputs

Every discovery replicate should retain:

- selected genome manifest and seed
- raw domain hits
- retained candidate proteins
- MMseqs cluster assignments
- filter decisions with explicit reason codes
- family HMMs
- reciprocal family matching to the reference catalog

