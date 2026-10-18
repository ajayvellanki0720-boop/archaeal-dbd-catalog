# Minimal worked example

This synthetic example exercises the reciprocal canonical-domain exclusion rule without downloading a database. The candidate domain overlaps a canonical hit and is excluded with a 5-bit tolerance but retained with a 2-bit tolerance.

Run from the repository root:

```bash
python -m archaeal_dbd.smoke
```

The command compares the computed rows with `expected_retention.csv` and exits nonzero if they differ.
