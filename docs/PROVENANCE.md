# Provenance and computational environment

## Versioned software

The portable Python analyses use Python 3.12 with NumPy 2.0.2, pandas 2.3.3, SciPy 1.13.1, Matplotlib 3.9.4, and pytest 8.3.5. External-tool versions specified by `environment.yml` are HMMER 3.4, MMseqs2 15.6f452, MAFFT 7.526, trimAl 1.4.1, IQ-TREE 2.3.6, and Foldseek 9.427df8a. `requirements.lock.txt`, `environment.yml`, and `Dockerfile` provide three installation routes.

## Database provenance

The archived analysis identifies Pfam clans CL0123 and CL0057, CDD version 3.21, AlphaFold Database structures, NCBI RefSeq/GenBank assemblies, GTDB-based lineage groupings, and GEO series GSE293312. Exact Pfam/GTDB releases, database download dates, per-record AlphaFold versions, and the full genome-download manifest were not preserved in the supplied archive. They must not be inferred retrospectively. Their absence is recorded as a release limitation and the required fields are specified in `ADDITIONAL_DATA_SCHEMA.md`.

## Catalog sequence recovery

On 19 October 2026, `scripts/recover_catalog_sequences.py` used the assembly FTP paths recorded in the deposited genome metadata, repaired archived paths to their current versioned NCBI assembly locations, and extracted the accession named by every membership row. The assembly protein files yielded 13,232 of 13,237 sequences. Five WP accessions absent from the current assembly files were retrieved individually through NCBI Protein E-utilities. `release/tables/catalog_sequence_retrieval.csv` records each source URL, SHA-256 checksum, byte count, date, target count, recovered count, and fallback count. `catalog_sequence_fallback.csv` records the five fallback accessions and URLs. The resulting deterministic gzip archive contains one sequence for every membership row, and automated tests verify identifiers and lengths.

## Seeds

- Original panel sampling and legacy permutation analyses: 42
- Within-panel sensitivity: 10,000 plus replicate index
- Independent projection sensitivity: 20,000 plus replicate index
- Genome-blocked motif validation: 30,000 plus family and replicate offsets
- Genus-blocked motif validation: 40,000 plus family and replicate offsets

Replicate-level output files retain the replicate identifier. The source code defines the exact offset rule.

## Commands

```bash
python -m pip install -e '.[test]'
make smoke
make test
make all
python scripts/recover_catalog_sequences.py
```

The first two commands provide a rapid installation and integrity check. `make all` regenerates every analysis supported by the deposited inputs. It does not reconstruct missing raw discovery or corrected 3,526-genome matrices.

## Integrity

`MANIFEST.sha256` records a SHA-256 digest for every deposited file except the manifest itself and transient version-control/cache directories. Regenerate it only after all release files are finalized:

```bash
python scripts/build_manifest.py
shasum -a 256 -c MANIFEST.sha256
```
