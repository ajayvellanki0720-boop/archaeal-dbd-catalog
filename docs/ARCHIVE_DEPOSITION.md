# Immutable archive deposition

1. [x] Make the repository public only after checking that every released sequence-derived file may be redistributed under its source terms.
2. [x] Create a versioned GitHub release from the tested commit and attach the source archive plus `MANIFEST.sha256`.
3. [x] Connect the repository to Zenodo; use the stable concept DOI `10.5281/zenodo.23461089` in version-independent citations.
4. [x] Insert the concept DOI in `CITATION.cff`, the manuscript Data and Code Availability statement, and the cover letter.
5. [x] Confirm that the Zenodo draft lists both authors, the MIT software license, database-source caveats, and the manuscript title.
6. [x] Publish the prepared Zenodo record so that the DOI is registered.
7. [x] Download the published archive independently and verify `shasum -a 256 -c MANIFEST.sha256`.

Release metadata are recorded in `.zenodo.json`, `codemeta.json`, and `CITATION.cff`. The published archive SHA-256 is `aa78941b644c240d966f84c97628c3b30e094156daa66f13816d9649bafc905a`.

Version 0.5.0 adds the complete 13,237-sequence catalog archive and provenance tables. Its version-specific Zenodo record and archive checksum should be recorded here after publication; until then, the concept DOI remains the correct stable citation target.
