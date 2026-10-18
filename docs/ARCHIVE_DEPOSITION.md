# Immutable archive deposition

1. [x] Make the repository public only after checking that every released sequence-derived file may be redistributed under its source terms.
2. [x] Create a versioned GitHub release from the tested commit and attach the source archive plus `MANIFEST.sha256`.
3. [x] Connect the repository to Zenodo and reserve DOI `10.5281/zenodo.23461090` for version 0.4.0.
4. [x] Insert the reserved DOI in `CITATION.cff`, the manuscript Data and Code Availability statement, and the cover letter.
5. [x] Confirm that the Zenodo draft lists both authors, the MIT software license, database-source caveats, and the manuscript title.
6. [ ] Publish the prepared Zenodo record so that the DOI is registered.
7. [ ] Download the published archive independently and verify `shasum -a 256 -c MANIFEST.sha256`.

Release metadata are pre-populated in `.zenodo.json`, `codemeta.json`, and `CITATION.cff`. The DOI is currently reserved; it becomes registered when the Zenodo draft is published.
