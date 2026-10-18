# Immutable archive deposition

1. Make the repository public only after checking that every released sequence-derived file may be redistributed under its source terms.
2. Create a versioned GitHub release from the tested commit and attach the source archive plus `MANIFEST.sha256`.
3. Connect the repository to Zenodo, archive the release, and reserve or mint the DOI.
4. Replace the DOI placeholder in `CITATION.cff`, the manuscript Data and Code Availability statement, and the cover letter.
5. Confirm that the Zenodo record lists both authors, the MIT software license, database-source caveats, and the manuscript title.
6. Download the archive from Zenodo independently and verify `shasum -a 256 -c MANIFEST.sha256`.

No DOI is included in the current package because a DOI cannot be truthfully assigned before deposit.
