# GitHub release procedure for v0.4.0

1. Push this release candidate to `hbatebi/GPCR-SelectivityMap` on a review branch.
2. Confirm the README, `CITATION.cff`, `CHANGELOG.md`, `RELEASE_NOTES_v0.4.0.md`, and revision source tables.
3. Let GitHub Actions run the test matrix.
4. Merge to the default branch after the checks pass.
5. Create the annotated tag `v0.4.0`.
6. Create a GitHub Release from `v0.4.0` using `RELEASE_NOTES_v0.4.0.md`.
7. If you distribute built wheel/source archives, build them in a network-enabled environment and attach them to the release.
8. Confirm the Zenodo-GitHub integration archives the release.
9. Add the minted version-specific Zenodo DOI to `CITATION.cff` and to the manuscript Data Availability statement. Keep `10.5281/zenodo.21820928` as the concept DOI.
10. Deposit any provenance-sensitive or large Supporting Source Data files in the same versioned Zenodo record if they are not appropriate for GitHub.

## Recommended repository settings

- Description: `Static GPCR profiling and phylogeny-aware benchmarking of captured-class information`
- Topics: `gpcr`, `structural-bioinformatics`, `phylogeny`, `machine-learning`, `protein-structure`, `g-protein`, `reproducible-research`
- Issues enabled; Discussions optional.
