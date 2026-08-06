# GitHub release procedure

1. Create the public repository `hbatebi/GPCR-SelectivityMap` without an auto-generated README, license, or `.gitignore`.
2. Push the prepared repository and the annotated tag `v0.3.1`.
3. Create a GitHub release from tag `v0.3.1` using `RELEASE_NOTES_v0.3.1.md`.
4. Attach the wheel, source distribution, release ZIP, and SHA-256 checksum file from `release-assets/`.
5. Confirm that the GitHub Actions test matrix passes on Python 3.10, 3.11, and 3.12.
6. Enable the Zenodo GitHub integration, create a new release if needed, and add the minted DOI to `CITATION.cff` in the next patch release.

## Recommended repository settings

- Description: `Static GPCR profiling and phylogeny-aware benchmarking of endpoint-structure information`
- Website: leave blank until a documentation or DOI landing page is available
- Topics: `gpcr`, `structural-bioinformatics`, `phylogeny`, `machine-learning`, `protein-structure`, `g-protein`, `reproducible-research`
- Features: Issues enabled; Discussions optional; Wiki disabled unless maintained
