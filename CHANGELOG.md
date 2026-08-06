# Changelog

## 0.3.1 - 2026-08-06

- Added the manuscript software workflow as editable SVG, PDF, and high-resolution PNG.
- Repositioned the public release around conserved activation, evolutionary organization, receptor-side compatibility, partner-conditioned accommodation, and receptor-only transfer.
- Moved interface-complementarity and receptor-compatibility commands to an explicitly experimental section.
- Clarified that the software separates association, evolutionary organisation, and transfer, and does not report calibrated functional coupling probabilities.
- Updated release metadata, manuscript title, figure assets, and public documentation without changing numerical analysis code.

## 0.3.0

- Added receptor–Gα physical contact extraction.
- Added Kyte–Doolittle hydropathy, charge, side-chain-volume, aromatic, and combined interface-complementarity scores.
- Added complex cognate-versus-decoy chemistry audit.
- Added training-only receptor-alone Gα-family templates.
- Added sequence-only and structure-coverage receptor modes.
- Added receptor-grouped and low-homology validation with identity-kNN baselines.
- Added receptor-position and contact-map/Gα-chemistry scrambling controls.
- Added ligand-permissive and polymer-only receptor-only transfer workflows.
- Added resumable Fritz quick and full jobs.

## 0.2.0 - 2026-08-05

Manuscript-aligned, phylogeny-aware release.

### Added

- `cluster-benchmark` for receptor-grouped and 30%, 40%, and 50% sequence-cluster transfer
- global sequence-identity kNN baseline
- coarse taxonomy and fine receptor-family baselines
- predefined intracellular-interface sequence model
- `interface-audit` with non-interface, matched-position, and scrambling controls
- `provenance-audit` for direct, consensus, curated, and homology-inferred labels
- `aggregate` for receptor and receptor-by-transducer feature aggregation
- interpretation gate separating receptor-grouped association from low-homology transfer
- manuscript sequence-phylogeny audit figures and source tables
- exact Fritz sequence-interface audit workflow

### Changed

- repositioned the package from coupling prediction to representation profiling and benchmarking
- revised documentation to state that receptor-grouped performance can reflect evolutionary proximity
- demoted all probability language that could be mistaken for validated functional coupling prediction

### Scientific interpretation represented by this release

Sequence and interface representations can be predictive within the sampled evolutionary neighbourhood, but their advantage weakens under stringent homology removal and is matched by global sequence proximity. Endpoint geometry can retain more information at intermediate divergence, while neither modality establishes a universal receptor-intrinsic coupling code at the strictest threshold.

## 0.1.0 - 2026-08-05

Initial GitHub-ready release with static profiling, batch extraction, grouped benchmarking, documentation, and manuscript source data.
