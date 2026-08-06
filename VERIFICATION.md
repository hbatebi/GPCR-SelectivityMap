# Verification

Release 0.3.1 was checked on 6 August 2026.

## Public-release checks

- All Python files in `src/`, `tests/`, and `manuscript/figure_scripts/` compiled successfully.
- Thirteen automated tests passed under the available Python 3.12 environment.
- The package installed in editable mode without downloading dependencies.
- A wheel, `gpcr_selectivitymap-0.3.1-py3-none-any.whl`, built successfully with local build tools.
- The installed package reported version `0.3.1`, and the command-line help completed successfully.
- All four main-figure scripts ran from the deposited source tables and regenerated SVG, PDF, and PNG outputs.
- Every included shell and Slurm script passed `bash -n` syntax validation.
- YAML and JSON release metadata were parsed, and the repository was scanned for common secret patterns and private key material; none were found.
- Personal HPC home paths and account identifiers were removed from the public workflow templates.

## Scientific checks retained from earlier project verification

The real-structure checks performed during development remain applicable because the underlying profiling implementation is unchanged:

- experimental beta2 adrenergic receptor complex PDB 3SN6, receptor chain R;
- contact, surface, approximate electrostatic, and elastic-network susceptibility blocks returned `ok`;
- TSV, JSON, SVG, and HTML profile outputs were generated.

The real structure and generic-number table used for that check are not redistributed in the GitHub repository.

The phylogeny-aware synthetic and manuscript-workflow checks retained in this release include:

- sequence-identity calculation and connected-component clustering;
- receptor-grouped and sequence-cluster benchmarking;
- full sequence, interface sequence, non-interface sequence, endpoint geometry, and global identity nearest-neighbour models;
- paired receptor-bootstrap output and interpretation manifests;
- matched non-interface and positional scrambling controls;
- label-provenance auditing;
- receptor-only entity definitions and static-feature extensions.

## Reproduction boundary

The compact source tables and figure-generation assets are included. Full reproduction additionally requires the external mmCIF archive, GPCRdb mappings, final analysis frames, complete out-of-fold predictions, and bootstrap-level outputs described in `docs/manuscript_reproduction.md`. Those larger inputs should be archived in a versioned research-data record such as Zenodo.
