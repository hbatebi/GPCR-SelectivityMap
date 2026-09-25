# BJP revision materials (v0.4.0)

This directory contains the compact, auditable analysis and figure material added during the 2026 BJP peer-review revision.

- `analysis/revision_bjp_2026/`: collected revision run, manifests, checksums, outputs and run specification.
- `functional_repertoire/A05_Hauser2022_finalization/`: finalized Hauser et al. (2022) common-coupling-map mapping by exact human UniProt accession.
- `followup_fourfold/`: post-full-run four-fold class-coverage sensitivity files.
- `source_tables/`: high-value reviewer-facing source tables, including the 92-feature dictionary, 408-comparison table, G12/13 census and compact SI Tables S22-S26.
- `figures/`: revised five-figure manuscript assets.

Important distinction: `coordinate-level receptor-only` means no transducer in the deposited receptor coordinate model under the entity screen. `experimentally partner free` additionally requires a transducer-free/stabilizer-free experimental specimen under the specimen-level provenance audit.

### Complete Supporting Source Data S14

`source_data_S14/` contains the complete manuscript-facing provenance and reviewer-revision source package, including the row-level specimen audit and reproducible `build_provenance_audit.py` script, compact exports for SI Tables S22-S26, finalized Hauser et al. mapping, G12/13 census, descriptor dictionary, same-receptor 408-comparison table and related sensitivity outputs.
