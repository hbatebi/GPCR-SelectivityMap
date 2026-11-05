# GPCR SelectivityMap v0.5.2 — Nature Communications v5 source data and exploratory ligand case study

**Release candidate, 2026-11-05. Confirm scientific/metadata sign-off before publishing.**

This release adds a complete Nature Communications v5 reproducibility archive with Source Data 1–29 and an exploratory ligand-pocket versus mature Gα-interface case study. It preserves earlier receptor-based benchmark workflows and the existing BJP revision records, without changes to the core software APIs.

## New materials

- `source_data/natcomms_v5/`: 29 source-data archives, original filenames, `Supplementary_Figure_12.tif`, per-file SHA-256 manifest in `docs/`.
- `manuscript_workflow/ligand_interface_v5/`: inspectable Python scripts, lightweight inputs, sensitivity tables and figures; `bash RUN_FROM_SOURCE_DATA_29.sh` reruns full inputs from Source Data 29.
- `docs/NATCOMMS_V5_LIGAND_DATA_RELEASE.md`: scientific scope, source provenance and limits.

## Scope / limitations

Five 5-HT1A–Gi1 structures are compared in ten non-independent ligand pairs. Four same-Gα-construct-category comparisons retain a contrast between ligand-contact and mature Gα-contact fingerprints at tested cutoffs. **This is descriptive and cannot demonstrate a necessary ligand property, activation kinetics, or general convergence of mature GPCR endpoints.** Only three of five within-study ligands have verified numerical Gi1 Emax annotations.

## Versioning

Zenodo concept DOI: `10.5281/zenodo.21820928`. A *new version-specific DOI* should be minted by Zenodo after this GitHub release is published; the previous record is not overwritten. Preserve the Zenodo integration and confirm the new release is listed under the same concept record. No version DOI is preassigned here.
