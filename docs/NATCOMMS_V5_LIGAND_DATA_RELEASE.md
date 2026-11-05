# Nature Communications v5: full source-data release + exploratory ligand-interface case study

**Data cut:** 2026-11-05. **Status:** manuscript candidate, pending author's scientific and copyright sign-off.

## What this adds to GPCR-SelectivityMap

1. Full matching set of **Source Data 1–29** for the Nature Communications v5 candidate, in `source_data/natcomms_v5/`, with exact source filenames and SHA-256 manifest.
2. `Source_Data_29.zip`: underlying analysis code, complete five-ligand 5-HT1A–Gi1 coordinate series and four-receptor comparator panel, pharmacology evidence, computed contacts, microswitches, Gα-interface geometry, controls, source hashes, and scientific caveats.
3. An inspectable workflow (`manuscript_workflow/ligand_interface_v5/`) with scripts, lightweight input tables, results, and figures. `RUN_FROM_SOURCE_DATA_29.sh` re-extracts full inputs and executes the archived pipeline.
4. The companion image `source_data/natcomms_v5/Supplementary_Figure_12.tif`.

## Claim boundary

Four same-Gα-construct-category comparisons in one five-structure 5-HT1A case series exhibit higher receptor–Gα contact-set overlap than receptor–ligand contact overlap across the tested contact distance cutoffs. This is an **exploratory, descriptive observation**. Five experimental structures and four dependent pairwise comparisons cannot establish the existence of a ligand property necessary for productive coupling, a cross-GPCR principle, true endpoint equivalence, or ligand-efficacy causality. The original Warren et al. (2024) study already established ligand-dependent pharmacology and the binding-site context.

Only three of the five Warren ligands have verified quantitative same-assay Gi1 Emax values in the included evidence input, spanning 93.4–102.8% relative to serotonin. Do not impute missing measurements or treat ten pairwise comparisons as independent replicates. Different DNG constructs, ligand sizes and stabilizing conditions are explicitly flagged.

## Provenance and responsibilities

- Original coordinates: RCSB PDB; preserve accession IDs, original authorship and data citations.
- Structure catalog: GPCRdb 2026 export and GPCR-SelectivityMap v5 source data, where applicable.
- Quantitative pharmacology: Warren et al., *Nature* (2024), DOI 10.1038/s41586-024-07403-2, qualified evidence table supplied in Source Data 29.
- The `Source_Data_1–28` ZIPs are preserved **byte-for-byte** from the integrated manuscript candidate; `Source_Data_29.zip` is the final case-study archive, not the separately labelled, **provisional** archive-wide ligand analysis.
- This deposit deliberately omits unpublished cover letters, editable manuscript drafts, internal notes and acceptance-related language.

## Release/DOI

Repository: https://github.com/hbatebi/GPCR-SelectivityMap

Previous Zenodo concept DOI: https://doi.org/10.5281/zenodo.21820928. The version-specific DOI for this release does **not** exist until the corresponding GitHub Release has been published and ingested by Zenodo. Do not cite a guessed DOI.

Release target: `v0.5.2` (new data/workflow archive; no change to the existing GPCR-SelectivityMap package APIs). Existing code and `source_data/bjp_revision_2026/` must be kept intact.
