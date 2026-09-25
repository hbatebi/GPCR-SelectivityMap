# GPCR SelectivityMap v0.4.0

## Purpose

This release accompanies the revised manuscript **“Mature GPCR complexes carry distributed information about the captured G protein class.”** It extends v0.3.1 with the analyses and audit trail added during peer review.

## Main scientific clarification

The release distinguishes two operational levels of a receptor-only reference:

1. **Coordinate-level receptor-only:** the deposited receptor coordinate model contains no transducer under the entity-based screen.
2. **Specimen-level experimentally partner free:** the experimental specimen itself contains no intracellular transducer or stabilising partner, and the primary publication does not describe the receptor as inactive.

Focused or local refinement can produce receptor-only coordinate models from transducer-containing cryo-EM specimens. Under the stricter specimen-level definition, only three receptors currently retain active matched references and Gq/11 is not represented. The three-class matched model is therefore not fitted under this stricter definition. The earlier coordinate-level analysis is preserved only as a transparent audit trail.

## Added in v0.4.0

- revision analysis suite and validation manifests;
- complete specimen-provenance audit with row-level deposited metadata, shared-particle evidence, corrected matched-set membership, feasibility summaries, input hashes and reproducible audit script;
- exact-UniProt Hauser et al. (2022) functional-repertoire annotation;
- G12/13 structural census;
- 92-feature descriptor dictionary;
- same-receptor 408-comparison table;
- follow-up four-fold class-coverage sensitivity outputs;
- revised manuscript figure assets and figure organization;
- complete Supporting Source Data S14 package aligned to SI Tables S22-S26.

## Claim boundary

Captured structural class is the G protein family physically present in a deposited mature complex. It is not treated as the receptor's unique or complete functional coupling repertoire. This package reports structural associations and validation results, not calibrated functional coupling probabilities.

## Release state

This archive is a **local release candidate**. Assign the public release DOI only after pushing/tagging v0.4.0 and creating the Zenodo-backed release.
