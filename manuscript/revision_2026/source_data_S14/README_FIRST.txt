BJP Supporting Source Data S14 — complete revision package
=============================================================

Purpose
-------
This package contains the peer-review revision source data supporting SI Tables S22-S26, the specimen-level provenance audit, the finalized functional-repertoire annotation, descriptor definitions, same-receptor controls and related sensitivity analyses.

Key manuscript mappings
-----------------------
- SI Table S22: provenance_compact/Table_S22_specimen_provenance_audit.tsv; full row-level evidence in specimen_provenance_audit/audit_output/*.
- SI Table S23: provenance_compact/Table_S23_matched_structure_ids.tsv; corrected membership tables in specimen_provenance_audit/audit_output/corrected_*_reference_membership.tsv.
- SI Table S24: provenance_compact/Table_S24_matched_feasibility.tsv; source feasibility_summary.tsv/json in specimen_provenance_audit/audit_output/.
- SI Table S25: provenance_compact/Table_S25_provenance_stratified_margins.tsv; underlying classifier/margin outputs in expanded_matched_audit/ and followup_fourfold/.
- SI Table S26: provenance_compact/Table_S26_functional_repertoire_context.tsv; finalized Hauser et al. (2022) mapping in functional_repertoire/ and G12/13 census in g12_g13/.

Additional included material
----------------------------
- matched_92_feature_dictionary.tsv and the generating descriptor script;
- all 408 same-receptor feature-by-class-pair comparisons requested by Reviewer 1;
- finalized exact-human-UniProt Hauser et al. (2022) common-coupling-map annotation;
- G12/13 structural census;
- post-full-run four-fold class-coverage technical sensitivity;
- relevant revision scripts, prespecification notes and manifests.

Provenance terminology
----------------------
The revised manuscript distinguishes coordinate-level receptor-only entries from specimen-level experimentally partner-free structures. The specimen-level audit is supplied here in full so that the classification can be reproduced from deposited metadata and the recorded evidence.
