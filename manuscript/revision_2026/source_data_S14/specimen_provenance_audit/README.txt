BJP specimen-level provenance source data
========================================

This directory contains the row-level evidence and reproducible script supporting SI Tables S22-S24 and the provenance classification used in Table S25.

Operational distinction used in the revised manuscript
------------------------------------------------------
- Coordinate-level receptor-only: the deposited receptor coordinate model contains no transducer under the entity-based screen.
- Specimen-level experimentally partner free: the experimental specimen itself contains no intracellular transducer or stabilising partner, and the primary publication does not describe the receptor as inactive.

Key files
---------
- build_provenance_audit.py: reproducible audit script (no model fitting).
- local_cif_metadata.tsv: deposited metadata assembled for the audit.
- audit_output/strict_51_provenance_audit.tsv: all 51 original strict receptor-only coordinate models.
- audit_output/terminal_P1_56_provenance_audit.tsv: terminal-fusion sensitivity inventory.
- audit_output/ligand_permissive_58_provenance_audit.tsv: ligand-permissive sensitivity inventory.
- audit_output/terminal_ligand_permissive_P1_59_provenance_audit.tsv: combined sensitivity inventory.
- audit_output/shared_particle_evidence.tsv: direct shared-particle evidence where available.
- audit_output/corrected_*: provenance-aware eligible structures, matched bound units and reference membership.
- audit_output/feasibility_summary.tsv/json: class-coverage summary under each definition.
- audit_output/INPUT_SHA256SUMS.txt: input hashes used by the audit.

The compact manuscript-facing exports for Tables S22-S26 are in ../provenance_compact/.
