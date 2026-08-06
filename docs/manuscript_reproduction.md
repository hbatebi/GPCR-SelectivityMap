# Manuscript reproduction

The repository preserves two project-specific overlays.

## Static receptor-only extension

`manuscript_workflow/gpcr_icl2_overlay/` contains the exact feature extraction and Fritz launchers for contact networks, surface and electrostatic descriptors, and elastic-network susceptibility.

## Sequence-interface and phylogeny audit

`manuscript_workflow/sequence_interface_audit_v1/` contains the exact Fritz driver and Slurm launchers used to test:

- full sequence;
- predefined interface sequence;
- non-interface sequence;
- taxonomy baselines;
- global sequence-identity kNN;
- endpoint geometry;
- sequence plus geometry;
- matched non-interface and scrambling controls;
- granularity and label-provenance controls.

Compact completed outputs are under `manuscript/sequence_phylogeny_audit/`.

## External inputs required

Full reproduction additionally requires:

- the experimental mmCIF structure collection;
- GPCRdb generic-number mapping;
- the bound receptor-transducer representative table;
- the strict receptor-only audit;
- functional annotation tables with provenance;
- the final analysis configuration;
- the complete out-of-fold prediction and bootstrap archives.

These large or provenance-sensitive inputs should be retrieved from the manuscript's versioned Zenodo record.
