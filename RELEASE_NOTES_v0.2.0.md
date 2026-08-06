# GPCR SelectivityMap v0.2.0

This release aligns the software with the phylogeny-revised manuscript.

## Scientific change

The package no longer frames receptor-grouped classification as evidence for a transferable coupling code. It now explicitly tests whether sequence and static structure generalize beyond close homologues and whether predefined interface residues outperform global sequence proximity and matched non-interface positions.

## Main additions

- low-homology cluster holdouts at user-defined identity thresholds;
- global sequence-identity kNN baseline;
- taxonomy-only baselines;
- predefined intracellular-interface sequence model;
- matched non-interface and positional scrambling controls;
- receptor and receptor-by-transducer aggregation;
- functional-label provenance audit;
- automated, claim-bounded interpretation gate;
- exact Fritz manuscript audit workflow and compact completed outputs.

## Current manuscript-aligned interpretation

The deposited audit supports a phylogeny-dominated sequence signal rather than a universal interface code. Endpoint geometry retains more information at intermediate sequence divergence, while neither sequence nor static endpoint geometry demonstrates reliable transfer at the strictest homology threshold. This is absence of evidence for universal transfer, not proof that interface chemistry is irrelevant.

## Deposit recommendation

Upload the repository to GitHub as release `v0.2.0`. Attach the wheel and checksum to the release. Archive the same release in Zenodo and place the Zenodo DOI in `CITATION.cff` after minting. Do not commit the large structure collection, full checkpoints, or full bootstrap archives to GitHub.
