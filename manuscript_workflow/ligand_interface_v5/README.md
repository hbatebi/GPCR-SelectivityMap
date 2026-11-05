# Exploratory ligand-contact and mature Gα-interface case study

This workflow documents the **supporting, descriptive** ligand analysis added as a candidate Results subsection to the author's Nature Communications v5 manuscript.

Scope: four receptor examples (5-HT1A, D1, α2A, MC1R), and the five 5-HT1A–Gi1 Warren et al. (2024) structures (8FY8, 8FYT, 8FYE, 8FYL, 8FYX); ten dependent pairwise comparisons; four same-construct-category comparisons; contact cutoffs 4.0/4.5/5.0 Å; generic-numbered microswitch distances; α5 geometry; assay-evidence provenance.

**The five Warren structures and their pipeline outputs are contained in `../../source_data/natcomms_v5/Source_Data_29.zip`.** They are not duplicated in GitHub's inspectable workflow directory. Full RCSB-derived mmCIF files and large mapping inputs are included in that ZIP.

## One-command reproduction (from a cloned repository)

```bash
cd manuscript_workflow/ligand_interface_v5
bash RUN_FROM_SOURCE_DATA_29.sh
```

This expands Source Data 29 into a temporary work directory and reruns its original `RUN_ALL.sh`. Requires Python 3.10+ and internet access only if input coordinates are absent; use `SKIP_PIP_INSTALL=1` if Python dependencies are already installed. `RUN_ALL.sh` may create `GPCR_Ligand_Function_v5_COMPLETE_RESULTS.zip` in the extracted working directory. Run an offline smoke test with `RUN_OFFLINE_TEST.sh` from that directory if desired.

## Interpretation boundaries

This is a *case study* of receptor-side ligand-binding fingerprints versus mature Gα-contact fingerprints, not an experimental test of an activation probability or a cross-receptor rule. Numeric Gi1 Emax values for only three Warren ligands (93.4%, 97.4%, 102.8%) were verified; LSD and 5-MeO-DMT numerical Emax values are not imputed. Do not fit an efficacy–geometry correlation from the ten nonindependent pairwise comparisons.

PDB/mmCIF coordinate source: RCSB PDB. GtoPdb 2026.3.1 records are used for literature cross-reference. Source Data 29 contains underlying code and files with original input hashes.
