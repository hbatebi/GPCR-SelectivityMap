# One-command GPCR ligand–function analysis: results and editorial audit

**Run mode:** FULL. **Status:** complete 15-coordinate panel.
Four-receptor primary structural pairs: 4; sensitivity/control pairs: 3.
Within-Warren structural models processed: 5/5.
Within-Warren pairwise comparisons: 10/10; matching Gα construct group: 4.
Warren 2024 assay: 3 of 5 five-structure ligands have currently verified numeric Emax in the evidence input; the 2021 serotonin and aripiprazole structures have Emax measurements from the same 2024 functional assay, but they are not same-study specimens.

## Four primary receptor comparisons (4.5 Å heavy-atom contacts)
| Receptor | Ligand-contact Jaccard | Gα-contact Jaccard | α5-tip shift (Å) |
|---|---:|---:|---:|
| 5HT1A | 0.579 | 0.708 | 0.842 |
| DRD1 | 0.808 | 0.806 | 0.701 |
| ADRA2A | 0.923 | 0.789 | 1.075 |
| MC1R | 0.647 | 0.871 | 0.415 |

## Within-publication Warren five-ligand series
Pairwise comparisons use the same deposited structure multiple times; they are not independent experimental replicates.

| Pair | Ligand Jaccard | Gα Jaccard | α5 shift (Å) | Gα construct |
|---|---:|---:|---:|---|
| 5-MeO-DMT / 4-F,5-MeO-PyrT | 0.938 | 0.833 | 0.149 | mixed_Ga_construct |
| 5-MeO-DMT / vilazodone | 0.636 | 0.909 | 0.099 | mixed_Ga_construct |
| 5-MeO-DMT / LSD | 0.600 | 0.955 | 0.144 | mixed_Ga_construct |
| 5-MeO-DMT / buspirone | 0.520 | 0.955 | 0.478 | both_DNG |
| 4-F,5-MeO-PyrT / vilazodone | 0.667 | 0.826 | 0.174 | both_nonDNG |
| 4-F,5-MeO-PyrT / LSD | 0.632 | 0.870 | 0.272 | both_nonDNG |
| 4-F,5-MeO-PyrT / buspirone | 0.542 | 0.870 | 0.540 | mixed_Ga_construct |
| vilazodone / LSD | 0.565 | 0.952 | 0.141 | both_nonDNG |
| vilazodone / buspirone | 0.500 | 0.952 | 0.497 | mixed_Ga_construct |
| LSD / buspirone | 0.727 | 1.000 | 0.391 | mixed_Ga_construct |

## Prespecified sensitivity and evidence limitations
- Distance thresholds 4.0, 4.5, 5.0 Å and both Jaccard and overlap coefficients are reported in `outputs/contact_cutoff_and_ligand_size_sensitivity.tsv`.
- Primary α5-tip shifts also have an alternative receptor-core alignment sensitivity.
- All DNG-versus-non-DNG comparisons are identified. Do not interpret construct-dependent differences as pure ligand effects.
- Most ligand Gi1 Emax values cluster near 100%; there is too little range for a meaningful within-Warren regression.
- The 77.1% aripiprazole Gi1 Emax is assay-dependent; a different laboratory detected no quantifiable Gi1 signal in a different platform. Do not treat undetectable as zero.
- Original studies already described ligand recognition; the novelty here must concern *limitations of endpoint geometry* and must be qualified.
- No kinetic measurements, state probabilities, equivalence test or mechanistic ligand requirement are established by these atomic models.

## Editorial recommendation
A short **descriptive** subsection or supplementary case study can be justified after the full series is complete, provided the conclusion is phrased as examples of distinct orthosteric interactions coexisting with mature interfaces. The data do **not** establish why ligands are full/partial agonists, and do **not** demonstrate universal pocket-to-interface convergence. The v5 manuscript was not changed.

## Required author checks before submission
1. Inspect the newly downloaded deposited ligand conformations, full occupancy and construct records for 8FY8, 8FYT, 8FYE.
2. Verify numerical 5-MeO-DMT/LSD Emax and error from original source-data files if available; otherwise leave blank.
3. Cross-check all receptor–ligand/Gα author chain labels against PDB structures and papers.
4. Cite all original structural and functional mechanisms explicitly and avoid claiming newly discovered efficacy determinants.
5. If the construct-matched and cutoff sensitivity shows no consistent trend, describe cases individually, not as a general causal mechanism.
