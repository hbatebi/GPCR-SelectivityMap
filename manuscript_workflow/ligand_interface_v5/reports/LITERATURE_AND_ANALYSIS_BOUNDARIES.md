# Primary sources and scientific boundaries

**Study baseline:** Batebi Nature Communications manuscript v5 (uploaded by author). The original manuscript is unchanged. Embedded `v5_inputs/Source_Data_1_from_v5.zip` is an input copy, not a new publication source.

## Primary references
- 5-HT1A serotonin/aripiprazole complexes: Nature (2021), DOI 10.1038/s41586-021-03376-8, PDB 7E2Y, 7E2Z.
- Warren et al. (2024), *Structural pharmacology and therapeutic potential of 5-methoxytryptamines*, Nature 630:237–246, DOI 10.1038/s41586-024-07403-2. PDB 8FY8, 8FYT, 8FYE, 8FYL, 8FYX. The paper reports 5-HT1A Gi1 BRET Emax of 93.4% for buspirone, 97.4% for vilazodone, 102.8% for 4-F,5-MeO-PyrT and 77.1% for aripiprazole (5-HT normalized to 100%). LSD and 5-MeO-DMT described as full/near-full, not numerically transcribed here.
- D1 structures: DOI 10.1038/s41467-022-30929-w, PDB 7X2C, 7X2D, and a PAM-confounded control 7X2F.
- Alpha2A structures: DOI 10.1038/s12276-024-01296-x, PDB 9CBL, 9CBM.
- MC1R structures: DOI 10.1038/s41422-021-00557-y, PDB 7F4D, 7F4H, 7F4I.

## Confirmed interpretation boundaries
- Deposited mature complexes are conditionally selected experimental structures, not measurements of a ligand's probability of reaching the complex.
- No direct kinetic conclusions can be inferred from contact overlaps or Gα α5 pose.
- Five Warren structures have similar high Gi1 efficacies; numerical efficacy exists for only three of the five in this package. The serotonin and aripiprazole Emax values are from the same 2024 functional assay but their coordinate structures were published in 2021.
- Buspirone and 5-MeO-DMT use dominant-negative Gαi1 constructs in Warren et al. (2024), while the other three are labelled non-DNG. Matched-construct comparisons are reported separately.
- Ligand-contact Jaccard is not independent of ligand size, resolution or structure preparation. Sensitivity analyses use multiple distance thresholds, overlap coefficients and alignment windows.
- Biological independent units are receptor/experimental ligand perturbations, not pairwise combinations of five structures. Ten pairwise comparisons are not ten independent experimental replicates.
- Quantitative full vs partial labels from heterogeneous GtoPdb assays are not substituted for same-assay Emax values.
- Previously published binding mechanisms must be cited; this re-analysis may supply a cross-system comparison, not claim to rediscover them.
- Stronger conclusions require broader efficacy range, matched constructs and direct G protein activation assays or perturbations.
