#!/usr/bin/env python3
"""Produce pretransparent scientific digest and v5 conditional draft, no invented statistics."""
from pathlib import Path
import argparse, pandas as pd, numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs';REPORT=ROOT/'reports';F=ROOT/'figures'
a=argparse.ArgumentParser();a.add_argument('--mode',choices=['full','smoke'],required=True);args=a.parse_args()
q=pd.read_csv(OUT/'four_receptor_atomic/corrected_pair_comparisons.tsv',sep='\t');f=pd.read_csv(OUT/'atomic_five/corrected_pair_comparisons.tsv',sep='\t')
w=pd.read_csv(OUT/'five_ligand_structural_readouts.tsv',sep='\t')
cf=pd.read_csv(OUT/'contact_cutoff_and_ligand_size_sensitivity.tsv',sep='\t')
full=args.mode=='full'
if full and (len(f)!=10 or w['TM3.50_TM6.34_CA'].notna().sum()!=5):
    raise SystemExit('STOP: Full analysis requested, but five-ligand results are incomplete')
base=q[q.role=='primary'].copy();base.to_csv(OUT/'four_primary_case_summary.tsv',sep='\t',index=False)
ctrl=q[q.role!='primary'].copy();ctrl.to_csv(OUT/'four_receptor_controls.tsv',sep='\t',index=False)
# exact-interval sensitivity provided as per-pair, not an inferred population effect
all_pairs=pd.concat([base.assign(source='four_receptor_panel'),f.assign(source='Warren2024_series')],ignore_index=True)
all_pairs.to_csv(OUT/'primary_plus_Warren_pairwise_descriptive.tsv',sep='\t',index=False)
assay=pd.read_csv(ROOT/'input/5HT1A_functional_evidence.tsv',sep='\t')
assay.to_csv(OUT/'functional_assay_evidence_with_source.tsv',sep='\t',index=False)
match=f[f.role.ne('mixed_Ga_construct')]
numeric=assay[assay.Emax_percent_of_serotonin.notna()]
cut_base=cf[cf.cutoff_A.eq(4.5)]
text=['# One-command GPCR ligand–function analysis: results and editorial audit','',
      f'**Run mode:** {args.mode.upper()}. **Status:** {"complete 15-coordinate panel" if full else "offline smoke / partial five-structure series"}.',
      f'Four-receptor primary structural pairs: {len(base)}; sensitivity/control pairs: {len(ctrl)}.',
      f'Within-Warren structural models processed: {w.TM3_50_TM6_34_CA.notna().sum() if "TM3_50_TM6_34_CA" in w else int(w["TM3.50_TM6.34_CA"].notna().sum())}/5.',
      f'Within-Warren pairwise comparisons: {len(f)}/10; matching Gα construct group: {len(match)}.',
      f'Warren 2024 assay: {int(w.Gi1_BRET_Emax_percent_5HT.notna().sum())} of 5 five-structure ligands have currently verified numeric Emax in the evidence input; the 2021 serotonin and aripiprazole structures have Emax measurements from the same 2024 functional assay, but they are not same-study specimens.',
      '', '## Four primary receptor comparisons (4.5 Å heavy-atom contacts)',
      '| Receptor | Ligand-contact Jaccard | Gα-contact Jaccard | α5-tip shift (Å) |',
      '|---|---:|---:|---:|']
for r in base.itertuples(index=False):text.append(f'| {r.case_id} | {r.ligand_contact_jaccard:.3f} | {r.ga_contact_jaccard:.3f} | {r.alpha5_tip_shift_A:.3f} |')
text+=['','## Within-publication Warren five-ligand series',
       'Pairwise comparisons use the same deposited structure multiple times; they are not independent experimental replicates.']
if not f.empty:
    text+=['', '| Pair | Ligand Jaccard | Gα Jaccard | α5 shift (Å) | Gα construct |', '|---|---:|---:|---:|---|']
    for r in f.itertuples(index=False):
        text.append(f'| {r.first_ligand} / {r.second_ligand} | {r.ligand_contact_jaccard:.3f} | {r.ga_contact_jaccard:.3f} | {r.alpha5_tip_shift_A:.3f} | {r.role} |')
text+=['', '## Prespecified sensitivity and evidence limitations',
       '- Distance thresholds 4.0, 4.5, 5.0 Å and both Jaccard and overlap coefficients are reported in `outputs/contact_cutoff_and_ligand_size_sensitivity.tsv`.',
       '- Primary α5-tip shifts also have an alternative receptor-core alignment sensitivity.',
       '- All DNG-versus-non-DNG comparisons are identified. Do not interpret construct-dependent differences as pure ligand effects.',
       '- Most ligand Gi1 Emax values cluster near 100%; there is too little range for a meaningful within-Warren regression.',
       '- The 77.1% aripiprazole Gi1 Emax is assay-dependent; a different laboratory detected no quantifiable Gi1 signal in a different platform. Do not treat undetectable as zero.',
       '- Original studies already described ligand recognition; the novelty here must concern *limitations of endpoint geometry* and must be qualified.',
       '- No kinetic measurements, state probabilities, equivalence test or mechanistic ligand requirement are established by these atomic models.',
       '', '## Editorial recommendation',
       'A short **descriptive** subsection or supplementary case study can be justified after the full series is complete, provided the conclusion is phrased as examples of distinct orthosteric interactions coexisting with mature interfaces. The data do **not** establish why ligands are full/partial agonists, and do **not** demonstrate universal pocket-to-interface convergence. The v5 manuscript was not changed.',
       '', '## Required author checks before submission',
       '1. Inspect the newly downloaded deposited ligand conformations, full occupancy and construct records for 8FY8, 8FYT, 8FYE.',
       '2. Verify numerical 5-MeO-DMT/LSD Emax and error from original source-data files if available; otherwise leave blank.',
       '3. Cross-check all receptor–ligand/Gα author chain labels against PDB structures and papers.',
       '4. Cite all original structural and functional mechanisms explicitly and avoid claiming newly discovered efficacy determinants.',
       '5. If the construct-matched and cutoff sensitivity shows no consistent trend, describe cases individually, not as a general causal mechanism.',
       '']
(OUT/'INTEGRATED_SCIENTIFIC_REPORT.md').write_text('\n'.join(text))
# Manuscript-ready wording with cautious placeholders when the 5 are not complete
paragraph=['# Candidate short Results subsection for Nature Communications v5','',
           '**Proposed heading:** Mature G-protein interfaces incompletely capture ligand-dependent pharmacology', '',
           '**STATUS:** '+('DESCRIPTIVE DRAFT AFTER COORDINATE ANALYSIS; REQUIRES AUTHOR REVIEW.' if full else 'NOT FOR MANUSCRIPT: full series not yet run.'),'',
           'To examine whether mature complex geometry distinguishes ligand-associated signalling, we compared deposited class A receptor–G-protein structures within four ligand-paired receptor systems and a five-ligand 5-HT1A–Gi1 series from Warren et al. (2024). Within-receptor structural alignment and atom-level contact fingerprints quantified the orthosteric receptor surface and the receptor-facing Gα interface. The four receptor pairs showed heterogeneous ligand-to-interface similarity rather than a uniform convergence gradient. For example, MC1R bound to α-MSH or SHU9119 exhibited lower ligand-contact similarity than Gα-interface similarity, whereas α2A showed the reverse pattern. In the 5-HT1A–Gi1 series, the reported same-assay maximal responses of buspirone, vilazodone and 4-F,5-MeO-PyrT were 93.4%, 97.4% and 102.8% of serotonin, respectively. Structural differences across the series are presented descriptively and stratified by dominant-negative Gαi1 construct. These results illustrate that ligand-dependent orthosteric interactions can coexist with comparable deposited mature interfaces in particular systems, but do not establish a universal relationship between endpoint geometry and agonist efficacy or the kinetics of productive coupling.', '',
           'Before integration: replace the word "comparable" only after examining exact completed five-structure measurements; if mixed-construct comparisons are driving the effect, omit that clause. This text is an editorial candidate, not a confirmed new biological discovery.', '',
           'Sources: DOI 10.1038/s41586-021-03376-8; 10.1038/s41586-024-07403-2; 10.1038/s41467-022-30929-w; 10.1038/s12276-024-01296-x; 10.1038/s41422-021-00557-y.']
(REPORT/'CANDIDATE_RESULTS_SUBSECTION_REVIEW_REQUIRED.md').write_text('\n'.join(paragraph)+'\n')
print('Integrated report written:',OUT/'INTEGRATED_SCIENTIFIC_REPORT.md')
