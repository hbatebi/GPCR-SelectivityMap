#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs';F=ROOT/'figures';F.mkdir(exist_ok=True)
d=pd.read_csv(OUT/'five_ligand_structural_readouts.tsv',sep='\t')
fig,ax=plt.subplots(figsize=(7.4,4.6))
s=d.dropna(subset=['TM3.50_TM6.34_CA']);x=np.arange(len(s));bars=ax.bar(x,s['TM3.50_TM6.34_CA'])
ax.set_xticks(x,s.ligand,rotation=25,ha='right');ax.set_ylabel('TM3.50–TM6.34 Cα distance (Å)')
for rect,flag in zip(bars,s.Ga_construct_group):
 if flag=='DNG':rect.set_hatch('//')
ax.set_title('Warren 2024: 5-HT₁A–Gi₁ structural endpoints (descriptive)')
fig.tight_layout();fig.savefig(F/'five_ligand_tm3_tm6_CA.png',dpi=180);plt.close(fig)
q=OUT/'ten_pairwise_comparisons.tsv'
if q.is_file():
 p=pd.read_csv(q,sep='\t')
 fig,ax=plt.subplots(figsize=(6.7,4.4))
 for role,g in p.groupby('role'):
  ax.scatter(g.ligand_contact_jaccard,g.ga_contact_jaccard,label=role,s=70)
 ax.set(xlabel='Receptor ligand-contact Jaccard',ylabel='Receptor Gα-contact Jaccard',xlim=(0,1.02),ylim=(0,1.02))
 ax.plot([0,1],[0,1],linestyle='--',linewidth=1,color='gray');ax.legend(fontsize=8)
 fig.tight_layout();fig.savefig(F/'five_ligand_contact_overlap.png',dpi=180);plt.close(fig)
 n=len(p); nco=p.within_same_Ga_construct_group.sum()
else: n=nco=0
lines=[
'# Warren et al. (2024) five-ligand structural audit', '',
'STATUS: Exploratory descriptive comparison. NOT a validation of an efficacy mechanism.', '',
f'Coordinates processed: {int(d["TM3.50_TM6.34_CA"].notna().sum())}/5.',
f'Pairwise structural comparisons: {n}/10; same Gα construct group: {nco}.',
f'Numerically labelled Emax measurements: {d.Gi1_BRET_Emax_percent_5HT.notna().sum()}/5.', '',
'The baseline 93.4%, 97.4%, and 102.8% Emax measurements are from one Gi1 BRET assay in Warren et al., Nature (2024), doi:10.1038/s41586-024-07403-2. Only numbers marked verified_number in the evidence input are used. Missing values are not imputed.', '',
'Controls: 8FY8 and 8FYX are labeled DNG, whereas 8FYT, 8FYE, 8FYL are non-DNG in the publication structural table. Receptor structure L125(3.41)W and other engineered/stabilizing conditions mean the deposited structure is not a WT dynamic ensemble.', '',
'Ten pairwise distances from five unique structures DO NOT constitute ten independent biological replicates.',
'Contact Jaccard is sensitive to ligand size and the resolved atoms; does not measure activation kinetics.',
'Independent experimental data here are insufficient to infer which ligand contacts determine Gi1 efficacy.', '',
'RECOMMENDATION: do not change Nature Communications v5 main results on this basis alone. A substantially wider efficacy range or perturbation experiments would be needed.',
]
(OUT/'RESULTS_AND_EDITORIAL_VERDICT.md').write_text('\n'.join(lines)+'\n')
print('Figures and summary written:',F,OUT/'RESULTS_AND_EDITORIAL_VERDICT.md')
