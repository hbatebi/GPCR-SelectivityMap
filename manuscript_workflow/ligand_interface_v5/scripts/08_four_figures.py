from pathlib import Path
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'figures';out.mkdir(exist_ok=True)
d=pd.read_csv(ROOT/'outputs/four_receptor_atomic/corrected_pair_comparisons.tsv',sep='\t');p=d[d.role=='primary'].copy();p['name']=p.case_id.map({'5HT1A':'5-HT1A','DRD1':'D1','ADRA2A':'α2A','MC1R':'MC1R'})
x=np.arange(len(p));w=.36
fig,ax=plt.subplots(figsize=(8,4.5));ax.bar(x-w/2,p.ligand_contact_jaccard,w,label='Ligand contacts');ax.bar(x+w/2,p.ga_contact_jaccard,w,label='Gα contacts');ax.set_xticks(x,p.name);ax.set_ylabel('Receptor-residue Jaccard similarity, 4.5 Å');ax.set_ylim(0,1.1);ax.legend();ax.set_title('Four matched receptor–transducer case studies');fig.tight_layout();fig.savefig(out/'four_receptor_contact_comparisons.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.5));ax.plot(x,p.alpha5_tip_shift_A,'o-',label='Default receptor-core alignment');ax.plot(x,p.alpha5_tip_shift_alt_alignment_A,'s--',label='Alternative core');ax.set_xticks(x,p.name);ax.set_ylabel('α5 tip displacement (Å)');ax.legend();fig.tight_layout();fig.savefig(out/'four_receptor_alignment_sensitivity.png',dpi=180);plt.close(fig)
