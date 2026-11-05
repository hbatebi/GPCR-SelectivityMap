#!/usr/bin/env python3
"""Add GPCR generic-number microswitch geometry and same-assay measured Emax.
C-alpha distances only, to avoid residue-dependent sidechain missingness. Descriptive,
no claim of functional causal effect or statistics with only n=3 numerical efficacies.
"""
from pathlib import Path
import numpy as np, pandas as pd, gemmi, itertools
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
M=pd.read_csv(ROOT/'input'/'generic_number_mapping_1360.csv')
assay=pd.read_csv(ROOT/'input'/'warren2024_five_ligand_manifest.tsv',sep='\t')
# Update-only route: the dedicated evidence table is authoritative, so the author
# can add subsequently *verified* same-assay numbers without receiving a new ZIP.
evidence=pd.read_csv(ROOT/'input'/'5HT1A_functional_evidence.tsv',sep='\t')
verified=evidence[evidence.quantitative_status.isin(['verified_number','reference'])].set_index('pdb_id')['Emax_percent_of_serotonin']
assay['Gi1_BRET_Emax_percent_5HT']=assay.pdb_id.map(verified)

M['generic_number']=pd.to_numeric(M['generic_number'],errors='coerce')
P=[('TM3.50_TM6.34_CA',3.50,6.34),('PIF_TM5.50_TM6.44_CA',5.50,6.44),
   ('TM3.50_TM7.53_CA',3.50,7.53),('TM6.48_TM7.53_CA',6.48,7.53),
   ('TM3.40_TM5.50_CA',3.40,5.50)]
rows=[]
for r in assay.itertuples(index=False):
 cif=ROOT/'structures'/f'{r.pdb_id}.cif'
 if not cif.is_file(): continue
 s=gemmi.read_structure(str(cif));model=s[0];rec=model.find_chain('R')
 if rec is None: raise RuntimeError('Missing receptor chain R for '+r.pdb_id)
 by_num={int(x.seqid.num):x for x in rec if x.seqid.num is not None and any(a.name=='CA' for a in x)}
 sub=M[M.pdb_id.eq(r.pdb_id)]
 pos={round(float(x.generic_number),2):int(x.author_residue_number) for x in sub.itertuples(index=False)
      if pd.notna(x.generic_number) and pd.notna(x.author_residue_number)}
 dd={'pdb_id':r.pdb_id}
 for name,u,v in P:
  a,b=by_num.get(pos.get(u,-1)),by_num.get(pos.get(v,-1))
  if a is None or b is None: dd[name]=float('nan');continue
  c1=next((x for x in a if x.name=='CA'),None);c2=next((x for x in b if x.name=='CA'),None)
  dd[name]=float(np.linalg.norm(np.array([c1.pos.x-c2.pos.x,c1.pos.y-c2.pos.y,c1.pos.z-c2.pos.z]))) if c1 and c2 else float('nan')
 rows.append(dd)
geo=pd.DataFrame(rows)
joined=assay.merge(geo,on='pdb_id',how='left',validate='one_to_one')
atom=ROOT/'outputs'/'atomic_five'/'corrected_structure_metrics.tsv'
if atom.is_file():
 t=pd.read_csv(atom,sep='\t');joined=joined.merge(t[['pdb_id','n_ligand_heavy_atoms','n_receptor_ligand_contact_residues','n_receptor_ga_contact_residues','n_receptor_alpha5_contact_residues']],on='pdb_id',how='left',validate='one_to_one')
joined.to_csv(ROOT/'outputs'/'five_ligand_structural_readouts.tsv',sep='\t',index=False,float_format='%.7f')
# Pairwise structural differences and efficacy differences. NO regression: 10 pairwise comparisons are dependent.
paired=ROOT/'outputs'/'atomic_five'/'corrected_pair_comparisons.tsv'
if paired.is_file():
 p=pd.read_csv(paired,sep='\t')
 e=assay.set_index('pdb_id')['Gi1_BRET_Emax_percent_5HT'].to_dict()
 p['first_Emax']=p.first_pdb.map(e);p['second_Emax']=p.second_pdb.map(e)
 p['absolute_Emax_difference_pp']=(p.first_Emax-p.second_Emax).abs()
 p['within_same_Ga_construct_group']=p.role.ne('mixed_Ga_construct')
 for name,_,_ in P:
  x=geo.set_index('pdb_id')[name].to_dict()
  p['absolute_diff_'+name]=(p.first_pdb.map(x)-p.second_pdb.map(x)).abs()
 p.to_csv(ROOT/'outputs'/'ten_pairwise_comparisons.tsv',sep='\t',index=False,float_format='%.7f')
print(joined[['pdb_id','ligand','Gi1_BRET_Emax_percent_5HT']+[q[0] for q in P]].to_string(index=False))
print(f'NOTE: {assay.Gi1_BRET_Emax_percent_5HT.notna().sum()} of 5 same-assay numerical Emax values available; no correlation inference from pairwise distances.')
