#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial.distance import cosine, pdist
from bjp_revision_common import *

def analyse(label,paired,free_structures,features,out):
    dcols=[f'delta__{f}' for f in features if f'delta__{f}' in paired]
    # unit-level standardized deformation vectors
    D=paired[dcols].apply(pd.to_numeric,errors='coerce')
    med=D.median(); sd=D.std(ddof=0).replace(0,np.nan); Z=(D-med)/sd; Z=Z.fillna(0.0)
    pair_rows=[]
    for i in range(len(paired)):
        for j in range(i+1,len(paired)):
            vi=Z.iloc[i].to_numpy(float); vj=Z.iloc[j].to_numpy(float)
            ni=np.linalg.norm(vi); nj=np.linalg.norm(vj)
            cs=np.nan if ni==0 or nj==0 else float(np.dot(vi,vj)/(ni*nj))
            pair_rows.append({'i':i,'j':j,'receptor_i':paired.iloc[i].receptor_name,'receptor_j':paired.iloc[j].receptor_name,'class_i':paired.iloc[i].transducer_family,'class_j':paired.iloc[j].transducer_family,'same_class':paired.iloc[i].transducer_family==paired.iloc[j].transducer_family,'cosine_similarity':cs})
    pairtab=pd.DataFrame(pair_rows); write_tsv(pairtab,out/f'{label}_deformation_cosines.tsv')
    cos_summary=pairtab.groupby('same_class').cosine_similarity.agg(['count','mean','median']).reset_index().to_dict('records') if len(pairtab) else []
    # magnitude per unit/receptor
    mag=np.linalg.norm(Z.to_numpy(float),axis=1); mtab=paired[['receptor_name','transducer_family']].copy(); mtab['deformation_magnitude_z']=mag
    # partner-free within-receptor structural variability on same features, standardized with free global SD
    f=free_structures.copy(); available=[x for x in features if x in f]
    for x in available: f[x]=pd.to_numeric(f[x],errors='coerce')
    f_sd=f[available].std(ddof=0).replace(0,np.nan)
    varrows=[]
    for rec,g in f.groupby('receptor_name'):
        if len(g)<2: continue
        X=((g[available]-g[available].median())/f_sd).fillna(0.0).to_numpy(float)
        norms=np.linalg.norm(X,axis=1)
        varrows.append({'receptor_name':rec,'n_partner_free_structures':len(g),'partner_free_internal_variability_median':float(np.median(norms)),'partner_free_internal_variability_max':float(np.max(norms))})
    vartab=pd.DataFrame(varrows); mrec=mtab.groupby('receptor_name',as_index=False).deformation_magnitude_z.mean().merge(vartab,on='receptor_name',how='left')
    mrec['deformation_to_free_variability_ratio']=mrec.deformation_magnitude_z/mrec.partner_free_internal_variability_median.replace(0,np.nan); write_tsv(mrec,out/f'{label}_deformation_vs_free_variability.tsv')
    # within-class dispersion of bound versus partner-free endpoints in standardized 92D space
    disp=[]
    for endpoint,prefix in [('bound','bound__'),('partner_free','receptor_only__')]:
        cols=[prefix+x for x in features if prefix+x in paired]; X=paired[cols].apply(pd.to_numeric,errors='coerce'); X=(X-X.median())/X.std(ddof=0).replace(0,np.nan); X=X.fillna(0.0)
        for cls,idx in paired.groupby('transducer_family').groups.items():
            A=X.loc[idx].to_numpy(float); val=float(pdist(A,'euclidean').mean()) if len(A)>1 else np.nan
            disp.append({'population':label,'endpoint':endpoint,'captured_class':cls,'n_units':len(A),'mean_pairwise_distance':val})
    disptab=pd.DataFrame(disp); write_tsv(disptab,out/f'{label}_within_class_dispersion.tsv')
    summary={'population':label,'n_units':len(paired),'n_receptors':paired.receptor_name.nunique(),'cosine_summary':cos_summary,'median_deformation_magnitude_z':float(np.median(mag)),'n_receptors_with_free_variability':int(mrec.partner_free_internal_variability_median.notna().sum())}
    write_json(summary,out/f'{label}_summary.json')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A07_deformation_diagnostics')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['p0_matched','strict_frame']); feats=p0_base_features(paths)
    analyse('P0_strict',pd.read_csv(paths['p0_matched']),pd.read_csv(paths['strict_frame']),feats,out)
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/paired_units.tsv'; f1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict_free_inventory.tsv'
    if p1.exists():
        # Recover exactly the partner-free PDBs used in P1, not every active/no-transducer
        # structure from those receptors.
        p1df=pd.read_csv(p1,sep='\t')
        ids=set()
        for value in p1df.receptor_only_pdb_ids.dropna().astype(str):
            ids.update(x.strip().upper() for x in value.split(';') if x.strip())
        active=pd.read_csv(paths['active_none']); strict=pd.read_csv(paths['strict_frame'])
        if 'chain_id' not in active.columns: active['chain_id']=active.get('preferred_chain','')
        pool=pd.concat([strict,active],ignore_index=True,sort=False)
        pool['pdb_id']=pool.pdb_id.astype(str).str.upper()
        free=pool.loc[pool.pdb_id.isin(ids)].drop_duplicates(['pdb_id','chain_id'] if 'chain_id' in pool else ['pdb_id'])
        analyse('P1_terminal_strict',p1df,free,feats,out)
    write_json(input_manifest(root,paths,['p0_matched','strict_frame','active_none'],'A07_deformation_diagnostics'),out/'manifest.json')
if __name__=='__main__': main()
