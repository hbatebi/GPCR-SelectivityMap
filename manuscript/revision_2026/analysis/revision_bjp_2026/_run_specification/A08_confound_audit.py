#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from bjp_revision_common import *

def aggregate_meta(df):
    def agg(g):
        res=pd.to_numeric(g.get('resolution'),errors='coerce')
        return pd.Series({'pdb_chains':unique_join((g.pdb_id.astype(str)+':'+g.get('chain_id',pd.Series('',index=g.index)).fillna('').astype(str)).tolist()),'methods':unique_join(g.get('experimental_method',pd.Series(dtype=str))),'resolution_median':res.median(),'species':unique_join(g.get('species',pd.Series(dtype=str))),'ligand_types':unique_join(g.get('ligand_type',pd.Series(dtype=str))),'any_fusion':bool(norm_bool_series(g.get('fusion_protein',pd.Series(False,index=g.index))).any()),'fusion_names':unique_join(g.get('fusion_names',pd.Series(dtype=str)))})
    rows=[]
    for rec,g in df.groupby('receptor_name'):
        z=agg(g).to_dict(); z['receptor_name']=rec; rows.append(z)
    return pd.DataFrame(rows)

def analyse(label,paired,pred,bound,free,out):
    rec,summary=receptor_margin_table(pred); b=bound[bound.pdb_id.astype(str).isin(set(paired.pdb_id.astype(str)))].copy(); bm=aggregate_meta(b).add_prefix('bound_').rename(columns={'bound_receptor_name':'receptor_name'})
    # select exact PDBs used by partner-free references
    ids=set();
    for x in paired.receptor_only_pdb_ids.dropna().astype(str): ids.update(x.split(';'))
    fm=aggregate_meta(free[free.pdb_id.astype(str).isin(ids)].copy()).add_prefix('free_').rename(columns={'free_receptor_name':'receptor_name'})
    tab=rec.merge(bm,on='receptor_name',how='left').merge(fm,on='receptor_name',how='left'); tab['resolution_difference_bound_minus_free']=tab.bound_resolution_median-tab.free_resolution_median; tab['absolute_resolution_difference']=tab.resolution_difference_bound_minus_free.abs(); tab['same_method']=tab.apply(lambda r: bool(set(str(r.bound_methods).split(';')) & set(str(r.free_methods).split(';'))),axis=1); tab['species_overlap']=tab.apply(lambda r: bool(set(str(r.bound_species).split(';')) & set(str(r.free_species).split(';'))),axis=1); tab['ligand_type_overlap']=tab.apply(lambda r: bool(set(str(r.bound_ligand_types).split(';')) & set(str(r.free_ligand_types).split(';'))),axis=1)
    write_tsv(tab,out/f'{label}_metadata_and_margin.tsv')
    assoc=[]
    for c in ['resolution_difference_bound_minus_free','absolute_resolution_difference']:
        rho,p,n=spearman_safe(tab[c],tab.margin_gain); assoc.append({'population':label,'variable':c,'n':n,'spearman_rho':rho,'nominal_p':p})
    write_tsv(pd.DataFrame(assoc),out/f'{label}_continuous_associations.tsv')
    sens=[]
    for name,mask in [('all',pd.Series(True,index=tab.index)),('same_method',tab.same_method),('no_bound_fusion',~tab.bound_any_fusion.fillna(False)),('ligand_type_overlap',tab.ligand_type_overlap)]:
        v=tab.loc[mask,'margin_gain'].dropna().to_numpy(float); sens.append({'population':label,'subset':name,'n':len(v),'positive':int((v>0).sum()),'mean_margin_gain':v.mean() if len(v) else np.nan,'signflip_p':exact_or_mc_signflip(v) if len(v) else np.nan})
    write_tsv(pd.DataFrame(sens),out/f'{label}_sensitivity_subsets.tsv'); summary.update({'same_method_n':int(tab.same_method.sum()),'species_overlap_n':int(tab.species_overlap.sum()),'ligand_overlap_n':int(tab.ligand_type_overlap.sum())}); write_json(summary,out/f'{label}_summary.json')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A08_confound_audit')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['p0_matched','p0_predictions','bound_frame','strict_frame','active_none']); bound=pd.read_csv(paths['bound_frame']); strict=pd.read_csv(paths['strict_frame']); active=pd.read_csv(paths['active_none']); freeall=pd.concat([strict,active],ignore_index=True,sort=False).drop_duplicates(['pdb_id','chain_id'] if 'chain_id' in active else ['pdb_id'])
    analyse('P0_strict',pd.read_csv(paths['p0_matched']),pd.read_csv(paths['p0_predictions']),bound,freeall,out)
    pdir=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict'
    if (pdir/'paired_units.tsv').exists() and (pdir/'predictions.tsv').exists(): analyse('P1_terminal_strict',pd.read_csv(pdir/'paired_units.tsv',sep='\t'),pd.read_csv(pdir/'predictions.tsv',sep='\t'),bound,freeall,out)
    write_json(input_manifest(root,paths,['p0_matched','p0_predictions','bound_frame','strict_frame','active_none'],'A08_confound_audit'),out/'manifest.json')
if __name__=='__main__': main()
