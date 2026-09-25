#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from bjp_revision_common import *

FEATURES=['tm6_r350_634_distance','tm6_aligned_displacement','tm5_r350_561_distance','tm7_displacement','dry_r350_630_min_distance','pif_550_644_distance','y753_displacement','icl2_helix_fraction','icl2_centroid_to_tm6_634','intracellular_cavity_core_volume','intracellular_cavity_core_mouth_area']

def analyse(label,paired,out):
    rows=[]; long=[]
    for f in FEATURES:
        dc=f'delta__{f}'
        if dc not in paired: continue
        # receptor means prevent multi-class receptors from receiving extra weight in pooled descriptives
        rr=paired.groupby('receptor_name',as_index=False).agg(delta=(dc,'mean'))
        vals=pd.to_numeric(rr.delta,errors='coerce').dropna(); pos=int((vals>0).sum()); neg=int((vals<0).sum())
        rows.append({'population':label,'feature':f,'n_receptors':len(vals),'median_delta':vals.median(),'q1':vals.quantile(.25),'q3':vals.quantile(.75),'mean_delta':vals.mean(),'positive':pos,'negative':neg,'ties':int((vals==0).sum()),'sign_test_p':sign_test_two_sided(pos,neg)})
        for cls,g in paired.groupby('transducer_family'):
            v=pd.to_numeric(g[dc],errors='coerce').dropna(); rows.append({'population':label,'feature':f,'captured_class':cls,'n_units':len(v),'median_delta':v.median(),'q1':v.quantile(.25),'q3':v.quantile(.75),'mean_delta':v.mean(),'positive':int((v>0).sum()),'negative':int((v<0).sum())})
        z=paired[[c for c in ['receptor_name','transducer_family',f'bound__{f}',f'receptor_only__{f}',dc] if c in paired]].copy(); z['feature']=f; long.append(z)
    tab=pd.DataFrame(rows); pooled=tab[tab.captured_class.isna()] if 'captured_class' in tab else tab
    q=bh_fdr(pooled.sign_test_p if 'sign_test_p' in pooled else []); tab.loc[pooled.index,'sign_test_q']=q
    write_tsv(tab,out/f'{label}_feature_summary.tsv'); write_tsv(pd.concat(long,ignore_index=True) if long else pd.DataFrame(),out/f'{label}_per_unit_values.tsv')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A04_structural_descriptives')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['p0_matched'])
    analyse('P0_strict',pd.read_csv(paths['p0_matched']),out)
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/paired_units.tsv'
    if p1.exists(): analyse('P1_terminal_strict',pd.read_csv(p1,sep='\t'),out)
    write_json(input_manifest(root,paths,['p0_matched'],'A04_structural_descriptives'),out/'manifest.json')
if __name__=='__main__': main()
