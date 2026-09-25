#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from bjp_revision_common import *

FOUR=['tm6_r350_634_distance','tm3_tm7_distance','pif_550_644_distance','y753_displacement']
CORE10=FOUR+['tm6_aligned_displacement','dry_r350_630_min_distance','icl2_helix_fraction','icl2_centroid_to_tm6_634','tm7_h8_angle','intracellular_cavity_core_mouth_area']

def analyse(root,out,label,paired,seed,quick):
    base=p0_base_features(discover_inputs(root)); rows=[]
    for i,(name,features) in enumerate([('four_feature_activation_code',FOUR),('core10_interpretable',CORE10)]):
        miss=[f for f in features if f not in base]
        if miss: rows.append({'population':label,'feature_set':name,'status':'missing_features','missing':';'.join(miss)}); continue
        cv=run_representation_cv(paired,features,root,seed+i*31,quick,nested=False,solver_override='lbfgs',names=(f'{name}__bound',f'{name}__free'))
        d=out/label/name; d.mkdir(parents=True,exist_ok=True)
        for k in ('predictions','performance','audit','coefficients'): write_tsv(cv[k],d/f'{k}.tsv')
        rec,s=receptor_margin_table(cv['predictions'],f'{name}__bound',f'{name}__free'); write_tsv(rec,d/'receptor_margins.tsv')
        s.update(paired_auc_from_predictions(cv['predictions'],f'{name}__bound',f'{name}__free')); s.update({'population':label,'feature_set':name,'features':features}); write_json(s,d/'summary.json'); rows.append(s)
    for k in (3,5):
        cv=run_pca_cv(paired,base,k,root,seed+200+k,quick); d=out/label/f'pca{k}'; d.mkdir(parents=True,exist_ok=True)
        for key in ('predictions','performance','audit'): write_tsv(cv[key],d/f'{key}.tsv')
        rec,s=receptor_margin_table(cv['predictions'],f'bound_pca{k}',f'receptor_only_pca{k}'); write_tsv(rec,d/'receptor_margins.tsv')
        s.update(paired_auc_from_predictions(cv['predictions'],f'bound_pca{k}',f'receptor_only_pca{k}')); s.update({'population':label,'feature_set':f'PCA{k}','k':k}); write_json(s,d/'summary.json'); rows.append(s)
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A02_low_dimensional')); ap.add_argument('--seed',type=int,default=20291729); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True)
    paths=discover_inputs(root); require_paths(paths,['p0_matched']); rows=[]
    p0=pd.read_csv(paths['p0_matched']); rows+=analyse(root,out,'P0_strict',p0,args.seed,args.quick)
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/paired_units.tsv'
    if p1.exists(): rows+=analyse(root,out,'P1_terminal_strict',pd.read_csv(p1,sep='\t'),args.seed+1000,args.quick)
    write_tsv(pd.DataFrame(rows),out/'A02_summary.tsv'); write_json(input_manifest(root,paths,['p0_matched'],'A02_low_dimensional_sensitivity'),out/'manifest.json')
if __name__=='__main__': main()
