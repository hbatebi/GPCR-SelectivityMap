#!/usr/bin/env python3
from __future__ import annotations
import argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.exceptions import ConvergenceWarning
from bjp_revision_common import *


def seq_identity(a,b,cols):
    aa=a[cols].astype(str).to_numpy(); bb=b[cols].astype(str).to_numpy(); valid=(aa!='nan')&(bb!='nan')&(aa!='')&(bb!='')
    return float((aa[valid]==bb[valid]).mean()) if valid.any() else np.nan

def receptor_seq_table(frame):
    cols=[c for c in frame.columns if c.startswith('residue_') and c.endswith('_identity')]
    rows=[]
    for r,g in frame.groupby('receptor_name'):
        row={'receptor_name':r}
        for c in cols:
            m=g[c].dropna().astype(str)
            row[c]=m.mode().iloc[0] if len(m) else np.nan
        rows.append(row)
    return pd.DataFrame(rows),cols

def fit_and_score(root,train,paired,features,variant,out,seed,quick):
    ProjectConfig,_analysis_settings,*_=ensure_project_imports(root)
    from gpcr_icl2.microswitches.publication_controls import _pipeline,_tune
    config=ProjectConfig.from_dir(root/'config'); settings=_analysis_settings(config,quick)
    if quick: settings['solver']='lbfgs'
    y=train.transducer_family.astype(str).to_numpy(); groups=train.receptor_name.astype(str).to_numpy()
    X=train[list(features)].apply(pd.to_numeric,errors='coerce')
    base=_pipeline(features,(),settings,seed)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',ConvergenceWarning)
        if quick: model=base.fit(X,y); params={'quick_fixed':True}
        else: model,params=_tune(base,X,y,groups,None,settings,seed)
    classes=tuple(str(x) for x in model.named_steps['model'].classes_)
    predparts=[]
    for exp,prefix in [('archive_scored_bound','bound__'),('archive_scored_free','receptor_only__')]:
        XX=paired[[prefix+f for f in features]].copy(); XX.columns=features; XX=XX.apply(pd.to_numeric,errors='coerce')
        pr=model.predict_proba(XX)
        p=paired[[c for c in ('pdb_id','chain_id','receptor_name','transducer_family') if c in paired]].copy(); p['experiment']=exp
        for j,c in enumerate(classes): p[f'probability_{c}']=pr[:,j]
        predparts.append(p)
    pred=pd.concat(predparts,ignore_index=True); rec,s=receptor_margin_table(pred,'archive_scored_bound','archive_scored_free'); s.update(paired_auc_from_predictions(pred,'archive_scored_bound','archive_scored_free')); s.update({'variant':variant,'n_training_units':len(train),'n_training_receptors':train.receptor_name.nunique(),'best_parameters':params})
    out.mkdir(parents=True,exist_ok=True); write_tsv(pred,out/'predictions.tsv'); write_tsv(rec,out/'receptor_margins.tsv'); write_json(s,out/'summary.json')
    return s

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A03_archive_trained_scoring')); ap.add_argument('--seed',type=int,default=20291729); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['bound_frame','p0_matched'])
    frame=pd.read_csv(paths['bound_frame']); features=p0_base_features(paths); frame=frame[frame.transducer_family.isin(PRIMARY_CLASSES)].copy()
    pops={'P0_strict':pd.read_csv(paths['p0_matched'])}
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/paired_units.tsv'
    if p1.exists(): pops['P1_terminal_strict']=pd.read_csv(p1,sep='\t')
    seqtab,seqcols=receptor_seq_table(frame); summaries=[]; close_rows=[]
    for pi,(label,paired) in enumerate(pops.items()):
        held=set(paired.receptor_name.astype(str)); train=frame[~frame.receptor_name.astype(str).isin(held)].copy()
        summaries.append(fit_and_score(root,train,paired,features,f'{label}__exclude_matched',out/label/'exclude_matched',args.seed+pi*100,args.quick))
        # second variant excludes any training receptor with >=50% identity to any held receptor, using generic-position sequence fingerprints.
        heldseq=seqtab[seqtab.receptor_name.isin(held)]; trseq=seqtab[~seqtab.receptor_name.isin(held)]
        exclude=set()
        for _,r in trseq.iterrows():
            vals=[seq_identity(r,h,seqcols) for _,h in heldseq.iterrows()]; mx=np.nanmax(vals) if np.isfinite(vals).any() else np.nan
            close_rows.append({'population':label,'training_receptor':r.receptor_name,'max_identity_to_matched':mx,'exclude_at_0.50':bool(np.isfinite(mx) and mx>=0.50)})
            if np.isfinite(mx) and mx>=0.50: exclude.add(r.receptor_name)
        train2=train[~train.receptor_name.isin(exclude)].copy()
        # only run if all 3 classes survive
        if set(train2.transducer_family.unique())==set(PRIMARY_CLASSES): summaries.append(fit_and_score(root,train2,paired,features,f'{label}__exclude_50pct_homologues',out/label/'exclude_50pct_homologues',args.seed+pi*100+17,args.quick))
    write_tsv(pd.DataFrame(close_rows),out/'close_homologue_audit.tsv'); write_tsv(pd.DataFrame(summaries),out/'A03_summary.tsv'); write_json(input_manifest(root,paths,['bound_frame','p0_matched'],'A03_archive_trained_scoring'),out/'manifest.json')
if __name__=='__main__': main()
