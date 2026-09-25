#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from bjp_revision_common import *

def seq_identity_rows(a,b,cols):
    aa=np.array([str(a.get(c,'')) for c in cols],dtype=object); bb=np.array([str(b.get(c,'')) for c in cols],dtype=object); valid=(aa!='')&(bb!='')&(aa!='nan')&(bb!='nan')
    return float((aa[valid]==bb[valid]).mean()) if valid.any() else np.nan

def repseq(frame,cols):
    rows=[]
    for rec,g in frame.groupby('receptor_name'):
        r={'receptor_name':rec}
        for c in cols:
            v=g[c].dropna().astype(str); r[c]=v.mode().iloc[0] if len(v) else np.nan
        rows.append(r)
    return pd.DataFrame(rows).set_index('receptor_name')

def analyse(label,pred,coefs,bound,out):
    free=pred[pred.experiment.eq('receptor_only_geometry')].copy(); free['classification_margin']=free.apply(classification_margin,axis=1); pcols=[f'probability_{c}' for c in PRIMARY_CLASSES]; free['predicted_class']=free[pcols].idxmax(axis=1).str.replace('probability_','',regex=False); free['correct']=free.predicted_class.eq(free.transducer_family)
    seqcols=[c for c in bound if c.startswith('residue_') and c.endswith('_identity')]; st=repseq(bound,seqcols)
    population=set(free.receptor_name.astype(str)); units=bound[bound.receptor_name.astype(str).isin(population)][['receptor_name','transducer_family']].drop_duplicates()
    rows=[]
    for _,r in free.iterrows():
        rec=str(r.receptor_name); fold=int(r.outer_fold); test_recs=set(free.loc[free.outer_fold.eq(fold),'receptor_name'].astype(str)); train_units=units[~units.receptor_name.astype(str).isin(test_recs)]
        if rec not in st.index: continue
        vals=[]
        for _,u in train_units.iterrows():
            tr=str(u.receptor_name)
            if tr not in st.index: continue
            vals.append((seq_identity_rows(st.loc[rec],st.loc[tr],seqcols),tr,str(u.transducer_family)))
        vals=[x for x in vals if np.isfinite(x[0])]
        same=[x for x in vals if x[2]==str(r.transducer_family)]
        best=max(vals,key=lambda x:x[0]) if vals else (np.nan,'',''); bestsame=max(same,key=lambda x:x[0]) if same else (np.nan,'','')
        rows.append({'population':label,'receptor_name':rec,'captured_class':r.transducer_family,'outer_fold':fold,'partner_free_margin':r.classification_margin,'partner_free_correct':r.correct,'nearest_any_identity':best[0],'nearest_any_receptor':best[1],'nearest_any_class':best[2],'nearest_same_class_identity':bestsame[0],'nearest_same_class_receptor':bestsame[1]})
    tab=pd.DataFrame(rows); write_tsv(tab,out/f'{label}_relatedness.tsv')
    rho,p,n=spearman_safe(tab.nearest_same_class_identity,tab.partner_free_margin); summary={'population':label,'n':len(tab),'spearman_margin_vs_nearest_same_class_identity':rho,'nominal_p':p,'n_for_correlation':n,'mean_nearest_same_class_identity_correct':float(tab.loc[tab.partner_free_correct,'nearest_same_class_identity'].mean()),'mean_nearest_same_class_identity_incorrect':float(tab.loc[~tab.partner_free_correct,'nearest_same_class_identity'].mean())}
    if coefs is not None and len(coefs):
        c=coefs[coefs.experiment.eq('receptor_only_geometry')].copy(); top=c.groupby('transformed_feature').absolute_coefficient.mean().sort_values(ascending=False).head(30).reset_index(); write_tsv(top,out/f'{label}_partner_free_top_coefficients.tsv')
    write_json(summary,out/f'{label}_summary.json')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A10_partner_free_relatedness')); args=ap.parse_args(); root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['p0_predictions','bound_frame']); bound=pd.read_csv(paths['bound_frame'])
    c0=paths['p0_imprint_dir']/'transducer_imprint_model_coefficients.csv'; analyse('P0_strict',pd.read_csv(paths['p0_predictions']),pd.read_csv(c0) if c0.exists() else None,bound,out)
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict'
    if (p1/'predictions.tsv').exists(): analyse('P1_terminal_strict',pd.read_csv(p1/'predictions.tsv',sep='\t'),pd.read_csv(p1/'coefficients.tsv',sep='\t') if (p1/'coefficients.tsv').exists() else None,bound,out)
    write_json(input_manifest(root,paths,['p0_predictions','bound_frame'],'A10_partner_free_relatedness'),out/'manifest.json')
if __name__=='__main__': main()
