#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from bjp_revision_common import *


def summarize_predictions(pred:pd.DataFrame, analysis:str)->pd.DataFrame:
    if 'outer_fold' not in pred: return pd.DataFrame()
    rows=[]
    for exp,g in pred.groupby('experiment'):
        for fold,te in g[g.outer_fold>=0].groupby('outer_fold'):
            test_recs=set(te.receptor_name.astype(str)); all_exp=g[g.outer_fold>=0]
            train=all_exp[~all_exp.receptor_name.astype(str).isin(test_recs)]
            row={'analysis':analysis,'experiment':exp,'outer_fold':int(fold),'test_units':len(te),'test_receptors':te.receptor_name.nunique(),'train_units_approx':len(train),'train_receptors_approx':train.receptor_name.nunique(),'test_receptor_names':unique_join(te.receptor_name)}
            for c in PRIMARY_CLASSES:
                row[f'test_{c}']=int((te.transducer_family.astype(str)==c).sum()) if 'transducer_family' in te else None
                row[f'train_{c}_approx']=int((train.transducer_family.astype(str)==c).sum()) if 'transducer_family' in train else None
            rows.append(row)
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/D02_fold_inventory')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root)
    blocks=[]
    if paths.get('p0_predictions') and paths['p0_predictions'].exists(): blocks.append(summarize_predictions(pd.read_csv(paths['p0_predictions']),'P0_strict_matched'))
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/predictions.tsv'
    if p1.exists(): blocks.append(summarize_predictions(pd.read_csv(p1,sep='\t'),'P1_terminal_strict_matched'))
    # Existing 207-unit sequence/interface audit predictions if present.
    seq=root/'results/microswitches/sequence_interface_audit_v1/predictions/out_of_fold_predictions.tsv'
    if seq.exists():
        q=pd.read_csv(seq,sep='\t')
        if 'outer_fold' in q.columns and 'experiment' in q.columns: blocks.append(summarize_predictions(q,'sequence_interface_audit'))
    result=pd.concat([x for x in blocks if len(x)],ignore_index=True) if any(len(x) for x in blocks) else pd.DataFrame()
    write_tsv(result,out/'fold_inventory.tsv')
    write_json({'analysis':'D02_fold_inventory','rows':len(result),'warning':'train counts marked approx are reconstructed from OOF fold membership. Use project audit tables as authoritative when a dedicated training-count field is available.'},out/'manifest.json')

if __name__=='__main__': main()
