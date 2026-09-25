#!/usr/bin/env python3
from __future__ import annotations
import argparse
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd
from bjp_revision_common import *


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--project-root',type=Path,default=Path.cwd())
    ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/R14_same_receptor_408'))
    args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True)
    paths=discover_inputs(root); require_paths(paths,['same_receptor','bound_frame'])
    comp=pd.read_csv(paths['same_receptor']); bound=pd.read_csv(paths['bound_frame'])
    comp=comp.copy(); comp['comparison']=comp['class_a'].astype(str)+' vs '+comp['class_b'].astype(str)
    write_tsv(comp,out/'all_408_feature_by_class_pair_comparisons.tsv')
    summary=pd.DataFrame([{
        'n_comparisons':len(comp),
        'n_features':comp.feature.nunique(),
        'n_class_pairs':comp.comparison.nunique(),
        'minimum_nominal_p':pd.to_numeric(comp.p_value,errors='coerce').min(),
        'minimum_fdr_q':pd.to_numeric(comp.p_fdr,errors='coerce').min(),
        'n_fdr_below_0_05':int((pd.to_numeric(comp.p_fdr,errors='coerce')<0.05).sum()),
        'paired_receptor_n_min':int(pd.to_numeric(comp.n_paired_receptors,errors='coerce').min()),
        'paired_receptor_n_max':int(pd.to_numeric(comp.n_paired_receptors,errors='coerce').max()),
    }])
    write_tsv(summary,out/'comparison_summary.tsv')

    # Reconstruct which receptors can contribute to each class-pair comparison.
    unit=bound[['receptor_name','transducer_family','pdb_id']+[c for c in ['preferred_chain','chain_id'] if c in bound]].drop_duplicates()
    inv=[]
    for rec,g in unit.groupby('receptor_name'):
        classes=sorted(set(g.transducer_family.astype(str)) & set(PRIMARY_CLASSES))
        if len(classes)<2: continue
        for a,b in combinations(classes,2):
            ga=g[g.transducer_family.astype(str).eq(a)]; gb=g[g.transducer_family.astype(str).eq(b)]
            chaincol='chain_id' if 'chain_id' in g else ('preferred_chain' if 'preferred_chain' in g else None)
            inv.append({
                'receptor_name':rec,'class_a':a,'class_b':b,'comparison':a+' vs '+b,
                'class_a_pdbs':unique_join(ga.pdb_id),'class_b_pdbs':unique_join(gb.pdb_id),
                'class_a_chains':unique_join(ga[chaincol]) if chaincol else '',
                'class_b_chains':unique_join(gb[chaincol]) if chaincol else '',
            })
    inv=pd.DataFrame(inv)
    write_tsv(inv,out/'multitransducer_receptor_class_pair_inventory.tsv')
    if len(inv):
        counts=inv.groupby('comparison',as_index=False).agg(n_receptors=('receptor_name','nunique'),receptors=('receptor_name',unique_join))
        write_tsv(counts,out/'receptors_per_class_pair.tsv')
    write_json(input_manifest(root,paths,['same_receptor','bound_frame'],'R14_same_receptor_408',extra={'n_rows':len(comp)}),out/'manifest.json')

if __name__=='__main__': main()
