#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from bjp_revision_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A06_g12_g13_census')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['active_transducers','strict_frame'])
    a=pd.read_csv(paths['active_transducers']); strict=pd.read_csv(paths['strict_frame']);
    fam='transducer_family' if 'transducer_family' in a else 'transducer_type'; g=a[a[fam].astype(str).eq('G12/13')].copy()
    strict_receptors=set(strict.receptor_name.astype(str)); g['has_strict_partner_free_active']=g.receptor_name.astype(str).isin(strict_receptors)
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict_free_inventory.tsv'
    p1receptors=set(pd.read_csv(p1,sep='\t').receptor_name.astype(str)) if p1.exists() else set()
    g['matchable_in_P1_terminal_strict']=g.receptor_name.astype(str).isin(p1receptors)
    cols=[c for c in ['pdb_id','chain_id','preferred_chain','receptor_name','species','experimental_method','resolution','year_of_deposition',fam,'has_strict_partner_free_active','matchable_in_P1_terminal_strict'] if c in g]
    write_tsv(g[cols].sort_values(['receptor_name','pdb_id']),out/'g12_g13_structural_census.tsv')
    summary={'n_structures':len(g),'n_receptors':g.receptor_name.nunique(),'receptors':sorted(g.receptor_name.astype(str).unique()),'n_with_strict_partner_free':int(g.loc[g.has_strict_partner_free_active,'receptor_name'].nunique()),'n_matchable_P1_terminal_strict':int(g.loc[g.matchable_in_P1_terminal_strict,'receptor_name'].nunique())}
    write_json(summary,out/'A06_summary.json'); write_json(input_manifest(root,paths,['active_transducers','strict_frame'],'A06_g12_g13_census'),out/'manifest.json'); print(summary)
if __name__=='__main__': main()
