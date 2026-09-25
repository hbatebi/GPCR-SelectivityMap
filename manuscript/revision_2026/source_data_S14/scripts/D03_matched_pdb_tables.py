#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np,pandas as pd
from bjp_revision_common import *

META=['experimental_method','resolution','species','ligand_type','fusion_protein','fusion_names','stabilising_mutations','engineered_construct']

def meta_for_ids(pool:pd.DataFrame,ids:set[str],prefix:str)->dict:
    g=pool[pool.pdb_id.astype(str).str.upper().isin(ids)].copy()
    out={prefix+'pdb_ids':unique_join(g.pdb_id.astype(str).str.upper())}
    chain='chain_id' if 'chain_id' in g else ('preferred_chain' if 'preferred_chain' in g else None)
    if chain: out[prefix+'pdb_chains']=unique_join((g.pdb_id.astype(str).str.upper()+':'+g[chain].fillna('').astype(str)).tolist())
    for c in META:
        if c not in g: continue
        if c=='resolution': out[prefix+'resolution_values']=unique_join(pd.to_numeric(g[c],errors='coerce').round(3))
        else: out[prefix+c]=unique_join(g[c])
    return out

def build(label,paired,bound,freepool,out):
    rows=[]
    for _,r in paired.iterrows():
        bpid=str(r.get('pdb_id','')).upper(); fids={x.strip().upper() for x in str(r.get('receptor_only_pdb_ids','')).split(';') if x.strip() and x!='nan'}
        row={'population':label,'receptor_name':r.receptor_name,'captured_transducer_class':r.transducer_family,'bound_pdb_id':bpid,'bound_receptor_chain':r.get('chain_id',r.get('preferred_chain','')),'partner_free_reference_definition':'median feature vector across all eligible structures listed in partner_free_pdb_ids','n_partner_free_structures':r.get('receptor_only_structure_count',len(fids))}
        row.update(meta_for_ids(bound,{bpid},'bound_')); row.update(meta_for_ids(freepool,fids,'partner_free_')); rows.append(row)
    tab=pd.DataFrame(rows); write_tsv(tab,out/f'{label}_matched_structure_table.tsv')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/D03_matched_pdb_tables')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['p0_matched','bound_frame','strict_frame','active_none'])
    bound=pd.read_csv(paths['bound_frame']); strict=pd.read_csv(paths['strict_frame']); active=pd.read_csv(paths['active_none'])
    for d in (bound,strict,active):
        d['pdb_id']=d.pdb_id.astype(str).str.upper()
        if 'chain_id' not in d: d['chain_id']=d.get('preferred_chain','')
    pool=pd.concat([strict,active],ignore_index=True,sort=False).drop_duplicates(['pdb_id','chain_id'])
    build('P0_strict',pd.read_csv(paths['p0_matched']),bound,pool,out)
    p1=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/paired_units.tsv'
    if p1.exists(): build('P1_terminal_strict',pd.read_csv(p1,sep='\t'),bound,pool,out)
    write_json(input_manifest(root,paths,['p0_matched','bound_frame','strict_frame','active_none'],'D03_matched_pdb_tables',extra={'reference_definition':'median across all eligible partner-free structures per receptor'}),out/'manifest.json')

if __name__=='__main__': main()
