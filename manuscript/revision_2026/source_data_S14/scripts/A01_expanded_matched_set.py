#!/usr/bin/env python3
from __future__ import annotations

import argparse, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from bjp_revision_common import *

FUSION_TERMS=("lysozyme","endolysin","b562","bril","rubredoxin","flavodoxin","thioredoxin","apocytochrome")
FUSION_ACCESSIONS={"P00720","P0ABE7","D9IEF7","P0A6H1"}
BAD_EXTRA_TERMS=("nanobody","antibody","fab","scfv","g protein","guanine nucleotide-binding","heavy chain","light chain","arrestin","grk")


def aslist(v):
    if v is None: return []
    return v if isinstance(v,list) else [v]


def classify_topology(cif_path: Path, chain: str, uniprot: str, desc: str="") -> tuple[str,str,float]:
    """Classify fusion location from mmCIF reference-sequence spans.

    Returns topology, evidence string, confidence in [0,1]. High-confidence N/C
    terminal calls require either explicit receptor+fusion reference spans or a
    fusion span at a chain terminus when the receptor is not database-mapped.
    """
    if not cif_path.exists(): return "missing_cif","coordinate file missing",0.0
    try: d=MMCIF2Dict(str(cif_path))
    except Exception as e: return "parse_error",str(e),0.0
    ids=aslist(d.get('_struct_ref.id')); codes=aslist(d.get('_struct_ref.db_code')); accs=aslist(d.get('_struct_ref.pdbx_db_accession')); dbn=aslist(d.get('_struct_ref.db_name'))
    ref={str(i):{'code':str(c),'acc':str(a),'db':str(db)} for i,c,a,db in zip(ids,codes,accs,dbn)}
    refids=aslist(d.get('_struct_ref_seq.ref_id')); strands=aslist(d.get('_struct_ref_seq.pdbx_strand_id')); begs=aslist(d.get('_struct_ref_seq.seq_align_beg')); ends=aslist(d.get('_struct_ref_seq.seq_align_end')); aacc=aslist(d.get('_struct_ref_seq.pdbx_db_accession'))
    seg=[]
    for ri,s,b,e,a in zip(refids,strands,begs,ends,aacc):
        if chain not in [x.strip() for x in str(s).split(',')]: continue
        try: b=int(float(b)); e=int(float(e))
        except: continue
        meta=ref.get(str(ri),{}); txt=(meta.get('code','')+' '+meta.get('acc','')+' '+str(a)).lower()
        kind='other'
        if str(uniprot) and str(uniprot) not in ('nan','None') and str(uniprot) in {meta.get('acc'),str(a)}: kind='receptor'
        if any(t in txt for t in FUSION_TERMS) or str(a) in FUSION_ACCESSIONS or meta.get('acc') in FUSION_ACCESSIONS: kind='fusion'
        seg.append((b,e,kind,txt))
    # chain sequence length
    L=None
    chains=aslist(d.get('_entity_poly.pdbx_strand_id')); seqs=aslist(d.get('_entity_poly.pdbx_seq_one_letter_code_can'))
    for ch,s in zip(chains,seqs):
        if chain in [x.strip() for x in str(ch).split(',')]:
            L=len(re.sub(r'\s','',str(s))); break
    fs=[x for x in seg if x[2]=='fusion']; rs=[x for x in seg if x[2]=='receptor']
    evidence=f"segments={seg}; chain_length={L}; description={desc}"
    if not fs:
        # `cand` contains only records already independently classified as
        # fusion constructs by the entity audit.  mmCIF often names the
        # inserted protein with an abbreviated database code (for example
        # RUBR_CLOPA or FLAV_DESVH), so keyword matching alone can miss the
        # fusion span.  When a non-receptor reference span is geometrically
        # unambiguous relative to receptor reference spans, use that topology
        # information without consulting any model outcome.
        other=[x for x in seg if x[2]=='other']
        if rs and other:
            for ob,oe,_,_ in other:
                before=any(x[1] < ob for x in rs)
                after=any(x[0] > oe for x in rs)
                if before and after:
                    return "internal",evidence+"; inferred_from_nonreceptor_span=internal",0.95
            first_r=min(x[0] for x in rs); last_r=max(x[1] for x in rs)
            if all(x[1] < first_r for x in other):
                return "N_terminal",evidence+"; inferred_from_nonreceptor_span=N_terminal",0.90
            if all(x[0] > last_r for x in other):
                return "C_terminal",evidence+"; inferred_from_nonreceptor_span=C_terminal",0.90
        return "manual_review_no_fusion_span",evidence,0.0
    fb,fe=min(x[0] for x in fs),max(x[1] for x in fs)
    if rs:
        rb,re_=min(x[0] for x in rs),max(x[1] for x in rs)
        before=any(x[1] < fb for x in rs); after=any(x[0] > fe for x in rs)
        if before and after: return "internal",evidence,1.0
        if fe < rb: return "N_terminal",evidence,1.0
        if fb > re_: return "C_terminal",evidence,1.0
        return "manual_review_overlap",evidence,0.2
    if L:
        if fb <= max(12,int(0.2*L)) and fe < 0.48*L: return "N_terminal",evidence,0.75
        if fe >= min(L-12,int(0.8*L)) and fb > 0.52*L: return "C_terminal",evidence,0.75
    return "manual_review_no_receptor_span",evidence,0.0


def extra_polymer_policy(desc) -> tuple[bool,str]:
    if pd.isna(desc) or not str(desc).strip(): return True,"none"
    s=str(desc).lower()
    if any(t in s for t in BAD_EXTRA_TERMS): return False,"excluded_auxiliary_or_transducer"
    return True,"extracellular_polymer_ligand_permitted_only_in_ligand_permissive_variant"


def run_population(root: Path, out: Path, name: str, bound: pd.DataFrame, free: pd.DataFrame, features: list[str], seed: int, quick: bool, original_receptors:set[str]):
    paired=make_paired(bound,free,features)
    popdir=out/name; popdir.mkdir(parents=True,exist_ok=True)
    coverage=[]
    for f in features:
        b=pd.to_numeric(paired.get('bound__'+f),errors='coerce')
        a=pd.to_numeric(paired.get('receptor_only__'+f),errors='coerce')
        coverage.append({'feature':f,'bound_nonmissing_fraction':float(b.notna().mean()),'partner_free_nonmissing_fraction':float(a.notna().mean()),'both_nonmissing_fraction':float((b.notna()&a.notna()).mean()),'complete_all_units':bool((b.notna()&a.notna()).all())})
    cov=pd.DataFrame(coverage); write_tsv(cov,popdir/'feature_coverage.tsv')
    write_tsv(paired,popdir/'paired_units.tsv')
    cv=run_representation_cv(paired,features,root,seed,quick,nested=not quick)
    for key in ('predictions','performance','audit','coefficients'): write_tsv(cv[key],popdir/f'{key}.tsv')
    rec,summary=receptor_margin_table(cv['predictions'])
    rec['original_P0_receptor']=rec.receptor_name.isin(original_receptors)
    write_tsv(rec,popdir/'receptor_margins.tsv')
    summary.update(paired_auc_from_predictions(cv['predictions'],'bound_geometry','receptor_only_geometry'))
    summary.update({"population":name,"n_units":int(len(paired)),"n_common_features":len(features),"n_features_complete_all_units":int(cov.complete_all_units.sum()),"n_free_structures":int(free.pdb_id.nunique()),"n_new_receptors":int(len(set(paired.receptor_name)-original_receptors))})
    # Report the new-receptor subset two ways: as a subset of the expanded OOF predictions,
    # and, where class coverage permits, as a fully independent CV analysis fit only on
    # receptors that were not present in P0. The latter is the prespecified replication.
    new=rec.loc[~rec.original_P0_receptor].copy()
    if len(new):
        vals=new.margin_gain.to_numpy(float); lo,hi=bootstrap_mean_ci(vals,200 if quick else 2000,seed+77)
        summary['new_receptors_subset_of_expanded_fit']={"n":len(new),"positive":int((vals>0).sum()),"mean":float(vals.mean()),"ci_low":lo,"ci_high":hi,"signflip_p":exact_or_mc_signflip(vals,seed+88)}
        new_paired=paired.loc[~paired.receptor_name.astype(str).isin(original_receptors)].copy()
        repdir=popdir/'new_receptors_only_independent'
        repdir.mkdir(parents=True,exist_ok=True)
        try:
            if new_paired.receptor_name.nunique() < 6 or set(new_paired.transducer_family.astype(str)) != set(PRIMARY_CLASSES):
                raise RuntimeError('insufficient receptor/class coverage for independent three-class CV')
            repcv=run_representation_cv(new_paired,features,root,seed+707,quick,nested=not quick)
            for key in ('predictions','performance','audit','coefficients'): write_tsv(repcv[key],repdir/f'{key}.tsv')
            reprec,reps=receptor_margin_table(repcv['predictions']); write_tsv(reprec,repdir/'receptor_margins.tsv')
            reps.update(paired_auc_from_predictions(repcv['predictions'],'bound_geometry','receptor_only_geometry'))
            reps.update({'n_units':int(len(new_paired)),'n_receptors':int(new_paired.receptor_name.nunique()),'status':'ok'})
            write_json(reps,repdir/'summary.json'); summary['new_receptors_only_independent']=reps
        except Exception as e:
            reps={'status':'not_estimable','n_units':int(len(new_paired)),'n_receptors':int(new_paired.receptor_name.nunique()),'reason':f'{type(e).__name__}: {e}'}
            write_json(reps,repdir/'summary.json'); summary['new_receptors_only_independent']=reps
    write_json(summary,popdir/'summary.json')
    return paired,summary


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A01_expanded_matched_set')); ap.add_argument('--seed',type=int,default=20291729); ap.add_argument('--topology-overrides',type=Path,default=None); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True)
    paths=discover_inputs(root); require_paths(paths,['bound_frame','active_none','strict_frame','entity_audit','p0_matched'])
    bound=pd.read_csv(paths['bound_frame']); active=pd.read_csv(paths['active_none']); strict=pd.read_csv(paths['strict_frame']); audit=pd.read_csv(paths['entity_audit']); p0=pd.read_csv(paths['p0_matched'])
    features=p0_base_features(paths); original=set(p0.receptor_name.astype(str))
    for d in (bound,active,strict,audit):
        if 'pdb_id' in d: d['pdb_id']=d.pdb_id.astype(str).str.upper()
        if 'receptor_name' in d: d['receptor_name']=d.receptor_name.astype(str)
    if 'transducer_family' not in bound and 'transducer_type' in bound: bound['transducer_family']=bound.transducer_type
    # merge audit with active feature frame
    keep=[c for c in ['pdb_id','receptor_name','chain_id','preferred_chain','uniprot','activation_state','fusion_protein','fusion_names','experimental_method','resolution','ligand_type','species']+features if c in active]
    m=audit.merge(active[keep].drop_duplicates(['pdb_id','receptor_name']),on=['pdb_id','receptor_name'],how='left',suffixes=('','_feature'))
    cand=m[(m['final_audit_category'].astype(str)=='fusion_construct') & (m['activation_state'].astype(str).str.lower().isin(['active','active_like'])) & m.receptor_name.isin(set(bound.receptor_name))].copy()
    rows=[]
    for _,r in cand.iterrows():
        chain=str(r.get('chain_id') if pd.notna(r.get('chain_id')) else r.get('preferred_chain'))
        cif=Path(str(r.get('cif_path',''))); cif = cif if cif.is_absolute() else root/cif
        top,evidence,conf=classify_topology(cif,chain,str(r.get('uniprot','')),str(r.get('receptor_entity_descriptions','')))
        extra_ok,extra_reason=extra_polymer_policy(r.get('extra_polymer_descriptions'))
        rows.append({"pdb_id":r.pdb_id,"receptor_name":r.receptor_name,"chain_id":chain,"fusion_names":r.get('fusion_names'),"topology":top,"topology_confidence":conf,"extra_polymer_descriptions":r.get('extra_polymer_descriptions'),"extra_ok_ligand_permissive":extra_ok,"extra_policy_reason":extra_reason,"evidence":evidence})
    topo=pd.DataFrame(rows)
    if len(topo):
        topo['topology_source']='automatic_mmcif'
        override=args.topology_overrides if args.topology_overrides and args.topology_overrides.is_absolute() else (root/args.topology_overrides if args.topology_overrides else root/'revision_inputs/fusion_topology_overrides.tsv')
        if override.exists():
            ov=pd.read_csv(override,sep='\t' if override.suffix.lower() in ['.tsv','.txt'] else ',')
            required={'pdb_id','receptor_name','topology'}
            if not required.issubset(ov.columns): raise ValueError(f'Fusion override file must contain {sorted(required)}')
            ov['pdb_id']=ov.pdb_id.astype(str).str.upper(); ov['receptor_name']=ov.receptor_name.astype(str)
            allowed={'N_terminal','C_terminal','internal','exclude','unresolved'}
            bad=set(ov.topology.astype(str))-allowed
            if bad: raise ValueError(f'Invalid topology override values: {sorted(bad)}')
            for _,rr in ov.iterrows():
                mask=topo.pdb_id.astype(str).str.upper().eq(rr.pdb_id) & topo.receptor_name.astype(str).eq(rr.receptor_name)
                if mask.any():
                    topo.loc[mask,'topology']=str(rr.topology); topo.loc[mask,'topology_confidence']=1.0 if str(rr.topology) in {'N_terminal','C_terminal','internal','exclude'} else 0.0; topo.loc[mask,'topology_source']='manual_override'
                    if 'rationale' in ov: topo.loc[mask,'override_rationale']=str(rr.get('rationale',''))
        write_tsv(topo,out/'fusion_topology_audit.tsv')
    else:
        write_tsv(topo,out/'fusion_topology_audit.tsv')
    if len(topo):
        manual=topo.loc[~topo.topology.isin(['N_terminal','C_terminal','internal']) | (pd.to_numeric(topo.topology_confidence,errors='coerce')<0.70)].copy()
        write_tsv(manual,out/'fusion_topology_manual_review.tsv')
    active2=active.merge(topo[['pdb_id','receptor_name','topology','topology_confidence','extra_ok_ligand_permissive']],on=['pdb_id','receptor_name'],how='left')
    if 'chain_id' not in active2.columns:
        active2['chain_id']=active2.get('preferred_chain',pd.Series('',index=active2.index))
    terminal=active2.topology.isin(['N_terminal','C_terminal']) & (pd.to_numeric(active2.topology_confidence,errors='coerce')>=0.70)
    # exact planned population: only fusion criterion relaxed, so no other polymer entities are allowed.
    audit_extra=audit[['pdb_id','receptor_name','extra_polymer_entity_ids']].copy()
    active2=active2.merge(audit_extra,on=['pdb_id','receptor_name'],how='left')
    no_extra=active2.extra_polymer_entity_ids.isna() | active2.extra_polymer_entity_ids.astype(str).isin(['','nan'])
    free_main=pd.concat([strict,active2.loc[terminal & no_extra,strict.columns.intersection(active2.columns)]],ignore_index=True,sort=False).drop_duplicates(['pdb_id','chain_id'] if 'chain_id' in strict else ['pdb_id'])
    # more permissive variant: terminal fusion plus extracellular polymer ligand, but never binder/transducer.
    free_lp=pd.concat([strict,active2.loc[terminal & active2.extra_ok_ligand_permissive.fillna(False),strict.columns.intersection(active2.columns)]],ignore_index=True,sort=False).drop_duplicates(['pdb_id','chain_id'] if 'chain_id' in strict else ['pdb_id'])
    # N-terminal-only sensitivity using same strict extra-polymer rule
    nterm=active2.topology.eq('N_terminal') & (pd.to_numeric(active2.topology_confidence,errors='coerce')>=0.70) & no_extra
    free_n=pd.concat([strict,active2.loc[nterm,strict.columns.intersection(active2.columns)]],ignore_index=True,sort=False).drop_duplicates(['pdb_id','chain_id'] if 'chain_id' in strict else ['pdb_id'])
    # Exploratory all-fusion set: no auxiliary/transducer extra polymer; may include internal fusions. Not for primary claims.
    allfus=active2['topology'].isin(['N_terminal','C_terminal','internal']) & active2.extra_ok_ligand_permissive.fillna(False)
    free_all=pd.concat([strict,active2.loc[allfus,strict.columns.intersection(active2.columns)]],ignore_index=True,sort=False).drop_duplicates(['pdb_id','chain_id'] if 'chain_id' in strict else ['pdb_id'])

    summaries={}
    for i,(name,free) in enumerate([('P1_terminal_strict',free_main),('P1_terminal_ligand_permissive',free_lp),('P1_N_terminal_strict',free_n),('P1_all_fusion_exploratory',free_all)]):
        try:
            _,s=run_population(root,out,name,bound,free,features,args.seed+i*101,args.quick,original); summaries[name]=s
        except Exception as e:
            summaries[name]={"status":"failed","error":f"{type(e).__name__}: {e}"}
    # save exact free inventories used
    for name,free in [('P1_terminal_strict',free_main),('P1_terminal_ligand_permissive',free_lp),('P1_N_terminal_strict',free_n),('P1_all_fusion_exploratory',free_all)]:
        cols=[c for c in ['pdb_id','chain_id','receptor_name','experimental_method','resolution','ligand_type','fusion_protein','fusion_names'] if c in free]
        write_tsv(free[cols].sort_values(['receptor_name','pdb_id']),out/f'{name}_free_inventory.tsv')
    write_json({"analysis":"A01_expanded_matched_set","P0_n_receptors":len(original),"topology_counts":topo.topology.value_counts(dropna=False).to_dict(),"populations":summaries,"note":"P1_terminal_strict is the prespecified sensitivity. Ligand-permissive and all-fusion populations are secondary/exploratory and must not replace P0."},out/'A01_summary.json')
    write_json(input_manifest(root,paths,['bound_frame','active_none','strict_frame','entity_audit','p0_matched'],'A01_expanded_matched_set'),out/'manifest.json')
    print((out/'A01_summary.json').read_text())
if __name__=='__main__': main()
