#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from bjp_revision_common import *

FAMS=['Gs','Gi/o','Gq/11','G12/13']

def load_annotations(root,paths,external):
    parts=[]
    # Existing curated receptor-only audit is a valid local snapshot for receptors represented there.
    if paths.get('entity_audit') and paths['entity_audit'].exists():
        a=pd.read_csv(paths['entity_audit']); cols=[c for c in ['receptor_name']+FAMS+['annotation_source','annotation_notes'] if c in a]
        a=a[cols].drop_duplicates('receptor_name'); parts.append(a)
    if paths.get('functional_provenance') and paths['functional_provenance'].exists():
        p=pd.read_csv(paths['functional_provenance'],sep='\t'); cols=[c for c in ['receptor_name']+FAMS+['annotation_source','notes'] if c in p]; parts.append(p[cols].drop_duplicates('receptor_name'))
    ext=external if external and external.is_absolute() else (root/external if external else root/'revision_inputs/functional_coupling_annotations.tsv')
    if ext.exists():
        e=pd.read_csv(ext,sep='\t' if ext.suffix.lower() in ['.tsv','.txt'] else ','); parts.append(e)
    if not parts: return pd.DataFrame(columns=['receptor_name']+FAMS)
    # Later sources override earlier ones cell-by-cell.
    allrec=sorted(set().union(*(set(x.receptor_name.astype(str)) for x in parts)))
    out=pd.DataFrame({'receptor_name':allrec}).set_index('receptor_name')
    for p in parts:
        p=p.copy(); p['receptor_name']=p.receptor_name.astype(str); p=p.set_index('receptor_name')
        for c in p.columns:
            # Pandas >=2.2/3.x no longer permits assigning strings into a
            # float64 column created from np.nan.  Keep the merge layer
            # dtype-neutral; family flags are normalised to numeric below.
            if c not in out:
                out[c] = pd.Series(pd.NA, index=out.index, dtype="object")
            else:
                out[c] = out[c].astype("object")
            vals = p[c].astype("object")
            out.loc[p.index, c] = vals.to_numpy()
    return out.reset_index()

def run_geometry_subset(root,frame,features,name,out,seed,quick):
    ProjectConfig,_analysis_settings,ExperimentSpec,FeatureSpec,_,_,cross_validate_experiment=ensure_project_imports(root)
    config=ProjectConfig.from_dir(root/'config'); settings=_analysis_settings(config,quick)
    ex=ExperimentSpec(name,FeatureSpec(name,tuple(features)),nested_tuning=not quick)
    res=cross_validate_experiment(frame,ex,settings,seed,minimum_fold_coverage=0.0,bootstrap_iterations=100 if quick else 2000)
    write_tsv(res['predictions'],out/f'{name}_predictions.tsv'); write_tsv(res['audit'],out/f'{name}_fold_audit.tsv'); write_json(dict(res['performance']),out/f'{name}_performance.json')
    return dict(res['performance'])

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A05_functional_coupling_repertoire')); ap.add_argument('--annotation-file',type=Path,default=None); ap.add_argument('--seed',type=int,default=20291729); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True); paths=discover_inputs(root); require_paths(paths,['bound_frame','p0_matched'])
    bound=pd.read_csv(paths['bound_frame']); ann=load_annotations(root,paths,args.annotation_file)
    # Normalize family flags: 1/true = documented, 0/false = documented absent, blank = unknown.
    for c in FAMS:
        if c not in ann: ann[c]=np.nan
        raw=ann[c]
        numeric=pd.to_numeric(raw,errors='coerce')
        # accept common boolean/string exports without converting unknowns to negatives
        txt=raw.astype(str).str.strip().str.lower()
        numeric=numeric.mask(numeric.isna() & txt.isin(['true','yes','y']),1.0)
        numeric=numeric.mask(numeric.isna() & txt.isin(['false','no','n']),0.0)
        ann[c]=numeric
    ann['n_known_families']=ann[FAMS].notna().sum(axis=1)
    ann['n_documented_positive_families']=ann[FAMS].eq(1).sum(axis=1).astype(float)
    ann.loc[ann.n_known_families.eq(0),'n_documented_positive_families']=np.nan
    analysed=['Gs','Gi/o','Gq/11']
    ann['n_known_analyzed_classes']=ann[analysed].notna().sum(axis=1)
    ann['n_documented_positive_analyzed_classes']=ann[analysed].eq(1).sum(axis=1).astype(float)
    ann.loc[ann.n_known_analyzed_classes.eq(0),'n_documented_positive_analyzed_classes']=np.nan
    write_tsv(ann,out/'functional_annotation_snapshot.tsv')
    x=bound.merge(ann,on='receptor_name',how='left',suffixes=('','__functional'))
    x['captured_family_documented']=x.apply(lambda r: r.get(str(r.transducer_family),np.nan)==1 if str(r.transducer_family) in FAMS else np.nan,axis=1)
    x['functional_promiscuous_2plus']=x.n_documented_positive_families>=2
    write_tsv(x[[c for c in ['pdb_id','chain_id','receptor_name','transducer_family','n_documented_positive_families','n_known_families','captured_family_documented','functional_promiscuous_2plus','primary_family','annotation_source','annotation_notes','notes'] if c in x]],out/'bound_units_functional_annotation.tsv')
    unit_coverage=x.n_known_families.gt(0)
    complete3=x.n_known_analyzed_classes.eq(3)
    complete4=x.n_known_families.eq(4)
    descript={
      'n_bound_units':len(x),'n_bound_receptors':x.receptor_name.nunique(),
      'n_units_with_any_functional_annotation':int(unit_coverage.sum()),'n_receptors_with_any_functional_annotation':int(x.loc[unit_coverage,'receptor_name'].nunique()),
      'n_units_complete_for_three_analyzed_classes':int(complete3.sum()),
      'n_receptors_complete_for_three_analyzed_classes':int(x.loc[complete3,'receptor_name'].nunique()),
      'complete3_units_from_receptors_with_1_analyzed_class':int((complete3 & x.n_documented_positive_analyzed_classes.eq(1)).sum()),
      'complete3_units_from_receptors_with_2_analyzed_classes':int((complete3 & x.n_documented_positive_analyzed_classes.eq(2)).sum()),
      'complete3_units_from_receptors_with_3_analyzed_classes':int((complete3 & x.n_documented_positive_analyzed_classes.eq(3)).sum()),
      'complete3_units_from_promiscuous_receptors_2plus':int((complete3 & x.n_documented_positive_analyzed_classes.ge(2)).sum()),
      'n_units_complete_for_all_four_families':int(complete4.sum()),
      'complete4_units_from_receptors_with_1_family':int((complete4 & x.n_documented_positive_families.eq(1)).sum()),
      'complete4_units_from_receptors_with_2_families':int((complete4 & x.n_documented_positive_families.eq(2)).sum()),
      'complete4_units_from_receptors_with_3_families':int((complete4 & x.n_documented_positive_families.eq(3)).sum()),
      'complete4_units_from_receptors_with_4_families':int((complete4 & x.n_documented_positive_families.eq(4)).sum()),
    }
    # P0/P1 matched descriptives
    for label,p in [('P0_strict',paths['p0_matched']),('P1_terminal_strict',root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict/paired_units.tsv')]:
        if not p.exists(): continue
        m=pd.read_csv(p,sep='\t' if p.suffix=='.tsv' else ',').merge(ann,on='receptor_name',how='left')
        rr=m[['receptor_name','n_known_analyzed_classes','n_documented_positive_analyzed_classes']].drop_duplicates()
        complete=rr.n_known_analyzed_classes.eq(3)
        descript[label]={'n_receptors':rr.receptor_name.nunique(),'receptors_with_any_analyzed_class_annotation':int(rr.n_known_analyzed_classes.gt(0).sum()),'receptors_complete_for_three_analyzed_classes':int(complete.sum()),'complete_receptors_promiscuous_2plus':int((complete & rr.n_documented_positive_analyzed_classes.ge(2)).sum())}
        write_tsv(m[[c for c in ['receptor_name','transducer_family','n_documented_positive_families','n_documented_positive_analyzed_classes','n_known_analyzed_classes']+FAMS+['primary_family','annotation_source'] if c in m]],out/f'{label}_functional_annotation.tsv')
    # Sensitivity model on receptors documented as single-family among the analysed three.
    features=p0_base_features(paths); perf=[]
    single=x[(x[['Gs','Gi/o','Gq/11']].eq(1).sum(axis=1)==1) & x[['Gs','Gi/o','Gq/11']].notna().all(axis=1)].copy()
    if single.receptor_name.nunique()>=12 and set(single.transducer_family.unique())==set(PRIMARY_CLASSES): perf.append(run_geometry_subset(root,single,features,'single_family_geometry',out,args.seed,args.quick))
    if 'primary_family' in x:
        primary=x[x.primary_family.astype(str).eq(x.transducer_family.astype(str))].copy()
        if primary.receptor_name.nunique()>=12 and set(primary.transducer_family.unique())==set(PRIMARY_CLASSES): perf.append(run_geometry_subset(root,primary,features,'captured_is_primary_geometry',out,args.seed+17,args.quick))
    descript['sensitivity_models']=perf
    descript['annotation_completeness_warning']='Promiscuity counts are denominator-qualified. Counts across the three modeled families use only receptors with all three family flags known; four-family counts use only receptors with all four flags known. Unknown values are never converted to non-coupling.'
    write_json(descript,out/'A05_summary.json')
    # Missing receptor template for later curation, if needed.
    missing=sorted(set(bound.receptor_name.astype(str))-set(ann.loc[ann.n_known_families.gt(0),'receptor_name'].astype(str)))
    templ=pd.DataFrame({'receptor_name':missing,'Gs':np.nan,'Gi/o':np.nan,'Gq/11':np.nan,'G12/13':np.nan,'primary_family':'','annotation_source':'','evidence_note':''}); write_tsv(templ,out/'NEEDS_FUNCTIONAL_ANNOTATION.tsv')
    write_json(input_manifest(root,paths,['bound_frame','p0_matched','entity_audit'],'A05_functional_coupling_repertoire'),out/'manifest.json')
    print((out/'A05_summary.json').read_text())
if __name__=='__main__': main()
