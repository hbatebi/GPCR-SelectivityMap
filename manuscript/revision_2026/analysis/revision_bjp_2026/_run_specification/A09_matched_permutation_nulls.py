#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from bjp_revision_common import *


def analyse(root,label,paired,features,observed_performance,observed_predictions,out,iterations,seed,quick):
    ProjectConfig,_analysis_settings,ExperimentSpec,FeatureSpec,*_=ensure_project_imports(root)
    from gpcr_icl2.microswitches.publication_controls import grouped_permutation_test
    config=ProjectConfig.from_dir(root/'config'); settings=_analysis_settings(config,quick)
    # The submitted full transducer-imprint run explicitly used saga in full mode.
    if not quick:
        settings['solver']='saga'; settings['max_iter']=max(150000,int(settings.get('max_iter',50000))); settings['tolerance']=min(1e-4,float(settings.get('tolerance',1e-4)))
    ex=[
        ExperimentSpec('bound_geometry',FeatureSpec('bound_geometry',tuple('bound__'+f for f in features)),nested_tuning=False),
        ExperimentSpec('receptor_only_geometry',FeatureSpec('receptor_only_geometry',tuple('receptor_only__'+f for f in features)),nested_tuning=False),
    ]
    obsperf=observed_performance[observed_performance.experiment.isin(['bound_geometry','receptor_only_geometry'])].copy()
    nperm=min(iterations,20) if quick else iterations
    # Each permutation refits a fixed-parameter model; observed performance comes from
    # the already-generated P0/P1 model run so A09 does not repeat expensive nested CV.
    null,summary=grouped_permutation_test(paired,ex,obsperf,settings,nperm,seed+900)
    d=out/label; d.mkdir(parents=True,exist_ok=True)
    write_tsv(obsperf,d/'observed_performance.tsv'); write_tsv(summary,d/'permutation_model_summary.tsv'); write_tsv(null,d/'permutation_null.tsv')
    if not null.empty:
        piv=null.pivot(index='iteration',columns='experiment',values='macro_roc_auc').dropna()
        piv['delta_auc_bound_minus_free']=piv['bound_geometry']-piv['receptor_only_geometry']; write_tsv(piv.reset_index(),d/'paired_delta_auc_null.tsv')
        obs=paired_auc_from_predictions(observed_predictions,'bound_geometry','receptor_only_geometry'); delta=obs['delta_auc']
        pv=float((1+(piv.delta_auc_bound_minus_free>=delta).sum())/(1+len(piv))) if np.isfinite(delta) else np.nan
        pv2=float((1+(piv.delta_auc_bound_minus_free.abs()>=abs(delta)).sum())/(1+len(piv))) if np.isfinite(delta) else np.nan
        result={'population':label,**obs,'paired_delta_permutation_p_upper':pv,'paired_delta_permutation_p_two_sided_abs':pv2,'null_delta_mean':float(piv.delta_auc_bound_minus_free.mean()),'iterations':len(piv)}
    else:
        result={'population':label,'status':'no_null'}
    write_json(result,d/'summary.json'); return result


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/A09_matched_permutation_nulls')); ap.add_argument('--iterations',type=int,default=1000); ap.add_argument('--seed',type=int,default=20291729); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True)
    paths=discover_inputs(root); require_paths(paths,['p0_matched','p0_predictions','p0_performance']); feats=p0_base_features(paths); rows=[]
    rows.append(analyse(root,'P0_strict',pd.read_csv(paths['p0_matched']),feats,pd.read_csv(paths['p0_performance']),pd.read_csv(paths['p0_predictions']),out,args.iterations,args.seed,args.quick))
    p1dir=root/'revision_bjp_2026/A01_expanded_matched_set/P1_terminal_strict'
    if (p1dir/'paired_units.tsv').exists() and (p1dir/'performance.tsv').exists() and (p1dir/'predictions.tsv').exists():
        rows.append(analyse(root,'P1_terminal_strict',pd.read_csv(p1dir/'paired_units.tsv',sep='\t'),feats,pd.read_csv(p1dir/'performance.tsv',sep='\t'),pd.read_csv(p1dir/'predictions.tsv',sep='\t'),out,args.iterations,args.seed+1000,args.quick))
    write_tsv(pd.DataFrame(rows),out/'A09_summary.tsv'); write_json(input_manifest(root,paths,['p0_matched','p0_predictions','p0_performance'],'A09_matched_permutation_nulls'),out/'manifest.json')

if __name__=='__main__': main()
