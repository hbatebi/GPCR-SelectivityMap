#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd


def loadj(p):
    try: return json.loads(p.read_text())
    except Exception as e: return {'_read_error':str(e)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--expected-permutations',type=int,default=None); args=ap.parse_args()
    root=args.project_root.resolve(); out=root/'revision_bjp_2026'; out.mkdir(exist_ok=True)
    required=[
      '00_preflight/preflight.json','A01_expanded_matched_set/A01_summary.json','A02_low_dimensional/A02_summary.tsv',
      'A03_archive_trained_scoring/A03_summary.tsv','A04_structural_descriptives/P0_strict_feature_summary.tsv',
      'A05_functional_coupling_repertoire/A05_summary.json','A06_g12_g13_census/A06_summary.json',
      'A07_deformation_diagnostics/P0_strict_summary.json','A08_confound_audit/P0_strict_summary.json',
      'A09_matched_permutation_nulls/P0_strict/summary.json','A10_partner_free_relatedness/P0_strict_summary.json',
      'R14_same_receptor_408/all_408_feature_by_class_pair_comparisons.tsv','D01_feature_dictionary/matched_92_feature_dictionary.tsv',
      'D02_fold_inventory/fold_inventory.tsv','D03_matched_pdb_tables/P0_strict_matched_structure_table.tsv']
    rows=[]; ok=True
    for rel in required:
        p=out/rel; exists=p.exists() and p.stat().st_size>0; rows.append({'item':rel,'exists_nonempty':exists}); ok &= exists
    pre=loadj(out/'00_preflight/preflight.json') if (out/'00_preflight/preflight.json').exists() else {}
    frozen=bool(pre.get('frozen_result_matches_expected',False)); ok &= frozen
    a1=loadj(out/'A01_expanded_matched_set/A01_summary.json') if (out/'A01_expanded_matched_set/A01_summary.json').exists() else {}
    p1=(a1.get('populations') or {}).get('P1_terminal_strict',{})
    a5=loadj(out/'A05_functional_coupling_repertoire/A05_summary.json') if (out/'A05_functional_coupling_repertoire/A05_summary.json').exists() else {}
    a9=loadj(out/'A09_matched_permutation_nulls/P0_strict/summary.json') if (out/'A09_matched_permutation_nulls/P0_strict/summary.json').exists() else {}
    if args.expected_permutations is not None:
        perm_ok=int(a9.get('iterations',-1))>=args.expected_permutations; ok &= perm_ok
    else: perm_ok=None
    index={
      'overall_file_validation_passed':bool(ok),
      'frozen_P0_reproduced':frozen,
      'P0_positive_receptors':pre.get('p0_reproduced_positive_receptors'),
      'P0_mean_margin_gain':pre.get('p0_reproduced_mean_margin_gain'),
      'P1_status':p1.get('status'),'P1_n_receptors':p1.get('n_receptors'),'P1_n_units':p1.get('n_units'),'P1_positive_receptors':p1.get('n_positive'),'P1_mean_margin_gain':p1.get('mean_margin_gain'),
      'A05_bound_receptors':a5.get('n_bound_receptors'),'A05_receptors_complete_three_classes':a5.get('n_receptors_complete_for_three_analyzed_classes'),
      'A09_P0_iterations':a9.get('iterations'),'A09_expected_permutations_passed':perm_ok,
      'manual_fusion_topology_rows':None,
    }
    m=out/'A01_expanded_matched_set/fusion_topology_manual_review.tsv'
    if m.exists():
        try: index['manual_fusion_topology_rows']=len(pd.read_csv(m,sep='\t'))
        except: pass
    pd.DataFrame(rows).to_csv(out/'VALIDATION_CHECKLIST.tsv',sep='\t',index=False)
    (out/'REVISION_ANALYSIS_INDEX.json').write_text(json.dumps(index,indent=2,sort_keys=True,default=str))
    print(json.dumps(index,indent=2,sort_keys=True,default=str))
    if not ok: raise SystemExit('Revision output validation failed; inspect VALIDATION_CHECKLIST.tsv and logs.')

if __name__=='__main__': main()
