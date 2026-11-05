#!/usr/bin/env python3
"""Fail-fast checks against independently audited original baseline coordinates."""
from pathlib import Path
import argparse,json,math,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--mode',choices=['full','smoke'],required=True);args=a.parse_args()
f=pd.read_csv(ROOT/'outputs/four_receptor_atomic/corrected_pair_comparisons.tsv',sep='\t')
ref={'5HT1A':(0.5789473684,0.7083333333,0.841914),
     'DRD1':(0.8076923077,0.8055555556,0.701018),
     'ADRA2A':(0.9230769231,0.7894736842,1.074543),
     'MC1R':(0.6470588235,0.8709677419,0.414942)}
for label,(lig,ga,tip) in ref.items():
    r=f[(f.case_id==label)&f.role.eq('primary')]
    if len(r)!=1:raise SystemExit(f'REGRESSION FAILURE: unexpected {label} rows')
    row=r.iloc[0]
    for nm,exp in [('ligand_contact_jaccard',lig),('ga_contact_jaccard',ga),('alpha5_tip_shift_A',tip)]:
        if abs(float(row[nm])-exp)>3e-5:raise SystemExit(f'REGRESSION FAILURE: {label} {nm} = {row[nm]}, expected {exp}')
q=pd.read_csv(ROOT/'outputs/atomic_five/corrected_pair_comparisons.tsv',sep='\t')
r=q[((q.first_pdb=='8FYL')&(q.second_pdb=='8FYX')) | ((q.first_pdb=='8FYX')&(q.second_pdb=='8FYL'))]
if len(r)!=1:raise SystemExit('REGRESSION FAILURE: Vilazodone/buspirone pair not unique')
row=r.iloc[0]
if abs(row.alpha5_tip_shift_A-0.497222)>3e-5 or abs(row.ligand_contact_jaccard-0.5)>1e-6 or abs(row.ga_contact_jaccard-0.952381)>3e-5:
    raise SystemExit('REGRESSION FAILURE: Vilazodone/buspirone original 2024 comparison not reproduced')
if a.parse_args().mode=='full' and len(q)!=10:raise SystemExit('REGRESSION FAILURE: incomplete Warren series')
status={'baseline_four_primary':'PASS','baseline_Warren_8FYL_8FYX':'PASS','Warren_pairs':len(q),'mode':a.parse_args().mode}
(ROOT/'audit/REPRODUCIBILITY_TEST.json').write_text(json.dumps(status,indent=2))
print('Reproducibility baseline PASS:',status)
