#!/usr/bin/env python3
"""Create one independently auditable results archive; fail if a full run is incomplete."""
from pathlib import Path
import argparse,zipfile,hashlib,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--mode',choices=['full','smoke'],required=True);args=a.parse_args()
manifest=pd.read_csv(ROOT/'input/pdb_manifest.tsv',sep='\t')
if args.mode=='full':
    missing=[p for p in manifest.pdb_id if not (ROOT/'structures'/f'{p}.cif').is_file()]
    f=pd.read_csv(ROOT/'outputs/atomic_five/corrected_pair_comparisons.tsv',sep='\t')
    if missing or len(f)!=10:raise SystemExit(f'STOP: Not a complete run; missing={missing}, comparisons={len(f)}/10')
zipname=('GPCR_Ligand_Function_v5_COMPLETE_RESULTS.zip' if args.mode=='full' else 'GPCR_Ligand_Function_v5_OFFLINE_TEST_RESULTS.zip')
out=ROOT/zipname
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for sub in ['outputs','figures','reports','audit','input','scripts','structures']:
        for p in sorted((ROOT/sub).rglob('*')):
            if p.is_file():z.write(p,p.relative_to(ROOT))
    for f in ['README.md','requirements.txt','RUN_ALL.sh','RUN_OFFLINE_TEST.sh']:
        z.write(ROOT/f,f)
with zipfile.ZipFile(out) as z:
    if z.testzip() is not None:raise SystemExit('ZIP integrity test failed')
print('SUCCESS: '+str(out));print(f'Size: {out.stat().st_size:,} bytes')
