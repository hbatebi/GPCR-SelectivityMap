#!/usr/bin/env python3
"""Input gate: download missing Warren structures, validate all coordinates and chains."""
from pathlib import Path
import argparse, hashlib, json, subprocess
import pandas as pd
import gemmi
ROOT=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--download',action='store_true');a.add_argument('--require-all',action='store_true');args=a.parse_args()
meta=pd.read_csv(ROOT/'input/pdb_manifest.tsv',sep='\t')
missing=[]
for pdb in meta.pdb_id:
    path=ROOT/'structures'/f'{pdb}.cif'
    if (not path.is_file() or path.stat().st_size<50000) and args.download:
        print(f'Downloading deposited RCSB coordinate {pdb} ...',flush=True)
        tmp=path.with_suffix('.cif.tmp')
        subprocess.run(['curl','-fLSsS','--retry','3','--connect-timeout','30',
                        f'https://files.rcsb.org/download/{pdb}.cif','-o',str(tmp)],check=True)
        tmp.replace(path)
    if not path.is_file() or path.stat().st_size<50000:missing.append(pdb)
if missing and args.require_all:
    raise SystemExit('STOP: Missing coordinate files: '+', '.join(missing)+'. Check internet access to files.rcsb.org.')
# Mandatory expected ligand components for the small-molecule Warren series, fail if missed.
expected={row.pdb_id:row.expected_component_id for row in pd.read_csv(ROOT/'input/warren2024_five_ligand_manifest.tsv',sep='\t').itertuples(index=False)}
chains={'7X2C':('F','A'),'7X2D':('F','A'),'7X2F':('F','A')}
reports=[]
for pdb in meta.pdb_id:
    p=ROOT/'structures'/f'{pdb}.cif'
    if not p.exists():continue
    try:
        st=gemmi.read_structure(str(p));mo=st[0]
        rec,ga=chains.get(pdb,('R','A'))
        if mo.find_chain(rec) is None or mo.find_chain(ga) is None:
            raise ValueError(f'receptor/Gα author chain {rec}/{ga} missing')
        if pdb in expected:
            ligand_comp=expected[pdb]
            if not any(res.name==ligand_comp for ch in mo for res in ch):
                raise ValueError(f'expected ligand {ligand_comp} not present')
        reports.append(dict(pdb_id=pdb,size_bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),status='PASS'))
    except Exception as e:raise SystemExit(f'STOP: Bad structure {pdb}: {e}')
pd.DataFrame(reports).to_csv(ROOT/'audit/COORDINATE_SHA256_CURRENT.tsv',sep='\t',index=False)
(ROOT/'audit/coordinate_status.json').write_text(json.dumps({'available':len(reports),'expected':len(meta),'missing':missing,'complete':len(reports)==len(meta)},indent=2))
print(f'Validated {len(reports)}/{len(meta)} structures. Missing: {missing or "none"}')
