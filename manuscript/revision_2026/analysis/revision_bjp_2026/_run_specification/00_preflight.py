#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys
from pathlib import Path
import pandas as pd
from bjp_revision_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',type=Path,default=Path.cwd()); ap.add_argument('--outdir',type=Path,default=Path('revision_bjp_2026/00_preflight')); args=ap.parse_args()
    root=args.project_root.resolve(); out=args.outdir if args.outdir.is_absolute() else root/args.outdir; out.mkdir(parents=True,exist_ok=True)
    paths=discover_inputs(root); required=['bound_frame','active_none','strict_frame','entity_audit','active_transducers','same_receptor','p0_matched','p0_predictions']; require_paths(paths,required)
    p0=pd.read_csv(paths['p0_matched']); pred=pd.read_csv(paths['p0_predictions']); feats=p0_base_features(paths); rec,s=receptor_margin_table(pred)
    checks={
      'python':sys.version,
      'p0_units':len(p0),'p0_receptors':p0.receptor_name.nunique(),'p0_features':len(feats),
      'p0_reproduced_positive_receptors':s.get('n_positive'),'p0_reproduced_mean_margin_gain':s.get('mean_margin_gain'),'p0_reproduced_signflip_p':s.get('signflip_p_two_sided'),
      'strict_partner_free_structures':pd.read_csv(paths['strict_frame']).pdb_id.nunique(),
      'bound_units':len(pd.read_csv(paths['bound_frame'])),'same_receptor_comparisons':len(pd.read_csv(paths['same_receptor'])),
      'structure_cif_count':len(list((root/'data/structures').glob('*.cif'))) if (root/'data/structures').exists() else 0,
    }
    checks['frozen_result_matches_expected']=bool(checks['p0_units']==15 and checks['p0_receptors']==12 and checks['p0_reproduced_positive_receptors']==11 and abs(checks['p0_reproduced_mean_margin_gain']-0.3887713288847121)<1e-8)
    write_tsv(rec,out/'P0_receptor_margin_reproduction.tsv'); write_json(checks,out/'preflight.json'); write_json(input_manifest(root,paths,required,'00_preflight'),out/'manifest.json')
    print(pd.Series(checks).to_string())
    if not checks['frozen_result_matches_expected']: raise SystemExit('Preflight failed: frozen P0 result did not reproduce exactly.')
if __name__=='__main__': main()
