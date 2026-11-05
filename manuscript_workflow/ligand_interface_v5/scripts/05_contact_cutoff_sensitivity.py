#!/usr/bin/env python3
"""Sensitivity of residue contact overlap to heavy-atom cutoff and ligand size.

Residue-contact overlap is *not* an efficacy measure and contact-region Jaccards
should not be directly equated across regions with very different residue counts.
"""
from pathlib import Path
import ast, argparse, itertools, numpy as np, pandas as pd, gemmi
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs';OUT.mkdir(exist_ok=True)
a=argparse.ArgumentParser();a.add_argument('--available-only',action='store_true');args=a.parse_args()
def literals(path,attribute):
    tree=ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==attribute for t in node.targets):
            try:return ast.literal_eval(node.value)
            except ValueError:raise SystemExit('Unable to parse '+attribute)
    raise SystemExit('Missing '+attribute)
FOUR=literals(ROOT/'scripts/01_four_receptor_atomic.py','STRUCT')
FIVE=literals(ROOT/'scripts/02_five_ligand_atomic.py','STRUCT')
FOUR_PAIRS=literals(ROOT/'scripts/01_four_receptor_atomic.py','PAIRS')
STRUCT={**FOUR,**FIVE}
PAIRS=[('FOUR_'+name,p,q,role) for name,p,q,role in FOUR_PAIRS]
DNG={'8FY8','8FYX'}
PAIRS += [('WARREN2024_5HT1A',a,b,'same_Ga_construct' if ((a in DNG)==(b in DNG)) else 'mixed_Ga_construct') for a,b in itertools.combinations(sorted(FIVE),2)]
def heavy(res):return np.array([(at.pos.x,at.pos.y,at.pos.z) for at in res if at.element.name not in ('H','D')],dtype=float).reshape(-1,3)
def key(r):return (int(r.seqid.num),str(r.seqid.icode).strip())
load={}
for pdb,(rc,gc,ligns,name) in STRUCT.items():
    p=ROOT/'structures'/f'{pdb}.cif'
    if not p.is_file():
        if args.available_only:continue
        raise SystemExit(f'Missing {pdb}.cif, run bash RUN_ALL.sh to download')
    st=gemmi.read_structure(str(p))[0];rec=st.find_chain(rc);ga=st.find_chain(gc)
    res={key(r):r for r in rec if r.het_flag!='H' and any(a.name=='CA' for a in r)}
    ligand=[]
    for kind,ch,name_comp,num in ligns:
        ligand_chain=st.find_chain(ch)
        if kind=='peptide':ligand += [r for r in ligand_chain if r.name not in ('HOH','WAT')]
        else:
            match=[r for r in ligand_chain if r.name==name_comp and (num is None or r.seqid.num==num)]
            if len(match)!=1:raise SystemExit(f'{pdb}: expected one component {name_comp}, found {len(match)}')
            ligand+=match
    coords_lig=np.vstack([heavy(r) for r in ligand if len(heavy(r))]);
    coords_ga=np.vstack([heavy(r) for r in ga if any(a.name=='CA' for a in r) and len(heavy(r))]);
    if not len(coords_lig) or not len(coords_ga):raise SystemExit('no target heavy atoms for '+pdb)
    tr_l, tr_g=cKDTree(coords_lig),cKDTree(coords_ga)
    rec_dist={k:(float(tr_l.query(heavy(r))[0].min()),float(tr_g.query(heavy(r))[0].min())) for k,r in res.items() if len(heavy(r))}
    load[pdb]=dict(dist=rec_dist,n_ligand_atoms=len(coords_lig),n_ga_atoms=len(coords_ga))
def scores(a,b):
    inter=len(a&b);n1=len(a);n2=len(b)
    union=len(a|b)
    return inter/union if union else float('nan'),2*inter/(n1+n2) if n1+n2 else float('nan'), inter/min(n1,n2) if min(n1,n2) else float('nan')
rows=[]
for name,p,q,role in PAIRS:
    if p not in load or q not in load:continue
    for cut in [4.0,4.5,5.0]:
        Lp={k for k,(dl,dg) in load[p]['dist'].items() if dl<cut}
        Lq={k for k,(dl,dg) in load[q]['dist'].items() if dl<cut}
        Gp={k for k,(dl,dg) in load[p]['dist'].items() if dg<cut}
        Gq={k for k,(dl,dg) in load[q]['dist'].items() if dg<cut}
        ls=scores(Lp,Lq);gs=scores(Gp,Gq)
        rows.append(dict(group=name,pdb_a=p,pdb_b=q,comparison_role=role,cutoff_A=cut,
            n_ligand_atoms_a=load[p]['n_ligand_atoms'],n_ligand_atoms_b=load[q]['n_ligand_atoms'],
            ligand_contact_count_a=len(Lp),ligand_contact_count_b=len(Lq),
            ga_contact_count_a=len(Gp),ga_contact_count_b=len(Gq),
            ligand_jaccard=ls[0],ligand_dice=ls[1],ligand_overlap_coefficient=ls[2],
            ga_jaccard=gs[0],ga_dice=gs[1],ga_overlap_coefficient=gs[2],
            ga_minus_ligand_jaccard=gs[0]-ls[0]))
res=pd.DataFrame(rows);res.to_csv(OUT/'contact_cutoff_and_ligand_size_sensitivity.tsv',sep='\t',index=False,float_format='%.8f')
print(f'Cutoff sensitivity: {len(rows)} comparisons across {len(rows)//3} PDB pairs. Available PDBs {len(load)}.')
