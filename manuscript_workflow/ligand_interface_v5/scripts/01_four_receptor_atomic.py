#!/usr/bin/env python3
"""Recalculate comparative geometry from verified deposited receptor/ligand chains.

Designed for an exploratory, same-receptor, same-family structural comparison.
Not an efficacy predictor; no kinetics or inferred causal mechanism. Gemmi + NumPy +
SciPy + pandas. Input PDBx/mmCIF structures are provided in structures/.

Explicit orthosteric assignments avoid the peptide-detection and fenoldopam-copy
errors in the earlier automatic extraction. Validate against each paper prior to
submission. All atom distances are from coordinates, without inference on dynamics.
"""
from pathlib import Path
import numpy as np, pandas as pd, gemmi
from scipy.spatial import cKDTree
import hashlib

ROOT=Path(__file__).resolve().parents[1]
CIFS=ROOT/'structures'; OUT=ROOT/'outputs'/'four_receptor_atomic';OUT.mkdir(exist_ok=True)
# These are deposited auth chains. Ligand assignments checked against published studies
# and deposited chain/residue names. 7X2C has 2 fenoldopam molecules in OBP and EBP.
STRUCT={
 '7E2Y':('R','A',[('component','R','SRO',501)],'serotonin'),
 '7E2Z':('R','A',[('component','R','9SC',501)],'aripiprazole'),
 '7X2C':('F','A',[('component','F','G3C',502),('component','F','G3C',503)],'fenoldopam (two copies)'),
 '7X2D':('F','A',[('component','F','86W',502)],'tavapadon'),
 '7X2F':('F','A',[('component','F','LDP',504)],'dopamine (+mevidalen; nonmatched control)'),
 '9CBL':('R','A',[('component','R','ALE',501)],'adrenaline'),
 '9CBM':('R','A',[('component','R','CZX',501)],'dexmedetomidine'),
 '7F4D':('R','A',[('peptide','M',None,None)],'alpha-MSH'),
 '7F4H':('R','A',[('peptide','F',None,None)],'afamelanotide'),
 '7F4I':('R','A',[('peptide','U',None,None)],'SHU9119')
}
PAIRS=[('5HT1A','7E2Y','7E2Z','primary'),('DRD1','7X2C','7X2D','primary'),
       ('ADRA2A','9CBL','9CBM','primary'),('MC1R','7F4D','7F4I','primary'),
       ('MC1R','7F4H','7F4I','sensitivity'),('MC1R','7F4D','7F4H','full_full_control'),
       ('DRD1','7X2F','7X2D','PAM_confounded_control')]
CORE=[(2.43,2.54),(3.42,3.48),(4.43,4.56)]
ALTERNATE=[(1.46,1.60),(2.39,2.54),(3.39,3.55),(4.38,4.56)]
CUT=4.5

MAPPING_CSV=ROOT/'input'/'generic_number_mapping_1360.csv'
mapdf=pd.read_csv(MAPPING_CSV,low_memory=False)

def h_atoms(r): return [a for a in r if a.element.name not in ('H','D')]
def xyz(a):return np.array([a.pos.x,a.pos.y,a.pos.z],dtype=float)
def coords(res):
 a=h_atoms(res)
 return np.array([xyz(x) for x in a]) if a else np.empty((0,3),float)
def allxyz(reses):
 aa=[coords(r) for r in reses]; aa=[a for a in aa if len(a)]
 return np.vstack(aa) if aa else np.empty((0,3),float)
def reskey(r):return (int(r.seqid.num),str(r.seqid.icode).strip())
def ca(r):
 return next((xyz(a) for a in r if a.name=='CA'),None)
def chain(s,name):
 x=s[0].find_chain(name)
 if x is None:raise ValueError('No chain '+name)
 return x

def mapped_pos(st,pdb,rchain):
 rc=chain(st,rchain);lookup={reskey(r):r for r in rc if r.seqid.num is not None};ans={}
 for z in mapdf.loc[mapdf.pdb_id==pdb].itertuples(index=False):
  try:gn=round(float(z.generic_number),2);number=int(z.author_residue_number)
  except (ValueError,TypeError):continue
  ic='' if pd.isna(z.author_insertion_code) else str(z.author_insertion_code).strip()
  if (number,ic) in lookup:ans[gn]=lookup[number,ic]
 return ans

def ligand_reses(st,pdb):
 ret=[]
 for kind,cname,rname,rnum in STRUCT[pdb][2]:
  cc=chain(st,cname)
  if kind=='peptide':
   rr=[r for r in cc if r.name not in ('HOH','WAT')]
   if not 4<=sum(ca(r) is not None for r in rr)<=30:raise ValueError(f'Not a short peptide {pdb} {cname}')
   ret+=rr
  else:
   match=[r for r in cc if r.name==rname and r.seqid.num==rnum]
   if len(match)!=1:raise ValueError(f'Component not unique {pdb} {cname} {rname}:{rnum}')
   ret+=match
 return ret

def receptor_residues(st,pdb):
 rec=chain(st,STRUCT[pdb][0]);return {reskey(r):r for r in rec if r.het_flag!='H' and ca(r) is not None}
def contact_set(resdict,target_xyz,cut=CUT):
 tree=cKDTree(target_xyz)
 out={}
 for key,res in resdict.items():
  pts=coords(res)
  if not len(pts):continue
  d=float(tree.query(pts,k=1)[0].min())
  if d<cut:out[key]=(res.name,d)
 return out

def alpha5(al):
 rs=[r for r in al if ca(r) is not None]
 return rs[-17:] if len(rs)>=17 else []
def alpha_pose(al):
 rs=alpha5(al)
 if len(rs)<17: return None,None
 pts=np.array([ca(r) for r in rs])
 axis=pts[-3:].mean(axis=0)-pts[:3].mean(axis=0)
 axis=axis/np.linalg.norm(axis)
 return axis,pts[-3:].mean(axis=0)

def align(d,e,windows):
 a,b=[],[]
 for gn in sorted(set(d['gp'])&set(e['gp'])):
  if not any(x<=gn<=y for x,y in windows):continue
  p,q=ca(d['gp'][gn]),ca(e['gp'][gn])
  if p is not None and q is not None:a.append(p);b.append(q)
 if len(a)<6:return None,None,len(a),float('nan')
 a,b=np.array(a),np.array(b)
 ac,bc=a.mean(axis=0),b.mean(axis=0)
 u,s,v=np.linalg.svd((b-bc).T@(a-ac));f=np.eye(3);f[-1,-1]=np.linalg.det(u@v)
 rot=u@f@v;shift=ac-bc@rot
 rmsd=float(np.sqrt(np.mean(np.sum(((b@rot+shift)-a)**2,axis=1))))
 return rot,shift,len(a),rmsd

def jacc(a,b):
 return float(len(a&b)/len(a|b)) if a|b else float('nan')

def resnamekeys(contact): return ';'.join(f'{v[0]}{k[0]}' for k,v in sorted(contact.items()))

def displacement_by_res(d,e,rot,shift,keys,include_backbone=False):
 moved=[]
 for k in sorted(keys):
  if k not in d['rec'] or k not in e['rec']:continue
  a={a.name:xyz(a) for a in h_atoms(d['rec'][k]) if include_backbone or a.name not in ('N','CA','C','O')}
  b={a.name:xyz(a) for a in h_atoms(e['rec'][k]) if include_backbone or a.name not in ('N','CA','C','O')}
  shared=sorted(a.keys()&b.keys())
  if len(shared)<2:continue
  dist=np.array([np.linalg.norm(a[t]-(b[t]@rot+shift)) for t in shared])
  moved.append((k,d['rec'][k].name,len(shared),float(np.sqrt(np.mean(dist**2)))))
 return moved

load={}; srows=[]; contactrows=[]
for pdb,meta in STRUCT.items():
 st=gemmi.read_structure(str(CIFS/(pdb+'.cif')))
 rec=receptor_residues(st,pdb);lig=ligand_reses(st,pdb);li=allxyz(lig)
 if len(li)==0:raise RuntimeError('No ligand atoms: '+pdb)
 al=chain(st,meta[1]);ac=allxyz([r for r in al if ca(r) is not None]);alpha5c=allxyz(alpha5(al))
 if len(alpha5c)==0:raise RuntimeError('No alpha5: '+pdb)
 lcont=contact_set(rec,li)
 gcont=contact_set(rec,ac)
 hcont=contact_set(rec,alpha5c)
 load[pdb]={'pdb':pdb,'rec':rec,'lig':lig,'gp':mapped_pos(st,pdb,meta[0]),'alpha':al,
            'lcont':lcont,'gcont':gcont,'hcont':hcont}
 for k,(name,d) in lcont.items():
  contactrows.append(dict(pdb_id=pdb,ligand_name=meta[3],author_residue_number=k[0],author_insertion_code=k[1],amino_acid=name,min_heavy_distance_A=d))
 srows.append(dict(pdb_id=pdb,ligand_name=meta[3],receptor_chain=meta[0],ga_chain=meta[1],
                   assigned_ligand_chains=';'.join(sorted(set(x[1] for x in meta[2]))),
                   ligand_component_identifiers=';'.join(str(x) for x in meta[2]),
                   n_ligand_heavy_atoms=len(li), n_receptor_ligand_contact_residues=len(lcont),
                   n_receptor_ga_contact_residues=len(gcont),n_receptor_alpha5_contact_residues=len(hcont),
                   receptor_ligand_contact_residues=resnamekeys(lcont),
                   n_mapped_generic_positions=len(load[pdb]['gp']),
                   ga_c_terminal_17_residues=f'{alpha5(al)[0].seqid.num}-{alpha5(al)[-1].seqid.num}',
                   cif_sha256=hashlib.sha256((CIFS/(pdb+'.cif')).read_bytes()).hexdigest()))

pd.DataFrame(srows).to_csv(OUT/'corrected_structure_metrics.tsv',sep='\t',index=False,float_format='%.9f')
pd.DataFrame(contactrows).to_csv(OUT/'all_receptor_ligand_contacts_4p5A.tsv',sep='\t',index=False,float_format='%.9f')

rows=[];move=[]
for label,full,part,role in PAIRS:
 d,e=load[full],load[part]
 rot,shift,n,core_rmsd=align(d,e,CORE)
 rot2,shift2,n2,core_rmsd2=align(d,e,ALTERNATE)
 if rot is None:raise RuntimeError(f'Cannot align {label}')
 first,last=alpha_pose(d['alpha']),alpha_pose(e['alpha'])
 ax1,tip1=first;ax2,tip2=last
 axis_angle=float(np.degrees(np.arccos(np.clip(np.dot(ax1,ax2@rot),-1,1))))
 tip=float(np.linalg.norm(tip1-(tip2@rot+shift)))
 tip2_val=float(np.linalg.norm(tip1-(tip2@rot2+shift2))) if rot2 is not None else float('nan')
 # Pocket: union of receptor-ligand 4.5 A contacts, all resolved receptor-side atoms.
 pocket=d['lcont'].keys()|e['lcont'].keys()
 intr=d['hcont'].keys()|e['hcont'].keys()
 allga=d['gcont'].keys()|e['gcont'].keys()
 side=displacement_by_res(d,e,rot,shift,pocket,include_backbone=False)
 alpha_side=displacement_by_res(d,e,rot,shift,intr,include_backbone=False)
 # Only residues mapped by generic numbering as robust canonical key in supplementary columns
 inv={reskey(r):p for p,r in d['gp'].items()}
 def namelist(keys):
  return ';'.join(f'{d["rec"].get(k,e["rec"].get(k)).name}{k[0]}' for k in sorted(keys))
 for k,name,atoms,disp in side:
  move.append(dict(case_id=label,comparison_role=role,full_pdb=full,partial_pdb=part,
                   residue_number=k[0],amino_acid=name, generic_number=inv.get(k,float('nan')),
                   common_sidechain_atoms=atoms,aligned_sidechain_rmsd_A=disp,
                   ligand_contact_full=k in d['lcont'],ligand_contact_partial=k in e['lcont']))
 rows.append(dict(case_id=label,role=role,first_pdb=full,second_pdb=part,
  first_ligand=STRUCT[full][3],second_ligand=STRUCT[part][3],
  n_receptor_core_alignment_CA=n,core_alignment_rmsd_A=core_rmsd,
  n_alt_alignment_CA=n2,alt_alignment_rmsd_A=core_rmsd2,
  alpha5_tip_shift_A=tip,alpha5_tip_shift_alt_alignment_A=tip2_val,
  alpha5_axis_rotation_deg=axis_angle,
  ligand_contact_residues_first=len(d['lcont']),ligand_contact_residues_second=len(e['lcont']),
  ligand_contact_jaccard=jacc(d['lcont'].keys(),e['lcont'].keys()),
  ligand_contact_shared_count=len(d['lcont'].keys() & e['lcont'].keys()),
  ligand_contact_first_only=namelist(d['lcont'].keys()-e['lcont'].keys()),
  ligand_contact_second_only=namelist(e['lcont'].keys()-d['lcont'].keys()),
  n_comparable_pocket_sidechains=len(side),
  mean_pocket_sidechain_rmsd_A=float(np.mean([v for _,_,_,v in side])) if side else float('nan'),
  median_pocket_sidechain_rmsd_A=float(np.median([v for _,_,_,v in side])) if side else float('nan'),
  top_pocket_motion=';'.join(f'{name}{k[0]}:{v:.2f}' for k,name,atoms,v in sorted(side,key=lambda x:-x[3])[:7]),
  ga_contact_residues_first=len(d['gcont']),ga_contact_residues_second=len(e['gcont']),
  ga_contact_jaccard=jacc(d['gcont'].keys(),e['gcont'].keys()),
  ga_contact_shared_count=len(d['gcont'].keys()&e['gcont'].keys()),
  ga_contact_first_only=namelist(d['gcont'].keys()-e['gcont'].keys()),
  ga_contact_second_only=namelist(e['gcont'].keys()-d['gcont'].keys()),
  ga_alpha5_contact_residues_first=len(d['hcont']),ga_alpha5_contact_residues_second=len(e['hcont']),
  ga_alpha5_contact_jaccard=jacc(d['hcont'].keys(),e['hcont'].keys()),
  n_comparable_ga_alpha5_interface_sidechains=len(alpha_side),
  mean_ga_alpha5_interface_sidechain_rmsd_A=float(np.mean([v for _,_,_,v in alpha_side])) if alpha_side else float('nan')))
pd.DataFrame(rows).to_csv(OUT/'corrected_pair_comparisons.tsv',sep='\t',index=False,float_format='%.9f')
pd.DataFrame(move).to_csv(OUT/'pocket_sidechain_motion.tsv',sep='\t',index=False,float_format='%.9f')
print(pd.DataFrame(rows)[['case_id','role','first_pdb','second_pdb','n_receptor_core_alignment_CA','alpha5_tip_shift_A','alpha5_axis_rotation_deg','ligand_contact_residues_first','ligand_contact_residues_second','ligand_contact_jaccard','ga_contact_jaccard','mean_pocket_sidechain_rmsd_A']].to_string(index=False))
print('\nCreated:',OUT)
