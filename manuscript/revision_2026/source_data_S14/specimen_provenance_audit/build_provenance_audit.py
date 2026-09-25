#!/usr/bin/env python3
"""Build provenance-aware partner-free audits for the BJP revision.

This script does NOT fit any predictive model. It reconstructs the strict P0,
terminal-fusion P1, and ligand-permissive inventories and applies a frozen
provenance adjudication table. The purpose is to determine whether the
same-receptor three-class matched analysis remains estimable after excluding
receptor-focused/local reconstructions obtained from transducer-containing
experimental specimens and correcting state misannotations.
"""
from __future__ import annotations
import argparse, json, hashlib
from pathlib import Path
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# Frozen manual adjudications, based on deposited mmCIF metadata and cited primary/PDB sources.
# Anything not listed defaults to eligible only if it is an X-ray strict entry with no contrary evidence;
# for EM entries, absence from this table is an error to prevent silent acceptance.
DECISIONS = {
    # Strict X-ray genuine transducer-free comparators
    '2YDO': ('eligible_genuine_partner_free', 'active_or_active_like', 'Thermostabilised agonist-bound A2A receptor crystal; no intracellular transducer in specimen/model.'),
    '2YDV': ('eligible_genuine_partner_free', 'active_or_active_like', 'Thermostabilised agonist-bound A2A receptor crystal; no intracellular transducer in specimen/model.'),
    '4UG2': ('eligible_genuine_partner_free', 'active_or_active_like', 'Thermostabilised agonist-bound A2A receptor crystal; no intracellular transducer in specimen/model.'),
    '4UHR': ('eligible_genuine_partner_free', 'active_or_active_like', 'Thermostabilised agonist-bound A2A receptor crystal; no intracellular transducer in specimen/model.'),
    '3CAP': ('eligible_genuine_partner_free', 'active_or_active_like', 'Native ligand-free opsin crystal; no intracellular transducer.'),
    '3PXO': ('eligible_genuine_partner_free', 'active_or_active_like', 'Metarhodopsin II-alone crystal; related peptide-bound structure is a separate PDB entry.'),
    '4BEZ': ('eligible_genuine_partner_free', 'active_or_active_like', 'G90D rhodopsin active-conformation crystal; no intracellular transducer.'),
    '5DYS': ('eligible_genuine_partner_free', 'active_or_active_like', 'T94I rhodopsin crystal; no intracellular transducer.'),
    '5TE3': ('eligible_genuine_partner_free', 'active_or_active_like', 'Opsin crystal; no intracellular transducer.'),
    '6FK6': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FK7': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FK8': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FK9': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FKA': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FKB': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FKC': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6FKD': ('eligible_genuine_partner_free', 'active_or_active_like', 'Pharmacologically stabilised opsin crystal; small-molecule ligand, no intracellular transducer.'),
    '6PT3': ('eligible_genuine_partner_free', 'active', 'Active DOP crystal with small-molecule agonist DPI-287; no intracellular transducer.'),

    # Strict EM entries that are receptor-focused/local products of transducer-containing specimens
    '7VV6': ('exclude_transducer_conditioned_focus', 'active', 'MRGPRX2 receptor local reconstruction from a receptor-complex cryo-EM study; not an independent partner-free specimen.'),
    '7YU4': ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly identifies LPA1-Gi complex focused on receptor.'),
    '8EA0': ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly identifies miniGq-coupled M3 receptor local refinement.'),
    '8IWE': ('exclude_transducer_conditioned_focus', 'active', 'Receptor-focused TAAR9 reconstruction; particle count matches SPE-mTAAR9-Gs full-complex entry 8IW4.'),
    '8IWM': ('exclude_transducer_conditioned_focus', 'active', 'Receptor-focused TAAR9 reconstruction; particle count matches PEA-mTAAR9-Gs full-complex entry 8IW7.'),
    '8K5B': ('exclude_transducer_conditioned_focus', 'active', 'PDB title says local refinement; source paper reports three ligand-bound HCAR2-Gi signaling complexes.'),
    '8K5C': ('exclude_transducer_conditioned_focus', 'active', 'Local HCAR2 refinement; particle count matches full HCAR2-Gi entry 8I7V.'),
    '8K5D': ('exclude_transducer_conditioned_focus', 'active', 'Local HCAR2 refinement; particle count matches full HCAR2-Gi entry 8I7W.'),
    '8W8S': ('exclude_transducer_conditioned_focus', 'active', 'Locally refined AA14-GPR101 receptor model; particle count matches AA14-GPR101-Gs entry 8W8R.'),
    '8WCC': ('exclude_transducer_conditioned_focus', 'active', 'CHA-mTAAR1 receptor reconstruction has the same 296,951 particles as CHA-mTAAR1-Gq entry 8WCB.'),
    '8XBG': ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly identifies GPR34-Gi complex receptor-focused; particle count matches full complex 8XBE.'),
    '8XBI': ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly identifies GPR34-Gi complex receptor-focused; particle count matches full complex 8XBH.'),
    '8IYX': ('exclude_inactive_state', 'inactive', 'Antagonist YL-365-bound GPR34 structure; primary publication explicitly describes it as inactive.'),

    # β2AR transition series: all 20 are locally refined receptor maps from Gs-containing particles
    **{pid: ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly states locally refined receptor from beta2AR-Gs complex transition intermediate.')
       for pid in ['8GGI','8GGJ','8GGK','8GGL','8GGM','8GGN','8GGO','8GGP','8GGQ','8GGR','8GGS','8GGT','8GGU','8GGV','8GGW','8GGX','8GGY','8GGZ','8GH0','8GH1']},

    # P1 terminal-fusion additions
    '8UGW': ('eligible_partner_free_but_terminal_fusion', 'active_or_active_like', 'A2A X-ray monomer with engineered endolysin fusion; no transducer; designed highly signaling-active receptor.'),
    '7ZL9': ('exclude_inactive_state', 'inactive', 'HCA2-b562 X-ray structure from a study explicitly describing the crystal structures as inactive.'),
    '8IHK': ('exclude_transducer_conditioned_focus', 'active', 'HCA3-Gi complex local reconstruction; particle count matches full HCA3-Gi entry 8IHJ.'),
    '8KEX': ('exclude_transducer_conditioned_focus', 'active', 'Gq-coupled MRGPRX4 local reconstruction; particle count matches full complex 8K4S.'),
    '6AK3': ('eligible_partner_free_but_terminal_fusion', 'active', 'Agonist-bound EP3 X-ray receptor with b562 terminal fusion; no G protein/transducer.'),

    # P1 terminal + extracellular-ligand additions
    '6LW5': ('eligible_partner_free_but_terminal_fusion', 'active', 'Active-conformation FPR2 X-ray structure with b562 fusion and extracellular peptide agonist; no intracellular transducer.'),
    '4XEE': ('eligible_partner_free_but_terminal_fusion', 'active_like', 'Active-like rat NTSR1 X-ray structure with endolysin fusion and extracellular neurotensin peptide; no intracellular transducer.'),
    '4XES': ('eligible_partner_free_but_terminal_fusion', 'active_like', 'Active-like rat NTSR1 X-ray structure with endolysin fusion and extracellular neurotensin peptide; no intracellular transducer.'),

    # Ligand-permissive additions
    '4XT3': ('eligible_genuine_partner_free_ligand_permissive', 'active_like', 'US28 crystal bound to extracellular chemokine CX3CL1; active-state-like receptor, no intracellular transducer.'),
    '6PT2': ('eligible_genuine_partner_free_ligand_permissive', 'active', 'Active DOP crystal with extracellular peptide agonist KGCHM07; no intracellular transducer.'),
    '7VV0': ('exclude_transducer_conditioned_focus', 'active', 'MRGPRX2 local receptor reconstruction from receptor-complex cryo-EM study.'),
    '7VV4': ('exclude_transducer_conditioned_focus', 'active', 'MRGPRX2 local receptor reconstruction; particle count matches full complex entry 7VV3.'),
    '8IY6': ('exclude_transducer_conditioned_focus', 'active', 'PDB/EMDB sample is Endothelin-1-ETB-Gi complex, focused on receptor.'),
    '8JPF': ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly says focused refinement of NTSR1 in NTSR1-GRK2-Galpha(q) complexes.'),
    '8ZRT': ('exclude_transducer_conditioned_focus', 'active', 'PDB title explicitly says receptor-focused ET-1-bound ETBR-DNGI complex.'),
}

EVIDENCE = {
 '2YDO': 'https://www.rcsb.org/structure/2YDO', '2YDV': 'https://www.rcsb.org/structure/2YDV',
 '4UG2': 'https://www.rcsb.org/structure/4UG2', '4UHR': 'https://www.rcsb.org/structure/4UHR',
 '3CAP': 'https://www.rcsb.org/structure/3CAP', '3PXO': 'https://www.rcsb.org/structure/3PXO',
 '4BEZ': 'https://www.rcsb.org/structure/4BEZ', '5DYS': 'https://www.rcsb.org/structure/5DYS',
 '5TE3': 'https://www.rcsb.org/structure/5TE3', '6PT3': 'https://www.rcsb.org/structure/6PT3',
 '6PT2': 'https://www.rcsb.org/structure/6PT2', '4XT3': 'https://www.rcsb.org/structure/4XT3',
 '7YU4': 'https://www.rcsb.org/structure/7YU4', '8EA0': 'https://www.rcsb.org/structure/8EA0',
 '8IYX': 'https://doi.org/10.1073/pnas.2308435120', '8K5B': 'https://doi.org/10.1038/s41467-023-42764-8',
 '8K5C': 'https://doi.org/10.1038/s41467-023-42764-8', '8K5D': 'https://doi.org/10.1038/s41467-023-42764-8',
 '8W8S': 'https://doi.org/10.1038/s41589-023-01456-6', '8WCC': 'https://www.rcsb.org/structure/8WCC',
 '8IWE': 'https://doi.org/10.1038/s41586-023-06106-4', '8IWM': 'https://doi.org/10.1038/s41586-023-06106-4',
 '8XBG': 'https://www.rcsb.org/structure/8XBG', '8XBI': 'https://www.rcsb.org/structure/8XBI',
 '7VV6': 'https://www.rcsb.org/structure/7VV6', '7VV0': 'https://www.rcsb.org/structure/7VV0', '7VV4': 'https://www.rcsb.org/structure/7VV4',
 '8IY6': 'https://www.rcsb.org/structure/8IY6', '8JPF': 'https://www.rcsb.org/structure/8JPF', '8ZRT': 'https://www.rcsb.org/structure/8ZRT',
 '7ZL9': 'https://doi.org/10.1038/s41467-023-37177-6', '8IHK': 'https://www.rcsb.org/structure/8IHK',
 '8KEX': 'https://www.rcsb.org/structure/8KEX', '8UGW': 'https://www.rcsb.org/structure/8UGW', '6AK3': 'https://www.rcsb.org/structure/6AK3',
 '6LW5': 'https://www.rcsb.org/structure/6LW5', '4XEE': 'https://www.rcsb.org/structure/4XEE', '4XES': 'https://www.rcsb.org/structure/4XES',
}
for pid in ['6FK6','6FK7','6FK8','6FK9','6FKA','6FKB','6FKC','6FKD']:
    EVIDENCE[pid] = f'https://www.rcsb.org/structure/{pid}'
for pid in ['8GGI','8GGJ','8GGK','8GGL','8GGM','8GGN','8GGO','8GGP','8GGQ','8GGR','8GGS','8GGT','8GGU','8GGV','8GGW','8GGX','8GGY','8GGZ','8GH0','8GH1']:
    EVIDENCE[pid] = f'https://www.rcsb.org/structure/{pid}'

SHARED_PARTICLE_PAIRS = [
    ('8W8S','8W8R','GPR101 receptor-focused vs Gs complex'),
    ('8WCC','8WCB','TAAR1 receptor-focused vs Gq complex'),
    ('8IWE','8IW4','TAAR9 receptor-focused vs Gs complex'),
    ('8IWM','8IW7','TAAR9 receptor-focused vs Gs complex'),
    ('8K5C','8I7V','HCAR2 local vs Gi complex'),
    ('8K5D','8I7W','HCAR2 local vs Gi complex'),
    ('8XBG','8XBE','GPR34 receptor-focused vs Gi complex'),
    ('8XBI','8XBH','GPR34 receptor-focused vs Gi complex'),
    ('8IHK','8IHJ','HCAR3 local vs Gi complex'),
    ('8KEX','8K4S','MRGPRX4 local vs Gq complex'),
    ('7VV4','7VV3','MRGPRX2 local vs Gi complex'),
]
PAIR_MAP = {a:(b,n) for a,b,n in SHARED_PARTICLE_PAIRS}

def sval(d, key):
    x=d.get(key, '')
    if isinstance(x, list): return '; '.join(str(i) for i in x)
    return str(x)

def cif_meta(base: Path, pid: str):
    f=base/'data'/'structures'/f'{pid.upper()}.cif'
    if not f.exists():
        return dict(cif_title='', emdb_ids='', em_sample_name='', em_particles='', citation_title='', citation_doi='')
    d=MMCIF2Dict(str(f))
    dbid=d.get('_database_2.database_id',[]); code=d.get('_database_2.database_code',[])
    if not isinstance(dbid,list): dbid=[dbid]
    if not isinstance(code,list): code=[code]
    emdb=';'.join(c for a,c in zip(dbid,code) if str(a).upper()=='EMDB')
    return dict(
        cif_title=sval(d,'_struct.title'), emdb_ids=emdb,
        em_sample_name=sval(d,'_em_entity_assembly.name'), em_particles=sval(d,'_em_3d_reconstruction.num_particles'),
        citation_title=sval(d,'_citation.title'), citation_doi=sval(d,'_citation.pdbx_database_id_DOI'))

def audit_frame(df, base, population):
    rows=[]
    for _,r in df.iterrows():
        pid=str(r['pdb_id']).upper()
        if pid not in DECISIONS:
            raise RuntimeError(f'No frozen provenance decision for {pid} ({population}); explicit adjudication required.')
        status,state,reason=DECISIONS[pid]
        meta=cif_meta(base,pid)
        eligible=status.startswith('eligible_')
        row={
            'population':population,'pdb_id':pid,'receptor_name':r.get('receptor_name',''),
            'chain_id':r.get('chain_id',r.get('preferred_chain','')),
            'experimental_method':r.get('experimental_method',''),'resolution':r.get('resolution',''),
            'original_activation_state':r.get('activation_state',''),'original_ligand_type':r.get('ligand_type',''),
            'fusion_protein':r.get('fusion_protein',''),'fusion_names':r.get('fusion_names',''),
            **meta,
            'provenance_status':status,'state_status':state,'final_eligibility':bool(eligible),
            'reason':reason,'evidence_url':EVIDENCE.get(pid,f'https://www.rcsb.org/structure/{pid}'),
            'local_evidence':'Deposited mmCIF _struct.title/_em_entity_assembly.name/_em_3d_reconstruction.num_particles where available'
        }
        if pid in PAIR_MAP:
            full_pid, pair_note = PAIR_MAP[pid]
            full_meta = cif_meta(base, full_pid)
            row.update({
                'paired_full_complex_pdb': full_pid,
                'paired_full_complex_particles': full_meta['em_particles'],
                'exact_particle_count_match_to_full_complex': bool(meta['em_particles'] and meta['em_particles']==full_meta['em_particles']),
                'paired_full_complex_title': full_meta['cif_title'],
                'pair_evidence_note': pair_note,
            })
        else:
            row.update({'paired_full_complex_pdb':'','paired_full_complex_particles':'',
                        'exact_particle_count_match_to_full_complex':'',
                        'paired_full_complex_title':'','pair_evidence_note':''})
        rows.append(row)
    return pd.DataFrame(rows)

def sha256(p):
    h=hashlib.sha256();
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--bundle-root', default='/mnt/data/provenance_audit/revision_bundle')
    ap.add_argument('--revision-root', default='/mnt/data/provenance_audit/revision_bjp_2026')
    ap.add_argument('--out', default='/mnt/data/provenance_audit/audit_output')
    a=ap.parse_args(); base=Path(a.bundle_root); rev=Path(a.revision_root); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    strict=pd.read_csv(base/'results/microswitches/publication_controls/strict_receptor_only_structures.csv')
    p1=pd.read_csv(rev/'A01_expanded_matched_set/P1_terminal_strict_free_inventory.tsv',sep='\t')
    lig=pd.read_csv(base/'data/processed/receptor_only_audit/receptor_only_ligand_permissive_structures.csv')
    p1_tlig=pd.read_csv(rev/'A01_expanded_matched_set/P1_terminal_ligand_permissive_free_inventory.tsv',sep='\t')
    audits={
      'strict_51':audit_frame(strict,base,'strict_51'),
      'terminal_P1_56':audit_frame(p1,base,'terminal_P1_56'),
      'terminal_ligand_permissive_P1_59':audit_frame(p1_tlig,base,'terminal_ligand_permissive_P1_59'),
      'ligand_permissive_58':audit_frame(lig,base,'ligand_permissive_58')}
    for k,v in audits.items(): v.to_csv(out/f'{k}_provenance_audit.tsv',sep='\t',index=False)

    # Shared-particle evidence directly from deposited mmCIF metadata.
    spr=[]
    for local,full,note in SHARED_PARTICLE_PAIRS:
        a1=cif_meta(base,local); a2=cif_meta(base,full)
        spr.append({'focused_or_local_pdb':local,'full_complex_pdb':full,
                    'focused_particles':a1['em_particles'],'full_complex_particles':a2['em_particles'],
                    'exact_particle_count_match':a1['em_particles']!='' and a1['em_particles']==a2['em_particles'],
                    'focused_title':a1['cif_title'],'full_title':a2['cif_title'],'note':note})
    pd.DataFrame(spr).to_csv(out/'shared_particle_evidence.tsv',sep='\t',index=False)

    reps=pd.read_csv(base/'results/microswitches/coupling_analysis_representatives.csv')
    reps['pdb_id']=reps.pdb_id.astype(str).str.upper()
    scenarios={
      'corrected_strict_P0':audits['strict_51'],
      'corrected_terminal_P1':audits['terminal_P1_56'],
      'corrected_terminal_ligand_permissive_P1':audits['terminal_ligand_permissive_P1_59'],
      'corrected_ligand_permissive':audits['ligand_permissive_58'],
    }
    summaries=[]
    for name,adf in scenarios.items():
        eligible=adf[adf.final_eligibility].copy()
        eligible.to_csv(out/f'{name}_eligible_structures.tsv',sep='\t',index=False)
        membership=(eligible.groupby('receptor_name',sort=True)
                    .agg(n_partner_free_structures=('pdb_id','size'), partner_free_pdbs=('pdb_id',lambda x:';'.join(x))))
        membership.to_csv(out/f'{name}_reference_membership.tsv',sep='\t')
        recs=sorted(eligible.receptor_name.unique())
        units=reps[reps.receptor_name.isin(recs) & reps.transducer_family.isin(['Gs','Gi/o','Gq/11'])].copy()
        units=units[['receptor_name','transducer_family','pdb_id','chain_id']].drop_duplicates()
        units.to_csv(out/f'{name}_matched_bound_units.tsv',sep='\t',index=False)
        class_rec={c:int(units.loc[units.transducer_family.eq(c),'receptor_name'].nunique()) for c in ['Gs','Gi/o','Gq/11']}
        three_classes=all(class_rec[c]>0 for c in class_rec)
        grouped_train_test_possible=all(class_rec[c]>=2 for c in class_rec) # minimum necessary condition
        summaries.append({
          'scenario':name,'eligible_partner_free_structures':len(eligible),'eligible_receptors':len(recs),
          'matched_receptor_class_units':len(units),'Gs_unique_receptors':class_rec['Gs'],'Gi_o_unique_receptors':class_rec['Gi/o'],
          'Gq_11_unique_receptors':class_rec['Gq/11'],'all_three_classes_present':three_classes,
          'minimum_grouped_train_test_class_coverage_possible':grouped_train_test_possible,
          'model_status':'NOT ESTIMABLE as original receptor-grouped 3-class analysis' if not grouped_train_test_possible else 'coverage prerequisite met; further fold audit required'})
    sdf=pd.DataFrame(summaries)
    sdf.to_csv(out/'feasibility_summary.tsv',sep='\t',index=False)
    with open(out/'feasibility_summary.json','w') as f: json.dump(summaries,f,indent=2)

    # Decision counts and receptor-level audit summary
    allsum=[]
    for k,v in audits.items():
        for status,n in v.provenance_status.value_counts().items():
            allsum.append({'population':k,'provenance_status':status,'n_structures':int(n),'n_receptors':int(v.loc[v.provenance_status.eq(status),'receptor_name'].nunique())})
    pd.DataFrame(allsum).to_csv(out/'audit_status_counts.tsv',sep='\t',index=False)

    inputs=[base/'results/microswitches/publication_controls/strict_receptor_only_structures.csv',
            rev/'A01_expanded_matched_set/P1_terminal_strict_free_inventory.tsv',
            rev/'A01_expanded_matched_set/P1_terminal_ligand_permissive_free_inventory.tsv',
            base/'data/processed/receptor_only_audit/receptor_only_ligand_permissive_structures.csv']
    with open(out/'INPUT_SHA256SUMS.txt','w') as f:
        for p in inputs: f.write(f'{sha256(p)}  {p}\n')

if __name__=='__main__': main()
