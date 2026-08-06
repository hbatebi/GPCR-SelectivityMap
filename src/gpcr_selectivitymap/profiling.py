from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any
import pandas as pd
from .features import StructureTask, extract_structure_features, feature_dictionary
from .io import read_table, write_json, write_table
from .reporting import html_report, profile_svg

def profile_structure(*, structure: str | Path, mapping: str | Path, pdb_id: str,
                      receptor_name: str, chain_id: str, output: str | Path,
                      dataset_role: str = "profile", transducer_family: str = "",
                      sasa_points: int = 80, contact_cutoff: float = 8.0,
                      anm_cutoff: float = 15.0, anm_modes: int = 12) -> dict[str, Any]:
    structure = Path(structure).resolve()
    if not structure.exists(): raise FileNotFoundError(structure)
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    structures = out / "_structure_input"; structures.mkdir(exist_ok=True)
    target = structures / f"{pdb_id.upper()}.cif"
    if target.resolve() != structure:
        target.write_bytes(structure.read_bytes())
    task = StructureTask(pdb_id.upper(), receptor_name, chain_id, dataset_role, transducer_family)
    row = extract_structure_features(task, mapping, structures, sasa_points, contact_cutoff, anm_cutoff, anm_modes)
    frame = pd.DataFrame([row])
    write_table(frame, out / "receptor_profile.tsv")
    write_json(out / "receptor_profile.json", row)
    write_table(feature_dictionary(frame), out / "feature_dictionary.tsv")
    profile_svg(row, out / "profile.svg")
    html_report(row, out / "report.html", "profile.svg")
    return row

def _batch_worker(payload: dict[str, Any]) -> dict[str, Any]:
    task = StructureTask(payload["pdb_id"], payload["receptor_name"], payload.get("chain_id", ""),
                         payload.get("dataset_role", "dataset"), payload.get("transducer_family", ""))
    try:
        row = extract_structure_features(task, payload["mapping"], payload["structures_dir"],
                                         payload["sasa_points"], payload["contact_cutoff"],
                                         payload["anm_cutoff"], payload["anm_modes"])
        row["extraction_status"] = "ok"; row["extraction_error"] = ""
        return row
    except Exception as error:
        return {"pdb_id": task.pdb_id, "receptor_name": task.receptor_name, "chain_id": task.chain_id,
                "dataset_role": task.dataset_role, "transducer_family": task.transducer_family,
                "extraction_status": "failed", "extraction_error": repr(error)}

def batch_profile(*, manifest: str | Path, mapping: str | Path, structures_dir: str | Path,
                  output: str | Path, workers: int = 1, sasa_points: int = 80,
                  contact_cutoff: float = 8.0, anm_cutoff: float = 15.0,
                  anm_modes: int = 12) -> pd.DataFrame:
    mf = read_table(manifest)
    required = {"pdb_id", "receptor_name"}
    missing = required.difference(mf.columns)
    if missing: raise ValueError(f"Manifest missing columns: {sorted(missing)}")
    payloads=[]
    for r in mf.fillna("").to_dict(orient="records"):
        payloads.append({**r, "pdb_id": str(r["pdb_id"]).upper(), "mapping": str(mapping),
                         "structures_dir": str(structures_dir), "sasa_points": sasa_points,
                         "contact_cutoff": contact_cutoff, "anm_cutoff": anm_cutoff, "anm_modes": anm_modes})
    rows=[]
    if workers <= 1:
        rows=[_batch_worker(p) for p in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(_batch_worker,p) for p in payloads]
            for f in as_completed(futures): rows.append(f.result())
    frame=pd.DataFrame(rows).sort_values(["receptor_name","pdb_id"]).reset_index(drop=True)
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    write_table(frame,out/"receptor_level_features.tsv")
    write_table(feature_dictionary(frame),out/"feature_dictionary.tsv")
    audit=frame[[c for c in ["pdb_id","receptor_name","chain_id","extraction_status","extraction_error","quality_contact","quality_surface","quality_electro","quality_mechanical"] if c in frame.columns]]
    write_table(audit,out/"extraction_audit.tsv")
    return frame
