from __future__ import annotations

"""GPCR–Gα interface chemistry and receptor-alone compatibility analysis.

This module implements two related but explicitly separated analyses:

1. Complex-side audit: extract physical receptor–Gα contacts from experimental
   complexes and compare the observed Gα chemistry with family-decoy chemistry.
2. Receptor-alone mode: derive Gα-family contact/chemistry templates strictly
   from training complexes and score a receptor sequence or receptor-only
   structure against those templates without using a test G protein.

Scores are compatibility descriptors, not calibrated functional coupling
probabilities.  Validation is receptor-aware and can additionally hold out
30%, 40%, and 50% receptor-sequence identity clusters.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence
import hashlib
import json
import math
import os
import re
import time
import warnings

import numpy as np
import pandas as pd
from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .io import read_table, write_json, write_table
from .sequence_audit import DEFAULT_INTERFACE_POSITIONS, pair_identity, sequence_columns

FAMILIES: tuple[str, ...] = ("Gs", "Gi/o", "Gq/11")

# Human GNAS short isoform.  The project already used this sequence to place
# heterogeneous Gα chains onto a common reference coordinate system.
GALPHA_S_REFERENCE = (
    "MGCLGNSKTEDQRNEEKAQREANKKIEKQLQKDKQVYRATHRLLLLGAGESGKSTIVKQMRILH"
    "VNGFNGEGGEEDPQAARSNSDGEKATKVQDIKNNLKEAIETIVAAMSNLVPPVELANPENQFRV"
    "DYILSVMNVPDFDFPPEFYEHAKALWEDEGVRACYERSNEYQLIDCAQYFLDKIDVIKQADYVP"
    "SDQDLLRCRVLTSGIFETKFQVDKVNFHMFDVGGQRDERRKWIQCFNDVTAIIFVVASSSYNMV"
    "IREDNQTNRLQEALNLFKSIWNNRWLRTISVILFLNKQDLLAEKVLAGKSKIEDYFPEFARYTT"
    "PEDATPEPGEDPRVTRAKYFIRDEFLRISTASGDGRHYCYPHFTCAVDTENIRRVFNDCRDIIQ"
    "RMHLRQYELL"
)

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {aa: 0.0 for aa in KD}
CHARGE.update({"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1})
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}
AROMATIC = {aa: float(aa in {"F", "W", "Y", "H"}) for aa in KD}
VOL_RANGE = max(VOLUME.values()) - min(VOLUME.values())

MODEL_NAMES: tuple[str, ...] = (
    "hydropathy",
    "charge",
    "hydropathy_charge",
    "volume",
    "combined_physics",
)


def _clean_aa(value: object) -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip().upper()
    if value in KD:
        return value
    if len(value) == 3:
        try:
            converted = seq1(value, custom_map={"MSE": "M"}).upper()
            return converted if converted in KD else ""
        except Exception:
            return ""
    return ""


def normalize_family(value: object) -> str:
    text = str(value).strip().lower().replace(" ", "")
    if text in {"gs", "g_s"}:
        return "Gs"
    if text in {"gi/o", "gio", "gi", "go", "gi1", "gi2", "gi3"}:
        return "Gi/o"
    if text in {"gq/11", "gq11", "gq", "g11"}:
        return "Gq/11"
    return str(value).strip()


def normalize_generic_number(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (int, float, np.integer, np.floating)):
        number = float(value)
        helix = int(number)
        position = int(round((number - helix) * 100))
        return f"{helix}x{position:02d}"
    text = str(value).strip().lower().replace(".", "x")
    match = re.fullmatch(r"(\d+)x(\d+)", text)
    if match:
        digits = match.group(2)
        position = int(digits.ljust(2, "0")[:2]) if len(digits) < 2 else int(digits)
        return f"{int(match.group(1))}x{position:02d}"
    try:
        number = float(value)
        helix = int(number)
        position = int(round((number - helix) * 100))
        return f"{helix}x{position:02d}"
    except Exception:
        return text


def galpha_element(reference_position: int) -> str:
    p = int(reference_position)
    if 8 <= p <= 45:
        return "alphaN"
    if 200 <= p <= 218:
        return "b2b3"
    if 236 <= p <= 255:
        return "a4b6"
    if 360 <= p <= 394:
        return "alpha5"
    return "other"


def aa_descriptor(aa: str) -> dict[str, float]:
    aa = _clean_aa(aa)
    if not aa:
        return {"kd": np.nan, "charge": np.nan, "volume": np.nan, "aromatic": np.nan}
    return {
        "kd": float(KD[aa]),
        "charge": float(CHARGE[aa]),
        "volume": float(VOLUME[aa]),
        "aromatic": float(AROMATIC[aa]),
    }


def pair_components(receptor_aa: str, galpha_aa: str) -> dict[str, float]:
    r = aa_descriptor(receptor_aa)
    g = aa_descriptor(galpha_aa)
    if not np.isfinite(r["kd"]) or not np.isfinite(g["kd"]):
        return {name: np.nan for name in (
            "hydropathy_similarity", "hydrophobic_packing", "charge_complementarity",
            "volume_similarity", "aromatic_pair", "hydropathy", "charge",
            "hydropathy_charge", "volume", "combined_physics",
        )}
    hydro_similarity = 1.0 - abs(r["kd"] - g["kd"]) / 9.0
    hydro_packing = max(r["kd"], 0.0) * max(g["kd"], 0.0) / (4.5 * 4.5)
    charge_comp = -r["charge"] * g["charge"]
    volume_similarity = 1.0 - abs(r["volume"] - g["volume"]) / VOL_RANGE
    aromatic_pair = r["aromatic"] * g["aromatic"]
    hydro = 0.60 * hydro_similarity + 0.40 * hydro_packing
    charge_scaled = (charge_comp + 1.0) / 2.0
    hydro_charge = 0.65 * hydro + 0.35 * charge_scaled
    combined = (
        0.30 * hydro_similarity
        + 0.20 * hydro_packing
        + 0.25 * charge_scaled
        + 0.20 * volume_similarity
        + 0.05 * aromatic_pair
    )
    return {
        "hydropathy_similarity": float(hydro_similarity),
        "hydrophobic_packing": float(hydro_packing),
        "charge_complementarity": float(charge_comp),
        "volume_similarity": float(volume_similarity),
        "aromatic_pair": float(aromatic_pair),
        "hydropathy": float(hydro),
        "charge": float(charge_scaled),
        "hydropathy_charge": float(hydro_charge),
        "volume": float(volume_similarity),
        "combined_physics": float(combined),
    }


def _chain_sequence(chain) -> tuple[str, list]:
    residues, letters = [], []
    for residue in chain:
        if residue.id[0] != " " or "CA" not in residue:
            continue
        aa = _clean_aa(residue.resname)
        letters.append(aa or "X")
        residues.append(residue)
    return "".join(letters), residues


def _align_chain_to_galpha(chain) -> tuple[dict[tuple, int], float, float, float]:
    sequence, residues = _chain_sequence(chain)
    if not sequence:
        return {}, 0.0, 0.0, -np.inf
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    try:
        alignment = aligner.align(sequence, GALPHA_S_REFERENCE)[0]
    except Exception:
        return {}, 0.0, 0.0, -np.inf
    mapping: dict[tuple, int] = {}
    matches = 0
    aligned = 0
    query_blocks, target_blocks = alignment.aligned
    for (qs, qe), (ts, te) in zip(query_blocks, target_blocks):
        for offset in range(int(qe - qs)):
            q = int(qs + offset)
            t = int(ts + offset)
            mapping[residues[q].id] = t + 1
            aligned += 1
            matches += int(sequence[q] == GALPHA_S_REFERENCE[t])
    identity = matches / aligned if aligned else 0.0
    coverage = aligned / len(sequence) if sequence else 0.0
    return mapping, float(identity), float(coverage), float(alignment.score)


def _heavy_coordinates(chain) -> np.ndarray:
    coords = []
    for residue in chain:
        for atom in residue:
            element = str(getattr(atom, "element", "")).upper()
            if element == "H" or atom.name.upper().startswith("H"):
                continue
            coords.append(np.asarray(atom.coord, dtype=float))
    return np.asarray(coords, dtype=float) if coords else np.empty((0, 3), dtype=float)


def _minimum_chain_distance(left, right) -> float:
    a = _heavy_coordinates(left)
    b = _heavy_coordinates(right)
    if len(a) == 0 or len(b) == 0:
        return np.inf
    tree = cKDTree(b)
    distance, _ = tree.query(a, k=1)
    return float(np.min(distance))


def _mapping_for_structure(mapping: pd.DataFrame, pdb_id: str, chain_id: str) -> dict[tuple[int, str], str]:
    part = mapping[
        mapping["pdb_id"].astype(str).str.upper().eq(str(pdb_id).upper())
        & mapping["chain_id"].astype(str).eq(str(chain_id))
    ]
    out: dict[tuple[int, str], str] = {}
    for row in part.itertuples(index=False):
        author = getattr(row, "author_residue_number", np.nan)
        if pd.isna(author):
            continue
        insertion = getattr(row, "author_insertion_code", " ")
        insertion = " " if pd.isna(insertion) else str(insertion).strip() or " "
        generic = normalize_generic_number(getattr(row, "generic_number", ""))
        if generic:
            out[(int(float(author)), insertion)] = generic
    return out


def _select_receptor_chain(model, row: Mapping[str, object], mapping: pd.DataFrame) -> str:
    available = {chain.id for chain in model}
    for key in ("chain_id", "preferred_chain"):
        candidate = str(row.get(key, "")).strip()
        if candidate and candidate in available:
            return candidate
    pdb = str(row.get("pdb_id", "")).upper()
    part = mapping[mapping["pdb_id"].astype(str).str.upper().eq(pdb)]
    if not part.empty:
        counts = part.groupby("chain_id").size().sort_values(ascending=False)
        for candidate in counts.index.astype(str):
            if candidate in available:
                return candidate
    raise ValueError(f"Could not resolve receptor chain for {pdb}")


def _select_galpha_chain(model, receptor_chain_id: str, preferred: str = "") -> tuple[str, dict, dict]:
    candidates = []
    for chain in model:
        if chain.id == receptor_chain_id:
            continue
        sequence, _ = _chain_sequence(chain)
        if len(sequence) < 120:
            continue
        mapping, identity, coverage, score = _align_chain_to_galpha(chain)
        # Engineered mini-G and chimeric G proteins can be only ~40–50% identical
        # to GNAS while still aligning strongly over the full Ras-like domain.
        if coverage < 0.45 or identity < 0.38 or score < 250.0:
            continue
        distance = _minimum_chain_distance(model[receptor_chain_id], chain)
        candidates.append({
            "chain": chain.id,
            "mapping": mapping,
            "identity": identity,
            "coverage": coverage,
            "alignment_score": score,
            "receptor_distance": distance,
            "preferred": int(str(chain.id) == str(preferred)),
        })
    if not candidates:
        raise ValueError("No G alpha-like chain found")
    candidates.sort(key=lambda x: (x["receptor_distance"], -x["preferred"], -x["identity"]))
    chosen = candidates[0]
    audit = {k: v for k, v in chosen.items() if k not in {"mapping"}}
    audit["n_galpha_candidates"] = len(candidates)
    audit["candidate_chains"] = ";".join(str(c["chain"]) for c in candidates)
    return str(chosen["chain"]), chosen["mapping"], audit


def _extract_one_structure(task: dict) -> dict:
    pdb_id = str(task["pdb_id"]).upper()
    structure_path = Path(task["structure_path"])
    mapping_path = Path(task["mapping_path"])
    row = task["row"]
    cutoff = float(task.get("cutoff", 4.5))
    preferred_galpha = str(task.get("preferred_galpha", ""))
    mapping = pd.read_csv(mapping_path)
    started = time.time()
    result = {"pdb_id": pdb_id, "status": "failed", "contacts": [], "galpha_profile": []}
    try:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(structure_path))
        model = next(structure.get_models())
        receptor_chain = _select_receptor_chain(model, row, mapping)
        receptor_map = _mapping_for_structure(mapping, pdb_id, receptor_chain)
        if not receptor_map:
            raise ValueError(f"No generic-number mapping for receptor chain {receptor_chain}")
        galpha_chain, galpha_map, chain_audit = _select_galpha_chain(
            model, receptor_chain, preferred=preferred_galpha
        )
        receptor_atoms = []
        for residue in model[receptor_chain]:
            insertion = str(residue.id[2]).strip() or " "
            key = (int(residue.id[1]), insertion)
            generic = receptor_map.get(key)
            if not generic:
                continue
            aa = _clean_aa(residue.resname)
            if not aa:
                continue
            for atom in residue:
                element = str(getattr(atom, "element", "")).upper()
                if element == "H" or atom.name.upper().startswith("H"):
                    continue
                receptor_atoms.append((np.asarray(atom.coord, dtype=float), residue, generic, aa))
        galpha_atoms = []
        profile_rows = []
        for residue in model[galpha_chain]:
            ref_position = galpha_map.get(residue.id)
            if ref_position is None:
                continue
            aa = _clean_aa(residue.resname)
            if not aa:
                continue
            descriptor = aa_descriptor(aa)
            profile_rows.append({
                "pdb_id": pdb_id,
                "transducer_family": normalize_family(row.get("transducer_family", row.get("transducer_type", ""))),
                "galpha_chain": galpha_chain,
                "galpha_author_residue": int(residue.id[1]),
                "galpha_reference_position": int(ref_position),
                "galpha_element": galpha_element(int(ref_position)),
                "galpha_aa": aa,
                "galpha_kd": descriptor["kd"],
                "galpha_charge": descriptor["charge"],
                "galpha_volume": descriptor["volume"],
                "galpha_aromatic": descriptor["aromatic"],
            })
            for atom in residue:
                element = str(getattr(atom, "element", "")).upper()
                if element == "H" or atom.name.upper().startswith("H"):
                    continue
                galpha_atoms.append((np.asarray(atom.coord, dtype=float), residue, int(ref_position), aa))
        if not receptor_atoms or not galpha_atoms:
            raise ValueError("Empty receptor or G alpha atom set")
        tree = cKDTree(np.asarray([entry[0] for entry in galpha_atoms], dtype=float))
        pairs: dict[tuple, dict] = {}
        for coordinate, receptor_residue, receptor_position, receptor_aa in receptor_atoms:
            for neighbour in tree.query_ball_point(coordinate, cutoff):
                g_coordinate, g_residue, g_position, g_aa = galpha_atoms[int(neighbour)]
                distance = float(np.linalg.norm(coordinate - g_coordinate))
                key = (
                    receptor_position,
                    int(receptor_residue.id[1]),
                    int(g_position),
                    int(g_residue.id[1]),
                )
                if key not in pairs:
                    pairs[key] = {
                        "min_distance": distance,
                        "atom_contact_count": 1,
                        "receptor_aa": receptor_aa,
                        "galpha_aa": g_aa,
                    }
                else:
                    pairs[key]["min_distance"] = min(pairs[key]["min_distance"], distance)
                    pairs[key]["atom_contact_count"] += 1
        family = normalize_family(row.get("transducer_family", row.get("transducer_type", "")))
        contacts = []
        for key, values in pairs.items():
            rpos, r_author, gpos, g_author = key
            components = pair_components(values["receptor_aa"], values["galpha_aa"])
            weight = (1.0 + math.log1p(values["atom_contact_count"])) * max(
                0.05, 1.0 - max(values["min_distance"] - 2.5, 0.0) / max(cutoff - 2.5, 0.1)
            )
            contacts.append({
                "pdb_id": pdb_id,
                "receptor_name": str(row.get("receptor_name", "")),
                "transducer_family": family,
                "receptor_chain": receptor_chain,
                "galpha_chain": galpha_chain,
                "receptor_generic_position": rpos,
                "receptor_author_residue": r_author,
                "receptor_aa": values["receptor_aa"],
                "galpha_reference_position": gpos,
                "galpha_author_residue": g_author,
                "galpha_element": galpha_element(gpos),
                "galpha_aa": values["galpha_aa"],
                "min_distance": values["min_distance"],
                "atom_contact_count": values["atom_contact_count"],
                "contact_weight": weight,
                **components,
            })
        if not contacts:
            raise ValueError(f"No receptor–G alpha contacts within {cutoff:.2f} Å")
        result.update({
            "status": "success",
            "contacts": contacts,
            "galpha_profile": profile_rows,
            "receptor_chain": receptor_chain,
            "galpha_chain": galpha_chain,
            "n_contacts": len(contacts),
            "n_receptor_positions": len({r["receptor_generic_position"] for r in contacts}),
            "n_galpha_positions": len({r["galpha_reference_position"] for r in contacts}),
            **chain_audit,
        })
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    result["elapsed_seconds"] = time.time() - started
    return result


def extract_contact_dataset(
    *,
    frame: pd.DataFrame,
    mapping_path: str | Path,
    structures_dir: str | Path,
    gprotein_contacts: pd.DataFrame | None,
    output: str | Path,
    workers: int = 1,
    cutoff: float = 4.5,
    quick_limit: int | None = None,
    heartbeat_seconds: int = 60,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = Path(output)
    checkpoint_dir = out / "checkpoints" / "contacts"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    representatives = frame.copy()
    if "representative_for_receptor_transducer" in representatives.columns:
        mask = representatives["representative_for_receptor_transducer"].fillna(False).astype(bool)
        representatives = representatives[mask]
    family_col = "transducer_family" if "transducer_family" in representatives.columns else "transducer_type"
    representatives[family_col] = representatives[family_col].map(normalize_family)
    representatives = representatives[representatives[family_col].isin(FAMILIES)].copy()
    representatives = representatives.drop_duplicates("pdb_id").sort_values("pdb_id")
    if quick_limit:
        sampled = []
        per_family = max(2, int(math.ceil(quick_limit / len(FAMILIES))))
        for family in FAMILIES:
            sampled.append(representatives[representatives[family_col].eq(family)].head(per_family))
        representatives = pd.concat(sampled, ignore_index=True).head(int(quick_limit))
    preferred_map = {}
    if gprotein_contacts is not None and not gprotein_contacts.empty:
        id_col = "id" if "id" in gprotein_contacts.columns else "pdb_id"
        if "transducer_chain" in gprotein_contacts.columns:
            preferred_map = dict(zip(
                gprotein_contacts[id_col].astype(str).str.upper(),
                gprotein_contacts["transducer_chain"].fillna("").astype(str),
            ))
    tasks, cached = [], []
    structures_dir = Path(structures_dir)
    for row in representatives.to_dict(orient="records"):
        pdb = str(row["pdb_id"]).upper()
        checkpoint = checkpoint_dir / f"{pdb}.json"
        if checkpoint.exists():
            try:
                cached.append(json.loads(checkpoint.read_text()))
                continue
            except Exception:
                checkpoint.unlink(missing_ok=True)
        structure_path = structures_dir / f"{pdb}.cif"
        if not structure_path.exists():
            structure_path = structures_dir / f"{pdb.lower()}.cif"
        tasks.append({
            "pdb_id": pdb,
            "structure_path": str(structure_path),
            "mapping_path": str(mapping_path),
            "row": row,
            "preferred_galpha": preferred_map.get(pdb, ""),
            "cutoff": cutoff,
        })
    results = list(cached)
    total = len(tasks)
    started = time.time()
    last_heartbeat = started
    if total:
        with ProcessPoolExecutor(max_workers=max(1, int(workers))) as executor:
            futures = {executor.submit(_extract_one_structure, task): task for task in tasks}
            done = 0
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                done += 1
                checkpoint = checkpoint_dir / f"{result['pdb_id']}.json"
                checkpoint.write_text(json.dumps(result, indent=2, allow_nan=True))
                now = time.time()
                if now - last_heartbeat >= heartbeat_seconds or done == total:
                    rate = done / max(now - started, 1e-9)
                    eta = (total - done) / rate if rate > 0 else np.nan
                    print(
                        f"[contacts] {done}/{total} new tasks; "
                        f"success={sum(r.get('status') == 'success' for r in results)}; "
                        f"elapsed={(now-started)/60:.1f} min; ETA={eta/60:.1f} min",
                        flush=True,
                    )
                    last_heartbeat = now
    contacts = pd.DataFrame([
        contact for result in results if result.get("status") == "success"
        for contact in result.get("contacts", [])
    ])
    profiles = pd.DataFrame([
        residue for result in results if result.get("status") == "success"
        for residue in result.get("galpha_profile", [])
    ])
    audit = pd.DataFrame([
        {k: v for k, v in result.items() if k not in {"contacts", "galpha_profile"}}
        for result in results
    ])
    write_table(contacts, out / "extraction" / "contact_pairs.tsv")
    write_table(profiles, out / "extraction" / "galpha_profiles.tsv")
    write_table(audit, out / "extraction" / "contact_extraction_audit.tsv")
    return contacts, profiles, audit


def receptor_profile_from_row(
    row: Mapping[str, object],
    positions: Iterable[str] | None = None,
    *,
    profile_mode: str = "sequence",
) -> dict[str, str]:
    """Return receptor interface residues.

    profile_mode="sequence" uses expected receptor sequence and therefore works
    without a structure.  "observed" uses only residues modelled in the supplied
    structure.  "observed_then_sequence" prefers modelled identities but falls
    back to expected sequence.  The observed mode is a structure-coverage mode;
    it does not claim to encode side-chain exposure or dynamics.
    """
    positions = tuple(positions or DEFAULT_INTERFACE_POSITIONS)
    profile = {}
    for position in positions:
        expected = (
            f"gpcrdb_{position}_expected_aa",
            f"gpcrdb_{position.replace('x', '.')}_expected_aa",
        )
        observed = (f"residue_{position}_identity",)
        if profile_mode == "observed":
            candidates = observed
        elif profile_mode == "observed_then_sequence":
            candidates = observed + expected
        else:
            candidates = expected + observed
        aa = ""
        for column in candidates:
            if column in row:
                aa = _clean_aa(row[column])
                if aa:
                    break
        if aa:
            profile[position] = aa
    return profile


def _profile_table(frame: pd.DataFrame, positions: Iterable[str], profile_mode: str = "sequence") -> dict[str, dict[str, str]]:
    return {
        str(row["pdb_id"]).upper(): receptor_profile_from_row(row, positions, profile_mode=profile_mode)
        for row in frame.to_dict(orient="records")
    }


def build_family_templates(
    contacts: pd.DataFrame,
    galpha_profiles: pd.DataFrame,
    receptor_profiles: Mapping[str, Mapping[str, str]],
    train_pdbs: Iterable[str],
    *,
    minimum_edge_count: int = 2,
    minimum_edge_frequency: float = 0.03,
) -> dict[str, pd.DataFrame]:
    train_set = {str(p).upper() for p in train_pdbs}
    c = contacts[contacts["pdb_id"].astype(str).str.upper().isin(train_set)].copy()
    g = galpha_profiles[galpha_profiles["pdb_id"].astype(str).str.upper().isin(train_set)].copy()
    templates: dict[str, pd.DataFrame] = {}
    for family in FAMILIES:
        cf = c[c["transducer_family"].map(normalize_family).eq(family)].copy()
        gf = g[g["transducer_family"].map(normalize_family).eq(family)].copy()
        family_pdbs = sorted(cf["pdb_id"].astype(str).str.upper().unique())
        n_family = len(family_pdbs)
        if n_family == 0:
            templates[family] = pd.DataFrame()
            continue
        edges = (
            cf.groupby(["receptor_generic_position", "galpha_reference_position"], as_index=False)
            .agg(
                n_structures=("pdb_id", "nunique"),
                mean_distance=("min_distance", "mean"),
                mean_atom_contacts=("atom_contact_count", "mean"),
                mean_contact_weight=("contact_weight", "mean"),
            )
        )
        edges["edge_frequency"] = edges["n_structures"] / n_family
        edges = edges[
            (edges["n_structures"] >= int(minimum_edge_count))
            & (edges["edge_frequency"] >= float(minimum_edge_frequency))
        ].copy()
        if edges.empty:
            templates[family] = edges
            continue
        gstats = (
            gf.groupby("galpha_reference_position", as_index=False)
            .agg(
                galpha_kd=("galpha_kd", "mean"),
                galpha_charge=("galpha_charge", "mean"),
                galpha_volume=("galpha_volume", "mean"),
                galpha_aromatic=("galpha_aromatic", "mean"),
                galpha_position_coverage=("pdb_id", "nunique"),
            )
        )
        edges = edges.merge(gstats, on="galpha_reference_position", how="left")
        # Receptor residue frequency is retained only as a diagnostic upper comparator.
        receptor_rows = []
        for pdb in family_pdbs:
            for position, aa in receptor_profiles.get(pdb, {}).items():
                receptor_rows.append((position, aa))
        if receptor_rows:
            rr = pd.DataFrame(receptor_rows, columns=["receptor_generic_position", "receptor_aa"])
            modal = rr.groupby("receptor_generic_position")["receptor_aa"].agg(
                lambda values: values.mode().iloc[0] if len(values.mode()) else ""
            ).rename("family_modal_receptor_aa").reset_index()
            edges = edges.merge(modal, on="receptor_generic_position", how="left")
        edges["template_weight"] = (
            edges["edge_frequency"]
            * (1.0 + np.log1p(edges["mean_atom_contacts"].clip(lower=0)))
            * np.exp(-np.maximum(edges["mean_distance"] - 2.5, 0.0) / 2.0)
        )
        edges["transducer_family"] = family
        edges["n_training_complexes"] = n_family
        templates[family] = edges.sort_values(
            ["receptor_generic_position", "galpha_reference_position"]
        ).reset_index(drop=True)
    return templates


def _components_against_mean(receptor_aa: str, row: Mapping[str, object]) -> dict[str, float]:
    r = aa_descriptor(receptor_aa)
    if not np.isfinite(r["kd"]):
        return {name: np.nan for name in MODEL_NAMES}
    g_kd = float(row.get("galpha_kd", np.nan))
    g_charge = float(row.get("galpha_charge", np.nan))
    g_volume = float(row.get("galpha_volume", np.nan))
    g_aromatic = float(row.get("galpha_aromatic", np.nan))
    if not np.isfinite(g_kd):
        return {name: np.nan for name in MODEL_NAMES}
    hydro_similarity = 1.0 - abs(r["kd"] - g_kd) / 9.0
    hydro_packing = max(r["kd"], 0.0) * max(g_kd, 0.0) / (4.5 * 4.5)
    charge_comp = -r["charge"] * g_charge
    charge_scaled = (charge_comp + 1.0) / 2.0
    volume_similarity = 1.0 - abs(r["volume"] - g_volume) / VOL_RANGE
    aromatic_pair = r["aromatic"] * g_aromatic
    hydro = 0.60 * hydro_similarity + 0.40 * hydro_packing
    return {
        "hydropathy": float(hydro),
        "charge": float(charge_scaled),
        "hydropathy_charge": float(0.65 * hydro + 0.35 * charge_scaled),
        "volume": float(volume_similarity),
        "combined_physics": float(
            0.30 * hydro_similarity
            + 0.20 * hydro_packing
            + 0.25 * charge_scaled
            + 0.20 * volume_similarity
            + 0.05 * aromatic_pair
        ),
    }


def score_receptor_against_template(
    receptor_profile: Mapping[str, str], template: pd.DataFrame
) -> dict[str, float]:
    if template is None or template.empty:
        return {**{name: np.nan for name in MODEL_NAMES}, "coverage": 0.0, "n_edges": 0}
    totals = {name: 0.0 for name in MODEL_NAMES}
    weights = {name: 0.0 for name in MODEL_NAMES}
    available_edges = 0
    all_weight = float(template["template_weight"].sum())
    for row in template.to_dict(orient="records"):
        position = str(row["receptor_generic_position"])
        aa = _clean_aa(receptor_profile.get(position, ""))
        if not aa:
            continue
        components = _components_against_mean(aa, row)
        weight = float(row.get("template_weight", 1.0))
        available_edges += 1
        for name, value in components.items():
            if np.isfinite(value):
                totals[name] += weight * value
                weights[name] += weight
    scores = {
        name: totals[name] / weights[name] if weights[name] > 0 else np.nan
        for name in MODEL_NAMES
    }
    used_weight = max(weights.values()) if weights else 0.0
    scores["coverage"] = used_weight / all_weight if all_weight > 0 else 0.0
    scores["n_edges"] = available_edges
    return scores


def _softmax(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).any():
        return np.full(len(values), 1.0 / len(values))
    fill = np.nanmin(values[np.isfinite(values)]) - 1.0
    values = np.where(np.isfinite(values), values, fill)
    shifted = values - np.max(values)
    exp = np.exp(shifted / 0.08)
    return exp / exp.sum()


def score_receptor_rows(
    frame: pd.DataFrame,
    templates: Mapping[str, pd.DataFrame],
    *,
    positions: Iterable[str],
    true_family_column: str | None = "transducer_family",
    mode: str = "receptor_alone",
    profile_mode: str = "sequence",
) -> pd.DataFrame:
    rows = []
    for source_row in frame.to_dict(orient="records"):
        profile = receptor_profile_from_row(source_row, positions, profile_mode=profile_mode)
        true_family = normalize_family(source_row.get(true_family_column, "")) if true_family_column else ""
        identity = {
            "pdb_id": str(source_row.get("pdb_id", "")).upper(),
            "receptor_name": str(source_row.get("receptor_name", "")),
            "true_family": true_family,
            "mode": mode,
        }
        family_scores: dict[str, dict[str, float]] = {
            family: score_receptor_against_template(profile, templates.get(family, pd.DataFrame()))
            for family in FAMILIES
        }
        for model in MODEL_NAMES:
            values = np.array([family_scores[family][model] for family in FAMILIES], dtype=float)
            probabilities = _softmax(values)
            order = np.argsort(-np.nan_to_num(values, nan=-np.inf))
            predicted = FAMILIES[int(order[0])] if len(order) else ""
            true_index = FAMILIES.index(true_family) if true_family in FAMILIES else None
            rank = int(np.where(order == true_index)[0][0] + 1) if true_index is not None else np.nan
            true_score = values[true_index] if true_index is not None else np.nan
            decoys = np.delete(values, true_index) if true_index is not None else np.array([])
            margin = true_score - np.nanmax(decoys) if len(decoys) and np.isfinite(decoys).any() else np.nan
            for family, score, probability in zip(FAMILIES, values, probabilities):
                rows.append({
                    **identity,
                    "model": model,
                    "candidate_family": family,
                    "compatibility_score": score,
                    "probability_like_score": probability,
                    "predicted_family": predicted,
                    "true_family_rank": rank,
                    "true_minus_best_decoy": margin,
                    "profile_positions": len(profile),
                    "template_coverage": family_scores[family]["coverage"],
                    "template_edges_used": family_scores[family]["n_edges"],
                })
    return pd.DataFrame(rows)


def score_complex_contact_map_decoys(
    test_pdbs: Iterable[str],
    contacts: pd.DataFrame,
    galpha_profiles: pd.DataFrame,
    train_pdbs: Iterable[str],
) -> pd.DataFrame:
    """Score actual chemistry and family-decoy Gα chemistry on each test contact map.

    The test complex contact graph is held fixed.  The observed Gα residues are
    compared with training-family mean residues at the same common Gα positions.
    This is a complex-side mechanistic audit and is not receptor-alone prediction.
    """
    train_set = {str(p).upper() for p in train_pdbs}
    test_set = {str(p).upper() for p in test_pdbs}
    train_g = galpha_profiles[galpha_profiles["pdb_id"].astype(str).str.upper().isin(train_set)].copy()
    means = (
        train_g.groupby(["transducer_family", "galpha_reference_position"], as_index=False)
        .agg(
            galpha_kd=("galpha_kd", "mean"),
            galpha_charge=("galpha_charge", "mean"),
            galpha_volume=("galpha_volume", "mean"),
            galpha_aromatic=("galpha_aromatic", "mean"),
        )
    )
    mean_lookup = {
        (normalize_family(row.transducer_family), int(row.galpha_reference_position)): row._asdict()
        for row in means.itertuples(index=False)
    }
    rows = []
    for pdb, part in contacts[contacts["pdb_id"].astype(str).str.upper().isin(test_set)].groupby("pdb_id"):
        true_family = normalize_family(part["transducer_family"].iloc[0])
        receptor_name = str(part["receptor_name"].iloc[0])
        observed_by_model = {}
        for model in MODEL_NAMES:
            observed_by_model[model] = np.average(
                part[model].astype(float), weights=part["contact_weight"].astype(float)
            )
        candidate_scores = {family: {model: [0.0, 0.0] for model in MODEL_NAMES} for family in FAMILIES}
        for contact in part.to_dict(orient="records"):
            raa = contact["receptor_aa"]
            gpos = int(contact["galpha_reference_position"])
            weight = float(contact["contact_weight"])
            for family in FAMILIES:
                mean = mean_lookup.get((family, gpos))
                if mean is None:
                    continue
                components = _components_against_mean(raa, mean)
                for model, value in components.items():
                    if np.isfinite(value):
                        candidate_scores[family][model][0] += weight * value
                        candidate_scores[family][model][1] += weight
        for model in MODEL_NAMES:
            values = np.array([
                candidate_scores[family][model][0] / candidate_scores[family][model][1]
                if candidate_scores[family][model][1] > 0 else np.nan
                for family in FAMILIES
            ])
            probabilities = _softmax(values)
            order = np.argsort(-np.nan_to_num(values, nan=-np.inf))
            true_index = FAMILIES.index(true_family)
            rank = int(np.where(order == true_index)[0][0] + 1)
            margin = values[true_index] - np.nanmax(np.delete(values, true_index))
            for family, score, probability in zip(FAMILIES, values, probabilities):
                rows.append({
                    "pdb_id": str(pdb).upper(),
                    "receptor_name": receptor_name,
                    "true_family": true_family,
                    "mode": "complex_contact_map_decoy",
                    "model": model,
                    "candidate_family": family,
                    "compatibility_score": score,
                    "probability_like_score": probability,
                    "predicted_family": FAMILIES[int(order[0])],
                    "true_family_rank": rank,
                    "true_minus_best_decoy": margin,
                    "observed_actual_pair_score": observed_by_model[model],
                })
    return pd.DataFrame(rows)


def _macro_metrics(predictions: pd.DataFrame) -> dict[str, float]:
    if predictions.empty:
        return {"macro_auc": np.nan, "log_loss": np.nan, "balanced_accuracy": np.nan, "top1_accuracy": np.nan}
    wide = predictions.pivot_table(
        index=["pdb_id", "receptor_name", "true_family"],
        columns="candidate_family",
        values="probability_like_score",
        aggfunc="first",
    ).reset_index()
    if not set(FAMILIES).issubset(wide.columns):
        return {"macro_auc": np.nan, "log_loss": np.nan, "balanced_accuracy": np.nan, "top1_accuracy": np.nan}
    y = wide["true_family"].astype(str).to_numpy()
    metric_classes = tuple(sorted(FAMILIES))
    prob = wide[list(metric_classes)].astype(float).to_numpy()
    prob = prob / np.maximum(prob.sum(axis=1, keepdims=True), 1e-12)
    predicted = np.asarray(metric_classes)[np.argmax(prob, axis=1)]
    if set(metric_classes).issubset(set(y)):
        try:
            auc = roc_auc_score(y, prob, labels=list(metric_classes), multi_class="ovr", average="macro")
        except Exception:
            auc = np.nan
    else:
        auc = np.nan
    try:
        ll = log_loss(y, prob, labels=list(metric_classes))
    except Exception:
        ll = np.nan
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bal = balanced_accuracy_score(y, predicted)
    except Exception:
        bal = np.nan
    return {
        "macro_auc": float(auc),
        "log_loss": float(ll),
        "balanced_accuracy": float(bal),
        "top1_accuracy": float(np.mean(predicted == y)),
        "mean_true_rank": float(
            predictions.drop_duplicates(["pdb_id", "model"])["true_family_rank"].mean()
        ),
        "mean_true_minus_decoy": float(
            predictions.drop_duplicates(["pdb_id", "model"])["true_minus_best_decoy"].mean()
        ),
        "n_units": int(len(wide)),
        "n_receptors": int(wide["receptor_name"].nunique()),
    }


def _bootstrap_metric(
    predictions: pd.DataFrame,
    *,
    bootstraps: int,
    seed: int,
) -> tuple[float, float]:
    if predictions.empty or bootstraps <= 0:
        return np.nan, np.nan
    metric_classes = tuple(sorted(FAMILIES))
    wide = predictions.pivot_table(
        index=["pdb_id", "receptor_name", "true_family"],
        columns="candidate_family",
        values="probability_like_score",
        aggfunc="first",
    ).reset_index()
    if not set(metric_classes).issubset(wide.columns):
        return np.nan, np.nan
    receptors = np.array(sorted(wide["receptor_name"].astype(str).unique()))
    if len(receptors) < 3:
        return np.nan, np.nan
    receptor_indices = {
        receptor: np.flatnonzero(wide["receptor_name"].astype(str).to_numpy() == receptor)
        for receptor in receptors
    }
    y_all = wide["true_family"].astype(str).to_numpy()
    prob_all = wide[list(metric_classes)].astype(float).to_numpy()
    prob_all = prob_all / np.maximum(prob_all.sum(axis=1, keepdims=True), 1e-12)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(int(bootstraps)):
        sampled = rng.choice(receptors, size=len(receptors), replace=True)
        indices = np.concatenate([receptor_indices[str(receptor)] for receptor in sampled])
        y = y_all[indices]
        if not set(metric_classes).issubset(set(y)):
            continue
        try:
            value = roc_auc_score(
                y, prob_all[indices], labels=list(metric_classes),
                multi_class="ovr", average="macro",
            )
        except Exception:
            continue
        if np.isfinite(value):
            values.append(float(value))
    if not values:
        return np.nan, np.nan
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def _identity_knn_predictions(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    sequence_cols: Sequence[str],
    k: int = 5,
) -> pd.DataFrame:
    train_receptors = {}
    for receptor, part in train.groupby("receptor_name"):
        signature = []
        for column in sequence_cols:
            values = part[column].map(_clean_aa)
            values = values[values.ne("")]
            signature.append(values.mode().iloc[0] if len(values) else "")
        labels = part["transducer_family"].map(normalize_family).value_counts()
        train_receptors[str(receptor)] = (tuple(signature), labels)
    rows = []
    for source in test.to_dict(orient="records"):
        signature = tuple(_clean_aa(source.get(column, "")) for column in sequence_cols)
        neighbours = []
        for receptor, (candidate, labels) in train_receptors.items():
            neighbours.append((pair_identity(signature, candidate), receptor, labels))
        neighbours.sort(key=lambda value: (-value[0], value[1]))
        selected = neighbours[: max(1, int(k))]
        votes = {family: 1e-6 for family in FAMILIES}
        for identity, _, labels in selected:
            weight = max(float(identity), 1e-3)
            for family in FAMILIES:
                votes[family] += weight * float(labels.get(family, 0.0))
        total = sum(votes.values())
        probability = {family: votes[family] / total for family in FAMILIES}
        predicted = max(FAMILIES, key=lambda family: probability[family])
        true = normalize_family(source["transducer_family"])
        order = sorted(FAMILIES, key=lambda family: probability[family], reverse=True)
        true_rank = order.index(true) + 1
        margin = probability[true] - max(probability[f] for f in FAMILIES if f != true)
        for family in FAMILIES:
            rows.append({
                "pdb_id": str(source["pdb_id"]).upper(),
                "receptor_name": str(source["receptor_name"]),
                "true_family": true,
                "mode": "baseline",
                "model": "global_identity_knn",
                "candidate_family": family,
                "compatibility_score": probability[family],
                "probability_like_score": probability[family],
                "predicted_family": predicted,
                "true_family_rank": true_rank,
                "true_minus_best_decoy": margin,
            })
    return pd.DataFrame(rows)


def _interface_logistic_predictions(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    positions: Sequence[str],
) -> pd.DataFrame:
    columns = [f"gpcrdb_{position}_expected_aa" for position in positions]
    columns = [column for column in columns if column in train.columns and column in test.columns]
    if not columns:
        return pd.DataFrame()
    transformer = ColumnTransformer([
        ("aa", Pipeline([
            ("impute", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore")),
        ]), columns)
    ])
    pipeline = Pipeline([
        ("preprocess", transformer),
        ("model", LogisticRegression(max_iter=30000, solver="lbfgs", C=1.0)),
    ])
    try:
        pipeline.fit(train[columns], train["transducer_family"].map(normalize_family))
        raw = pipeline.predict_proba(test[columns])
    except Exception:
        return pd.DataFrame()
    model_classes = list(pipeline.named_steps["model"].classes_)
    aligned = np.zeros((len(test), len(FAMILIES)), dtype=float)
    for i, family in enumerate(FAMILIES):
        if family in model_classes:
            aligned[:, i] = raw[:, model_classes.index(family)]
    aligned = aligned / np.maximum(aligned.sum(axis=1, keepdims=True), 1e-12)
    rows = []
    for source, probabilities in zip(test.to_dict(orient="records"), aligned):
        true = normalize_family(source["transducer_family"])
        predicted = FAMILIES[int(np.argmax(probabilities))]
        order = np.argsort(-probabilities)
        true_index = FAMILIES.index(true)
        rank = int(np.where(order == true_index)[0][0] + 1)
        margin = probabilities[true_index] - np.max(np.delete(probabilities, true_index))
        for family, probability in zip(FAMILIES, probabilities):
            rows.append({
                "pdb_id": str(source["pdb_id"]).upper(),
                "receptor_name": str(source["receptor_name"]),
                "true_family": true,
                "mode": "baseline",
                "model": "interface_sequence_logistic",
                "candidate_family": family,
                "compatibility_score": probability,
                "probability_like_score": probability,
                "predicted_family": predicted,
                "true_family_rank": rank,
                "true_minus_best_decoy": margin,
            })
    return pd.DataFrame(rows)


def _condition_splits(frame: pd.DataFrame, condition: str, folds: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    group_column = "receptor_name" if condition == "receptor_grouped" else condition
    if group_column not in frame.columns:
        return []
    groups = frame[group_column].astype(str).to_numpy()
    y = frame["transducer_family"].map(normalize_family).to_numpy()
    n_splits = min(int(folds), int(pd.Series(groups).nunique()))
    if n_splits < 2:
        return []
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    valid = []
    classes = set(FAMILIES)
    for train, test in splitter.split(frame, y, groups):
        if not classes.issubset(set(y[train])):
            continue
        if len(set(y[test])) < 2:
            continue
        valid.append((train, test))
    return valid


def _scramble_profile(profile: Mapping[str, str], rng: np.random.Generator) -> dict[str, str]:
    keys = list(profile)
    values = [profile[key] for key in keys]
    rng.shuffle(values)
    return dict(zip(keys, values))


def _scramble_template_contact_map(template: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    scrambled = template.copy()
    if len(scrambled) > 1:
        values = scrambled["galpha_reference_position"].to_numpy(copy=True)
        rng.shuffle(values)
        scrambled["galpha_reference_position"] = values
        for column in ("galpha_kd", "galpha_charge", "galpha_volume", "galpha_aromatic"):
            vals = scrambled[column].to_numpy(copy=True)
            rng.shuffle(vals)
            scrambled[column] = vals
    return scrambled


def _control_distribution(
    test: pd.DataFrame,
    templates: Mapping[str, pd.DataFrame],
    *,
    positions: Sequence[str],
    repeats: int,
    seed: int,
) -> pd.DataFrame:
    if repeats <= 0:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rows = []
    source_rows = test.to_dict(orient="records")
    for repeat in range(int(repeats)):
        # receptor-position scramble
        pred_rows = []
        for source in source_rows:
            profile = _scramble_profile(receptor_profile_from_row(source, positions, profile_mode="sequence"), rng)
            true = normalize_family(source["transducer_family"])
            scores = {
                family: score_receptor_against_template(profile, templates[family])["combined_physics"]
                for family in FAMILIES
            }
            probs = _softmax(np.array([scores[f] for f in FAMILIES]))
            for family, probability in zip(FAMILIES, probs):
                pred_rows.append({
                    "pdb_id": str(source["pdb_id"]).upper(),
                    "receptor_name": str(source["receptor_name"]),
                    "true_family": true,
                    "candidate_family": family,
                    "probability_like_score": probability,
                    "model": "receptor_position_scramble",
                    "true_family_rank": np.nan,
                    "true_minus_best_decoy": np.nan,
                })
        rows.append({
            "control": "receptor_position_scramble",
            "repeat": repeat,
            "macro_auc": _macro_metrics(pd.DataFrame(pred_rows))["macro_auc"],
        })
        # contact-map and G-alpha chemistry scramble
        scrambled_templates = {
            family: _scramble_template_contact_map(templates[family], rng)
            for family in FAMILIES
        }
        scored = score_receptor_rows(
            test, scrambled_templates, positions=positions, mode="contact_map_scramble"
        )
        scored = scored[scored["model"].eq("combined_physics")]
        rows.append({
            "control": "contact_map_and_galpha_scramble",
            "repeat": repeat,
            "macro_auc": _macro_metrics(scored)["macro_auc"],
        })
    return pd.DataFrame(rows)


def run_cross_validated_complementarity(
    *,
    frame: pd.DataFrame,
    contacts: pd.DataFrame,
    galpha_profiles: pd.DataFrame,
    output: str | Path,
    positions: Sequence[str] = DEFAULT_INTERFACE_POSITIONS,
    folds: int = 5,
    seed: int = 20272729,
    bootstraps: int = 500,
    scrambles: int = 50,
    k_neighbors: int = 5,
    quick: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = Path(output)
    checkpoints = out / "checkpoints" / "cross_validation"
    checkpoints.mkdir(parents=True, exist_ok=True)
    work = frame.copy()
    if "representative_for_receptor_transducer" in work.columns:
        work = work[work["representative_for_receptor_transducer"].fillna(False).astype(bool)]
    family_col = "transducer_family" if "transducer_family" in work.columns else "transducer_type"
    work["transducer_family"] = work[family_col].map(normalize_family)
    work = work[work["transducer_family"].isin(FAMILIES)].drop_duplicates("pdb_id").reset_index(drop=True)
    extracted = set(contacts["pdb_id"].astype(str).str.upper())
    work = work[work["pdb_id"].astype(str).str.upper().isin(extracted)].reset_index(drop=True)
    profile_lookup = _profile_table(work, positions)
    seq_cols = sequence_columns(work)
    conditions = ["receptor_grouped", "sequence_cluster_50", "sequence_cluster_40", "sequence_cluster_30"]
    all_predictions, all_controls, fold_audit = [], [], []
    for condition_index, condition in enumerate(conditions):
        splits = _condition_splits(work, condition, folds, seed + condition_index)
        print(f"[validation] {condition}: {len(splits)} valid folds", flush=True)
        for fold, (train_idx, test_idx) in enumerate(splits):
            checkpoint = checkpoints / f"{condition}_fold{fold}.tsv"
            control_checkpoint = checkpoints / f"{condition}_fold{fold}_controls.tsv"
            if checkpoint.exists():
                prediction = pd.read_csv(checkpoint, sep="\t")
                all_predictions.append(prediction)
                if control_checkpoint.exists():
                    all_controls.append(pd.read_csv(control_checkpoint, sep="\t"))
                fold_audit.append({"condition": condition, "fold": fold, "status": "cached"})
                continue
            train = work.iloc[train_idx].copy()
            test = work.iloc[test_idx].copy()
            train_pdbs = train["pdb_id"].astype(str).str.upper()
            test_pdbs = test["pdb_id"].astype(str).str.upper()
            templates = build_family_templates(
                contacts, galpha_profiles, profile_lookup, train_pdbs,
                minimum_edge_count=1 if quick else 2,
                minimum_edge_frequency=0.02 if quick else 0.03,
            )
            receptor_predictions = score_receptor_rows(
                test, templates, positions=positions, mode="receptor_alone_template",
                profile_mode="sequence",
            )
            structure_covered_predictions = score_receptor_rows(
                test, templates, positions=positions, mode="receptor_structure_covered_template",
                profile_mode="observed",
            )
            complex_predictions = score_complex_contact_map_decoys(
                test_pdbs, contacts, galpha_profiles, train_pdbs
            )
            identity_predictions = _identity_knn_predictions(
                train, test, sequence_cols=seq_cols, k=k_neighbors
            )
            interface_predictions = _interface_logistic_predictions(
                train, test, positions=positions
            )
            prediction = pd.concat(
                [receptor_predictions, structure_covered_predictions, complex_predictions, identity_predictions, interface_predictions],
                ignore_index=True,
            )
            prediction["condition"] = condition
            prediction["fold"] = fold
            write_table(prediction, checkpoint)
            all_predictions.append(prediction)
            controls = _control_distribution(
                test, templates, positions=positions,
                repeats=max(2, scrambles // max(len(splits), 1)),
                seed=seed + 1000 * condition_index + fold,
            )
            if not controls.empty:
                controls["condition"] = condition
                controls["fold"] = fold
                write_table(controls, control_checkpoint)
                all_controls.append(controls)
            for family, template in templates.items():
                write_table(template, out / "templates" / condition / f"fold_{fold}_{family.replace('/', '_')}.tsv")
            fold_audit.append({
                "condition": condition,
                "fold": fold,
                "status": "success",
                "n_train": len(train),
                "n_test": len(test),
                "n_train_receptors": train["receptor_name"].nunique(),
                "n_test_receptors": test["receptor_name"].nunique(),
                **{f"template_edges_{family}": len(templates[family]) for family in FAMILIES},
            })
    predictions = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame()
    controls = pd.concat(all_controls, ignore_index=True) if all_controls else pd.DataFrame()
    audit = pd.DataFrame(fold_audit)
    write_table(predictions, out / "predictions" / "cross_validated_predictions.tsv")
    write_table(controls, out / "controls" / "scramble_controls.tsv")
    write_table(audit, out / "audit" / "fold_audit.tsv")
    metric_rows = []
    if not predictions.empty:
        for (condition, mode, model), part in predictions.groupby(["condition", "mode", "model"]):
            metric = _macro_metrics(part)
            lo, hi = _bootstrap_metric(part, bootstraps=bootstraps, seed=seed + len(metric_rows))
            metric_rows.append({
                "condition": condition,
                "mode": mode,
                "model": model,
                **metric,
                "auc_ci_low": lo,
                "auc_ci_high": hi,
            })
    metrics = pd.DataFrame(metric_rows)
    write_table(metrics, out / "metrics" / "cross_validated_metrics.tsv")
    if not controls.empty and not metrics.empty:
        control_summary = []
        observed_lookup = metrics[
            (metrics["mode"].eq("receptor_alone_template"))
            & (metrics["model"].eq("combined_physics"))
        ].set_index("condition")["macro_auc"].to_dict()
        for (condition, control), part in controls.groupby(["condition", "control"]):
            observed = observed_lookup.get(condition, np.nan)
            values = part["macro_auc"].dropna().to_numpy()
            p = (1 + np.sum(values >= observed)) / (1 + len(values)) if len(values) and np.isfinite(observed) else np.nan
            control_summary.append({
                "condition": condition,
                "control": control,
                "observed_combined_auc": observed,
                "control_mean_auc": float(np.mean(values)) if len(values) else np.nan,
                "control_ci_low": float(np.quantile(values, 0.025)) if len(values) else np.nan,
                "control_ci_high": float(np.quantile(values, 0.975)) if len(values) else np.nan,
                "empirical_p": p,
                "n_repeats": len(values),
            })
        write_table(pd.DataFrame(control_summary), out / "controls" / "scramble_control_summary.tsv")
    return predictions, metrics, controls


def score_external_receptor_only(
    *,
    bound_frame: pd.DataFrame,
    contacts: pd.DataFrame,
    galpha_profiles: pd.DataFrame,
    receptor_only_table: pd.DataFrame,
    labels: pd.DataFrame | None,
    output: str | Path,
    definition: str,
    positions: Sequence[str] = DEFAULT_INTERFACE_POSITIONS,
    bootstraps: int = 1000,
    seed: int = 20272729,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = Path(output)
    family_col = "transducer_family" if "transducer_family" in bound_frame.columns else "transducer_type"
    bound = bound_frame.copy()
    if "representative_for_receptor_transducer" in bound.columns:
        bound = bound[bound["representative_for_receptor_transducer"].fillna(False).astype(bool)]
    bound["transducer_family"] = bound[family_col].map(normalize_family)
    bound = bound[bound["transducer_family"].isin(FAMILIES)].drop_duplicates("pdb_id")
    bound_profiles = _profile_table(bound, positions)
    templates = build_family_templates(
        contacts, galpha_profiles, bound_profiles, bound["pdb_id"],
        minimum_edge_count=2, minimum_edge_frequency=0.03,
    )
    sequence_predictions = score_receptor_rows(
        receptor_only_table, templates, positions=positions,
        true_family_column=None, mode=f"receptor_only_{definition}_sequence",
        profile_mode="sequence",
    )
    observed_predictions = score_receptor_rows(
        receptor_only_table, templates, positions=positions,
        true_family_column=None, mode=f"receptor_only_{definition}_structure_covered",
        profile_mode="observed",
    )
    predictions = pd.concat([sequence_predictions, observed_predictions], ignore_index=True)
    # Aggregate structure predictions to independent receptors.
    aggregate = (
        predictions.groupby(["receptor_name", "mode", "model", "candidate_family"], as_index=False)
        .agg(
            compatibility_score=("compatibility_score", "median"),
            probability_like_score=("probability_like_score", "median"),
            n_structures=("pdb_id", "nunique"),
            mean_template_coverage=("template_coverage", "mean"),
        )
    )
    if labels is not None and not labels.empty:
        keep = [column for column in ["receptor_name", *FAMILIES, "annotation_source", "notes"] if column in labels.columns]
        aggregate = aggregate.merge(labels[keep].drop_duplicates("receptor_name"), on="receptor_name", how="left")
    metric_rows = []
    if labels is not None and not labels.empty:
        for (mode, model, family), part in aggregate.groupby(["mode", "model", "candidate_family"]):
            if family not in part.columns:
                continue
            valid = part[family].notna() & part["probability_like_score"].notna()
            y = part.loc[valid, family].astype(float).to_numpy()
            score = part.loc[valid, "probability_like_score"].astype(float).to_numpy()
            if len(np.unique(y)) < 2:
                auc = np.nan
            else:
                auc = float(roc_auc_score(y, score))
            bootstrap_values = []
            if len(y) >= 4 and len(np.unique(y)) == 2 and bootstraps > 0:
                rng = np.random.default_rng(seed + len(metric_rows))
                indices = np.arange(len(y))
                for _ in range(int(bootstraps)):
                    sampled = rng.choice(indices, size=len(indices), replace=True)
                    ys = y[sampled]
                    if len(np.unique(ys)) < 2:
                        continue
                    bootstrap_values.append(float(roc_auc_score(ys, score[sampled])))
            metric_rows.append({
                "definition": definition,
                "mode": mode,
                "model": model,
                "family": family,
                "auc": auc,
                "auc_ci_low": float(np.quantile(bootstrap_values, 0.025)) if bootstrap_values else np.nan,
                "auc_ci_high": float(np.quantile(bootstrap_values, 0.975)) if bootstrap_values else np.nan,
                "n_bootstrap_valid": len(bootstrap_values),
                "n_receptors": int(valid.sum()),
                "n_positive": int(np.sum(y == 1)) if len(y) else 0,
                "n_negative": int(np.sum(y == 0)) if len(y) else 0,
            })
    metrics = pd.DataFrame(metric_rows)
    write_table(predictions, out / "predictions" / f"{definition}_structure_predictions.tsv")
    write_table(aggregate, out / "predictions" / f"{definition}_receptor_predictions.tsv")
    write_table(metrics, out / "metrics" / f"{definition}_functional_metrics.tsv")
    for family, template in templates.items():
        write_table(template, out / "templates" / "all_bound" / f"{family.replace('/', '_')}.tsv")
    return aggregate, metrics


def _make_summary_figures(metrics: pd.DataFrame, external: pd.DataFrame, output: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    condition_order = ["receptor_grouped", "sequence_cluster_50", "sequence_cluster_40", "sequence_cluster_30"]
    labels = ["Receptor grouped", "50%", "40%", "30%"]
    selected = [
        ("receptor_alone_template", "combined_physics", "Receptor–Gα chemistry"),
        ("baseline", "global_identity_knn", "Global identity kNN"),
        ("baseline", "interface_sequence_logistic", "Interface sequence"),
        ("receptor_structure_covered_template", "combined_physics", "Structure-covered chemistry"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    markers = ["o", "s", "^", "D"]
    for (mode, model, label), marker in zip(selected, markers):
        values, lows, highs = [], [], []
        for condition in condition_order:
            part = metrics[
                metrics["condition"].eq(condition)
                & metrics["mode"].eq(mode)
                & metrics["model"].eq(model)
            ]
            if len(part):
                row = part.iloc[0]
                values.append(float(row["macro_auc"]))
                lows.append(float(row["auc_ci_low"]))
                highs.append(float(row["auc_ci_high"]))
            else:
                values.append(np.nan); lows.append(np.nan); highs.append(np.nan)
        x = np.arange(len(condition_order))
        yerr = np.array([
            [v - lo if np.isfinite(v) and np.isfinite(lo) else np.nan for v, lo in zip(values, lows)],
            [hi - v if np.isfinite(v) and np.isfinite(hi) else np.nan for v, hi in zip(values, highs)],
        ])
        ax.errorbar(x, values, yerr=yerr, marker=marker, linewidth=1.5, capsize=2.5, label=label)
    ax.axhline(0.5, color="0.65", linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(labels)), labels)
    ax.set_ylabel("Macro AUC")
    ax.set_ylim(0.35, 1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(figure_dir / "complementarity_validation.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "complementarity_validation.png", dpi=400, bbox_inches="tight")
    plt.close(fig)

    if external is None or external.empty:
        return
    part = external[external["model"].eq("combined_physics")].copy()
    if part.empty:
        return
    def _external_label(row):
        definition = "Ligand-permissive" if row["definition"] == "ligand_permissive" else "Polymer-only"
        mode = "sequence" if str(row["mode"]).endswith("_sequence") else "structure-covered"
        return f"{definition}: {mode}"
    part["label"] = part.apply(_external_label, axis=1)
    labels_unique = list(dict.fromkeys(part["label"]))
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    x = np.arange(len(FAMILIES))
    width = 0.8 / max(len(labels_unique), 1)
    for i, label in enumerate(labels_unique):
        p = part[part["label"].eq(label)].set_index("family")
        values = [float(p.loc[f, "auc"]) if f in p.index else np.nan for f in FAMILIES]
        ax.bar(x - 0.4 + width / 2 + i * width, values, width=width, label=label)
    ax.axhline(0.5, color="0.65", linestyle="--", linewidth=1)
    ax.set_xticks(x, FAMILIES)
    ax.set_ylabel("Functional one-vs-rest AUC")
    ax.set_ylim(0, 1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(figure_dir / "receptor_only_compatibility.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "receptor_only_compatibility.png", dpi=400, bbox_inches="tight")
    plt.close(fig)

def _interpretation_gate(metrics: pd.DataFrame, external: pd.DataFrame | None = None) -> dict:
    def value(condition: str, mode: str, model: str) -> float:
        if metrics is None or metrics.empty or not {"condition", "mode", "model", "macro_auc"}.issubset(metrics.columns):
            return np.nan
        part = metrics[
            metrics["condition"].eq(condition)
            & metrics["mode"].eq(mode)
            & metrics["model"].eq(model)
        ]
        return float(part["macro_auc"].iloc[0]) if len(part) else np.nan
    rg = value("receptor_grouped", "receptor_alone_template", "combined_physics")
    id_rg = value("receptor_grouped", "baseline", "global_identity_knn")
    s50 = value("sequence_cluster_50", "receptor_alone_template", "combined_physics")
    id50 = value("sequence_cluster_50", "baseline", "global_identity_knn")
    s40 = value("sequence_cluster_40", "receptor_alone_template", "combined_physics")
    id40 = value("sequence_cluster_40", "baseline", "global_identity_knn")
    s30 = value("sequence_cluster_30", "receptor_alone_template", "combined_physics")
    if np.isfinite(s40) and np.isfinite(id40) and s40 >= 0.65 and s40 - id40 >= 0.03:
        outcome = "transferable_interface_complementarity_candidate"
    elif np.isfinite(rg) and rg >= 0.65 and (not np.isfinite(s40) or s40 < 0.62):
        outcome = "family_constrained_complementarity"
    elif np.isfinite(rg) and rg >= 0.58:
        outcome = "limited_support"
    else:
        outcome = "no_support"
    return {
        "outcome": outcome,
        "receptor_grouped_combined_auc": rg,
        "receptor_grouped_identity_knn_auc": id_rg,
        "sequence_cluster_50_combined_auc": s50,
        "sequence_cluster_50_identity_knn_auc": id50,
        "sequence_cluster_40_combined_auc": s40,
        "sequence_cluster_40_identity_knn_auc": id40,
        "sequence_cluster_30_combined_auc": s30,
        "interpretation_limits": [
            "Compatibility scores are not calibrated functional coupling probabilities.",
            "A positive bound-complex result requires independent receptor-only functional transfer.",
            "Low-identity analyses may have few valid folds and must be interpreted with their uncertainty.",
            "The external receptor-only set is small and multi-label; family-wise AUCs are exploratory.",
        ],
    }


def run_complementarity_audit(
    *,
    frame: str | Path,
    mapping: str | Path,
    structures: str | Path,
    output: str | Path,
    gprotein_contacts: str | Path | None = None,
    receptor_only_ligand: str | Path | None = None,
    receptor_only_polymer: str | Path | None = None,
    receptor_only_labels: str | Path | None = None,
    interface_positions: str | Path | None = None,
    workers: int = 1,
    folds: int = 5,
    seed: int = 20272729,
    bootstraps: int = 500,
    scrambles: int = 50,
    contact_cutoff: float = 4.5,
    quick: bool = False,
    quick_limit: int = 36,
    k_neighbors: int = 5,
) -> dict:
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    frame_df = read_table(frame)
    gp_contacts = read_table(gprotein_contacts) if gprotein_contacts and Path(gprotein_contacts).exists() else None
    if interface_positions and Path(interface_positions).exists():
        positions_df = read_table(interface_positions)
        column = next((c for c in ("generic_position", "position") if c in positions_df.columns), None)
        positions = tuple(positions_df[column].dropna().astype(str)) if column else DEFAULT_INTERFACE_POSITIONS
    else:
        positions = DEFAULT_INTERFACE_POSITIONS
    contacts, galpha_profiles, extraction_audit = extract_contact_dataset(
        frame=frame_df,
        mapping_path=mapping,
        structures_dir=structures,
        gprotein_contacts=gp_contacts,
        output=out,
        workers=workers,
        cutoff=contact_cutoff,
        quick_limit=quick_limit if quick else None,
    )
    predictions, metrics, controls = run_cross_validated_complementarity(
        frame=frame_df,
        contacts=contacts,
        galpha_profiles=galpha_profiles,
        output=out,
        positions=positions,
        folds=folds,
        seed=seed,
        bootstraps=min(bootstraps, 100) if quick else bootstraps,
        scrambles=min(scrambles, 12) if quick else scrambles,
        k_neighbors=k_neighbors,
        quick=quick,
    )
    external_metrics = []
    labels_df = read_table(receptor_only_labels) if receptor_only_labels and Path(receptor_only_labels).exists() else None
    for definition, path in (("ligand_permissive", receptor_only_ligand), ("polymer_only", receptor_only_polymer)):
        if path and Path(path).exists():
            table = read_table(path)
            _, ext_metrics = score_external_receptor_only(
                bound_frame=frame_df,
                contacts=contacts,
                galpha_profiles=galpha_profiles,
                receptor_only_table=table,
                labels=labels_df,
                output=out / "external_receptor_only",
                definition=definition,
                positions=positions,
                bootstraps=min(bootstraps, 100) if quick else bootstraps,
                seed=seed,
            )
            if not ext_metrics.empty:
                external_metrics.append(ext_metrics)
    external = pd.concat(external_metrics, ignore_index=True) if external_metrics else pd.DataFrame()
    if not external.empty:
        write_table(external, out / "external_receptor_only" / "metrics" / "all_functional_metrics.tsv")
    _make_summary_figures(metrics, external, out)
    gate = _interpretation_gate(metrics, external)
    write_json(out / "interpretation_gate.json", gate)
    manifest = {
        "status": "success" if not metrics.empty else "failed",
        "analysis": "gpcr_galpha_interface_complementarity",
        "version": "1.0.0",
        "quick": bool(quick),
        "seed": int(seed),
        "contact_cutoff_angstrom": float(contact_cutoff),
        "n_input_rows": int(len(frame_df)),
        "n_contact_structures_success": int(extraction_audit["status"].eq("success").sum()) if not extraction_audit.empty else 0,
        "n_contact_structures_failed": int(extraction_audit["status"].ne("success").sum()) if not extraction_audit.empty else 0,
        "n_contact_pairs": int(len(contacts)),
        "n_galpha_profile_rows": int(len(galpha_profiles)),
        "n_prediction_rows": int(len(predictions)),
        "n_metric_rows": int(len(metrics)),
        "n_control_rows": int(len(controls)),
        "n_external_metric_rows": int(len(external)),
        "elapsed_seconds": time.time() - started,
        "interpretation_outcome": gate["outcome"],
        "inputs": {
            "frame": str(frame),
            "mapping": str(mapping),
            "structures": str(structures),
            "gprotein_contacts": str(gprotein_contacts or ""),
            "receptor_only_ligand": str(receptor_only_ligand or ""),
            "receptor_only_polymer": str(receptor_only_polymer or ""),
            "receptor_only_labels": str(receptor_only_labels or ""),
        },
    }
    write_json(out / "RUN_MANIFEST.json", manifest)
    summary = [
        "# GPCR–Gα interface complementarity audit",
        "",
        f"Outcome: **{gate['outcome']}**",
        "",
        "The package evaluates actual complex contact chemistry and a separate receptor-alone mode built from training-only Gα-family templates.",
        "Scores are compatibility descriptors and are not calibrated functional coupling probabilities.",
        "",
        "## Key validation values",
        "",
        f"- Receptor-grouped combined receptor-alone AUC: {gate['receptor_grouped_combined_auc']:.3f}" if np.isfinite(gate['receptor_grouped_combined_auc']) else "- Receptor-grouped combined receptor-alone AUC: unavailable",
        f"- Receptor-grouped identity-kNN AUC: {gate['receptor_grouped_identity_knn_auc']:.3f}" if np.isfinite(gate['receptor_grouped_identity_knn_auc']) else "- Receptor-grouped identity-kNN AUC: unavailable",
        f"- 50% cluster combined receptor-alone AUC: {gate['sequence_cluster_50_combined_auc']:.3f}" if np.isfinite(gate['sequence_cluster_50_combined_auc']) else "- 50% cluster combined receptor-alone AUC: unavailable",
        f"- 40% cluster combined receptor-alone AUC: {gate['sequence_cluster_40_combined_auc']:.3f}" if np.isfinite(gate['sequence_cluster_40_combined_auc']) else "- 40% cluster combined receptor-alone AUC: unavailable",
        f"- 30% cluster combined receptor-alone AUC: {gate['sequence_cluster_30_combined_auc']:.3f}" if np.isfinite(gate['sequence_cluster_30_combined_auc']) else "- 30% cluster combined receptor-alone AUC: unavailable",
        "",
        "See `metrics/`, `controls/`, `predictions/`, and `external_receptor_only/` for complete results.",
    ]
    (out / "run_summary.md").write_text("\n".join(summary) + "\n")
    return manifest


def score_receptor_compatibility_from_templates(
    *,
    receptor_table: str | Path,
    template_dir: str | Path,
    output: str | Path,
    interface_positions: str | Path | None = None,
    profile_mode: str = "sequence",
) -> pd.DataFrame:
    frame = read_table(receptor_table)
    if interface_positions and Path(interface_positions).exists():
        positions_df = read_table(interface_positions)
        column = "generic_position" if "generic_position" in positions_df.columns else "position"
        positions = tuple(positions_df[column].dropna().astype(str))
    else:
        positions = DEFAULT_INTERFACE_POSITIONS
    template_dir = Path(template_dir)
    templates = {}
    for family in FAMILIES:
        candidates = [
            template_dir / f"{family.replace('/', '_')}.tsv",
            template_dir / f"{family}.tsv",
        ]
        path = next((candidate for candidate in candidates if candidate.exists()), None)
        if path is None:
            raise FileNotFoundError(f"Missing template for {family} in {template_dir}")
        templates[family] = read_table(path)
    predictions = score_receptor_rows(
        frame, templates, positions=positions, true_family_column=None,
        mode="standalone_receptor_compatibility", profile_mode=profile_mode,
    )
    if predictions.empty or predictions["profile_positions"].fillna(0).max() <= 0:
        raise ValueError(
            "The receptor table contains none of the required GPCR generic-position residue columns. "
            "Use a publication-control or receptor-only structure table with gpcrdb_<position>_expected_aa "
            "and/or residue_<position>_identity columns."
        )
    write_table(predictions, Path(output) / "receptor_compatibility.tsv")
    write_json(Path(output) / "receptor_compatibility_manifest.json", {
        "warning": "Compatibility scores are similarities to bound-complex-derived templates, not validated functional coupling probabilities.",
        "n_input_rows": len(frame),
        "n_output_rows": len(predictions),
    })
    return predictions
