from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
import yaml

from .features import StructureTask, extract_structure_features, feature_dictionary
from .modelling import (
    PRIMARY_CLASSES,
    build_model_specs,
    cross_validate_model,
    grouped_prediction_permutation,
    paired_model_difference,
    strict_receptor_transfer,
)


DEFAULT_PATHS = {
    "bound_frame": "results/microswitches/publication_controls/publication_control_analysis_frame.csv",
    "mapping": "results/microswitches/measurements/generic_number_mapping.csv",
    "structures": "data/structures",
    "receptor_only_primary": "data/manual/receptor_only_ligand_permissive_structures.csv",
    "receptor_only_conservative": "data/manual/strict_polymer_only_structures.csv",
    "receptor_only_labels": "data/manual/receptor_only_receptor_level_labels.csv",
}

PRIMARY_COMPARISONS = (
    ("missingness_only", "contact_network"),
    ("missingness_only", "surface_electrostatics"),
    ("missingness_only", "mechanical_susceptibility"),
    ("missingness_only", "all_new_static"),
    ("existing_geometry", "all_new_static"),
    ("sequence_only", "all_new_static"),
    ("sequence_only", "sequence_plus_new_static"),
    ("existing_geometry", "existing_geometry_plus_new_static"),
)


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    temporary.replace(path)


def atomic_csv(path: Path, frame: pd.DataFrame, sep: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, sep=sep)
    temporary.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, low_memory=False)


def _path(project: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project / path


def _resolve_bound_frame(project: Path, configured: Path) -> Path:
    """Resolve the bound analysis frame across supported repository layouts."""
    candidates = [
        configured,
        project / "results/microswitches/publication_controls/publication_control_analysis_frame.csv",
        project / "results/microswitches/coupling_analysis_representatives.csv",
        project / "results/microswitches/coupling_specificity/data/coupling_analysis_representatives.csv",
    ]
    seen: set[Path] = set()
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        if candidate.is_file():
            return candidate
    attempted = "\n  ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        "Could not locate a bound receptor analysis frame. Tried:\n  " + attempted
    )


def load_config(project: Path, config_path: Path | None) -> dict[str, Any]:
    path = config_path or project / "config/receptor_only_dynamic_selectivity.yaml"
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    payload["_config_path"] = str(path)
    return payload


def select_quick_bound(frame: pd.DataFrame, per_class: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    selected = []
    for label in PRIMARY_CLASSES:
        part = frame.loc[frame["transducer_family"].astype(str).eq(label)].copy()
        receptors = part["receptor_name"].astype(str).unique()
        rng.shuffle(receptors)
        chosen_receptors = set(receptors[:per_class])
        selected.append(part.loc[part["receptor_name"].astype(str).isin(chosen_receptors)])
    result = pd.concat(selected, ignore_index=True).drop_duplicates("pdb_id")
    return result


def build_manifests(project: Path, config: Mapping[str, Any], mode: str, seed: int) -> dict[str, pd.DataFrame]:
    paths = {key: _path(project, config.get("paths", {}).get(key, default)) for key, default in DEFAULT_PATHS.items()}
    paths["bound_frame"] = _resolve_bound_frame(project, paths["bound_frame"])
    bound = _read_table(paths["bound_frame"])
    bound = bound.loc[bound["transducer_family"].astype(str).isin(PRIMARY_CLASSES)].copy()
    primary = _read_table(paths["receptor_only_primary"])
    conservative = _read_table(paths["receptor_only_conservative"])
    labels = _read_table(paths["receptor_only_labels"])

    for frame in (bound, primary, conservative):
        frame["pdb_id"] = frame["pdb_id"].astype(str).str.upper()
        if "preferred_chain" not in frame and "chain_id" in frame:
            frame["preferred_chain"] = frame["chain_id"]
        if "chain_id" not in frame and "preferred_chain" in frame:
            frame["chain_id"] = frame["preferred_chain"]

    if mode == "quick":
        per_class = int(config.get("quick", {}).get("bound_receptors_per_class", 8))
        bound = select_quick_bound(bound, per_class, seed)
        quick_receptors = int(config.get("quick", {}).get("receptor_only_receptors", 8))
        receptor_names = sorted(primary["receptor_name"].astype(str).unique())[:quick_receptors]
        primary = primary.loc[primary["receptor_name"].astype(str).isin(receptor_names)].copy()
        conservative = conservative.loc[conservative["receptor_name"].astype(str).isin(receptor_names)].copy()

    return {
        "bound": bound.reset_index(drop=True),
        "primary": primary.reset_index(drop=True),
        "conservative": conservative.reset_index(drop=True),
        "labels": labels.reset_index(drop=True),
        "mapping_path": pd.DataFrame({"path": [str(paths["mapping"])]}),
        "structures_path": pd.DataFrame({"path": [str(paths["structures"])]}),
    }


def run_input_audit(project: Path, output_root: Path, config: Mapping[str, Any], manifests: Mapping[str, pd.DataFrame]) -> dict[str, Any]:
    audit_dir = output_root / "input_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    mapping_path = Path(manifests["mapping_path"].iloc[0, 0])
    structures_dir = Path(manifests["structures_path"].iloc[0, 0])
    mapping = _read_table(mapping_path)
    structure_files = {path.stem.upper() for path in structures_dir.glob("*.cif")}

    trajectory_extensions = {".xtc", ".trr", ".dcd", ".nc", ".netcdf"}
    topology_extensions = {".tpr", ".gro", ".psf", ".prmtop", ".parm7"}
    trajectories = [path for path in project.rglob("*") if path.is_file() and path.suffix.lower() in trajectory_extensions]
    topologies = [path for path in project.rglob("*") if path.is_file() and path.suffix.lower() in topology_extensions]

    inventory_rows = []
    for role in ("bound", "primary", "conservative"):
        frame = manifests[role]
        for row in frame.itertuples(index=False):
            pdb_id = str(row.pdb_id).upper()
            part = mapping.loc[mapping["pdb_id"].astype(str).str.upper().eq(pdb_id)]
            inventory_rows.append(
                {
                    "dataset_role": role,
                    "pdb_id": pdb_id,
                    "receptor_name": str(row.receptor_name),
                    "chain_id": str(getattr(row, "preferred_chain", "") or getattr(row, "chain_id", "") or ""),
                    "transducer_family": str(getattr(row, "transducer_family", "") or ""),
                    "structure_file_present": pdb_id in structure_files,
                    "generic_mapping_rows": int(len(part)),
                    "generic_positions": int(part["generic_number"].dropna().astype(str).nunique()) if not part.empty else 0,
                    "functional_annotation_available": bool(getattr(row, "functional_annotation_available", False)) if role != "bound" else False,
                }
            )
    inventory = pd.DataFrame(inventory_rows)
    atomic_csv(audit_dir / "structure_inventory.tsv", inventory)

    feasibility = pd.DataFrame(
        [
            {
                "analysis": "contact_and_sidechain_networks",
                "status": "immediately_feasible",
                "reason": "Experimental mmCIF receptor chains and generic-position mappings are present.",
            },
            {
                "analysis": "intracellular_surface_chemistry",
                "status": "immediately_feasible",
                "reason": "Biopython Shrake-Rupley receptor-chain SASA can be calculated without transducer atoms.",
            },
            {
                "analysis": "approximate_coulombic_electrostatics",
                "status": "immediately_feasible",
                "reason": "Transparent charged-residue Coulombic descriptors can be computed; these are not PB potentials.",
            },
            {
                "analysis": "anisotropic_network_mechanical_susceptibility",
                "status": "immediately_feasible",
                "reason": "Generic-position receptor CA coordinates support a coarse-grained ANM.",
            },
            {
                "analysis": "poisson_boltzmann_electrostatics",
                "status": "not_in_primary_implementation",
                "reason": "Consistent protonation and APBS/PDB2PQR preparation were not available for all structures.",
            },
            {
                "analysis": "trajectory_ensemble_accessibility",
                "status": "requires_new_simulations",
                "reason": f"Found {len(trajectories)} trajectory and {len(topologies)} topology files.",
            },
            {
                "analysis": "hydration_kinetics_and_water_wires",
                "status": "requires_new_simulations",
                "reason": "No validated receptor-only trajectory ensemble is present.",
            },
            {
                "analysis": "markov_state_models_and_mfpt",
                "status": "requires_new_simulations",
                "reason": "No validated receptor-only trajectory ensemble is present.",
            },
        ]
    )
    atomic_csv(audit_dir / "feature_feasibility.tsv", feasibility)

    label_coverage = []
    for role in ("primary", "conservative"):
        frame = manifests[role]
        for label in PRIMARY_CLASSES:
            values = pd.to_numeric(frame.get(label, pd.Series(dtype=float)), errors="coerce")
            receptor_frame = frame.assign(_label=values).drop_duplicates("receptor_name")
            known = receptor_frame["_label"].notna()
            label_coverage.append(
                {
                    "definition": role,
                    "target": label,
                    "structures": int(len(frame)),
                    "receptors": int(frame["receptor_name"].nunique()),
                    "known_receptors": int(known.sum()),
                    "positive_receptors": int((receptor_frame.loc[known, "_label"] == 1).sum()),
                    "negative_receptors": int((receptor_frame.loc[known, "_label"] == 0).sum()),
                }
            )
    atomic_csv(audit_dir / "label_coverage.tsv", pd.DataFrame(label_coverage))

    source_root = Path(__file__).resolve().parent
    source_hashes = {
        path.name: sha256_file(path)
        for path in sorted(source_root.glob("*.py"))
    }
    payload = {
        "analysis": "receptor_only_dynamic_selectivity_static_first_pass",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_dir": str(project),
        "config_path": config.get("_config_path"),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "bound_structures": int(len(manifests["bound"])),
        "bound_receptors": int(manifests["bound"]["receptor_name"].nunique()),
        "receptor_only_primary_structures": int(len(manifests["primary"])),
        "receptor_only_primary_receptors": int(manifests["primary"]["receptor_name"].nunique()),
        "receptor_only_conservative_structures": int(len(manifests["conservative"])),
        "receptor_only_conservative_receptors": int(manifests["conservative"]["receptor_name"].nunique()),
        "mapping_rows": int(len(mapping)),
        "structure_files": int(len(structure_files)),
        "trajectory_files": [str(path.relative_to(project)) for path in trajectories],
        "topology_files": [str(path.relative_to(project)) for path in topologies],
        "source_hashes": source_hashes,
        "critical_definition": "Receptor-only permits ligands but excludes transducer, auxiliary binder and fusion protein entities.",
        "transducer_atom_rule": "Every new feature is extracted from the selected receptor chain only.",
    }
    atomic_json(audit_dir / "input_audit.json", payload)
    (audit_dir / "software_versions.txt").write_text(
        "\n".join(
            [
                f"python={sys.version.split()[0]}",
                f"numpy={np.__version__}",
                f"pandas={pd.__version__}",
                f"scipy={scipy.__version__}",
                f"scikit-learn={sklearn.__version__}",
                f"matplotlib={matplotlib.__version__}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    audit_rows = pd.DataFrame([{"item": key, "value": json.dumps(value) if isinstance(value, (dict, list)) else value} for key, value in payload.items()])
    atomic_csv(audit_dir / "input_audit.tsv", audit_rows)
    markdown = [
        "# Receptor-only dynamic-selectivity input audit",
        "",
        f"- Bound structures: {payload['bound_structures']} from {payload['bound_receptors']} receptors",
        f"- Primary receptor-only structures: {payload['receptor_only_primary_structures']} from {payload['receptor_only_primary_receptors']} receptors",
        f"- Conservative polymer-only structures: {payload['receptor_only_conservative_structures']} from {payload['receptor_only_conservative_receptors']} receptors",
        f"- Trajectory files: {len(trajectories)}",
        f"- Topology files: {len(topologies)}",
        "",
        "The executable first pass is restricted to receptor-chain contact networks, surface chemistry, approximate Coulombic electrostatics and generic-position anisotropic-network susceptibility.",
        "Trajectory ensembles, hydration kinetics, MSMs and accessibility free-energy proxies are not executed because the required inputs are absent.",
    ]
    (audit_dir / "input_audit.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")
    return payload


def _extract_worker(payload: Mapping[str, Any]) -> dict[str, Any]:
    task = StructureTask(**payload["task"])
    started = time.time()
    try:
        row = extract_structure_features(
            task,
            payload["mapping_path"],
            payload["structures_dir"],
            sasa_points=int(payload["sasa_points"]),
            contact_cutoff=float(payload["contact_cutoff"]),
            anm_cutoff=float(payload["anm_cutoff"]),
            anm_modes=int(payload["anm_modes"]),
        )
        row["status"] = "ok"
        row["runtime_seconds"] = time.time() - started
        return row
    except Exception as error:  # noqa: BLE001
        return {
            **payload["task"],
            "status": "failed",
            "runtime_seconds": time.time() - started,
            "error": repr(error),
            "traceback": traceback.format_exc(),
        }


def extract_all_features(
    output_root: Path,
    manifests: Mapping[str, pd.DataFrame],
    config: Mapping[str, Any],
    workers: int,
    resume: bool,
    mode: str,
) -> pd.DataFrame:
    feature_dir = output_root / "features"
    checkpoint_dir = feature_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    mapping_path = str(manifests["mapping_path"].iloc[0, 0])
    structures_dir = str(manifests["structures_path"].iloc[0, 0])

    union = pd.concat(
        [
            manifests["bound"].assign(dataset_role="bound"),
            manifests["primary"].assign(dataset_role="receptor_only"),
        ],
        ignore_index=True,
    )
    union = union.drop_duplicates("pdb_id", keep="first")
    task_rows = []
    for row in union.itertuples(index=False):
        task = {
            "pdb_id": str(row.pdb_id).upper(),
            "receptor_name": str(row.receptor_name),
            "chain_id": str(getattr(row, "preferred_chain", "") or getattr(row, "chain_id", "") or ""),
            "dataset_role": str(getattr(row, "dataset_role", "")),
            "transducer_family": str(getattr(row, "transducer_family", "") or ""),
        }
        checkpoint = checkpoint_dir / f"{task['pdb_id']}_{task['chain_id'] or 'NA'}.json"
        if resume and checkpoint.exists():
            try:
                existing = json.loads(checkpoint.read_text(encoding="utf-8"))
                if existing.get("status") == "ok":
                    continue
            except json.JSONDecodeError:
                pass
        extraction = config.get("feature_extraction", {})
        task_rows.append(
            {
                "task": task,
                "checkpoint": str(checkpoint),
                "mapping_path": mapping_path,
                "structures_dir": structures_dir,
                "sasa_points": extraction.get("sasa_points_quick" if mode == "quick" else "sasa_points_full", 60 if mode == "quick" else 100),
                "contact_cutoff": extraction.get("contact_cutoff_angstrom", 8.0),
                "anm_cutoff": extraction.get("anm_cutoff_angstrom", 15.0),
                "anm_modes": extraction.get("anm_modes_quick" if mode == "quick" else "anm_modes_full", 8 if mode == "quick" else 12),
            }
        )

    if task_rows:
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(workers, len(task_rows))) as executor:
            futures = {executor.submit(_extract_worker, payload): payload for payload in task_rows}
            for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                result = future.result()
                payload = futures[future]
                atomic_json(Path(payload["checkpoint"]), result)
                print(
                    f"[features {index}/{len(task_rows)}] {result.get('pdb_id')} "
                    f"status={result.get('status')} runtime_s={result.get('runtime_seconds', np.nan):.1f}",
                    flush=True,
                )

    rows = []
    for path in sorted(checkpoint_dir.glob("*.json")):
        try:
            rows.append(json.loads(path.read_text(encoding="utf-8")))
        except json.JSONDecodeError:
            continue
    audit = pd.DataFrame(rows)
    atomic_csv(output_root / "logs/excluded_cases.tsv", audit.loc[audit.get("status", pd.Series(dtype=str)).ne("ok")])
    features = audit.loc[audit.get("status", pd.Series(dtype=str)).eq("ok")].copy()
    if features.empty:
        raise RuntimeError("No structure produced a complete feature checkpoint")
    drop_columns = [column for column in ("status", "traceback") if column in features]
    features = features.drop(columns=drop_columns)
    atomic_csv(feature_dir / "receptor_level_features.tsv", features)
    atomic_csv(
        feature_dir / "trajectory_level_features.tsv",
        pd.DataFrame(
            [
                {
                    "status": "not_run_missing_validated_receptor_only_trajectories",
                    "reason": "No supported trajectory/topology pairs were present in the supplied repository.",
                }
            ]
        ),
    )
    atomic_csv(feature_dir / "feature_dictionary.tsv", feature_dictionary(features))
    return features


def merge_features(manifest: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    left = manifest.copy()
    left["pdb_id"] = left["pdb_id"].astype(str).str.upper()
    right = features.drop(columns=[column for column in ("receptor_name", "chain_id", "transducer_family", "dataset_role") if column in features]).copy()
    right["pdb_id"] = right["pdb_id"].astype(str).str.upper()
    return left.merge(right.drop_duplicates("pdb_id"), on="pdb_id", how="left")


def run_models(
    project: Path,
    output_root: Path,
    config: Mapping[str, Any],
    manifests: Mapping[str, pd.DataFrame],
    features: pd.DataFrame,
    mode: str,
    seed: int,
    resume: bool,
) -> dict[str, Any]:
    model_dir = output_root / "metrics"
    prediction_dir = output_root / "predictions"
    coefficient_dir = output_root / "coefficients"
    bootstrap_dir = output_root / "bootstraps"
    permutation_dir = output_root / "permutations"
    fold_dir = output_root / "fold_manifests"
    model_checkpoint_dir = output_root / "model_checkpoints"
    for directory in (model_dir, prediction_dir, coefficient_dir, bootstrap_dir, permutation_dir, fold_dir, model_checkpoint_dir):
        directory.mkdir(parents=True, exist_ok=True)

    bound = merge_features(manifests["bound"], features)
    primary = merge_features(manifests["primary"], features)
    conservative = merge_features(manifests["conservative"], features)
    prepared, specs, availability = build_model_specs(bound, project / "config")
    atomic_csv(output_root / "feature_dictionary/model_feature_availability.tsv", availability)

    validation = config.get("validation", {})
    quick = config.get("quick", {})
    full = config.get("full", {})
    settings = quick if mode == "quick" else full
    outer_folds = int(settings.get("outer_folds", 3 if mode == "quick" else 5))
    inner_folds = int(settings.get("inner_folds", 2 if mode == "quick" else 3))
    bootstrap_iterations = int(settings.get("bootstrap_iterations", 100 if mode == "quick" else 1000))
    permutation_iterations = int(settings.get("permutation_iterations", 50 if mode == "quick" else 1000))
    c_grid = [float(value) for value in validation.get("C_grid", [0.03, 0.1, 0.3, 1.0, 3.0, 10.0])]
    minimum_coverage = float(validation.get("minimum_training_feature_coverage", 0.60))
    max_iter = int(validation.get("max_iter", 50000))

    model_names = list(
        settings.get(
            "models",
            [
                "existing_geometry",
                "sequence_only",
                "missingness_only",
                "contact_network",
                "surface_electrostatics",
                "mechanical_susceptibility",
                "all_new_static",
                "sequence_plus_new_static",
                "existing_geometry_plus_new_static",
            ],
        )
    )

    performance_rows = []
    prediction_frames = []
    audit_frames = []
    coefficient_frames = []
    bootstrap_frames = []
    permutation_rows = []

    for index, model_name in enumerate(model_names, start=1):
        print(f"[model {index}/{len(model_names)}] {model_name}", flush=True)
        checkpoint = model_checkpoint_dir / model_name
        performance_path = checkpoint / "performance.json"
        predictions_path = checkpoint / "predictions.tsv"
        audit_path = checkpoint / "audit.tsv"
        coefficients_path = checkpoint / "coefficients.tsv"
        bootstrap_path = checkpoint / "bootstrap.tsv"
        permutation_path = checkpoint / "permutation.json"
        loaded = False
        if resume and all(path.exists() for path in (performance_path, predictions_path, audit_path, coefficients_path, permutation_path)):
            try:
                result = {
                    "performance": json.loads(performance_path.read_text(encoding="utf-8")),
                    "predictions": pd.read_csv(predictions_path, sep="\t", low_memory=False),
                    "audit": pd.read_csv(audit_path, sep="\t", low_memory=False),
                    "coefficients": pd.read_csv(coefficients_path, sep="\t", low_memory=False),
                    "bootstrap": pd.read_csv(bootstrap_path, sep="\t", low_memory=False) if bootstrap_path.exists() and bootstrap_path.stat().st_size > 1 else pd.DataFrame(),
                }
                permutation_record = json.loads(permutation_path.read_text(encoding="utf-8"))
                loaded = True
                print(f"[model {model_name}] reused checkpoint", flush=True)
            except Exception:
                loaded = False
        if not loaded:
            result = cross_validate_model(
                prepared,
                specs[model_name],
                seed=seed,
                outer_folds=outer_folds,
                inner_folds=inner_folds,
                c_grid=c_grid,
                bootstrap_iterations=bootstrap_iterations,
                minimum_coverage=minimum_coverage,
                max_iter=max_iter,
            )
            permutation_record = {
                "experiment": model_name,
                **grouped_prediction_permutation(
                    result["predictions"],
                    iterations=permutation_iterations,
                    seed=seed + 30000 + index,
                ),
            }
            checkpoint.mkdir(parents=True, exist_ok=True)
            atomic_json(performance_path, result["performance"])
            atomic_csv(predictions_path, result["predictions"])
            atomic_csv(audit_path, result["audit"])
            atomic_csv(coefficients_path, result["coefficients"])
            atomic_csv(bootstrap_path, result["bootstrap"].assign(experiment=model_name) if not result["bootstrap"].empty else pd.DataFrame())
            atomic_json(permutation_path, permutation_record)
        performance_rows.append(result["performance"])
        prediction_frames.append(result["predictions"])
        audit_frames.append(result["audit"])
        coefficient_frames.append(result["coefficients"])
        if not result["bootstrap"].empty:
            bootstrap_frames.append(result["bootstrap"].assign(experiment=model_name) if "experiment" not in result["bootstrap"] else result["bootstrap"])
        permutation_rows.append(permutation_record)

    performance = pd.DataFrame(performance_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    audits = pd.concat(audit_frames, ignore_index=True) if audit_frames else pd.DataFrame()
    coefficients = pd.concat(coefficient_frames, ignore_index=True) if coefficient_frames else pd.DataFrame()
    bootstraps = pd.concat(bootstrap_frames, ignore_index=True) if bootstrap_frames else pd.DataFrame()
    permutations = pd.DataFrame(permutation_rows)

    atomic_csv(model_dir / "model_metrics.tsv", performance)
    atomic_csv(prediction_dir / "out_of_fold_predictions.tsv", predictions)
    atomic_csv(fold_dir / "fold_convergence_and_leakage_audit.tsv", audits)
    atomic_csv(coefficient_dir / "fold_coefficients.tsv", coefficients)
    stability_dir = output_root / "feature_stability"
    stability_dir.mkdir(parents=True, exist_ok=True)
    if not coefficients.empty:
        stability = (
            coefficients.assign(
                coefficient_sign=lambda frame: np.sign(frame["coefficient"]),
                nonzero=lambda frame: frame["absolute_coefficient"].gt(1e-8).astype(float),
            )
            .groupby(["experiment", "class_label", "transformed_feature"], as_index=False)
            .agg(
                folds_observed=("fold", "nunique"),
                mean_coefficient=("coefficient", "mean"),
                mean_absolute_coefficient=("absolute_coefficient", "mean"),
                nonzero_fraction=("nonzero", "mean"),
                mean_sign=("coefficient_sign", "mean"),
            )
        )
        stability["sign_stability"] = stability["mean_sign"].abs()
    else:
        stability = pd.DataFrame()
    atomic_csv(stability_dir / "feature_selection_stability.tsv", stability)
    fold_assignments = (
        predictions.loc[predictions["experiment"].eq(model_names[0]), ["sample_id", "pdb_id", "receptor_name", "transducer_family", "outer_fold"]]
        .drop_duplicates()
        .sort_values(["outer_fold", "receptor_name", "pdb_id"])
    )
    atomic_csv(fold_dir / "fold_assignments.tsv", fold_assignments)
    atomic_csv(bootstrap_dir / "receptor_cluster_bootstraps.tsv", bootstraps)
    atomic_csv(permutation_dir / "grouped_permutation_results.tsv", permutations)

    paired_rows = []
    for index, (baseline, candidate) in enumerate(PRIMARY_COMPARISONS):
        if baseline not in model_names or candidate not in model_names:
            continue
        paired_rows.append(
            paired_model_difference(
                predictions,
                baseline,
                candidate,
                iterations=bootstrap_iterations,
                seed=seed + 50000 + index,
            )
        )
    paired = pd.DataFrame(paired_rows)
    atomic_csv(model_dir / "paired_model_differences.tsv", paired)

    transfer_names = list(
        settings.get(
            "strict_transfer_models",
            [
                "existing_geometry",
                "sequence_only",
                "contact_network",
                "surface_electrostatics",
                "mechanical_susceptibility",
                "all_new_static",
                "sequence_plus_new_static",
                "existing_geometry_plus_new_static",
            ],
        )
    )
    strict_dir = output_root / "strict_receptor_only"
    strict_dir.mkdir(parents=True, exist_ok=True)
    transfer_outputs = []
    transfer_prediction_frames = []
    transfer_receptor_frames = []
    transfer_audit_frames = []
    for definition, frame in (("ligand_permissive", primary), ("polymer_only", conservative)):
        definition_dir = strict_dir / "checkpoints" / definition
        paths = {
            "metrics": definition_dir / "metrics.tsv",
            "structure_predictions": definition_dir / "structure_predictions.tsv",
            "receptor_predictions": definition_dir / "receptor_predictions.tsv",
            "audit": definition_dir / "audit.tsv",
        }
        loaded = False
        if resume and all(path.exists() for path in paths.values()):
            try:
                transfer = {key: pd.read_csv(path, sep="\t", low_memory=False) for key, path in paths.items()}
                loaded = True
                print(f"[strict transfer {definition}] reused checkpoint", flush=True)
            except Exception:
                loaded = False
        if not loaded:
            transfer = strict_receptor_transfer(
                prepared,
                frame,
                specs,
                transfer_names,
                seed=seed,
                inner_folds=inner_folds,
                c_grid=c_grid,
                max_iter=max_iter,
                minimum_coverage=minimum_coverage,
                bootstrap_iterations=bootstrap_iterations,
                definition=definition,
                fixed_c=float(validation.get("strict_transfer_C", 1.0)),
            )
            definition_dir.mkdir(parents=True, exist_ok=True)
            for key, path in paths.items():
                atomic_csv(path, transfer[key])
        if not transfer["metrics"].empty:
            transfer_outputs.append(transfer["metrics"])
        if not transfer["structure_predictions"].empty:
            transfer_prediction_frames.append(transfer["structure_predictions"])
        if not transfer["receptor_predictions"].empty:
            transfer_receptor_frames.append(transfer["receptor_predictions"])
        if not transfer["audit"].empty:
            transfer_audit_frames.append(transfer["audit"])

    transfer_metrics = pd.concat(transfer_outputs, ignore_index=True) if transfer_outputs else pd.DataFrame()
    transfer_structure = pd.concat(transfer_prediction_frames, ignore_index=True) if transfer_prediction_frames else pd.DataFrame()
    transfer_receptor = pd.concat(transfer_receptor_frames, ignore_index=True) if transfer_receptor_frames else pd.DataFrame()
    transfer_audit = pd.concat(transfer_audit_frames, ignore_index=True) if transfer_audit_frames else pd.DataFrame()
    atomic_csv(strict_dir / "functional_transfer_metrics.tsv", transfer_metrics)
    atomic_csv(strict_dir / "structure_predictions.tsv", transfer_structure)
    atomic_csv(strict_dir / "receptor_predictions.tsv", transfer_receptor)
    atomic_csv(strict_dir / "transfer_training_audit.tsv", transfer_audit)

    return {
        "bound": prepared,
        "primary": primary,
        "conservative": conservative,
        "specs": specs,
        "performance": performance,
        "predictions": predictions,
        "audits": audits,
        "paired": paired,
        "permutations": permutations,
        "transfer_metrics": transfer_metrics,
        "transfer_receptor_predictions": transfer_receptor,
    }


def make_figures(output_root: Path, results: Mapping[str, Any], dpi: int = 600) -> None:
    figure_dir = output_root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    performance = results["performance"].copy()
    performance = performance.loc[performance["experiment"].ne("new_feature_availability_only")]
    performance = performance.sort_values("macro_roc_auc", ascending=True)

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    y = np.arange(len(performance))
    values = performance["macro_roc_auc"].to_numpy(float)
    low = performance["macro_roc_auc_ci_low"].to_numpy(float)
    high = performance["macro_roc_auc_ci_high"].to_numpy(float)
    ax.errorbar(values, y, xerr=np.vstack([values - low, high - values]), fmt="o", capsize=4)
    ax.axvline(0.5, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(performance["experiment"].str.replace("_", " "))
    ax.set_xlabel("Macro ROC–AUC")
    ax.set_title("Matched receptor-grouped predictive comparison")
    fig.tight_layout()
    for suffix in ("png", "svg"):
        fig.savefig(figure_dir / f"matched_predictive_comparison.{suffix}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    atomic_csv(figure_dir / "matched_predictive_comparison_source.tsv", performance)

    paired = results["paired"].copy()
    paired = paired.loc[paired.get("status", pd.Series(dtype=str)).eq("ok")]
    if not paired.empty:
        paired["label"] = paired["candidate"] + " − " + paired["baseline"]
        paired = paired.sort_values("delta_auc_candidate_minus_baseline")
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        y = np.arange(len(paired))
        values = paired["delta_auc_candidate_minus_baseline"].to_numpy(float)
        low = paired["delta_auc_ci_low"].to_numpy(float)
        high = paired["delta_auc_ci_high"].to_numpy(float)
        ax.errorbar(values, y, xerr=np.vstack([values - low, high - values]), fmt="o", capsize=4)
        ax.axvline(0.0, linestyle="--", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(paired["label"].str.replace("_", " "))
        ax.set_xlabel("Paired ΔAUC")
        ax.set_title("Incremental value of new receptor-only feature blocks")
        fig.tight_layout()
        for suffix in ("png", "svg"):
            fig.savefig(figure_dir / f"incremental_value_forest.{suffix}", dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        atomic_csv(figure_dir / "incremental_value_forest_source.tsv", paired)

    transfer = results["transfer_metrics"].copy()
    transfer = transfer.loc[transfer.get("target", pd.Series(dtype=str)).eq("macro_functional")]
    if not transfer.empty:
        pivot = transfer.pivot(index="experiment", columns="definition", values="roc_auc")
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        pivot.plot(kind="bar", ax=ax)
        ax.axhline(0.5, linestyle="--", linewidth=1)
        ax.set_ylabel("Macro functional ROC–AUC")
        ax.set_xlabel("")
        ax.set_title("Held-receptor-out transfer to receptor-only functional labels")
        ax.legend(title="Definition", frameon=False)
        fig.tight_layout()
        for suffix in ("png", "svg"):
            fig.savefig(figure_dir / f"strict_receptor_only_transfer.{suffix}", dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        atomic_csv(figure_dir / "strict_receptor_only_transfer_source.tsv", transfer)


def interpretation_gate(results: Mapping[str, Any]) -> dict[str, Any]:
    paired = results["paired"].copy()
    lookup = {
        (row.baseline, row.candidate): row
        for row in paired.itertuples(index=False)
        if getattr(row, "status", "") == "ok"
    }
    new_vs_missing = lookup.get(("missingness_only", "all_new_static"))
    new_vs_existing = lookup.get(("existing_geometry", "all_new_static"))
    new_vs_sequence = lookup.get(("sequence_only", "all_new_static"))
    seq_increment = lookup.get(("sequence_only", "sequence_plus_new_static"))

    def supported(row: Any | None) -> bool:
        return bool(row is not None and getattr(row, "auc_improvement_supported", False))

    transfer = results["transfer_metrics"]
    transfer_primary = transfer.loc[
        transfer.get("definition", pd.Series(dtype=str)).eq("ligand_permissive")
        & transfer.get("target", pd.Series(dtype=str)).eq("macro_functional")
        & transfer.get("experiment", pd.Series(dtype=str)).isin(["all_new_static", "sequence_plus_new_static"])
    ]
    transfer_above_chance = bool(
        not transfer_primary.empty
        and pd.to_numeric(transfer_primary["bootstrap_ci_low"], errors="coerce").gt(0.5).any()
    )

    if supported(new_vs_missing) and supported(new_vs_existing) and supported(seq_increment) and transfer_above_chance:
        outcome = "strong_support"
        wording = "Receptor-only physicochemical or mechanical features contain transferable coupling-preference information beyond sequence and static endpoint geometry."
    elif supported(new_vs_missing) or supported(new_vs_existing):
        outcome = "limited_support"
        wording = "New receptor-chain-only features show coupling-associated information, but a transferable receptor-intrinsic code remains unproven."
    else:
        outcome = "no_support"
        wording = "The tested receptor-only static physicochemical and mechanical descriptors do not reveal a transferable coupling-preference signal in the currently available sample."

    return {
        "outcome": outcome,
        "automated_cautious_wording": wording,
        "new_static_exceeds_missingness": supported(new_vs_missing),
        "new_static_exceeds_existing_geometry": supported(new_vs_existing),
        "new_static_exceeds_sequence": supported(new_vs_sequence),
        "sequence_plus_new_exceeds_sequence": supported(seq_increment),
        "strict_transfer_ci_above_chance": transfer_above_chance,
        "critical_limit": "The strict receptor-only functional evaluation contains few independent receptors and is exploratory.",
        "prohibited_inference": "Do not conclude that biological coupling preference is absent when this analysis is negative.",
    }


def write_summary(output_root: Path, audit: Mapping[str, Any], results: Mapping[str, Any], gate: Mapping[str, Any]) -> None:
    performance = results["performance"].sort_values("macro_roc_auc", ascending=False)
    lines = [
        "# Receptor-only static selectivity run summary",
        "",
        f"Automated outcome: **{gate['outcome']}**",
        "",
        gate["automated_cautious_wording"],
        "",
        "## Scope",
        "",
        "This first pass tests receptor-chain contact networks, intracellular surface chemistry, approximate Coulombic electrostatics and generic-position anisotropic-network susceptibility. It does not run trajectory, hydration-kinetic, Markov-state or free-energy analyses because validated receptor-only trajectories were not present.",
        "",
        "## Matched bound-structure performance",
        "",
        "| Model | Macro ROC–AUC | 95% receptor-bootstrap interval | Converged folds |",
        "|---|---:|---:|---:|",
    ]
    for row in performance.itertuples(index=False):
        lines.append(
            f"| {row.experiment} | {row.macro_roc_auc:.3f} | {row.macro_roc_auc_ci_low:.3f}–{row.macro_roc_auc_ci_high:.3f} | {row.valid_folds} |"
        )
    lines += [
        "",
        "## Interpretation limits",
        "",
        f"- {gate['critical_limit']}",
        f"- {gate['prohibited_inference']}",
        "- Approximate electrostatic descriptors are not Poisson–Boltzmann potentials.",
        "- Normal-mode descriptors represent mechanical susceptibility, not free energy.",
        "- Every new feature was extracted from the selected receptor chain only.",
    ]
    (output_root / "run_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(
    project_dir: str | Path,
    output_root: str | Path,
    *,
    config_path: str | Path | None = None,
    mode: str = "quick",
    workers: int = 8,
    resume: bool = True,
    seed: int = 20272729,
    stage: str = "all",
) -> dict[str, Any]:
    project = Path(project_dir).expanduser().resolve()
    output = Path(output_root).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = load_config(project, Path(config_path).resolve() if config_path else None)
    manifests = build_manifests(project, config, mode, seed)

    audit = run_input_audit(project, output, config, manifests)
    if stage == "audit":
        manifest = {"success": True, "stage": "audit", "output_root": str(output)}
        atomic_json(output / "analysis_manifest.json", manifest)
        return manifest

    features = extract_all_features(output, manifests, config, workers, resume, mode)
    if stage == "features":
        manifest = {"success": True, "stage": "features", "output_root": str(output), "feature_rows": len(features)}
        atomic_json(output / "analysis_manifest.json", manifest)
        return manifest

    results = run_models(project, output, config, manifests, features, mode, seed, resume)
    make_figures(output, results, dpi=int(config.get("figures", {}).get("dpi", 600)))
    gate = interpretation_gate(results)
    atomic_json(output / "interpretation_gate.json", gate)
    write_summary(output, audit, results, gate)

    audit_frame = results["audits"]
    leakage_ok = bool(
        not audit_frame.empty
        and pd.to_numeric(audit_frame.get("group_overlap", pd.Series(dtype=float)), errors="coerce").fillna(0).eq(0).all()
    )
    convergence_ok = bool(
        not results["performance"].empty
        and results["performance"]["all_folds_converged"].fillna(False).astype(bool).all()
    )
    manifest = {
        "analysis": "receptor_only_dynamic_selectivity_static_first_pass",
        "version": "0.1.1",
        "mode": mode,
        "stage": stage,
        "success": leakage_ok and convergence_ok,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "output_root": str(output),
        "bound_structures": int(len(manifests["bound"])),
        "bound_receptors": int(manifests["bound"]["receptor_name"].nunique()),
        "feature_rows": int(len(features)),
        "leakage_audit_passed": leakage_ok,
        "all_model_folds_converged": convergence_ok,
        "interpretation_outcome": gate["outcome"],
        "config_sha256": sha256_file(Path(config["_config_path"])),
    }
    atomic_json(output / "analysis_manifest.json", manifest)
    return manifest


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Run the static receptor-only dynamic-selectivity extension.")
    command.add_argument("--project-dir", required=True)
    command.add_argument("--output-root", required=True)
    command.add_argument("--config")
    command.add_argument("--mode", choices=("quick", "full"), default="quick")
    command.add_argument("--workers", type=int, default=8)
    command.add_argument("--seed", type=int, default=20272729)
    command.add_argument("--stage", choices=("audit", "features", "all"), default="all")
    command.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    return command


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    manifest = run_pipeline(
        args.project_dir,
        args.output_root,
        config_path=args.config,
        mode=args.mode,
        workers=args.workers,
        resume=args.resume,
        seed=args.seed,
        stage=args.stage,
    )
    print(json.dumps(manifest, indent=2), flush=True)
    return 0 if manifest.get("success", False) else 1


if __name__ == "__main__":
    raise SystemExit(main())
