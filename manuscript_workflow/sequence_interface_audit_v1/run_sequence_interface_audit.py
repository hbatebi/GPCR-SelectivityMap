#!/usr/bin/env python3
"""Sequence-interface, phylogeny, granularity and label-provenance audit.

This analysis is designed as the final explanatory control for the GPCR
endpoint-selectivity manuscript.  It tests whether the sequence advantage is
consistent with transferable coupling-interface chemistry, broad phylogeny,
or differences in the measurement granularity of sequence and structure.

The script reuses the v0.8.1 receptor-grouped modelling utilities where
possible.  It adds:

* full-sequence, physically defined interface-sequence, non-interface,
  taxonomy, interface-plus-taxonomy, geometry and sequence-plus-geometry
  models on identical receptor and sequence-cluster folds;
* a global-sequence-identity k-nearest-neighbour phylogeny baseline;
* matched non-interface position controls;
* within-receptor positional scrambling and across-receptor position-wise
  scrambling controls;
* receptor-by-transducer-family granularity audit and a single-transducer
  receptor sensitivity analysis;
* bound-label and receptor-only functional-label provenance audits;
* paired receptor-bootstrap comparisons, coefficient localisation and
  publication-quality SVG/PDF/PNG figures.

Important interpretation limit
------------------------------
The physically defined interface region is fixed before fitting.  The
coefficient localisation analysis is descriptive.  The script does not claim
causality and does not treat database-consensus functional labels as directly
measured unless that provenance is explicitly present in the input table.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import inspect
import json
import math
import multiprocessing as mp
import os
import re
import sys
import time
import traceback
import warnings
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import entropy as scipy_entropy
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    log_loss,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer

PRIMARY_CLASSES = ("Gi/o", "Gq/11", "Gs")

# This region is defined from the receptor side of the intracellular coupling
# interface, not from fitted coefficients.  The shorter top-weighted list is
# reported only as a descriptive, post-hoc localisation output.
INTERFACE_POSITIONS = (
    # Cytoplasmic TM3 and DRY-proximal positions
    "3x49", "3x50", "3x53", "3x55",
    # Cytoplasmic TM5 extension
    "5x58", "5x61", "5x64", "5x68", "5x72",
    # Cytoplasmic TM6 / ICL3-proximal positions
    "6x23", "6x30", "6x34", "6x37", "6x40",
    # Cytoplasmic TM7
    "7x53", "7x54", "7x55", "7x56",
    # Helix 8
    "8x47", "8x50", "8x53", "8x59",
    # ICL2 block
    "34x50", "34x51", "34x52", "34x53",
    "34x54", "34x55", "34x56", "34x57",
)

# Positions cited after inspection of the original full-sequence coefficients.
# These are explicitly descriptive and are not the primary feature definition.
POSTHOC_TOP_INTERFACE_POSITIONS = (
    "34x51", "34x57", "34x54", "34x52", "34x53", "34x56",
    "5x68", "5x64", "5x61", "6x34", "6x30", "8x53",
    "7x54", "3x53", "3x55",
)

MODEL_NAMES = (
    "full_sequence",
    "interface_sequence",
    "noninterface_sequence",
    "coarse_taxonomy",
    "fine_family",
    "interface_plus_taxonomy",
    "endpoint_geometry",
    "sequence_plus_geometry",
)

VALIDATION_MAP = {
    "receptor": ("receptor_grouped_5fold", "receptor_name"),
    "seq30": ("sequence_cluster_30_grouped_5fold", "sequence_cluster_30"),
    "seq40": ("sequence_cluster_40_grouped_5fold", "sequence_cluster_40"),
    "seq50": ("sequence_cluster_50_grouped_5fold", "sequence_cluster_50"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GPCR sequence-interface and phylogeny controls."
    )
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--mode", choices=("smoke", "quick", "full"), default="full"
    )
    parser.add_argument("--workers", type=int, default=48)
    parser.add_argument("--seed", type=int, default=20272729)
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--control-sets", type=int, default=50)
    parser.add_argument("--scrambles", type=int, default=50)
    parser.add_argument("--coefficient-permutations", type=int, default=10000)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--resume", action="store_true", default=True)
    return parser.parse_args()


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    tmp.replace(path)


def atomic_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False)
    tmp.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalise_position(column: str) -> str:
    match = re.search(r"gpcrdb_([^_]+)_expected_aa", str(column))
    return match.group(1) if match else str(column)


def position_column(position: str) -> str:
    return f"gpcrdb_{position}_expected_aa"


def sequence_columns(frame: pd.DataFrame) -> list[str]:
    return sorted(
        column
        for column in frame.columns
        if column.startswith("gpcrdb_") and column.endswith("_expected_aa")
    )


def coerce_sequence(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        out[column] = (
            out[column]
            .where(out[column].notna(), np.nan)
            .astype("object")
        )
        mask = out[column].notna()
        out.loc[mask, column] = (
            out.loc[mask, column]
            .astype(str)
            .str.strip()
            .str.upper()
            .replace({"": np.nan, "NAN": np.nan, "NONE": np.nan})
        )
    return out


def receptor_signatures(
    frame: pd.DataFrame, columns: Sequence[str]
) -> dict[str, tuple[str, ...]]:
    signatures: dict[str, tuple[str, ...]] = {}
    for receptor, part in frame.groupby("receptor_name", sort=True):
        values = []
        for column in columns:
            observed = part[column].dropna().astype(str)
            values.append(observed.mode().iloc[0] if len(observed) else "")
        signatures[str(receptor)] = tuple(values)
    return signatures


def pair_identity(a: Sequence[str], b: Sequence[str]) -> float:
    matches = []
    for x, y in zip(a, b):
        if x and y:
            matches.append(x == y)
    return float(np.mean(matches)) if matches else 0.0


def entropy_and_coverage(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    receptor_frame = frame.drop_duplicates("receptor_name")
    rows = []
    for column in columns:
        values = receptor_frame[column].dropna().astype(str)
        frequencies = values.value_counts(normalize=True)
        rows.append(
            {
                "feature": column,
                "generic_position": normalise_position(column),
                "coverage": float(receptor_frame[column].notna().mean()),
                "n_unique": int(values.nunique()),
                "entropy_nats": float(scipy_entropy(frequencies.to_numpy(float)))
                if len(frequencies)
                else np.nan,
                "is_interface": normalise_position(column) in INTERFACE_POSITIONS,
                "is_posthoc_top_interface": normalise_position(column)
                in POSTHOC_TOP_INTERFACE_POSITIONS,
            }
        )
    return pd.DataFrame(rows)


def make_settings(config: Any, *, stable_l2: bool = False) -> dict[str, Any]:
    settings = dict(config.analysis.get("coupling_specificity", {}))
    settings["inner_folds"] = int(settings.get("nested_inner_folds", 3))
    settings["outer_folds"] = int(settings.get("outer_folds", 5))
    settings["max_iter"] = max(50000, int(settings.get("max_iter", 50000)))
    settings["tolerance"] = float(settings.get("tolerance", 1e-4))
    if stable_l2:
        settings.update(
            {
                "solver": "lbfgs",
                "elastic_net_grid": {
                    "C": [0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
                    "l1_ratio": [0.0],
                },
            }
        )
    return settings


def feature_specs(frame: pd.DataFrame, config: Any) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    from gpcr_icl2.microswitches.publication_controls import (
        FeatureSpec,
        build_feature_catalog,
    )

    work, catalog, availability = build_feature_catalog(
        frame, config, minimum_feature_coverage=0.75
    )
    seq = sequence_columns(work)
    interface = [
        position_column(position)
        for position in INTERFACE_POSITIONS
        if position_column(position) in seq
    ]
    noninterface = [column for column in seq if column not in interface]

    coarse = [
        column
        for column in ("receptor_ligand_class", "receptor_subfamily")
        if column in work and work[column].dropna().astype(str).nunique() > 1
    ]
    fine = [
        column
        for column in ("receptor_family",)
        if column in work and work[column].dropna().astype(str).nunique() > 1
    ]

    specs = {
        "full_sequence": FeatureSpec(
            "full_sequence", (), tuple(seq), "All generic-position sequence features"
        ),
        "interface_sequence": FeatureSpec(
            "interface_sequence",
            (),
            tuple(interface),
            "Physically defined intracellular coupling-interface sequence",
        ),
        "noninterface_sequence": FeatureSpec(
            "noninterface_sequence",
            (),
            tuple(noninterface),
            "All available non-interface sequence positions",
        ),
        "coarse_taxonomy": FeatureSpec(
            "coarse_taxonomy", (), tuple(coarse), "Ligand class and receptor subfamily"
        ),
        "fine_family": FeatureSpec(
            "fine_family", (), tuple(fine), "Fine receptor-family identifier"
        ),
        "interface_plus_taxonomy": FeatureSpec(
            "interface_plus_taxonomy",
            (),
            tuple(dict.fromkeys(interface + coarse)),
            "Interface sequence plus coarse taxonomy",
        ),
        "endpoint_geometry": catalog["SCVcoreE_full_adjusted"],
        "sequence_plus_geometry": catalog["SCVcoreEQ_full_sequence_adjusted"],
    }
    return work, specs, availability


def run_existing_model_task(task: Mapping[str, Any]) -> dict[str, Any]:
    from gpcr_icl2.microswitches.config import ProjectConfig
    from gpcr_icl2.microswitches.publication_controls import (
        ExperimentSpec,
        FeatureSpec,
        cross_validate_experiment,
    )

    started = time.time()
    checkpoint = Path(task["checkpoint"])
    try:
        frame = pd.read_pickle(task["frame_pickle"])
        config = ProjectConfig.from_dir(task["config_dir"])
        spec_data = task["feature_spec"]
        spec = FeatureSpec(
            name=spec_data["name"],
            numeric=tuple(spec_data.get("numeric", [])),
            categorical=tuple(spec_data.get("categorical", [])),
            description=spec_data.get("description", ""),
        )
        experiment = ExperimentSpec(
            name=task["model"],
            feature_spec=spec,
            validation=task["validation"],
            group_column=task["group_column"],
            frame_subset=task.get("frame_subset", "all"),
            nested_tuning=True,
        )
        settings = make_settings(config, stable_l2=True)
        result = cross_validate_experiment(
            frame,
            experiment,
            settings,
            int(task["seed"]),
            minimum_fold_coverage=0.0,
            bootstrap_iterations=0,
        )
        performance = dict(result["performance"])
        predictions = result["predictions"]
        audit = result["audit"]
        coefficients = result.get("coefficients", pd.DataFrame())

        payload = {
            "status": "ok",
            "task_id": task["task_id"],
            "model": task["model"],
            "validation": task["validation"],
            "group_column": task["group_column"],
            "frame_subset": task.get("frame_subset", "all"),
            "runtime_seconds": time.time() - started,
            "model_specification": "nested_L2_lbfgs",
            "performance": performance,
            "predictions": predictions.to_dict(orient="records"),
            "audit": audit.to_dict(orient="records"),
            "coefficients": coefficients.to_dict(orient="records")
            if isinstance(coefficients, pd.DataFrame)
            else [],
        }
        atomic_json(checkpoint, payload)
        return payload
    except Exception as error:
        payload = {
            "status": "failed",
            "task_id": task["task_id"],
            "model": task.get("model"),
            "validation": task.get("validation"),
            "runtime_seconds": time.time() - started,
            "error": repr(error),
            "traceback": traceback.format_exc(),
        }
        atomic_json(checkpoint, payload)
        return payload


def model_metrics(
    y_true: np.ndarray, probabilities: np.ndarray, classes: Sequence[str]
) -> dict[str, float]:
    valid = np.isfinite(probabilities).all(axis=1)
    if not valid.any():
        return {
            "macro_roc_auc": np.nan,
            "macro_pr_auc": np.nan,
            "balanced_accuracy": np.nan,
            "log_loss": np.nan,
            "matthews_correlation": np.nan,
        }
    y = y_true[valid]
    p = probabilities[valid]
    predicted = np.asarray(classes)[np.argmax(p, axis=1)]
    try:
        auc = float(
            roc_auc_score(
                y, p, labels=list(classes), multi_class="ovr", average="macro"
            )
        )
    except ValueError:
        auc = np.nan
    pr = []
    for index, label in enumerate(classes):
        target = (y == label).astype(int)
        if target.min() == target.max():
            continue
        pr.append(float(average_precision_score(target, p[:, index])))
    return {
        "macro_roc_auc": auc,
        "macro_pr_auc": float(np.mean(pr)) if pr else np.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
        "log_loss": float(log_loss(y, p, labels=list(classes))),
        "matthews_correlation": float(matthews_corrcoef(y, predicted)),
    }


def dynamic_splits(
    y: np.ndarray, groups: np.ndarray, requested: int, seed: int
) -> list[tuple[np.ndarray, np.ndarray]]:
    unique_groups = np.unique(groups)
    class_group_counts = [len(np.unique(groups[y == label])) for label in np.unique(y)]
    n_splits = min(requested, len(unique_groups), *(class_group_counts or [len(unique_groups)]))
    if n_splits < 2:
        return []
    splitter = StratifiedGroupKFold(
        n_splits=n_splits, shuffle=True, random_state=seed
    )
    return list(splitter.split(np.zeros(len(y)), y, groups))


def receptor_label_profiles(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    classes = list(PRIMARY_CLASSES)
    profiles: dict[str, np.ndarray] = {}
    for receptor, part in frame.groupby("receptor_name"):
        labels = sorted(set(part["transducer_family"].astype(str)))
        vector = np.zeros(len(classes), dtype=float)
        for label in labels:
            if label in classes:
                vector[classes.index(label)] = 1.0
        if vector.sum() > 0:
            vector /= vector.sum()
        profiles[str(receptor)] = vector
    return profiles


def tune_knn_k(
    train_frame: pd.DataFrame,
    signatures: Mapping[str, Sequence[str]],
    candidate_k: Sequence[int],
    seed: int,
) -> int:
    receptors = train_frame["receptor_name"].astype(str).to_numpy()
    y = train_frame["transducer_family"].astype(str).to_numpy()
    splits = dynamic_splits(y, receptors, 3, seed)
    if len(splits) < 2:
        return 5
    scores = defaultdict(list)
    for inner_train, inner_test in splits:
        fit = train_frame.iloc[inner_train]
        test = train_frame.iloc[inner_test]
        profiles = receptor_label_profiles(fit)
        train_receptors = sorted(profiles)
        for k in candidate_k:
            probabilities = []
            truth = []
            for row in test.itertuples(index=False):
                receptor = str(row.receptor_name)
                neighbours = []
                for candidate in train_receptors:
                    identity = pair_identity(signatures[receptor], signatures[candidate])
                    neighbours.append((identity, candidate))
                neighbours.sort(reverse=True)
                selected = neighbours[: min(k, len(neighbours))]
                if not selected:
                    probability = np.ones(len(PRIMARY_CLASSES)) / len(PRIMARY_CLASSES)
                else:
                    weights = np.asarray([max(value, 1e-6) ** 2 for value, _ in selected])
                    vectors = np.asarray([profiles[name] for _, name in selected])
                    probability = np.average(vectors, axis=0, weights=weights)
                    if probability.sum() <= 0:
                        probability = np.ones(len(PRIMARY_CLASSES)) / len(PRIMARY_CLASSES)
                    else:
                        probability = (probability + 0.01) / (probability.sum() + 0.01 * len(PRIMARY_CLASSES))
                probabilities.append(probability)
                truth.append(str(row.transducer_family))
            score = model_metrics(
                np.asarray(truth), np.asarray(probabilities), PRIMARY_CLASSES
            )["macro_roc_auc"]
            if np.isfinite(score):
                scores[k].append(float(score))
    if not scores:
        return 5
    return max(scores, key=lambda value: np.mean(scores[value]))


def run_knn_task(task: Mapping[str, Any]) -> dict[str, Any]:
    started = time.time()
    checkpoint = Path(task["checkpoint"])
    try:
        frame = pd.read_pickle(task["frame_pickle"]).reset_index(drop=True)
        seq_cols = sequence_columns(frame)
        signatures = receptor_signatures(frame, seq_cols)
        classes = tuple(PRIMARY_CLASSES)
        y = frame["transducer_family"].astype(str).to_numpy()
        groups = frame[task["group_column"]].fillna("unknown").astype(str).to_numpy()
        splits = dynamic_splits(y, groups, 5, int(task["seed"]))
        probabilities = np.full((len(frame), len(classes)), np.nan)
        fold_ids = np.full(len(frame), -1, dtype=int)
        audit_rows = []
        neighbour_rows = []

        for fold, (train_idx, test_idx) in enumerate(splits):
            train = frame.iloc[train_idx].copy()
            test = frame.iloc[test_idx].copy()
            if set(train["transducer_family"].astype(str)) != set(classes):
                audit_rows.append(
                    {
                        "fold": fold,
                        "status": "skipped_training_missing_class",
                        "train_n": len(train),
                        "test_n": len(test),
                    }
                )
                continue
            chosen_k = tune_knn_k(
                train,
                signatures,
                candidate_k=(1, 3, 5, 9, 15),
                seed=int(task["seed"]) + fold + 100,
            )
            profiles = receptor_label_profiles(train)
            train_receptors = sorted(profiles)
            for index in test_idx:
                receptor = str(frame.loc[index, "receptor_name"])
                neighbours = []
                for candidate in train_receptors:
                    identity = pair_identity(signatures[receptor], signatures[candidate])
                    neighbours.append((identity, candidate))
                neighbours.sort(reverse=True)
                selected = neighbours[: min(chosen_k, len(neighbours))]
                weights = np.asarray([max(value, 1e-6) ** 2 for value, _ in selected])
                vectors = np.asarray([profiles[name] for _, name in selected])
                probability = (
                    np.average(vectors, axis=0, weights=weights)
                    if len(selected)
                    else np.ones(len(classes)) / len(classes)
                )
                probability = (probability + 0.01) / (probability.sum() + 0.01 * len(classes))
                probabilities[index] = probability
                fold_ids[index] = fold
                neighbour_rows.append(
                    {
                        "analysis_row_id": frame.loc[index, "analysis_row_id"],
                        "fold": fold,
                        "receptor_name": receptor,
                        "chosen_k": chosen_k,
                        "nearest_training_identity": selected[0][0] if selected else np.nan,
                        "nearest_training_receptor": selected[0][1] if selected else "",
                    }
                )
            audit_rows.append(
                {
                    "fold": fold,
                    "status": "ok",
                    "train_n": len(train_idx),
                    "test_n": len(test_idx),
                    "train_groups": len(np.unique(groups[train_idx])),
                    "test_groups": len(np.unique(groups[test_idx])),
                    "group_overlap": len(
                        set(groups[train_idx]).intersection(groups[test_idx])
                    ),
                    "chosen_k": chosen_k,
                }
            )

        valid = np.isfinite(probabilities).all(axis=1)
        predicted = np.full(len(frame), "", dtype=object)
        predicted[valid] = np.asarray(classes)[np.argmax(probabilities[valid], axis=1)]
        predictions = frame[
            [
                "analysis_row_id",
                "pdb_id",
                "chain_id",
                "receptor_name",
                "transducer_family",
            ]
        ].copy()
        predictions["experiment"] = "global_identity_knn"
        predictions["validation"] = task["validation"]
        predictions["outer_fold"] = fold_ids
        for index, label in enumerate(classes):
            predictions[f"probability_{label}"] = probabilities[:, index]
        predictions["predicted_class"] = predicted
        metrics = model_metrics(y, probabilities, classes)
        performance = {
            "experiment": "global_identity_knn",
            "validation": task["validation"],
            "n": len(frame),
            "n_receptors": int(frame["receptor_name"].nunique()),
            "n_predictions": int(valid.sum()),
            "valid_folds": int(
                sum(row.get("status") == "ok" for row in audit_rows)
            ),
            "all_folds_converged": True,
            **metrics,
        }
        payload = {
            "status": "ok",
            "task_id": task["task_id"],
            "model": "global_identity_knn",
            "validation": task["validation"],
            "runtime_seconds": time.time() - started,
            "performance": performance,
            "predictions": predictions.to_dict(orient="records"),
            "audit": audit_rows,
            "neighbours": neighbour_rows,
        }
        atomic_json(checkpoint, payload)
        return payload
    except Exception as error:
        payload = {
            "status": "failed",
            "task_id": task["task_id"],
            "model": "global_identity_knn",
            "validation": task.get("validation"),
            "runtime_seconds": time.time() - started,
            "error": repr(error),
            "traceback": traceback.format_exc(),
        }
        atomic_json(checkpoint, payload)
        return payload


def fixed_l2_predictions(
    frame: pd.DataFrame,
    columns: Sequence[str],
    group_column: str,
    seed: int,
    transform: str = "none",
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Fast fixed-C L2 control using fold-local sequence scrambling."""
    work = frame.reset_index(drop=True).copy()
    classes = tuple(PRIMARY_CLASSES)
    y = work["transducer_family"].astype(str).to_numpy()
    groups = work[group_column].fillna("unknown").astype(str).to_numpy()
    splits = dynamic_splits(y, groups, 5, seed)
    probabilities = np.full((len(work), len(classes)), np.nan)
    folds = np.full(len(work), -1, dtype=int)
    audits = []

    for fold, (train_idx, test_idx) in enumerate(splits):
        train = work.iloc[train_idx].copy()
        test = work.iloc[test_idx].copy()
        rng = np.random.default_rng(seed + 1000 + fold)

        if transform == "within_receptor_position_shuffle":
            # Sequence is constant within receptor.  Generate one shuffled
            # signature per receptor, then apply it to all rows of that receptor.
            for part in (train, test):
                for receptor, indices in part.groupby("receptor_name").groups.items():
                    row = part.loc[next(iter(indices)), list(columns)].to_numpy(object)
                    shuffled = row.copy()
                    rng.shuffle(shuffled)
                    part.loc[list(indices), list(columns)] = shuffled
        elif transform == "position_shuffle_across_receptors":
            # Shuffle receptor identities separately at each position within
            # train and test.  This preserves position-wise amino-acid frequency
            # but destroys receptor-specific combinations.  Test labels are not
            # used and train/test data are never mixed.
            for part in (train, test):
                unique = part.drop_duplicates("receptor_name").set_index("receptor_name")
                shuffled_map: dict[str, dict[str, Any]] = {
                    receptor: {} for receptor in unique.index.astype(str)
                }
                receptors = unique.index.astype(str).to_numpy()
                for column in columns:
                    values = unique[column].to_numpy(object).copy()
                    rng.shuffle(values)
                    for receptor, value in zip(receptors, values):
                        shuffled_map[receptor][column] = value
                for receptor, indices in part.groupby("receptor_name").groups.items():
                    for column in columns:
                        part.loc[list(indices), column] = shuffled_map[str(receptor)][column]
        elif transform != "none":
            raise ValueError(f"Unknown transform: {transform}")

        preprocess = ColumnTransformer(
            [
                (
                    "categorical",
                    Pipeline(
                        [
                            (
                                "imputer",
                                SimpleImputer(
                                    strategy="most_frequent", keep_empty_features=True
                                ),
                            ),
                            (
                                "onehot",
                                OneHotEncoder(
                                    handle_unknown="ignore", sparse_output=False
                                ),
                            ),
                        ]
                    ),
                    list(columns),
                )
            ],
            remainder="drop",
        )
        model = LogisticRegression(
            solver="lbfgs",
            penalty="l2",
            C=1.0,
            class_weight="balanced",
            max_iter=50000,
            tol=1e-4,
            random_state=seed + fold,
        )
        pipe = Pipeline([("preprocess", preprocess), ("model", model)])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            pipe.fit(train[list(columns)], train["transducer_family"].astype(str))
        local = pipe.predict_proba(test[list(columns)])
        for index, label in enumerate(pipe.named_steps["model"].classes_):
            probabilities[test_idx, classes.index(str(label))] = local[:, index]
        folds[test_idx] = fold
        n_iter = int(np.max(pipe.named_steps["model"].n_iter_))
        audits.append(
            {
                "fold": fold,
                "status": "ok",
                "train_n": len(train_idx),
                "test_n": len(test_idx),
                "group_overlap": len(
                    set(groups[train_idx]).intersection(groups[test_idx])
                ),
                "n_iter": n_iter,
                "converged": not any(
                    issubclass(warning.category, ConvergenceWarning)
                    for warning in caught
                )
                and n_iter < 50000,
            }
        )

    valid = np.isfinite(probabilities).all(axis=1)
    predicted = np.full(len(work), "", dtype=object)
    predicted[valid] = np.asarray(classes)[np.argmax(probabilities[valid], axis=1)]
    predictions = work[
        [
            "analysis_row_id",
            "pdb_id",
            "chain_id",
            "receptor_name",
            "transducer_family",
        ]
    ].copy()
    predictions["outer_fold"] = folds
    for index, label in enumerate(classes):
        predictions[f"probability_{label}"] = probabilities[:, index]
    predictions["predicted_class"] = predicted
    return predictions, model_metrics(y, probabilities, classes), pd.DataFrame(audits)


def run_fast_control_task(task: Mapping[str, Any]) -> dict[str, Any]:
    started = time.time()
    checkpoint = Path(task["checkpoint"])
    try:
        frame = pd.read_pickle(task["frame_pickle"])
        predictions, metrics, audit = fixed_l2_predictions(
            frame,
            task["columns"],
            task["group_column"],
            int(task["seed"]),
            task.get("transform", "none"),
        )
        predictions["experiment"] = task["model"]
        predictions["validation"] = task["validation"]
        performance = {
            "experiment": task["model"],
            "validation": task["validation"],
            "control_iteration": task.get("iteration", -1),
            "n": len(frame),
            "n_receptors": int(frame["receptor_name"].nunique()),
            "n_predictions": int(predictions["predicted_class"].ne("").sum()),
            "valid_folds": int(len(audit)),
            "all_folds_converged": bool(audit["converged"].all())
            if len(audit)
            else False,
            **metrics,
        }
        payload = {
            "status": "ok",
            "task_id": task["task_id"],
            "model": task["model"],
            "validation": task["validation"],
            "iteration": task.get("iteration", -1),
            "columns": list(task["columns"]),
            "runtime_seconds": time.time() - started,
            "performance": performance,
            "predictions": predictions.to_dict(orient="records"),
            "audit": audit.to_dict(orient="records"),
        }
        atomic_json(checkpoint, payload)
        return payload
    except Exception as error:
        payload = {
            "status": "failed",
            "task_id": task["task_id"],
            "model": task.get("model"),
            "validation": task.get("validation"),
            "iteration": task.get("iteration", -1),
            "runtime_seconds": time.time() - started,
            "error": repr(error),
            "traceback": traceback.format_exc(),
        }
        atomic_json(checkpoint, payload)
        return payload


def matched_noninterface_sets(
    frame: pd.DataFrame,
    interface_columns: Sequence[str],
    noninterface_columns: Sequence[str],
    n_sets: int,
    seed: int,
) -> tuple[list[list[str]], pd.DataFrame]:
    """Generate unsupervised entropy/coverage matched non-interface controls."""
    stats = entropy_and_coverage(frame, list(interface_columns) + list(noninterface_columns))
    stats = stats.set_index("feature")
    rng = np.random.default_rng(seed)
    sets = []
    rows = []
    target_size = min(len(interface_columns), len(noninterface_columns))

    for iteration in range(n_sets):
        available = set(noninterface_columns)
        selected = []
        order = list(interface_columns)
        rng.shuffle(order)
        for interface_column in order[:target_size]:
            if not available:
                break
            target = stats.loc[interface_column]
            candidates = []
            for candidate in available:
                record = stats.loc[candidate]
                distance = math.sqrt(
                    (float(target["coverage"]) - float(record["coverage"])) ** 2
                    + (
                        float(target["entropy_nats"])
                        - float(record["entropy_nats"])
                    )
                    ** 2
                )
                candidates.append((distance, float(rng.random()), candidate))
            candidates.sort()
            chosen = candidates[0][2]
            available.remove(chosen)
            selected.append(chosen)
            rows.append(
                {
                    "control_iteration": iteration,
                    "interface_feature": interface_column,
                    "interface_position": normalise_position(interface_column),
                    "matched_noninterface_feature": chosen,
                    "matched_noninterface_position": normalise_position(chosen),
                    "distance": candidates[0][0],
                }
            )
        sets.append(sorted(selected))
    return sets, pd.DataFrame(rows)


def parse_checkpoint(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if payload.get("status") == "ok" else None


def prediction_frame(payload: Mapping[str, Any]) -> pd.DataFrame:
    frame = pd.DataFrame(payload.get("predictions", []))
    if frame.empty:
        return frame
    if "analysis_row_id" not in frame:
        # Existing cross_validate_experiment does not preserve this custom
        # column.  Reconstruct it from stable identifiers.
        frame["analysis_row_id"] = (
            frame["pdb_id"].astype(str)
            + "|"
            + frame.get("chain_id", pd.Series("", index=frame.index)).astype(str)
            + "|"
            + frame["receptor_name"].astype(str)
            + "|"
            + frame["transducer_family"].astype(str)
        )
    return frame


def receptor_bootstrap_metrics(
    predictions: pd.DataFrame,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    probability_columns = [f"probability_{label}" for label in PRIMARY_CLASSES]
    valid = predictions.dropna(subset=probability_columns).copy().reset_index(drop=True)
    if valid.empty:
        return {
            "bootstrap_valid_iterations": 0,
            "macro_roc_auc_ci_low": np.nan,
            "macro_roc_auc_ci_high": np.nan,
            "log_loss_ci_low": np.nan,
            "log_loss_ci_high": np.nan,
        }
    receptors = sorted(valid["receptor_name"].astype(str).unique())
    by_receptor = {
        receptor: valid.index[valid["receptor_name"].astype(str).eq(receptor)].to_numpy()
        for receptor in receptors
    }
    y = valid["transducer_family"].astype(str).to_numpy()
    probabilities = valid[probability_columns].to_numpy(float)
    rng = np.random.default_rng(seed)
    auc_values = []
    loss_values = []
    for _ in range(iterations):
        sampled = rng.choice(receptors, size=len(receptors), replace=True)
        indices = np.concatenate([by_receptor[receptor] for receptor in sampled])
        metrics = model_metrics(y[indices], probabilities[indices], PRIMARY_CLASSES)
        if np.isfinite(metrics["macro_roc_auc"]):
            auc_values.append(metrics["macro_roc_auc"])
        if np.isfinite(metrics["log_loss"]):
            loss_values.append(metrics["log_loss"])
    return {
        "bootstrap_valid_iterations": int(min(len(auc_values), len(loss_values))),
        "macro_roc_auc_ci_low": float(np.quantile(auc_values, 0.025)) if auc_values else np.nan,
        "macro_roc_auc_ci_high": float(np.quantile(auc_values, 0.975)) if auc_values else np.nan,
        "log_loss_ci_low": float(np.quantile(loss_values, 0.025)) if loss_values else np.nan,
        "log_loss_ci_high": float(np.quantile(loss_values, 0.975)) if loss_values else np.nan,
    }


def receptor_bootstrap_metric_difference(
    candidate: pd.DataFrame,
    reference: pd.DataFrame,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    keys = ["analysis_row_id"]
    probability_columns = [f"probability_{label}" for label in PRIMARY_CLASSES]
    left = candidate[
        keys
        + ["receptor_name", "transducer_family"]
        + probability_columns
    ].rename(columns={column: f"candidate__{column}" for column in probability_columns})
    right = reference[keys + probability_columns].rename(
        columns={column: f"reference__{column}" for column in probability_columns}
    )
    merged = left.merge(right, on=keys, how="inner")
    merged = merged.dropna()
    if merged.empty:
        return {
            "n_rows": 0,
            "n_receptors": 0,
            "auc_difference": np.nan,
            "auc_difference_ci_low": np.nan,
            "auc_difference_ci_high": np.nan,
            "log_loss_difference": np.nan,
            "log_loss_difference_ci_low": np.nan,
            "log_loss_difference_ci_high": np.nan,
        }
    receptors = sorted(merged["receptor_name"].astype(str).unique())
    by_receptor = {
        receptor: merged.index[merged["receptor_name"].astype(str).eq(receptor)].to_numpy()
        for receptor in receptors
    }
    y = merged["transducer_family"].astype(str).to_numpy()
    candidate_prob = merged[
        [f"candidate__probability_{label}" for label in PRIMARY_CLASSES]
    ].to_numpy(float)
    reference_prob = merged[
        [f"reference__probability_{label}" for label in PRIMARY_CLASSES]
    ].to_numpy(float)
    observed_candidate = model_metrics(y, candidate_prob, PRIMARY_CLASSES)
    observed_reference = model_metrics(y, reference_prob, PRIMARY_CLASSES)
    rng = np.random.default_rng(seed)
    auc_differences = []
    loss_differences = []
    for _ in range(iterations):
        sampled = rng.choice(receptors, size=len(receptors), replace=True)
        indices = np.concatenate([by_receptor[receptor] for receptor in sampled])
        candidate_metrics = model_metrics(
            y[indices], candidate_prob[indices], PRIMARY_CLASSES
        )
        reference_metrics = model_metrics(
            y[indices], reference_prob[indices], PRIMARY_CLASSES
        )
        if np.isfinite(candidate_metrics["macro_roc_auc"]) and np.isfinite(
            reference_metrics["macro_roc_auc"]
        ):
            auc_differences.append(
                candidate_metrics["macro_roc_auc"]
                - reference_metrics["macro_roc_auc"]
            )
        if np.isfinite(candidate_metrics["log_loss"]) and np.isfinite(
            reference_metrics["log_loss"]
        ):
            loss_differences.append(
                candidate_metrics["log_loss"] - reference_metrics["log_loss"]
            )
    return {
        "n_rows": int(len(merged)),
        "n_receptors": int(len(receptors)),
        "auc_difference": observed_candidate["macro_roc_auc"]
        - observed_reference["macro_roc_auc"],
        "auc_difference_ci_low": float(np.quantile(auc_differences, 0.025))
        if auc_differences
        else np.nan,
        "auc_difference_ci_high": float(np.quantile(auc_differences, 0.975))
        if auc_differences
        else np.nan,
        "log_loss_difference": observed_candidate["log_loss"]
        - observed_reference["log_loss"],
        "log_loss_difference_ci_low": float(
            np.quantile(loss_differences, 0.025)
        )
        if loss_differences
        else np.nan,
        "log_loss_difference_ci_high": float(
            np.quantile(loss_differences, 0.975)
        )
        if loss_differences
        else np.nan,
    }


def coefficient_position_summary(coefficients: pd.DataFrame) -> pd.DataFrame:
    if coefficients.empty:
        return pd.DataFrame()
    work = coefficients.copy()
    work["generic_position"] = work["transformed_feature"].map(normalise_position)
    work = work[work["generic_position"].str.contains("x", regex=False)]
    work["is_interface"] = work["generic_position"].isin(INTERFACE_POSITIONS)
    work["is_posthoc_top_interface"] = work["generic_position"].isin(
        POSTHOC_TOP_INTERFACE_POSITIONS
    )
    summary = (
        work.groupby(
            ["experiment", "validation", "generic_position", "is_interface", "is_posthoc_top_interface"],
            as_index=False,
        )
        .agg(
            mean_absolute_coefficient=("absolute_coefficient", "mean"),
            median_absolute_coefficient=("absolute_coefficient", "median"),
            maximum_absolute_coefficient=("absolute_coefficient", "max"),
            coefficient_entries=("absolute_coefficient", "size"),
            nonzero_fraction=(
                "absolute_coefficient",
                lambda values: float(np.mean(np.asarray(values) > 1e-8)),
            ),
        )
    )
    summary["rank_within_model_validation"] = summary.groupby(
        ["experiment", "validation"]
    )["mean_absolute_coefficient"].rank(method="min", ascending=False)
    return summary.sort_values(
        ["experiment", "validation", "rank_within_model_validation"]
    )


def interface_coefficient_enrichment(
    summary: pd.DataFrame, iterations: int, seed: int
) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rows = []
    for (experiment, validation), part in summary.groupby(
        ["experiment", "validation"]
    ):
        values = part.set_index("generic_position")["mean_absolute_coefficient"]
        interface = [position for position in INTERFACE_POSITIONS if position in values]
        if not interface or len(interface) >= len(values):
            continue
        observed = float(values.loc[interface].mean() - values.drop(interface).mean())
        all_positions = values.index.to_numpy()
        null = []
        for _ in range(iterations):
            selected = rng.choice(all_positions, size=len(interface), replace=False)
            other = np.setdiff1d(all_positions, selected)
            null.append(float(values.loc[selected].mean() - values.loc[other].mean()))
        p_value = (1 + int(np.sum(np.asarray(null) >= observed))) / (iterations + 1)
        rows.append(
            {
                "experiment": experiment,
                "validation": validation,
                "n_positions": int(len(values)),
                "n_interface_positions": int(len(interface)),
                "observed_interface_minus_noninterface_mean_abs_coefficient": observed,
                "permutation_p_value_one_sided": p_value,
                "null_mean": float(np.mean(null)),
                "null_ci_low": float(np.quantile(null, 0.025)),
                "null_ci_high": float(np.quantile(null, 0.975)),
                "interpretation": "descriptive_posthoc_coefficient_localisation",
            }
        )
    return pd.DataFrame(rows)


def build_granularity_audit(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    unit = (
        frame.groupby(["receptor_name", "transducer_family"], as_index=False)
        .agg(
            n_structure_rows=("pdb_id", "size"),
            pdb_ids=("pdb_id", lambda values: ";".join(sorted(values.astype(str)))),
        )
    )
    receptor = (
        frame.groupby("receptor_name", as_index=False)
        .agg(
            n_structure_rows=("pdb_id", "size"),
            n_transducer_families=("transducer_family", "nunique"),
            transducer_families=(
                "transducer_family",
                lambda values: ";".join(sorted(set(values.astype(str)))),
            ),
        )
    )
    receptor["single_transducer_receptor"] = receptor["n_transducer_families"].eq(1)
    return unit, receptor


def build_label_provenance(project: Path, frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    bound = frame[
        [
            column
            for column in (
                "pdb_id",
                "chain_id",
                "receptor_name",
                "transducer_family",
                "publication_reference",
            )
            if column in frame
        ]
    ].copy()
    bound["label_target"] = "bound_transducer_family"
    bound["label_provenance_category"] = "direct_structural_complex_assignment"
    bound["homology_propagation_risk"] = False
    bound["provenance_note"] = (
        "The class label is the transducer family physically present in the deposited complex."
    )

    candidates = [
        project / "data/manual/receptor_only_receptor_level_labels.csv",
        project / "data/manual/functional_coupling_annotations_audit.csv",
        project / "functional_coupling_annotations_audit.csv",
    ]
    functional_path = next((path for path in candidates if path.exists()), None)
    if functional_path is None:
        functional = pd.DataFrame(
            columns=[
                "receptor_name",
                "Gs",
                "Gi/o",
                "Gq/11",
                "annotation_source",
                "notes",
                "label_provenance_category",
                "direct_vs_inferred_resolved",
            ]
        )
        note = "No receptor-only functional annotation table was found."
    else:
        functional = pd.read_csv(functional_path, low_memory=False)
        functional["label_target"] = "receptor_only_functional_coupling"
        functional["label_provenance_category"] = (
            "database_consensus_plus_manual_curation"
        )
        functional["direct_vs_inferred_resolved"] = False
        functional["provenance_note"] = (
            "The current table does not provide source-level experimental versus homology-inferred provenance for each family label."
        )
        note = (
            f"Functional labels were loaded from {functional_path}.  The table reports consensus/manual curation but does not resolve every label into direct measurement versus homology inference."
        )
    return bound, functional, note


def save_figures(
    performance: pd.DataFrame,
    paired: pd.DataFrame,
    controls: pd.DataFrame,
    positions: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    import matplotlib.pyplot as plt

    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    if not performance.empty:
        display_models = [
            "full_sequence",
            "interface_sequence",
            "noninterface_sequence",
            "coarse_taxonomy",
            "global_identity_knn",
            "endpoint_geometry",
        ]
        part = performance[
            performance["experiment"].isin(display_models)
        ].copy()
        order_validations = [
            "receptor_grouped_5fold",
            "sequence_cluster_30_grouped_5fold",
            "sequence_cluster_40_grouped_5fold",
            "sequence_cluster_50_grouped_5fold",
        ]
        fig, ax = plt.subplots(figsize=(9.0, 5.2))
        x = np.arange(len(order_validations))
        width = 0.12
        for index, model in enumerate(display_models):
            values = []
            for validation in order_validations:
                row = part[
                    part["experiment"].eq(model)
                    & part["validation"].eq(validation)
                ]
                values.append(
                    float(row["macro_roc_auc"].iloc[0]) if len(row) else np.nan
                )
            ax.plot(x, values, marker="o", linewidth=1.8, label=model.replace("_", " "))
        ax.axhline(0.5, linestyle="--", linewidth=1)
        ax.set_xticks(x, ["Receptor", "30% cluster", "40% cluster", "50% cluster"])
        ax.set_ylabel("Macro ROC–AUC")
        ax.set_ylim(0.35, 0.9)
        ax.legend(frameon=False, ncol=2, fontsize=8)
        fig.tight_layout()
        for suffix in ("svg", "pdf", "png"):
            fig.savefig(
                figures / f"sequence_interface_transfer.{suffix}",
                dpi=dpi if suffix == "png" else None,
                bbox_inches="tight",
            )
        plt.close(fig)

    if not paired.empty:
        key = paired[
            paired["comparison"].isin(
                [
                    "interface_sequence_minus_coarse_taxonomy",
                    "interface_sequence_minus_global_identity_knn",
                    "interface_sequence_minus_noninterface_sequence",
                    "full_sequence_minus_endpoint_geometry",
                ]
            )
        ].copy()
        if not key.empty:
            key = key.sort_values(["validation", "comparison"])
            fig, ax = plt.subplots(figsize=(8.5, max(4.0, 0.38 * len(key) + 1.5)))
            y = np.arange(len(key))
            value = key["auc_difference"].to_numpy(float)
            low = key["auc_difference_ci_low"].to_numpy(float)
            high = key["auc_difference_ci_high"].to_numpy(float)
            ax.errorbar(
                value,
                y,
                xerr=np.vstack([value - low, high - value]),
                fmt="o",
                capsize=3,
            )
            labels = [
                f"{row.comparison.replace('_', ' ')} | {row.validation.replace('_grouped_5fold', '')}"
                for row in key.itertuples()
            ]
            ax.set_yticks(y, labels)
            ax.axvline(0, linestyle="--", linewidth=1)
            ax.set_xlabel("Paired ΔAUC")
            fig.tight_layout()
            for suffix in ("svg", "pdf", "png"):
                fig.savefig(
                    figures / f"paired_incremental_value.{suffix}",
                    dpi=dpi if suffix == "png" else None,
                    bbox_inches="tight",
                )
            plt.close(fig)

    if not controls.empty:
        part = controls[
            controls["validation"].eq("sequence_cluster_30_grouped_5fold")
        ].copy()
        if not part.empty:
            fig, ax = plt.subplots(figsize=(8.0, 4.8))
            groups = [
                "matched_noninterface_control",
                "within_receptor_position_shuffle",
                "position_shuffle_across_receptors",
            ]
            data = [
                part.loc[part["experiment"].eq(group), "macro_roc_auc"].dropna().to_numpy(float)
                for group in groups
            ]
            ax.boxplot(data, tick_labels=[g.replace("_", " ") for g in groups], showfliers=False)
            interface_row = controls[
                controls["experiment"].eq("interface_fixed_l2")
                & controls["validation"].eq("sequence_cluster_30_grouped_5fold")
            ]
            if len(interface_row):
                ax.axhline(
                    float(interface_row["macro_roc_auc"].iloc[0]),
                    linewidth=1.5,
                    label="Observed interface sequence",
                )
            ax.axhline(0.5, linestyle="--", linewidth=1)
            ax.set_ylabel("Macro ROC–AUC")
            ax.tick_params(axis="x", rotation=15)
            ax.legend(frameon=False)
            fig.tight_layout()
            for suffix in ("svg", "pdf", "png"):
                fig.savefig(
                    figures / f"interface_sequence_controls.{suffix}",
                    dpi=dpi if suffix == "png" else None,
                    bbox_inches="tight",
                )
            plt.close(fig)

    if not positions.empty:
        part = positions[
            positions["experiment"].eq("full_sequence")
            & positions["validation"].eq("receptor_grouped_5fold")
        ].nsmallest(30, "rank_within_model_validation")
        if not part.empty:
            part = part.sort_values("mean_absolute_coefficient")
            fig, ax = plt.subplots(figsize=(7.5, 7.0))
            y = np.arange(len(part))
            ax.hlines(y, 0, part["mean_absolute_coefficient"], linewidth=1)
            ax.plot(part["mean_absolute_coefficient"], y, "o")
            ax.set_yticks(y, part["generic_position"])
            ax.set_xlabel("Mean absolute coefficient across folds and classes")
            ax.set_title("Full-sequence coefficient localisation")
            fig.tight_layout()
            for suffix in ("svg", "pdf", "png"):
                fig.savefig(
                    figures / f"sequence_position_localisation.{suffix}",
                    dpi=dpi if suffix == "png" else None,
                    bbox_inches="tight",
                )
            plt.close(fig)


def aggregate_results(
    output_dir: Path,
    bootstraps: int,
    coefficient_permutations: int,
    seed: int,
) -> dict[str, pd.DataFrame]:
    payloads = []
    for path in sorted((output_dir / "checkpoints").glob("*.json")):
        payload = parse_checkpoint(path)
        if payload is not None:
            payloads.append(payload)

    primary_performance = []
    primary_predictions = []
    audits = []
    coefficients = []
    neighbours = []
    control_performance = []
    control_predictions = []

    for payload in payloads:
        model = payload.get("model", "")
        if "performance" in payload:
            record = dict(payload["performance"])
            record["task_id"] = payload.get("task_id")
            record["runtime_seconds"] = payload.get("runtime_seconds")
            if model in MODEL_NAMES or model == "global_identity_knn" or model.endswith("_single_transducer"):
                primary_performance.append(record)
            else:
                control_performance.append(record)
        pred = prediction_frame(payload)
        if not pred.empty:
            if model in MODEL_NAMES or model == "global_identity_knn" or model.endswith("_single_transducer"):
                primary_predictions.append(pred)
            else:
                control_predictions.append(pred)
        if payload.get("audit"):
            local = pd.DataFrame(payload["audit"])
            local["experiment"] = model
            local["validation"] = payload.get("validation")
            local["task_id"] = payload.get("task_id")
            audits.append(local)
        if payload.get("coefficients"):
            local = pd.DataFrame(payload["coefficients"])
            coefficients.append(local)
        if payload.get("neighbours"):
            local = pd.DataFrame(payload["neighbours"])
            local["validation"] = payload.get("validation")
            neighbours.append(local)

    performance = pd.DataFrame(primary_performance)
    predictions = pd.concat(primary_predictions, ignore_index=True) if primary_predictions else pd.DataFrame()
    audit = pd.concat(audits, ignore_index=True) if audits else pd.DataFrame()
    coefficient_frame = pd.concat(coefficients, ignore_index=True) if coefficients else pd.DataFrame()
    neighbour_frame = pd.concat(neighbours, ignore_index=True) if neighbours else pd.DataFrame()
    controls = pd.DataFrame(control_performance)
    controls_predictions = pd.concat(control_predictions, ignore_index=True) if control_predictions else pd.DataFrame()

    if not performance.empty and not predictions.empty:
        ci_rows = []
        for (experiment, validation), part in predictions.groupby(["experiment", "validation"], sort=True):
            ci = receptor_bootstrap_metrics(
                part,
                bootstraps,
                seed + int(hashlib.sha256(f"{experiment}|{validation}|model_ci".encode()).hexdigest()[:8], 16) % 1_000_000,
            )
            ci_rows.append({"experiment": experiment, "validation": validation, **ci})
        ci_frame = pd.DataFrame(ci_rows)
        performance = performance.drop(
            columns=[
                column
                for column in (
                    "macro_roc_auc_ci_low", "macro_roc_auc_ci_high",
                    "log_loss_ci_low", "log_loss_ci_high",
                    "bootstrap_valid_iterations",
                )
                if column in performance
            ],
            errors="ignore",
        ).merge(ci_frame, on=["experiment", "validation"], how="left")
    atomic_tsv(output_dir / "metrics/model_metrics.tsv", performance)
    atomic_tsv(output_dir / "predictions/out_of_fold_predictions.tsv", predictions)
    atomic_tsv(output_dir / "fold_manifests/convergence_and_leakage_audit.tsv", audit)
    atomic_tsv(output_dir / "coefficients/model_coefficients.tsv", coefficient_frame)
    atomic_tsv(output_dir / "phylogeny/global_identity_neighbours.tsv", neighbour_frame)
    atomic_tsv(output_dir / "controls/control_model_metrics.tsv", controls)
    atomic_tsv(output_dir / "controls/control_predictions.tsv", controls_predictions)

    position_summary = coefficient_position_summary(coefficient_frame)
    position_enrichment = interface_coefficient_enrichment(
        position_summary, coefficient_permutations, seed + 9000
    )
    atomic_tsv(output_dir / "coefficients/position_coefficient_summary.tsv", position_summary)
    atomic_tsv(output_dir / "coefficients/interface_coefficient_enrichment.tsv", position_enrichment)

    paired_rows = []
    comparisons = (
        ("interface_sequence", "coarse_taxonomy"),
        ("interface_sequence", "global_identity_knn"),
        ("interface_sequence", "noninterface_sequence"),
        ("interface_plus_taxonomy", "coarse_taxonomy"),
        ("full_sequence", "endpoint_geometry"),
        ("sequence_plus_geometry", "full_sequence"),
        ("full_sequence", "interface_sequence"),
    )
    if not predictions.empty:
        for validation in sorted(predictions["validation"].dropna().unique()):
            for candidate, reference in comparisons:
                candidate_frame = predictions[
                    predictions["experiment"].eq(candidate)
                    & predictions["validation"].eq(validation)
                ]
                reference_frame = predictions[
                    predictions["experiment"].eq(reference)
                    & predictions["validation"].eq(validation)
                ]
                if candidate_frame.empty or reference_frame.empty:
                    continue
                result = receptor_bootstrap_metric_difference(
                    candidate_frame,
                    reference_frame,
                    bootstraps,
                    seed
                    + int(hashlib.sha256(f'{validation}|{candidate}|{reference}'.encode()).hexdigest()[:8], 16) % 1_000_000,
                )
                paired_rows.append(
                    {
                        "comparison": f"{candidate}_minus_{reference}",
                        "candidate_model": candidate,
                        "reference_model": reference,
                        "validation": validation,
                        **result,
                    }
                )
    paired = pd.DataFrame(paired_rows)
    atomic_tsv(output_dir / "metrics/paired_model_differences.tsv", paired)

    control_summary_rows = []
    if not controls.empty:
        for (experiment, validation), part in controls.groupby(
            ["experiment", "validation"]
        ):
            values = pd.to_numeric(part["macro_roc_auc"], errors="coerce").dropna()
            control_summary_rows.append(
                {
                    "experiment": experiment,
                    "validation": validation,
                    "valid_iterations": int(len(values)),
                    "mean_auc": float(values.mean()) if len(values) else np.nan,
                    "median_auc": float(values.median()) if len(values) else np.nan,
                    "auc_ci_low": float(values.quantile(0.025)) if len(values) else np.nan,
                    "auc_ci_high": float(values.quantile(0.975)) if len(values) else np.nan,
                }
            )
    control_summary = pd.DataFrame(control_summary_rows)
    atomic_tsv(output_dir / "controls/control_summary.tsv", control_summary)

    control_test_rows = []
    if not controls.empty:
        for validation in sorted(controls["validation"].dropna().unique()):
            observed_row = controls[
                controls["experiment"].eq("interface_fixed_l2")
                & controls["validation"].eq(validation)
            ]
            if observed_row.empty:
                continue
            observed = float(observed_row["macro_roc_auc"].iloc[0])
            for null_model in (
                "matched_noninterface_control",
                "within_receptor_position_shuffle",
                "position_shuffle_across_receptors",
            ):
                null = pd.to_numeric(
                    controls.loc[
                        controls["experiment"].eq(null_model)
                        & controls["validation"].eq(validation),
                        "macro_roc_auc",
                    ],
                    errors="coerce",
                ).dropna()
                if len(null):
                    p_value = (1 + int((null >= observed).sum())) / (len(null) + 1)
                    percentile = float((null < observed).mean())
                else:
                    p_value = np.nan
                    percentile = np.nan
                control_test_rows.append(
                    {
                        "validation": validation,
                        "observed_model": "interface_fixed_l2",
                        "observed_auc": observed,
                        "null_model": null_model,
                        "null_iterations": int(len(null)),
                        "null_mean_auc": float(null.mean()) if len(null) else np.nan,
                        "null_ci_low": float(null.quantile(0.025)) if len(null) else np.nan,
                        "null_ci_high": float(null.quantile(0.975)) if len(null) else np.nan,
                        "empirical_p_value_one_sided": p_value,
                        "observed_percentile_in_null": percentile,
                    }
                )
    control_tests = pd.DataFrame(control_test_rows)
    atomic_tsv(output_dir / "controls/interface_control_tests.tsv", control_tests)

    return {
        "performance": performance,
        "predictions": predictions,
        "audit": audit,
        "coefficients": coefficient_frame,
        "position_summary": position_summary,
        "position_enrichment": position_enrichment,
        "paired": paired,
        "controls": controls,
        "control_summary": control_summary,
        "control_tests": control_tests,
    }


def interpretation_gate(
    performance: pd.DataFrame,
    paired: pd.DataFrame,
    controls: pd.DataFrame,
) -> dict[str, Any]:
    def metric(model: str, validation: str) -> float:
        row = performance[
            performance["experiment"].eq(model)
            & performance["validation"].eq(validation)
        ]
        return float(row["macro_roc_auc"].iloc[0]) if len(row) else np.nan

    receptor_full = metric("full_sequence", "receptor_grouped_5fold")
    cluster_full = metric("full_sequence", "sequence_cluster_30_grouped_5fold")
    cluster_interface = metric(
        "interface_sequence", "sequence_cluster_30_grouped_5fold"
    )
    cluster_knn = metric(
        "global_identity_knn", "sequence_cluster_30_grouped_5fold"
    )
    cluster_taxonomy = metric(
        "coarse_taxonomy", "sequence_cluster_30_grouped_5fold"
    )

    def paired_record(name: str, validation: str) -> dict[str, Any] | None:
        row = paired[
            paired["comparison"].eq(name)
            & paired["validation"].eq(validation)
        ]
        return row.iloc[0].to_dict() if len(row) else None

    interface_vs_knn = paired_record(
        "interface_sequence_minus_global_identity_knn",
        "sequence_cluster_30_grouped_5fold",
    )
    interface_vs_noninterface = paired_record(
        "interface_sequence_minus_noninterface_sequence",
        "sequence_cluster_30_grouped_5fold",
    )
    full_vs_geometry = paired_record(
        "full_sequence_minus_endpoint_geometry", "receptor_grouped_5fold"
    )

    strong_interface_support = bool(
        interface_vs_knn
        and interface_vs_noninterface
        and np.isfinite(cluster_interface)
        and cluster_interface > 0.60
        and interface_vs_knn.get("auc_difference_ci_low", -np.inf) > 0
        and interface_vs_noninterface.get("auc_difference_ci_low", -np.inf) > 0
    )
    phylogeny_dominated = bool(
        np.isfinite(cluster_full)
        and cluster_full < 0.58
        and np.isfinite(cluster_knn)
        and abs(cluster_full - cluster_knn) < 0.04
    )
    mixed_support = bool(
        not strong_interface_support
        and not phylogeny_dominated
        and np.isfinite(cluster_interface)
        and cluster_interface > 0.55
    )

    if strong_interface_support:
        outcome = "strong_interface_support"
        wording = (
            "Coupling-interface sequence carries transferable bound-class information beyond global sequence-identity and non-interface controls."
        )
    elif phylogeny_dominated:
        outcome = "phylogeny_dominated"
        wording = (
            "The sequence advantage weakens under low-identity transfer and is largely consistent with evolutionary organisation rather than a general interface code."
        )
    elif mixed_support:
        outcome = "mixed_support"
        wording = (
            "Evolutionary relatedness explains part of sequence performance, with limited residual information localised to the coupling interface."
        )
    else:
        outcome = "no_resolved_interface_increment"
        wording = (
            "The present controls do not resolve coupling-interface chemistry as an incremental transferable signal beyond phylogeny and non-interface sequence."
        )

    return {
        "outcome": outcome,
        "recommended_wording": wording,
        "key_metrics": {
            "full_sequence_receptor_grouped_auc": receptor_full,
            "full_sequence_30pct_cluster_auc": cluster_full,
            "interface_sequence_30pct_cluster_auc": cluster_interface,
            "global_identity_knn_30pct_cluster_auc": cluster_knn,
            "coarse_taxonomy_30pct_cluster_auc": cluster_taxonomy,
        },
        "key_paired_results": {
            "interface_vs_global_identity_knn_30pct": interface_vs_knn,
            "interface_vs_noninterface_30pct": interface_vs_noninterface,
            "full_sequence_vs_geometry_receptor_grouped": full_vs_geometry,
        },
        "interpretation_limits": [
            "The interface region is physically defined, but coefficient localisation remains descriptive.",
            "Bound-complex labels are structural observations, whereas receptor-only functional labels require separate provenance review.",
            "The current representative frame already contains one structure per receptor-by-transducer-family unit; promiscuous receptors cannot be collapsed to one multiclass row without changing the target to multilabel.",
            "A negative low-identity result does not show that interface chemistry is biologically irrelevant.",
        ],
    }


def main() -> int:
    args = parse_args()
    project = args.project_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    analysis_frame_path = (
        project
        / "results/microswitches/publication_controls/publication_control_analysis_frame.csv"
    )
    config_dir = project / "config"
    if not analysis_frame_path.exists():
        raise FileNotFoundError(analysis_frame_path)
    if not config_dir.exists():
        raise FileNotFoundError(config_dir)

    from gpcr_icl2.microswitches.config import ProjectConfig
    from gpcr_icl2.microswitches.publication_controls import (
        ExperimentSpec,
        FeatureSpec,
        add_generic_sequence_clusters,
    )
    import gpcr_icl2.microswitches.publication_controls as publication_controls

    config = ProjectConfig.from_dir(config_dir)
    source_path = Path(inspect.getfile(publication_controls)).resolve()
    frame = pd.read_csv(analysis_frame_path, low_memory=False)
    frame = frame[frame["transducer_family"].isin(PRIMARY_CLASSES)].copy()
    frame["pdb_id"] = frame["pdb_id"].astype(str).str.upper()
    frame["analysis_row_id"] = (
        frame["pdb_id"].astype(str)
        + "|"
        + frame.get("chain_id", pd.Series("", index=frame.index)).astype(str)
        + "|"
        + frame["receptor_name"].astype(str)
        + "|"
        + frame["transducer_family"].astype(str)
    )
    seq_cols = sequence_columns(frame)
    frame = coerce_sequence(frame, seq_cols)
    if not all(column in frame for column in ("sequence_cluster_30", "sequence_cluster_40", "sequence_cluster_50")):
        frame, cluster_assignments = add_generic_sequence_clusters(frame)
    else:
        cluster_assignments = frame[
            [
                "receptor_name",
                "sequence_cluster_30",
                "sequence_cluster_40",
                "sequence_cluster_50",
            ]
        ].drop_duplicates()
    atomic_tsv(output_dir / "audit/sequence_cluster_assignments.tsv", cluster_assignments)

    frame, specs, availability = feature_specs(frame, config)
    frame_pickle = output_dir / "analysis_frame.pkl"
    frame.to_pickle(frame_pickle)
    atomic_tsv(output_dir / "audit/feature_availability.tsv", availability)
    seq_stats = entropy_and_coverage(frame, sequence_columns(frame))
    atomic_tsv(output_dir / "audit/sequence_position_statistics.tsv", seq_stats)

    interface_columns = list(specs["interface_sequence"].categorical)
    noninterface_columns = list(specs["noninterface_sequence"].categorical)
    control_sets, control_match_table = matched_noninterface_sets(
        frame,
        interface_columns,
        noninterface_columns,
        args.control_sets if args.mode == "full" else (5 if args.mode == "quick" else 1),
        args.seed,
    )
    atomic_tsv(output_dir / "controls/matched_noninterface_position_sets.tsv", control_match_table)

    unit_audit, receptor_audit = build_granularity_audit(frame)
    atomic_tsv(output_dir / "granularity/receptor_transducer_unit_multiplicity.tsv", unit_audit)
    atomic_tsv(output_dir / "granularity/receptor_level_multiplicity.tsv", receptor_audit)

    bound_provenance, functional_provenance, provenance_note = build_label_provenance(project, frame)
    atomic_tsv(output_dir / "label_provenance/bound_label_provenance.tsv", bound_provenance)
    atomic_tsv(output_dir / "label_provenance/functional_label_provenance.tsv", functional_provenance)
    (output_dir / "label_provenance/README.md").write_text(
        "# Label provenance audit\n\n" + provenance_note + "\n",
        encoding="utf-8",
    )

    if args.mode == "smoke":
        # Retain enough receptors for all three classes while keeping the test short.
        selected = []
        for label in PRIMARY_CLASSES:
            receptors = sorted(
                frame.loc[frame["transducer_family"].eq(label), "receptor_name"].unique()
            )[:20]
            selected.extend(receptors)
        frame = frame[frame["receptor_name"].isin(sorted(set(selected)))].copy()
        frame.to_pickle(frame_pickle)
        validations = {"receptor": VALIDATION_MAP["receptor"]}
        selected_models = (
            "full_sequence",
            "interface_sequence",
            "coarse_taxonomy",
            "endpoint_geometry",
        )
        scrambles = 1
    elif args.mode == "quick":
        validations = {
            key: VALIDATION_MAP[key] for key in ("receptor", "seq30")
        }
        selected_models = MODEL_NAMES
        scrambles = min(5, args.scrambles)
    else:
        validations = VALIDATION_MAP
        selected_models = MODEL_NAMES
        scrambles = args.scrambles

    run_config = {
        "analysis": "sequence_interface_phylogeny_granularity_audit_v1",
        "mode": args.mode,
        "project_dir": str(project),
        "analysis_frame": str(analysis_frame_path),
        "publication_controls_source": str(source_path),
        "publication_controls_sha256": sha256_file(source_path),
        "analysis_frame_sha256": sha256_file(analysis_frame_path),
        "config_digest": config.digest,
        "seed": args.seed,
        "workers": args.workers,
        "bootstraps": args.bootstraps,
        "control_sets": len(control_sets),
        "scrambles": scrambles,
        "interface_positions": list(INTERFACE_POSITIONS),
        "posthoc_top_interface_positions": list(POSTHOC_TOP_INTERFACE_POSITIONS),
        "validations": {
            key: {"validation": value[0], "group_column": value[1]}
            for key, value in validations.items()
        },
        "models": list(selected_models) + ["global_identity_knn"],
        "model_specification": "nested L2 logistic regression with lbfgs and training-fold C tuning",
        "primary_label_definition": "transducer family physically present in deposited bound complex",
        "granularity_note": "The current frame has one row per receptor-by-transducer-family unit.",
    }
    config_path = output_dir / "run_config.json"
    if config_path.exists():
        existing = json.loads(config_path.read_text(encoding="utf-8"))
        immutable = (
            "publication_controls_sha256",
            "analysis_frame_sha256",
            "seed",
            "interface_positions",
            "validations",
            "models",
        )
        mismatches = [
            key for key in immutable if existing.get(key) != run_config.get(key)
        ]
        if mismatches:
            raise RuntimeError(
                "Existing output directory has incompatible settings: "
                + ", ".join(mismatches)
            )
    else:
        atomic_json(config_path, run_config)

    tasks: list[dict[str, Any]] = []
    for validation_key, (validation, group_column) in validations.items():
        for model in selected_models:
            identifier = f"primary__{model}__{validation_key}"
            checkpoint = checkpoint_dir / f"{identifier}.json"
            if parse_checkpoint(checkpoint) is not None:
                continue
            spec = specs[model]
            tasks.append(
                {
                    "kind": "primary",
                    "task_id": identifier,
                    "model": model,
                    "validation": validation,
                    "group_column": group_column,
                    "feature_spec": {
                        "name": spec.name,
                        "numeric": list(spec.numeric),
                        "categorical": list(spec.categorical),
                        "description": spec.description,
                    },
                    "frame_pickle": str(frame_pickle),
                    "config_dir": str(config_dir),
                    "seed": args.seed,
                    "checkpoint": str(checkpoint),
                }
            )
        identifier = f"knn__{validation_key}"
        checkpoint = checkpoint_dir / f"{identifier}.json"
        if parse_checkpoint(checkpoint) is None:
            tasks.append(
                {
                    "kind": "knn",
                    "task_id": identifier,
                    "validation": validation,
                    "group_column": group_column,
                    "frame_pickle": str(frame_pickle),
                    "seed": args.seed,
                    "checkpoint": str(checkpoint),
                }
            )

    # Natural-granularity sensitivity: receptors observed with only one class.
    single_receptors = receptor_audit.loc[
        receptor_audit["single_transducer_receptor"], "receptor_name"
    ].astype(str)
    single_frame = frame[frame["receptor_name"].astype(str).isin(single_receptors)].copy()
    single_pickle = output_dir / "single_transducer_analysis_frame.pkl"
    single_frame.to_pickle(single_pickle)
    for base_model in ("full_sequence", "interface_sequence", "endpoint_geometry"):
        identifier = f"single__{base_model}"
        checkpoint = checkpoint_dir / f"{identifier}.json"
        if parse_checkpoint(checkpoint) is not None:
            continue
        spec = specs[base_model]
        tasks.append(
            {
                "kind": "primary",
                "task_id": identifier,
                "model": f"{base_model}_single_transducer",
                "validation": "single_transducer_receptor_grouped_5fold",
                "group_column": "receptor_name",
                "frame_subset": "single_transducer_receptors",
                "feature_spec": {
                    "name": spec.name,
                    "numeric": list(spec.numeric),
                    "categorical": list(spec.categorical),
                    "description": spec.description,
                },
                "frame_pickle": str(single_pickle),
                "config_dir": str(config_dir),
                "seed": args.seed,
                "checkpoint": str(checkpoint),
            }
        )

    # Fast, repeated controls are restricted to receptor-grouped and the most
    # stringent 30% cluster transfer.
    control_validations = [
        ("receptor",) + VALIDATION_MAP["receptor"],
        ("seq30",) + VALIDATION_MAP["seq30"],
    ]
    for validation_key, validation, group_column in control_validations:
        identifier = f"interface_fixed_l2__{validation_key}"
        checkpoint = checkpoint_dir / f"{identifier}.json"
        if parse_checkpoint(checkpoint) is None:
            tasks.append(
                {
                    "kind": "control",
                    "task_id": identifier,
                    "model": "interface_fixed_l2",
                    "validation": validation,
                    "group_column": group_column,
                    "columns": interface_columns,
                    "iteration": -1,
                    "transform": "none",
                    "frame_pickle": str(frame_pickle),
                    "seed": args.seed,
                    "checkpoint": str(checkpoint),
                }
            )
    for iteration, columns in enumerate(control_sets):
        for validation_key, validation, group_column in control_validations:
            identifier = f"control_noninterface__{iteration:03d}__{validation_key}"
            checkpoint = checkpoint_dir / f"{identifier}.json"
            if parse_checkpoint(checkpoint) is None:
                tasks.append(
                    {
                        "kind": "control",
                        "task_id": identifier,
                        "model": "matched_noninterface_control",
                        "validation": validation,
                        "group_column": group_column,
                        "columns": columns,
                        "iteration": iteration,
                        "transform": "none",
                        "frame_pickle": str(frame_pickle),
                        "seed": args.seed + iteration,
                        "checkpoint": str(checkpoint),
                    }
                )
    for iteration in range(scrambles):
        for transform, model in (
            ("within_receptor_position_shuffle", "within_receptor_position_shuffle"),
            ("position_shuffle_across_receptors", "position_shuffle_across_receptors"),
        ):
            for validation_key, validation, group_column in control_validations:
                identifier = f"{model}__{iteration:03d}__{validation_key}"
                checkpoint = checkpoint_dir / f"{identifier}.json"
                if parse_checkpoint(checkpoint) is None:
                    tasks.append(
                        {
                            "kind": "control",
                            "task_id": identifier,
                            "model": model,
                            "validation": validation,
                            "group_column": group_column,
                            "columns": interface_columns,
                            "iteration": iteration,
                            "transform": transform,
                            "frame_pickle": str(frame_pickle),
                            "seed": args.seed + 10000 + iteration,
                            "checkpoint": str(checkpoint),
                        }
                    )

    print(
        json.dumps(
            {
                "mode": args.mode,
                "rows": len(frame),
                "receptors": int(frame["receptor_name"].nunique()),
                "sequence_positions": len(sequence_columns(frame)),
                "interface_positions_available": len(interface_columns),
                "noninterface_positions_available": len(noninterface_columns),
                "remaining_tasks": len(tasks),
                "workers": min(args.workers, max(1, len(tasks))),
            },
            indent=2,
        ),
        flush=True,
    )

    failures = []
    if tasks:
        context = mp.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(args.workers, len(tasks)), mp_context=context
        ) as executor:
            future_map = {}
            for task in tasks:
                if task["kind"] == "primary":
                    future = executor.submit(run_existing_model_task, task)
                elif task["kind"] == "knn":
                    future = executor.submit(run_knn_task, task)
                else:
                    future = executor.submit(run_fast_control_task, task)
                future_map[future] = task
            for index, future in enumerate(
                concurrent.futures.as_completed(future_map), start=1
            ):
                task = future_map[future]
                try:
                    payload = future.result()
                except Exception as error:
                    payload = {
                        "status": "failed",
                        "task_id": task["task_id"],
                        "error": repr(error),
                    }
                if payload.get("status") != "ok":
                    failures.append(payload)
                print(
                    f"[{index}/{len(tasks)}] {payload.get('task_id')} "
                    f"status={payload.get('status')} "
                    f"runtime_s={payload.get('runtime_seconds', float('nan')):.1f} "
                    f"error={payload.get('error', '')}",
                    flush=True,
                )

    results = aggregate_results(
        output_dir,
        args.bootstraps if args.mode == "full" else min(300, args.bootstraps),
        args.coefficient_permutations
        if args.mode == "full"
        else min(1000, args.coefficient_permutations),
        args.seed,
    )
    save_figures(
        results["performance"],
        results["paired"],
        results["controls"],
        results["position_summary"],
        output_dir,
        args.dpi,
    )
    gate = interpretation_gate(
        results["performance"], results["paired"], results["controls"]
    )
    atomic_json(output_dir / "interpretation_gate.json", gate)
    summary_lines = [
        "# Sequence-interface and phylogeny audit",
        "",
        f"Outcome: **{gate['outcome']}**",
        "",
        gate["recommended_wording"],
        "",
        "## Design",
        "",
        f"* {len(frame)} receptor-by-transducer-family rows from {frame['receptor_name'].nunique()} receptors.",
        f"* {len(interface_columns)} physically defined interface positions and {len(noninterface_columns)} non-interface positions.",
        "* Receptor-grouped and 30%, 40% and 50% sequence-cluster transfer tests in full mode.",
        "* Nested L2 logistic regression with fold-local C tuning for the primary model comparison.",
        "* Global sequence-identity kNN and coarse taxonomy baselines.",
        "* Matched non-interface, positional scrambling and receptor-specific barcode scrambling controls.",
        "",
        "## Granularity conclusion",
        "",
        "The publication-control frame already contains one structure per receptor-by-transducer-family unit.  Receptors observed with more than one family are represented by different bound complexes, whereas sequence is identical across those rows.  A separate single-transducer-receptor sensitivity analysis is included.",
        "",
        "## Label provenance",
        "",
        "Bound-class labels are direct assignments from the transducer physically present in each deposited complex.  Receptor-only functional labels are audited separately, and the current consensus table does not resolve direct measurement versus homology inference for every family label.",
        "",
        "## Interpretation limits",
        "",
    ]
    summary_lines.extend(f"* {item}" for item in gate["interpretation_limits"])
    (output_dir / "run_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    failed_checkpoints = []
    for path in checkpoint_dir.glob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if payload.get("status") != "ok":
            failed_checkpoints.append(path.name)

    all_primary_converged = True
    audit = results["audit"]
    if not audit.empty and "converged" in audit:
        primary_audit = audit[audit["experiment"].isin(list(selected_models))]
        if not primary_audit.empty:
            converged_values = primary_audit["converged"].map(
                lambda value: bool(value) if pd.notna(value) else False
            )
            all_primary_converged = bool(converged_values.all())

    manifest = {
        "success": not failed_checkpoints and all_primary_converged,
        "completed_utc": pd.Timestamp.utcnow().isoformat(),
        "output_dir": str(output_dir),
        "mode": args.mode,
        "failed_checkpoint_count": len(failed_checkpoints),
        "failed_checkpoints": failed_checkpoints,
        "all_primary_model_folds_converged": all_primary_converged,
        "interpretation_outcome": gate["outcome"],
        "key_output_files": [
            "metrics/model_metrics.tsv",
            "metrics/paired_model_differences.tsv",
            "predictions/out_of_fold_predictions.tsv",
            "coefficients/position_coefficient_summary.tsv",
            "coefficients/interface_coefficient_enrichment.tsv",
            "controls/control_summary.tsv",
            "controls/interface_control_tests.tsv",
            "granularity/receptor_transducer_unit_multiplicity.tsv",
            "label_provenance/bound_label_provenance.tsv",
            "label_provenance/functional_label_provenance.tsv",
            "interpretation_gate.json",
            "run_summary.md",
        ],
    }
    atomic_json(output_dir / "RUN_MANIFEST.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)
    return 0 if manifest["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
