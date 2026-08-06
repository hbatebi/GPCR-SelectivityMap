from __future__ import annotations

"""Phylogeny-aware evaluation of GPCR sequence and static structure representations.

The routines in this module are intentionally conservative. They distinguish
association under receptor-grouped validation from transfer to low-homology
sequence clusters and include a global sequence-identity kNN baseline. Outputs
are representation benchmarks, not validated functional coupling probabilities.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
import json
import math
import re

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .io import read_table, write_json, write_table

DEFAULT_INTERFACE_POSITIONS: tuple[str, ...] = (
    "3x49", "3x50", "3x53", "3x55",
    "5x58", "5x61", "5x64", "5x68", "5x72",
    "6x23", "6x30", "6x34", "6x37", "6x40",
    "7x53", "7x54", "7x55", "7x56",
    "8x47", "8x50", "8x53", "8x59",
    "34x50", "34x51", "34x52", "34x53", "34x54", "34x55", "34x56", "34x57",
)

DEFAULT_C_GRID: tuple[float, ...] = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0)


@dataclass(frozen=True)
class RepresentationSpec:
    name: str
    columns: tuple[str, ...]
    kind: str = "logistic"  # logistic or identity_knn


class UnionFind:
    def __init__(self, items: Iterable[str]):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x: str) -> str:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def normalise_position(column: str) -> str:
    match = re.search(r"gpcrdb_([^_]+)_expected_aa", str(column))
    return match.group(1) if match else str(column)


def position_column(position: str) -> str:
    return f"gpcrdb_{position}_expected_aa"


def sequence_columns(frame: pd.DataFrame) -> list[str]:
    return sorted(
        column for column in frame.columns
        if column.startswith("gpcrdb_") and column.endswith("_expected_aa")
    )


def load_interface_positions(path: str | Path | None = None) -> tuple[str, ...]:
    if path is None:
        return DEFAULT_INTERFACE_POSITIONS
    table = read_table(path)
    for candidate in ("generic_position", "position", "gpcrdb_position"):
        if candidate in table.columns:
            values = table[candidate].dropna().astype(str).str.strip()
            return tuple(dict.fromkeys(v for v in values if v))
    raise ValueError("Interface-position table needs a generic_position column")


def _clean_aa(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().upper()
    return "" if text in {"", "NAN", "NONE", "NA"} else text


def receptor_signatures(
    frame: pd.DataFrame,
    *,
    receptor_column: str,
    columns: Sequence[str],
) -> dict[str, tuple[str, ...]]:
    signatures: dict[str, tuple[str, ...]] = {}
    for receptor, part in frame.groupby(receptor_column, sort=True):
        values: list[str] = []
        for column in columns:
            observed = part[column].map(_clean_aa)
            observed = observed[observed.ne("")]
            values.append(observed.mode().iloc[0] if len(observed) else "")
        signatures[str(receptor)] = tuple(values)
    return signatures


def pair_identity(a: Sequence[str], b: Sequence[str]) -> float:
    comparable = [(x, y) for x, y in zip(a, b) if x and y]
    if not comparable:
        return 0.0
    return float(np.mean([x == y for x, y in comparable]))


def sequence_identity_matrix(signatures: dict[str, tuple[str, ...]]) -> pd.DataFrame:
    receptors = sorted(signatures)
    matrix = np.eye(len(receptors), dtype=float)
    for i, left in enumerate(receptors):
        for j in range(i + 1, len(receptors)):
            value = pair_identity(signatures[left], signatures[receptors[j]])
            matrix[i, j] = matrix[j, i] = value
    return pd.DataFrame(matrix, index=receptors, columns=receptors)


def assign_identity_clusters(
    identity: pd.DataFrame,
    thresholds: Sequence[float],
) -> pd.DataFrame:
    receptors = list(identity.index.astype(str))
    output = pd.DataFrame({"receptor_name": receptors})
    for threshold in thresholds:
        uf = UnionFind(receptors)
        for i, left in enumerate(receptors):
            for j in range(i + 1, len(receptors)):
                right = receptors[j]
                if float(identity.loc[left, right]) >= float(threshold):
                    uf.union(left, right)
        roots = {receptor: uf.find(receptor) for receptor in receptors}
        root_order = {root: idx for idx, root in enumerate(sorted(set(roots.values())))}
        key = int(round(float(threshold) * 100))
        output[f"sequence_cluster_{key}"] = [
            f"seq{key}_{root_order[roots[receptor]]:04d}" for receptor in receptors
        ]
    return output


def _macro_auc(y: np.ndarray, prob: np.ndarray, classes: np.ndarray) -> float:
    return float(roc_auc_score(y, prob, labels=classes, multi_class="ovr", average="macro"))


def _build_pipeline(frame: pd.DataFrame, columns: Sequence[str], c_value: float) -> Pipeline:
    numeric = [c for c in columns if pd.api.types.is_numeric_dtype(frame[c])]
    categorical = [c for c in columns if c not in numeric]
    transformers = []
    if numeric:
        transformers.append((
            "numeric",
            Pipeline([
                ("impute", SimpleImputer(strategy="median", add_indicator=True)),
                ("scale", StandardScaler()),
            ]),
            numeric,
        ))
    if categorical:
        transformers.append((
            "categorical",
            Pipeline([
                ("impute", SimpleImputer(strategy="most_frequent")),
                ("onehot", OneHotEncoder(handle_unknown="ignore")),
            ]),
            categorical,
        ))
    if not transformers:
        raise ValueError("No usable columns for model")
    return Pipeline([
        ("preprocess", ColumnTransformer(transformers, remainder="drop")),
        ("model", LogisticRegression(max_iter=30000, solver="lbfgs", C=float(c_value))),
    ])


def _valid_split(y_train: np.ndarray, y_test: np.ndarray, classes: np.ndarray) -> bool:
    return set(classes).issubset(set(y_train)) and len(set(y_test)) >= 2


def _outer_splits(
    frame: pd.DataFrame,
    *,
    target: str,
    group_column: str,
    folds: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    y = frame[target].astype(str).to_numpy()
    groups = frame[group_column].astype(str).to_numpy()
    n_groups = pd.Series(groups).nunique()
    n_splits = min(int(folds), int(n_groups))
    if n_splits < 2:
        raise ValueError(f"Need at least two unique groups in {group_column}")
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return [(train, test) for train, test in splitter.split(frame, y, groups)]


def _choose_c(
    frame: pd.DataFrame,
    columns: Sequence[str],
    *,
    target: str,
    receptor_column: str,
    train_indices: np.ndarray,
    c_grid: Sequence[float],
    seed: int,
) -> float:
    train_frame = frame.iloc[train_indices].reset_index(drop=True)
    y = train_frame[target].astype(str).to_numpy()
    groups = train_frame[receptor_column].astype(str).to_numpy()
    classes = np.array(sorted(np.unique(y)))
    n_groups = pd.Series(groups).nunique()
    n_splits = min(3, int(n_groups))
    if n_splits < 2 or len(classes) < 2:
        return 1.0
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    best_c, best_score = 1.0, -np.inf
    for c_value in c_grid:
        predictions: list[dict[str, object]] = []
        for inner_fold, (inner_train, inner_test) in enumerate(splitter.split(train_frame, y, groups)):
            if not _valid_split(y[inner_train], y[inner_test], classes):
                continue
            pipe = _build_pipeline(train_frame, columns, float(c_value))
            try:
                pipe.fit(train_frame.iloc[inner_train][list(columns)], y[inner_train])
            except Exception:
                continue
            raw = pipe.predict_proba(train_frame.iloc[inner_test][list(columns)])
            model_classes = pipe.named_steps["model"].classes_
            aligned = np.zeros((len(inner_test), len(classes)), dtype=float)
            for j, label in enumerate(model_classes):
                if label in classes:
                    aligned[:, np.where(classes == label)[0][0]] = raw[:, j]
            for idx, probs in zip(inner_test, aligned):
                predictions.append({"index": int(idx), "actual": y[idx], "probs": probs})
        if not predictions:
            continue
        actual = np.array([r["actual"] for r in predictions], dtype=str)
        probs = np.vstack([r["probs"] for r in predictions])
        try:
            score = _macro_auc(actual, probs, classes)
        except ValueError:
            continue
        if score > best_score + 1e-12:
            best_score, best_c = score, float(c_value)
    return best_c


def evaluate_logistic_representation(
    frame: pd.DataFrame,
    spec: RepresentationSpec,
    *,
    target: str,
    receptor_column: str,
    group_column: str,
    validation_name: str,
    folds: int,
    seed: int,
    c_grid: Sequence[float] = DEFAULT_C_GRID,
    tune_c: bool = True,
    splits: list[tuple[np.ndarray, np.ndarray]] | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    columns = [c for c in spec.columns if c in frame.columns]
    if not columns:
        raise ValueError(f"No columns found for {spec.name}")
    classes = np.array(sorted(frame[target].dropna().astype(str).unique()))
    y = frame[target].astype(str).to_numpy()
    outer = splits or _outer_splits(frame, target=target, group_column=group_column, folds=folds, seed=seed)
    records: list[dict[str, object]] = []
    convergence = True
    selected_cs: list[float] = []
    valid_folds = 0
    for fold, (train, test) in enumerate(outer):
        if not _valid_split(y[train], y[test], classes):
            continue
        c_value = _choose_c(
            frame, columns, target=target, receptor_column=receptor_column,
            train_indices=train, c_grid=c_grid, seed=seed + fold,
        ) if tune_c else float(c_grid[0])
        selected_cs.append(c_value)
        pipe = _build_pipeline(frame, columns, c_value)
        try:
            pipe.fit(frame.iloc[train][columns], y[train])
            raw = pipe.predict_proba(frame.iloc[test][columns])
        except Exception:
            convergence = False
            continue
        model_classes = pipe.named_steps["model"].classes_
        aligned = np.zeros((len(test), len(classes)), dtype=float)
        for j, label in enumerate(model_classes):
            if label in classes:
                aligned[:, np.where(classes == label)[0][0]] = raw[:, j]
        predicted = classes[np.argmax(aligned, axis=1)]
        valid_folds += 1
        for row_index, actual, pred, probs in zip(test, y[test], predicted, aligned):
            record: dict[str, object] = {
                "model": spec.name,
                "validation": validation_name,
                "row_index": int(row_index),
                "fold": int(fold),
                "receptor_name": str(frame.iloc[row_index][receptor_column]),
                "actual": str(actual),
                "predicted": str(pred),
            }
            record.update({f"probability_{label}": float(probs[j]) for j, label in enumerate(classes)})
            records.append(record)
    predictions = pd.DataFrame(records)
    if predictions.empty:
        raise ValueError(f"No valid folds for {spec.name} under {validation_name}")
    pcols = [f"probability_{c}" for c in classes]
    metrics = {
        "experiment": spec.name,
        "validation": validation_name,
        "n": int(len(frame)),
        "n_predictions": int(len(predictions)),
        "n_receptors": int(frame[receptor_column].nunique()),
        "n_features": int(len(columns)),
        "valid_folds": int(valid_folds),
        "all_folds_converged": bool(convergence and valid_folds == len(outer)),
        "selected_C_median": float(np.median(selected_cs)) if selected_cs else np.nan,
        "macro_roc_auc": _macro_auc(predictions.actual.to_numpy(), predictions[pcols].to_numpy(), classes),
        "log_loss": float(log_loss(predictions.actual, predictions[pcols], labels=classes)),
        "balanced_accuracy": float(balanced_accuracy_score(predictions.actual, predictions.predicted)),
    }
    return predictions, metrics


def evaluate_identity_knn(
    frame: pd.DataFrame,
    identity: pd.DataFrame,
    *,
    target: str,
    receptor_column: str,
    group_column: str,
    validation_name: str,
    folds: int,
    seed: int,
    k: int = 5,
    splits: list[tuple[np.ndarray, np.ndarray]] | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    classes = np.array(sorted(frame[target].dropna().astype(str).unique()))
    y = frame[target].astype(str).to_numpy()
    outer = splits or _outer_splits(frame, target=target, group_column=group_column, folds=folds, seed=seed)
    records: list[dict[str, object]] = []
    valid_folds = 0
    for fold, (train, test) in enumerate(outer):
        if not _valid_split(y[train], y[test], classes):
            continue
        train_frame = frame.iloc[train]
        receptor_labels = train_frame.groupby(receptor_column)[target].apply(lambda x: sorted(set(x.astype(str))))
        train_receptors = sorted(receptor_labels.index.astype(str))
        if not train_receptors:
            continue
        valid_folds += 1
        for row_index in test:
            receptor = str(frame.iloc[row_index][receptor_column])
            neighbours = sorted(
                ((other, float(identity.loc[receptor, other])) for other in train_receptors if other != receptor),
                key=lambda item: item[1], reverse=True,
            )[: max(1, int(k))]
            votes = {label: 0.0 for label in classes}
            if neighbours:
                for neighbour, similarity in neighbours:
                    labels = receptor_labels.loc[neighbour]
                    weight = max(similarity, 1e-8) / max(1, len(labels))
                    for label in labels:
                        if label in votes:
                            votes[label] += weight
            total = sum(votes.values())
            probs = np.array([votes[c] / total if total else 1.0 / len(classes) for c in classes])
            predicted = str(classes[int(np.argmax(probs))])
            record: dict[str, object] = {
                "model": "global_identity_knn",
                "validation": validation_name,
                "row_index": int(row_index),
                "fold": int(fold),
                "receptor_name": receptor,
                "actual": str(y[row_index]),
                "predicted": predicted,
            }
            record.update({f"probability_{label}": float(probs[j]) for j, label in enumerate(classes)})
            records.append(record)
    predictions = pd.DataFrame(records)
    if predictions.empty:
        raise ValueError(f"No valid folds for identity kNN under {validation_name}")
    pcols = [f"probability_{c}" for c in classes]
    metrics = {
        "experiment": "global_identity_knn",
        "validation": validation_name,
        "n": int(len(frame)),
        "n_predictions": int(len(predictions)),
        "n_receptors": int(frame[receptor_column].nunique()),
        "n_features": int(identity.shape[1]),
        "valid_folds": int(valid_folds),
        "all_folds_converged": True,
        "selected_C_median": np.nan,
        "macro_roc_auc": _macro_auc(predictions.actual.to_numpy(), predictions[pcols].to_numpy(), classes),
        "log_loss": float(log_loss(predictions.actual, predictions[pcols], labels=classes)),
        "balanced_accuracy": float(balanced_accuracy_score(predictions.actual, predictions.predicted)),
    }
    return predictions, metrics


def receptor_bootstrap_interval(
    predictions: pd.DataFrame,
    *,
    iterations: int,
    seed: int,
) -> tuple[float, float, int]:
    classes = np.array(sorted(predictions.actual.astype(str).unique()))
    pcols = [f"probability_{c}" for c in classes]
    receptors = predictions.receptor_name.astype(str).unique()
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(int(iterations)):
        chosen = rng.choice(receptors, size=len(receptors), replace=True)
        parts = [predictions[predictions.receptor_name.astype(str).eq(receptor)] for receptor in chosen]
        sample = pd.concat(parts, ignore_index=True)
        try:
            values.append(_macro_auc(sample.actual.to_numpy(), sample[pcols].to_numpy(), classes))
        except ValueError:
            continue
    if not values:
        return np.nan, np.nan, 0
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)), len(values)


def paired_bootstrap_difference(
    candidate: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    iterations: int,
    seed: int,
) -> dict[str, float | int]:
    key = ["validation", "row_index", "receptor_name", "actual"]
    merged = candidate.merge(reference, on=key, suffixes=("_candidate", "_reference"))
    classes = np.array(sorted(merged.actual.astype(str).unique()))
    cand_cols = [f"probability_{c}_candidate" for c in classes]
    ref_cols = [f"probability_{c}_reference" for c in classes]
    observed = _macro_auc(merged.actual.to_numpy(), merged[cand_cols].to_numpy(), classes) - _macro_auc(
        merged.actual.to_numpy(), merged[ref_cols].to_numpy(), classes
    )
    receptors = merged.receptor_name.astype(str).unique()
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(int(iterations)):
        chosen = rng.choice(receptors, size=len(receptors), replace=True)
        sample = pd.concat([merged[merged.receptor_name.astype(str).eq(r)] for r in chosen], ignore_index=True)
        try:
            values.append(
                _macro_auc(sample.actual.to_numpy(), sample[cand_cols].to_numpy(), classes)
                - _macro_auc(sample.actual.to_numpy(), sample[ref_cols].to_numpy(), classes)
            )
        except ValueError:
            continue
    return {
        "auc_difference": float(observed),
        "auc_difference_ci_low": float(np.quantile(values, 0.025)) if values else np.nan,
        "auc_difference_ci_high": float(np.quantile(values, 0.975)) if values else np.nan,
        "bootstrap_valid_iterations": int(len(values)),
        "n_rows": int(len(merged)),
        "n_receptors": int(merged.receptor_name.nunique()),
    }


def _resolve_columns(frame: pd.DataFrame, prefixes: Sequence[str], explicit: Sequence[str] = ()) -> list[str]:
    columns = [c for c in explicit if c in frame.columns]
    for prefix in prefixes:
        columns.extend([c for c in frame.columns if c.startswith(prefix)])
    return sorted(set(columns))


def build_default_representations(
    frame: pd.DataFrame,
    *,
    interface_positions: Sequence[str],
    geometry_prefixes: Sequence[str],
    coarse_taxonomy_columns: Sequence[str],
    fine_taxonomy_columns: Sequence[str],
) -> dict[str, RepresentationSpec]:
    seq = sequence_columns(frame)
    interface = [position_column(p) for p in interface_positions if position_column(p) in seq]
    noninterface = [c for c in seq if c not in interface]
    geometry = _resolve_columns(frame, geometry_prefixes)
    coarse = [c for c in coarse_taxonomy_columns if c in frame.columns]
    fine = [c for c in fine_taxonomy_columns if c in frame.columns]
    specs: dict[str, RepresentationSpec] = {
        "full_sequence": RepresentationSpec("full_sequence", tuple(seq)),
        "interface_sequence": RepresentationSpec("interface_sequence", tuple(interface)),
        "noninterface_sequence": RepresentationSpec("noninterface_sequence", tuple(noninterface)),
    }
    if geometry:
        specs["endpoint_geometry"] = RepresentationSpec("endpoint_geometry", tuple(geometry))
        specs["sequence_plus_geometry"] = RepresentationSpec("sequence_plus_geometry", tuple(seq + geometry))
    if coarse:
        specs["coarse_taxonomy"] = RepresentationSpec("coarse_taxonomy", tuple(coarse))
        if interface:
            specs["interface_plus_taxonomy"] = RepresentationSpec("interface_plus_taxonomy", tuple(interface + coarse))
    if fine:
        specs["fine_family"] = RepresentationSpec("fine_family", tuple(fine))
    specs["global_identity_knn"] = RepresentationSpec("global_identity_knn", tuple(seq), kind="identity_knn")
    return specs


def _interpret(metrics: pd.DataFrame, paired: pd.DataFrame) -> dict[str, object]:
    def metric(model: str, validation: str) -> float:
        rows = metrics[(metrics.experiment == model) & (metrics.validation == validation)]
        return float(rows.iloc[0].macro_roc_auc) if len(rows) else np.nan

    seq30 = metric("full_sequence", "sequence_cluster_30_grouped")
    interface30 = metric("interface_sequence", "sequence_cluster_30_grouped")
    knn30 = metric("global_identity_knn", "sequence_cluster_30_grouped")
    outcome = "insufficient_evidence"
    wording = "The available analyses do not yet distinguish transferable representation signal from evolutionary organisation."
    rows_a = paired[paired.comparison.eq("interface_sequence_minus_noninterface_sequence") & paired.validation.eq("sequence_cluster_30_grouped")]
    rows_b = paired[paired.comparison.eq("interface_sequence_minus_global_identity_knn") & paired.validation.eq("sequence_cluster_30_grouped")]
    null_a = len(rows_a) and float(rows_a.iloc[0].auc_difference_ci_low) <= 0 <= float(rows_a.iloc[0].auc_difference_ci_high)
    null_b = len(rows_b) and float(rows_b.iloc[0].auc_difference_ci_low) <= 0 <= float(rows_b.iloc[0].auc_difference_ci_high)
    if np.isfinite(seq30) and seq30 < 0.60 and null_a and null_b:
        outcome = "phylogeny_dominated"
        wording = (
            "The sequence advantage weakens under low-identity transfer and is largely consistent "
            "with evolutionary organisation rather than a general interface code."
        )
    elif np.isfinite(interface30) and interface30 >= 0.65 and np.isfinite(knn30) and interface30 > knn30 + 0.05:
        outcome = "transferable_interface_signal"
        wording = "Interface sequence retains information beyond global sequence proximity under low-homology transfer."
    return {
        "outcome": outcome,
        "recommended_wording": wording,
        "interpretation_limits": [
            "Above-chance receptor-grouped performance is not evidence of low-homology transfer.",
            "A negative low-identity result does not show that interface chemistry is biologically irrelevant.",
            "Bound-complex labels and receptor-only functional labels require separate provenance treatment.",
            "Scores are representation benchmarks and not validated functional coupling probabilities.",
        ],
    }



def plot_transfer_summary(metrics: pd.DataFrame, output_dir: str | Path) -> None:
    import matplotlib.pyplot as plt

    order = ["receptor_grouped", "sequence_cluster_50_grouped", "sequence_cluster_40_grouped", "sequence_cluster_30_grouped"]
    labels = {
        "receptor_grouped": "Receptor grouped",
        "sequence_cluster_50_grouped": "50% identity",
        "sequence_cluster_40_grouped": "40% identity",
        "sequence_cluster_30_grouped": "30% identity",
    }
    selected = ["full_sequence", "interface_sequence", "endpoint_geometry", "global_identity_knn"]
    available_order = [v for v in order if v in set(metrics.validation)]
    if not available_order:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for model in selected:
        part = metrics[metrics.experiment.eq(model)].set_index("validation")
        x, y, low, high = [], [], [], []
        for idx, validation in enumerate(available_order):
            if validation not in part.index:
                continue
            row = part.loc[validation]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            x.append(idx)
            y.append(float(row.macro_roc_auc))
            low.append(float(row.get("macro_roc_auc_ci_low", np.nan)))
            high.append(float(row.get("macro_roc_auc_ci_high", np.nan)))
        if not x:
            continue
        yerr = None
        if all(np.isfinite(low)) and all(np.isfinite(high)):
            yerr = np.vstack([np.array(y) - np.array(low), np.array(high) - np.array(y)])
        ax.errorbar(x, y, yerr=yerr, marker="o", capsize=3, label=model.replace("_", " "))
    ax.axhline(0.5, linewidth=1, linestyle="--")
    ax.set_xticks(range(len(available_order)), [labels[v] for v in available_order])
    ax.set_ylabel("Macro ROC AUC")
    ax.set_ylim(0.35, 1.0)
    ax.set_title("Transfer of GPCR sequence and static structure representations")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "representation_transfer.svg", bbox_inches="tight")
    fig.savefig(out / "representation_transfer.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def run_cluster_benchmark(
    *,
    features: str | Path,
    target: str,
    receptor_column: str,
    output: str | Path,
    thresholds: Sequence[float] = (0.30, 0.40, 0.50),
    folds: int = 5,
    seed: int = 20272729,
    bootstraps: int = 500,
    interface_positions: str | Path | None = None,
    geometry_prefixes: Sequence[str] = ("scv_", "geometry_", "cavity_", "microswitch_", "intracellular_"),
    coarse_taxonomy_columns: Sequence[str] = ("receptor_ligand_class", "receptor_subfamily"),
    fine_taxonomy_columns: Sequence[str] = ("receptor_family",),
    models: Sequence[str] | None = None,
    k_neighbors: int = 5,
) -> pd.DataFrame:
    frame = read_table(features).dropna(subset=[target, receptor_column]).reset_index(drop=True)
    if frame.empty:
        raise ValueError("No rows remain after removing missing target and receptor values")
    seq = sequence_columns(frame)
    if not seq:
        raise ValueError("No gpcrdb_*_expected_aa sequence columns were found")
    positions = load_interface_positions(interface_positions)
    signatures = receptor_signatures(frame, receptor_column=receptor_column, columns=seq)
    identity = sequence_identity_matrix(signatures)
    clusters = assign_identity_clusters(identity, thresholds).rename(columns={"receptor_name": receptor_column})
    cluster_columns = [c for c in clusters.columns if c != receptor_column]
    frame = frame.drop(columns=[c for c in cluster_columns if c in frame.columns], errors="ignore")
    frame = frame.merge(clusters, on=receptor_column, how="left", validate="many_to_one")
    specs = build_default_representations(
        frame,
        interface_positions=positions,
        geometry_prefixes=geometry_prefixes,
        coarse_taxonomy_columns=coarse_taxonomy_columns,
        fine_taxonomy_columns=fine_taxonomy_columns,
    )
    if models:
        requested = set(models)
        specs = {name: spec for name, spec in specs.items() if name in requested}
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    write_table(clusters, out / "audit" / "sequence_cluster_assignments.tsv")
    write_table(identity.reset_index(names=receptor_column), out / "audit" / "sequence_identity_matrix.tsv")

    validations: list[tuple[str, str]] = [("receptor_grouped", receptor_column)]
    for threshold in thresholds:
        key = int(round(float(threshold) * 100))
        validations.append((f"sequence_cluster_{key}_grouped", f"sequence_cluster_{key}"))

    all_predictions: list[pd.DataFrame] = []
    all_metrics: list[dict[str, object]] = []
    for validation_name, group_column in validations:
        splits = _outer_splits(frame, target=target, group_column=group_column, folds=folds, seed=seed)
        for name, spec in specs.items():
            try:
                if spec.kind == "identity_knn":
                    pred, metric_row = evaluate_identity_knn(
                        frame, identity, target=target, receptor_column=receptor_column,
                        group_column=group_column, validation_name=validation_name,
                        folds=folds, seed=seed, k=k_neighbors, splits=splits,
                    )
                else:
                    pred, metric_row = evaluate_logistic_representation(
                        frame, spec, target=target, receptor_column=receptor_column,
                        group_column=group_column, validation_name=validation_name,
                        folds=folds, seed=seed, splits=splits,
                    )
            except ValueError:
                continue
            lo, hi, valid = receptor_bootstrap_interval(pred, iterations=bootstraps, seed=seed)
            metric_row.update({
                "macro_roc_auc_ci_low": lo,
                "macro_roc_auc_ci_high": hi,
                "bootstrap_valid_iterations": valid,
            })
            all_predictions.append(pred)
            all_metrics.append(metric_row)

    predictions = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame()
    metrics = pd.DataFrame(all_metrics)
    if metrics.empty:
        raise RuntimeError("No benchmark model completed")
    pairs = [
        ("full_sequence", "endpoint_geometry"),
        ("interface_sequence", "noninterface_sequence"),
        ("interface_sequence", "global_identity_knn"),
        ("sequence_plus_geometry", "full_sequence"),
    ]
    paired_rows: list[dict[str, object]] = []
    for validation_name, _ in validations:
        for candidate, reference in pairs:
            cand = predictions[(predictions.model == candidate) & (predictions.validation == validation_name)]
            ref = predictions[(predictions.model == reference) & (predictions.validation == validation_name)]
            if cand.empty or ref.empty:
                continue
            result = paired_bootstrap_difference(cand, ref, iterations=bootstraps, seed=seed + 101)
            result.update({
                "comparison": f"{candidate}_minus_{reference}",
                "candidate_model": candidate,
                "reference_model": reference,
                "validation": validation_name,
            })
            paired_rows.append(result)
    paired = pd.DataFrame(paired_rows)
    gate = _interpret(metrics, paired)

    write_table(metrics.sort_values(["validation", "macro_roc_auc"], ascending=[True, False]), out / "metrics" / "model_metrics.tsv")
    write_table(predictions, out / "predictions" / "out_of_fold_predictions.tsv")
    write_table(paired, out / "metrics" / "paired_model_differences.tsv")
    write_json(out / "interpretation_gate.json", gate)
    plot_transfer_summary(metrics, out / "figures")
    write_json(out / "analysis_manifest.json", {
        "analysis": "phylogeny_aware_representation_benchmark",
        "target": target,
        "receptor_column": receptor_column,
        "thresholds": list(map(float, thresholds)),
        "folds": int(folds),
        "seed": int(seed),
        "bootstraps": int(bootstraps),
        "interface_positions": list(positions),
        "models": sorted(specs),
        "scope": "Representation benchmarking only; outputs are not validated functional coupling probabilities.",
    })
    return metrics


def run_phylogeny_baseline(**kwargs: object) -> pd.DataFrame:
    kwargs = dict(kwargs)
    kwargs["models"] = ("global_identity_knn", "coarse_taxonomy", "fine_family")
    return run_cluster_benchmark(**kwargs)
