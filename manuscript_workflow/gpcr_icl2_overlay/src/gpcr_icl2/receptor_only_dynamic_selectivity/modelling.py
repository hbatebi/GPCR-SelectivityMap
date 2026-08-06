from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
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
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from gpcr_icl2.microswitches.config import ProjectConfig
from gpcr_icl2.microswitches.publication_controls import FeatureSpec, build_feature_catalog

PRIMARY_CLASSES = ("Gs", "Gi/o", "Gq/11")


@dataclass(frozen=True)
class ModelSpec:
    name: str
    numeric: tuple[str, ...] = ()
    categorical: tuple[str, ...] = ()
    description: str = ""


def _coerce_boolean_columns(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    for column in output.columns:
        series = output[column]
        if str(series.dtype) == "boolean" or series.dtype == bool:
            output[column] = series.astype("Float64")
            continue
        if series.dtype == object:
            values = set(series.dropna().astype(str).str.lower().unique())
            if values and values.issubset({"true", "false"}):
                output[column] = series.map(
                    lambda value: np.nan
                    if pd.isna(value)
                    else float(str(value).lower() == "true")
                )
    return output


def _numeric(frame: pd.DataFrame, columns: Sequence[str]) -> tuple[str, ...]:
    selected = []
    for column in columns:
        if column not in frame:
            continue
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.notna().any() and values.dropna().nunique() > 1:
            selected.append(column)
    return tuple(dict.fromkeys(selected))


def _categorical(frame: pd.DataFrame, columns: Sequence[str]) -> tuple[str, ...]:
    selected = []
    for column in columns:
        if column not in frame:
            continue
        values = frame[column].dropna().astype(str)
        if values.nunique() > 1:
            selected.append(column)
    return tuple(dict.fromkeys(selected))


def build_model_specs(frame: pd.DataFrame, config_dir: str | Path) -> tuple[pd.DataFrame, dict[str, ModelSpec], pd.DataFrame]:
    config = ProjectConfig.from_dir(Path(config_dir))
    prepared, catalog, availability = build_feature_catalog(
        _coerce_boolean_columns(frame),
        config,
        minimum_feature_coverage=0.75,
    )

    contact = _numeric(prepared, [column for column in prepared if column.startswith("contact_")])
    surface_electro = _numeric(
        prepared,
        [
            column
            for column in prepared
            if column.startswith("surface_") or column.startswith("electro_")
        ],
    )
    mechanical = _numeric(prepared, [column for column in prepared if column.startswith("mechanical_")])
    availability_controls = _numeric(prepared, [column for column in prepared if column.startswith("availability_")])
    all_new = tuple(dict.fromkeys(contact + surface_electro + mechanical))

    def from_feature_spec(name: str, spec: FeatureSpec, description: str | None = None) -> ModelSpec:
        return ModelSpec(
            name=name,
            numeric=tuple(spec.numeric),
            categorical=tuple(spec.categorical),
            description=description or spec.description,
        )

    existing = from_feature_spec(
        "existing_geometry",
        catalog["SCVcoreE_full_adjusted"],
        "Existing activation, intracellular geometry, core cavity and adjusted metadata",
    )
    sequence = from_feature_spec(
        "sequence_only",
        catalog["Q_sequence_fingerprint"],
        "Generic-position receptor sequence fingerprint",
    )
    missingness = from_feature_spec(
        "missingness_only",
        catalog["M_missingness_only"],
        "Existing coordinate availability and measurement completeness controls",
    )

    specs = {
        existing.name: existing,
        sequence.name: sequence,
        missingness.name: missingness,
        "contact_network": ModelSpec(
            "contact_network", contact, (), "Receptor-chain contact and side-chain network descriptors"
        ),
        "surface_electrostatics": ModelSpec(
            "surface_electrostatics",
            surface_electro,
            (),
            "Intracellular surface chemistry and approximate Coulombic electrostatics",
        ),
        "mechanical_susceptibility": ModelSpec(
            "mechanical_susceptibility",
            mechanical,
            (),
            "Generic-position anisotropic-network mechanical susceptibility",
        ),
        "all_new_static": ModelSpec(
            "all_new_static",
            all_new,
            (),
            "All new receptor-chain-only static physicochemical and mechanical features",
        ),
        "sequence_plus_new_static": ModelSpec(
            "sequence_plus_new_static",
            all_new,
            sequence.categorical,
            "Sequence plus all new receptor-chain-only static features",
        ),
        "existing_geometry_plus_new_static": ModelSpec(
            "existing_geometry_plus_new_static",
            tuple(dict.fromkeys(existing.numeric + all_new)),
            existing.categorical,
            "Existing adjusted endpoint structure plus all new static features",
        ),
        "new_feature_availability_only": ModelSpec(
            "new_feature_availability_only",
            availability_controls,
            (),
            "Availability of the new feature blocks only",
        ),
    }

    extra_availability = []
    for model_name, spec in specs.items():
        for feature in spec.numeric + spec.categorical:
            extra_availability.append(
                {
                    "feature_block": model_name,
                    "feature": feature,
                    "feature_type": "numeric" if feature in spec.numeric else "categorical",
                    "n_available": int(prepared[feature].notna().sum()) if feature in prepared else 0,
                    "fraction_available": float(prepared[feature].notna().mean()) if feature in prepared else 0.0,
                    "n_unique": int(prepared[feature].dropna().astype(str).nunique()) if feature in prepared else 0,
                }
            )
    availability = pd.concat([availability, pd.DataFrame(extra_availability)], ignore_index=True)
    return prepared, specs, availability


def _preprocessor(numeric: Sequence[str], categorical: Sequence[str]) -> ColumnTransformer:
    transformers = []
    if numeric:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
                        ("scaler", StandardScaler()),
                    ]
                ),
                list(numeric),
            )
        )
    if categorical:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent", keep_empty_features=True)),
                        ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                    ]
                ),
                list(categorical),
            )
        )
    if not transformers:
        raise ValueError("No usable features")
    return ColumnTransformer(transformers, remainder="drop", verbose_feature_names_out=True)


def _pipeline(numeric: Sequence[str], categorical: Sequence[str], c_value: float, seed: int, max_iter: int) -> Pipeline:
    return Pipeline(
        [
            ("preprocess", _preprocessor(numeric, categorical)),
            (
                "model",
                LogisticRegression(
                    solver="lbfgs",
                    l1_ratio=0.0,
                    C=float(c_value),
                    class_weight="balanced",
                    max_iter=int(max_iter),
                    tol=1e-4,
                    random_state=int(seed),
                ),
            ),
        ]
    )


def _dynamic_splits(y: np.ndarray, groups: np.ndarray, requested: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    unique_groups = np.unique(groups)
    class_group_counts = [len(np.unique(groups[y == label])) for label in np.unique(y)]
    n_splits = min(requested, len(unique_groups), *(class_group_counts or [len(unique_groups)]))
    if n_splits < 2:
        return []
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(splitter.split(np.zeros(len(y)), y, groups))


def _fold_features(train: pd.DataFrame, spec: ModelSpec, minimum_coverage: float) -> tuple[tuple[str, ...], tuple[str, ...]]:
    numeric = []
    for column in spec.numeric:
        if column not in train:
            continue
        values = pd.to_numeric(train[column], errors="coerce")
        if values.notna().mean() >= minimum_coverage and values.dropna().nunique() > 1:
            numeric.append(column)
    categorical = []
    for column in spec.categorical:
        if column not in train:
            continue
        values = train[column].dropna().astype(str)
        if train[column].notna().mean() >= minimum_coverage and values.nunique() > 1:
            categorical.append(column)
    return tuple(numeric), tuple(categorical)


def _metric_record(y_true: np.ndarray, probabilities: np.ndarray, classes: Sequence[str]) -> dict[str, float]:
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
    auc_values = []
    pr_values = []
    for index, label in enumerate(classes):
        target = (y == label).astype(int)
        if target.min() != target.max():
            auc_values.append(float(roc_auc_score(target, p[:, index])))
            pr_values.append(float(average_precision_score(target, p[:, index])))
    class_index = {label: index for index, label in enumerate(classes)}
    true_probability = np.asarray([p[row, class_index[label]] for row, label in enumerate(y)], dtype=float)
    return {
        "macro_roc_auc": float(np.mean(auc_values)) if auc_values else np.nan,
        "macro_pr_auc": float(np.mean(pr_values)) if pr_values else np.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
        "log_loss": float(-np.mean(np.log(np.clip(true_probability, 1e-15, 1.0)))),
        "matthews_correlation": float(matthews_corrcoef(y, predicted)),
    }


def _tune_c(
    x: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    numeric: Sequence[str],
    categorical: Sequence[str],
    c_grid: Sequence[float],
    seed: int,
    inner_folds: int,
    max_iter: int,
    sample_weights: np.ndarray | None,
) -> tuple[Pipeline, dict[str, Any]]:
    splits = _dynamic_splits(y, groups, inner_folds, seed)
    classes = set(np.unique(y))
    splits = [
        (train, test)
        for train, test in splits
        if set(np.unique(y[train])) == classes and set(np.unique(y[test])) == classes
    ]
    best_c = float(c_grid[len(c_grid) // 2])
    best_score = -np.inf
    if len(splits) >= 2:
        for c_value in c_grid:
            scores = []
            for train_index, test_index in splits:
                candidate = _pipeline(numeric, categorical, c_value, seed, max_iter)
                kwargs = {}
                if sample_weights is not None:
                    kwargs["model__sample_weight"] = sample_weights[train_index]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ConvergenceWarning)
                    candidate.fit(x.iloc[train_index], y[train_index], **kwargs)
                probabilities = candidate.predict_proba(x.iloc[test_index])
                try:
                    score = roc_auc_score(
                        y[test_index],
                        probabilities,
                        labels=candidate.classes_,
                        multi_class="ovr",
                        average="macro",
                    )
                except ValueError:
                    continue
                scores.append(float(score))
            if scores and float(np.mean(scores)) > best_score:
                best_score = float(np.mean(scores))
                best_c = float(c_value)
    fitted = _pipeline(numeric, categorical, best_c, seed, max_iter)
    kwargs = {"model__sample_weight": sample_weights} if sample_weights is not None else {}
    fitted.fit(x, y, **kwargs)
    return fitted, {"C": best_c, "inner_macro_auc": best_score, "inner_valid_folds": len(splits)}


def cross_validate_model(
    frame: pd.DataFrame,
    spec: ModelSpec,
    *,
    seed: int,
    outer_folds: int,
    inner_folds: int,
    c_grid: Sequence[float],
    bootstrap_iterations: int,
    minimum_coverage: float = 0.60,
    max_iter: int = 50000,
) -> dict[str, Any]:
    work = _coerce_boolean_columns(frame).reset_index(drop=True)
    work = work.loc[work["transducer_family"].astype(str).isin(PRIMARY_CLASSES)].reset_index(drop=True)
    y = work["transducer_family"].astype(str).to_numpy()
    groups = work["receptor_name"].astype(str).to_numpy()
    classes = tuple(sorted(np.unique(y)))
    splits = _dynamic_splits(y, groups, outer_folds, seed)
    probabilities = np.full((len(work), len(classes)), np.nan)
    fold_ids = np.full(len(work), -1, dtype=int)
    audit_rows = []
    coefficient_rows = []
    receptor_counts = work.groupby("receptor_name")["transducer_family"].size()
    weights = work["receptor_name"].map(lambda value: 1.0 / receptor_counts.get(value, 1)).to_numpy(float)

    for fold, (train_index, test_index) in enumerate(splits):
        train = work.iloc[train_index]
        test = work.iloc[test_index]
        numeric, categorical = _fold_features(train, spec, minimum_coverage)
        if not numeric and not categorical:
            audit_rows.append({"experiment": spec.name, "fold": fold, "status": "no_usable_features"})
            continue
        columns = list(numeric) + list(categorical)
        fitted, parameters = _tune_c(
            train[columns],
            y[train_index],
            groups[train_index],
            numeric,
            categorical,
            c_grid,
            seed + fold,
            inner_folds,
            max_iter,
            weights[train_index],
        )
        local = fitted.predict_proba(test[columns])
        for local_index, label in enumerate(fitted.classes_):
            probabilities[test_index, classes.index(str(label))] = local[:, local_index]
        fold_ids[test_index] = fold
        model = fitted.named_steps["model"]
        n_iter = int(np.max(np.asarray(model.n_iter_)))
        converged = n_iter < model.max_iter
        names = fitted.named_steps["preprocess"].get_feature_names_out()
        for class_label, coefficients in zip(model.classes_, np.asarray(model.coef_)):
            for name, coefficient in zip(names, coefficients):
                coefficient_rows.append(
                    {
                        "experiment": spec.name,
                        "fold": fold,
                        "class_label": str(class_label),
                        "transformed_feature": str(name),
                        "coefficient": float(coefficient),
                        "absolute_coefficient": abs(float(coefficient)),
                    }
                )
        audit_rows.append(
            {
                "experiment": spec.name,
                "fold": fold,
                "status": "ok",
                "train_n": len(train_index),
                "test_n": len(test_index),
                "train_receptors": int(train["receptor_name"].nunique()),
                "test_receptors": int(test["receptor_name"].nunique()),
                "group_overlap": int(len(set(groups[train_index]).intersection(groups[test_index]))),
                "n_numeric_features": len(numeric),
                "n_categorical_features": len(categorical),
                "n_iter": n_iter,
                "max_iter": model.max_iter,
                "converged": converged,
                "best_parameters": json.dumps(parameters, sort_keys=True),
                "selected_numeric_features": ";".join(numeric),
                "selected_categorical_features": ";".join(categorical),
            }
        )

    valid = np.isfinite(probabilities).all(axis=1)
    predictions = work[
        [
            column
            for column in (
                "pdb_id",
                "chain_id",
                "preferred_chain",
                "receptor_name",
                "receptor_family",
                "receptor_ligand_class",
                "transducer_family",
            )
            if column in work
        ]
    ].copy()
    predictions["sample_id"] = (
        predictions["pdb_id"].astype(str)
        + "|"
        + predictions["receptor_name"].astype(str)
        + "|"
        + predictions["transducer_family"].astype(str)
    )
    predictions["experiment"] = spec.name
    predictions["outer_fold"] = fold_ids
    for index, label in enumerate(classes):
        predictions[f"probability_{label}"] = probabilities[:, index]
    predicted = np.full(len(work), "", dtype=object)
    predicted[valid] = np.asarray(classes)[np.argmax(probabilities[valid], axis=1)]
    predictions["predicted_class"] = predicted

    metrics = _metric_record(y, probabilities, classes)
    audit = pd.DataFrame(audit_rows)
    ok = audit.loc[audit.get("status", pd.Series(dtype=str)).eq("ok")] if not audit.empty else pd.DataFrame()
    performance = {
        "experiment": spec.name,
        "description": spec.description,
        "n": int(len(work)),
        "n_receptors": int(work["receptor_name"].nunique()),
        "n_predictions": int(valid.sum()),
        "valid_folds": int(len(ok)),
        "all_folds_converged": bool(not ok.empty and ok["converged"].fillna(False).astype(bool).all()),
        **metrics,
    }

    bootstrap = receptor_cluster_bootstrap(predictions, bootstrap_iterations, seed + 10000)
    for metric in (
        "macro_roc_auc",
        "macro_pr_auc",
        "balanced_accuracy",
        "log_loss",
        "matthews_correlation",
    ):
        values = pd.to_numeric(bootstrap.get(metric, pd.Series(dtype=float)), errors="coerce").dropna()
        performance[f"{metric}_ci_low"] = float(np.percentile(values, 2.5)) if len(values) else np.nan
        performance[f"{metric}_ci_high"] = float(np.percentile(values, 97.5)) if len(values) else np.nan

    return {
        "performance": performance,
        "predictions": predictions,
        "audit": audit,
        "coefficients": pd.DataFrame(coefficient_rows),
        "bootstrap": bootstrap,
    }


def receptor_cluster_bootstrap(predictions: pd.DataFrame, iterations: int, seed: int) -> pd.DataFrame:
    if iterations <= 0:
        return pd.DataFrame()
    valid = predictions.loc[predictions["predicted_class"].ne("")].reset_index(drop=True)
    if valid.empty:
        return pd.DataFrame()
    classes = tuple(sorted(valid["transducer_family"].astype(str).unique()))
    pcols = [f"probability_{label}" for label in classes]
    receptors = valid["receptor_name"].astype(str).unique()
    by_receptor = {
        receptor: valid.index[valid["receptor_name"].astype(str).eq(receptor)].to_numpy()
        for receptor in receptors
    }
    rng = np.random.default_rng(seed)
    rows = []
    for iteration in range(iterations):
        sample = rng.choice(receptors, size=len(receptors), replace=True)
        indices = np.concatenate([by_receptor[receptor] for receptor in sample])
        subset = valid.iloc[indices]
        rows.append(
            {
                "iteration": iteration,
                **_metric_record(
                    subset["transducer_family"].astype(str).to_numpy(),
                    subset[pcols].to_numpy(float),
                    classes,
                ),
            }
        )
    return pd.DataFrame(rows)


def paired_model_difference(
    predictions: pd.DataFrame,
    baseline: str,
    candidate: str,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    classes = PRIMARY_CLASSES
    pcols = [f"probability_{label}" for label in classes]
    keys = ["sample_id", "pdb_id", "receptor_name", "transducer_family"]
    left = predictions.loc[predictions["experiment"].eq(baseline), keys + pcols].copy()
    right = predictions.loc[predictions["experiment"].eq(candidate), keys + pcols].copy()
    merged = left.merge(right, on=keys, suffixes=("_baseline", "_candidate"))
    if merged.empty:
        return {"baseline": baseline, "candidate": candidate, "status": "no_paired_rows"}

    def evaluate(subset: pd.DataFrame) -> tuple[float, float]:
        y = subset["transducer_family"].astype(str).to_numpy()
        p0 = np.column_stack([subset[f"probability_{label}_baseline"] for label in classes])
        p1 = np.column_stack([subset[f"probability_{label}_candidate"] for label in classes])
        valid = np.isfinite(p0).all(axis=1) & np.isfinite(p1).all(axis=1)
        y, p0, p1 = y[valid], p0[valid], p1[valid]
        if len(y) == 0:
            return np.nan, np.nan
        metric0 = _metric_record(y, p0, classes)
        metric1 = _metric_record(y, p1, classes)
        delta_auc = float(metric1["macro_roc_auc"] - metric0["macro_roc_auc"])
        delta_log_loss = float(metric1["log_loss"] - metric0["log_loss"])
        return delta_auc, delta_log_loss

    observed_auc, observed_log_loss = evaluate(merged)
    receptors = merged["receptor_name"].astype(str).unique()
    rng = np.random.default_rng(seed)
    boot_auc = []
    boot_log = []
    for _ in range(iterations):
        sample = rng.choice(receptors, size=len(receptors), replace=True)
        subset = pd.concat(
            [merged.loc[merged["receptor_name"].astype(str).eq(receptor)] for receptor in sample],
            ignore_index=True,
        )
        delta_auc, delta_log = evaluate(subset)
        if np.isfinite(delta_auc):
            boot_auc.append(delta_auc)
        if np.isfinite(delta_log):
            boot_log.append(delta_log)
    return {
        "baseline": baseline,
        "candidate": candidate,
        "status": "ok",
        "n": int(len(merged)),
        "n_receptors": int(len(receptors)),
        "delta_auc_candidate_minus_baseline": observed_auc,
        "delta_auc_ci_low": float(np.percentile(boot_auc, 2.5)) if boot_auc else np.nan,
        "delta_auc_ci_high": float(np.percentile(boot_auc, 97.5)) if boot_auc else np.nan,
        "delta_log_loss_candidate_minus_baseline": observed_log_loss,
        "delta_log_loss_ci_low": float(np.percentile(boot_log, 2.5)) if boot_log else np.nan,
        "delta_log_loss_ci_high": float(np.percentile(boot_log, 97.5)) if boot_log else np.nan,
        "auc_improvement_supported": bool(boot_auc and np.percentile(boot_auc, 2.5) > 0),
        "log_loss_improvement_supported": bool(boot_log and np.percentile(boot_log, 97.5) < 0),
    }


def grouped_prediction_permutation(
    predictions: pd.DataFrame,
    iterations: int,
    seed: int,
) -> dict[str, Any]:
    valid = predictions.loc[predictions["predicted_class"].ne("")].reset_index(drop=True)
    if valid.empty or iterations <= 0:
        return {"status": "not_run", "iterations": 0}
    classes = PRIMARY_CLASSES
    pcols = [f"probability_{label}" for label in classes]
    probabilities = valid[pcols].to_numpy(float)
    observed = _metric_record(valid["transducer_family"].astype(str).to_numpy(), probabilities, classes)["macro_roc_auc"]
    groups = []
    for receptor, part in valid.groupby("receptor_name", sort=True):
        ordered = part.sort_values("sample_id")
        groups.append((str(receptor), ordered.index.to_numpy(), ordered["transducer_family"].astype(str).to_numpy()))
    strata: dict[int, list[tuple[str, np.ndarray, np.ndarray]]] = {}
    for item in groups:
        strata.setdefault(len(item[1]), []).append(item)
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(iterations):
        permuted = valid["transducer_family"].astype(str).to_numpy().copy()
        for size, items in strata.items():
            donor_order = rng.permutation(len(items))
            for target_index, donor_index in enumerate(donor_order):
                target_rows = items[target_index][1]
                donor_labels = items[int(donor_index)][2]
                permuted[target_rows] = donor_labels
        score = _metric_record(permuted, probabilities, classes)["macro_roc_auc"]
        if np.isfinite(score):
            null.append(score)
    p_value = (1 + sum(value >= observed for value in null)) / (1 + len(null)) if null else np.nan
    return {
        "status": "ok" if null else "not_testable",
        "iterations": int(len(null)),
        "observed_macro_roc_auc": observed,
        "null_mean_macro_roc_auc": float(np.mean(null)) if null else np.nan,
        "null_sd_macro_roc_auc": float(np.std(null, ddof=1)) if len(null) > 1 else np.nan,
        "grouped_permutation_p_value": float(p_value),
        "permutation_unit": "receptor label vector within equal-row-count strata",
    }


def _align_test_columns(test: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    output = test.copy()
    for column in columns:
        if column not in output:
            output[column] = np.nan
    return output


def strict_receptor_transfer(
    bound_frame: pd.DataFrame,
    receptor_only_frame: pd.DataFrame,
    specs: Mapping[str, ModelSpec],
    model_names: Sequence[str],
    *,
    seed: int,
    inner_folds: int,
    c_grid: Sequence[float],
    max_iter: int,
    minimum_coverage: float,
    bootstrap_iterations: int,
    definition: str,
    fixed_c: float = 1.0,
) -> dict[str, pd.DataFrame]:
    bound = _coerce_boolean_columns(bound_frame).copy()
    receptor_only = _coerce_boolean_columns(receptor_only_frame).copy()
    label_columns = [label for label in PRIMARY_CLASSES if label in receptor_only]
    if not label_columns:
        raise ValueError("Receptor-only frame lacks functional Gs/Gi/o/Gq/11 annotations")
    annotated = receptor_only.loc[
        receptor_only[label_columns].apply(pd.to_numeric, errors="coerce").notna().any(axis=1)
    ].copy()
    predictions_rows = []
    audit_rows = []

    for receptor in sorted(annotated["receptor_name"].astype(str).unique()):
        train = bound.loc[bound["receptor_name"].astype(str).ne(receptor)].copy()
        test = annotated.loc[annotated["receptor_name"].astype(str).eq(receptor)].copy()
        y = train["transducer_family"].astype(str).to_numpy()
        groups = train["receptor_name"].astype(str).to_numpy()
        counts = train.groupby("receptor_name")["transducer_family"].size()
        weights = train["receptor_name"].map(lambda value: 1.0 / counts.get(value, 1)).to_numpy(float)

        for model_name in model_names:
            spec = specs[model_name]
            numeric, categorical = _fold_features(train, spec, minimum_coverage)
            columns = list(numeric) + list(categorical)
            if not columns:
                audit_rows.append(
                    {
                        "definition": definition,
                        "receptor_name": receptor,
                        "experiment": model_name,
                        "status": "no_usable_features",
                    }
                )
                continue
            aligned_test = _align_test_columns(test, columns)
            fitted = _pipeline(
                numeric,
                categorical,
                float(fixed_c),
                seed + sum(ord(char) for char in receptor + model_name),
                max_iter,
            )
            fitted.fit(train[columns], y, model__sample_weight=weights)
            parameters = {
                "C": float(fixed_c),
                "selection": "prespecified_fixed_for_held_receptor_transfer",
            }
            probabilities = fitted.predict_proba(aligned_test[columns])
            class_to_index = {str(label): index for index, label in enumerate(fitted.classes_)}
            model = fitted.named_steps["model"]
            audit_rows.append(
                {
                    "definition": definition,
                    "receptor_name": receptor,
                    "experiment": model_name,
                    "status": "ok",
                    "training_receptors": int(train["receptor_name"].nunique()),
                    "held_out_receptor_absent_from_training": True,
                    "test_structures": int(len(test)),
                    "n_numeric_features": len(numeric),
                    "n_categorical_features": len(categorical),
                    "n_iter": int(np.max(model.n_iter_)),
                    "max_iter": int(model.max_iter),
                    "converged": bool(int(np.max(model.n_iter_)) < int(model.max_iter)),
                    "best_parameters": json.dumps(parameters, sort_keys=True),
                }
            )
            for row_index, (_, row) in enumerate(test.iterrows()):
                output = {
                    "definition": definition,
                    "experiment": model_name,
                    "pdb_id": str(row["pdb_id"]),
                    "receptor_name": receptor,
                }
                for label in PRIMARY_CLASSES:
                    output[f"probability_{label}"] = (
                        float(probabilities[row_index, class_to_index[label]])
                        if label in class_to_index
                        else np.nan
                    )
                    output[label] = pd.to_numeric(pd.Series([row.get(label)]), errors="coerce").iloc[0]
                predictions_rows.append(output)

    structure_predictions = pd.DataFrame(predictions_rows)
    if structure_predictions.empty:
        return {
            "structure_predictions": structure_predictions,
            "receptor_predictions": pd.DataFrame(),
            "metrics": pd.DataFrame(),
            "audit": pd.DataFrame(audit_rows),
        }

    aggregation = {
        **{f"probability_{label}": "mean" for label in PRIMARY_CLASSES},
        **{label: "first" for label in PRIMARY_CLASSES},
        "pdb_id": lambda values: ";".join(sorted(set(values.astype(str)))),
    }
    receptor_predictions = (
        structure_predictions.groupby(["definition", "experiment", "receptor_name"], as_index=False)
        .agg(aggregation)
    )

    metric_rows = []
    rng = np.random.default_rng(seed + 80000)
    for (current_definition, experiment), part in receptor_predictions.groupby(["definition", "experiment"]):
        family_aucs = []
        family_prs = []
        for label in PRIMARY_CLASSES:
            truth = pd.to_numeric(part[label], errors="coerce")
            score = pd.to_numeric(part[f"probability_{label}"], errors="coerce")
            valid = truth.notna() & score.notna()
            target = truth[valid].astype(int).to_numpy()
            values = score[valid].to_numpy(float)
            if len(target) >= 3 and len(np.unique(target)) == 2:
                auc = float(roc_auc_score(target, values))
                pr = float(average_precision_score(target, values))
                family_aucs.append(auc)
                family_prs.append(pr)
            else:
                auc = np.nan
                pr = np.nan
            metric_rows.append(
                {
                    "definition": current_definition,
                    "experiment": experiment,
                    "target": label,
                    "n_receptors": int(valid.sum()),
                    "n_positive": int(target.sum()) if len(target) else 0,
                    "n_negative": int(len(target) - target.sum()) if len(target) else 0,
                    "roc_auc": auc,
                    "pr_auc": pr,
                    "bootstrap_ci_low": np.nan,
                    "bootstrap_ci_high": np.nan,
                    "valid_bootstraps": 0,
                }
            )

        # Receptor bootstrap of macro family AUC. Unknown labels remain missing.
        boot = []
        receptor_values = part["receptor_name"].astype(str).to_numpy()
        for _ in range(bootstrap_iterations):
            indices = rng.integers(0, len(part), size=len(part))
            sample = part.iloc[indices]
            values = []
            for label in PRIMARY_CLASSES:
                truth = pd.to_numeric(sample[label], errors="coerce")
                score = pd.to_numeric(sample[f"probability_{label}"], errors="coerce")
                valid = truth.notna() & score.notna()
                target = truth[valid].astype(int).to_numpy()
                if len(target) >= 3 and len(np.unique(target)) == 2:
                    values.append(float(roc_auc_score(target, score[valid].to_numpy(float))))
            if values:
                boot.append(float(np.mean(values)))
        metric_rows.append(
            {
                "definition": current_definition,
                "experiment": experiment,
                "target": "macro_functional",
                "n_receptors": int(part["receptor_name"].nunique()),
                "n_positive": np.nan,
                "n_negative": np.nan,
                "roc_auc": float(np.mean(family_aucs)) if family_aucs else np.nan,
                "pr_auc": float(np.mean(family_prs)) if family_prs else np.nan,
                "bootstrap_ci_low": float(np.percentile(boot, 2.5)) if len(boot) >= 50 else np.nan,
                "bootstrap_ci_high": float(np.percentile(boot, 97.5)) if len(boot) >= 50 else np.nan,
                "valid_bootstraps": int(len(boot)),
            }
        )

    return {
        "structure_predictions": structure_predictions,
        "receptor_predictions": receptor_predictions,
        "metrics": pd.DataFrame(metric_rows),
        "audit": pd.DataFrame(audit_rows),
    }
