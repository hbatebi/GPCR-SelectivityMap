from __future__ import annotations

"""Interface versus non-interface and scrambling controls."""

from pathlib import Path
from typing import Sequence
import numpy as np
import pandas as pd

from .io import read_table, write_json, write_table
from .sequence_audit import (
    DEFAULT_INTERFACE_POSITIONS,
    RepresentationSpec,
    assign_identity_clusters,
    evaluate_logistic_representation,
    load_interface_positions,
    normalise_position,
    pair_identity,
    paired_bootstrap_difference,
    position_column,
    receptor_signatures,
    sequence_columns,
    sequence_identity_matrix,
    _outer_splits,
)


def _position_statistics(frame: pd.DataFrame, columns: Sequence[str], receptor_column: str) -> pd.DataFrame:
    receptor_frame = frame.drop_duplicates(receptor_column)
    rows = []
    for column in columns:
        values = receptor_frame[column].dropna().astype(str)
        frequencies = values.value_counts(normalize=True)
        entropy = float(-(frequencies * np.log(frequencies)).sum()) if len(frequencies) else np.nan
        rows.append({
            "feature": column,
            "generic_position": normalise_position(column),
            "coverage": float(receptor_frame[column].notna().mean()),
            "n_unique": int(values.nunique()),
            "entropy_nats": entropy,
        })
    return pd.DataFrame(rows)


def _matched_noninterface_sets(
    stats: pd.DataFrame,
    interface_columns: Sequence[str],
    noninterface_columns: Sequence[str],
    *,
    n_sets: int,
    seed: int,
) -> list[list[str]]:
    rng = np.random.default_rng(seed)
    lookup = stats.set_index("feature")
    controls: list[list[str]] = []
    for _ in range(int(n_sets)):
        available = list(noninterface_columns)
        selected: list[str] = []
        order = list(interface_columns)
        rng.shuffle(order)
        for interface in order:
            if not available:
                break
            target = lookup.loc[interface]
            distances = []
            for candidate in available:
                row = lookup.loc[candidate]
                distance = abs(float(row.coverage) - float(target.coverage)) + 0.35 * abs(
                    float(row.entropy_nats) - float(target.entropy_nats)
                )
                distance += float(rng.uniform(0, 0.02))
                distances.append((distance, candidate))
            _, chosen = min(distances)
            selected.append(chosen)
            available.remove(chosen)
        controls.append(sorted(selected))
    return controls


def _shuffle_within_rows(frame: pd.DataFrame, columns: Sequence[str], rng: np.random.Generator) -> pd.DataFrame:
    out = frame.copy()
    values = out[list(columns)].to_numpy(dtype=object)
    for row in values:
        rng.shuffle(row)
    out.loc[:, list(columns)] = values
    return out


def _shuffle_positions_across_rows(frame: pd.DataFrame, columns: Sequence[str], rng: np.random.Generator) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        values = out[column].to_numpy(dtype=object).copy()
        rng.shuffle(values)
        out[column] = values
    return out



def _plot_interface_controls(summary: pd.DataFrame, output_dir: str | Path) -> None:
    import matplotlib.pyplot as plt

    if summary.empty:
        return
    validations = list(dict.fromkeys(summary.validation.astype(str)))
    control_types = list(dict.fromkeys(summary.null_model.astype(str)))
    width = 0.8 / max(1, len(control_types))
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for j, control in enumerate(control_types):
        values = []
        for validation in validations:
            rows = summary[(summary.validation == validation) & (summary.null_model == control)]
            values.append(float(rows.iloc[0].null_mean_auc) if len(rows) else np.nan)
        positions = np.arange(len(validations)) - 0.4 + width / 2 + j * width
        ax.bar(positions, values, width=width, label=control.replace("_", " "))
    observed = []
    for validation in validations:
        rows = summary[summary.validation == validation]
        observed.append(float(rows.iloc[0].observed_auc) if len(rows) else np.nan)
    ax.plot(np.arange(len(validations)), observed, marker="o", linewidth=2, label="observed interface")
    ax.axhline(0.5, linewidth=1, linestyle="--")
    ax.set_xticks(np.arange(len(validations)), [v.replace("_grouped", "").replace("_", " ") for v in validations])
    ax.set_ylabel("Macro ROC AUC")
    ax.set_ylim(0.35, 1.0)
    ax.set_title("Predefined interface sequence versus control representations")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "interface_controls.svg", bbox_inches="tight")
    fig.savefig(out / "interface_controls.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def run_interface_audit(
    *,
    features: str | Path,
    target: str,
    receptor_column: str,
    output: str | Path,
    interface_positions: str | Path | None = None,
    identity_threshold: float = 0.30,
    folds: int = 5,
    seed: int = 20272729,
    bootstraps: int = 500,
    matched_controls: int = 50,
    scrambles: int = 50,
) -> pd.DataFrame:
    frame = read_table(features).dropna(subset=[target, receptor_column]).reset_index(drop=True)
    seq = sequence_columns(frame)
    positions = load_interface_positions(interface_positions)
    interface = [position_column(p) for p in positions if position_column(p) in seq]
    noninterface = [c for c in seq if c not in interface]
    if not interface or not noninterface:
        raise ValueError("Both interface and non-interface sequence columns are required")

    signatures = receptor_signatures(frame, receptor_column=receptor_column, columns=seq)
    identity = sequence_identity_matrix(signatures)
    clusters = assign_identity_clusters(identity, [identity_threshold]).rename(columns={"receptor_name": receptor_column})
    cluster_columns = [c for c in clusters.columns if c != receptor_column]
    frame = frame.drop(columns=[c for c in cluster_columns if c in frame.columns], errors="ignore")
    frame = frame.merge(clusters, on=receptor_column, how="left", validate="many_to_one")
    cluster_column = f"sequence_cluster_{int(round(identity_threshold * 100))}"
    validations = [("receptor_grouped", receptor_column), (f"sequence_cluster_{int(round(identity_threshold * 100))}_grouped", cluster_column)]
    stats = _position_statistics(frame, seq, receptor_column)
    stats["is_interface"] = stats.feature.isin(interface)
    control_sets = _matched_noninterface_sets(stats, interface, noninterface, n_sets=matched_controls, seed=seed)

    output_rows: list[dict[str, object]] = []
    control_definitions: list[dict[str, object]] = []
    direct_predictions: dict[tuple[str, str], pd.DataFrame] = {}
    rng = np.random.default_rng(seed)

    for validation_name, group_column in validations:
        splits = _outer_splits(frame, target=target, group_column=group_column, folds=folds, seed=seed)
        for name, columns in (("interface_fixed_l2", interface), ("noninterface_fixed_l2", noninterface)):
            pred, metric = evaluate_logistic_representation(
                frame, RepresentationSpec(name, tuple(columns)), target=target,
                receptor_column=receptor_column, group_column=group_column,
                validation_name=validation_name, folds=folds, seed=seed,
                c_grid=(1.0,), tune_c=False, splits=splits,
            )
            metric["control_type"] = "observed"
            output_rows.append(metric)
            direct_predictions[(validation_name, name)] = pred

        for index, columns in enumerate(control_sets):
            pred, metric = evaluate_logistic_representation(
                frame, RepresentationSpec(f"matched_noninterface_{index:03d}", tuple(columns)),
                target=target, receptor_column=receptor_column, group_column=group_column,
                validation_name=validation_name, folds=folds, seed=seed,
                c_grid=(1.0,), tune_c=False, splits=splits,
            )
            metric["control_type"] = "matched_noninterface_control"
            output_rows.append(metric)
            control_definitions.append({"control_set": index, "validation": validation_name, "positions": ";".join(map(normalise_position, columns))})

        for index in range(int(scrambles)):
            shuffled = _shuffle_within_rows(frame, interface, rng)
            _, metric = evaluate_logistic_representation(
                shuffled, RepresentationSpec(f"within_receptor_position_shuffle_{index:03d}", tuple(interface)),
                target=target, receptor_column=receptor_column, group_column=group_column,
                validation_name=validation_name, folds=folds, seed=seed,
                c_grid=(1.0,), tune_c=False, splits=splits,
            )
            metric["control_type"] = "within_receptor_position_shuffle"
            output_rows.append(metric)

            shuffled = _shuffle_positions_across_rows(frame, interface, rng)
            _, metric = evaluate_logistic_representation(
                shuffled, RepresentationSpec(f"position_shuffle_across_receptors_{index:03d}", tuple(interface)),
                target=target, receptor_column=receptor_column, group_column=group_column,
                validation_name=validation_name, folds=folds, seed=seed,
                c_grid=(1.0,), tune_c=False, splits=splits,
            )
            metric["control_type"] = "position_shuffle_across_receptors"
            output_rows.append(metric)

    metrics = pd.DataFrame(output_rows)
    summary_rows = []
    for validation_name, _ in validations:
        observed = metrics[(metrics.validation == validation_name) & (metrics.experiment == "interface_fixed_l2")].iloc[0]
        for control_type in ("matched_noninterface_control", "within_receptor_position_shuffle", "position_shuffle_across_receptors"):
            null = metrics[(metrics.validation == validation_name) & (metrics.control_type == control_type)].macro_roc_auc.to_numpy(float)
            pvalue = (1 + int(np.sum(null >= float(observed.macro_roc_auc)))) / (1 + len(null))
            summary_rows.append({
                "validation": validation_name,
                "observed_model": "interface_fixed_l2",
                "observed_auc": float(observed.macro_roc_auc),
                "null_model": control_type,
                "null_iterations": int(len(null)),
                "null_mean_auc": float(np.mean(null)) if len(null) else np.nan,
                "null_ci_low": float(np.quantile(null, 0.025)) if len(null) else np.nan,
                "null_ci_high": float(np.quantile(null, 0.975)) if len(null) else np.nan,
                "empirical_p_value_one_sided": float(pvalue),
            })
    summary = pd.DataFrame(summary_rows)

    paired_rows = []
    for validation_name, _ in validations:
        result = paired_bootstrap_difference(
            direct_predictions[(validation_name, "interface_fixed_l2")],
            direct_predictions[(validation_name, "noninterface_fixed_l2")],
            iterations=bootstraps, seed=seed + 501,
        )
        result.update({"validation": validation_name, "comparison": "interface_minus_noninterface"})
        paired_rows.append(result)

    out = Path(output)
    write_table(stats, out / "audit" / "sequence_position_statistics.tsv")
    write_table(pd.DataFrame(control_definitions), out / "controls" / "matched_noninterface_position_sets.tsv")
    write_table(metrics, out / "controls" / "control_model_metrics.tsv")
    write_table(summary, out / "controls" / "interface_control_tests.tsv")
    write_table(pd.DataFrame(paired_rows), out / "controls" / "paired_interface_noninterface.tsv")
    _plot_interface_controls(summary, out / "figures")
    write_json(out / "analysis_manifest.json", {
        "analysis": "interface_specificity_audit",
        "target": target,
        "receptor_column": receptor_column,
        "identity_threshold": float(identity_threshold),
        "interface_positions": list(positions),
        "matched_controls": int(matched_controls),
        "scrambles": int(scrambles),
        "scope": "Interface localisation is descriptive unless it exceeds matched controls under low-homology transfer.",
    })
    return summary
