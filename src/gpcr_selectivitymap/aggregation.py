from __future__ import annotations

from pathlib import Path
from typing import Sequence
import pandas as pd
from .io import read_table, write_json, write_table


def _mode_or_first(series: pd.Series):
    nonmissing = series.dropna()
    if nonmissing.empty:
        return pd.NA
    mode = nonmissing.mode()
    return mode.iloc[0] if len(mode) else nonmissing.iloc[0]


def aggregate_feature_table(
    *,
    features: str | Path,
    units: Sequence[str],
    output: str | Path,
    numeric_method: str = "median",
) -> pd.DataFrame:
    frame = read_table(features)
    missing = [c for c in units if c not in frame.columns]
    if missing:
        raise ValueError(f"Missing aggregation units: {', '.join(missing)}")
    numeric = [c for c in frame.select_dtypes(include="number").columns if c not in units]
    categorical = [c for c in frame.columns if c not in set(units) | set(numeric)]
    if numeric_method not in {"median", "mean"}:
        raise ValueError("numeric_method must be median or mean")
    agg: dict[str, object] = {c: numeric_method for c in numeric}
    agg.update({c: _mode_or_first for c in categorical})
    result = frame.groupby(list(units), dropna=False, as_index=False).agg(agg)
    out = Path(output)
    if out.suffix.lower() not in {".csv", ".tsv", ".tab"}:
        out.mkdir(parents=True, exist_ok=True)
        table_path = out / "aggregated_features.tsv"
        manifest_path = out / "aggregation_manifest.json"
    else:
        table_path = out
        manifest_path = out.with_suffix(".manifest.json")
    write_table(result, table_path)
    write_json(manifest_path, {
        "input_rows": int(len(frame)),
        "output_rows": int(len(result)),
        "units": list(units),
        "numeric_method": numeric_method,
        "categorical_method": "mode_or_first",
    })
    return result
