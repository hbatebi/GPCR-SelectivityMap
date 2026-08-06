from __future__ import annotations

"""Label provenance validation for bound-complex and functional annotations."""

from pathlib import Path
import pandas as pd
from .io import read_table, write_json, write_table

PROVENANCE_COLUMNS = (
    "label_source",
    "label_type",
    "experimental_support",
    "curation_status",
    "homology_inferred",
    "reference",
)


def run_provenance_audit(
    *,
    labels: str | Path,
    output: str | Path,
    receptor_column: str = "receptor_name",
    target: str = "transducer_family",
) -> pd.DataFrame:
    frame = read_table(labels)
    missing_required = [c for c in (receptor_column, target) if c not in frame.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {', '.join(missing_required)}")
    missing_provenance = [c for c in PROVENANCE_COLUMNS if c not in frame.columns]
    summary_rows = [{
        "metric": "n_rows",
        "value": int(len(frame)),
        "detail": "All label rows",
    }, {
        "metric": "n_receptors",
        "value": int(frame[receptor_column].nunique()),
        "detail": "Independent receptor identifiers",
    }, {
        "metric": "missing_provenance_columns",
        "value": int(len(missing_provenance)),
        "detail": ";".join(missing_provenance),
    }]
    for column in PROVENANCE_COLUMNS:
        if column not in frame.columns:
            continue
        values = frame[column].fillna("unknown").astype(str).str.strip().replace({"": "unknown"})
        for value, count in values.value_counts(dropna=False).items():
            summary_rows.append({
                "metric": f"{column}_count",
                "value": int(count),
                "detail": str(value),
            })
    if "homology_inferred" in frame.columns:
        normalized = frame.homology_inferred.fillna("unknown").astype(str).str.lower().str.strip()
        potentially_inferred = normalized.isin({"true", "1", "yes", "y", "homology", "inferred"})
        summary_rows.append({
            "metric": "potentially_homology_inferred_rows",
            "value": int(potentially_inferred.sum()),
            "detail": "Rows explicitly marked as inferred from homology",
        })
    summary = pd.DataFrame(summary_rows)
    out = Path(output)
    write_table(summary, out / "label_provenance_summary.tsv")
    write_table(frame, out / "label_provenance_audited.tsv")
    write_json(out / "label_provenance_manifest.json", {
        "labels": str(labels),
        "receptor_column": receptor_column,
        "target": target,
        "missing_provenance_columns": missing_provenance,
        "complete": not missing_provenance,
        "interpretation": (
            "Bound-complex labels may be direct structural observations. Functional labels should not be described as directly measured "
            "unless provenance fields support that statement."
        ),
    })
    return summary
