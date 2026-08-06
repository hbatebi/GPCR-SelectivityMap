from pathlib import Path
import pandas as pd

from gpcr_selectivitymap.interface_audit import run_interface_audit
from gpcr_selectivitymap.provenance import run_provenance_audit


def _frame() -> pd.DataFrame:
    labels = ["Gs", "Gi/o", "Gq/11"]
    amino = list("ACDEFGHIKLMNPQRSTVWY")
    rows = []
    for i in range(18):
        rows.append({
            "receptor_name": f"r{i}",
            "transducer_family": labels[i % 3],
            "gpcrdb_3x49_expected_aa": amino[i % 20],
            "gpcrdb_5x61_expected_aa": amino[(i + 2) % 20],
            "gpcrdb_1x46_expected_aa": amino[(i + 4) % 20],
            "gpcrdb_2x43_expected_aa": amino[(i + 6) % 20],
            "gpcrdb_4x45_expected_aa": amino[(i + 8) % 20],
        })
    return pd.DataFrame(rows)


def test_interface_audit_smoke(tmp_path: Path) -> None:
    table = tmp_path / "features.tsv"
    _frame().to_csv(table, sep="\t", index=False)
    summary = run_interface_audit(
        features=table,
        target="transducer_family",
        receptor_column="receptor_name",
        output=tmp_path / "audit",
        identity_threshold=1.0,
        folds=3,
        bootstraps=10,
        matched_controls=2,
        scrambles=2,
        seed=19,
    )
    assert len(summary) == 6
    assert (tmp_path / "audit/controls/interface_control_tests.tsv").exists()


def test_provenance_audit_flags_missing_columns(tmp_path: Path) -> None:
    labels = tmp_path / "labels.tsv"
    pd.DataFrame({
        "receptor_name": ["r1", "r2"],
        "transducer_family": ["Gs", "Gi/o"],
        "label_source": ["paper", "database"],
    }).to_csv(labels, sep="\t", index=False)
    result = run_provenance_audit(labels=labels, output=tmp_path / "prov")
    missing = result[result.metric.eq("missing_provenance_columns")].iloc[0]
    assert int(missing.value) > 0
    assert (tmp_path / "prov/label_provenance_manifest.json").exists()
