from pathlib import Path
import numpy as np
import pandas as pd

from gpcr_selectivitymap.sequence_audit import (
    assign_identity_clusters,
    pair_identity,
    run_cluster_benchmark,
    sequence_identity_matrix,
)


def _synthetic_frame() -> pd.DataFrame:
    labels = ["Gs", "Gi/o", "Gq/11"]
    amino = list("ACDEFGHIKLMNPQRSTVWY")
    rows = []
    for i in range(18):
        label = labels[i % 3]
        seq = [amino[(i + shift * 3) % len(amino)] for shift in range(6)]
        for rep in range(2):
            rows.append({
                "receptor_name": f"r{i:02d}",
                "transducer_family": label,
                "receptor_subfamily": f"sub{i % 6}",
                "receptor_family": f"fam{i:02d}",
                "gpcrdb_3x49_expected_aa": seq[0],
                "gpcrdb_5x61_expected_aa": seq[1],
                "gpcrdb_1x46_expected_aa": seq[2],
                "gpcrdb_2x43_expected_aa": seq[3],
                "gpcrdb_4x45_expected_aa": seq[4],
                "gpcrdb_7x51_expected_aa": seq[5],
                "geometry_signal": float(i % 3) + 0.03 * rep,
            })
    return pd.DataFrame(rows)


def test_pair_identity_and_clusters() -> None:
    assert pair_identity(("A", "C", ""), ("A", "D", "G")) == 0.5
    signatures = {"a": ("A", "C"), "b": ("A", "C"), "c": ("D", "E")}
    matrix = sequence_identity_matrix(signatures)
    clusters = assign_identity_clusters(matrix, [1.0])
    assert clusters.loc[clusters.receptor_name.eq("a"), "sequence_cluster_100"].iloc[0] == clusters.loc[
        clusters.receptor_name.eq("b"), "sequence_cluster_100"
    ].iloc[0]


def test_cluster_benchmark_smoke(tmp_path: Path) -> None:
    table = tmp_path / "features.tsv"
    _synthetic_frame().to_csv(table, sep="\t", index=False)
    metrics = run_cluster_benchmark(
        features=table,
        target="transducer_family",
        receptor_column="receptor_name",
        output=tmp_path / "out",
        thresholds=(1.0,),
        folds=3,
        bootstraps=10,
        models=("full_sequence", "interface_sequence", "noninterface_sequence", "endpoint_geometry", "global_identity_knn"),
        seed=11,
    )
    assert not metrics.empty
    assert (tmp_path / "out/metrics/model_metrics.tsv").exists()
    assert (tmp_path / "out/predictions/out_of_fold_predictions.tsv").exists()
    assert (tmp_path / "out/analysis_manifest.json").exists()
