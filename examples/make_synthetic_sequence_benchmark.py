#!/usr/bin/env python3
"""Create a clearly synthetic table for testing the phylogeny-aware CLI."""
from pathlib import Path
import pandas as pd


def main() -> None:
    labels = ["Gs", "Gi/o", "Gq/11"]
    amino = list("ACDEFGHIKLMNPQRSTVWY")
    rows = []
    for i in range(24):
        sequence = [amino[(i + 3 * j) % len(amino)] for j in range(6)]
        rows.append({
            "receptor_name": f"synthetic_r{i:02d}",
            "transducer_family": labels[i % 3],
            "receptor_subfamily": f"synthetic_subfamily_{i % 6}",
            "receptor_family": f"synthetic_family_{i:02d}",
            "gpcrdb_3x49_expected_aa": sequence[0],
            "gpcrdb_5x61_expected_aa": sequence[1],
            "gpcrdb_1x46_expected_aa": sequence[2],
            "gpcrdb_2x43_expected_aa": sequence[3],
            "gpcrdb_4x45_expected_aa": sequence[4],
            "gpcrdb_7x51_expected_aa": sequence[5],
            "geometry_signal": float(i % 3),
        })
    path = Path("synthetic_sequence_benchmark.tsv")
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    print(path.resolve())


if __name__ == "__main__":
    main()
