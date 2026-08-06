from __future__ import annotations

import numpy as np
import pandas as pd

from gpcr_selectivitymap.structure import ResidueKey, ResidueRecord
from gpcr_selectivitymap.features import (
    build_receptor_frame,
    contact_features,
    mechanical_features,
    normalise_generic,
)


def residue(number: int, coord: tuple[float, float, float], name: str = "ALA") -> ResidueRecord:
    ca = np.asarray(coord, dtype=float)
    atoms = {"CA": ca, "CB": ca + np.array([0.5, 0.2, 0.1])}
    return ResidueRecord(
        chain_id="R",
        key=ResidueKey(number, ""),
        resname=name,
        atoms=atoms,
        occupancies={key: 1.0 for key in atoms},
        altlocs={key: "" for key in atoms},
    )


def synthetic_generic() -> dict[str, ResidueRecord]:
    positions = [
        "2.39", "2.43", "3.39", "3.40", "3.43", "3.46", "3.49", "3.50", "3.53",
        "34.50", "34.51", "34.52", "34.53", "34.54", "34.55",
        "4.38", "4.45", "4.50", "5.42", "5.46", "5.50", "5.58", "5.61", "5.64",
        "6.23", "6.30", "6.34", "6.40", "6.44", "6.48", "7.45", "7.50", "7.53", "7.55",
        "8.47", "8.50", "8.53",
    ]
    result = {}
    for index, position in enumerate(positions, start=1):
        major = float(position.split(".")[0])
        minor = float(position.split(".")[1])
        coord = (
            np.cos(index * 0.43) * (8.0 + major * 0.2),
            np.sin(index * 0.43) * (8.0 + major * 0.2),
            (minor - 45.0) * 0.35 + major * 0.5,
        )
        result[position] = residue(index, coord, "LYS" if index % 11 == 0 else "GLU" if index % 13 == 0 else "ALA")
    # Force a well-defined receptor frame.
    result["3.50"] = residue(100, (-5.0, 0.0, 0.0))
    result["5.58"] = residue(101, (0.0, 5.0, 0.0))
    result["6.34"] = residue(102, (5.0, 0.0, 0.0))
    result["7.53"] = residue(103, (0.0, -5.0, 0.0))
    result["3.40"] = residue(104, (-5.0, 0.0, 12.0))
    result["5.46"] = residue(105, (0.0, 5.0, 12.0))
    result["6.44"] = residue(106, (5.0, 0.0, 12.0))
    result["7.45"] = residue(107, (0.0, -5.0, 12.0))
    return result


def test_generic_number_normalisation() -> None:
    assert normalise_generic("3x50") == "3.50"
    assert normalise_generic(34.5) == "34.50"


def test_receptor_frame_is_orthonormal() -> None:
    origin, x_axis, y_axis, z_axis = build_receptor_frame(synthetic_generic())
    assert np.isfinite(origin).all()
    matrix = np.stack([x_axis, y_axis, z_axis])
    assert np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-6)


def test_contact_and_mechanical_features_are_finite() -> None:
    generic = synthetic_generic()
    contacts = contact_features(generic)
    mechanics = mechanical_features(generic, cutoff=18.0, modes_requested=6)
    assert contacts["contact_network_nodes"] >= 20
    assert np.isfinite(contacts["contact_network_density"])
    assert mechanics["mechanical_mode_count"] >= 3
    assert np.isfinite(mechanics["mechanical_total_softness"])
