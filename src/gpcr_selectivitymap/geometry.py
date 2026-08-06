from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence

import numpy as np

BACKBONE = {"N", "CA", "C", "O", "OXT"}
CHI_ATOMS: dict[str, tuple[tuple[str, str, str, str], ...]] = {
    "ARG": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")),
    "ASN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")),
    "ASP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")),
    "CYS": (("N", "CA", "CB", "SG"),),
    "GLN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")),
    "GLU": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")),
    "HIS": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")),
    "ILE": (("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")),
    "LEU": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "LYS": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")),
    "MET": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD")),
    "PHE": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "PRO": (("N", "CA", "CB", "CG"),),
    "SER": (("N", "CA", "CB", "OG"),),
    "THR": (("N", "CA", "CB", "OG1"),),
    "TRP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "TYR": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "VAL": (("N", "CA", "CB", "CG1"),),
}
CHARGED_ATOMS = {
    "ARG": {"NH1", "NH2", "NE"}, "LYS": {"NZ"}, "HIS": {"ND1", "NE2"},
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
}
AROMATIC_RING_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}


def as_coord(value: object) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise ValueError(f"Expected a finite 3-vector, got shape={arr.shape}")
    return arr


def distance(a: object, b: object) -> float:
    return float(np.linalg.norm(as_coord(a) - as_coord(b)))


def pairwise_min_distance(atoms_a: Mapping[str, np.ndarray], atoms_b: Mapping[str, np.ndarray],
                          names_a: Iterable[str] | None = None,
                          names_b: Iterable[str] | None = None) -> float | None:
    keys_a = [k for k in (names_a or atoms_a.keys()) if k in atoms_a and not k.startswith("H")]
    keys_b = [k for k in (names_b or atoms_b.keys()) if k in atoms_b and not k.startswith("H")]
    if not keys_a or not keys_b:
        return None
    aa = np.stack([atoms_a[k] for k in keys_a])
    bb = np.stack([atoms_b[k] for k in keys_b])
    d = np.linalg.norm(aa[:, None, :] - bb[None, :, :], axis=2)
    return float(np.nanmin(d))


def centroid(atoms: Mapping[str, np.ndarray], names: Iterable[str] | None = None) -> np.ndarray | None:
    keys = [k for k in (names or atoms.keys()) if k in atoms and not k.startswith("H")]
    if not keys:
        return None
    return np.mean(np.stack([atoms[k] for k in keys]), axis=0)


def sidechain_centroid(resname: str, atoms: Mapping[str, np.ndarray]) -> np.ndarray | None:
    names = [k for k in atoms if k not in BACKBONE and not k.startswith("H")]
    if not names:
        return atoms.get("CA")
    return centroid(atoms, names)


def aromatic_centroid(resname: str, atoms: Mapping[str, np.ndarray]) -> np.ndarray | None:
    ring = AROMATIC_RING_ATOMS.get(resname.upper())
    if ring:
        out = centroid(atoms, ring)
        if out is not None:
            return out
    return sidechain_centroid(resname, atoms)


def dihedral(p0: object, p1: object, p2: object, p3: object) -> float:
    p0, p1, p2, p3 = map(as_coord, (p0, p1, p2, p3))
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    norm = np.linalg.norm(b1)
    if norm == 0:
        return float("nan")
    b1 /= norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) == 0 or np.linalg.norm(w) == 0:
        return float("nan")
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def chi_angles(resname: str, atoms: Mapping[str, np.ndarray]) -> tuple[float | None, float | None]:
    defs = CHI_ATOMS.get(resname.upper(), ())
    values: list[float | None] = []
    for names in defs[:2]:
        if all(name in atoms for name in names):
            val = dihedral(*(atoms[name] for name in names))
            values.append(None if not np.isfinite(val) else val)
        else:
            values.append(None)
    while len(values) < 2:
        values.append(None)
    return values[0], values[1]


def angle_between(v1: object, v2: object) -> float | None:
    a = as_coord(v1)
    b = as_coord(v2)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return None
    cosine = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return rotation, translation and RMSD mapping mobile onto target."""
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    if mobile.shape != target.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("mobile and target must be matching (n, 3) arrays")
    if len(mobile) < 3:
        raise ValueError("At least three points are required for alignment")
    mob_center = mobile.mean(axis=0)
    tar_center = target.mean(axis=0)
    x = mobile - mob_center
    y = target - tar_center
    u, _s, vt = np.linalg.svd(x.T @ y)
    det = np.linalg.det(u @ vt)
    correction = np.eye(3)
    correction[-1, -1] = -1.0 if det < 0 else 1.0
    rotation = u @ correction @ vt
    translation = tar_center - mob_center @ rotation
    fitted = mobile @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def transform(coords: np.ndarray | Sequence[float], rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    arr = np.asarray(coords, dtype=float)
    return arr @ rotation + translation


def circular_difference(a: float | None, b: float | None) -> float | None:
    if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
        return None
    return float(abs((a - b + 180.0) % 360.0 - 180.0))


def pocket_volume(points: Sequence[np.ndarray]) -> float | None:
    pts = np.asarray(points, dtype=float)
    if len(pts) < 4 or np.linalg.matrix_rank(pts - pts.mean(axis=0)) < 3:
        return None
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(pts).volume)
    except Exception:
        return None


def rms_spread(points: Sequence[np.ndarray]) -> float | None:
    pts = np.asarray(points, dtype=float)
    if len(pts) < 2:
        return None
    center = pts.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((pts - center) ** 2, axis=1))))


def safe_round(value: float | int | None, digits: int = 4):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if not math.isfinite(value) else round(value, digits)
