from __future__ import annotations

from collections import deque
from functools import lru_cache
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from scipy.linalg import eigh
from scipy.spatial import cKDTree, distance_matrix

from .structure import ChainRecord, ResidueKey, ResidueRecord, load_chains

# The panel is intentionally restricted to generic positions already mapped by
# the parent pipeline. No atom from a G protein or other non-receptor chain is
# ever read into the feature matrix.
INTRACELLULAR_POSITIONS = (
    "2.39", "2.43", "2.46", "2.47", "2.50", "2.53",
    "3.43", "3.46", "3.49", "3.50", "3.53", "3.55",
    "34.50", "34.51", "34.52", "34.53", "34.54", "34.55", "34.56", "34.57",
    "4.38", "4.45", "4.46", "4.50", "4.53",
    "5.50", "5.53", "5.58", "5.61", "5.64", "5.68", "5.72",
    "6.23", "6.30", "6.34", "6.37", "6.40", "6.44", "6.47", "6.48", "6.50",
    "7.45", "7.46", "7.49", "7.50", "7.51", "7.52", "7.53", "7.54", "7.55", "7.56",
    "8.47", "8.50", "8.53", "8.59",
)

ORTHOSTERIC_POSITIONS = ("3.39", "5.42", "5.46", "5.50", "6.44", "6.48", "7.45", "7.46")
INTERFACE_POSITIONS = ("3.50", "34.51", "34.54", "5.61", "5.64", "6.30", "6.34", "7.53", "8.50")
FRAME_INNER = ("3.50", "5.58", "6.34", "7.53")
FRAME_OUTER = ("3.40", "5.46", "6.44", "7.45")

REGION_PREFIXES: Mapping[str, tuple[str, ...]] = {
    "tm3": ("3.",),
    "icl2": ("34.",),
    "tm5": ("5.",),
    "tm6": ("6.",),
    "tm7": ("7.",),
    "h8": ("8.",),
}

POSITIVE = {"ARG", "LYS"}
NEGATIVE = {"ASP", "GLU"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
CHARGE = {"ARG": 1.0, "LYS": 1.0, "HIS": 0.1, "ASP": -1.0, "GLU": -1.0}
DONOR_RESIDUES = {"ARG", "LYS", "HIS", "SER", "THR", "ASN", "GLN", "TYR", "TRP", "CYS"}
ACCEPTOR_RESIDUES = {"ASP", "GLU", "HIS", "SER", "THR", "ASN", "GLN", "TYR", "CYS"}


@dataclass(frozen=True)
class StructureTask:
    pdb_id: str
    receptor_name: str
    chain_id: str
    dataset_role: str
    transducer_family: str = ""


def normalise_generic(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip().lower().replace("x", ".")
    try:
        number = float(text)
    except ValueError:
        return text
    major = int(np.floor(number + 1e-9))
    minor = int(round((number - major) * 100))
    return f"{major}.{minor:02d}"


def _unit(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm < 1e-8:
        raise ValueError("Cannot normalise a zero or non-finite vector")
    return vector / norm


def _centroid(values: Iterable[np.ndarray]) -> np.ndarray | None:
    arrays = [np.asarray(value, dtype=float) for value in values if value is not None]
    if not arrays:
        return None
    return np.mean(np.stack(arrays), axis=0)


def build_receptor_frame(generic_residues: Mapping[str, ResidueRecord | None]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    inner = _centroid(
        generic_residues[position].ca
        for position in FRAME_INNER
        if position in generic_residues and generic_residues[position] is not None
    )
    outer = _centroid(
        generic_residues[position].ca
        for position in FRAME_OUTER
        if position in generic_residues and generic_residues[position] is not None
    )
    if inner is None or outer is None:
        raise ValueError("Insufficient generic positions to define the receptor-centred frame")
    z_axis = _unit(outer - inner)

    tm3 = generic_residues.get("3.50")
    tm6 = generic_residues.get("6.34") or generic_residues.get("6.30")
    if tm3 is None or tm6 is None or tm3.ca is None or tm6.ca is None:
        candidates = [
            residue.ca
            for position, residue in generic_residues.items()
            if residue is not None and residue.ca is not None and position in INTRACELLULAR_POSITIONS
        ]
        if len(candidates) < 3:
            raise ValueError("Insufficient intracellular coordinates for frame x axis")
        matrix = np.stack(candidates) - inner
        _, _, vh = np.linalg.svd(matrix, full_matrices=False)
        raw_x = vh[0]
    else:
        raw_x = tm6.ca - tm3.ca
    raw_x = raw_x - z_axis * float(np.dot(raw_x, z_axis))
    x_axis = _unit(raw_x)
    y_axis = _unit(np.cross(z_axis, x_axis))
    # Re-orthogonalise x to guard against numerical drift.
    x_axis = _unit(np.cross(y_axis, z_axis))
    return inner, x_axis, y_axis, z_axis


def _mapping_for_task(mapping: pd.DataFrame, task: StructureTask) -> dict[str, ResidueKey]:
    part = mapping.loc[mapping["pdb_id"].astype(str).str.upper().eq(task.pdb_id)].copy()
    if task.chain_id:
        chain_part = part.loc[part["chain_id"].astype(str).eq(task.chain_id)]
        if not chain_part.empty:
            part = chain_part
    result: dict[str, ResidueKey] = {}
    for row in part.itertuples(index=False):
        position = normalise_generic(getattr(row, "generic_number", ""))
        number = getattr(row, "author_residue_number", np.nan)
        if not position or pd.isna(number):
            continue
        raw_insertion = getattr(row, "author_insertion_code", "")
        insertion = "" if pd.isna(raw_insertion) else str(raw_insertion).strip()
        result[position] = ResidueKey(int(float(number)), insertion)
    return result


def _select_chain(chains: Mapping[str, ChainRecord], requested: str) -> ChainRecord:
    if requested and requested in chains:
        return chains[requested]
    if len(chains) == 1:
        return next(iter(chains.values()))
    # The mapping-derived chain is preferable; if unavailable, choose the chain
    # with the largest amino-acid content and record the fallback in quality.
    return max(chains.values(), key=lambda chain: len(chain.residues))


def _mapped_residues(chain: ChainRecord, keys: Mapping[str, ResidueKey]) -> dict[str, ResidueRecord | None]:
    by_key = chain.by_key()
    return {position: by_key.get(key) for position, key in keys.items()}


def _representative_coordinate(residue: ResidueRecord) -> np.ndarray | None:
    return residue.sidechain_center if residue.sidechain_center is not None else residue.ca


def _minimum_atom_distance(left: ResidueRecord, right: ResidueRecord, left_atoms: Sequence[str] | None = None, right_atoms: Sequence[str] | None = None) -> float:
    lnames = left_atoms if left_atoms is not None else tuple(left.atoms)
    rnames = right_atoms if right_atoms is not None else tuple(right.atoms)
    lcoords = [left.atoms[name] for name in lnames if name in left.atoms and not name.startswith("H")]
    rcoords = [right.atoms[name] for name in rnames if name in right.atoms and not name.startswith("H")]
    if not lcoords or not rcoords:
        return float("nan")
    return float(np.min(distance_matrix(np.stack(lcoords), np.stack(rcoords))))


def _sidechain_hetero_atoms(residue: ResidueRecord) -> list[np.ndarray]:
    coords = []
    for name, coord in residue.atoms.items():
        if name in {"N", "CA", "C", "O", "OXT"} or name.startswith("H"):
            continue
        if name.startswith(("N", "O", "S")):
            coords.append(coord)
    return coords


def _graph_metrics(adjacency: np.ndarray, positions: Sequence[str]) -> dict[str, float]:
    n = int(adjacency.shape[0])
    if n == 0:
        return {"contact_network_nodes": 0.0}
    adjacency = (adjacency > 0).astype(int)
    np.fill_diagonal(adjacency, 0)
    degrees = adjacency.sum(axis=1).astype(float)
    edges = int(adjacency.sum() // 2)
    density = float(2 * edges / (n * (n - 1))) if n > 1 else 0.0

    clustering = []
    for index in range(n):
        neighbours = np.flatnonzero(adjacency[index])
        k = len(neighbours)
        if k < 2:
            clustering.append(0.0)
            continue
        links = int(adjacency[np.ix_(neighbours, neighbours)].sum() // 2)
        clustering.append(2.0 * links / (k * (k - 1)))

    # Brandes betweenness for the small generic-position graph.
    betweenness = np.zeros(n, dtype=float)
    for source in range(n):
        stack: list[int] = []
        predecessors: list[list[int]] = [[] for _ in range(n)]
        sigma = np.zeros(n, dtype=float)
        sigma[source] = 1.0
        distance = np.full(n, -1, dtype=int)
        distance[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            vertex = queue.popleft()
            stack.append(vertex)
            for neighbour in np.flatnonzero(adjacency[vertex]):
                if distance[neighbour] < 0:
                    queue.append(int(neighbour))
                    distance[neighbour] = distance[vertex] + 1
                if distance[neighbour] == distance[vertex] + 1:
                    sigma[neighbour] += sigma[vertex]
                    predecessors[neighbour].append(vertex)
        dependency = np.zeros(n, dtype=float)
        while stack:
            vertex = stack.pop()
            if sigma[vertex] > 0:
                factor = (1.0 + dependency[vertex]) / sigma[vertex]
                for predecessor in predecessors[vertex]:
                    dependency[predecessor] += sigma[predecessor] * factor
            if vertex != source:
                betweenness[vertex] += dependency[vertex]
    if n > 2:
        betweenness /= ((n - 1) * (n - 2))

    # Floyd-Warshall is inexpensive for ~50 generic-position nodes.
    distances = np.full((n, n), np.inf, dtype=float)
    distances[adjacency > 0] = 1.0
    np.fill_diagonal(distances, 0.0)
    for k in range(n):
        distances = np.minimum(distances, distances[:, [k]] + distances[[k], :])
    finite = distances[np.isfinite(distances) & (distances > 0)]

    laplacian = np.diag(degrees) - adjacency
    eigenvalues = np.linalg.eigvalsh(laplacian) if n > 1 else np.array([0.0])
    algebraic = float(eigenvalues[1]) if len(eigenvalues) > 1 else 0.0

    position_to_index = {position: index for index, position in enumerate(positions)}
    source_indices = [position_to_index[p] for p in ORTHOSTERIC_POSITIONS if p in position_to_index]
    target_indices = [position_to_index[p] for p in INTERFACE_POSITIONS if p in position_to_index]
    path_values = [distances[i, j] for i in source_indices for j in target_indices if np.isfinite(distances[i, j])]

    return {
        "contact_network_nodes": float(n),
        "contact_network_edges": float(edges),
        "contact_network_density": density,
        "contact_network_degree_mean": float(np.mean(degrees)),
        "contact_network_degree_max": float(np.max(degrees)),
        "contact_network_clustering_mean": float(np.mean(clustering)),
        "contact_network_betweenness_mean": float(np.mean(betweenness)),
        "contact_network_betweenness_max": float(np.max(betweenness)),
        "contact_network_algebraic_connectivity": algebraic,
        "contact_network_average_path_length": float(np.mean(finite)) if finite.size else float("nan"),
        "contact_network_orthosteric_interface_path_min": float(np.min(path_values)) if path_values else float("nan"),
        "contact_network_orthosteric_interface_path_mean": float(np.mean(path_values)) if path_values else float("nan"),
    }


def contact_features(generic_residues: Mapping[str, ResidueRecord | None], cutoff: float = 8.0) -> dict[str, float]:
    nodes = [
        (position, residue, _representative_coordinate(residue))
        for position, residue in sorted(generic_residues.items())
        if residue is not None and _representative_coordinate(residue) is not None
    ]
    positions = [item[0] for item in nodes]
    residues = [item[1] for item in nodes]
    coords = np.stack([item[2] for item in nodes]) if nodes else np.empty((0, 3))
    adjacency = np.zeros((len(nodes), len(nodes)), dtype=int)
    if len(nodes) > 1:
        adjacency = (distance_matrix(coords, coords) <= cutoff).astype(int)
        np.fill_diagonal(adjacency, 0)
    features = _graph_metrics(adjacency, positions)

    salt_bridges = 0
    aromatic_contacts = 0
    cation_pi_contacts = 0
    polar_close_contacts = 0
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            left, right = residues[i], residues[j]
            if left.resname in POSITIVE and right.resname in NEGATIVE or left.resname in NEGATIVE and right.resname in POSITIVE:
                charged_atoms = {
                    "ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",),
                    "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
                }
                if _minimum_atom_distance(left, right, charged_atoms.get(left.resname), charged_atoms.get(right.resname)) <= 4.0:
                    salt_bridges += 1
            center_distance = float(np.linalg.norm(coords[i] - coords[j]))
            if left.resname in AROMATIC and right.resname in AROMATIC and center_distance <= 6.0:
                aromatic_contacts += 1
            if ((left.resname in POSITIVE and right.resname in AROMATIC) or (right.resname in POSITIVE and left.resname in AROMATIC)) and center_distance <= 6.0:
                cation_pi_contacts += 1
            left_hetero = _sidechain_hetero_atoms(left)
            right_hetero = _sidechain_hetero_atoms(right)
            if left_hetero and right_hetero and float(np.min(distance_matrix(np.stack(left_hetero), np.stack(right_hetero)))) <= 3.5:
                polar_close_contacts += 1

    features.update(
        {
            "contact_salt_bridge_count": float(salt_bridges),
            "contact_aromatic_pair_count": float(aromatic_contacts),
            "contact_cation_pi_count": float(cation_pi_contacts),
            "contact_polar_close_count": float(polar_close_contacts),
        }
    )

    def min_between(left_prefixes: Sequence[str], right_prefixes: Sequence[str]) -> float:
        left_indices = [i for i, p in enumerate(positions) if p.startswith(tuple(left_prefixes))]
        right_indices = [i for i, p in enumerate(positions) if p.startswith(tuple(right_prefixes))]
        values = [float(np.linalg.norm(coords[i] - coords[j])) for i in left_indices for j in right_indices if i != j]
        return float(np.min(values)) if values else float("nan")

    features.update(
        {
            "contact_icl2_tm3_min_distance": min_between(("34.",), ("3.",)),
            "contact_icl2_tm5_min_distance": min_between(("34.",), ("5.",)),
            "contact_icl2_tm6_min_distance": min_between(("34.",), ("6.",)),
            "contact_tm3_tm6_min_distance": min_between(("3.",), ("6.",)),
            "contact_tm7_h8_min_distance": min_between(("7.",), ("8.",)),
        }
    )
    return features


def _bio_chain_sasa(structure_path: Path, chain_id: str, n_points: int) -> dict[ResidueKey, float]:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(structure_path.stem, str(structure_path))
    model = next(structure.get_models())
    if chain_id not in model:
        chains = list(model.get_chains())
        if len(chains) == 1:
            chain = chains[0]
        else:
            chain = max(chains, key=lambda item: sum(1 for _ in item.get_residues()))
    else:
        chain = model[chain_id]
    ShrakeRupley(n_points=n_points).compute(chain, level="R")
    result: dict[ResidueKey, float] = {}
    for residue in chain:
        _hetflag, number, insertion = residue.id
        if hasattr(residue, "sasa"):
            result[ResidueKey(int(number), str(insertion).strip())] = float(residue.sasa)
    return result


def surface_features(
    structure_path: Path,
    chain: ChainRecord,
    generic_residues: Mapping[str, ResidueRecord | None],
    n_points: int = 80,
) -> dict[str, float]:
    sasa_by_key = _bio_chain_sasa(structure_path, chain.chain_id, n_points)
    rows = []
    for position in INTRACELLULAR_POSITIONS:
        residue = generic_residues.get(position)
        if residue is None:
            continue
        sasa = sasa_by_key.get(residue.key, np.nan)
        rows.append((position, residue, float(sasa)))
    available = [(p, r, s) for p, r, s in rows if np.isfinite(s)]
    total = float(sum(s for _, _, s in available))
    features: dict[str, float] = {
        "surface_residue_count": float(len(available)),
        "surface_total_sasa": total,
        "surface_mean_sasa": float(np.mean([s for _, _, s in available])) if available else float("nan"),
    }
    property_sets = {
        "positive": POSITIVE,
        "negative": NEGATIVE,
        "aromatic": AROMATIC,
        "hydrophobic": HYDROPHOBIC,
        "polar": POLAR,
        "donor": DONOR_RESIDUES,
        "acceptor": ACCEPTOR_RESIDUES,
    }
    for name, residues in property_sets.items():
        area = float(sum(s for _, residue, s in available if residue.resname in residues))
        features[f"surface_{name}_sasa"] = area
        features[f"surface_{name}_fraction"] = area / total if total > 0 else float("nan")

    # Patch continuity proxy: among exposed residues of each physicochemical
    # class, report the fraction in the largest <=8 A connected component.
    for name, residueset in {
        "positive": POSITIVE,
        "negative": NEGATIVE,
        "hydrophobic": HYDROPHOBIC,
        "polar": POLAR,
    }.items():
        selected = [
            (residue, _representative_coordinate(residue))
            for _, residue, sasa in available
            if residue.resname in residueset and sasa >= 10.0 and _representative_coordinate(residue) is not None
        ]
        if not selected:
            features[f"surface_{name}_largest_patch_fraction"] = float("nan")
            continue
        coordinates = np.stack([coord for _, coord in selected])
        graph = (distance_matrix(coordinates, coordinates) <= 8.0).astype(int)
        np.fill_diagonal(graph, 0)
        seen: set[int] = set()
        sizes = []
        for start in range(len(selected)):
            if start in seen:
                continue
            queue = [start]
            seen.add(start)
            count = 0
            while queue:
                node = queue.pop()
                count += 1
                for neighbour in np.flatnonzero(graph[node]):
                    neighbour = int(neighbour)
                    if neighbour not in seen:
                        seen.add(neighbour)
                        queue.append(neighbour)
            sizes.append(count)
        features[f"surface_{name}_largest_patch_fraction"] = max(sizes) / len(selected)
    return features


def electrostatic_features(chain: ChainRecord, generic_residues: Mapping[str, ResidueRecord | None]) -> dict[str, float]:
    origin, x_axis, y_axis, z_axis = build_receptor_frame(generic_residues)
    charges = []
    for residue in chain.residues:
        charge = CHARGE.get(residue.resname, 0.0)
        coordinate = _representative_coordinate(residue)
        if charge and coordinate is not None and np.linalg.norm(coordinate - origin) <= 35.0:
            charges.append((charge, coordinate))

    sample_points = {
        "z0": origin,
        "z2": origin + 2.0 * z_axis,
        "z4": origin + 4.0 * z_axis,
        "z6": origin + 6.0 * z_axis,
        "xp": origin + 2.0 * z_axis + 4.0 * x_axis,
        "xm": origin + 2.0 * z_axis - 4.0 * x_axis,
        "yp": origin + 2.0 * z_axis + 4.0 * y_axis,
        "ym": origin + 2.0 * z_axis - 4.0 * y_axis,
    }

    def potential(point: np.ndarray) -> float:
        return float(sum(charge / max(float(np.linalg.norm(point - coordinate)), 1.0) for charge, coordinate in charges))

    potentials = {name: potential(point) for name, point in sample_points.items()}
    field = np.zeros(3, dtype=float)
    for charge, coordinate in charges:
        delta = origin - coordinate
        distance = max(float(np.linalg.norm(delta)), 1.0)
        field += charge * delta / distance**3

    features: dict[str, float] = {
        "electro_charge_site_count": float(len(charges)),
        "electro_potential_mean": float(np.mean(list(potentials.values()))),
        "electro_potential_std": float(np.std(list(potentials.values()))),
        "electro_potential_min": float(np.min(list(potentials.values()))),
        "electro_potential_max": float(np.max(list(potentials.values()))),
        "electro_potential_z_gradient": potentials["z6"] - potentials["z0"],
        "electro_potential_x_asymmetry": potentials["xp"] - potentials["xm"],
        "electro_potential_y_asymmetry": potentials["yp"] - potentials["ym"],
        "electro_field_magnitude": float(np.linalg.norm(field)),
        "electro_field_z_projection": float(np.dot(field, z_axis)),
        "electro_field_x_projection": float(np.dot(field, x_axis)),
        "electro_field_y_projection": float(np.dot(field, y_axis)),
    }
    for radius in (10.0, 15.0, 20.0):
        features[f"electro_net_charge_within_{int(radius)}A"] = float(
            sum(charge for charge, coordinate in charges if np.linalg.norm(coordinate - origin) <= radius)
        )
    for name, value in potentials.items():
        features[f"electro_potential_{name}"] = value
    return features


def _anm_hessian(coords: np.ndarray, cutoff: float) -> np.ndarray:
    n = len(coords)
    hessian = np.zeros((3 * n, 3 * n), dtype=float)
    tree = cKDTree(coords)
    for i, j in tree.query_pairs(cutoff):
        delta = coords[j] - coords[i]
        distance2 = float(np.dot(delta, delta))
        if distance2 < 1e-8:
            continue
        block = np.outer(delta, delta) / distance2
        si = slice(3 * i, 3 * i + 3)
        sj = slice(3 * j, 3 * j + 3)
        hessian[si, si] += block
        hessian[sj, sj] += block
        hessian[si, sj] -= block
        hessian[sj, si] -= block
    return hessian


def mechanical_features(
    generic_residues: Mapping[str, ResidueRecord | None],
    cutoff: float = 15.0,
    modes_requested: int = 12,
) -> dict[str, float]:
    nodes = [
        (position, residue.ca)
        for position, residue in sorted(generic_residues.items())
        if residue is not None and residue.ca is not None
    ]
    if len(nodes) < 15:
        raise ValueError("Fewer than 15 mapped CA atoms are available for the ANM")
    positions = [position for position, _ in nodes]
    coords = np.stack([coordinate for _, coordinate in nodes])
    hessian = _anm_hessian(coords, cutoff)
    upper = min(hessian.shape[0] - 1, 6 + max(modes_requested * 2, 24))
    values, vectors = eigh(hessian, subset_by_index=(0, upper), check_finite=False)
    positive = np.flatnonzero(values > 1e-7)
    if len(positive) < 3:
        raise ValueError("ANM yielded fewer than three non-zero modes")
    selected = positive[:modes_requested]
    eigenvalues = values[selected]
    modes = vectors[:, selected]
    inverse = 1.0 / eigenvalues

    mode_array = modes.reshape(len(nodes), 3, len(selected))
    fluctuations = np.sum(mode_array**2 * inverse[None, None, :], axis=(1, 2))
    position_to_index = {position: index for index, position in enumerate(positions)}
    origin, x_axis, _y_axis, _z_axis = build_receptor_frame(generic_residues)

    features: dict[str, float] = {
        "mechanical_node_count": float(len(nodes)),
        "mechanical_mode_count": float(len(selected)),
        "mechanical_eigenvalue_1": float(eigenvalues[0]),
        "mechanical_eigenvalue_2": float(eigenvalues[1]) if len(eigenvalues) > 1 else float("nan"),
        "mechanical_eigenvalue_3": float(eigenvalues[2]) if len(eigenvalues) > 2 else float("nan"),
        "mechanical_total_softness": float(np.sum(inverse)),
        "mechanical_low_mode_gap": float(eigenvalues[1] - eigenvalues[0]) if len(eigenvalues) > 1 else float("nan"),
        "mechanical_fluctuation_mean": float(np.mean(fluctuations)),
        "mechanical_fluctuation_max": float(np.max(fluctuations)),
    }

    for region, prefixes in REGION_PREFIXES.items():
        indices = [index for index, position in enumerate(positions) if position.startswith(prefixes)]
        features[f"mechanical_{region}_fluctuation_mean"] = float(np.mean(fluctuations[indices])) if indices else float("nan")

    for label, position in {
        "tm6": "6.34",
        "tm5": "5.64",
        "tm7": "7.53",
        "icl2": "34.54",
    }.items():
        index = position_to_index.get(position)
        if index is None:
            features[f"mechanical_{label}_opening_softness"] = float("nan")
            features[f"mechanical_{label}_max_mode_overlap"] = float("nan")
            continue
        node_modes = mode_array[index]  # 3 x modes
        projections = x_axis @ node_modes
        features[f"mechanical_{label}_opening_softness"] = float(np.sum((projections**2) * inverse))
        features[f"mechanical_{label}_max_mode_overlap"] = float(np.max(np.abs(projections)))

    source_indices = [position_to_index[p] for p in ORTHOSTERIC_POSITIONS if p in position_to_index]
    target_indices = [position_to_index[p] for p in INTERFACE_POSITIONS if p in position_to_index]
    coupling_values = []
    for i in source_indices:
        ui = mode_array[i]
        auto_i = float(np.sum(ui**2 * inverse[None, :]))
        for j in target_indices:
            uj = mode_array[j]
            auto_j = float(np.sum(uj**2 * inverse[None, :]))
            covariance = (ui * inverse[None, :]) @ uj.T
            denominator = max(np.sqrt(auto_i * auto_j), 1e-12)
            coupling_values.append(float(np.linalg.norm(covariance, ord="fro") / denominator))
    features["mechanical_orthosteric_interface_coupling_mean"] = float(np.mean(coupling_values)) if coupling_values else float("nan")
    features["mechanical_orthosteric_interface_coupling_max"] = float(np.max(coupling_values)) if coupling_values else float("nan")

    tm6_index = position_to_index.get("6.34")
    if tm6_index is not None:
        radial = coords[tm6_index] - origin
        radial = radial - np.dot(radial, x_axis) * x_axis
        features["mechanical_tm6_radial_distance"] = float(np.linalg.norm(coords[tm6_index] - origin))
    else:
        features["mechanical_tm6_radial_distance"] = float("nan")
    return features


@lru_cache(maxsize=4)
def _load_mapping_cached(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    frame["pdb_id"] = frame["pdb_id"].astype(str).str.upper()
    return frame


def extract_structure_features(
    task: StructureTask,
    mapping_path: str | Path,
    structures_dir: str | Path,
    sasa_points: int = 80,
    contact_cutoff: float = 8.0,
    anm_cutoff: float = 15.0,
    anm_modes: int = 12,
) -> dict[str, Any]:
    structure_path = Path(structures_dir) / f"{task.pdb_id}.cif"
    if not structure_path.exists():
        structure_path = Path(structures_dir) / f"{task.pdb_id.lower()}.cif"
    if not structure_path.exists():
        raise FileNotFoundError(f"Structure file not found for {task.pdb_id}")

    mapping = _load_mapping_cached(str(Path(mapping_path).resolve()))
    keys = _mapping_for_task(mapping, task)
    if len(keys) < 15:
        raise ValueError(f"Only {len(keys)} mapped generic positions are available")

    chains = load_chains(structure_path)
    chain = _select_chain(chains, task.chain_id)
    generic_residues = _mapped_residues(chain, keys)

    row: dict[str, Any] = {
        "pdb_id": task.pdb_id,
        "receptor_name": task.receptor_name,
        "chain_id": chain.chain_id,
        "dataset_role": task.dataset_role,
        "transducer_family": task.transducer_family,
        "feature_receptor_chain_only": True,
        "feature_mapped_generic_count": int(sum(residue is not None for residue in generic_residues.values())),
        "feature_chain_residue_count": int(len(chain.residues)),
    }

    blocks = {
        "contact": lambda: contact_features(generic_residues, contact_cutoff),
        "surface": lambda: surface_features(structure_path, chain, generic_residues, sasa_points),
        "electro": lambda: electrostatic_features(chain, generic_residues),
        "mechanical": lambda: mechanical_features(generic_residues, anm_cutoff, anm_modes),
    }
    for name, function in blocks.items():
        try:
            values = function()
            row.update(values)
            row[f"availability_{name}"] = 1.0
            row[f"quality_{name}"] = "ok"
            row[f"error_{name}"] = ""
        except Exception as error:  # noqa: BLE001 - every failed feature block is audited
            row[f"availability_{name}"] = 0.0
            row[f"quality_{name}"] = f"error:{type(error).__name__}"
            row[f"error_{name}"] = str(error)
    return row


def feature_dictionary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for column in frame.columns:
        if column.startswith("contact_"):
            block = "contact_network"
        elif column.startswith("surface_"):
            block = "surface_chemistry"
        elif column.startswith("electro_"):
            block = "approximate_coulombic_electrostatics"
        elif column.startswith("mechanical_"):
            block = "anisotropic_network_mechanical_susceptibility"
        elif column.startswith("availability_"):
            block = "availability_control"
        else:
            continue
        rows.append(
            {
                "feature": column,
                "feature_block": block,
                "type": "numeric",
                "uses_transducer_atoms": False,
                "interpretation": (
                    "Approximate Coulombic descriptor; not Poisson-Boltzmann potential"
                    if column.startswith("electro_")
                    else "Mechanical susceptibility; not a free-energy estimate"
                    if column.startswith("mechanical_")
                    else "Receptor-chain-only static descriptor"
                ),
            }
        )
    return pd.DataFrame(rows)
