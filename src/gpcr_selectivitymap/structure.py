from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa

from .geometry import BACKBONE, sidechain_centroid

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


@dataclass(frozen=True, order=True)
class ResidueKey:
    resseq: int
    icode: str = ""

    def label(self) -> str:
        return f"{self.resseq}{self.icode.strip()}"


@dataclass
class ResidueRecord:
    chain_id: str
    key: ResidueKey
    resname: str
    atoms: dict[str, np.ndarray]
    occupancies: dict[str, float]
    altlocs: dict[str, str]

    @property
    def one_letter(self) -> str:
        return AA3_TO_1.get(self.resname.upper(), "X")

    @property
    def ca(self) -> np.ndarray | None:
        return self.atoms.get("CA")

    @property
    def sidechain_center(self) -> np.ndarray | None:
        return sidechain_centroid(self.resname, self.atoms)

    @property
    def missing_sidechain(self) -> bool:
        if self.resname in {"GLY", "ALA"}:
            return False
        return not any(name not in BACKBONE and not name.startswith("H") for name in self.atoms)


@dataclass
class ChainRecord:
    chain_id: str
    residues: list[ResidueRecord]

    @property
    def sequence(self) -> str:
        return "".join(r.one_letter for r in self.residues)

    def by_key(self) -> dict[ResidueKey, ResidueRecord]:
        return {r.key: r for r in self.residues}


@dataclass
class AlignmentResult:
    chain_id: str
    uniprot_to_key: dict[int, ResidueKey]
    identity: float
    receptor_coverage: float
    observed_coverage: float
    aligned_count: int
    score: float
    mismatches: list[tuple[int, str, str]]


def _best_atom(atom) -> tuple[np.ndarray, float, str]:
    if atom.is_disordered() and hasattr(atom, "child_dict"):
        choices = list(atom.child_dict.values())
        atom = max(choices, key=lambda a: ((a.get_occupancy() or 0.0), str(a.get_altloc())))
    occ = float(atom.get_occupancy() or 0.0)
    return np.asarray(atom.get_coord(), dtype=float), occ, str(atom.get_altloc() or "").strip()


def load_chains(path: str | Path) -> dict[str, ChainRecord]:
    path = Path(path)
    suffix = path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models())
    chains: dict[str, ChainRecord] = {}
    for chain in model:
        residues: list[ResidueRecord] = []
        for residue in chain:
            if not is_aa(residue, standard=False):
                continue
            hetflag, resseq, icode = residue.id
            if str(hetflag).strip() not in {"", " "} and residue.get_resname().upper() not in {"MSE"}:
                continue
            atoms: dict[str, np.ndarray] = {}
            occupancies: dict[str, float] = {}
            altlocs: dict[str, str] = {}
            for atom in residue:
                name = str(atom.get_name()).strip()
                coord, occ, alt = _best_atom(atom)
                if name not in atoms or occ > occupancies[name]:
                    atoms[name] = coord
                    occupancies[name] = occ
                    altlocs[name] = alt
            if "CA" not in atoms:
                continue
            residues.append(ResidueRecord(
                chain_id=str(chain.id).strip(),
                key=ResidueKey(int(resseq), str(icode).strip()),
                resname=str(residue.get_resname()).upper(),
                atoms=atoms,
                occupancies=occupancies,
                altlocs=altlocs,
            ))
        if residues:
            chains[str(chain.id).strip()] = ChainRecord(str(chain.id).strip(), residues)
    if not chains:
        raise ValueError(f"No amino-acid chains with CA atoms found in {path}")
    return chains


def align_receptor_to_chain(uniprot_seq: Mapping[int, str], chain: ChainRecord) -> AlignmentResult:
    unp_numbers = sorted(int(x) for x in uniprot_seq)
    target = "".join(str(uniprot_seq[n]).upper()[0] for n in unp_numbers)
    query = chain.sequence
    if not target or not query:
        return AlignmentResult(chain.chain_id, {}, 0.0, 0.0, 0.0, 0, float("-inf"), [])
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.open_end_gap_score = 0
    aligner.extend_end_gap_score = 0
    alignment = aligner.align(target, query)[0]
    mapping: dict[int, ResidueKey] = {}
    mismatches: list[tuple[int, str, str]] = []
    identities = 0
    aligned = 0
    for (ts, te), (qs, qe) in zip(*alignment.aligned):
        length = min(te - ts, qe - qs)
        for offset in range(length):
            ti = ts + offset
            qi = qs + offset
            number = unp_numbers[ti]
            observed = chain.residues[qi].one_letter
            expected = target[ti]
            mapping[number] = chain.residues[qi].key
            aligned += 1
            if expected == observed:
                identities += 1
            else:
                mismatches.append((number, expected, observed))
    identity = identities / aligned if aligned else 0.0
    receptor_coverage = aligned / len(target) if target else 0.0
    observed_coverage = aligned / len(query) if query else 0.0
    composite = identity * 0.65 + receptor_coverage * 0.20 + observed_coverage * 0.15
    return AlignmentResult(
        chain.chain_id, mapping, identity, receptor_coverage, observed_coverage,
        aligned, composite, mismatches,
    )


def select_receptor_chain(chains: Mapping[str, ChainRecord], uniprot_seq: Mapping[int, str],
                          preferred_chain: str | None = None,
                          minimum_identity: float = 0.35,
                          minimum_aligned: int = 70) -> tuple[ChainRecord, AlignmentResult]:
    results = [align_receptor_to_chain(uniprot_seq, chain) for chain in chains.values()]
    results.sort(
        key=lambda r: (
            r.score,
            r.aligned_count,
            r.chain_id == str(preferred_chain or "").strip(),
            r.chain_id,
        ), reverse=True,
    )
    best = results[0]
    if best.identity < minimum_identity or best.aligned_count < minimum_aligned:
        raise ValueError(
            f"No credible receptor chain: best={best.chain_id}, identity={best.identity:.3f}, "
            f"aligned={best.aligned_count}"
        )
    return chains[best.chain_id], best


def mapped_residues(chain: ChainRecord, uniprot_to_key: Mapping[int, ResidueKey],
                    generic_to_uniprot: Mapping[str, int | None]) -> dict[str, ResidueRecord | None]:
    by_key = chain.by_key()
    out: dict[str, ResidueRecord | None] = {}
    for generic, sequence_number in generic_to_uniprot.items():
        if sequence_number is None:
            out[generic] = None
            continue
        key = uniprot_to_key.get(int(sequence_number))
        out[generic] = by_key.get(key) if key is not None else None
    return out


def observed_ions(path: str | Path, ion_names: Iterable[str] = ("NA", "LI", "K", "CA", "MG")) -> list[tuple[str, np.ndarray]]:
    path = Path(path)
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models())
    wanted = {x.upper() for x in ion_names}
    ions: list[tuple[str, np.ndarray]] = []
    for residue in model.get_residues():
        name = str(residue.get_resname()).upper().strip()
        if name not in wanted:
            continue
        for atom in residue:
            coord, _occ, _alt = _best_atom(atom)
            ions.append((name, coord))
            break
    return ions
