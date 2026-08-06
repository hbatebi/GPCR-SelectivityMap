#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run the GPCR–Galpha interface complementarity audit")
    p.add_argument("--project-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--mode", choices=("quick", "full"), default="full")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--seed", type=int, default=20272729)
    p.add_argument("--bootstraps", type=int, default=2000)
    p.add_argument("--scrambles", type=int, default=500)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--contact-cutoff", type=float, default=4.5)
    p.add_argument("--quick-limit", type=int, default=36)
    return p


def first_existing(*paths: Path) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def main() -> int:
    args = parser().parse_args()
    project = Path(args.project_dir).resolve()
    software_src = project / "software" / "gpcr-selectivitymap-v0.3.1" / "src"
    if not software_src.exists():
        raise SystemExit(f"Installed software source not found: {software_src}")
    sys.path.insert(0, str(software_src))

    from gpcr_selectivitymap.complementarity import run_complementarity_audit

    frame = project / "results/microswitches/publication_controls/publication_control_analysis_frame.csv"
    mapping = project / "results/microswitches/measurements/generic_number_mapping.csv"
    structures = project / "data/structures"
    gprotein_contacts = project / "data/processed/gprotein_contact.csv"
    audit_dir = project / "data/processed/receptor_only_audit"
    interface_positions = project / "config/interface_complementarity/intracellular_interface_positions.tsv"

    required = [frame, mapping, structures, gprotein_contacts]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing required inputs:\n" + "\n".join(missing))

    ligand = first_existing(
        audit_dir / "receptor_only_ligand_permissive_structures.csv",
        project / "results/microswitches/receptor_only_audit/receptor_only_ligand_permissive_structures.csv",
    )
    polymer = first_existing(
        audit_dir / "strict_polymer_only_structures.csv",
        project / "results/microswitches/receptor_only_audit/strict_polymer_only_structures.csv",
    )
    labels = first_existing(
        audit_dir / "receptor_only_receptor_level_labels.csv",
        project / "results/microswitches/receptor_only_audit/receptor_only_receptor_level_labels.csv",
    )

    manifest = run_complementarity_audit(
        frame=frame,
        mapping=mapping,
        structures=structures,
        output=Path(args.output_dir),
        gprotein_contacts=gprotein_contacts,
        receptor_only_ligand=ligand,
        receptor_only_polymer=polymer,
        receptor_only_labels=labels,
        interface_positions=interface_positions if interface_positions.exists() else None,
        workers=args.workers,
        folds=args.folds,
        seed=args.seed,
        bootstraps=args.bootstraps,
        scrambles=args.scrambles,
        contact_cutoff=args.contact_cutoff,
        quick=args.mode == "quick",
        quick_limit=args.quick_limit,
    )
    print("Run status:", manifest.get("status"), flush=True)
    print("Interpretation:", manifest.get("interpretation_outcome"), flush=True)
    return 0 if manifest.get("status") == "success" else 2


if __name__ == "__main__":
    raise SystemExit(main())
