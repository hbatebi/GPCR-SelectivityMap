from __future__ import annotations

import argparse

from .aggregation import aggregate_feature_table
from .complementarity import run_complementarity_audit, score_receptor_compatibility_from_templates
from .benchmark import run_benchmark
from .interface_audit import run_interface_audit
from .profiling import batch_profile, profile_structure
from .provenance import run_provenance_audit
from .sequence_audit import run_cluster_benchmark, run_phylogeny_baseline


def _csv_list(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gpcr-selectivitymap",
        description=(
            "Static GPCR receptor-chain profiling and phylogeny-aware, leakage-controlled "
            "representation benchmarking"
        ),
    )
    sub = p.add_subparsers(dest="command", required=True)

    one = sub.add_parser("profile", help="Profile one GPCR structure")
    one.add_argument("--structure", required=True)
    one.add_argument("--mapping", required=True)
    one.add_argument("--pdb-id", required=True)
    one.add_argument("--receptor-name", required=True)
    one.add_argument("--chain", default="")
    one.add_argument("--output", required=True)
    one.add_argument("--sasa-points", type=int, default=80)
    one.add_argument("--anm-modes", type=int, default=12)

    batch = sub.add_parser("batch", help="Profile a manifest of GPCR structures")
    batch.add_argument("--manifest", required=True)
    batch.add_argument("--mapping", required=True)
    batch.add_argument("--structures", required=True)
    batch.add_argument("--output", required=True)
    batch.add_argument("--workers", type=int, default=1)

    bench = sub.add_parser("benchmark", help="Run receptor-grouped classification on a feature table")
    bench.add_argument("--features", required=True)
    bench.add_argument("--target", required=True)
    bench.add_argument("--group", default="receptor_name")
    bench.add_argument("--config")
    bench.add_argument("--output", required=True)
    bench.add_argument("--folds", type=int, default=5)
    bench.add_argument("--seed", type=int, default=20272729)
    bench.add_argument("--bootstraps", type=int, default=500)
    bench.add_argument("--permutations", type=int, default=100)

    cluster = sub.add_parser(
        "cluster-benchmark",
        help="Compare sequence and static structure under receptor and low-homology cluster holdouts",
    )
    cluster.add_argument("--features", required=True)
    cluster.add_argument("--target", required=True)
    cluster.add_argument("--receptor-column", default="receptor_name")
    cluster.add_argument("--output", required=True)
    cluster.add_argument("--identity-thresholds", default="0.30,0.40,0.50")
    cluster.add_argument("--interface-positions")
    cluster.add_argument("--geometry-prefixes", default="scv_,geometry_,cavity_,microswitch_,intracellular_")
    cluster.add_argument("--coarse-taxonomy-columns", default="receptor_ligand_class,receptor_subfamily")
    cluster.add_argument("--fine-taxonomy-columns", default="receptor_family")
    cluster.add_argument("--models", default="")
    cluster.add_argument("--folds", type=int, default=5)
    cluster.add_argument("--seed", type=int, default=20272729)
    cluster.add_argument("--bootstraps", type=int, default=500)
    cluster.add_argument("--k-neighbors", type=int, default=5)

    phylogeny = sub.add_parser(
        "phylogeny-baseline",
        help="Evaluate global sequence-identity kNN and taxonomy-only baselines",
    )
    phylogeny.add_argument("--features", required=True)
    phylogeny.add_argument("--target", required=True)
    phylogeny.add_argument("--receptor-column", default="receptor_name")
    phylogeny.add_argument("--output", required=True)
    phylogeny.add_argument("--identity-thresholds", default="0.30,0.40,0.50")
    phylogeny.add_argument("--coarse-taxonomy-columns", default="receptor_ligand_class,receptor_subfamily")
    phylogeny.add_argument("--fine-taxonomy-columns", default="receptor_family")
    phylogeny.add_argument("--folds", type=int, default=5)
    phylogeny.add_argument("--seed", type=int, default=20272729)
    phylogeny.add_argument("--bootstraps", type=int, default=500)
    phylogeny.add_argument("--k-neighbors", type=int, default=5)

    interface = sub.add_parser(
        "interface-audit",
        help="Test predefined interface positions against non-interface and scrambling controls",
    )
    interface.add_argument("--features", required=True)
    interface.add_argument("--target", required=True)
    interface.add_argument("--receptor-column", default="receptor_name")
    interface.add_argument("--output", required=True)
    interface.add_argument("--interface-positions")
    interface.add_argument("--identity-threshold", type=float, default=0.30)
    interface.add_argument("--folds", type=int, default=5)
    interface.add_argument("--seed", type=int, default=20272729)
    interface.add_argument("--bootstraps", type=int, default=500)
    interface.add_argument("--matched-controls", type=int, default=50)
    interface.add_argument("--scrambles", type=int, default=50)

    provenance = sub.add_parser(
        "provenance-audit",
        help="Audit whether functional labels are directly measured, consensus, or homology inferred",
    )
    provenance.add_argument("--labels", required=True)
    provenance.add_argument("--output", required=True)
    provenance.add_argument("--receptor-column", default="receptor_name")
    provenance.add_argument("--target", default="transducer_family")

    aggregate = sub.add_parser(
        "aggregate",
        help="Aggregate structure-level features to receptor or receptor-by-transducer units",
    )
    aggregate.add_argument("--features", required=True)
    aggregate.add_argument("--units", required=True, help="Comma-separated grouping columns")
    aggregate.add_argument("--numeric-method", choices=("median", "mean"), default="median")
    aggregate.add_argument("--output", required=True)


    complementarity = sub.add_parser(
        "complementarity-audit",
        help="Experimental audit of GPCR–Galpha interface chemistry and receptor-alone compatibility",
    )
    complementarity.add_argument("--frame", required=True)
    complementarity.add_argument("--mapping", required=True)
    complementarity.add_argument("--structures", required=True)
    complementarity.add_argument("--output", required=True)
    complementarity.add_argument("--gprotein-contacts")
    complementarity.add_argument("--receptor-only-ligand")
    complementarity.add_argument("--receptor-only-polymer")
    complementarity.add_argument("--receptor-only-labels")
    complementarity.add_argument("--interface-positions")
    complementarity.add_argument("--workers", type=int, default=1)
    complementarity.add_argument("--folds", type=int, default=5)
    complementarity.add_argument("--seed", type=int, default=20272729)
    complementarity.add_argument("--bootstraps", type=int, default=500)
    complementarity.add_argument("--scrambles", type=int, default=50)
    complementarity.add_argument("--contact-cutoff", type=float, default=4.5)
    complementarity.add_argument("--quick", action="store_true")
    complementarity.add_argument("--quick-limit", type=int, default=36)
    complementarity.add_argument("--k-neighbors", type=int, default=5)

    compatibility = sub.add_parser(
        "receptor-compatibility",
        help="Experimentally score receptor profiles against frozen Galpha-family templates; not calibrated probabilities",
    )
    compatibility.add_argument("--receptor-table", required=True)
    compatibility.add_argument("--template-dir", required=True)
    compatibility.add_argument("--output", required=True)
    compatibility.add_argument("--interface-positions")
    compatibility.add_argument("--profile-mode", choices=("sequence", "observed", "observed_then_sequence"), default="sequence")

    return p


def main(argv: list[str] | None = None) -> int:
    a = parser().parse_args(argv)
    if a.command == "profile":
        profile_structure(
            structure=a.structure,
            mapping=a.mapping,
            pdb_id=a.pdb_id,
            receptor_name=a.receptor_name,
            chain_id=a.chain,
            output=a.output,
            sasa_points=a.sasa_points,
            anm_modes=a.anm_modes,
        )
    elif a.command == "batch":
        batch_profile(
            manifest=a.manifest,
            mapping=a.mapping,
            structures_dir=a.structures,
            output=a.output,
            workers=a.workers,
        )
    elif a.command == "benchmark":
        run_benchmark(
            features=a.features,
            target=a.target,
            group=a.group,
            config=a.config,
            output=a.output,
            folds=a.folds,
            seed=a.seed,
            bootstraps=a.bootstraps,
            permutations=a.permutations,
        )
    elif a.command in {"cluster-benchmark", "phylogeny-baseline"}:
        kwargs = dict(
            features=a.features,
            target=a.target,
            receptor_column=a.receptor_column,
            output=a.output,
            thresholds=tuple(float(v) for v in _csv_list(a.identity_thresholds)),
            folds=a.folds,
            seed=a.seed,
            bootstraps=a.bootstraps,
            coarse_taxonomy_columns=_csv_list(a.coarse_taxonomy_columns),
            fine_taxonomy_columns=_csv_list(a.fine_taxonomy_columns),
            k_neighbors=a.k_neighbors,
        )
        if a.command == "cluster-benchmark":
            kwargs.update(
                interface_positions=a.interface_positions,
                geometry_prefixes=_csv_list(a.geometry_prefixes),
                models=_csv_list(a.models) or None,
            )
            run_cluster_benchmark(**kwargs)
        else:
            run_phylogeny_baseline(**kwargs)
    elif a.command == "interface-audit":
        run_interface_audit(
            features=a.features,
            target=a.target,
            receptor_column=a.receptor_column,
            output=a.output,
            interface_positions=a.interface_positions,
            identity_threshold=a.identity_threshold,
            folds=a.folds,
            seed=a.seed,
            bootstraps=a.bootstraps,
            matched_controls=a.matched_controls,
            scrambles=a.scrambles,
        )
    elif a.command == "provenance-audit":
        run_provenance_audit(
            labels=a.labels,
            output=a.output,
            receptor_column=a.receptor_column,
            target=a.target,
        )
    elif a.command == "aggregate":
        aggregate_feature_table(
            features=a.features,
            units=_csv_list(a.units),
            output=a.output,
            numeric_method=a.numeric_method,
        )

    elif a.command == "complementarity-audit":
        run_complementarity_audit(
            frame=a.frame,
            mapping=a.mapping,
            structures=a.structures,
            output=a.output,
            gprotein_contacts=a.gprotein_contacts,
            receptor_only_ligand=a.receptor_only_ligand,
            receptor_only_polymer=a.receptor_only_polymer,
            receptor_only_labels=a.receptor_only_labels,
            interface_positions=a.interface_positions,
            workers=a.workers,
            folds=a.folds,
            seed=a.seed,
            bootstraps=a.bootstraps,
            scrambles=a.scrambles,
            contact_cutoff=a.contact_cutoff,
            quick=a.quick,
            quick_limit=a.quick_limit,
            k_neighbors=a.k_neighbors,
        )
    elif a.command == "receptor-compatibility":
        score_receptor_compatibility_from_templates(
            receptor_table=a.receptor_table,
            template_dir=a.template_dir,
            output=a.output,
            interface_positions=a.interface_positions,
            profile_mode=a.profile_mode,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
