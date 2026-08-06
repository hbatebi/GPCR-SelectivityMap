# Exact manuscript sequence-interface audit

> **Cluster configuration.** The included Slurm files use generic defaults. Run them from the project directory or set `PROJECT_DIR` and `PYTHON_BIN`; add your centre-specific account and partition options at submission time.


This directory preserves the exact Fritz/NHR workflow used for the manuscript's sequence-interface, phylogeny, granularity, and label-provenance audit.

The generic package commands in `src/gpcr_selectivitymap/` are intended for reusable datasets. The driver here is project-specific and reuses the original `gpcr_icl2` publication-control feature catalogue and modelling utilities.

## Fritz project

```text
/path/to/gpcr-project
```

## Full run

```bash
sbatch run_sequence_interface_audit_fritz.sbatch
```

## Quick run

```bash
sbatch run_sequence_interface_audit_quick_fritz.sbatch
```

## Scope

The workflow evaluates full sequence, predefined interface sequence, non-interface sequence, taxonomy, global identity kNN, endpoint geometry, and sequence plus geometry under receptor and sequence-cluster holdouts. It includes matched non-interface controls, positional scrambling, granularity checks, label-provenance audits, paired receptor bootstrap intervals, and publication figures.

The deposited manuscript result was classified as `phylogeny_dominated`. That classification means the sequence advantage weakened under low-identity transfer and was consistent with evolutionary organisation rather than a universal interface code. It does not imply that interface chemistry is biologically irrelevant.
