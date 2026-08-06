# Fritz interface-complementarity workflow

> **Cluster configuration.** The included Slurm files use generic defaults. Run them from the project directory or set `PROJECT_DIR` and `PYTHON_BIN`; add your centre-specific account and partition options at submission time.


This workflow performs the bounded analysis discussed for the manuscript.

## Included analyses

### Experimental complex audit

- extracts receptor–Gα heavy-atom contacts at 4.5 Å
- maps receptor residues to GPCR generic positions
- maps Gα residues to a common GNAS reference coordinate system
- calculates Kyte–Doolittle hydropathy similarity and hydrophobic packing
- calculates charge, side-chain volume, aromatic, and fixed combined scores
- compares observed complex chemistry with Gα-family decoys on the same contact map

### Receptor-alone analysis

- derives Gs, Gi/o, and Gq/11 templates inside each training fold
- scores a receptor against all candidate Gα-family templates without a test G protein
- reports a sequence-only mode
- reports a modelled-residue coverage sensitivity mode
- tests ligand-permissive and strict polymer-only receptor-only structures

### Required controls

- receptor-grouped validation
- 50%, 40%, and 30% sequence-cluster holdouts
- global sequence-identity kNN
- interface-sequence logistic baseline
- receptor-position scrambling
- contact-map and Gα-chemistry scrambling
- receptor bootstrap confidence intervals
- explicit extraction, coverage, and failed-structure audits

## Installation on Fritz

Copy the complete workflow directory or the release ZIP to the project, then run:

```bash
bash manuscript_workflow/interface_complementarity_v1/install_overlay.sh
```

The default project path is:

```text
/path/to/gpcr-project
```

## Quick run

```bash
cd /path/to/gpcr-project
sbatch run_interface_complementarity_quick_fritz.sbatch
```

Quick output:

```text
results/microswitches/interface_complementarity_v1_quick
```

## Full run

```bash
sbatch run_interface_complementarity_full_fritz.sbatch
```

Full output:

```text
results/microswitches/interface_complementarity_v1_full
```

## Monitor

```bash
squeue -u "$USER"
tail -f gpcr_complementarity_JOBID.out
```

## Main outputs

```text
RUN_MANIFEST.json
interpretation_gate.json
run_summary.md
extraction/contact_pairs.tsv
extraction/galpha_profiles.tsv
extraction/contact_extraction_audit.tsv
metrics/cross_validated_metrics.tsv
predictions/cross_validated_predictions.tsv
controls/scramble_control_summary.tsv
external_receptor_only/metrics/all_functional_metrics.tsv
external_receptor_only/predictions/*_receptor_predictions.tsv
templates/
```

## Restart behaviour

Contact extraction and cross-validation folds are checkpointed. Re-submitting the same job to the same output directory reuses completed checkpoints.

## Interpretation boundary

The outputs are chemical compatibility scores. They are not validated coupling probabilities. The analysis is static and does not overlap with MD trajectory, rigid-body motion, or kinetic work.
