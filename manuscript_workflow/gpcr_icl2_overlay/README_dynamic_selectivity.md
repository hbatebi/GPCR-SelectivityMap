# Receptor-only static selectivity extension v0.1.1

> **Cluster configuration.** The included Slurm files use generic defaults. Run them from the project directory or set `PROJECT_DIR` and `PYTHON_BIN`; add your centre-specific account and partition options at submission time.


This overlay adds a bounded, static-first extension to the existing GPCR endpoint-selectivity repository. It tests whether richer receptor-chain-only representations contain G-protein-class information beyond the existing endpoint geometry and sequence baselines.

## Implemented feature blocks

1. Generic-position contact and side-chain network descriptors.
2. Intracellular receptor-chain surface chemistry using Shrake–Rupley SASA.
3. Approximate charged-residue Coulombic descriptors in a conserved receptor-centred frame.
4. Generic-position anisotropic-network mechanical susceptibility.

Every feature is extracted from the selected receptor chain only. No G-protein, arrestin, nanobody, Fab, fusion-partner or other non-receptor atom is passed to the feature extractor.

Approximate electrostatic descriptors are **not** Poisson–Boltzmann potentials. Normal-mode descriptors are **mechanical susceptibility**, not free energy.

## Not executed

No validated receptor-only trajectories or topology files were found in the supplied repository. The package therefore does not calculate trajectory ensemble populations, water residence, water wires, kinetic networks, MSMs, MFPTs or accessibility free-energy proxies.

## Receptor-only definitions bundled with this overlay

- `receptor_only_ligand_permissive_structures.csv`: primary definition; ligands may remain.
- `strict_polymer_only_structures.csv`: conservative sensitivity definition.
- `receptor_only_receptor_level_labels.csv`: independent functional Gs, Gi/o and Gq/11 annotations.

## Installation on Fritz

Unpack the bundle on your laptop, transfer the unpacked directory, then run:

```bash
bash install_overlay.sh /path/to/gpcr-project
```

The installer only adds new files. It does not overwrite existing v0.8.1 result directories.

## Quick run

```bash
cd /path/to/gpcr-project
sbatch run_dynamic_selectivity_quick_fritz.sbatch
```

The quick run is a technical and directional gate, not manuscript inference. It uses a receptor-balanced subset, three outer folds, 100 receptor bootstraps and 100 grouped prediction-label permutations.

## Full run

Submit only after checking the quick manifest, feature ranges, exclusions, convergence and leakage audit:

```bash
sbatch run_dynamic_selectivity_full_fritz.sbatch
```

The full output directory is checkpointed and resumable:

```text
results/receptor_only_dynamic_selectivity/v1_full/
```

Submitting the same full batch file again reuses successful structure-feature checkpoints.

## Main outputs

- `input_audit/input_audit.{json,tsv,md}`
- `input_audit/feature_feasibility.tsv`
- `features/receptor_level_features.tsv`
- `features/feature_dictionary.tsv`
- `metrics/model_metrics.tsv`
- `metrics/paired_model_differences.tsv`
- `predictions/out_of_fold_predictions.tsv`
- `permutations/grouped_permutation_results.tsv`
- `fold_manifests/fold_convergence_and_leakage_audit.tsv`
- `strict_receptor_only/functional_transfer_metrics.tsv`
- `strict_receptor_only/transfer_training_audit.tsv`
- `interpretation_gate.json`
- `analysis_manifest.json`
- `run_summary.md`
- SVG and PNG figures with source tables.

## Validation and leakage controls

- Receptor-grouped outer folds.
- Training-fold-only preprocessing and feature availability selection.
- Nested training-fold tuning of L2 regularization.
- Receptor-clustered bootstrap intervals.
- Receptor label-vector permutations within equal-row-count strata.
- Held-receptor-out bound-to-receptor-only transfer: a receptor's own bound complex is excluded when its receptor-only structures are predicted.
- Explicit fold group-overlap and convergence audit.
- Unknown functional labels remain missing.

## Interpretation

The automated gate allows three outcomes:

- **Strong support:** new features beat missingness and existing geometry, add beyond sequence, and transfer to receptor-only functional labels.
- **Limited support:** new features carry coupling-associated information but do not establish a transferable receptor-intrinsic code.
- **No support:** the tested representation does not reveal a transferable signal at current coverage.

A negative result never means that biological coupling preference is absent.
