# GPCR SelectivityMap

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21820928.svg)](https://doi.org/10.5281/zenodo.21820928)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


**GPCR SelectivityMap** is an open, auditable toolkit for static receptor profiling and phylogeny-aware benchmarking of class A G protein-coupled receptors. It integrates GPCR generic-position sequence, endpoint geometry, contact networks, intracellular surface chemistry, approximate electrostatics, and elastic-network susceptibility, then tests whether apparent signals survive receptor grouping, low-homology transfer, missingness controls, and interface-specificity audits.

The package is designed to separate **bound-class association**, **family-constrained signal**, and **transferable evidence** before biological interpretation. It does not report calibrated functional G-protein coupling probabilities.

## Experimental interface-complementarity module

The optional `complementarity-audit` and `receptor-compatibility` commands evaluate receptor–Gα hydropathy, charge, side-chain-volume, aromatic, and combined physicochemical scores. These commands are retained as transparent exploratory analyses. They are not part of the validated manuscript claim and their outputs must not be interpreted as functional coupling probabilities.

```bash
gpcr-selectivitymap complementarity-audit \
  --frame publication_control_analysis_frame.csv \
  --mapping generic_number_mapping.csv \
  --structures data/structures \
  --gprotein-contacts data/processed/gprotein_contact.csv \
  --output results/interface_complementarity
```

After frozen templates are created, receptor-alone scoring is available through:

```bash
gpcr-selectivitymap receptor-compatibility \
  --receptor-table receptors.csv \
  --template-dir templates/all_bound \
  --profile-mode sequence \
  --output receptor_compatibility
```

These are compatibility scores, not calibrated functional coupling probabilities.

## What changed in v0.2.0

Version 0.2.0 adds the controls required to distinguish transferable representation signal from receptor homology:

- 30%, 40%, and 50% sequence-identity cluster holdouts
- global sequence-identity k-nearest-neighbour baseline
- coarse taxonomy and receptor-family baselines
- predefined intracellular interface sequence model
- interface versus non-interface comparison
- matched non-interface and positional scrambling controls
- receptor-level or receptor-by-transducer aggregation
- label-provenance auditing
- automated interpretation gates that separate receptor-grouped association from low-homology transfer

The associated manuscript audit found that sequence and interface representations are predictive under conventional receptor-grouped validation but lose their advantage under stringent homology removal. Global sequence proximity matches the interface model at 30% identity, while endpoint geometry retains more information at intermediate divergence. The software therefore emphasizes representation benchmarking rather than coupling prediction.

## Scientific scope and claim boundary

The software reports static receptor descriptors and tests whether representations generalize beyond receptor identity and sequence homology. It does **not** report validated functional G-protein coupling probabilities.

Elastic-network susceptibility is a low-frequency static proxy. It is not molecular dynamics, a time-dependent rigid-body trajectory, a free-energy surface, or a kinetic rate.

A model that performs above chance under receptor-grouped validation may still be learning evolutionary proximity. Use `cluster-benchmark`, `phylogeny-baseline`, and `interface-audit` before interpreting a signal as transferable chemistry.

## Installation

```bash
python -m pip install .

# Development installation
python -m pip install -e '.[test]'
pytest
```

## Profile one receptor

A GPCRdb-derived mapping table is required. It must contain `pdb_id`, `chain_id`, `generic_number`, `author_residue_number`, and optionally `author_insertion_code`.

```bash
gpcr-selectivitymap profile \
  --structure receptor.cif \
  --mapping generic_number_mapping.csv \
  --pdb-id 3SN6 \
  --receptor-name adrb2_human \
  --chain R \
  --output example_profile
```

Outputs include tabular and JSON features, an SVG profile, and an HTML report.

## Profile a structure dataset

```bash
gpcr-selectivitymap batch \
  --manifest examples/manifest_example.tsv \
  --mapping generic_number_mapping.csv \
  --structures data/structures \
  --workers 8 \
  --output dataset_features
```

## Standard receptor-grouped benchmark

```bash
gpcr-selectivitymap benchmark \
  --features dataset_features/receptor_level_features.tsv \
  --target transducer_family \
  --group receptor_name \
  --config config/benchmark_models.example.yaml \
  --folds 5 \
  --bootstraps 500 \
  --permutations 100 \
  --output benchmark_results
```

All preprocessing is fit inside the training fold. Outputs are not functional coupling probabilities.

## Phylogeny-aware sequence versus structure benchmark

The input feature table should contain:

- one or more rows per receptor or receptor-by-transducer unit
- a receptor identifier
- the target label
- sequence columns named `gpcrdb_<position>_expected_aa`
- optional geometry columns matching the supplied prefixes
- optional taxonomy columns

```bash
gpcr-selectivitymap cluster-benchmark \
  --features publication_control_analysis_frame.csv \
  --target transducer_family \
  --receptor-column receptor_name \
  --identity-thresholds 0.30,0.40,0.50 \
  --interface-positions config/intracellular_interface_positions.tsv \
  --geometry-prefixes scv_,geometry_,cavity_,microswitch_,intracellular_ \
  --folds 5 \
  --bootstraps 2000 \
  --output results/cluster_benchmark
```

Key outputs:

```text
results/cluster_benchmark/
├── analysis_manifest.json
├── interpretation_gate.json
├── audit/
│   ├── sequence_cluster_assignments.tsv
│   └── sequence_identity_matrix.tsv
├── metrics/
│   ├── model_metrics.tsv
│   └── paired_model_differences.tsv
├── predictions/
│   └── out_of_fold_predictions.tsv
└── figures/
    ├── representation_transfer.svg
    └── representation_transfer.png
```

## Phylogeny-only baseline

```bash
gpcr-selectivitymap phylogeny-baseline \
  --features publication_control_analysis_frame.csv \
  --target transducer_family \
  --receptor-column receptor_name \
  --identity-thresholds 0.30,0.40,0.50 \
  --output results/phylogeny_baseline
```

This evaluates global sequence-identity kNN together with available taxonomy-only baselines.

## Interface specificity audit

The interface set is physically predefined in `config/intracellular_interface_positions.tsv`; it is not selected from fitted coefficients.

```bash
gpcr-selectivitymap interface-audit \
  --features publication_control_analysis_frame.csv \
  --target transducer_family \
  --receptor-column receptor_name \
  --interface-positions config/intracellular_interface_positions.tsv \
  --identity-threshold 0.30 \
  --matched-controls 50 \
  --scrambles 50 \
  --bootstraps 2000 \
  --output results/interface_audit
```

The audit compares the interface with non-interface positions and with matched and scrambled controls. Interface localization is descriptive unless it survives low-homology transfer and exceeds these controls.

## Aggregate structures to the natural target granularity

```bash
gpcr-selectivitymap aggregate \
  --features structure_level_features.tsv \
  --units receptor_name,transducer_family \
  --numeric-method median \
  --output results/receptor_transducer_aggregated
```

Use `--units receptor_name` for a receptor-level target and `--units receptor_name,transducer_family` for bound-class units.

## Label-provenance audit

```bash
gpcr-selectivitymap provenance-audit \
  --labels functional_coupling_labels.tsv \
  --receptor-column receptor_name \
  --target transducer_family \
  --output results/provenance
```

The recommended provenance fields are provided in `config/label_provenance_template.tsv`.

## Included feature modules

- structure and generic-number audit
- endpoint geometry
- generic-position contact networks
- intracellular surface chemistry
- transparent approximate Coulombic descriptors
- anisotropic-network susceptibility
- receptor-grouped benchmarking
- sequence-cluster transfer
- sequence-identity kNN and taxonomy baselines
- interface and scrambling controls
- label-provenance audit

## Deliberately excluded

- validated G-protein coupling probability claims
- molecular-dynamics trajectory analysis
- time-dependent rigid-body decomposition
- hinge or screw-axis trajectories
- Markov-state models and transition rates
- free-energy or kinetic inference from elastic-network modes

These exclusions keep the package scientifically honest and distinct from trajectory-based rigid-body and kinetic work.

## Software workflow figure

The editable workflow overview used in the manuscript is provided as:

```text
manuscript/figures/Figure_4_GPCR_SelectivityMap_workflow.svg
```

PNG and PDF exports are included beside it.

## Manuscript reproducibility

- `manuscript/` contains figure source data and editable SVG figures.
- `manuscript/sequence_phylogeny_audit/` contains the decisive phylogeny-audit summary tables and figures.
- `manuscript_workflow/gpcr_icl2_overlay/` preserves the static receptor-only extension.
- `manuscript_workflow/sequence_interface_audit_v1/` preserves the exact Fritz analysis driver and batch files used for the manuscript audit.

Large structure archives, full prediction tables, and bootstrap distributions should be deposited in a versioned Zenodo record rather than committed to GitHub.

## Citation

Please cite the associated manuscript and archived software release. See `CITATION.cff`.

## Associated manuscript

Batebi H. *Bound G proteins leave a readable imprint on GPCR structures* (submitted, 2026).
Source tables and fold assignments underlying the main figures are in `manuscript/source_data/`
and `manuscript/sequence_phylogeny_audit/source_data/`.

## Use of generative AI

See [AI_USE.md](AI_USE.md). Generative AI was used for prose editing, code assistance, and consistency
checking only. No figure or image content was AI-generated, and AI was not used to generate
results or conclusions.

## Citation

See [CITATION.cff](CITATION.cff), or use the GitHub "Cite this repository" button.

Archived at Zenodo. Cite the concept DOI [10.5281/zenodo.21820928](https://doi.org/10.5281/zenodo.21820928) for the software in general, or [10.5281/zenodo.21820929](https://doi.org/10.5281/zenodo.21820929) for the exact v0.3.1 snapshot used in the manuscript.
