# Receptor–Gα interface complementarity

Version 0.3.0 adds a predefined, two-sided interface chemistry audit and a separate receptor-alone compatibility mode.

## Scientific question

The analysis asks whether the chemistry of the receptor interface is more compatible with the experimentally observed Gα family than with alternative Gα-family templates, and whether that compatibility transfers beyond receptor phylogeny.

It is not a search over arbitrary residue scales. The primary descriptors are fixed before evaluation:

- Kyte–Doolittle hydropathy similarity
- hydrophobic packing
- charge complementarity
- side-chain volume similarity
- aromatic pairing
- a fixed combined-physics score

## Complex-side mode

Physical receptor–Gα residue contacts are extracted at a 4.5 Å heavy-atom cutoff. Gα chains are aligned to a common GNAS reference coordinate system. Each experimental contact map is scored with the observed Gα chemistry and with training-family decoy chemistry.

This mode uses the test complex contact graph and is therefore a mechanistic audit, not receptor-alone prediction.

## Receptor-alone mode

Gα-family templates are derived from training complexes only. A test receptor is scored against Gs, Gi/o, and Gq/11 templates without using a test G protein.

Two receptor-alone profiles are reported:

1. `sequence`: expected amino acids at predefined intracellular generic positions. This requires no receptor structure.
2. `observed`: only residues modelled in a receptor structure. This is a structure-coverage sensitivity analysis and does not encode dynamics or side-chain exposure.

The output values are compatibility scores. They are not calibrated functional coupling probabilities.

## Validation

All templates and baselines are reconstructed inside each outer training fold. The package reports:

- receptor-grouped validation
- 50%, 40%, and 30% sequence-cluster holdouts
- global sequence-identity kNN
- interface-sequence logistic baseline
- receptor-position scrambling
- contact-map and Gα-chemistry scrambling
- ligand-permissive and strict polymer-only receptor-only transfer

## Scope boundary

The module is static. It does not analyse trajectories, time-dependent helix motion, rigid-body decomposition, transition kinetics, or dynamic coupling pathways.
