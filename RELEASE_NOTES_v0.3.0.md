# GPCR SelectivityMap v0.3.0

Version 0.3.0 adds a predefined receptor–Gα interface-complementarity workflow and receptor-alone compatibility mode.

## New

- physical receptor–Gα contact extraction from mmCIF structures
- Gα mapping to a common reference position system
- Kyte–Doolittle hydropathy similarity and hydrophobic-packing scores
- charge, side-chain-volume, aromatic, and fixed combined-physics scores
- cognate-versus-decoy complex-side chemistry audit
- training-only Gs, Gi/o, and Gq/11 templates
- receptor-alone sequence compatibility
- receptor-structure coverage sensitivity mode
- receptor-grouped and 30%, 40%, and 50% cluster validation
- identity-kNN and interface-sequence baselines
- receptor-position and contact-map/Gα-chemistry scrambling controls
- external ligand-permissive and polymer-only receptor-only evaluation
- resumable Fritz quick and full workflows

## Interpretation

Compatibility scores are not calibrated functional coupling probabilities. A receptor-alone claim requires performance beyond sequence identity and transfer to independently audited receptor-only functional labels.

## Scope

The release remains static. It does not include MD trajectories, rigid-body decomposition, time-dependent motions, or kinetics.
