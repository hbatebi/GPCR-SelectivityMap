# BJP 2026 major-revision analysis specification

Frozen for the full HPC revision run: 24 September 2026.

## Status and provenance

This is a **revision analysis specification**, not a prospective preregistration of the original study. The submitted manuscript and its frozen P0 results predate this document. During development of the revision code, quick/debug runs were used on an uploaded project bundle to catch software and data-integration errors. Those engineering runs are not used for manuscript inference. The full Alex2 run executed from this frozen code package is the analysis record to be reviewed before any scientific revision text is frozen.

No analysis below may be silently dropped because of an unfavourable result. If an analysis is not estimable, the output must state why. Unknown annotations are not converted into negative evidence.

## Populations

**P0 strict primary population.** The submitted definition is unchanged: active partner free structures for which the receptor is the sole polymer entity, matched by receptor to Gs, Gi/o or Gq/11 bound receptor-transducer units. The partner free feature vector is the within-receptor median across all eligible strict structures, not an arbitrarily selected single PDB. P0 remains the primary analysis regardless of sensitivity results.

**P1 terminal-fusion sensitivity population.** Relax only the fusion exclusion: add active partner free constructs carrying a high-confidence N-terminal or C-terminal fusion according to the mmCIF topology audit, while retaining the other strict entity exclusions. Internal-loop fusions are excluded. The actual number of additional receptors is determined by the audit and is reported as observed; no target N (including the projected 28 in the planning document) is forced.

Secondary A1 populations are: (i) N-terminal-fusion-only under the same strict entity rules; (ii) terminal fusion plus ligand-permissive extracellular polymeric ligand; and (iii) all-fusion exploratory, which may include internal fusions and is explicitly not used to replace P0 or P1.

## Shared rules

1. G protein labels are the **captured transducer class** in deposited complexes, not the receptor's complete functional coupling repertoire or intrinsic preference.
2. Gs, Gi/o and Gq/11 are the three modeled classes unless an analysis is explicitly a G12/13 census.
3. Receptor grouping prevents the same receptor entering both train and test sets.
4. Preprocessing, imputation, scaling, dimensionality reduction, feature filtering and hyperparameter selection occur inside training folds where applicable.
5. The submitted 92-feature representation is frozen for P0. P1 uses the same 92 named features for direct comparability and reports feature availability explicitly.
6. Receptor-level classification margin is probability of the captured class minus the largest competing-class probability. Multi-unit receptors are averaged before the receptor-level paired test.
7. Primary paired inference uses an exact sign-flip test when enumeration is feasible and 100,000 Monte Carlo sign flips otherwise. Bootstrap intervals use receptor-level resampling.
8. All matched endpoint language is associational. None of the archive comparisons is interpreted as a causal effect of G protein binding.
9. Outputs are written to a new `revision_bjp_2026/` tree. Submitted/frozen results are read-only.
10. Matched-model revision analyses use seed 20291729, which reproduces the submitted transducer-imprint outer folds (submitted driver seed 20272729 plus the analysis's internal 19000 offset). New P1 folds use the same seed policy.

## A01 — Expanded matched set with terminal fusions (R2.1; informs R1)

**Question.** Does the direction of the P0 bound-versus-partner-free result persist when the partner-free definition is relaxed only to admit terminal fusion constructs?

**Main sensitivity.** P1 terminal-fusion strict population above. Same 92-feature pipeline and model family as P0. Receptor-grouped outer folds with class coverage. Report actual receptor/unit counts, fold composition, bound and partner-free macro ROC AUC, receptor-level margin differences, positive/negative receptor counts, bootstrap CI, sign-flip P, and class-specific composition. Also report feature coverage.

**Independent replication subset.** Refit/evaluate the same pipeline using only receptors newly admitted by P1, if all three classes and sufficient receptor coverage permit valid grouped CV. If not estimable, report that explicitly. A simple subset of P1 OOF predictions is also recorded but is not called an independent replication.

**Secondary subsets.** N-terminal-only strict; terminal ligand-permissive; all-fusion exploratory.

**Decision rule.** Positive mean receptor margin difference plus a majority of positive receptors = direction preserved. Attenuation is reported as attenuation. Direction reversal triggers scientific reassessment before manuscript rewriting. P0 remains primary in every case.

## A02 — Low-dimensional matched sensitivity (R2.1)

Run P0 and P1 when available using: (a) the submitted four-feature activation geometry (`tm6_r350_634_distance`, `tm3_tm7_distance`, `pif_550_644_distance`, `y753_displacement`); (b) a frozen 10-feature interpretable set adding TM6 displacement, DRY geometry, ICL2 helicity/position, TM7-H8 angle and cavity mouth area; (c) PCA with K=3 and K=5 fitted inside training folds. Use a simple L2 multinomial model for the fixed low-dimensional sets. Report the same paired receptor-level metrics and macro AUCs.

## A03 — Archive-trained held-receptor scoring (R2.1)

Train on the 207-unit bound archive after excluding all matched receptors, then score the matched bound and partner-free endpoints with the same fitted model. Repeat after excluding training receptors with >=50% generic-position sequence identity to any matched receptor when sequence information permits. Because training is on bound structures, this analysis is a supporting consistency check with explicit distribution-shift bias, not primary evidence.

## A04 — Descriptive structural comparison (R2.1)

For P0 and P1, report bound-minus-partner-free changes for named interpretable receptor features: TM6 opening/displacement, TM5 intracellular distance, TM7 displacement, DRY and PIF geometry, Y7.53 displacement, ICL2 helicity/position, and intracellular cavity volume/mouth area. Receptor-level pooled summaries report median, IQR, mean and sign counts. Sign tests with FDR correction are descriptive only. Class-stratified summaries are unit-level and are not used for strong class-specific inference at small N.

## A05 — Functional coupling repertoire (R1.2; informs R2 interpretation)

Merge the local frozen functional annotations with an optional curated revision snapshot at `revision_inputs/functional_coupling_annotations.tsv`. Record Gs, Gi/o, Gq/11 and G12/13 as 1 documented positive, 0 documented negative, blank unknown. Unknown is never treated as absence. Report annotation coverage before reporting percentages. Quantify receptors/units with one, two, three or four documented families when coverage supports it; annotate P0 and P1 separately. If `primary_family` is available, report whether the captured structural class is primary/secondary and run a bound-geometry sensitivity restricted to captured-is-primary units. Run a single-family sensitivity only when all three modeled classes and adequate receptor coverage remain.

If the local snapshot is incomplete, generate `NEEDS_FUNCTIONAL_ANNOTATION.tsv`; no complete-archive claim is made until that table is curated or an authoritative snapshot is supplied.

## A06 — G12/13 structural census (R1.3)

List every G12/13 structure in the frozen archive with receptor, PDB/chain, method/resolution/date where available. Flag receptors with a strict P0 partner-free active structure and separately receptors with an eligible P1 terminal-fusion partner-free structure. This is a census, not a modeled fourth-class analysis unless class counts later justify one.

## A07 — Deformation diagnostics (R2.6)

For P0 and P1: (i) compare cosine similarity of standardized bound-minus-partner-free deformation vectors within versus between captured classes; (ii) compare deformation magnitude with within-receptor variability among the exact partner-free structures used to construct each reference; (iii) compare within-class endpoint dispersion for bound versus partner-free coordinates. These diagnostics ask whether class information at the mature endpoint requires a single coherent class-specific deformation direction. They do not estimate causal deformation.

## A08 — Confound audit (R1.1; informs P1 interpretation)

For P0 and P1, tabulate method, resolution, species, ligand-type and fusion metadata for exact matched endpoints. Relate receptor-level margin gain descriptively to signed and absolute resolution difference; report same-method, no-bound-fusion and ligand-overlap restrictions. Small-N P values in restricted subsets are sensitivity descriptors, not pass/fail significance tests. Missing engineered-construct metadata remain unknown.

## A09 — Matched permutation nulls (R2.1/R2.2)

Run grouped receptor-label permutations with refitting for bound and partner-free representations separately, and derive the paired bound-minus-partner-free AUC-difference null from the same permutation iterations. Full run uses >=1000 permutations. Report model-level null summaries and both an upper-tail and absolute two-sided paired-difference permutation P. The receptor-level sign-flip result remains the main paired inferential result.

## A10 — Partner-free relatedness (R2.2)

Within the same outer-fold structure, calculate for each matched receptor the nearest training receptor sequence identity overall and to the same captured class using generic-position residue identities. Relate nearest same-class identity descriptively to partner-free classification margin/correctness. Export the highest-magnitude partner-free coefficients. This analysis is used to interpret, not dismiss, the non-zero partner-free signal.

## R14 — Explain the 408 same-receptor comparisons (R1.4)

Export all 408 feature-by-class-pair rows exactly as produced by the submitted analysis, summarize feature/class-pair counts and FDR results, and provide a receptor/class-pair/PDB inventory showing which multi-transducer receptors can contribute. No new significance threshold is introduced.

## D01 — Feature dictionary (R2.4)

Generate a machine-readable dictionary for the frozen 92 matched features with block, quantity type, unit, generic positions inferable from the feature name, missing-data handling and source-code hashes. The generated summaries are a starting table; exact atom selections and cavity definitions are checked against `measurement.py`/`cavity.py` before manuscript freezing.

## D02 — Fold inventory (R2.5)

Export fold membership and train/test receptor/unit counts for P0 and P1 and consolidate existing fold information where available. Dedicated project audit tables remain authoritative when they contain more exact training information than can be reconstructed from OOF predictions.

## D03 — Matched PDB/chain tables (R1.5)

For P0 and P1, list receptor, captured class, bound PDB/receptor chain and **all** partner-free PDB/chains that contribute to the within-receptor median reference, with method, resolution, species, ligand and construct/fusion metadata where available. Do not describe the median reference as a single selected partner-free PDB.

## Freeze rule after HPC execution

No manuscript, SI or point-by-point scientific response is finalized until all full-run outputs are reviewed together. Any follow-up analysis must be separately documented with its rationale and must not replace an unfavorable prespecified result.
