# Sequence-interface and phylogeny audit

Outcome: **phylogeny_dominated**

The sequence advantage weakens under low-identity transfer and is largely consistent with evolutionary organisation rather than a general interface code.

## Design

* 207 receptor-by-transducer-family rows from 169 receptors.
* 30 physically defined interface positions and 34 non-interface positions.
* Receptor-grouped and 30%, 40% and 50% sequence-cluster transfer tests in full mode.
* Nested L2 logistic regression with fold-local C tuning for the primary model comparison.
* Global sequence-identity kNN and coarse taxonomy baselines.
* Matched non-interface, positional scrambling and receptor-specific barcode scrambling controls.

## Granularity conclusion

The publication-control frame already contains one structure per receptor-by-transducer-family unit.  Receptors observed with more than one family are represented by different bound complexes, whereas sequence is identical across those rows.  A separate single-transducer-receptor sensitivity analysis is included.

## Label provenance

Bound-class labels are direct assignments from the transducer physically present in each deposited complex.  Receptor-only functional labels are audited separately, and the current consensus table does not resolve direct measurement versus homology inference for every family label.

## Interpretation limits

* The interface region is physically defined, but coefficient localisation remains descriptive.
* Bound-complex labels are structural observations, whereas receptor-only functional labels require separate provenance review.
* The current representative frame already contains one structure per receptor-by-transducer-family unit; promiscuous receptors cannot be collapsed to one multiclass row without changing the target to multilabel.
* A negative low-identity result does not show that interface chemistry is biologically irrelevant.
