# Phylogeny-aware benchmarking

## Why receptor-grouped validation is not sufficient

Holding out receptor identifiers prevents the same receptor from appearing in training and test folds. It does not prevent a close homologue from appearing in the training set. Sequence models can therefore perform well by transferring labels along evolutionary neighbourhoods.

`cluster-benchmark` constructs connected components at user-specified global sequence-identity thresholds and holds out entire components. This asks whether a representation transfers beyond close homologues.

## Recommended interpretation hierarchy

1. **Receptor-grouped association**: the representation is associated with the target in the sampled evolutionary neighbourhood.
2. **Low-homology transfer**: the association persists when close homologues are removed.
3. **Beyond sequence proximity**: the representation exceeds global identity kNN.
4. **Interface specificity**: predefined interface positions exceed matched non-interface positions under low-homology transfer.
5. **Functional transfer**: a frozen representation predicts independently curated receptor-only functional labels.

Only the later levels support strong receptor-intrinsic claims.

## Connected-component clustering

At each identity threshold, receptors are linked if their available generic-position sequence identity is at least the threshold. Transitive links form a cluster. This conservative procedure may produce fewer than five valid multiclass folds at stringent thresholds. The software reports valid fold counts and does not silently substitute receptor-grouped folds.

## Global identity kNN

The kNN baseline knows only how similar a test receptor is to training receptors. It has no explicit interface definition. If an interface model matches this baseline, the signal is consistent with evolutionary proximity rather than a universal interface code.

## Limitations

- Generic-position coverage determines which sites contribute to sequence identity.
- Connected-component cluster definitions depend on the threshold and available positions.
- Small numbers of clusters may yield wide intervals and fewer valid folds.
- Failure to transfer is absence of evidence, not proof that interface chemistry is biologically irrelevant.
