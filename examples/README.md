# Examples

## Structure profiling manifest

`manifest_example.tsv` shows the columns expected by the `batch` command. Real mmCIF structures and generic-number mappings are not redistributed here.

## Synthetic phylogeny benchmark

The generator below creates a small table solely for testing the command line. It is not biological data and must not be used for scientific interpretation.

```bash
cd examples
python make_synthetic_sequence_benchmark.py

gpcr-selectivitymap cluster-benchmark \
  --features synthetic_sequence_benchmark.tsv \
  --target transducer_family \
  --receptor-column receptor_name \
  --identity-thresholds 1.0 \
  --geometry-prefixes geometry_ \
  --folds 3 \
  --bootstraps 20 \
  --output synthetic_results
```
