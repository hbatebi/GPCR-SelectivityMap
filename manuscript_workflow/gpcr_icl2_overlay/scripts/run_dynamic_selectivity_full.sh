#!/usr/bin/env bash
set -Eeuo pipefail
PROJECT_DIR="${1:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTPUT_ROOT="$PROJECT_DIR/results/receptor_only_dynamic_selectivity/v1_full"
export PYTHONPATH="$PROJECT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$PYTHON_BIN" -u "$PROJECT_DIR/scripts/run_receptor_only_dynamic_selectivity.py" \
  --project-dir "$PROJECT_DIR" \
  --output-root "$OUTPUT_ROOT" \
  --mode full \
  --workers "${WORKERS:-32}" \
  --resume
printf 'Full results: %s\n' "$OUTPUT_ROOT"
