#!/usr/bin/env bash
set -Eeuo pipefail
PROJECT_DIR="${1:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_ROOT="$PROJECT_DIR/results/receptor_only_dynamic_selectivity/quick_${STAMP}"
export PYTHONPATH="$PROJECT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$PYTHON_BIN" -u "$PROJECT_DIR/scripts/run_receptor_only_dynamic_selectivity.py" \
  --project-dir "$PROJECT_DIR" \
  --output-root "$OUTPUT_ROOT" \
  --mode quick \
  --workers "${WORKERS:-8}" \
  --resume
printf 'Quick results: %s\n' "$OUTPUT_ROOT"
