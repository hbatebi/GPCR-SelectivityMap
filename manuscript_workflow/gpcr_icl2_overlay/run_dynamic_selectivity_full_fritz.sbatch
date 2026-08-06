#!/bin/bash -l
#SBATCH --job-name=gpcr-ro-dyn-full
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=23:50:00
#SBATCH --export=NONE
#SBATCH --output=gpcr_ro_dynamic_full_%j.out
#SBATCH --error=gpcr_ro_dynamic_full_%j.err

set -Eeuo pipefail
unset SLURM_EXPORT_ENV
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
OUTPUT_ROOT="$PROJECT_DIR/results/receptor_only_dynamic_selectivity/v1_full"
[[ -d "$PROJECT_DIR" ]] || { echo "Missing project: $PROJECT_DIR" >&2; exit 2; }
[[ -x "$PYTHON_BIN" ]] || { echo "Missing Python: $PYTHON_BIN" >&2; exit 2; }
cd "$PROJECT_DIR"
export PATH="$(dirname "$PYTHON_BIN"):$PATH"
export PYTHONPATH="$PROJECT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
"$PYTHON_BIN" -m pytest -q tests/test_dynamic_selectivity.py
srun --nodes=1 --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK:-32}" \
  "$PYTHON_BIN" -u scripts/run_receptor_only_dynamic_selectivity.py \
  --project-dir "$PROJECT_DIR" \
  --output-root "$OUTPUT_ROOT" \
  --mode full \
  --workers "${SLURM_CPUS_PER_TASK:-32}" \
  --resume
ARCHIVE="$PROJECT_DIR/receptor_only_dynamic_selectivity_v1_full_${SLURM_JOB_ID}.tar.gz"
tar -czf "$ARCHIVE" -C "$PROJECT_DIR" "results/receptor_only_dynamic_selectivity/v1_full"
echo "Results: $OUTPUT_ROOT"
echo "Archive: $ARCHIVE"
