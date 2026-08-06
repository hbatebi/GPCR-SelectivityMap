#!/bin/bash -l
#SBATCH --job-name=gpcr-comp-quick
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=03:50:00
#SBATCH --export=NONE
#SBATCH --output=gpcr_complementarity_quick_%j.out
#SBATCH --error=gpcr_complementarity_quick_%j.err

set -Eeuo pipefail
unset SLURM_EXPORT_ENV

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
DRIVER="$PROJECT_DIR/run_interface_complementarity_audit.py"
OUTPUT_DIR="$PROJECT_DIR/results/microswitches/interface_complementarity_v1_quick"
ARCHIVE="$PROJECT_DIR/interface_complementarity_v1_quick_${SLURM_JOB_ID}.tar.gz"
WORKERS="${SLURM_CPUS_PER_TASK:-24}"

cd "$PROJECT_DIR"
mkdir -p "$OUTPUT_DIR"
export PYTHONPATH="$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1/src:$PROJECT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

"$PYTHON_BIN" -m py_compile "$DRIVER" \
  "$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1/src/gpcr_selectivitymap/complementarity.py"

set +e
srun --nodes=1 --ntasks=1 --cpus-per-task="$WORKERS" \
  "$PYTHON_BIN" -u "$DRIVER" \
  --project-dir "$PROJECT_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --mode quick \
  --workers "$WORKERS" \
  --seed 20272729 \
  --folds 5 \
  --bootstraps 300 \
  --scrambles 20 \
  --contact-cutoff 4.5 \
  --quick-limit 36
STATUS=$?
set -e

tar -czf "$ARCHIVE" -C "$PROJECT_DIR" \
  "results/microswitches/interface_complementarity_v1_quick"

echo "Driver exit code: $STATUS"
echo "Archive: $ARCHIVE"
exit "$STATUS"
