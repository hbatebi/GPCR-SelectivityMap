#!/bin/bash -l
#SBATCH --job-name=gpcr-complementarity
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=23:50:00
#SBATCH --export=NONE
#SBATCH --output=gpcr_complementarity_%j.out
#SBATCH --error=gpcr_complementarity_%j.err

set -Eeuo pipefail
unset SLURM_EXPORT_ENV

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
DRIVER="$PROJECT_DIR/run_interface_complementarity_audit.py"
OUTPUT_DIR="$PROJECT_DIR/results/microswitches/interface_complementarity_v1_full"
ARCHIVE="$PROJECT_DIR/interface_complementarity_v1_full_${SLURM_JOB_ID}.tar.gz"
WORKERS="${SLURM_CPUS_PER_TASK:-64}"

fail() { echo "ERROR: $*" >&2; exit 2; }
[[ -d "$PROJECT_DIR" ]] || fail "Project directory not found: $PROJECT_DIR"
[[ -x "$PYTHON_BIN" ]] || fail "Python not found: $PYTHON_BIN"
[[ -f "$DRIVER" ]] || fail "Driver not found: $DRIVER"
[[ -d "$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1/src" ]] || fail "Software overlay is not installed"

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
ulimit -n 65535 2>/dev/null || true

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Started: $(date -Is)"
echo "Project: $PROJECT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Workers: $WORKERS"

"$PYTHON_BIN" -m py_compile "$DRIVER" \
  "$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1/src/gpcr_selectivitymap/complementarity.py"
"$PYTHON_BIN" - <<'PY'
import numpy, pandas, scipy, sklearn, Bio
from gpcr_selectivitymap.complementarity import FAMILIES, MODEL_NAMES
print("NumPy", numpy.__version__)
print("pandas", pandas.__version__)
print("SciPy", scipy.__version__)
print("scikit-learn", sklearn.__version__)
print("Biopython", Bio.__version__)
print("Families", FAMILIES)
print("Models", MODEL_NAMES)
print("Complementarity preflight passed")
PY

set +e
srun --nodes=1 --ntasks=1 --cpus-per-task="$WORKERS" \
  "$PYTHON_BIN" -u "$DRIVER" \
  --project-dir "$PROJECT_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --mode full \
  --workers "$WORKERS" \
  --seed 20272729 \
  --folds 5 \
  --bootstraps 2000 \
  --scrambles 500 \
  --contact-cutoff 4.5
STATUS=$?
set -e

tar -czf "$ARCHIVE" -C "$PROJECT_DIR" \
  "results/microswitches/interface_complementarity_v1_full"

echo "Finished: $(date -Is)"
echo "Driver exit code: $STATUS"
echo "Archive: $ARCHIVE"
exit "$STATUS"
