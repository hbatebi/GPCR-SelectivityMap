#!/bin/bash -l
#SBATCH --job-name=gpcr-seq-interface
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=23:50:00
#SBATCH --export=NONE
#SBATCH --output=gpcr_sequence_interface_%j.out
#SBATCH --error=gpcr_sequence_interface_%j.err

set -Eeuo pipefail
unset SLURM_EXPORT_ENV

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
DRIVER="$PROJECT_DIR/run_sequence_interface_audit.py"
OUTPUT_DIR="$PROJECT_DIR/results/microswitches/sequence_interface_audit_v1"
ARCHIVE="$PROJECT_DIR/sequence_interface_audit_v1_${SLURM_JOB_ID}.tar.gz"

SEED=20272729
BOOTSTRAPS=2000
CONTROL_SETS=50
SCRAMBLES=50
COEFFICIENT_PERMUTATIONS=10000
WORKERS="${SLURM_CPUS_PER_TASK:-64}"

fail() {
    echo "ERROR: $*" >&2
    exit 2
}

[[ -d "$PROJECT_DIR" ]] || fail "Project directory not found: $PROJECT_DIR"
[[ -x "$PYTHON_BIN" ]] || fail "Python executable not found: $PYTHON_BIN"
[[ -f "$DRIVER" ]] || fail "Driver not found: $DRIVER"
[[ -f "$PROJECT_DIR/results/microswitches/publication_controls/publication_control_analysis_frame.csv" ]] || \
    fail "Publication-control analysis frame is missing"

cd "$PROJECT_DIR"
mkdir -p "$OUTPUT_DIR"

export PATH="$(dirname "$PYTHON_BIN"):$PATH"
export PYTHONPATH="$PROJECT_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg

# Every worker is a separate Python process.  Prevent BLAS and OpenMP from
# creating hidden threads inside each worker.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

ulimit -n 65535 2>/dev/null || true

cleanup() {
    local rc=$?
    jobs -pr | xargs -r kill 2>/dev/null || true
    exit "$rc"
}
trap cleanup INT TERM

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Started: $(date -Is)"
echo "Project: $PROJECT_DIR"
echo "Python: $PYTHON_BIN"
echo "Workers: $WORKERS"
echo "Output: $OUTPUT_DIR"

"$PYTHON_BIN" --version
"$PYTHON_BIN" -m py_compile "$DRIVER"
"$PYTHON_BIN" - <<'PY'
import inspect
import numpy
import pandas
import scipy
import sklearn
from gpcr_icl2.microswitches import publication_controls

required = (
    "cross_validate_experiment",
    "build_feature_catalog",
    "add_generic_sequence_clusters",
)
for name in required:
    assert hasattr(publication_controls, name), name

print("NumPy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("SciPy:", scipy.__version__)
print("scikit-learn:", sklearn.__version__)
print("publication_controls:", inspect.getfile(publication_controls))
print("Sequence-interface audit preflight passed")
PY

set +e
srun --nodes=1 --ntasks=1 --cpus-per-task="$WORKERS" \
    "$PYTHON_BIN" -u "$DRIVER" \
    --project-dir "$PROJECT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --mode full \
    --workers "$WORKERS" \
    --seed "$SEED" \
    --bootstraps "$BOOTSTRAPS" \
    --control-sets "$CONTROL_SETS" \
    --scrambles "$SCRAMBLES" \
    --coefficient-permutations "$COEFFICIENT_PERMUTATIONS" \
    --dpi 600 \
    --resume
DRIVER_STATUS=$?
set -e

# Archive all completed checkpoints and summaries even when a convergence gate
# returns a non-zero status.  This preserves diagnostically useful results.
tar -czf "$ARCHIVE" -C "$PROJECT_DIR" \
    "results/microswitches/sequence_interface_audit_v1"

echo "Finished: $(date -Is)"
echo "Driver exit code: $DRIVER_STATUS"
echo "Results: $OUTPUT_DIR"
echo "Archive: $ARCHIVE"

exit "$DRIVER_STATUS"
