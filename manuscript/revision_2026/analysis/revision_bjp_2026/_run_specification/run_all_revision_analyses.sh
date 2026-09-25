#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
MODE="full"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) ROOT="$(cd "$2" && pwd)"; shift 2 ;;
    --quick) MODE="quick"; shift ;;
    --full) MODE="full"; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

SUITE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"
OUT="$ROOT/revision_bjp_2026"
mkdir -p "$OUT/logs" "$ROOT/revision_inputs"
# Install blank, editable input templates only when the user has not already curated them.
[[ -e "$ROOT/revision_inputs/fusion_topology_overrides.tsv" ]] || cp "$SUITE/revision_inputs/fusion_topology_overrides_TEMPLATE.tsv" "$ROOT/revision_inputs/fusion_topology_overrides.tsv"
[[ -e "$ROOT/revision_inputs/functional_coupling_annotations.tsv" ]] || cp "$SUITE/revision_inputs/functional_coupling_annotations_TEMPLATE.tsv" "$ROOT/revision_inputs/functional_coupling_annotations.tsv"

# Copy the exact specification/code snapshot used by this run into the output tree.
mkdir -p "$OUT/_run_specification"
cp -f "$SUITE/PRESPEC.md" "$OUT/_run_specification/PRESPEC.md"
cp -f "$SUITE"/*.py "$OUT/_run_specification/" 2>/dev/null || true
cp -f "$SUITE"/*.sh "$OUT/_run_specification/" 2>/dev/null || true

if [[ "$MODE" == "quick" ]]; then
  Q="--quick"
  PERM="20"
else
  Q=""
  PERM="1000"
fi

run_stage () {
  local name="$1"; shift
  echo "[$(date -Is)] START $name" | tee -a "$OUT/RUN_STATUS.log"
  "$@" 2>&1 | tee "$OUT/logs/${name}.log"
  echo "[$(date -Is)] DONE  $name" | tee -a "$OUT/RUN_STATUS.log"
}

cd "$ROOT"
run_stage 00_preflight "$PYTHON" "$SUITE/00_preflight.py" --project-root "$ROOT"
run_stage A01 "$PYTHON" "$SUITE/A01_expanded_matched_set.py" --project-root "$ROOT" $Q
run_stage A02 "$PYTHON" "$SUITE/A02_low_dimensional_sensitivity.py" --project-root "$ROOT" $Q
run_stage A03 "$PYTHON" "$SUITE/A03_archive_trained_scoring.py" --project-root "$ROOT" $Q
run_stage A04 "$PYTHON" "$SUITE/A04_structural_descriptives.py" --project-root "$ROOT"
run_stage A05 "$PYTHON" "$SUITE/A05_functional_coupling_repertoire.py" --project-root "$ROOT" $Q
run_stage A06 "$PYTHON" "$SUITE/A06_g12_g13_census.py" --project-root "$ROOT"
run_stage A07 "$PYTHON" "$SUITE/A07_deformation_diagnostics.py" --project-root "$ROOT"
run_stage A08 "$PYTHON" "$SUITE/A08_confound_audit.py" --project-root "$ROOT"
run_stage A09 "$PYTHON" "$SUITE/A09_matched_permutation_nulls.py" --project-root "$ROOT" --iterations "$PERM" $Q
run_stage A10 "$PYTHON" "$SUITE/A10_partner_free_relatedness.py" --project-root "$ROOT"
run_stage R14 "$PYTHON" "$SUITE/R14_same_receptor_408_table.py" --project-root "$ROOT"
run_stage D01 "$PYTHON" "$SUITE/D01_feature_dictionary.py" --project-root "$ROOT"
run_stage D02 "$PYTHON" "$SUITE/D02_fold_inventory.py" --project-root "$ROOT"
run_stage D03 "$PYTHON" "$SUITE/D03_matched_pdb_tables.py" --project-root "$ROOT"
run_stage 99_validate "$PYTHON" "$SUITE/99_validate_outputs.py" --project-root "$ROOT" --expected-permutations "$PERM"

(
  cd "$ROOT"
  find revision_bjp_2026 -type f -print0 | sort -z | xargs -0 sha256sum > revision_bjp_2026/SHA256SUMS.txt
)
cat > "$OUT/RUN_COMPLETE.txt" <<EOF
mode=$MODE
completed=$(date -Is)
project_root=$ROOT
python=$($PYTHON --version 2>&1)
suite=$SUITE
EOF

echo "Revision suite complete: $OUT"
