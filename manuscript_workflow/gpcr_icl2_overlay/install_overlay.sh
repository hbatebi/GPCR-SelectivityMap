#!/usr/bin/env bash
set -Eeuo pipefail
PROJECT_DIR="${1:-$(pwd)}"
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -d "$PROJECT_DIR/src/gpcr_icl2" ]] || { echo "Not a GPCR project: $PROJECT_DIR" >&2; exit 2; }
mkdir -p \
  "$PROJECT_DIR/src/gpcr_icl2/receptor_only_dynamic_selectivity" \
  "$PROJECT_DIR/scripts" \
  "$PROJECT_DIR/config" \
  "$PROJECT_DIR/tests" \
  "$PROJECT_DIR/data/manual"
cp -a "$SOURCE_DIR/src/gpcr_icl2/receptor_only_dynamic_selectivity/." "$PROJECT_DIR/src/gpcr_icl2/receptor_only_dynamic_selectivity/"
cp -a "$SOURCE_DIR/scripts/." "$PROJECT_DIR/scripts/"
cp -a "$SOURCE_DIR/config/receptor_only_dynamic_selectivity.yaml" "$PROJECT_DIR/config/"
cp -a "$SOURCE_DIR/tests/test_dynamic_selectivity.py" "$PROJECT_DIR/tests/"
cp -a "$SOURCE_DIR/data/manual/." "$PROJECT_DIR/data/manual/"
cp -a "$SOURCE_DIR/run_dynamic_selectivity_quick_fritz.sbatch" "$PROJECT_DIR/"
cp -a "$SOURCE_DIR/run_dynamic_selectivity_full_fritz.sbatch" "$PROJECT_DIR/"
chmod +x \
  "$PROJECT_DIR/scripts/run_receptor_only_dynamic_selectivity.py" \
  "$PROJECT_DIR/scripts/run_dynamic_selectivity_quick.sh" \
  "$PROJECT_DIR/scripts/run_dynamic_selectivity_full.sh" \
  "$PROJECT_DIR/run_dynamic_selectivity_quick_fritz.sbatch" \
  "$PROJECT_DIR/run_dynamic_selectivity_full_fritz.sbatch"
echo "Installed receptor-only dynamic-selectivity overlay into: $PROJECT_DIR"
