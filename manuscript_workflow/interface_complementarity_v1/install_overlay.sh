#!/usr/bin/env bash
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${1:-$(pwd)}"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

[[ -d "$PROJECT_DIR" ]] || { echo "Project directory not found: $PROJECT_DIR" >&2; exit 2; }
mkdir -p "$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1" \
         "$PROJECT_DIR/config/interface_complementarity" \
         "$PROJECT_DIR/data/processed/receptor_only_audit"

rm -rf "$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1/src"
cp -a "$REPO_ROOT/src" "$PROJECT_DIR/software/gpcr-selectivitymap-v0.3.1/"
cp "$HERE/run_interface_complementarity_audit.py" "$PROJECT_DIR/"
cp "$HERE/run_interface_complementarity_full_fritz.sbatch" "$PROJECT_DIR/"
cp "$HERE/run_interface_complementarity_quick_fritz.sbatch" "$PROJECT_DIR/"
cp "$HERE/intracellular_interface_positions.tsv" \
   "$PROJECT_DIR/config/interface_complementarity/"
cp "$REPO_ROOT/config/interface_complementarity.yaml" \
   "$PROJECT_DIR/config/interface_complementarity/"
cp "$HERE/data/"*.csv "$PROJECT_DIR/data/processed/receptor_only_audit/"
chmod +x "$PROJECT_DIR/run_interface_complementarity_audit.py" \
         "$PROJECT_DIR/run_interface_complementarity_full_fritz.sbatch" \
         "$PROJECT_DIR/run_interface_complementarity_quick_fritz.sbatch"

echo "Installed GPCR SelectivityMap v0.3.0 complementarity overlay into:"
echo "  $PROJECT_DIR"
echo "Quick submission:"
echo "  cd $PROJECT_DIR && sbatch run_interface_complementarity_quick_fritz.sbatch"
echo "Full submission:"
echo "  cd $PROJECT_DIR && sbatch run_interface_complementarity_full_fritz.sbatch"
