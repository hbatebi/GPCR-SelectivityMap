#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${SKIP_PIP_INSTALL:-0}" != "1" ]]; then
  python -m pip install -r requirements.txt
fi
python scripts/00_validate_coordinates.py
python scripts/01_four_receptor_atomic.py
python scripts/02_five_ligand_atomic.py --available-only
python scripts/03_five_ligand_microswitch.py
python scripts/04_five_ligand_figures.py
python scripts/05_contact_cutoff_sensitivity.py --available-only
python scripts/08_four_figures.py
python scripts/06_integrated_report.py --mode smoke
python scripts/09_regression_check.py --mode smoke
python scripts/07_export_results.py --mode smoke
printf '\nOFFLINE TEST SUCCESS; not a complete five-ligand result.\n'
