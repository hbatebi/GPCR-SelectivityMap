#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${SKIP_PIP_INSTALL:-0}" != "1" ]]; then
  python -m pip install -r requirements.txt
fi
python scripts/00_validate_coordinates.py --download --require-all
python scripts/01_four_receptor_atomic.py
python scripts/02_five_ligand_atomic.py
python scripts/03_five_ligand_microswitch.py
python scripts/04_five_ligand_figures.py
python scripts/05_contact_cutoff_sensitivity.py
python scripts/08_four_figures.py
python scripts/06_integrated_report.py --mode full
python scripts/09_regression_check.py --mode full
python scripts/07_export_results.py --mode full
printf '\nSUCCESS. Upload GPCR_Ligand_Function_v5_COMPLETE_RESULTS.zip (in this directory).\n'
