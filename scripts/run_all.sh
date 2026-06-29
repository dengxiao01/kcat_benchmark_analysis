#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Reproduce the analysis_results/ in the project root.
#
# Usage:
#     cd <project_root>
#     bash scripts/run_all.sh
#
# Steps:
#     1. Integrate 5 kcat prediction sources (JSON + 4 CSVs)
#     2. Convert integrated JSON to wide CSV
#     3. Add GPR column from eciML1515 model
#     4. Add database and fill_method columns
#     5. Run 10-point kcat analysis
# ----------------------------------------------------------------------------
set -euo pipefail

cd "$(dirname "$0")/.."

echo "▶ 01_integrate_kcat.py"
python3 scripts/01_integrate_kcat.py

echo
echo "▶ 02_create_kcat_csv.py"
python3 scripts/02_create_kcat_csv.py

echo
echo "▶ 03_filter_kcat_add_gpr.py"
python3 scripts/03_filter_kcat_add_gpr.py

echo
echo "▶ 04_add_database_fill.py"
python3 scripts/04_add_database_fill.py

echo
echo "▶ 05_comprehensive_analysis.py"
python3 scripts/05_comprehensive_analysis.py

echo
echo "✅ All steps complete. See analysis_results/ for outputs."
