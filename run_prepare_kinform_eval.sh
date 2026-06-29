#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${KINFORM_PREP_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"

cd "$ROOT"
"$PYTHON" src/21_prepare_kinform_eval.py
"$PYTHON" src/22_check_kinform_coverage.py
