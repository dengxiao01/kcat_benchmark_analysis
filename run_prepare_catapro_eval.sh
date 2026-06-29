#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${CATAPRO_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"

cd "$ROOT"
"$PYTHON" src/14_prepare_catapro_eval.py
"$PYTHON" src/16_filter_catapro_valid_smiles.py
