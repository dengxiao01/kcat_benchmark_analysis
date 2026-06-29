#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${KCATNET_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"

export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}"

cd "$ROOT"
"$PYTHON" src/24_prepare_kcatnet_eval.py
