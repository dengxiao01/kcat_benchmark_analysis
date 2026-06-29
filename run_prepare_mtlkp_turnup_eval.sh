#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${PYTHON:-python}"

cd "$ROOT"
"$PYTHON" src/43_prepare_mtlkp_turnup_eval.py
