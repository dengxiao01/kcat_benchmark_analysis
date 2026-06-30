#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${CATAPRO_PYTHON:-python}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"

cd "$ROOT"
"$PYTHON" src/14_prepare_catapro_eval.py
"$PYTHON" src/16_filter_catapro_valid_smiles.py
