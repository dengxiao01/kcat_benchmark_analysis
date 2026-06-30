#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PMAK_PYTHON:-python}"

export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"

cd "$ROOT"
"$PYTHON" src/18_prepare_pmak_eval.py
