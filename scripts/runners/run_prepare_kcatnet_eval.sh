#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${KCATNET_PYTHON:-python}"

export PYTHONNOUSERSITE=1

cd "$ROOT"
"$PYTHON" src/24_prepare_kcatnet_eval.py
