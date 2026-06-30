#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${DEKP_PYTHON:-python}"

export PYTHONNOUSERSITE=1

cd "$ROOT"
"$PYTHON" src/30_prepare_dekp_eval.py
