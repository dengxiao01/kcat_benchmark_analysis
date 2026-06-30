#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PRETKCAT_PYTHON:-python}"

export PYTHONNOUSERSITE=1

cd "$ROOT"
"$PYTHON" src/27_prepare_pretkcat_eval.py
