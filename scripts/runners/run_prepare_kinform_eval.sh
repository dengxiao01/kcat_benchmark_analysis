#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${KINFORM_PREP_PYTHON:-python}"

cd "$ROOT"
"$PYTHON" src/21_prepare_kinform_eval.py
"$PYTHON" src/22_check_kinform_coverage.py
