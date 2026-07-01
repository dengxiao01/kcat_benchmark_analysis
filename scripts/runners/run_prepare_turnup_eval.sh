#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PYTHON:-python}"

cd "$ROOT"
"$PYTHON" src/43_prepare_turnup_eval.py
