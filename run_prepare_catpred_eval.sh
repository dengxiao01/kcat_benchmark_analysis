#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 src/12_prepare_catpred_eval.py "$@"
