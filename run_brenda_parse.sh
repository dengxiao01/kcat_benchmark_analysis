#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 src/10_parse_brenda_kcat.py "$@"
