#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 src/08_fetch_sabiork_kcat.py "$@"
