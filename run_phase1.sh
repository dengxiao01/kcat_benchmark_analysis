#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 src/01_parse_models.py
python3 src/02_prepare_method_inputs.py
python3 src/03_match_experimental_kcat.py
python3 src/04_build_curation_queues.py

echo "Phase 1 local tables complete."
