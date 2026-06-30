#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

python3 src/07_fetch_pubchem_smiles.py "$@"
python3 src/02_prepare_method_inputs.py

echo "SMILES fetch complete."
