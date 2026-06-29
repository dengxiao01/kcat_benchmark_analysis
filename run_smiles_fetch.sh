#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 src/07_fetch_pubchem_smiles.py "$@"
python3 src/02_prepare_method_inputs.py

echo "SMILES fetch complete."
