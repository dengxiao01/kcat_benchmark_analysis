#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

python3 src/05_fetch_uniprot_sequences.py "$@"
python3 src/02_prepare_method_inputs.py

echo "UniProt sequence fetch complete."
