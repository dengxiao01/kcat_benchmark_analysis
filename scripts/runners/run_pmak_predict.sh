#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PMAK_PYTHON:-python}"
INPUT="${PMAK_INPUT:-$ROOT/data/final/pmak/pmak_kcat_input.csv}"
OUTPUT="${PMAK_OUTPUT:-$ROOT/data/final/pmak/pmak_kcat_input_output.csv}"
FEATURE_CACHE="${PMAK_FEATURE_CACHE:-$ROOT/data/final/pmak/pmak_feature_cache.pkl}"
DEVICE="${PMAK_DEVICE:-cuda:0}"

export HF_HOME="$ROOT/external_methods/PMAK/hf_cache"
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"

cd "$ROOT"
"$PYTHON" src/19_run_pmak_predictions.py \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --feature-cache "$FEATURE_CACHE" \
    --device "$DEVICE"

"$PYTHON" src/20_evaluate_pmak_predictions.py \
    --predictions "$OUTPUT" \
    --metadata "$ROOT/data/final/pmak/pmak_kcat_input_metadata.csv"

"$PYTHON" src/17_build_method_eval_summary.py
