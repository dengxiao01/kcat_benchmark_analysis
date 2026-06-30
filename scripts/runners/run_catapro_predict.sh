#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CATAPRO_DIR="$ROOT/external_methods/CataPro"
PYTHON="${CATAPRO_PYTHON:-python}"
INPUT="${CATAPRO_INPUT:-$ROOT/data/final/catapro/catapro_kcat_input_valid_smiles.csv}"
OUTPUT="${CATAPRO_OUTPUT:-$ROOT/data/final/catapro/catapro_kcat_input_output.csv}"
BATCH_SIZE="${CATAPRO_BATCH_SIZE:-64}"
DEVICE="${CATAPRO_DEVICE:-cuda:0}"

export HF_HOME="$CATAPRO_DIR/hf_cache"
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"

cd "$CATAPRO_DIR/inference"
"$PYTHON" predict.py \
    -inp_fpath "$INPUT" \
    -model_dpath "$CATAPRO_DIR/models" \
    -batch_size "$BATCH_SIZE" \
    -device "$DEVICE" \
    -out_fpath "$OUTPUT"

cd "$ROOT"
"$PYTHON" src/15_evaluate_catapro_predictions.py \
    --predictions "$OUTPUT" \
    --metadata "$ROOT/data/final/catapro/catapro_kcat_input_metadata.csv"
