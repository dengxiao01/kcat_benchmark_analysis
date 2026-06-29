#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
KINFORM_ROOT="$ROOT/external_methods/KinForm"
PYTHON="${KINFORM_PYTHON:-python}"
INPUT="${KINFORM_INPUT:-$ROOT/data/final/kinform/kinform_kcat_input_predictable.json}"
OUTPUT="${KINFORM_OUTPUT:-$ROOT/data/final/kinform/kinform_kcat_input_output.csv}"

export PYTHONPATH="$ROOT/external_methods/kinform_pydeps_exact:$ROOT/external_methods/kinform_pydeps:$KINFORM_ROOT/code:${PYTHONPATH:-}"
export TOKENIZERS_PARALLELISM=false
export CUDA_VISIBLE_DEVICES="${KINFORM_CUDA_VISIBLE_DEVICES:-}"

mkdir -p "$ROOT/data/final/kinform"
ln -sf "$KINFORM_ROOT/code/smiles_embeddings/smiles_transformer/vocab.pkl" "$KINFORM_ROOT/code/vocab.pkl"

cd "$KINFORM_ROOT/code"
"$PYTHON" main.py \
    --mode predict \
    --task kcat \
    --model_config KinForm-L \
    --data_path "$INPUT" \
    --save_results "$OUTPUT"

cd "$ROOT"
"$PYTHON" src/23_evaluate_kinform_predictions.py \
    --predictions "$OUTPUT" \
    --metadata "$ROOT/data/final/kinform/kinform_kcat_input_predictable_metadata.csv" \
    --all-metadata "$ROOT/data/final/kinform/kinform_feature_coverage.csv"

"$PYTHON" src/17_build_method_eval_summary.py
