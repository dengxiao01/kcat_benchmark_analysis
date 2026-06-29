#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${PRETKCAT_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"
INPUT="${PRETKCAT_INPUT:-$ROOT/data/final/pretkcat/pretkcat_kcat_input_valid_smiles.csv}"
OUTPUT="${PRETKCAT_OUTPUT:-$ROOT/data/final/pretkcat/pretkcat_kcat_input_output.csv}"
METADATA="${PRETKCAT_METADATA:-$ROOT/data/final/pretkcat/pretkcat_kcat_input_metadata.csv}"
FEATURE_CACHE="${PRETKCAT_FEATURE_CACHE:-$ROOT/data/final/pretkcat/pretkcat_feature_cache.pkl}"
MODEL_CACHE="${PRETKCAT_MODEL_CACHE:-$ROOT/data/final/pretkcat/pretkcat_extratrees_model.pkl}"
MOLGNET_MODEL="${PRETKCAT_MOLGNET_MODEL:-$ROOT/external_methods/PreTKcat/MolGNet.pt}"
PROTT5_MODEL="${PRETKCAT_PROTT5_MODEL:-$ROOT/external_methods/CataPro/models/prot_t5_xl_uniref50}"
DEVICE="${PRETKCAT_DEVICE:-cuda:0}"

export PYTHONNOUSERSITE=1
export TORCH_HOME="${TORCH_HOME:-$ROOT/external_methods/torch_cache}"
export TRANSFORMERS_OFFLINE="${TRANSFORMERS_OFFLINE:-1}"
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="$ROOT/external_methods/kcatnet_scatter_src:$ROOT/external_methods/catapro_pydeps:$ROOT/external_methods/PreTKcat/Pretrained_Model:${PYTHONPATH:-}"
export LD_LIBRARY_PATH="/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}"

cd "$ROOT"
"$PYTHON" src/28_run_pretkcat_predictions.py \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --feature-cache "$FEATURE_CACHE" \
    --model-cache "$MODEL_CACHE" \
    --molgnet-model "$MOLGNET_MODEL" \
    --prott5-model "$PROTT5_MODEL" \
    --device "$DEVICE" \
    --n-estimators "${PRETKCAT_N_ESTIMATORS:-100}" \
    --n-jobs "${PRETKCAT_N_JOBS:--1}"

"$PYTHON" src/29_evaluate_pretkcat_predictions.py \
    --predictions "$OUTPUT" \
    --metadata "$METADATA"

"$PYTHON" src/17_build_method_eval_summary.py
