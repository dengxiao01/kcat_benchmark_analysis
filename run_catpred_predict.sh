#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CATPRED_ROOT="${CATPRED_ROOT:-$ROOT/external_methods/CatPred}"
INPUT_PREFIX="${INPUT_PREFIX:-$ROOT/data/final/catpred/catpred_kcat_input}"
CATPRED_CAPSULE_DIR="${CATPRED_CAPSULE_DIR:-$ROOT/external_methods/CatPred_capsule}"
CHECKPOINT_DIR="${CHECKPOINT_DIR:-$CATPRED_CAPSULE_DIR/data/pretrained/production/kcat}"
CATPRED_CONDA_ENV="${CATPRED_CONDA_ENV:-LLM_4_enzymes_env}"
CATPRED_PYDEPS="${CATPRED_PYDEPS:-$ROOT/external_methods/catpred_pydeps}"
TORCH_HOME="${TORCH_HOME:-$ROOT/external_methods/torch_cache}"
CATPRED_CACHE_PATH="${CATPRED_CACHE_PATH:-$ROOT/external_methods/catpred_esm_cache}"
CATPRED_NO_CUDA="${CATPRED_NO_CUDA:-0}"

if [[ ! -d "$CATPRED_ROOT" ]]; then
  echo "CatPred source tree is missing: $CATPRED_ROOT" >&2
  echo "Expected it under external_methods/CatPred." >&2
  exit 2
fi

if [[ ! -f "${INPUT_PREFIX}.csv" ]]; then
  echo "CatPred input is missing: ${INPUT_PREFIX}.csv" >&2
  echo "Run: bash run_prepare_catpred_eval.sh" >&2
  exit 2
fi

if [[ ! -d "$CHECKPOINT_DIR" ]]; then
  echo "CatPred kcat checkpoint directory is missing: $CHECKPOINT_DIR" >&2
  echo "Download and extract CatPred capsule_data_update.tar.gz, then set CATPRED_CAPSULE_DIR or CHECKPOINT_DIR." >&2
  exit 2
fi

(
  cd "$CATPRED_ROOT"
  export PYTHONPATH="$CATPRED_PYDEPS:$CATPRED_ROOT${PYTHONPATH:+:$PYTHONPATH}"
  export TORCH_HOME
  export CATPRED_CACHE_PATH
  predict_args=(
    predict.py
    --test_path "${INPUT_PREFIX}.csv"
    --preds_path "${INPUT_PREFIX}_output.csv"
    --checkpoint_dir "$CHECKPOINT_DIR"
    --uncertainty_method mve
    --smiles_column SMILES
    --individual_ensemble_predictions
    --protein_records_path "${INPUT_PREFIX}.json.gz"
  )
  if [[ "$CATPRED_NO_CUDA" == "1" ]]; then
    predict_args+=(--no_cuda)
  fi
  conda run -n "$CATPRED_CONDA_ENV" python ./scripts/create_pdbrecords.py \
    --data_file "${INPUT_PREFIX}.csv" \
    --out_file "${INPUT_PREFIX}.json.gz"
  conda run -n "$CATPRED_CONDA_ENV" python "${predict_args[@]}"
)

python3 src/13_evaluate_catpred_predictions.py \
  --predictions "${INPUT_PREFIX}_output.csv" \
  --metadata "${INPUT_PREFIX}_metadata.csv"
