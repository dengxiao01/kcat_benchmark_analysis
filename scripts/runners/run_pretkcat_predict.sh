#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PRETKCAT_PYTHON:-python}"
INPUT="${PRETKCAT_INPUT:-$ROOT/data/final/pretkcat/pretkcat_kcat_input_valid_smiles.csv}"
OUTPUT="${PRETKCAT_OUTPUT:-$ROOT/data/final/pretkcat/pretkcat_kcat_input_output.csv}"
METADATA="${PRETKCAT_METADATA:-$ROOT/data/final/pretkcat/pretkcat_kcat_input_metadata.csv}"
FEATURE_CACHE="${PRETKCAT_FEATURE_CACHE:-$ROOT/data/final/pretkcat/pretkcat_feature_cache.pkl}"
MODEL_CACHE="${PRETKCAT_MODEL_CACHE:-$ROOT/data/final/pretkcat/pretkcat_extratrees_model.pkl}"
MOLGNET_MODEL="${PRETKCAT_MOLGNET_MODEL:-$ROOT/external_methods/PreTKcat/MolGNet.pt}"
PROTT5_MODEL="${PRETKCAT_PROTT5_MODEL:-$ROOT/external_methods/CataPro/models/prot_t5_xl_uniref50}"
DEVICE="${PRETKCAT_DEVICE:-cuda:0}"
DIAMOND="${PRETKCAT_DIAMOND:-diamond}"
OVERLAP_POLICY="${PRETKCAT_OVERLAP_POLICY:-exact-excluded}"
OVERLAP_AUDIT="${PRETKCAT_OVERLAP_AUDIT:-$ROOT/data/final/pretkcat/training_overlap_audit.json}"
EVAL_ROWS="${PRETKCAT_EVAL_ROWS:-$ROOT/data/final/pretkcat/pretkcat_kcat_predictions_evaluated.csv}"
EVAL_METRICS="${PRETKCAT_EVAL_METRICS:-$ROOT/reports/tables/pretkcat_eval_metrics.csv}"
EVAL_MISSING="${PRETKCAT_EVAL_MISSING:-$ROOT/data/final/pretkcat/pretkcat_invalid_or_unpredicted_rows.csv}"
EVAL_MISSING_SUMMARY="${PRETKCAT_EVAL_MISSING_SUMMARY:-$ROOT/reports/tables/pretkcat_invalid_or_unpredicted_summary.csv}"
REBUILD_GLOBAL_SUMMARY="${PRETKCAT_REBUILD_GLOBAL_SUMMARY:-1}"

export PYTHONNOUSERSITE=1
export TORCH_HOME="${TORCH_HOME:-$ROOT/external_methods/torch_cache}"
export TRANSFORMERS_OFFLINE="${TRANSFORMERS_OFFLINE:-1}"
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="$ROOT/external_methods/kcatnet_scatter_src:$ROOT/external_methods/catapro_pydeps:$ROOT/external_methods/PreTKcat/Pretrained_Model:${PYTHONPATH:-}"

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
    --n-jobs "${PRETKCAT_N_JOBS:--1}" \
    --overlap-policy "$OVERLAP_POLICY" \
    --diamond "$DIAMOND" \
    --overlap-audit-output "$OVERLAP_AUDIT"

"$PYTHON" src/29_evaluate_pretkcat_predictions.py \
    --predictions "$OUTPUT" \
    --metadata "$METADATA" \
    --out-rows "$EVAL_ROWS" \
    --out-metrics "$EVAL_METRICS" \
    --out-missing "$EVAL_MISSING" \
    --out-missing-summary "$EVAL_MISSING_SUMMARY"

if [[ "$REBUILD_GLOBAL_SUMMARY" == "1" ]]; then
    "$PYTHON" src/17_build_method_eval_summary.py
fi
