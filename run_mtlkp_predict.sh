#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${MTLKP_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"
INPUT="${MTLKP_INPUT:-$ROOT/data/final/mtlkp/mtlkp_kcat_input.csv}"
OUTPUT="${MTLKP_OUTPUT:-$ROOT/data/final/mtlkp/mtlkp_kcat_input_output.csv}"

export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

cd "$ROOT"
"$PYTHON" src/44_run_mtlkp_predictions.py \
  --input "$INPUT" \
  --output "$OUTPUT" \
  --task Kcat

"$PYTHON" src/41_evaluate_method_predictions.py \
  --predictions "$OUTPUT" \
  --metadata "$ROOT/data/final/mtlkp/mtlkp_kcat_input_metadata.csv" \
  --all-metadata "$ROOT/data/final/mtlkp/mtlkp_kcat_all_metadata.csv" \
  --out-rows "$ROOT/data/final/mtlkp/mtlkp_kcat_predictions_evaluated.csv" \
  --out-metrics "$ROOT/reports/tables/mtlkp_eval_metrics.csv" \
  --out-missing "$ROOT/data/final/mtlkp/mtlkp_invalid_or_unpredicted_rows.csv" \
  --out-missing-summary "$ROOT/reports/tables/mtlkp_missing_summary.csv" \
  --prediction-column prediction_log10

"$PYTHON" src/17_build_method_eval_summary.py
