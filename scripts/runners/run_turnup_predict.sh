#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${TURNUP_PYTHON:-python}"
INPUT="${TURNUP_INPUT:-$ROOT/data/final/turnup/turnup_kcat_input.csv}"
OUTPUT="${TURNUP_OUTPUT:-$ROOT/data/final/turnup/turnup_kcat_input_output.csv}"
CACHE="${TURNUP_ESM_CACHE:-$ROOT/data/final/turnup/turnup_esm1b_cache.pkl}"
DEVICE="${TURNUP_DEVICE:-cuda:0}"

export PYTHONNOUSERSITE=1

cd "$ROOT"
"$PYTHON" src/45_run_turnup_predictions.py \
  --input "$INPUT" \
  --output "$OUTPUT" \
  --esm-cache "$CACHE" \
  --device "$DEVICE"

"$PYTHON" src/41_evaluate_method_predictions.py \
  --predictions "$OUTPUT" \
  --metadata "$ROOT/data/final/turnup/turnup_kcat_input_metadata.csv" \
  --all-metadata "$ROOT/data/final/turnup/turnup_kcat_all_metadata.csv" \
  --out-rows "$ROOT/data/final/turnup/turnup_kcat_predictions_evaluated.csv" \
  --out-metrics "$ROOT/reports/tables/turnup_eval_metrics.csv" \
  --out-missing "$ROOT/data/final/turnup/turnup_invalid_or_unpredicted_rows.csv" \
  --out-missing-summary "$ROOT/reports/tables/turnup_missing_summary.csv" \
  --prediction-column prediction_log10

"$PYTHON" src/17_build_method_eval_summary.py
