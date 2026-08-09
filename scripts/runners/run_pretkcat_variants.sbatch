#!/usr/bin/env bash
#SBATCH -J pretkcat_variants
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/pretkcat_variants_%j.out
#SBATCH -e logs/pretkcat_variants_%j.err

set -euo pipefail

ROOT="${KCAT_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
PYTHON="${PRETKCAT_PYTHON:-python}"
cd "$ROOT"
mkdir -p logs data/final/pretkcat/variants reports/tables/pretkcat_variants

export PRETKCAT_PYTHON="$PYTHON"
export PRETKCAT_DIAMOND="${PRETKCAT_DIAMOND:-diamond}"
export PRETKCAT_REBUILD_GLOBAL_SUMMARY=0

bash "$ROOT/scripts/runners/run_prepare_pretkcat_eval.sh"

run_variant() {
    local variant="$1"
    local output="$2"
    local model="$3"
    local audit="$4"
    local rows="$5"
    local metrics="$6"
    local missing="$7"
    local missing_summary="$8"

    PRETKCAT_OVERLAP_POLICY="$variant" \
    PRETKCAT_OUTPUT="$output" \
    PRETKCAT_MODEL_CACHE="$model" \
    PRETKCAT_OVERLAP_AUDIT="$audit" \
    PRETKCAT_EVAL_ROWS="$rows" \
    PRETKCAT_EVAL_METRICS="$metrics" \
    PRETKCAT_EVAL_MISSING="$missing" \
    PRETKCAT_EVAL_MISSING_SUMMARY="$missing_summary" \
    bash "$ROOT/scripts/runners/run_pretkcat_predict.sh"
}

run_variant \
    exact-excluded \
    "$ROOT/data/final/pretkcat/pretkcat_kcat_input_output.csv" \
    "$ROOT/data/final/pretkcat/variants/exact-excluded/model.pkl" \
    "$ROOT/data/final/pretkcat/variants/exact-excluded/training_overlap_audit.json" \
    "$ROOT/data/final/pretkcat/pretkcat_kcat_predictions_evaluated.csv" \
    "$ROOT/reports/tables/pretkcat_eval_metrics.csv" \
    "$ROOT/data/final/pretkcat/pretkcat_invalid_or_unpredicted_rows.csv" \
    "$ROOT/reports/tables/pretkcat_invalid_or_unpredicted_summary.csv"

for variant in raw-public near-excluded; do
    variant_dir="$ROOT/data/final/pretkcat/variants/$variant"
    report_prefix="${variant//-/_}"
    mkdir -p "$variant_dir"
    run_variant \
        "$variant" \
        "$variant_dir/predictions.csv" \
        "$variant_dir/model.pkl" \
        "$variant_dir/training_overlap_audit.json" \
        "$variant_dir/predictions_evaluated.csv" \
        "$ROOT/reports/tables/pretkcat_variants/${report_prefix}_eval_metrics.csv" \
        "$variant_dir/missing_rows.csv" \
        "$ROOT/reports/tables/pretkcat_variants/${report_prefix}_missing_summary.csv"
done

"$PYTHON" src/17_build_method_eval_summary.py
