#!/usr/bin/env bash
#SBATCH -J catapro_smoke
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o logs/catapro_sample_%j.out
#SBATCH -e logs/catapro_sample_%j.err

set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
CATAPRO_DIR="$ROOT/external_methods/CataPro"
PYTHON="${CATAPRO_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"

export HF_HOME="$CATAPRO_DIR/hf_cache"
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"

cd "$ROOT"
mkdir -p logs
bash run_prepare_catapro_eval.sh

cd "$CATAPRO_DIR/inference"
"$PYTHON" predict.py \
    -inp_fpath "$ROOT/data/final/catapro/catapro_kcat_input_sample20.csv" \
    -model_dpath "$CATAPRO_DIR/models" \
    -batch_size 8 \
    -device cuda:0 \
    -out_fpath "$ROOT/data/final/catapro/catapro_kcat_input_sample20_output.csv"

cd "$ROOT"
"$PYTHON" src/15_evaluate_catapro_predictions.py \
    --predictions "$ROOT/data/final/catapro/catapro_kcat_input_sample20_output.csv" \
    --metadata "$ROOT/data/final/catapro/catapro_kcat_input_sample20_metadata.csv" \
    --out-rows "$ROOT/data/final/catapro/catapro_kcat_input_sample20_predictions_evaluated.csv" \
    --out-metrics "$ROOT/reports/tables/catapro_sample20_eval_metrics.csv" \
    --out-missing "$ROOT/data/final/catapro/catapro_sample20_invalid_or_unpredicted_rows.csv" \
    --out-missing-summary "$ROOT/reports/tables/catapro_sample20_invalid_or_unpredicted_summary.csv"
