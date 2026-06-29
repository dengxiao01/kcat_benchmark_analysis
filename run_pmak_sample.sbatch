#!/usr/bin/env bash
#SBATCH -J pmak_smoke
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -o logs/pmak_sample_%j.out
#SBATCH -e logs/pmak_sample_%j.err

set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${PMAK_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"

export HF_HOME="$ROOT/external_methods/PMAK/hf_cache"
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"
export LD_LIBRARY_PATH="/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}"

cd "$ROOT"
mkdir -p logs
bash run_prepare_pmak_eval.sh

"$PYTHON" src/19_run_pmak_predictions.py \
    --input "$ROOT/data/final/pmak/pmak_kcat_input_sample20.csv" \
    --output "$ROOT/data/final/pmak/pmak_kcat_input_sample20_output.csv" \
    --feature-cache "$ROOT/data/final/pmak/pmak_feature_cache_sample20.pkl" \
    --device cuda:0

"$PYTHON" src/20_evaluate_pmak_predictions.py \
    --predictions "$ROOT/data/final/pmak/pmak_kcat_input_sample20_output.csv" \
    --metadata "$ROOT/data/final/pmak/pmak_kcat_input_sample20_metadata.csv" \
    --out-rows "$ROOT/data/final/pmak/pmak_kcat_input_sample20_predictions_evaluated.csv" \
    --out-metrics "$ROOT/reports/tables/pmak_sample20_eval_metrics.csv" \
    --out-missing "$ROOT/data/final/pmak/pmak_sample20_invalid_or_unpredicted_rows.csv" \
    --out-missing-summary "$ROOT/reports/tables/pmak_sample20_invalid_or_unpredicted_summary.csv"
