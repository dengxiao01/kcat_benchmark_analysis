#!/usr/bin/env bash
#SBATCH -J pretkcat_smoke
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH -o logs/pretkcat_sample_%j.out
#SBATCH -e logs/pretkcat_sample_%j.err

set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${PRETKCAT_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"

export PYTHONNOUSERSITE=1
export TORCH_HOME="${TORCH_HOME:-$ROOT/external_methods/torch_cache}"
export TRANSFORMERS_OFFLINE="${TRANSFORMERS_OFFLINE:-1}"
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="$ROOT/external_methods/kcatnet_scatter_src:$ROOT/external_methods/catapro_pydeps:$ROOT/external_methods/PreTKcat/Pretrained_Model:${PYTHONPATH:-}"
export LD_LIBRARY_PATH="/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}"

cd "$ROOT"
mkdir -p logs
bash run_prepare_pretkcat_eval.sh

"$PYTHON" src/28_run_pretkcat_predictions.py \
    --input "$ROOT/data/final/pretkcat/pretkcat_kcat_input_sample20.csv" \
    --output "$ROOT/data/final/pretkcat/pretkcat_kcat_input_sample20_output.csv" \
    --feature-cache "$ROOT/data/final/pretkcat/pretkcat_feature_cache_sample20.pkl" \
    --model-cache "$ROOT/data/final/pretkcat/pretkcat_extratrees_model_sample20.pkl" \
    --molgnet-model "$ROOT/external_methods/PreTKcat/MolGNet.pt" \
    --prott5-model "$ROOT/external_methods/CataPro/models/prot_t5_xl_uniref50" \
    --device cuda:0 \
    --sample-train-size 500 \
    --n-estimators 50

"$PYTHON" src/29_evaluate_pretkcat_predictions.py \
    --predictions "$ROOT/data/final/pretkcat/pretkcat_kcat_input_sample20_output.csv" \
    --metadata "$ROOT/data/final/pretkcat/pretkcat_kcat_input_sample20_metadata.csv" \
    --out-rows "$ROOT/data/final/pretkcat/pretkcat_kcat_input_sample20_predictions_evaluated.csv" \
    --out-metrics "$ROOT/reports/tables/pretkcat_sample20_eval_metrics.csv" \
    --out-missing "$ROOT/data/final/pretkcat/pretkcat_sample20_invalid_or_unpredicted_rows.csv" \
    --out-missing-summary "$ROOT/reports/tables/pretkcat_sample20_invalid_or_unpredicted_summary.csv"
