#!/usr/bin/env bash
#SBATCH -J kcatnet_smoke
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -o logs/kcatnet_sample_%j.out
#SBATCH -e logs/kcatnet_sample_%j.err

set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${KCATNET_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"

export PYTHONNOUSERSITE=1
export KCATNET_PROTT5_DIR="${KCATNET_PROTT5_DIR:-$ROOT/external_methods/CataPro/models/prot_t5_xl_uniref50}"
export TORCH_HOME="${TORCH_HOME:-$ROOT/external_methods/torch_cache}"
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="$ROOT/external_methods/kcatnet_scatter_src:$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"
export LD_LIBRARY_PATH="/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}"

cd "$ROOT"
mkdir -p logs
bash run_prepare_kcatnet_eval.sh

"$PYTHON" src/25_run_kcatnet_predictions.py \
    --input "$ROOT/data/final/kcatnet/kcatnet_kcat_input_sample20.csv" \
    --output "$ROOT/data/final/kcatnet/kcatnet_kcat_input_sample20_output.csv" \
    --protein-cache "$ROOT/data/final/kcatnet/kcatnet_protein_cache_sample20.pkl" \
    --ligand-cache "$ROOT/data/final/kcatnet/kcatnet_ligand_cache_sample20.pkl" \
    --device cuda:0 \
    --batch-size 1

"$PYTHON" src/26_evaluate_kcatnet_predictions.py \
    --predictions "$ROOT/data/final/kcatnet/kcatnet_kcat_input_sample20_output.csv" \
    --metadata "$ROOT/data/final/kcatnet/kcatnet_kcat_input_sample20_metadata.csv" \
    --out-rows "$ROOT/data/final/kcatnet/kcatnet_kcat_input_sample20_predictions_evaluated.csv" \
    --out-metrics "$ROOT/reports/tables/kcatnet_sample20_eval_metrics.csv" \
    --out-missing "$ROOT/data/final/kcatnet/kcatnet_sample20_invalid_or_unpredicted_rows.csv" \
    --out-missing-summary "$ROOT/reports/tables/kcatnet_sample20_invalid_or_unpredicted_summary.csv"
