#!/usr/bin/env bash
#SBATCH --job-name=dekp_pub
#SBATCH --partition=qgpu_3090
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/dekp_public_retrained_%j.out
#SBATCH --error=logs/dekp_public_retrained_%j.err

set -euo pipefail

ROOT="${KCAT_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$ROOT"
mkdir -p logs

PYTHON="${DEKP_PYTHON:-python}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/external_methods/kcatnet_scatter_src:${PYTHONPATH:-}"

"$PYTHON" src/33_run_dekp_public_retrained.py \
  --epochs 20 \
  --batch-size 32 \
  --patience 5 \
  --device cuda:0
