#!/usr/bin/env bash
#SBATCH -J catpred_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -o logs/catpred_full_%j.out
#SBATCH -e logs/catpred_full_%j.err

set -euo pipefail

ROOT="${KCAT_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$ROOT"
mkdir -p logs

bash "$ROOT/scripts/runners/run_catpred_predict.sh"
