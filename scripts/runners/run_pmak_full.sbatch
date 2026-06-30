#!/usr/bin/env bash
#SBATCH -J pmak_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -o logs/pmak_full_%j.out
#SBATCH -e logs/pmak_full_%j.err

set -euo pipefail

ROOT="${KCAT_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$ROOT"
mkdir -p logs

bash "$ROOT/scripts/runners/run_prepare_pmak_eval.sh"
bash "$ROOT/scripts/runners/run_pmak_predict.sh"
