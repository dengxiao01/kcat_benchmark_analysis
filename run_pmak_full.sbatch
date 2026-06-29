#!/usr/bin/env bash
#SBATCH -J pmak_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -o logs/pmak_full_%j.out
#SBATCH -e logs/pmak_full_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

bash run_prepare_pmak_eval.sh
bash run_pmak_predict.sh
