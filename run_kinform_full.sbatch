#!/usr/bin/env bash
#SBATCH -J kinform_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -o logs/kinform_full_%j.out
#SBATCH -e logs/kinform_full_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

bash run_prepare_kinform_eval.sh
bash run_kinform_predict.sh
