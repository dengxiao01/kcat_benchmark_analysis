#!/usr/bin/env bash
#SBATCH -J kcatnet_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -o logs/kcatnet_full_%j.out
#SBATCH -e logs/kcatnet_full_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

bash run_prepare_kcatnet_eval.sh
bash run_kcatnet_predict.sh
