#!/usr/bin/env bash
#SBATCH -J catapro_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -o logs/catapro_full_%j.out
#SBATCH -e logs/catapro_full_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

bash run_prepare_catapro_eval.sh
bash run_catapro_predict.sh
