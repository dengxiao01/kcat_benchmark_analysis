#!/usr/bin/env bash
#SBATCH -J mtlkp_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -o logs/mtlkp_full_%j.out
#SBATCH -e logs/mtlkp_full_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

bash run_mtlkp_predict.sh
