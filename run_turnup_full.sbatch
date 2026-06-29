#!/usr/bin/env bash
#SBATCH -J turnup_kcat
#SBATCH -p qgpu_4090
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -o logs/turnup_full_%j.out
#SBATCH -e logs/turnup_full_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

bash run_turnup_predict.sh
