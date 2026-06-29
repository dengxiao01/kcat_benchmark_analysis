#!/usr/bin/env bash
#SBATCH --job-name=dekp_pub
#SBATCH --partition=qgpu_3090
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/dekp_public_retrained_%j.out
#SBATCH --error=logs/dekp_public_retrained_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

export LD_LIBRARY_PATH=/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}
export PYTHONNOUSERSITE=1
export PYTHONPATH=/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis/external_methods/kcatnet_scatter_src:${PYTHONPATH:-}

/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python src/33_run_dekp_public_retrained.py \
  --epochs 20 \
  --batch-size 32 \
  --patience 5 \
  --device cuda:0
