#!/usr/bin/env bash
#SBATCH --job-name=selfprot
#SBATCH --partition=qgpu_3090
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/selfprot_%j.out
#SBATCH --error=logs/selfprot_%j.err

set -euo pipefail

cd /hpcfs/fhome/dengxg/data/kcat_benchmark_analysis
mkdir -p logs

export LD_LIBRARY_PATH=/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}
export PYTHONNOUSERSITE=1
export TRANSFORMERS_OFFLINE=1

/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python src/36_run_selfprot_predictions.py \
  --input data/final/selfprot/selfprot_kcat_input_valid_smiles.csv \
  --output data/final/selfprot/selfprot_kcat_input_output.csv \
  --device cuda:0 \
  --batch-size 8
