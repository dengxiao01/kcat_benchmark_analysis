#!/usr/bin/env bash
set -euo pipefail

ROOT="/hpcfs/fhome/dengxg/data/kcat_benchmark_analysis"
PYTHON="${KCATNET_PYTHON:-/hpcfs/fhome/dengxg/.conda/envs/enyrnx/bin/python}"
INPUT="${KCATNET_INPUT:-$ROOT/data/final/kcatnet/kcatnet_kcat_input_valid_smiles.csv}"
OUTPUT="${KCATNET_OUTPUT:-$ROOT/data/final/kcatnet/kcatnet_kcat_input_output.csv}"
METADATA="${KCATNET_METADATA:-$ROOT/data/final/kcatnet/kcatnet_kcat_input_metadata.csv}"
PROTEIN_CACHE="${KCATNET_PROTEIN_CACHE:-$ROOT/data/final/kcatnet/kcatnet_protein_cache.pkl}"
LIGAND_CACHE="${KCATNET_LIGAND_CACHE:-$ROOT/data/final/kcatnet/kcatnet_ligand_cache.pkl}"
DEVICE="${KCATNET_DEVICE:-cuda:0}"

export PYTHONNOUSERSITE=1
export KCATNET_PROTT5_DIR="${KCATNET_PROTT5_DIR:-$ROOT/external_methods/CataPro/models/prot_t5_xl_uniref50}"
export TORCH_HOME="${TORCH_HOME:-$ROOT/external_methods/torch_cache}"
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="$ROOT/external_methods/kcatnet_scatter_src:$ROOT/external_methods/catapro_pydeps:${PYTHONPATH:-}"
export LD_LIBRARY_PATH="/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib:${LD_LIBRARY_PATH:-}"

cd "$ROOT"
"$PYTHON" src/25_run_kcatnet_predictions.py \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --protein-cache "$PROTEIN_CACHE" \
    --ligand-cache "$LIGAND_CACHE" \
    --device "$DEVICE" \
    --batch-size "${KCATNET_BATCH_SIZE:-1}"

"$PYTHON" src/26_evaluate_kcatnet_predictions.py \
    --predictions "$OUTPUT" \
    --metadata "$METADATA"

"$PYTHON" src/17_build_method_eval_summary.py
