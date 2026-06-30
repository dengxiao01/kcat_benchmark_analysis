# Runner Scripts

This directory contains the executable shell and Slurm entry points retained for
the current unified benchmark. Run commands from the repository root.

## Before running

1. Restore the required large assets from Zenodo:
   `python scripts/download_zenodo_assets.py --list`
2. Activate the Conda environment required by the selected method.
3. Create the log directory before submitting a Slurm job: `mkdir -p logs`.

## Script groups

- Benchmark construction: `run_phase1.sh`, `run_sequence_fetch.sh`,
  `run_smiles_fetch.sh`, `run_sabiork_fetch.sh`, `run_brenda_parse.sh`,
  and `run_finalize_benchmark_data.sh`.
- Method input preparation: scripts named `run_prepare_*_eval.sh`.
- Local prediction/evaluation: scripts named `run_*_predict.sh`.
- Full cluster jobs: scripts named `run_*_full.sbatch` plus the DEKP and
  SELFprot Slurm scripts.

## Examples

```bash
bash scripts/runners/run_prepare_catpred_eval.sh
bash scripts/runners/run_catpred_predict.sh
mkdir -p logs
sbatch scripts/runners/run_catpred_full.sbatch
```

Shell runners locate the repository from their own path. Slurm runners use the
submission directory by default; set `KCAT_BENCHMARK_ROOT` when submitting from
another directory. Method-specific Python interpreters can be selected with
variables such as `CATAPRO_PYTHON`, `PMAK_PYTHON`, `DEKP_PYTHON`, and
`SELFPROT_PYTHON`.
