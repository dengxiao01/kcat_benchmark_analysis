# kcat Benchmark Analysis

A reproducible benchmark of published kcat prediction methods on experimentally
supported enzyme-substrate records from *Escherichia coli* and *Saccharomyces
cerevisiae*. The project standardizes experimental truth, protein sequences,
substrate SMILES, method-specific inputs and outputs, and evaluation metrics on
the log10(kcat) scale.

<!-- ZENODO_START -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21024684.svg)](https://doi.org/10.5281/zenodo.21024684)

**Large assets:** [Zenodo record 21024684](https://zenodo.org/records/21024684)
(DOI: `10.5281/zenodo.21024684`). Checksums and restore locations for the
currently supported public assets are listed in `zenodo_assets_manifest.csv`.
<!-- ZENODO_END -->

## Benchmark summary

| Item | Value |
| --- | ---: |
| Experimental truth before model-input filtering | 1,072 rows |
| Unified sequence + substrate SMILES benchmark | 978 rows |
| *E. coli* | 513 rows, 451 reactions, 327 genes |
| *S. cerevisiae* | 465 rows, 224 reactions, 168 genes |
| Experimental sources | BRENDA and SABIO-RK |
| Evaluated methods | 12 |

The comparison is grouped by the information available to each method. Methods
with different input requirements or coverage should not be ranked as if they
were evaluated on identical samples.

| Comparison group | Methods | Evaluation scope |
| --- | --- | --- |
| Sequence + substrate SMILES | DLKcat, UniKP, CataPro, KcatNet, PreTKcat, SELFprot | 977 of 978 rows |
| Reaction-aware | TurNuP, PMAK | 780 rows with complete reaction SMILES |
| Method-specific subset | CatPred, KinForm | Rows accepted by each upstream pipeline |
| Public-data retraining | DEKP-public-retrained | Reproducible retraining, not the paper's private best checkpoint |
| GO functional assignment | GO-HKP | Non-AI baseline, 978 of 978 rows |

### Meaning of the `-official` suffix

A method name ending in `-official` means that the benchmark used the code and
checkpoint distributed by the method authors, rather than a local retraining or
an older overlap-only result. The suffix describes implementation provenance. It
does not mean "published", "accepted", or "endorsed by a journal."

## Canonical benchmark files

| File | Purpose | Use as model input? |
| --- | --- | --- |
| `data/final/experimental_kcat_truth.csv` | All experimental kcat values matched to model-derived enzyme-substrate entries | No; sequence and SMILES completeness is not guaranteed |
| `data/final/benchmark_ready_truth.csv` | The 978-row truth-only evaluation subset | No; it intentionally excludes predictor input fields |
| `data/final/benchmark_ready_catpred.csv` | The 978-row benchmark master table with sequence, SMILES, truth, and provenance | Yes, after extracting only the fields required by a method |

The `catpred` part of the third filename is historical: CatPred was the first
method integrated into the workflow. The file is now the method-independent
benchmark master table.

Method adapters create three logically separate files:

- `*_input.csv` contains only fields visible to the predictor.
- `*_metadata.csv` stores entry identifiers, species, reactions, provenance,
  and processing status.
- `*_truth.csv` is used only after inference for scoring.

This separation prevents experimental labels from leaking into predictor input.

## Data acquisition and benchmark construction

### 1. Candidate enzyme-reaction-substrate entries

The candidate universe is defined from two genome-scale metabolic models:

- `eciML1515.json` for *E. coli*.
- `yeast-GEM.xml` for *S. cerevisiae*.

`src/01_parse_models.py` parses reaction direction, GPR rules, EC numbers,
UniProt identifiers, metabolites, and database cross-references. In a GPR rule,
`or` branches are treated as alternative isozymes, while `and` branches are
retained as multi-subunit enzyme complexes. Candidate rows are generated at the
species + reaction + gene group + substrate level.

Non-cofactor reactants are preferred as candidate substrates. If none are
available, all reactants are retained. This fallback is why ATP, NAD(H), water,
CoA, and other currency-like compounds can still occur in the benchmark.

### 2. Protein sequences

UniProt accessions from the metabolic models are resolved with the UniProt REST
API. Results are cached in:

- `data/raw/uniprot_sequences.fasta`
- `data/interim/uniprot_sequences.csv`

Use `src/05_fetch_uniprot_sequences.py` or
`scripts/runners/run_sequence_fetch.sh`. Cached records are reused on later
runs.

### 3. Substrate structures

SMILES resolution follows a staged strategy:

1. Reuse structure annotations already present in the metabolic model.
2. Resolve BiGG, KEGG, ChEBI, and MetaNetX cross-references through the local CKB
   compound database.
3. Query PubChem PUG REST for unresolved names or registry identifiers.
4. Cache successful mappings and explicit failure reasons.

`src/06_fetch_metanetx_smiles.py` handles local cross-reference enrichment.
`src/07_fetch_pubchem_smiles.py` performs CKB and PubChem resolution. For an
offline-only pass, use:

```bash
python src/07_fetch_pubchem_smiles.py --ckb-only
```

For a small network test before a full query:

```bash
python src/07_fetch_pubchem_smiles.py --limit 20 --save-every 10
```

### 4. Reaction SMILES and protein structures

Complete reactant and product SMILES are prepared only for reaction-aware
methods such as TurNuP and PMAK. They are not required for entry into the common
sequence + substrate SMILES benchmark.

Protein structures are collected only for structure-aware workflows such as the
public-data DEKP retraining. AlphaFold files and local structure archives are
method-specific assets, not universal benchmark requirements.

### 5. Experimental kcat truth

Primary truth comes only from:

- BRENDA turnover-number records.
- SABIO-RK kcat records.

Values must be positive and are normalized to `s^-1`. BRENDA records marked as
mutants, variants, or mutations are excluded by default. Numeric ranges are
converted to their midpoint before aggregation.

The matching hierarchy is:

1. species + EC + UniProt + substrate database ID
2. species + EC + substrate database ID
3. species + EC + UniProt + normalized substrate name
4. species + EC + normalized substrate name

Only the highest available matching level is retained for an entry. Multiple
measurements at that level are aggregated by the median on the original kcat
scale, followed by log10 transformation. Median pH and temperature are retained
when available, together with source database, reference, and measurement count.

Use `src/08_fetch_sabiork_kcat.py` for SABIO-RK and
`src/10_parse_brenda_kcat.py` for BRENDA. Raw database exports are not
stored in Git and must be obtained under their source licenses.

### 6. Final benchmark filter

`src/11_finalize_benchmark_data.py` creates the final products:

- 1,072 matched experimental truth rows.
- 978 rows with a single-protein sequence and usable substrate SMILES.
- Method-ready input tables derived from the 978-row master table.

## Installation and assets

For table and report regeneration:

```bash
git clone https://github.com/dengxiao01/kcat_benchmark_analysis.git
cd kcat_benchmark_analysis
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

List and restore supported Zenodo bundles:

```bash
python scripts/download_zenodo_assets.py --list
python scripts/download_zenodo_assets.py \
  --asset kcat_benchmark_core_data_and_results.tar.gz \
  --restore
```

Third-party source locations, evaluated revisions, and upstream licenses are
listed in `external_methods/METHOD_SOURCES.md`. Deep-learning methods use
different Python, PyTorch, CUDA, and foundation-model versions. Activate the
method's upstream environment before running its predictor; do not force every
method into one environment.

## Script usage

Run commands from the repository root. Scripts with an argument parser expose
full options through:

```bash
python path/to/script.py -h
```

### Benchmark construction scripts

| Script | Purpose | Typical command |
| --- | --- | --- |
| `src/01_parse_models.py` | Parse GEM reactions, GPRs, ECs, UniProt IDs, and substrates | `python src/01_parse_models.py` |
| `src/02_prepare_method_inputs.py` | Build baseline sequence/SMILES input queues | `python src/02_prepare_method_inputs.py` |
| `src/03_match_experimental_kcat.py` | Maintain the experimental-truth schema and source-status table | `python src/03_match_experimental_kcat.py` |
| `src/04_build_curation_queues.py` | Build unresolved sequence, SMILES, and truth review queues | `python src/04_build_curation_queues.py` |
| `src/05_fetch_uniprot_sequences.py` | Fetch and cache UniProt sequences | `python src/05_fetch_uniprot_sequences.py` |
| `src/06_fetch_metanetx_smiles.py` | Resolve structures from local cross-references | `python src/06_fetch_metanetx_smiles.py` |
| `src/07_fetch_pubchem_smiles.py` | Resolve missing structures with CKB and PubChem | `python src/07_fetch_pubchem_smiles.py --limit 20` |
| `src/08_fetch_sabiork_kcat.py` | Query/cache SABIO-RK and match kcat records | `python src/08_fetch_sabiork_kcat.py --only-ready` |
| `src/10_parse_brenda_kcat.py` | Parse BRENDA JSON and match turnover numbers | `python src/10_parse_brenda_kcat.py` |
| `src/11_finalize_benchmark_data.py` | Build the canonical truth and benchmark files | `python src/11_finalize_benchmark_data.py` |

Convenience runners for this stage:

```bash
bash scripts/runners/run_phase1.sh
bash scripts/runners/run_sequence_fetch.sh
bash scripts/runners/run_smiles_fetch.sh
bash scripts/runners/run_sabiork_fetch.sh --only-ready
bash scripts/runners/run_brenda_parse.sh
bash scripts/runners/run_finalize_benchmark_data.sh
```

### Method workflows

| Method | Preparation, inference, and scoring scripts | Convenience entry point |
| --- | --- | --- |
| CatPred | `src/12_prepare_catpred_eval.py`, `src/13_evaluate_catpred_predictions.py`; upstream inference is wrapped by the runner | `bash scripts/runners/run_catpred_predict.sh` or `sbatch scripts/runners/run_catpred_full.sbatch` |
| CataPro | `src/14_prepare_catapro_eval.py`, `src/16_filter_catapro_valid_smiles.py`, `src/15_evaluate_catapro_predictions.py` | `sbatch scripts/runners/run_catapro_full.sbatch` |
| PMAK | `src/18_prepare_pmak_eval.py`, `src/19_run_pmak_predictions.py`, `src/20_evaluate_pmak_predictions.py` | `sbatch scripts/runners/run_pmak_full.sbatch` |
| KinForm | `src/21_prepare_kinform_eval.py`, `src/22_check_kinform_coverage.py`, `src/23_evaluate_kinform_predictions.py` | `sbatch scripts/runners/run_kinform_full.sbatch` |
| KcatNet | `src/24_prepare_kcatnet_eval.py`, `src/25_run_kcatnet_predictions.py`, `src/26_evaluate_kcatnet_predictions.py` | `sbatch scripts/runners/run_kcatnet_full.sbatch` |
| PreTKcat | `src/27_prepare_pretkcat_eval.py`, `src/28_run_pretkcat_predictions.py`, `src/29_evaluate_pretkcat_predictions.py` | `sbatch scripts/runners/run_pretkcat_full.sbatch` |
| DEKP-public-retrained | `src/30_prepare_dekp_eval.py`, `src/31_download_dekp_missing_structures.py`, `src/32_collect_dekp_structures.py`, `src/33_run_dekp_public_retrained.py`, `src/34_evaluate_dekp_predictions.py` | `sbatch scripts/runners/run_dekp_public_retrained.sbatch` |
| SELFprot | `src/35_prepare_selfprot_eval.py`, `src/36_run_selfprot_predictions.py`, `src/37_evaluate_selfprot_predictions.py` | `sbatch scripts/runners/run_selfprot_predictions.sbatch` |
| DLKcat | `src/38_run_dlkcat_official.py` followed by `src/41_evaluate_method_predictions.py` | Run the two Python scripts in the DLKcat environment |
| UniKP | `src/39_run_unikp_official_features.py`, `src/40_predict_unikp_official_py36.py`, then `src/41_evaluate_method_predictions.py` | Run feature extraction and prediction in their required environments |
| TurNuP | `src/43_prepare_turnup_eval.py`, `src/45_run_turnup_predictions.py`, then `src/41_evaluate_method_predictions.py` | `bash scripts/runners/run_prepare_turnup_eval.sh`, then `sbatch scripts/runners/run_turnup_full.sbatch` |
| GO-HKP | `src/48_prepare_go_hkp_eval.py` or `src/49_prepare_go_hkp_with_yeast_uniprot_go.py`, then `src/41_evaluate_method_predictions.py` | Run `src/49_*` for combined *E. coli* and yeast coverage |

Most method runners accept environment variables for paths, devices, and Python
interpreters. Examples include `CATAPRO_PYTHON`, `PMAK_PYTHON`,
`KINFORM_PYTHON`, `KCATNET_PYTHON`, `PRETKCAT_PYTHON`,
`TURNUP_PYTHON`, `DEKP_PYTHON`, and `SELFPROT_PYTHON`.

Shell runners locate the repository relative to their own path. Slurm runners
use `SLURM_SUBMIT_DIR` by default; set `KCAT_BENCHMARK_ROOT` when
submitting from another directory. Create `logs/` before `sbatch` because
Slurm resolves output paths at submission time.

### Shared evaluation and output scripts

| Script | Purpose | Command |
| --- | --- | --- |
| `src/17_build_method_eval_summary.py` | Combine the current 12 method metric files | `python src/17_build_method_eval_summary.py` |
| `src/41_evaluate_method_predictions.py` | Generic row alignment, missing-reason reporting, and metric calculation | `python src/41_evaluate_method_predictions.py -h` |
| `src/46_generate_benchmark_report.py` | Rebuild comparison tables, figures, and a local Markdown report | `python src/46_generate_benchmark_report.py` |
| `src/47_generate_dataset_method_context_report.py` | Rebuild dataset context, GO/KEGG tables, figures, and a local Markdown report | `python src/47_generate_dataset_method_context_report.py` |
| `src/50_export_report_tables.py` | Rebuild the standalone table export directory and manifest | `python src/50_export_report_tables.py` |

The local Markdown reports and chronological work notes are intentionally not
tracked in Git. Public tables and figures are stored under `reports/tables/`,
`reports/report_tables/`, and `reports/figures/`.

To rebuild the public summaries after method predictions are available:

```bash
python src/17_build_method_eval_summary.py
python src/46_generate_benchmark_report.py
python src/47_generate_dataset_method_context_report.py
python src/50_export_report_tables.py
```

### Release and Zenodo utilities

| Script | Purpose | Example |
| --- | --- | --- |
| `scripts/download_zenodo_assets.py` | List, download, verify, join, and restore public bundles | `python scripts/download_zenodo_assets.py --list` |
| `scripts/prepare_zenodo_release.py` | Build whitelisted public result/model bundles and checksums | `python scripts/prepare_zenodo_release.py --force` |
| `scripts/publish_zenodo_release.py` | Upload a prepared release using the local ignored token file | `python scripts/publish_zenodo_release.py --publish` |

## Evaluation metrics

All primary metrics are computed on log10(kcat):

- MAE: average absolute error. MAE = 1 means an average error of roughly one
  order of magnitude.
- RMSE: gives more weight to large errors.
- Pearson correlation: linear agreement.
- Spearman correlation: rank agreement.
- Bias: systematic overprediction or underprediction.
- Within-10-fold fraction: fraction with absolute log10 error <= 1.
- Coverage: number of scored rows divided by 978.

Always report coverage next to accuracy. A lower error on a restricted subset is
not automatically better than a slightly higher error on the near-full
benchmark.

## Repository layout

```text
configs/             Species and matching rules
data/final/          Three Git-tracked benchmark CSVs; method outputs are local
external_methods/    Only METHOD_SOURCES.md is tracked; source trees are local
reports/tables/      Public machine-readable summaries
reports/report_tables/ Standalone tables grouped by report purpose
reports/figures/     Public benchmark and dataset figures
scripts/runners/     Local shell and Slurm entry points
scripts/             Zenodo preparation, publishing, and download utilities
src/                 Benchmark construction, method adapters, and evaluation
release/             Local ignored Zenodo staging
logs/                Local ignored scheduler and run logs
```

Git stores code, configuration, the three canonical benchmark files, public
figures, and compact result tables. Raw databases, intermediate caches,
third-party repositories, model weights, structures, and large row-level method
outputs are excluded from Git.

## Data and licensing notes

- BRENDA and SABIO-RK remain the primary experimental sources.
- Raw databases must be obtained and used under their source licenses.
- Early inferred or database-filled kcat values are not part of benchmark truth.
- Training-set overlap can make performance optimistic. Publication analyses
  should report sequence, SMILES, and sequence-SMILES-pair overlap.
- The original project code is released under the MIT License in `LICENSE`.
- Third-party source code, checkpoints, and foundation models retain their
  upstream licenses. See `THIRD_PARTY_NOTICES.md` and
  `external_methods/METHOD_SOURCES.md`.

## Citation

The repository DOI is [10.5281/zenodo.21024684](https://doi.org/10.5281/zenodo.21024684).
When using the benchmark, also cite the original prediction methods and the
BRENDA, SABIO-RK, UniProt, PubChem, and metabolic-model resources used for the
relevant analysis.
