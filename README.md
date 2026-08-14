# kcat Benchmark Analysis

Benchmark data, evaluation code, supplementary tables, and figures for 12
published kcat prediction methods on *Escherichia coli* and *Saccharomyces
cerevisiae* enzyme-reaction-metabolite records.

| Release item | Value |
| --- | --- |
| Version | `v1.2.0` |
| Data revision | `1.2.0-r3` |
| Manuscript snapshot | `0814/V4` |
| Benchmark size | 1,246 records |
| Data freeze | 2026-07-17 |

<!-- ZENODO_START -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21859994.svg)](https://doi.org/10.5281/zenodo.21859994)

- Version DOI: [`10.5281/zenodo.21859994`](https://doi.org/10.5281/zenodo.21859994)
- Concept DOI: [`10.5281/zenodo.21024683`](https://doi.org/10.5281/zenodo.21024683)
- Asset manifest: [`zenodo_assets_manifest.csv`](zenodo_assets_manifest.csv)
<!-- ZENODO_END -->

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/dxg-9527/kcat_benchmark_analysis.git
cd kcat_benchmark_analysis
```

### 2. Install the benchmark environment

```bash
conda create -n kcat-benchmark python=3.11 -y
conda activate kcat-benchmark
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### 3. Restore the released assets

```bash
python scripts/download_zenodo_assets.py --list
python scripts/download_zenodo_assets.py \
  --asset kcat_benchmark_core_data_and_results.tar.gz \
  --restore
```

### 4. Validate the publication files

```bash
python scripts/generate_manuscript_figures.py
```

Expected output:

```text
0814/V4 snapshot OK: 4 composites, 20 panels, 25 checks
```

### 5. Rebuild the code-backed figure panels

```bash
python scripts/generate_manuscript_figures.py \
  --rebuild-code-panels /tmp/kcat_0814_rebuild
```

## Benchmark Files

| File | Rows | Content |
| --- | ---: | --- |
| `data/final/experimental_kcat_truth.csv` | 1,354 | Positive experimental matches before sequence and structure filtering |
| `data/final/benchmark_ready_truth.csv` | 1,246 | Experimental truth and provenance |
| `data/final/benchmark_ready_catpred.csv` | 1,246 | Master sequence, substrate SMILES, truth, and provenance table |
| `paper/Supplementary_tables.xlsx` | 24 sheets | `Index` and `Table S1` through `Table S23` |

Use `data/final/benchmark_ready_catpred.csv` to prepare model inputs. Select
only the columns required by the evaluated method.

### Supplementary workbook

| Sheet | Content |
| --- | --- |
| `Table S4` | 1,246-row wide benchmark and prediction matrix |
| `Table S22` | 1,246-row provenance, matching, dependence, measurement, role, direction, and overlap audit |
| `Table S23` | 14,952 method-record rows with method inputs, prediction status, outputs, and errors |

Full figure captions are available in `paper/FIGURE_CAPTIONS.md`. Final figures
are stored in `paper/figures/`. Plotting code and source data are stored in
`paper/figure_sources_0814/`.

## Benchmark Composition

| Item | Count |
| --- | ---: |
| Complete sequence and substrate-SMILES records | 1,246 |
| BRENDA substrate-supported records | 778 |
| SABIO-RK-only participant-ambiguous records | 468 |
| *E. coli* records | 781 |
| *S. cerevisiae* records | 465 |
| Model reactions | 773 |
| UniProt accessions | 518 |
| EC strings | 390 |
| Standardized sequence-substrate pairs | 871 |
| Label-assignment clusters | 840 |

All 1,246 substrate SMILES pass RDKit parsing. Quinate uses the neutral PubChem
CID 6508 structure.

Experimental kcat values are reported in `s^-1`. Primary evaluation metrics use
`log10(kcat)`.

## Evaluated Methods

| Scope | Methods | Scored records |
| --- | --- | ---: |
| Sequence and substrate checkpoint | DLKcat, UniKP, CataPro, KcatNet, SELFprot | 1,246 each |
| Public reconstruction | PreTKcat | 1,246 |
| Reaction-aware checkpoint | TurNuP, PMAK | 1,047 each |
| Method-specific checkpoint subset | CatPred | 1,156 |
| Method-specific checkpoint subset | KinForm-L | 729 |
| Structure-aware public reconstruction | DEKP | 1,246 |
| GO functional assignment | GO-HKP | 1,236 |

PreTKcat uses the exact-pair-excluded reconstruction as its primary result.
DEKP uses the public-data retrained model. Method sources, revisions, licenses,
and checkpoints are listed in `external_methods/METHOD_SOURCES.md`.

## Input and Output Convention

Method adapters use three input layers:

| Suffix | Content |
| --- | --- |
| `*_input.csv` | Predictor-visible features |
| `*_metadata.csv` | Entry identifiers, reactions, provenance, and processing status |
| `*_truth.csv` | Experimental labels used after inference |

Evaluated outputs are stored as `data/final/<method>/*_predictions_evaluated.csv`
after restoring the Zenodo core bundle or running the corresponding method.

## Benchmark Construction

Run commands from the repository root.

```bash
python src/01_parse_models.py
python src/02_prepare_method_inputs.py
python src/05_fetch_uniprot_sequences.py
python src/06_fetch_metanetx_smiles.py
python src/07_fetch_pubchem_smiles.py --ckb-only
python src/08_fetch_sabiork_kcat.py --only-ready
python src/10_parse_brenda_kcat.py --match-only
python src/11_finalize_benchmark_data.py
python src/12_finalize_substrate_roles.py
```

### Data sources

| Source | Use |
| --- | --- |
| `eciML1515.json` | *E. coli* reactions, genes, EC assignments, and metabolites |
| `yeast-GEM.xml` | Yeast9 reactions, genes, EC assignments, and metabolites |
| UniProt | Protein sequences and identifiers |
| BRENDA 2026.1 | Substrate-supported turnover numbers |
| SABIO-RK | Cached kcat records |
| MetaNetX, PubChem, ChEBI, KEGG, BiGG | Compound identifiers and structures |

Source checksums, release dates, query settings, and matching policies are in
`configs/benchmark_release.json`.

### Experimental matching hierarchy

1. Species + EC + UniProt + substrate database identifier
2. Species + EC + substrate database identifier
3. Species + EC + UniProt + normalized substrate name
4. Species + EC + normalized substrate name

Positive measurements at the highest available level are aggregated by the
median on the original kcat scale and then transformed to `log10(kcat)`.

## Method Workflows

| Method | Main entry point |
| --- | --- |
| CatPred | `bash scripts/runners/run_catpred_predict.sh` |
| CataPro | `sbatch scripts/runners/run_catapro_full.sbatch` |
| PMAK | `sbatch scripts/runners/run_pmak_full.sbatch` |
| KinForm-L | `sbatch scripts/runners/run_kinform_full.sbatch` |
| KcatNet | `sbatch scripts/runners/run_kcatnet_full.sbatch` |
| PreTKcat | `sbatch scripts/runners/run_pretkcat_full.sbatch` |
| DEKP | `sbatch scripts/runners/run_dekp_public_retrained.sbatch` |
| SELFprot | `sbatch scripts/runners/run_selfprot_predictions.sbatch` |
| TurNuP | `sbatch scripts/runners/run_turnup_full.sbatch` |

DLKcat, UniKP, and GO-HKP use their corresponding scripts under `src/`.
Deep-learning predictors require the upstream environments listed in
`external_methods/METHOD_SOURCES.md`.

## Evaluation and Audit

```bash
python src/17_build_method_eval_summary.py
python src/41_evaluate_method_predictions.py -h
python scripts/build_submission_audits.py
python scripts/rebuild_paper_tables.py
python scripts/recalculate_cluster_inference_v1_2.py
```

Generated tables, figures, and Markdown reports are written under `reports/`
and `analysis_results/`. Both directories are excluded from Git.

### Metrics

- MAE on `log10(kcat)`
- RMSE on `log10(kcat)`
- Pearson correlation
- Spearman correlation
- Mean signed error
- Fraction within 2-fold and 10-fold
- Prediction coverage

Compare accuracy together with coverage and evaluation scope.

## Repository Layout

```text
configs/             Release metadata and matching rules
data/final/          Canonical benchmark files and restored method outputs
external_methods/    Method source manifest and restored upstream projects
paper/               Supplementary workbook, captions, figures, and figure sources
scripts/             Audit, release, download, shell, and Slurm entry points
src/                 Benchmark construction, method adapters, and evaluation code
analysis_results/    Generated numerical audits
reports/             Generated tables, figures, and local reports
release/             Generated Zenodo staging files
```

## Licenses

- Project code: MIT License
- Raw BRENDA and SABIO-RK data: source terms apply
- Third-party code and checkpoints: upstream licenses apply
- Third-party notices: `THIRD_PARTY_NOTICES.md`

## Citation

Cite the benchmark version DOI:

```text
10.5281/zenodo.21859994
```

Use the concept DOI for the latest release:

```text
10.5281/zenodo.21024683
```

Citation metadata are available in `CITATION.cff`.
