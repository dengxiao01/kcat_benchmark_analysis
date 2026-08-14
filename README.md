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
python3 -m venv .venv
source .venv/bin/activate
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

| File | Rows | Content | Recommended use |
| --- | ---: | --- | --- |
| `data/final/experimental_kcat_truth.csv` | 1,354 | Positive experimental matches linked to model-derived enzyme-reactant records | Truth curation and coverage review |
| `data/final/benchmark_ready_truth.csv` | 1,246 | Final experimental truth and provenance without predictor-only fields | Label and provenance audit |
| `data/final/benchmark_ready_catpred.csv` | 1,246 | Master sequence, substrate SMILES, truth, and provenance table | Preparation of method-specific inputs |
| `data/derived/benchmark_context/` | 6 files | Row-level EC, reaction, and pathway context plus selected summaries | Dataset characterization and pathway analysis |
| `paper/Supplementary_tables.xlsx` | 24 sheets | `Index` and `Table S1` through `Table S23` | Publication tables and record-level audit |

Use `data/final/benchmark_ready_catpred.csv` for a new prediction method. Build
a separate input table containing only features accepted by that method. Join
experimental labels after inference by `entry_id`.

The `catpred` suffix records the first integrated workflow; the 1,246-row file
is method independent. All substrate SMILES pass RDKit parsing, including
Quinate represented by the neutral PubChem CID 6508 structure.

`data/derived/benchmark_context/benchmark_ready_catpred_enriched_context.csv`
adds EC classes, reaction cross-references, KEGG-like groups, direct Yeast9 KEGG
pathways, reaction names, and sequence/SMILES lengths. The companion summaries
cover the benchmark build funnel, EC distribution, reaction distribution, and
pathway distribution.

### Supplementary workbook

| Sheet | Content |
| --- | --- |
| `Table S4` | 1,246-row wide benchmark and prediction matrix |
| `Table S22` | 1,246-row provenance, matching, dependence, measurement, role, direction, and overlap audit |
| `Table S23` | 14,952 method-record rows with method inputs, prediction status, outputs, and errors |

Full figure captions are available in `paper/FIGURE_CAPTIONS.md`. Final figures
are stored in `paper/figures/`. Plotting code and source data are stored in
`paper/figure_sources_0814/`.

`Table S23` contains one row for every method-benchmark combination, including
unscored combinations. Its main fields cover:

- method name, inference regime, applicability, and prediction status;
- benchmark sequence, actual method sequence input, and sequence policy;
- model-encoded reaction equation and method reaction input;
- evaluated metabolite identifier, name, SMILES, stoichiometry, and chemical
  role;
- experimental substrate-support status;
- experimental and predicted kcat in linear and `log10` space;
- absolute and signed prediction errors.

`Table S22` contains source database, match level, reference and label clusters,
measurement dispersion, reaction direction, substrate role, structure mapping,
and training-proximity fields. Continuous training-neighbor values used by
Figure 4 are available in
`paper/figure_sources_0814/Figure4/data/training_proximity_record_audit.csv`.

## Benchmark Composition

| Item | Count |
| --- | ---: |
| Complete sequence and substrate-SMILES records | 1,246 |
| BRENDA substrate-supported records | 778 |
| SABIO-RK-only participant-ambiguous records | 468 |
| *E. coli* records / reactions / genes | 781 / 549 / 350 |
| *S. cerevisiae* records / reactions / genes | 465 / 224 / 168 |
| Model reactions | 773 |
| UniProt accessions | 518 |
| EC strings | 390 |
| Standardized sequence-substrate pairs | 871 |
| Label-assignment clusters | 840 |
| Evaluated methods | 12 |

All 1,246 substrate SMILES pass RDKit parsing. Quinate uses the neutral PubChem
CID 6508 structure.

Experimental kcat values are reported in `s^-1`. Primary evaluation metrics use
`log10(kcat)`.

Experimental matching support includes 459 species + EC + UniProt + substrate
identifier records, 784 species + EC + substrate identifier records, and three
normalized-name records. A total of 480 records inherit a weaker-evidence label
shared across multiple sequences. Use the protein, sequence-substrate pair,
reaction, reference, and label-assignment cluster fields in Table S22 for
dependence-aware analysis.

Chemical roles are assigned after experimental matching with
`configs/currency_cofactor_registry.csv`. The registry combines normalized
names, external identifiers, and standardized parent connectivity. The complete
resource contains 606 currency/cofactor records, 34 carrier-linked variable
records, and 606 other-reactant records. Figure 4 substrate-role comparisons
use the 778 BRENDA substrate-supported records.

Figure 3 reports three comparison designs:

- panel a: strict 1,047-record reaction-common comparison;
- panel b: available-case summaries within the CatPred and KinForm-L scopes,
  with the achieved sample size shown for each method;
- panel c: paired row-bootstrap differences on the strict 1,047-record set.

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

## Data Acquisition and Benchmark Construction

Run all commands from the repository root.

### 1. Candidate enzyme-reaction-metabolite records

| Model | Species | Use |
| --- | --- | --- |
| `eciML1515.json` | *E. coli* | Reactions, GPR rules, EC assignments, metabolites, and cross-references |
| `yeast-GEM.xml` | *S. cerevisiae* | Yeast9 reactions, GPR rules, EC assignments, metabolites, and cross-references |

`src/01_parse_models.py` parses reaction direction, GPR rules, EC numbers,
UniProt identifiers, metabolites, and database cross-references. `or` branches
in a GPR rule are treated as alternative isozymes; `and` branches are retained
as multi-subunit complexes. Candidate records are generated at the species,
reaction, gene group, and reactant level.

All model-encoded reactants enter experimental matching. ATP, NAD(P), water,
protons, CoA-linked compounds, and other common participants are assigned a
chemical role after matching rather than removed by model metabolite ID.

```bash
python src/01_parse_models.py
python src/02_prepare_method_inputs.py
```

### 2. Protein sequences

UniProt accessions from both metabolic models are resolved through the UniProt
REST API and cached locally.

```bash
python src/05_fetch_uniprot_sequences.py
# or
bash scripts/runners/run_sequence_fetch.sh
```

Generated files:

- `data/raw/uniprot_sequences.fasta`
- `data/interim/uniprot_sequences.csv`

### 3. Substrate structures

Structure mapping follows this order:

1. structure annotations in the metabolic model;
2. BiGG, KEGG, ChEBI, and MetaNetX cross-references;
3. local CKB identifier and compound mappings;
4. PubChem PUG REST queries for unresolved identifiers or names;
5. cached success and failure records.

```bash
python src/06_fetch_metanetx_smiles.py
python src/07_fetch_pubchem_smiles.py --ckb-only
python src/07_fetch_pubchem_smiles.py --limit 20 --save-every 10
```

The structure mapper rejects numeric pseudo-SMILES. Quinate is fixed to the
neutral PubChem CID 6508 structure (`C7H12O6`).

CKB is an internal compound-mapping snapshot. Table S22 provides the final
metabolite identifiers, mapping source, source identifier, and SMILES for every
released benchmark record. Rebuilding mappings for the full pre-benchmark
candidate universe requires an equivalent compound cross-reference resource.

### 4. Reaction SMILES and protein structures

TurNuP and PMAK use complete reactant and product structure representations.
Their common scope contains 1,047 records with complete reaction structures.

DEKP uses protein structures. The public-data reconstruction uses the collected
AlphaFold or local structure asset for each of the 1,246 benchmark records.

Reaction structures and protein structures are method assets; they are not
required for methods that use only protein sequence and one substrate SMILES.

### 5. Experimental kcat truth

| Source | Release used | Record handling |
| --- | --- | --- |
| BRENDA | 2026.1 local release | Turnover numbers with substrate-specific evidence; mutant records excluded |
| SABIO-RK | Cache frozen 2026-07-17 | Positive kcat records linked to reaction participants |

Values are normalized to `s^-1`. Numeric ranges are converted to their
midpoint. Matching follows this hierarchy:

1. species + EC + UniProt + substrate database identifier;
2. species + EC + substrate database identifier;
3. species + EC + UniProt + normalized substrate name;
4. species + EC + normalized substrate name.

The highest available match level is retained. Multiple measurements at that
level are aggregated by the median on the original kcat scale and then
transformed to `log10(kcat)`. Median pH and temperature, source database,
reference, and measurement count remain in the provenance fields.

```bash
python src/08_fetch_sabiork_kcat.py --only-ready
python src/10_parse_brenda_kcat.py --match-only
```

Raw BRENDA and SABIO-RK exports are not distributed through GitHub. Obtain and
use them under the source database licenses.

### 6. Final benchmark and substrate roles

```bash
python src/11_finalize_benchmark_data.py
python src/12_finalize_substrate_roles.py
```

`src/11_finalize_benchmark_data.py` creates:

- 1,354 positive experimental matches linked to model candidates;
- 1,246 records with a protein sequence and RDKit-parseable substrate SMILES;
- method-ready tables derived from the 1,246-row master benchmark.

The 1,246-record resource is used for coverage and complete-resource analyses.
The 778-record BRENDA substrate-supported set is used for substrate-supported
sensitivity analyses. The 468 SABIO-RK-only records retain the status
`participant_ambiguous`.

### 7. Release metadata

| File | Content |
| --- | --- |
| `VERSION` | Semantic benchmark version |
| `configs/benchmark_release.json` | Source dates, checksums, row counts, method policies, and statistical settings |
| `configs/matching_rules.yaml` | Experimental matching rules |
| `configs/currency_cofactor_registry.csv` | Currency/cofactor identifiers and structures |
| `configs/methods.yaml` | Method configuration |
| `configs/species.yaml` | Species and model configuration |
| `CHANGELOG.md` | Release changes |

### 8. Public reconstruction policies

| Method | Primary policy | Training rows |
| --- | --- | ---: |
| PreTKcat | Exact sequence-substrate pair overlap excluded | 16,001 |
| DEKP | Exact standardized benchmark pair overlap excluded | 13,171 |

PreTKcat also includes raw-public and joint near-neighbor-excluded sensitivity
models. The three policies fit 16,249, 16,001, and 15,493 training records,
respectively. Exact pair identity uses protein sequence and standardized
uncharged largest-fragment substrate connectivity. Near-neighbor exclusion also
uses sequence identity, bidirectional coverage, and chemical similarity.

## Script Reference

### Benchmark construction

| Script | Output |
| --- | --- |
| `src/01_parse_models.py` | Model reactions, GPRs, ECs, genes, and metabolites |
| `src/02_prepare_method_inputs.py` | Baseline sequence and substrate-SMILES inputs |
| `src/03_match_experimental_kcat.py` | Experimental-truth schema and source status |
| `src/04_build_curation_queues.py` | Unresolved sequence, structure, and truth queues |
| `src/05_fetch_uniprot_sequences.py` | UniProt sequence cache |
| `src/06_fetch_metanetx_smiles.py` | Local cross-reference structure mappings |
| `src/07_fetch_pubchem_smiles.py` | CKB and PubChem structure mappings |
| `src/08_fetch_sabiork_kcat.py` | SABIO-RK cache and matched kcat records |
| `src/10_parse_brenda_kcat.py` | Parsed BRENDA records and matched kcat truth |
| `src/11_finalize_benchmark_data.py` | Three canonical benchmark CSV files |
| `src/12_finalize_substrate_roles.py` | Substrate-role evidence and registry audit |

Convenience runners:

```bash
bash scripts/runners/run_phase1.sh
bash scripts/runners/run_sequence_fetch.sh
bash scripts/runners/run_smiles_fetch.sh
bash scripts/runners/run_sabiork_fetch.sh --only-ready
bash scripts/runners/run_brenda_parse.sh
bash scripts/runners/run_finalize_benchmark_data.sh
```

## Method Workflows

| Method | Preparation, inference, and scoring | Main entry point |
| --- | --- | --- |
| CatPred | `src/12_prepare_catpred_eval.py`, upstream CatPred inference, `src/13_evaluate_catpred_predictions.py` | `bash scripts/runners/run_catpred_predict.sh` or `sbatch scripts/runners/run_catpred_full.sbatch` |
| CataPro | `src/14_prepare_catapro_eval.py`, `src/16_filter_catapro_valid_smiles.py`, CataPro inference, `src/15_evaluate_catapro_predictions.py` | `sbatch scripts/runners/run_catapro_full.sbatch` |
| PMAK | `src/18_prepare_pmak_eval.py`, `src/19_run_pmak_predictions.py`, `src/20_evaluate_pmak_predictions.py` | `sbatch scripts/runners/run_pmak_full.sbatch` |
| KinForm-L | `src/21_prepare_kinform_eval.py`, `src/22_check_kinform_coverage.py`, `src/23_evaluate_kinform_predictions.py` | `sbatch scripts/runners/run_kinform_full.sbatch` |
| KcatNet | `src/24_prepare_kcatnet_eval.py`, `src/25_run_kcatnet_predictions.py`, `src/26_evaluate_kcatnet_predictions.py` | `sbatch scripts/runners/run_kcatnet_full.sbatch` |
| PreTKcat | `src/27_prepare_pretkcat_eval.py`, `src/28_run_pretkcat_predictions.py`, `src/29_evaluate_pretkcat_predictions.py` | `sbatch scripts/runners/run_pretkcat_full.sbatch` |
| DEKP | `src/30_prepare_dekp_eval.py` through `src/34_evaluate_dekp_predictions.py` | `sbatch scripts/runners/run_dekp_public_retrained.sbatch` |
| SELFprot | `src/35_prepare_selfprot_eval.py`, `src/36_run_selfprot_predictions.py`, `src/37_evaluate_selfprot_predictions.py` | `sbatch scripts/runners/run_selfprot_predictions.sbatch` |
| DLKcat | `src/38_run_dlkcat_official.py`, `src/41_evaluate_method_predictions.py` | Run both scripts in the DLKcat environment |
| UniKP | `src/39_run_unikp_official_features.py`, `src/40_predict_unikp_official_py36.py`, `src/41_evaluate_method_predictions.py` | Run feature extraction and prediction in the UniKP environments |
| TurNuP | `src/43_prepare_turnup_eval.py`, `src/45_run_turnup_predictions.py`, `src/41_evaluate_method_predictions.py` | `bash scripts/runners/run_prepare_turnup_eval.sh`, then `sbatch scripts/runners/run_turnup_full.sbatch` |
| GO-HKP | `src/49_prepare_go_hkp_with_yeast_uniprot_go.py`, `src/52_evaluate_go_hkp_predictions.py` | Run both scripts for combined *E. coli* and yeast coverage |

The benchmark environment supports data preparation, evaluation, statistics,
and figure generation. Predictor-specific Python, PyTorch, CUDA, foundation
model, and checkpoint requirements are listed in
`external_methods/METHOD_SOURCES.md`.

Method runners accept environment-specific Python paths. Supported variables
include `CATAPRO_PYTHON`, `PMAK_PYTHON`, `KINFORM_PYTHON`, `KCATNET_PYTHON`,
`PRETKCAT_PYTHON`, `TURNUP_PYTHON`, `DEKP_PYTHON`, and `SELFPROT_PYTHON`.
Set `KCAT_BENCHMARK_ROOT` for Slurm jobs submitted outside the repository root.
Create `logs/` before submitting jobs that write scheduler output there.

## Evaluation and Audit

### Shared evaluation scripts

| Script | Purpose |
| --- | --- |
| `src/17_build_method_eval_summary.py` | Combine method-level metrics |
| `src/41_evaluate_method_predictions.py` | Align predictions, preserve missing reasons, and calculate metrics |
| `src/46_generate_benchmark_report.py` | Generate benchmark comparison tables, figures, and a local report |
| `src/47_generate_dataset_method_context_report.py` | Generate species, EC, GO, KEGG-like, reaction, and method-context analyses |
| `src/50_export_report_tables.py` | Export report tables and their manifest |
| `src/51_rebuild_catpred_overlap_audit.py` | Rebuild CatPred reference-corpus overlap classes |
| `src/52_evaluate_go_hkp_predictions.py` | Evaluate GO-HKP predictions and missing assignments |

Generate local summaries after restoring or producing method predictions:

```bash
python src/17_build_method_eval_summary.py
python src/46_generate_benchmark_report.py
python src/47_generate_dataset_method_context_report.py
python src/50_export_report_tables.py
```

Generated tables, figures, and Markdown reports are written under `reports/`.
The directory is excluded from Git and can be regenerated from the released
method-level outputs.

### Submission audit

```bash
python scripts/build_submission_audits.py
python scripts/rebuild_paper_tables.py
python scripts/recalculate_cluster_inference_v1_2.py
```

Audit outputs are written under `analysis_results/paper_submission_audit/`.
The workflow covers record provenance, shared-label dependence, five cluster
definitions, sensitivity subsets, paired tests, measurement dispersion,
training proximity, reaction direction, substrate roles, and PreTKcat
reconstruction policies.

### Publication files

| Path | Content |
| --- | --- |
| `paper/Supplementary_tables.xlsx` | Formatted `Index` and `Table S1-S23` workbook |
| `paper/figures/` | Four composite figures and 20 standalone panels |
| `paper/FIGURE_CAPTIONS.md` | Full Figure 1-4 captions |
| `paper/figure_sources_0814/` | Plotting code, source data, and figure asset manifest |

```bash
python scripts/generate_manuscript_figures.py
python scripts/generate_manuscript_figures.py \
  --rebuild-code-panels /tmp/kcat_0814_rebuild
```

### Metrics

- MAE: mean absolute error on `log10(kcat)`; MAE = 1 corresponds to an average
  absolute error of one order of magnitude.
- RMSE: root mean squared error on `log10(kcat)`, with greater weight on large
  errors.
- Pearson correlation: linear agreement between predicted and observed values.
- Spearman correlation: rank agreement between predicted and observed values.
- Mean signed error: systematic overprediction or underprediction.
- Within-fold fraction: fraction within 2-fold or 10-fold of the observed kcat.
- Coverage: scored records divided by the defined evaluation scope.

Report the evaluation scope and scored record count with every accuracy metric.

### Release utilities

| Script | Command |
| --- | --- |
| List and restore Zenodo assets | `python scripts/download_zenodo_assets.py --list` |
| Prepare release bundles | `python scripts/prepare_zenodo_release.py --force` |
| Publish a prepared Zenodo release | `python scripts/publish_zenodo_release.py --publish` |

## Repository Layout

### Main files and directories

```text
.
├── configs/
│   ├── benchmark_release.json          Release counts, checksums, dates, and policies
│   ├── currency_cofactor_registry.csv  Substrate-role registry
│   ├── matching_rules.yaml             Experimental matching configuration
│   ├── methods.yaml                    Method configuration
│   └── species.yaml                    Species and model configuration
├── data/final/
│   ├── experimental_kcat_truth.csv     1,354 positive experimental matches
│   ├── benchmark_ready_truth.csv       1,246-row truth and provenance table
│   └── benchmark_ready_catpred.csv     1,246-row sequence/SMILES master table
├── data/derived/benchmark_context/      EC, reaction, and pathway context tables
├── external_methods/
│   └── METHOD_SOURCES.md               Upstream repositories, revisions, and licenses
├── paper/
│   ├── README.md                       Publication artifact guide
│   ├── FIGURE_CAPTIONS.md              Full figure captions
│   ├── Supplementary_tables.xlsx       Index and Table S1-S23
│   ├── figures/                        Final composite and standalone PNG files
│   └── figure_sources_0814/            Plotting code, data, and asset manifest
├── scripts/
│   ├── runners/                        Shell and Slurm method entry points
│   ├── build_submission_audits.py      Record-level and statistical audit
│   ├── rebuild_paper_tables.py         Numerical workbook reconstruction
│   ├── recalculate_cluster_inference_v1_2.py
│   ├── generate_manuscript_figures.py  Figure validation and panel rebuild
│   ├── download_zenodo_assets.py       Asset download, checksum, and restore
│   ├── prepare_zenodo_release.py       Release bundle preparation
│   └── publish_zenodo_release.py       Zenodo publication client
├── src/                                Data construction and method adapters
├── eciML1515.json                      E. coli metabolic model
├── yeast-GEM.xml                       Yeast9 metabolic model
├── requirements.txt                    Benchmark environment dependencies
├── zenodo_assets_manifest.csv          Released asset parts and checksums
├── VERSION
├── CHANGELOG.md
├── CITATION.cff
├── THIRD_PARTY_NOTICES.md
└── LICENSE
```

## Licenses

- Project code: MIT License
- Raw BRENDA and SABIO-RK data: source terms apply
- Third-party code and checkpoints: upstream licenses apply
- Third-party notices: `THIRD_PARTY_NOTICES.md`

Training-overlap classes are available only for methods with accessible
record-level training corpora. `no_joint_neighbor_under_thresholds` means that
no single training record reaches both sequence identity >=80% and chemical
similarity >=0.80; sequence-only or chemical-only overlap may still be present.
Methods without an accessible training corpus are marked `unknown`.

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
