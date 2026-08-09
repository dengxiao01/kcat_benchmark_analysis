# kcat Benchmark Analysis

A reproducible benchmark of published kcat prediction methods on experimentally
supported enzyme-substrate records from *Escherichia coli* and *Saccharomyces
cerevisiae*. The project standardizes experimental truth, protein sequences,
substrate SMILES, method-specific inputs and outputs, and evaluation metrics on
the log10(kcat) scale.

**Current release:** benchmark `v1.2.0`, released `2026-08-09`; source-data
freeze `2026-07-17`; artifact revision `1.2.0-r3` dated `2026-08-07`; table schema `1.2`.
The canonical resource contains 1,246 rows. Versioned release metadata are in
`configs/benchmark_release.json`, and human-readable changes are in
`CHANGELOG.md`.

<!-- ZENODO_START -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21859994.svg)](https://doi.org/10.5281/zenodo.21859994)

**Large assets:** [Zenodo record 21859994](https://zenodo.org/records/21859994) (version DOI: `10.5281/zenodo.21859994`). File checksums and restore paths are in `zenodo_assets_manifest.csv`.
All releases are linked by the concept DOI [`10.5281/zenodo.21024683`](https://doi.org/10.5281/zenodo.21024683).
<!-- ZENODO_END -->

## Benchmark summary

| Item | Value |
| --- | ---: |
| Positive experimental matches before sequence/structure filtering | 1,354 rows |
| Complete-resource sequence + substrate SMILES benchmark | 1,246 rows |
| BRENDA substrate-supported sensitivity set | 778 rows |
| SABIO-RK-only participant-ambiguous records | 468 rows |
| *E. coli* | 781 rows, 549 reactions, 350 genes |
| *S. cerevisiae* | 465 rows, 224 reactions, 168 genes |
| Unique model reactions / UniProt accessions / EC strings | 773 / 518 / 390 |
| Evaluated methods | 12 |

The benchmark retains all model-encoded reactants before experimental matching.
Candidate generation, experimental substrate support, and post hoc chemical-role
annotation are separate fields. BRENDA provides a substrate-specific field;
SABIO-RK-only cached records identify reaction participants but do not establish
which participant was varied in the assay.

Matching support consists of 459 species + EC + UniProt + substrate-ID rows,
784 species + EC + substrate-ID rows, and three normalized-name rows. There are
871 unique standardized sequence-substrate pairs and 840 label-assignment
clusters; 480 rows inherit a weaker-evidence label shared across multiple
sequences. Benchmark rows must therefore not be interpreted as independent
experiments. These fields are in `paper/tables_v1.2.0/S16_Label_audit.csv` and
`Record_audit.csv`.

Chemical roles are assigned after experimental matching with the versioned
registry in `configs/currency_cofactor_registry.csv`. It combines normalized
names, external database identifiers, and standardized parent connectivity and
does not rely on yeast-GEM `s_XXXX` identifiers. Across the complete resource,
606 rows are currency/cofactor, 34 carrier-linked variable, and 606 other
reactants. Role analyses in Figure 4 are restricted to the 778
substrate-supported rows.

Methods are grouped by inference-time information and implementation status:

| Comparison group | Methods | Evaluation scope |
| --- | --- | --- |
| Released sequence + substrate checkpoint | DLKcat, UniKP, CataPro, KcatNet, SELFprot | 1,246 of 1,246 rows |
| Temperature-conditioned public reconstruction | PreTKcat | 1,246 rows; exact-overlap-excluded primary model |
| Reaction-aware checkpoint | TurNuP, PMAK | 1,047 rows with complete reactant and product structures |
| Method-specific checkpoint subset | CatPred, KinForm-L | 1,156 and 729 rows |
| Structure-aware public reconstruction | DEKP | 1,246 rows; 1,246 valid structures |
| GO functional assignment | GO-HKP | 1,236 rows |

Figure 3a is an available-case summary and prints a separate `n` in every cell.
Panels 3b/c use the strict 1,047-record intersection for KcatNet, TurNuP, and
PMAK. Methods with different inputs or row coverage must not be pooled into one
leaderboard.

## Canonical benchmark files

| File | Purpose | Use as model input? |
| --- | --- | --- |
| `data/final/experimental_kcat_truth.csv` | 1,354 positive experimental matches linked to model-derived enzyme-reactant rows | No; sequence and SMILES completeness is not guaranteed |
| `data/final/benchmark_ready_truth.csv` | 1,246-row truth and provenance table | No; it intentionally excludes predictor-only fields |
| `data/final/benchmark_ready_catpred.csv` | 1,246-row master table with sequence, SMILES, truth, and provenance | Yes, after extracting only fields required by the method |

The `catpred` suffix in the master filename is historical: CatPred was the first
method integrated. The file is method independent. All 1,246 SMILES are
RDKit-parseable, including Quinate mapped to neutral PubChem CID 6508.

### Supplementary Data Table 0

`paper/tables_v1.2.0/Table0.csv` is the primary method-record table. It contains
14,952 rows (12 methods x 1,246 benchmark records), including explicit unscored
combinations. Important fields include:

- `method`, `inference_regime`, and `prediction_status`;
- `benchmark_sequence`, `method_sequence_input`, and `sequence_input_policy`;
- source-model forward reaction equations and whether complete reaction context
  is a direct method input;
- evaluated metabolite ID, name, SMILES, stoichiometry, chemical role, and
  experimental substrate-support status;
- experimental and predicted kcat in linear and log10 space.

`paper/tables_v1.2.0/Table0_wide.csv` is the companion 1,246-row matrix with one
prediction column per method. `paper/tables_v1.2.0/Record_audit.csv` contains
1,246 rows and 132 provenance, dependence, measurement, mapping, role, direction,
and training-proximity fields used to reproduce S16-S24.

Method adapters create three logically separate files:

- `*_input.csv` contains fields visible to the predictor;
- `*_metadata.csv` stores identifiers, reactions, provenance, and status;
- `*_truth.csv` is joined only after inference for scoring.

This separation prevents experimental labels from entering predictor input.

## Data acquisition and benchmark construction

### 1. Candidate enzyme-reaction-substrate entries

The candidate universe is defined from two genome-scale metabolic models:

- `eciML1515.json` for *E. coli* (`iML1515`, model version 1).
- `yeast-GEM.xml` for *S. cerevisiae* (`yeastGEM_v9.0.2`).

The exact model-file SHA256 checksums are stored in
`configs/benchmark_release.json`. They identify the local model snapshots used
to generate candidate records.

The same release file records the BRENDA archive and parsed-table checksums, the
SABIO-RK API/query template and cache checksums, and the local CKB compound
repository commit/database checksum. CKB is an internal snapshot for which no
public release was identified; it is documented for provenance and is not
presented as independently downloadable public input. `Record_audit.csv`
materializes each released row's model metabolite identifiers, final SMILES,
mapping source, and source identifier, so users can reconstruct the 1,246-row
benchmark evaluation without CKB. It does not reproduce mapping for the full
pre-benchmark candidate universe.

`src/01_parse_models.py` parses reaction direction, GPR rules, EC numbers,
UniProt identifiers, metabolites, and database cross-references. In a GPR rule,
`or` branches are treated as alternative isozymes, while `and` branches are
retained as multi-subunit enzyme complexes. Candidate rows are generated at the
species + reaction + gene group + substrate level.

Every model-encoded reactant is retained before experimental substrate
matching. ATP, NAD(P), water, protons, CoA-linked compounds, and other currency
participants are not removed by model ID. Experimental records select the
evaluated metabolite where substrate-specific evidence exists; chemical role is
annotated only after matching. Reaction side and stoichiometry are retained for
audit but are not treated as proof of the assay-variable substrate.

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

The mapper rejects pure numeric values before they can be treated as SMILES.
Benchmark v1.2.0 uses the neutral PubChem CID 6508 Quinate structure (`C7H12O6`)
and rejects numeric pseudo-SMILES before generating method inputs.

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

- BRENDA turnover-number records (local release 2026.1).
- SABIO-RK kcat records (cached snapshot frozen on 2026-07-17).

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
`src/10_parse_brenda_kcat.py` for BRENDA. After structures change, rebuild
truth from already parsed local records without reparsing the BRENDA archive:

```bash
python src/10_parse_brenda_kcat.py --match-only
python src/11_finalize_benchmark_data.py
```

Raw database exports are not stored in Git and must be obtained under their
source licenses.

### 6. Final benchmark filter

`src/11_finalize_benchmark_data.py` creates the final products:

- 1,354 positive experimental matches linked to model candidates;
- 1,246 rows with one protein sequence and an RDKit-parseable substrate SMILES;
- method-ready tables derived from the 1,246-row master table.

The complete 1,246-row resource is primary for coverage and traceability. The
778-row BRENDA substrate-supported set is the formal sensitivity analysis; the
468 SABIO-RK-only participant-ambiguous rows are not claimed to have a verified
assay-variable substrate.

## Release metadata and reconstruction policies

Release identity is recorded in three places:

- `VERSION`: short semantic version (`1.2.0`);
- `configs/benchmark_release.json`: dates, row counts, SHA256 values, source
  snapshots, method policies, and statistical settings;
- `CHANGELOG.md`: human-readable history.

PreTKcat is reported as a public reconstruction with three policies in S24:
raw-public fits 16,249 rows; the primary exact-excluded model removes 248 rows
and fits 16,001; near-excluded removes 756 rows in total and fits 15,493. DEKP
removes 230 exact-overlap training rows and fits 13,171 public rows. Exact pair
identity uses model sequence plus uncharged largest-fragment connectivity;
near exclusion additionally requires one training row to meet sequence identity,
bidirectional coverage, and chemical-similarity thresholds. Original SMILES
remain method features; standardized parent identity is used for audit and
exclusion.

## Installation and assets

For table and report regeneration:

```bash
git clone https://github.com/dxg-9527/kcat_benchmark_analysis.git
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
| `src/10_parse_brenda_kcat.py` | Parse BRENDA or rematch cached records after mapping changes | `python src/10_parse_brenda_kcat.py --match-only` |
| `src/11_finalize_benchmark_data.py` | Build the canonical truth and benchmark files | `python src/11_finalize_benchmark_data.py` |
| `src/12_finalize_substrate_roles.py` | Apply the auditable chemical-role registry and role evidence | `python src/12_finalize_substrate_roles.py` |

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
| KinForm-L | `src/21_prepare_kinform_eval.py`, `src/22_check_kinform_coverage.py`, `src/23_evaluate_kinform_predictions.py` | `sbatch scripts/runners/run_kinform_full.sbatch` |
| KcatNet | `src/24_prepare_kcatnet_eval.py`, `src/25_run_kcatnet_predictions.py`, `src/26_evaluate_kcatnet_predictions.py` | `sbatch scripts/runners/run_kcatnet_full.sbatch` |
| PreTKcat | `src/27_prepare_pretkcat_eval.py`, `src/28_run_pretkcat_predictions.py`, `src/29_evaluate_pretkcat_predictions.py` | `sbatch scripts/runners/run_pretkcat_full.sbatch` |
| DEKP-public-retrained | `src/30_prepare_dekp_eval.py`, `src/31_download_dekp_missing_structures.py`, `src/32_collect_dekp_structures.py`, `src/33_run_dekp_public_retrained.py`, `src/34_evaluate_dekp_predictions.py` | `sbatch scripts/runners/run_dekp_public_retrained.sbatch` |
| SELFprot | `src/35_prepare_selfprot_eval.py`, `src/36_run_selfprot_predictions.py`, `src/37_evaluate_selfprot_predictions.py` | `sbatch scripts/runners/run_selfprot_predictions.sbatch` |
| DLKcat | `src/38_run_dlkcat_official.py` followed by `src/41_evaluate_method_predictions.py` | Run the two Python scripts in the DLKcat environment |
| UniKP | `src/39_run_unikp_official_features.py`, `src/40_predict_unikp_official_py36.py`, then `src/41_evaluate_method_predictions.py` | Run feature extraction and prediction in their required environments |
| TurNuP | `src/43_prepare_turnup_eval.py`, `src/45_run_turnup_predictions.py`, then `src/41_evaluate_method_predictions.py` | `bash scripts/runners/run_prepare_turnup_eval.sh`, then `sbatch scripts/runners/run_turnup_full.sbatch` |
| GO-HKP | `src/49_prepare_go_hkp_with_yeast_uniprot_go.py`, then `src/52_evaluate_go_hkp_predictions.py` | Run both Python scripts for combined *E. coli* and yeast coverage |

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
| `src/51_rebuild_catpred_overlap_audit.py` | Rebuild CatPred reference-corpus overlap levels | `python src/51_rebuild_catpred_overlap_audit.py` |
| `src/52_evaluate_go_hkp_predictions.py` | Rebuild GO-HKP row-level output, metrics, and missing summary | `python src/52_evaluate_go_hkp_predictions.py` |

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

### Publication tables, figures, and audit

The public v1.2.0 tables and figures are rebuilt from current row-level outputs:

```bash
# build_submission_audits.py additionally requires RDKit, SciPy, DIAMOND,
# and locally available method training corpora.
python paper/build_submission_audits.py
python paper/build_table0.py
python paper/rebuild_paper_tables.py
python paper/generate_manuscript_figures.py
python paper/recalculate_cluster_inference_v1_2.py
```

`build_submission_audits.py` creates S16-S24 and the 132-field Record_audit:
matching and shared-label dependence, substrate-support and role sensitivities,
five cluster-bootstrap definitions, cluster-level paired tests, standardized
training proximity, experimental dispersion, mutation-screen stages, model
direction/reversibility, and the three PreTKcat reconstruction policies.

The v1.2.0 artifact audit validates canonical data, all method predictions,
Table0, Table0W, Record_audit, S1-S24, the consolidated XLSX, and all four
figures. The manuscript-matched build passed 299/299 encoded checks. The
submission-specific DOCX source and target-Word pagination remain part of the
journal workflow rather than the public code release.

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
- Coverage: number of scored rows divided by 1,246.

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
paper/               Versioned manuscript tables, figures, builders, and audit
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
  should report sequence, SMILES, and pair overlap. The v1.2.0 primary PreTKcat and DEKP public-data reconstructions exclude
  standardized benchmark pairs before fitting; PreTKcat raw and near-excluded
  variants are reported separately in S24. S20 reports exact, joint-threshold-neighbor, and
  `no_joint_neighbor_under_thresholds` counts and performance where record-level
  corpora are available. The last class only means that no single training
  record meets both sequence identity >=80% and chemical similarity >=0.80;
  separate sequence-only or chemical-only overlap may remain. Unavailable
  corpora are marked `unknown`, not zero.
- The original project code is released under the MIT License in `LICENSE`.
- Third-party source code, checkpoints, foundation models, and database-derived
  content retain their upstream licenses and terms. Inclusion in a Zenodo
  manifest or archive does not relicense those materials. See
  `THIRD_PARTY_NOTICES.md` and `external_methods/METHOD_SOURCES.md`.

## Citation

The all-versions repository DOI is
[10.5281/zenodo.21024683](https://doi.org/10.5281/zenodo.21024683). Benchmark
v1.2.0 is archived under the version DOI
[10.5281/zenodo.21859994](https://doi.org/10.5281/zenodo.21859994).
When using the benchmark, also cite the original prediction methods and the
BRENDA, SABIO-RK, UniProt, PubChem, and metabolic-model resources used for the
relevant analysis.
