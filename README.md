# kcat Benchmark Analysis

一个面向已发表 kcat 预测方法的统一、可追溯评测项目。当前标准集包含 978 条实验 kcat：*Escherichia coli* 513 条、*Saccharomyces cerevisiae* 465 条。项目统一了实验真值、蛋白序列、底物 SMILES、方法输入/输出和 log10(kcat) 评估指标，并保留每条记录的数据库来源与匹配层级。

An executable benchmark of published kcat prediction methods on 978 experimentally supported enzyme-substrate records from *E. coli* and *S. cerevisiae*.

<!-- ZENODO_START -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21024684.svg)](https://doi.org/10.5281/zenodo.21024684)

**Large assets:** [Zenodo record 21024684](https://zenodo.org/records/21024684) (DOI: `10.5281/zenodo.21024684`). File checksums and restore paths are in `zenodo_assets_manifest.csv`.
<!-- ZENODO_END -->

## Benchmark at a glance

| item | value |
| --- | ---: |
| experimental truth before model-input filtering | 1,072 rows |
| unified sequence+SMILES benchmark | 978 rows |
| *E. coli* | 513 rows / 451 reactions / 327 genes |
| *S. cerevisiae* | 465 rows / 224 reactions / 168 genes |
| experimental sources | BRENDA and SABIO-RK |
| evaluated methods | 13 current methods |

标准集的详细构建方法、反应分布、GO/KEGG 注释和目录分类见 [`reports/kcat_benchmark_dataset_and_method_context.md`](reports/kcat_benchmark_dataset_and_method_context.md)。完整结果解释见 [`reports/kcat_benchmark_analysis_report.md`](reports/kcat_benchmark_analysis_report.md)。

## What is compared

方法按实际输入和可评测范围分组，避免把不同口径的数字硬排在一起：

| group | methods | comparison scope |
| --- | --- | --- |
| sequence + substrate SMILES | DLKcat, UniKP, MTLKP, CataPro, KcatNet, PreTKcat, SELFprot | 977/978 near-full benchmark |
| reaction-aware | TurNuP, PMAK | 780 rows with complete reaction SMILES |
| method-specific subset | CatPred, KinForm | rows accepted by the official pipeline |
| public-data retraining | DEKP-public-retrained | reproducible retraining, not an official paper checkpoint |
| GO functional assignment | GO-HKP | non-AI baseline, 978/978 rows |

主指标在 log10(kcat) 空间计算。这里 MAE=1 表示平均约相差一个 10 倍数量级；同时报告 RMSE、Pearson、Spearman、bias、10 倍误差内比例和覆盖率。

## Canonical data files

| file | role |
| --- | --- |
| `data/final/experimental_kcat_truth.csv` | 匹配到代谢模型条目的实验真值全集，尚未要求模型输入完整 |
| `data/final/benchmark_ready_truth.csv` | 有单蛋白序列和底物 SMILES 的统一真值表 |
| `data/final/benchmark_ready_catpred.csv` | 978 条统一母表；`catpred` 只是历史文件名，并非 CatPred 专用输入 |

各方法的 `*_input.csv` 只放推理所需字段，`*_metadata.csv` 保存物种、反应、来源和处理状态，`*_truth.csv` 仅在预测完成后评分。这样做是为了避免把实验答案混入模型输入。

## Repository and Zenodo split

GitHub 保存可审阅的代码、配置、核心标准集、报告、图和小型统计表。以下内容不进入 Git 历史：第三方源码副本、模型权重、特征缓存、原始 BRENDA/compound 数据库、大体积方法输出、早期原型结果，以及 PPTX/HTML 展示文件。

Zenodo 资产清单记录在 `zenodo_assets_manifest.csv`，每个文件都有大小、SHA256、恢复位置和适用方法。为避免慢网络上传超时，大包按 128 MiB 分片；下载脚本会逐片校验、自动拼接并再次核对整包 SHA256。第三方方法源码从其官方仓库获取，固定版本和本地目录见 [`external_methods/METHOD_SOURCES.md`](external_methods/METHOD_SOURCES.md)。

## Quick start

```bash
git clone https://github.com/dengxiao01/kcat_benchmark_analysis.git
cd kcat_benchmark_analysis
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

查看、下载和恢复 Zenodo 大文件：

```bash
python scripts/download_zenodo_assets.py --list
python scripts/download_zenodo_assets.py --asset kcat_benchmark_core_data_and_results.tar.gz --restore
```

重新生成统一结果汇总、主报告、数据集/方法报告以及报告表格导出：

```bash
python src/41_evaluate_method_predictions.py
python src/17_build_method_eval_summary.py
python src/46_generate_benchmark_report.py
python src/47_generate_dataset_method_context_report.py
python src/50_export_report_tables.py
```

只阅读现有结果时无需安装各方法的深度学习环境。真正重跑某个预测器时，需要先按 `external_methods/METHOD_SOURCES.md` 放置上游源码，再恢复相应 Zenodo 任务权重，并从官方来源下载 ProtT5/ESM1b 等通用基础模型；不同论文代码依赖的 Python/PyTorch 版本不同，应使用各自上游环境，不建议强行合成一个 Conda 环境。

可执行的 shell 和 Slurm 入口统一放在 `scripts/runners/`。例如：

```bash
bash scripts/runners/run_prepare_catpred_eval.sh
sbatch scripts/runners/run_catpred_full.sbatch
```

## End-to-end workflow

1. `src/01_parse_models.py` 从 `eciML1515.json` 和 `yeast-GEM.xml` 生成酶-反应-候选底物条目。
2. `src/05_fetch_uniprot_sequences.py` 补蛋白序列；`src/06_fetch_metanetx_smiles.py`、`src/07_fetch_pubchem_smiles.py` 补底物结构。
3. `src/08_fetch_sabiork_kcat.py` 和 `src/10_parse_brenda_kcat.py` 整理实验 kcat；`src/03_match_experimental_kcat.py` 按物种、EC、UniProt、底物 ID/名称分层匹配。
4. `src/11_finalize_benchmark_data.py` 生成 1,072 条实验真值和 978 条 sequence+SMILES benchmark。
5. `src/12_*` 至 `src/49_*` 是当前方法的输入适配、推理、评分和 GO/结构补齐脚本。
6. `src/41_evaluate_method_predictions.py` 统一计算指标，`src/46_*`、`src/47_*` 生成报告和图。

联网脚本会把查询结果写入 `data/raw/` 或 `data/interim/` 作为本地缓存。原始数据库文件没有打包进 Git；如需从头构建，请自行下载相应版本并遵守数据库许可。

## Project layout

```text
configs/             species, method, and matching rules
src/                 benchmark construction and method adapters
scripts/             runners and release/download utilities
scripts/runners/     benchmark, method prediction, and Slurm entry points
data/final/          canonical benchmark tables; method outputs live here locally
reports/             manuscript-facing reports, tables, and figures
external_methods/    local third-party checkouts; only METHOD_SOURCES.md is tracked
docs/                chronological work log and reproducibility notes
release/             local Zenodo staging; never tracked by Git
```

更细的“benchmark 构建、GO、KEGG、MAE、species-level、method-level”等目录归属已写入方法学报告，并导出为 `reports/tables/project_directory_analysis_map.csv`。

## Data provenance and caveats

- BRENDA turnover number 和 SABIO-RK kcat 是主真值。只保留正值并统一为 `s^-1`；BRENDA 默认排除 mutant/variant 条目。
- 匹配优先级为 `species+EC+UniProt+substrate ID`、`species+EC+substrate ID`、`species+EC+UniProt+substrate name`、`species+EC+substrate name`。
- 同一 benchmark entry 先保留最高匹配层级，再在 kcat 原始尺度取中位数，最后计算 log10。
- GO-HKP 的 *E. coli* 部分使用 DeepGO-SE 反应级赋值，yeast 部分使用 UniProt GO 注释；两种来源在 metadata 中分别标注。
- 早期通过推断或数据库填充得到的 kcat 不进入统一 benchmark，也不在当前公开仓库中分发。
- 模型训练集重叠仍可能使结果偏乐观；发表时应继续报告 sequence、SMILES 和 sequence-SMILES pair 的去重/重叠分析。

## Licensing and citation

本项目原创代码按 [`LICENSE`](LICENSE) 中的 MIT License 发布。核心派生数据和报告应同时注明其来源数据库；BRENDA 数据当前为 CC BY 4.0。第三方源码、权重和预训练模型不由本项目重新许可，仍服从各自上游条款，详见 `external_methods/METHOD_SOURCES.md` 和 Zenodo 资产说明。

论文发表前请补充正式作者、文章题目和方法论文引用。Zenodo 发布后，可使用仓库顶部 DOI 引用本 benchmark 数据与结果。
