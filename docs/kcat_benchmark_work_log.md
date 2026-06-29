# kcat Benchmark 工作记录

> 本文件作为项目台账持续维护。以后新增数据、脚本、评测方法、失败原因和人工修订，都追加到这里，方便写文章时追溯。

## 0. 术语说明

- **kcat**：酶的转换数，可以通俗理解为“一个酶分子每秒能处理多少个底物分子”，常用单位是 `s^-1`。
- **SMILES**：一种把小分子结构写成字符串的格式，类似“化学结构的文本身份证”。CatPred、DLKcat 这类方法通常需要它来表示底物。
- **实验真值**：来自 BRENDA、SABIO-RK 等数据库的实验测定 kcat，用作评测预测方法好坏的参照。
- **log10(kcat)**：把 kcat 取 10 的对数。这样比较误差时，不会让特别大的 kcat 数值主导全部结果。
- **benchmark-ready row**：同时具备蛋白序列、底物 SMILES 和实验 kcat 的一行记录，可直接进入预测方法评测。

## 1. 当前总体目标

本项目要把已发表的 kcat 预测方法放到同一套数据和指标下比较。当前优先级是先把数据准备、方法输入、预测输出合并和指标评估走通，然后逐个接入方法。

第一条接入的方法定为 **CatPred**。

## 2. 已完成的数据准备

### 2.1 输入数据来源

| 数据 | 当前位置 | 用途 |
|---|---|---|
| 模型反应和基因条目 | `data/interim/enzyme_reaction_entries_with_sequence_smiles.csv` | 形成“物种-反应-基因-底物”的候选评测条目 |
| UniProt 序列 | 由 `src/05_fetch_uniprot_sequences.py` 获取 | 给预测方法提供蛋白序列 |
| MetaNetX SMILES | 由 `src/06_fetch_metanetx_smiles.py` 获取 | 补充底物结构 |
| PubChem SMILES | 由 `src/07_fetch_pubchem_smiles.py` 获取 | 通过 PubChem API 继续补结构 |
| CKB 本地化合物库 | `data/raw/compounds/ckb/compounds.sqlite` | 用 BiGG、MetaNetX、KEGG、ChEBI 和同义名映射底物 |
| SABIO-RK kcat | `data/raw/sabiork/sabiork_kcat_raw.csv` | 实验 kcat 真值来源之一 |
| BRENDA kcat | `data/raw/brenda/brenda_kcat_raw.csv` | 实验 kcat 真值来源之一 |
| BRENDA 原始包 | `data/raw/brenda/brenda_2026_1.json.tar.gz`, `data/raw/brenda/brenda_2026_1.txt.tar.gz` | 可重复解析和核查 |

### 2.2 关键脚本

| 脚本 | 作用 |
|---|---|
| `src/01_parse_models.py` | 解析模型，整理反应、基因、底物候选 |
| `src/02_prepare_method_inputs.py` | 生成预测方法通用输入 |
| `src/03_match_experimental_kcat.py` | 将实验 kcat 匹配到模型条目 |
| `src/04_build_curation_queues.py` | 生成需要人工或外部库补齐的队列 |
| `src/05_fetch_uniprot_sequences.py` | 获取 UniProt 蛋白序列 |
| `src/06_fetch_metanetx_smiles.py` | 获取 MetaNetX SMILES |
| `src/07_fetch_pubchem_smiles.py` | 获取 PubChem 和 CKB SMILES |
| `src/08_fetch_sabiork_kcat.py` | 获取 SABIO-RK kcat |
| `src/10_parse_brenda_kcat.py` | 解析 BRENDA kcat |
| `src/11_finalize_benchmark_data.py` | 生成最终 benchmark-ready 表和难映射清单 |

对应运行脚本包括：

```bash
bash run_phase1.sh
bash run_sequence_fetch.sh
bash run_smiles_fetch.sh
bash run_sabiork_fetch.sh
bash run_brenda_parse.sh
bash run_finalize_benchmark_data.sh
```

### 2.3 当前最终产物

| 文件 | 行数 | 说明 |
|---|---:|---|
| `data/final/experimental_kcat_truth.csv` | 1072 | 已匹配到模型条目的实验 kcat 真值 |
| `data/final/benchmark_ready_truth.csv` | 978 | 真值中同时有序列和 SMILES、可进入 SMILES 类预测器评测的条目 |
| `data/final/benchmark_ready_catpred.csv` | 978 | CatPred/CataPro 类方法可用的评测主表 |
| `data/interim/prediction_inputs/catpred_kcat_input.csv` | 11692 | 全量 CatPred 候选预测输入，不限是否已有真值 |
| `data/interim/unresolved_smiles_review.csv` | 471 | 暂时不适合自动补 SMILES 的底物 |
| `reports/tables/final_benchmark_readiness.csv` | 2 | 按物种汇总的最终可评测情况 |
| `reports/tables/unresolved_smiles_summary.csv` | 13 | 难映射底物分类汇总 |

`benchmark_ready_catpred.csv` 当前组成：

| 物种 | 可评测行数 |
|---|---:|
| ecoli | 513 |
| yeast | 465 |
| 合计 | 978 |

实验 kcat 来源：

| 来源 | 行数 |
|---|---:|
| BRENDA | 587 |
| SABIO-RK | 310 |
| BRENDA;SABIO-RK | 81 |

## 3. 难映射 SMILES 的处理原则

当前仍有 471 个底物没有自动补出可靠 SMILES。它们不是简单“查不到”，而是很多本身就不适合强行当作普通小分子处理。

主要类别：

| 类型 | 简单解释 | 当前处理 |
|---|---|---|
| protein_or_redox_carrier | 蛋白、电子载体或辅因子蛋白，比如 cytochrome、ferredoxin | 不强行给小分子 SMILES |
| nucleic_acid_polymer | DNA、RNA、tRNA 等核酸聚合物 | 需要专门结构表示或排除 |
| polymer_or_glycan | 多糖、糖链、肽聚糖等聚合物 | 需要确定聚合度和结构 |
| metal_cluster_or_cofactor_complex | 铁硫簇等金属簇 | 需要人工确定精确结构 |
| curated_lipid_needed | 模型里的脂质统称或可变脂质 | 需要专门脂质库映射 |
| ambiguous_or_variable_structure | 名称本身含糊或结构可变 | 需要人工指定 |
| no_database_mapping | 常规库未找到可靠映射 | 需要提供可信交叉引用或 SMILES |

如果后续希望继续扩大 benchmark，需要提供这些项之一：

- 精确 `SMILES` 或 `InChI`
- 可信数据库 ID，例如 KEGG、ChEBI、MetaNetX、BiGG、PubChem CID
- 对脂质、聚合物、金属簇等特殊底物的人工结构表

## 4. CatPred 作为第一个评测方法

### 4.1 源码下载状态

CatPred 已下载到：

```text
external_methods/CatPred
```

来源：

```text
https://github.com/maranasgroup/CatPred.git
```

当前 commit：

```text
8e72d324e9e6f7a9a24c3f8a720884c7c1740a9b
```

### 4.2 CatPred 输入格式

CatPred 的 kcat 批量预测输入至少需要：

| 列名 | 含义 |
|---|---|
| `SMILES` | 底物结构字符串 |
| `sequence` | 蛋白序列 |
| `pdbpath` | CatPred 用来索引蛋白记录的键 |

CatPred demo 里还包含 `Substrate`，用于记录底物名称。我们保留该列，便于检查。

CatPred 的 `scripts/create_pdbrecords.py` 会把 `pdbpath` 的 basename 当作蛋白键。如果同一个 `pdbpath` 对应不同序列，脚本会报错。因此我们现在把 `pdbpath` 改成由蛋白序列生成的稳定哈希，例如：

```text
seq_28f4584eabf40332
```

这样做的意思是：同一条序列永远得到同一个键，不同序列基本不会撞到同一个键。

### 4.3 我们和 CatPred 数据格式的对比

| 项目 | CatPred demo | 我们当前 benchmark |
|---|---|---|
| 主要用途 | 展示如何预测 | 系统评测，需要真值和来源追踪 |
| 输入列 | `Substrate`, `SMILES`, `sequence`, `pdbpath` | 另含 `entry_id`, `species`, `reaction_id`, `gene_id`, `uniprot_id`, `ec_number`, `true_kcat`, `source_database` 等 |
| 真值 | demo 输入不带实验真值 | 每行带实验 kcat 和 log10(kcat) |
| 物种 | demo 小样本 | 当前包括 ecoli 和 yeast |
| 可追溯性 | 主要面向预测演示 | 保留数据库来源、匹配层级、参考文献和测量次数 |
| 输出合并 | demo 直接生成预测表 | 我们用 metadata 按行顺序或 `entry_id` 合并预测和真值 |

可借鉴的地方：

- 保留 CatPred 原生输入格式，减少运行时兼容问题。
- 把 `pdbpath` 设计成每条唯一序列对应一个唯一键，避免 protein record 冲突。
- 将方法输入、评测元数据、实验真值分开保存。这样预测方法只看到它需要的列，评估脚本仍能追踪每个预测来自哪条反应和哪条实验证据。
- 评估输出要保留 CatPred 的不确定性列，例如 `SD_total`, `SD_aleatoric`, `SD_epistemic`。这些列后续可用于分析“模型对哪些样本不自信”。

### 4.4 已新增的 CatPred 评测脚本

| 文件 | 作用 |
|---|---|
| `src/12_prepare_catpred_eval.py` | 从 `benchmark_ready_catpred.csv` 生成 CatPred 输入、metadata 和 truth |
| `src/13_evaluate_catpred_predictions.py` | 读取 CatPred 输出，计算 MAE、RMSE、R2、Pearson、Spearman 等指标 |
| `run_prepare_catpred_eval.sh` | 一键生成 CatPred 评测输入 |
| `run_catpred_predict.sh` | 检查 CatPred checkpoint，运行预测，并调用评估脚本 |

已经生成的 CatPred 输入产物：

| 文件 | 说明 |
|---|---|
| `data/final/catpred/catpred_kcat_input.csv` | 正式 CatPred 输入，978 行 |
| `data/final/catpred/catpred_kcat_input_metadata.csv` | 评估用元数据，978 行 |
| `data/final/catpred/catpred_kcat_input_truth.csv` | 评估真值，978 行 |
| `data/final/catpred/catpred_kcat_input.json.gz` | CatPred protein records，491 条唯一序列 |
| `data/final/catpred/catpred_kcat_input_sample100.csv` | 小样本烟测输入，100 行 |
| `data/final/catpred/catpred_kcat_input_sample100.json.gz` | 小样本 protein records |
| `reports/tables/catpred_eval_readiness.csv` | CatPred 评测输入统计 |

当前 CatPred 输入统计：

| 分组 | 行数 | 唯一序列 | 唯一 SMILES | 唯一底物 |
|---|---:|---:|---:|---:|
| all | 978 | 491 | 284 | 347 |
| ecoli | 513 | 327 | 232 | 232 |
| yeast | 465 | 164 | 116 | 116 |

## 5. CatPred 评测流程

### 5.1 已跑通的步骤

生成输入：

```bash
bash run_prepare_catpred_eval.sh
```

验证 CatPred protein records 可生成：

```bash
python3 external_methods/CatPred/scripts/create_pdbrecords.py \
  --data_file data/final/catpred/catpred_kcat_input.csv \
  --out_file data/final/catpred/catpred_kcat_input.json.gz
```

评估脚本烟测：

```bash
python3 src/13_evaluate_catpred_predictions.py \
  --predictions /tmp/catpred_eval_smoke_predictions.csv \
  --metadata data/final/catpred/catpred_kcat_input_sample100_metadata.csv \
  --out-rows /tmp/catpred_eval_smoke_rows.csv \
  --out-metrics /tmp/catpred_eval_smoke_metrics.csv
```

烟测结果：

```text
Evaluated rows: 100
MAE log10: 0
RMSE log10: 0
```

这里使用的是临时构造的“预测值=真值”，目的只是确认合并和指标计算没有问题，不代表 CatPred 真实表现。

### 5.2 当前未完成的步骤

正式运行 CatPred 预测时，当前阻塞在 checkpoint 数据：

```text
CatPred kcat checkpoint directory is missing:
external_methods/CatPred_capsule/data/pretrained/production/kcat
```

CatPred README 要求下载：

```text
https://catpred.s3.us-east-1.amazonaws.com/capsule_data_update.tar.gz
```

已检查该文件大小约为：

```text
10,207,582,747 bytes
```

也就是约 10.2 GB。文件系统空间足够，但下载和解压需要明确确认。下载后还需要检查解压目录结构，并安装 CatPred 所需环境。

正式预测命令预期为：

```bash
bash run_catpred_predict.sh
```

如果 CatPred checkpoint 不放在默认目录，可以指定：

```bash
CHECKPOINT_DIR=/path/to/data/pretrained/production/kcat bash run_catpred_predict.sh
```

或指定 capsule 根目录：

```bash
CATPRED_CAPSULE_DIR=/path/to/CatPred_capsule bash run_catpred_predict.sh
```

正式预测完成后，评估结果会写到：

```text
data/final/catpred/catpred_kcat_predictions_evaluated.csv
reports/tables/catpred_eval_metrics.csv
```

## 6. 评估指标约定

主指标在 `log10(kcat)` 空间计算。

| 指标 | 通俗解释 |
|---|---|
| `MAE log10` | 平均绝对误差，越小越好 |
| `RMSE log10` | 均方根误差，对大错更敏感，越小越好 |
| `bias log10` | 平均偏高或偏低，正值表示整体预测偏大 |
| `R2 log10` | 解释真值变化的程度，越接近 1 越好 |
| `Pearson log10` | 线性相关性，越高说明整体趋势越一致 |
| `Spearman log10` | 排名相关性，越高说明大小排序越一致 |
| `within_0.3_log10_fraction` | 误差在约 2 倍以内的比例 |
| `within_1.0_log10_fraction` | 误差在 10 倍以内的比例 |

后续论文中建议同时报告整体指标、按物种指标、按数据来源指标。

## 7. 下一步

1. 确认是否下载 CatPred 的 10.2 GB capsule 数据包。
2. 解压并定位 `data/pretrained/production/kcat` checkpoint。
3. 安装或激活 CatPred 所需 conda 环境。
4. 运行 `bash run_catpred_predict.sh`。
5. 检查 `reports/tables/catpred_eval_metrics.csv`，形成 CatPred 第一版评测结果。
6. 若要扩大评测集，优先人工补充 `data/interim/unresolved_smiles_review.csv` 中可明确结构的小分子、脂质和特殊辅因子。

## 8. 工作日志

### 2026-06-11

- 整理了 BRENDA、SABIO-RK、UniProt、CKB、MetaNetX、PubChem 的当前产物状态。
- 确认最终可评测 CatPred/CataPro 类输入为 978 行，其中 ecoli 513 行，yeast 465 行。
- 确认仍有 471 个难映射底物，不宜强行映射为普通小分子 SMILES。
- 下载 CatPred 源码到 `external_methods/CatPred`。
- 对比 CatPred demo 输入与本项目 benchmark 表，决定采用“CatPred 原生输入 + 独立 metadata/truth”的组织方式。
- 新增 `src/12_prepare_catpred_eval.py`，生成 CatPred 输入、metadata、truth 和 100 行烟测样本。
- 新增 `src/13_evaluate_catpred_predictions.py`，用于读取 CatPred 输出并计算评估指标。
- 新增 `run_prepare_catpred_eval.sh` 和 `run_catpred_predict.sh`。
- 成功生成 CatPred protein records：正式集 491 条唯一序列，烟测集也通过。
- 评估脚本用临时真值预测文件烟测通过，MAE/RMSE 均为 0。
- 正式 CatPred 预测暂未运行，原因是缺少 10.2 GB CatPred capsule 数据包中的 kcat checkpoint。

### 2026-06-12

- 检查 `external_methods/CatPred_capsule`，发现目录中只有 `capsule_data_update.tar.gz.1` 等残片，没有解压出的 `data/pretrained/production/kcat`。
- 对残片运行 `tar -tzf` 时出现 `Unexpected EOF in archive`，说明压缩包没有下载完整。
- 确认官方 S3 文件支持 `Accept-Ranges: bytes`，完整大小应为 `10,207,582,747` bytes。
- 新增 `run_download_catpred_weights.sh`，使用 `curl --continue-at -` 安全断点续传，先写 `.part` 文件，完整后再改名和解压，避免把半成品当成完整模型包。
- 用户重新下载完整 `capsule_data_update.tar.gz` 后，解压 `data/pretrained/production/kcat`，确认有 10 个 `model.pt`，kcat checkpoint 大小约 102M。
- 原始 `mamba env create -f external_methods/CatPred/environment.yml` 失败，原因是公共 conda cache 无写权限；改用项目本地 conda package cache 后，完整环境求解仍然较慢。
- 采用轻量方案：复用已有 `LLM_4_enzymes_env`，并把缺失依赖安装到 `external_methods/catpred_pydeps`，不污染原环境。实际补充了 `typed-argument-parser`、`ipdb`、`rotary_embedding_torch==0.6.5` 等必要包。
- 固定 CatPred 运行缓存：
  - ESM2 模型缓存：`external_methods/torch_cache`
  - 蛋白 ESM 嵌入缓存：`external_methods/catpred_esm_cache`
- 对 CatPred 源码做了一个本地兼容补丁：读取 ESM 缓存时，如果当前机器无 GPU，就把 CUDA tensor 映射到 CPU，避免 CPU 节点读取 GPU 缓存时报错。
- 更新 `run_catpred_predict.sh`：默认使用 `LLM_4_enzymes_env` 和 `external_methods/catpred_pydeps`，自动使用 GPU；若需要强制 CPU，可设置 `CATPRED_NO_CUDA=1`。
- 新增 `run_catpred_full.sbatch`，在 `qgpu_4090` 分区申请 1 张 GPU、8 CPU、64G 内存运行全量 CatPred 评测。
- 删除失败的半成品 `catpred_eval` conda 克隆环境，并清理不再使用的 `external_methods/conda_pkgs` 临时缓存。
- 100 行样本真实预测跑通。CatPred 判定 12 行 SMILES 无效，实际评估 88 行；样本 MAE log10 为 0.8408，RMSE log10 为 1.172。
- 全量 Slurm 作业 `1635965` 在 `qgpu_4090` / `gnode18` 正常完成，状态 `COMPLETED 0:0`，GPU 推理阶段日志显示 CatPred elapsed time 为 35 秒。
- 全量输入 978 行，CatPred 输出 913 行有效预测；65 行未预测，已写入 `data/final/catpred/catpred_invalid_or_unpredicted_rows.csv`。
- 未预测项汇总写入 `reports/tables/catpred_invalid_or_unpredicted_summary.csv`：其中 yeast 的 `[H+]` 64 行，ecoli 的 Quinate 1 行，其 SMILES 当前为错误数值 `192.167`。
- 全量 CatPred 评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 913 | 0.8487 | 1.1854 | 0.4051 | 0.4044 | 0.6933 |
| ecoli | 512 | 0.9061 | 1.2604 | 0.3363 | 0.3124 | 0.6602 |
| yeast | 401 | 0.7755 | 1.0822 | 0.4865 | 0.5250 | 0.7357 |
| BRENDA | 586 | 0.8032 | 1.1653 | 0.4065 | 0.4146 | 0.7201 |
| BRENDA;SABIO-RK | 81 | 0.9286 | 1.1860 | 0.3235 | 0.2529 | 0.5802 |
| SABIO-RK | 246 | 0.9307 | 1.2320 | 0.4488 | 0.4703 | 0.6667 |

### 2026-06-12 CatPred-DB 数据组织对比

- 查阅 CatPred 论文和本地 `external_methods/CatPred_datas/catpred-db.tar.gz`。论文说明 CatPred-DB 来自 BRENDA 2022_2 和 SABIO-RK 2023-11；每条记录要求有 kinetic parameter、EC、物种、反应物/产物、UniProt 序列映射和 canonical SMILES。kcat 对相同“酶序列 + 底物/反应物 SMILES”保留最大值，Km/Ki 保留几何平均。CatPred 的 kcat 输入使用所有反应物拼接 SMILES。
- 本地 CatPred-DB 包结构：
  - `CatPred-DB/data/kcat/`：kcat 随机 train/val/test 和总表。
  - `CatPred-DB/data/km/`、`CatPred-DB/data/ki/`：Km/Ki 对应数据。
  - `CatPred-DB/data/splits_revision/iteration_*`：多轮 train/test 划分，以及按序列相似性 40/60/80/99 cluster 的外推测试集。
  - `CatPred-DB/data/processed_databases/brenda.csv`、`sabio.csv`：处理后的来源数据库表。
- CatPred-DB kcat 主表 `kcat-random_trainvaltest.csv`：23151 行，7177 条唯一序列，7226 个 UniProt，10853 个唯一 `reactant_smiles`，12266 个唯一 `reaction_smiles`，2653 个 EC，1684 个 taxonomy ID。
- CatPred-DB kcat 随机划分：
  - train：18751 行
  - val：2084 行
  - test：2316 行
- 已生成比较表：
  - `reports/tables/catpred_db_kcat_summary.csv`
  - `reports/tables/catpred_db_kcat_splits_summary.csv`
  - `reports/tables/catpred_db_vs_our_benchmark_overlap.csv`
  - `reports/tables/catpred_db_vs_our_benchmark_overlap_summary.csv`
- 与我们 978 行 benchmark 的粗略重叠：
  - 精确 `sequence + reactant_smiles` 重叠：3 行。
  - 精确 `uniprot + reactant_smiles` 重叠：3 行。
  - 同 UniProt 且我们的单底物 SMILES 出现在 CatPred 反应物组件中：24 行。
  - UniProt 重叠：438 行。
  - EC 重叠：907 行。
- 解释：精确重叠低，部分原因是 CatPred kcat 使用“所有反应物拼接 SMILES”，而我们当前 benchmark 多数是一条模型底物 SMILES；SMILES canonical 化规则也可能不同。但 UniProt/EC 重叠较高，说明 CatPred 已发表模型与我们的 BRENDA/SABIO 来源 benchmark 不是严格独立。后续论文需要把 CatPred 结果表述为“在我们模型映射 benchmark 上的表现”，并另外设计去 CatPred-DB 重叠的外部测试子集。
- 借鉴项：
  - 在我们最终表中增加 method-independent 的 `input_id`/`row_id`，所有方法输出都按它回连。
  - 保留 `sequence_id`、`substrate_id/SMILES_id` 和可选的 `sequence_cluster`、`substrate_cluster`，用于做 in-distribution / out-of-distribution 分层评估。
  - 对同一输入键的多条 kcat 真值明确聚合规则：kcat 可保留最大值，同时保留原始测量数和来源；也可以另存 median/geomean 做敏感性分析。
  - 对 kcat 需要决定“单底物 SMILES”还是“反应物拼接 SMILES”。为了兼容 CatPred/DLKcat 这类底物模型，保留单底物输入；若评测反应级方法，应新增 `reactant_smiles` 和 `reaction_smiles`。

### 2026-06-16 CataPro 评测接入

- 按 `kcat_benchmark_executable_plan(1).md` 的第二优先级接入 **CataPro**。
- 官方代码下载到 `external_methods/CataPro`，来源为 `https://github.com/zchwang/CataPro`，当前 commit 为 `cc89b2c81768665cf6fd76dfda607ce88691f601`。
- CataPro 官方输入需要 4 个核心字段：
  - `Enzyme_id`：酶-底物条目的 ID。
  - `type`：野生型或突变体；本 benchmark 全部按 `wild` 处理。
  - `sequence`：蛋白序列。
  - `smiles`：底物 SMILES。
- CataPro 输出的 `pred_log10[kcat(s^-1)]` 已经是 `log10(kcat)`，评估时直接与 `true_kcat_log10` 比较。
- CataPro 仓库自带 10 折任务权重：
  - `external_methods/CataPro/models/kcat_models`
  - `external_methods/CataPro/models/Km_models`
  - `external_methods/CataPro/models/act_models`
- 另外按官方 README 下载了两个特征提取模型：
  - `external_methods/CataPro/models/prot_t5_xl_uniref50`
  - `external_methods/CataPro/models/molt5-base-smiles2caption`
- 环境处理：
  - 默认 Python 中缺 RDKit，`LLM_4_enzymes_env` 虽然依赖齐全，但导入 transformers 时会触发 deepspeed 的 CUDA 编译检查，报 `CUDA_HOME does not exist`。
  - 最终采用 `enyrnx` conda 环境，并将缺失的 `sentencepiece` 安装到项目本地 `external_methods/catapro_pydeps`。
  - 运行脚本设置 `PYTHONNOUSERSITE=1`，避免用户级 `.local` 包干扰 `huggingface_hub` 等依赖。
- 对 CataPro 源码做了两个本地兼容/效率补丁：
  - `external_methods/CataPro/inference/act_model.py`：让 `ActivityModel` 内部的 kcat/Km 子模型跟随传入的 `device`，避免 CPU 路线误建 CUDA 模型。
  - `external_methods/CataPro/inference/predict.py`：相同蛋白序列和相同 SMILES 只抽一次 ProtT5/MolT5 特征，并让 `batch_size` 真正传入 DataLoader。
- 新增脚本：
  - `src/14_prepare_catapro_eval.py`：从 `data/final/benchmark_ready_catpred.csv` 生成 CataPro 输入、metadata 和 truth。
  - `src/15_evaluate_catapro_predictions.py`：读取 CataPro 输出并计算统一 log10 指标。
  - `src/16_filter_catapro_valid_smiles.py`：用 RDKit 过滤 CataPro 无法处理的坏 SMILES。
  - `src/17_build_method_eval_summary.py`：汇总已评测方法的 overall 指标。
  - `run_prepare_catapro_eval.sh`、`run_catapro_predict.sh`、`run_catapro_sample.sbatch`、`run_catapro_full.sbatch`。
- 输入准备结果：
  - `data/final/catapro/catapro_kcat_input.csv`：978 行全量 CataPro 输入。
  - `data/final/catapro/catapro_kcat_input_valid_smiles.csv`：977 行 RDKit 可解析输入。
  - `data/final/catapro/catapro_kcat_input_metadata.csv`：978 行全量 metadata。
  - `data/final/catapro/catapro_kcat_input_truth.csv`：978 行真值表。
  - `reports/tables/catapro_eval_readiness.csv`：CataPro 输入覆盖汇总。
  - `reports/tables/catapro_valid_smiles_summary.csv`：SMILES 可解析性汇总。
- 唯一无效 SMILES 为：
  - `ecoli / Quinate / 192.167`，写入 `data/final/catapro/catapro_invalid_smiles_rows.csv`。
  - `[H+]` 对 RDKit 是合法 SMILES，只会产生 warning，不影响 CataPro 运行。
- 样本烟测：
  - Slurm 作业 `1643698` 在 `qgpu_4090` / `gnode19` 跑通 20 行样本。
  - 样本评估写入 `reports/tables/catapro_sample20_eval_metrics.csv`，MAE log10 为 0.9864，RMSE log10 为 1.199。
- 全量评测：
  - Slurm 作业 `1643699` 在 `qgpu_4090` / `gnode19` 完成。
  - CataPro 输出 `data/final/catapro/catapro_kcat_input_output.csv`：977 行预测。
  - 合并评估表 `data/final/catapro/catapro_kcat_predictions_evaluated.csv`：977 行。
  - 未预测行 `data/final/catapro/catapro_invalid_or_unpredicted_rows.csv`：1 行，即 Quinate 的坏 SMILES。
  - 未预测汇总 `reports/tables/catapro_invalid_or_unpredicted_summary.csv`。
- CataPro 全量评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 977 | 0.7763 | 0.9977 | 0.5177 | 0.5158 | 0.7247 |
| ecoli | 512 | 0.8412 | 1.0714 | 0.4473 | 0.4043 | 0.6602 |
| yeast | 465 | 0.7048 | 0.9095 | 0.5840 | 0.6045 | 0.7957 |
| BRENDA | 586 | 0.8251 | 1.0769 | 0.4631 | 0.4403 | 0.6877 |
| BRENDA;SABIO-RK | 81 | 0.7830 | 0.9446 | 0.4735 | 0.4000 | 0.6790 |
| SABIO-RK | 310 | 0.6823 | 0.8436 | 0.6372 | 0.5857 | 0.8065 |

- 与 CatPred 的 overall 汇总写入 `reports/tables/method_eval_summary.csv`：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。

### 2026-06-17 PMAK 评测接入

- 按 `kcat_benchmark_executable_plan(1).md` 的第三优先级接入 **PMAK**。
- 官方代码下载到 `external_methods/PMAK`，来源为 `https://github.com/MrVincentCai/PMAK`，当前 commit 为 `1b1bea4580ef7bb908f893d3b13213a1486bbb98`。
- PMAK 与 CatPred/CataPro 的重要区别：
  - CatPred/CataPro 可以使用“蛋白序列 + 单底物 SMILES”。
  - PMAK 需要“蛋白序列 + 完整 reaction SMILES”。reaction SMILES 可以简单理解为把反应写成 `反应物SMILES>>产物SMILES`。
  - 因此 PMAK 不能直接吃 `benchmark_ready_catpred.csv` 中的单底物 SMILES，必须先为每个模型反应补齐所有反应物和产物的 SMILES。
- 代码和权重核查：
  - PMAK README 主要提供复现实验流程：`split_data.py`、`get_prot5.py`、`get_rxnfp.py`、`Add_representation.py`、`train_kcat.py`。
  - 仓库中自带 PMAK 模型权重，位于 `external_methods/PMAK/code/save_model/`。
  - 本次使用 `external_methods/PMAK/code/save_model/CV/Fold_0_reaction_cold.pth` 到 `Fold_4_reaction_cold.pth` 这 5 个 reaction-cold checkpoint，并对 5 个模型预测值取平均。
  - PMAK 的预测输出按 `log10(kcat)` 理解，与 `true_kcat_log10` 直接比较。
- 环境处理：
  - 复用 `enyrnx` 环境，其中已有 `rxnfp`、`torch`、`transformers`、`rdkit` 等关键依赖。
  - 复用 CataPro 阶段已下载的 ProtT5 模型：`external_methods/CataPro/models/prot_t5_xl_uniref50`。
  - 继续使用项目本地 `external_methods/catapro_pydeps` 中的 `sentencepiece`。
  - 运行脚本设置 `PYTHONNOUSERSITE=1` 和 `LD_LIBRARY_PATH=/hpcfs/fhome/dengxg/.conda/envs/enyrnx/lib`，避免用户级 Python 包和系统 `libstdc++` 干扰。
- 新增脚本：
  - `src/18_prepare_pmak_eval.py`：从已有 benchmark 表、模型反应表和底物 SMILES 映射中构造 PMAK 输入。
  - `src/19_run_pmak_predictions.py`：生成 ProtT5 逐残基蛋白特征、RXNFP 反应特征，并加载 5 个 PMAK reaction-cold 权重预测。
  - `src/20_evaluate_pmak_predictions.py`：读取 PMAK 输出并计算统一 log10 指标。
  - `run_prepare_pmak_eval.sh`、`run_pmak_predict.sh`、`run_pmak_sample.sbatch`、`run_pmak_full.sbatch`。
  - `src/17_build_method_eval_summary.py` 已扩展纳入 PMAK。
- 输入准备结果：
  - 原始 benchmark-ready 表：978 行。
  - 可构造完整 reaction SMILES 且 RDKit 可解析的 PMAK 输入：780 行。
  - 因缺少至少一个反应物/产物 SMILES 或结构无效而不能进入 PMAK 的行：198 行。
  - 按物种：
    - ecoli：513 行中 468 行可评测，45 行缺完整 reaction SMILES。
    - yeast：465 行中 312 行可评测，153 行缺完整 reaction SMILES。
  - 完整反应子集中包含 543 个唯一模型反应、361 个唯一 reaction SMILES、426 条唯一蛋白序列。
- 关键输入/输出：
  - `data/final/pmak/pmak_kcat_input.csv`：780 行 PMAK 实际输入。
  - `data/final/pmak/pmak_kcat_input_metadata.csv`：978 行全量 metadata，包含 PMAK 是否有完整 reaction SMILES 的标记。
  - `data/final/pmak/pmak_kcat_input_truth.csv`：978 行真值表。
  - `data/final/pmak/pmak_missing_reaction_smiles_rows.csv`：198 行不能构造完整 reaction SMILES 的条目。
  - `reports/tables/pmak_eval_readiness.csv`：PMAK 输入覆盖率汇总。
- 样本烟测：
  - Slurm 作业 `1644909` 在 `qgpu_4090` / `gnode19` 跑通 20 行样本。
  - 样本评估写入 `reports/tables/pmak_sample20_eval_metrics.csv`，MAE log10 为 0.8794，RMSE log10 为 1.308。
- 全量评测：
  - Slurm 作业 `1644910` 在 `qgpu_4090` / `gnode19` 完成。
  - PMAK 输出 `data/final/pmak/pmak_kcat_input_output.csv`：780 行预测。
  - 合并评估表 `data/final/pmak/pmak_kcat_predictions_evaluated.csv`：780 行。
  - 未预测行 `data/final/pmak/pmak_invalid_or_unpredicted_rows.csv`：198 行。
  - 未预测汇总 `reports/tables/pmak_invalid_or_unpredicted_summary.csv`：yeast 153 行、ecoli 45 行，原因均为 `missing_metabolite_smiles`。
  - 特征缓存 `data/final/pmak/pmak_feature_cache.pkl`：缓存 426 条蛋白特征和 361 条 reaction 特征，后续重跑可复用。
- PMAK 全量评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 780 | 0.7227 | 1.0091 | 0.3828 | 0.4303 | 0.7679 |
| ecoli | 468 | 0.7622 | 1.0844 | 0.3506 | 0.4092 | 0.7350 |
| yeast | 312 | 0.6633 | 0.8843 | 0.4410 | 0.4425 | 0.8173 |
| BRENDA | 452 | 0.7304 | 1.0458 | 0.3405 | 0.3693 | 0.7544 |
| BRENDA;SABIO-RK | 70 | 0.7372 | 0.9528 | 0.3373 | 0.2630 | 0.7571 |
| SABIO-RK | 258 | 0.7052 | 0.9571 | 0.4749 | 0.5386 | 0.7946 |

- 与 CatPred/CataPro/PMAK 的 overall 汇总已更新到 `reports/tables/method_eval_summary.csv`：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。
  - PMAK：780 行，MAE log10 0.7227，RMSE log10 1.0091，Pearson 0.3828，Spearman 0.4303，10 倍内比例 0.7679。
- 注意事项：
  - PMAK 的 MAE 较低，但评测覆盖只有 780/978 行，覆盖率与 CatPred/CataPro 不同。
  - PMAK 使用完整 reaction SMILES，缺反应物/产物结构的条目会被排除；论文比较时需要同时报告覆盖率，并最好增加“共同可预测子集”的横向比较。

### 2026-06-17 KinForm 评测接入

- 按 `kcat_benchmark_executable_plan(1).md` 继续接入 **KinForm**。
- 官方代码下载到 `external_methods/KinForm`，来源为 `https://github.com/Digital-Metabolic-Twin-Centre/KinForm`，当前 commit 为 `f7a70eb1cd6723ba3a8d606432e522ea2b0fa9fd`。
- KinForm 与前面方法的关键区别：
  - KinForm 使用“蛋白序列 + 底物 SMILES”预测 `kcat`，形式上接近 CatPred/CataPro。
  - 但 KinForm-L/H 还依赖预计算的蛋白 embedding 和 binding-site 加权特征；embedding 可以理解为把蛋白序列提前翻译成模型可读的数字向量。
  - 官方把训练权重、序列 ID 映射、binding-site 分数和预计算蛋白特征打包在 Zenodo 的 `results.tar.gz` 中，记录为 `https://zenodo.org/records/17433514`。
- 新增脚本：
  - `src/21_prepare_kinform_eval.py`：从 `data/final/benchmark_ready_catpred.csv` 生成 KinForm JSON 输入、metadata 和 truth。
  - `src/22_check_kinform_coverage.py`：检查我们的 benchmark 序列是否存在于 KinForm 的 `sequence_id_to_sequence.pkl`，以及预计算向量文件是否齐全；只把“可直接预测”的行写入预测输入。
  - `src/23_evaluate_kinform_predictions.py`：读取 KinForm 输出并按统一 log10 指标评估。
  - `run_prepare_kinform_eval.sh`、`run_kinform_predict.sh`、`run_kinform_full.sbatch`。
  - `src/17_build_method_eval_summary.py` 已扩展纳入 KinForm。
- 输入准备结果：
  - 原始 benchmark-ready 表：978 行。
  - RDKit 可解析 SMILES：977 行。
  - 无效 SMILES：1 行，仍为 `ecoli / Quinate / 192.167`，写入 `data/final/kinform/kinform_invalid_smiles_rows.csv`。
  - KinForm 官方代码对超过 1499 aa 的 kcat 序列会保留前 749 aa 和后 749 aa；准备脚本已记录 `kinform_model_sequence` 和 `kinform_sequence_truncated`，用于后续和官方逻辑对齐。
- 环境处理：
  - 默认 Python 的 pandas/numpy/sklearn/torch 可用，且 sklearn 版本更接近 KinForm 权重环境，但缺 RDKit 和 scikit-mol。
  - 由于 KinForm-L 实际使用 SMILES Transformer，不使用 RDKit/scikit-mol 指纹；为避免无关顶层 import 卡住预测，在 `external_methods/kinform_pydeps` 增加了轻量 import shim。若后续切换到 Morgan/MACCS/MinHash 等指纹模式，则需要安装真实 RDKit/scikit-mol。
  - KinForm 官方 `smiles_features.py` 按当前目录查找 `vocab.pkl`；运行脚本会在 `external_methods/KinForm/code/vocab.pkl` 建软链接指向官方自带词表，不直接修改外部方法源码。
  - KinForm 权重由 scikit-learn 1.7.2/joblib 1.4.2 保存；默认 Python 是 scikit-learn 1.8.0/joblib 1.5.3，会有模型反序列化版本 warning。为减少版本差异，已将 `scikit-learn==1.7.2` 和 `joblib==1.4.2` 安装到项目本地 `external_methods/kinform_pydeps_exact`，并在预测脚本中置于 `PYTHONPATH` 最前。
- Zenodo 资源包：
  - 下载文件：`external_methods/KinForm/results.tar.gz`。
  - MD5：`da0df5e0e819a3665193cf4e61677f17`，与 Zenodo 记录一致。
  - 解压目录：`external_methods/KinForm/results/`，总大小约 9.3G。
  - 关键文件检查通过：
    - `results/sequence_id_to_sequence.pkl`
    - `results/binding_sites/binding_sites_all.tsv`
    - `results/trained_models/kcat_KinForm-L/model.joblib`
    - `results/trained_models/kcat_KinForm-L/transformers.joblib`
    - `results/protein_embeddings/{prot_t5_last,prot_t5_layer_19,esmc_layer_24,esmc_layer_32,esm2_layer_26,esm2_layer_29}/{mean_vecs,weighted_vecs}`
- 覆盖率检查：
  - KinForm 资源齐全后，978 行 benchmark 中有 563 行可直接预测，415 行不可直接预测。
  - 可预测子集包含 287 条唯一模型序列。
  - 按物种：
    - ecoli：513 行中 326 行可预测。
    - yeast：465 行中 237 行可预测。
  - 不可预测原因：
    - `missing_sequence_lookup;missing_precomputed_features`：215 行，表示序列不在官方映射表里，因此也没有预计算向量。
    - `missing_precomputed_features`：199 行，表示序列能在映射表中找到，但至少缺一部分 KinForm 需要的预计算向量。
    - `invalid_smiles`：1 行，即 Quinate 的坏 SMILES。
  - 覆盖率汇总写入 `reports/tables/kinform_feature_coverage_summary.csv`，逐行覆盖信息写入 `data/final/kinform/kinform_feature_coverage.csv`。
- 全量可预测子集评测：
  - 预测输入：`data/final/kinform/kinform_kcat_input_predictable.json`，563 行。
  - KinForm 输出：`data/final/kinform/kinform_kcat_input_output.csv`，563 行，`y_pred` 为原始 kcat，不是 log10。
  - 合并评估表：`data/final/kinform/kinform_kcat_predictions_evaluated.csv`，563 行。
  - 未预测行：`data/final/kinform/kinform_invalid_or_unpredicted_rows.csv`，415 行。
  - 未预测汇总：`reports/tables/kinform_invalid_or_unpredicted_summary.csv`。
  - 输出对齐检查通过：预测输出中的 sequence 和 SMILES 与 metadata 均 0 mismatch。
- KinForm 全量可预测子集评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 563 | 0.7827 | 0.9894 | 0.6028 | 0.6297 | 0.7158 |
| ecoli | 326 | 0.7888 | 1.0425 | 0.4697 | 0.5216 | 0.7209 |
| yeast | 237 | 0.7744 | 0.9113 | 0.7766 | 0.7344 | 0.7089 |
| BRENDA | 338 | 0.8040 | 1.0426 | 0.5644 | 0.6406 | 0.7337 |
| BRENDA;SABIO-RK | 56 | 0.6137 | 0.7606 | 0.6635 | 0.5978 | 0.7857 |
| SABIO-RK | 169 | 0.7962 | 0.9460 | 0.6915 | 0.5505 | 0.6568 |

- 与 CatPred/CataPro/PMAK/KinForm 的 overall 汇总已更新到 `reports/tables/method_eval_summary.csv`：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。
  - PMAK：780 行，MAE log10 0.7227，RMSE log10 1.0091，Pearson 0.3828，Spearman 0.4303，10 倍内比例 0.7679。
  - KinForm：563 行，MAE log10 0.7827，RMSE log10 0.9894，Pearson 0.6028，Spearman 0.6297，10 倍内比例 0.7158。
- 注意事项：
  - 当前 KinForm 指标只代表“官方 Zenodo 包中已经有完整预计算特征”的 563 条样本，不代表 978 条全量 benchmark。
  - 若希望让 KinForm 覆盖剩余 414 条因 embedding 缺失的样本，需要额外配置 KinForm 官方 Path B 中的 ESM-2、ESM-C、ProtT5 和 Pseq2Sites 环境并为新序列计算 embedding/binding-site 特征；这一步计算量和环境复杂度都明显高于当前直接评测流程。
  - PyTorch 在登录节点查询 CUDA 时会打印驱动版本 warning，但本次 KinForm-L 预测使用预计算蛋白特征和 CPU 路径完成，预测输出已正常生成。

### 2026-06-18 KcatNet 评测接入

- 按 `kcat_benchmark_executable_plan(1).md` 继续接入 **KcatNet**。
- 官方代码下载到 `external_methods/KcatNet`，来源为 `https://github.com/BioColLab/KcatNet`，当前 commit 为 `7d370f69f9d1bbed517655d23d4d80bd76594321`。
- KcatNet 的输入是“蛋白序列 + 底物 SMILES”，输出是预测的 kcat。通俗地说，KcatNet 会先把蛋白序列和底物结构各自转成模型能读的数字特征，再用图神经网络预测 `log10(kcat)`。
- 官方仓库已自带 kcat 预训练权重：
  - `external_methods/KcatNet/RESULT/model_KcatNet.pt`
  - `external_methods/KcatNet/Dataset/degree.pt`
  - `external_methods/KcatNet/utils/trfm_12_23000.pkl`
  - `external_methods/KcatNet/utils/vocab.pkl`
- 新增脚本：
  - `src/24_prepare_kcatnet_eval.py`：从 `data/final/benchmark_ready_catpred.csv` 生成 KcatNet 输入、metadata、truth 和 sample20 文件。
  - `src/25_run_kcatnet_predictions.py`：包装 KcatNet 官方模型，读取 CSV 输入，缓存蛋白/SMILES 特征，输出 `prediction_log10` 和 `prediction_kcat`。
  - `src/26_evaluate_kcatnet_predictions.py`：读取 KcatNet 输出并按统一 log10 指标评估。
  - `run_prepare_kcatnet_eval.sh`、`run_kcatnet_predict.sh`、`run_kcatnet_sample.sbatch`、`run_kcatnet_full.sbatch`。
  - `src/17_build_method_eval_summary.py` 已扩展纳入 KcatNet。
- 输入准备结果：
  - 原始 benchmark-ready 表：978 行。
  - RDKit 可解析 SMILES：977 行。
  - 无效 SMILES：1 行，仍为 `ecoli / Quinate / 192.167`，写入 `data/final/kcatnet/kcatnet_invalid_smiles_rows.csv`。
  - KcatNet 对超过 1000 aa 的序列按官方逻辑截取前 1000 aa；准备脚本已记录 `kcatnet_model_sequence` 和 `kcatnet_sequence_truncated`。
  - 可预测子集包含 490 条唯一蛋白序列、283 个唯一 canonical SMILES。
- 环境处理：
  - 继续使用 `enyrnx` 环境，因为它已有 RDKit、fair-esm、BioPython、openpyxl、torch-geometric 等多数依赖。
  - `enyrnx` 缺 `sentencepiece`，已复用项目本地 `external_methods/catapro_pydeps/sentencepiece`。
  - `enyrnx` 缺 `torch_scatter`，这是 KcatNet 图神经网络层需要的底层算子库。PyG 官方 pip/conda wheel 与当前 conda-forge PyTorch 的 ABI 不兼容，表现为 `undefined symbol: parseSchemaOrName`。
  - 为解决 ABI 问题，新增本地编译配置 `external_methods/kcatnet_build_sysconfig/_sysconfigdata_kcatnet_enyrnx.py`，把 `enyrnx` 中残留的旧用户路径替换为当前路径；随后用 CUDA 12.1 从源码编译 `torch-scatter==2.1.2` 到 `external_methods/kcatnet_scatter_src`。运行脚本通过 `PYTHONPATH` 优先加载这个源码编译版本。
  - ProtT5 模型复用 CataPro 已下载的本地目录：`external_methods/CataPro/models/prot_t5_xl_uniref50`。
  - ESM2 权重使用已有本地缓存：`external_methods/torch_cache/hub/checkpoints/esm2_t33_650M_UR50D.pt` 和 `esm2_t33_650M_UR50D-contact-regression.pt`。
  - 对 KcatNet 官方 `utils/protein_init.py` 做了最小运行环境修补：把硬编码 ProtT5 路径改成 `KCATNET_PROTT5_DIR`，并让 ESM/ProtT5 明确服从 `KCATNET_DEVICE`。这不是算法改动，只是路径和设备选择改动。
- 样本烟测：
  - Slurm 作业 `1645387` 在 `qgpu_4090` / `gnode19` 跑通 20 行样本。
  - 样本评估写入 `reports/tables/kcatnet_sample20_eval_metrics.csv`，MAE log10 为 0.891，RMSE log10 为 1.233。
- 全量评测：
  - Slurm 作业 `1645388` 在 `qgpu_4090` / `gnode19` 完成。
  - KcatNet 输出：`data/final/kcatnet/kcatnet_kcat_input_output.csv`，977 行预测。
  - 合并评估表：`data/final/kcatnet/kcatnet_kcat_predictions_evaluated.csv`，977 行。
  - 未预测行：`data/final/kcatnet/kcatnet_invalid_or_unpredicted_rows.csv`，1 行，原因为 `invalid_smiles`。
  - 未预测汇总：`reports/tables/kcatnet_invalid_or_unpredicted_summary.csv`。
  - 特征缓存：
    - `data/final/kcatnet/kcatnet_protein_cache.pkl`：缓存 490 条蛋白特征，约 1.34G。
    - `data/final/kcatnet/kcatnet_ligand_cache.pkl`：缓存 283 个 SMILES 特征，约 4.8M。
  - 后续重跑 KcatNet 时会直接复用缓存，不再重复计算 ProtT5/ESM/SMILES Transformer 特征。
- KcatNet 全量评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 977 | 0.6841 | 0.9962 | 0.4738 | 0.5149 | 0.7789 |
| ecoli | 512 | 0.7759 | 1.1227 | 0.3611 | 0.4054 | 0.7227 |
| yeast | 465 | 0.5830 | 0.8349 | 0.6241 | 0.6542 | 0.8409 |
| BRENDA | 586 | 0.7066 | 1.0490 | 0.3966 | 0.4643 | 0.7679 |
| BRENDA;SABIO-RK | 81 | 0.7667 | 0.9778 | 0.4587 | 0.3328 | 0.7037 |
| SABIO-RK | 310 | 0.6198 | 0.8932 | 0.5925 | 0.5670 | 0.8194 |

- 与已有方法的 overall 汇总已更新到 `reports/tables/method_eval_summary.csv`：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。
  - PMAK：780 行，MAE log10 0.7227，RMSE log10 1.0091，Pearson 0.3828，Spearman 0.4303，10 倍内比例 0.7679。
  - KinForm：563 行，MAE log10 0.7827，RMSE log10 0.9894，Pearson 0.6028，Spearman 0.6297，10 倍内比例 0.7158。
  - KcatNet：977 行，MAE log10 0.6841，RMSE log10 0.9962，Pearson 0.4738，Spearman 0.5149，10 倍内比例 0.7789。
- 注意事项：
  - KcatNet 当前覆盖 977/978 行，与 CataPro 覆盖一致；唯一缺失项仍是 Quinate 的错误 SMILES `192.167`。
  - KcatNet 在当前 benchmark 上 MAE log10 最低，但 Pearson/Spearman 低于 KinForm 和 CataPro；写文章时需要同时报告误差指标、相关性指标和覆盖率。
  - KcatNet 使用 ESM2 contact map 和 ProtT5 residue embedding，缓存文件较大；不要随意删除 `kcatnet_protein_cache.pkl`，否则重跑会重新计算蛋白特征。

### 2026-06-22 PreTKcat 评测接入

- 按 `kcat_benchmark_executable_plan(1).md` 继续接入 **PreTKcat**。
- 官方代码下载到 `external_methods/PreTKcat`，来源为 `https://github.com/MrVincentCai/PreTKcat`，当前 commit 为 `b7bc0562a9b8555a201c5f6c72fc2a660dcdb76d`。
- PreTKcat 的输入是“蛋白序列 + 底物 SMILES + 温度”。通俗地说，它先把蛋白序列和底物分子结构分别转成模型能理解的数字向量，再把温度一起交给 ExtraTrees 回归器预测 `log10(kcat)`。
- 官方 README 说明：
  - 蛋白表示使用 ProtT5。
  - 底物分子图表示使用 MolGNet。
  - kcat 预测使用 ExtraTrees。
  - 需要从 MPG 项目下载 MolGNet 预训练模型。
- 本地资源处理：
  - PreTKcat 仓库自带公开训练表 `external_methods/PreTKcat/datasets/DLTKcat_data/kcat_merge_DLTKcat.csv`，共 16,249 行。
  - 官方仓库没有发布可直接加载的 kcat ExtraTrees 成品权重，也没有 GitHub release；因此本次评测不是“加载作者回归器权重”，而是“用作者公开训练集和特征流程重新训练 ExtraTrees 后预测我们的 benchmark”。
  - MolGNet 权重已从 MPG 官方 Google Drive 下载到 `external_methods/PreTKcat/MolGNet.pt`，文件大小约 204 MB。
  - ProtT5 继续复用 CataPro 阶段的本地模型目录：`external_methods/CataPro/models/prot_t5_xl_uniref50`。
- 兼容性处理：
  - 复用 `enyrnx` 环境、本地 `external_methods/catapro_pydeps` 中的 `sentencepiece`，以及本地编译的 `external_methods/kcatnet_scatter_src`。
  - PreTKcat 的 `MPG_util/graph_bert.py` 面向旧版 PyG；当前 PyG 中 `softmax` 参数顺序和 `MessagePassing` 默认 `node_dim` 不同。本次做了两处最小兼容修补：
    - `softmax(alpha, edge_index_i, size_i)` 改为 `softmax(alpha, edge_index_i, num_nodes=size_i)`。
    - `GraphAttentionConv` 显式设置 `aggr='add', node_dim=0`。
  - 这两处只解决新旧 PyG API 差异，不改变 PreTKcat 的模型结构或权重。
- 新增脚本：
  - `src/27_prepare_pretkcat_eval.py`：从 `data/final/benchmark_ready_catpred.csv` 生成 PreTKcat 输入、metadata、truth 和 sample20 文件。
  - `src/28_run_pretkcat_predictions.py`：使用 PreTKcat 公开训练集训练 ExtraTrees，并对 benchmark 预测。
  - `src/29_evaluate_pretkcat_predictions.py`：读取 PreTKcat 输出并按统一 log10 指标评估。
  - `run_prepare_pretkcat_eval.sh`、`run_pretkcat_predict.sh`、`run_pretkcat_sample.sbatch`、`run_pretkcat_full.sbatch`。
  - `src/17_build_method_eval_summary.py` 已扩展纳入 PreTKcat。
- 输入准备结果：
  - 原始 benchmark-ready 表：978 行。
  - RDKit 可解析 SMILES：977 行。
  - 无效 SMILES：1 行，仍为 `ecoli / Quinate / 192.167`，写入 `data/final/pretkcat/pretkcat_invalid_smiles_rows.csv`。
  - benchmark 中缺温度 148 行；由于 PreTKcat 需要温度，本次按作者公开训练集的中位温度 30 摄氏度补齐，并在 `pretkcat_temperature_imputed` 中标记。
  - 温度全部在作者训练数据范围 0 到 100 摄氏度内。
  - 公开训练集中与 benchmark 发生“模型可见的序列 + canonical SMILES”精确重叠的行有 26 行；这些行在 metadata 中用 `pretkcat_train_exact_pair_overlap` 标记。写文章时需要说明这可能让 PreTKcat 指标偏乐观。
  - 可预测子集包含 491 条唯一蛋白序列、284 个唯一 SMILES。
- 关键输入/输出：
  - `data/final/pretkcat/pretkcat_kcat_input.csv`：978 行 PreTKcat 输入。
  - `data/final/pretkcat/pretkcat_kcat_input_valid_smiles.csv`：977 行实际可预测输入。
  - `data/final/pretkcat/pretkcat_kcat_input_metadata.csv`：978 行全量 metadata，包含温度补值、序列截断和训练集重叠标记。
  - `data/final/pretkcat/pretkcat_kcat_input_truth.csv`：978 行真值表。
  - `reports/tables/pretkcat_eval_readiness.csv`：PreTKcat 输入覆盖率汇总。
- 样本烟测：
  - MolGNet 小分子前向检查已通过，输出向量维度为 768。
  - 已完成极小 CPU 烟测：使用 2 行公开训练集、预测 2 行 benchmark，验证 ProtT5、MolGNet、ExtraTrees 和输出写表流程可以串联跑通。输出为 `data/final/pretkcat/pretkcat_kcat_input_tiny2_output.csv`。
  - 初始 sample/full 作业 `1648364`、`1648365` 因 GPU 资源等待较久，已取消并降低内存请求后重提。
  - PreTKcat sample20 烟测作业 `1648367` 在 `qgpu_4090` 完成。
  - sample20 使用 500 行公开训练集和 50 棵 ExtraTrees，只用于验证流程，不作为正式指标。
  - sample20 输出 `data/final/pretkcat/pretkcat_kcat_input_sample20_output.csv`，评估结果写入 `reports/tables/pretkcat_sample20_eval_metrics.csv`；MAE log10 为 1.085，RMSE log10 为 1.416。
- 全量评测：
  - PreTKcat full 作业 `1648368` 在 `qgpu_4090` 完成。
  - 全量特征提取规模：8,299 条唯一模型序列、2,719 个唯一 SMILES。
  - 特征缓存 `data/final/pretkcat/pretkcat_feature_cache.pkl`，约 45 MB。
  - 重新训练的 ExtraTrees 模型 `data/final/pretkcat/pretkcat_extratrees_model.pkl`，约 217 MB。
  - PreTKcat 输出 `data/final/pretkcat/pretkcat_kcat_input_output.csv`：977 行预测。
  - 合并评估表 `data/final/pretkcat/pretkcat_kcat_predictions_evaluated.csv`：977 行。
  - 未预测行 `data/final/pretkcat/pretkcat_invalid_or_unpredicted_rows.csv`：1 行，原因为 `invalid_smiles`，仍是 `ecoli / Quinate / 192.167`。
  - 未预测汇总 `reports/tables/pretkcat_invalid_or_unpredicted_summary.csv`。
- PreTKcat 全量评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 977 | 0.8623 | 1.1105 | 0.4087 | 0.4581 | 0.6489 |
| ecoli | 512 | 0.8609 | 1.1547 | 0.3393 | 0.3439 | 0.6816 |
| yeast | 465 | 0.8639 | 1.0596 | 0.5363 | 0.6042 | 0.6129 |
| BRENDA | 586 | 0.8404 | 1.1469 | 0.3808 | 0.4162 | 0.6911 |
| BRENDA;SABIO-RK | 81 | 0.8522 | 1.0505 | 0.3435 | 0.3379 | 0.6543 |
| SABIO-RK | 310 | 0.9064 | 1.0542 | 0.4671 | 0.4499 | 0.5677 |

- 训练集精确重叠分组：

| 是否与公开训练集 exact pair 重叠 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| False | 951 | 0.8691 | 1.1042 | 0.4234 | 0.4592 | 0.6456 |
| True | 26 | 0.6134 | 1.3184 | 0.4277 | 0.4817 | 0.7692 |

- 与已有方法的 overall 汇总已更新到 `reports/tables/method_eval_summary.csv`：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。
  - PMAK：780 行，MAE log10 0.7227，RMSE log10 1.0091，Pearson 0.3828，Spearman 0.4303，10 倍内比例 0.7679。
  - KinForm：563 行，MAE log10 0.7827，RMSE log10 0.9894，Pearson 0.6028，Spearman 0.6297，10 倍内比例 0.7158。
  - KcatNet：977 行，MAE log10 0.6841，RMSE log10 0.9962，Pearson 0.4738，Spearman 0.5149，10 倍内比例 0.7789。
  - PreTKcat：977 行，MAE log10 0.8623，RMSE log10 1.1105，Pearson 0.4087，Spearman 0.4581，10 倍内比例 0.6489。
- 注意事项：
  - PreTKcat 官方仓库没有发布可直接加载的 kcat ExtraTrees 回归器权重，因此当前结果是“基于公开训练集重训 ExtraTrees”的评测结果；与 KcatNet/CatPred 这类直接加载官方模型权重的评测性质不同。
  - PreTKcat 覆盖 977/978 行，覆盖率与 CataPro/KcatNet 一致；唯一缺失项仍是 Quinate 的坏 SMILES。
  - PreTKcat 在当前 benchmark 上整体 MAE 和 RMSE 均弱于 KcatNet/CataPro，也弱于 PMAK 的可覆盖子集；但它使用了温度特征，后续写文章时可以额外讨论“温度补值”和“训练集重叠”对结果的影响。

### 2026-06-22 DEKP 评测接入准备

- 按 `kcat_benchmark_executable_plan(1).md` 继续接入 **DEKP**。
- 官方代码下载到 `external_methods/DEKP`，来源为 `https://github.com/wang-yi-zhen/DEKP`，当前 commit 为 `d2b8c1372b5c1855fd2de9aaadde19cf8cc7fa8d`。
- DEKP 论文信息：
  - Briefings in Bioinformatics 2025，题目为 *DEKP: a deep learning model for enzyme kinetic parameter prediction based on pretrained models and graph neural networks*。
  - 论文 DOI：`10.1093/bib/bbaf187`。
  - DEKP 同时预测 `kcat` 和 `Km`，核心特点是把蛋白序列、底物 SMILES、蛋白三维结构图等信息一起输入模型。
- 通俗说明：
  - CatPred、CataPro、KcatNet 主要吃“蛋白序列 + 底物 SMILES”。
  - DEKP 还要吃“蛋白结构”。这里的蛋白结构不是一个编号，而是实际 PDB/AlphaFold 坐标文件；模型会把残基之间的空间距离关系转成图结构，再交给图神经网络。
  - 因此 DEKP 的评测前置资产明显更重，不能只用 `benchmark_ready_catpred.csv` 直接预测。
- 官方 README 核查：
  - 预训练依赖包括 ProtT5-XL-U50、SMILES Transformer、PST、MolFormer。
  - 额外依赖包括 PyTorch Geometric。
  - README 只说明 `Encode` 文件夹用于提取不同特征，并提供 kcat/Km 数据集和蛋白结构数据链接。
  - 官方仓库没有发布 GitHub release，也没有随仓库提供可直接推理的 kcat 模型权重。
- 已下载的官方数据：
  - kcat 数据集：`external_methods/DEKP/datasets/kcat_dataset.csv`，来源为 README 的 Google Drive 链接，13,401 行。
  - 该 kcat 表为 TSV，字段为 `ECNumber, Organism, Smiles, Substrate, Sequence, Type, Label, Unit, UniprotID, CID`。
  - `Label` 是 `log10(kcat)`，范围约 -3 到 3。
  - 官方蛋白结构数据包：`external_methods/DEKP/datasets/protein_structure_datasets.rar`，来源为 Zenodo `https://zenodo.org/records/15081759`，约 671 MB。
  - 结构包可读取，包含约 5,849 个 PDB 文件；其中 kcat AlphaFold 结构去重后约 4,047 个 UniProt ID。
- 代码核查结果：
  - `DEKP/fine_tune.py` 默认需要 `Model/Pretrained_model_199_trfm,t5.pkl` 作为预训练权重。
  - 仓库中没有这个文件，也没有其它 `.pt/.pth/.ckpt` 或 final predictor 权重。
  - `fine_tune.py` 会加载 `feature/{trfm,t5,pst,dssp,molformer}.pkl` 和 `feature/pyg_graph.pkl`；即使只用默认 `trfm,t5` 特征，模型仍然需要 `pyg_graph.pkl` 中的蛋白结构图。
  - `pretrain.py` 可以从公开数据训练 DEKP 模型，但这属于“从头复现/重训”，不是 off-the-shelf 推理；若采用该路线，文章中必须标注为 `DEKP-public-retrained`，不能和加载官方成品权重的方法混为一谈。
- 环境核查：
  - 继续使用 `enyrnx` 环境，其中 RDKit、BioPython、torch、torch_geometric 和本地编译的 `torch_scatter` 可用。
  - 当前环境缺 `torch_cluster`，而 DEKP 原始 `extract_pdb_feature.py` 的 `radius_graph` 依赖它。
  - 后续若继续重训路线，可以有两种选择：
    - 安装/编译 `torch_cluster`。
    - 或写一个兼容版结构图提取脚本，用 `scipy.spatial.cKDTree` 生成半径邻接图，绕开 `torch_cluster`。
- 新增脚本：
  - `src/30_prepare_dekp_eval.py`：从 `data/final/benchmark_ready_catpred.csv` 生成 DEKP 输入、metadata、truth、结构覆盖率和训练集重叠标记。
  - `src/31_download_dekp_missing_structures.py`：针对缺结构的 benchmark UniProt 尝试从 AlphaFold 下载 PDB 文件。
  - `run_prepare_dekp_eval.sh`。
- 输入准备结果：
  - 原始 benchmark-ready 表：978 行。
  - RDKit 可解析 SMILES：977 行。
  - 无效 SMILES：1 行，仍为 `ecoli / Quinate / 192.167`。
  - benchmark 中唯一 UniProt：495 个。
  - 作者结构包覆盖 benchmark 唯一 UniProt：194/495 个。
  - 作者结构包覆盖 benchmark 行：388/978 行。
  - 同时满足“SMILES 可解析 + 结构可用”的当前 DEKP 可进入结构分支行数：387/978 行。
  - 缺结构行：590 行，涉及 301 个 UniProt。
  - 按物种：
    - ecoli：513 行中 243 行当前可进入 DEKP，269 行缺结构，1 行坏 SMILES。
    - yeast：465 行中 144 行当前可进入 DEKP，321 行缺结构。
  - 与 DEKP 官方 kcat 数据集的重叠：
    - benchmark 中 277 行的 UniProt 出现在 DEKP kcat 数据集中。
    - 17 行发生 “Sequence + canonical SMILES” exact pair 重叠。
- 关键输入/输出：
  - `data/final/dekp/dekp_kcat_input.csv`：978 行 DEKP 格式输入。
  - `data/final/dekp/dekp_kcat_input_valid_smiles.csv`：977 行 SMILES 可用输入。
  - `data/final/dekp/dekp_kcat_input_structure_available.csv`：388 行结构可用输入，其中 387 行同时 SMILES 可用。
  - `data/final/dekp/dekp_kcat_input_metadata.csv`：978 行全量 metadata，包含 SMILES、结构覆盖和训练集重叠标记。
  - `data/final/dekp/dekp_missing_structure_rows.csv`：590 行缺结构条目。
  - `reports/tables/dekp_eval_readiness.csv`：按物种汇总的 DEKP 输入覆盖率。
  - `reports/tables/dekp_asset_coverage.csv`：DEKP 所需资产覆盖汇总。
- AlphaFold 结构补齐尝试：
  - 已写好下载脚本 `src/31_download_dekp_missing_structures.py`。
  - 初始测试发现当前机器访问 AlphaFold 时会出现 SSL 证书校验问题，脚本已增加 `--no-check-certificate` 选项。
  - 单个缺失 UniProt `P00331` 测试返回 AlphaFold HTTP 404，说明部分缺失 ID 可能是旧 accession 或没有当前 AlphaFold 文件，需要做 UniProt accession 映射后再下载。
  - 由于当前网络/审批限制，批量 AlphaFold 下载没有继续完成。
- 当前结论：
  - **DEKP 还没有完成正式预测评估**。原因不是 benchmark 输入缺真值，而是 DEKP 的必要模型/结构资产不完整。
  - 若严格按 off-the-shelf benchmark，当前无法评估 DEKP，因为官方仓库没有提供成品 kcat 权重。
  - 若接受“公开数据重训版 DEKP”，下一步需要补齐：
    - 官方/作者的 kcat 预训练权重，最好是 `Model/Pretrained_model_199_trfm,t5.pkl` 或 final predictor；如果没有，就明确采用从头重训路线。
    - benchmark 余下 301 个 UniProt 的 PDB/AlphaFold 结构文件，放到 `external_methods/DEKP/structures/benchmark/AlphaFold/{UniProt}.pdb`。
    - 或提供 UniProt 旧 accession 到当前 accession 的映射，以便自动从 AlphaFold 重新下载。
    - 若使用完整 DEKP 特征，还需要 PST/MolFormer 预训练模型和对应 feature pickle；若只按作者 fine-tune 默认 `trfm,t5` 配置，则还至少需要生成 `trfm.pkl`、`t5.pkl` 和 `pyg_graph.pkl`。

### 2026-06-22 DEKP-public-retrained 正式评测

- 用户确认采用 **“公开数据重训版 DEKP”**，并要求继续补齐结构文件。
- 结构文件补齐：
  - 解压作者 Zenodo 结构包到 `external_methods/DEKP/structures/author_archive`。
  - 新增 `src/32_collect_dekp_structures.py`，把作者结构包中的 PDB 统一索引，并用符号链接整理到：
    - `external_methods/DEKP/structures/public_kcat/AlphaFold`
    - `external_methods/DEKP/structures/benchmark/AlphaFold`
  - 作者结构包共索引到 4,047 个唯一 PDB ID。
  - 公开 DEKP kcat 训练集 1,850 个 UniProt 全部从作者结构包链接成功。
  - benchmark 初始可从作者结构包链接 194 个 UniProt，仍缺 301 个 UniProt。
  - 修正 `src/31_download_dekp_missing_structures.py`：
    - AlphaFold DB 当前 API 返回 `model_v6`，旧脚本只尝试 v4/v3/v2，导致前期 404。
    - 新脚本优先请求 `https://alphafold.ebi.ac.uk/api/prediction/{uniprot}`，读取 `pdbUrl`，再回退 v6/v5/v4/v3/v2 URL。
  - 批量下载 benchmark 缺失 AlphaFold 结构：
    - 报告：`reports/tables/dekp_alphafold_download_report.csv`
    - 301 个缺失 UniProt 中，296 个新下载，5 个已存在，0 个失败。
  - 重新运行 `run_prepare_dekp_eval.sh` 后：
    - benchmark 总行数 978。
    - 结构可用行 978。
    - SMILES 可用且结构可用行 977。
    - 唯一不可评估行仍是坏 SMILES：`ecoli / Quinate / 192.167`。
- 新增正式训练/预测脚本：
  - `src/33_run_dekp_public_retrained.py`
  - `run_dekp_public_retrained.sbatch`
  - `src/34_evaluate_dekp_predictions.py`
- 训练实现说明：
  - 使用 DEKP 官方 `MetaDecoder` 模型结构和 PDB 几何图特征逻辑。
  - 由于官方仓库没有成品 kcat 权重，也没有完整预计算 feature pickle，本次标注为 `DEKP-public-retrained`。
  - 本次特征模式为 `trfm+sequence_cnn+structure_graph`：
    - `trfm`：复用项目中已有 UniKP SMILES Transformer 权重 `external_methods/CatPred/external/UniKP/trfm_12_23000.pkl` 和 `vocab.pkl`，生成 1024 维 SMILES 特征。
    - `sequence_cnn`：不依赖 ProtT5 tokenizer，使用氨基酸字符 tokenizer 和 DEKP 模型内置 CNN 序列分支训练。
    - `structure_graph`：从 PDB/AlphaFold 坐标生成残基空间图，输入 DEKP 图神经网络分支。
  - 通俗说明：
    - 这不是“作者已训练好模型拿来预测”，而是“用作者公开 kcat 数据重新训练一个 DEKP 架构模型”。
    - 因为没有官方最终权重和完整官方特征缓存，文章中必须写作 `DEKP-public-retrained`，不要与 CatPred/KcatNet 这类直接加载官方权重的方法混为一谈。
- 图结构生成：
  - `enyrnx` 环境没有 `torch_cluster`，原始 `radius_graph` 不能用。
  - 在 `src/33_run_dekp_public_retrained.py` 中用 `scipy.spatial.cKDTree` 生成半径邻接图，绕开 `torch_cluster`。
  - 初次全量作业发现 20 个 PDB 因残基原子不完整导致原始解析函数报 `inhomogeneous shape`。
  - 已实现稳健 PDB 解析：缺失原子用零坐标补齐，侧链缺失时用 CA 或零向量代替。
  - 重跑后图结构全部生成成功：
    - 图缓存：`data/final/dekp/dekp_public_retrained_graph_cache.pkl`，约 9.0 GB。
    - 成功图数：2,199。
    - 失败图数：0。
- 训练数据处理：
  - DEKP 公开 kcat 数据：13,401 行。
  - 清洗后仍为 13,401 行。
  - 与 benchmark 的 `Sequence + canonical SMILES` exact pair 重叠：16 行。
  - 为避免“把测试答案喂给模型”，本次训练排除了这些 exact pair 重叠行。
  - 实际用于训练/验证的数据：13,385 行。
  - 训练/验证划分：
    - train：12,046 行。
    - valid：1,339 行。
  - benchmark 预测行：977 行。
- Slurm 作业：
  - 初次全量作业 `1648705` 完成了特征与图缓存，并暴露出 20 个 PDB 解析失败问题。
  - 修复稳健 PDB 解析后，重新提交作业 `1648719` 到 `qgpu_3090`，在 `gnode8` 完成。
  - 最终训练 early stopping 于第 9 个 epoch。
  - 最佳 epoch：第 4 个。
  - 最佳验证集指标：
    - RMSE log10：1.2800。
    - MAE log10：1.0395。
    - Pearson：0.3217。
- 输出文件：
  - 预测输出：`data/final/dekp/dekp_public_retrained_kcat_input_output.csv`
  - 合并评估表：`data/final/dekp/dekp_public_retrained_kcat_predictions_evaluated.csv`
  - 未预测/无效行：`data/final/dekp/dekp_public_retrained_invalid_or_unpredicted_rows.csv`
  - 未预测汇总：`reports/tables/dekp_public_retrained_invalid_or_unpredicted_summary.csv`
  - 指标表：`reports/tables/dekp_public_retrained_eval_metrics.csv`
  - 运行报告：`reports/tables/dekp_public_retrained_run_report.csv`
  - 模型权重：`data/final/dekp/dekp_public_retrained_model.pt`
  - 跨方法总表：`reports/tables/method_eval_summary.csv`
- DEKP-public-retrained 全量评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 977 | 1.0454 | 1.2458 | 0.2112 | 0.2375 | 0.5056 |
| ecoli | 512 | 1.0425 | 1.2568 | 0.0921 | 0.0721 | 0.5430 |
| yeast | 465 | 1.0486 | 1.2336 | 0.3102 | 0.3675 | 0.4645 |
| BRENDA | 586 | 1.0060 | 1.2636 | 0.1318 | 0.0977 | 0.5904 |
| BRENDA;SABIO-RK | 81 | 0.9628 | 1.1327 | 0.2966 | 0.2681 | 0.5679 |
| SABIO-RK | 310 | 1.1416 | 1.2403 | 0.2699 | 0.3159 | 0.3290 |

- 训练集 exact pair 重叠分组：

| 是否与 DEKP 公开训练集 exact pair 重叠 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| False | 960 | 1.0477 | 1.2497 | 0.2109 | 0.2377 | 0.5021 |
| True | 17 | 0.9142 | 1.0076 | 0.2545 | 0.2663 | 0.7059 |

- 当前跨方法 overall 汇总：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。
  - PMAK：780 行，MAE log10 0.7227，RMSE log10 1.0091，Pearson 0.3828，Spearman 0.4303，10 倍内比例 0.7679。
  - KinForm：563 行，MAE log10 0.7827，RMSE log10 0.9894，Pearson 0.6028，Spearman 0.6297，10 倍内比例 0.7158。
  - KcatNet：977 行，MAE log10 0.6841，RMSE log10 0.9962，Pearson 0.4738，Spearman 0.5149，10 倍内比例 0.7789。
  - PreTKcat：977 行，MAE log10 0.8623，RMSE log10 1.1105，Pearson 0.4087，Spearman 0.4581，10 倍内比例 0.6489。
  - DEKP-public-retrained：977 行，MAE log10 1.0454，RMSE log10 1.2458，Pearson 0.2112，Spearman 0.2375，10 倍内比例 0.5056。
- 结论与注意事项：
  - DEKP-public-retrained 已走通评测流程，覆盖率与 CataPro/KcatNet/PreTKcat 一样为 977/978 行；唯一缺失仍是坏 SMILES。
  - 当前 DEKP-public-retrained 整体表现弱于 KcatNet、CataPro、PreTKcat，也弱于覆盖子集上的 PMAK/KinForm。
  - DEKP-public-retrained 有明显负偏差：bias log10 为 -0.6058，表示整体预测偏低。
  - 该结果不能代表“官方 DEKP 成品模型”的性能，只代表“基于公开 DEKP kcat 数据和当前可复现特征模式重训”的性能。

### 2026-06-22 日报：PreTKcat 与 DEKP 方法评测

今日围绕 kcat benchmark 继续完成 PreTKcat 和 DEKP 两个方法的接入、运行与统一评估。PreTKcat 方面，首先核查官方资源后确认仓库未提供可直接加载的 kcat ExtraTrees 成品回归器，因此采用公开训练数据重训路线；在已有 benchmark 输入基础上，整理生成 PreTKcat 所需的序列、SMILES、温度特征输入，复用本地 ProtT5 和 MolGNet 相关资源提取蛋白与小分子表示，并在 Slurm GPU 队列上完成全量训练和预测。PreTKcat 最终覆盖 977/978 行，唯一缺失项仍为 `ecoli / Quinate / 192.167` 的无效 SMILES；整体 MAE log10 为 0.8623，RMSE log10 为 1.1105，Pearson 为 0.4087，Spearman 为 0.4581，10 倍误差内比例为 0.6489。DEKP 方面，先核查论文、官方 GitHub 和 Zenodo 数据，确认官方未发布最终 kcat 模型权重，且 DEKP 需要额外的蛋白三维结构文件。为完成 `DEKP-public-retrained` 评测，先解压并索引作者结构包，将公开训练集和 benchmark 能命中的结构统一链接到本地结构目录；随后修正 AlphaFold 下载脚本，使其支持当前 AlphaFold DB 的 v6/API 下载规则，成功补齐 benchmark 缺失的 301 个 UniProt 结构文件。之后编写 DEKP 重训脚本，复用 DEKP 的 MetaDecoder 架构，采用 `trfm + sequence_cnn + structure_graph` 特征模式，其中 SMILES Transformer 特征来自本地 UniKP 权重，序列分支使用模型内置 CNN，结构分支由 PDB/AlphaFold 坐标构建残基空间图。由于环境缺少 `torch_cluster`，用 `scipy.spatial.cKDTree` 替代原始 `radius_graph` 生成图边；同时修复少数 PDB 残基原子不完整导致的解析失败问题。最终 DEKP-public-retrained 使用 13,385 行公开训练数据训练，排除了与 benchmark 完全相同的 `Sequence + canonical SMILES` 重叠样本，预测覆盖 977/978 行；整体 MAE log10 为 1.0454，RMSE log10 为 1.2458，Pearson 为 0.2112，Spearman 为 0.2375，10 倍误差内比例为 0.5056。两种方法的结果已写入统一汇总表 `reports/tables/method_eval_summary.csv`，相关预测、评估明细、模型和缓存文件均已保存在 `data/final/pretkcat` 与 `data/final/dekp` 目录中。需要注意的是，PreTKcat 和 DEKP 当前结果都属于“公开数据重训版”或“可复现重训版”，不是直接加载作者发布成品权重的 off-the-shelf 评测，后续论文中应明确标注这一点。

## 2026-06-23 SELFprot 方法评测记录

- 方法来源：
  - 论文：`SELFprot: Effective and Efficient Multitask Finetuning Methods for Protein Parameter Prediction`，JCIM 2025。
  - 官方仓库：`https://github.com/marltanwilson/SELFprot`。
  - 官方权重：Zenodo `https://zenodo.org/records/14266071`，其中 `SELFprot.zip` 包含默认 `models_fold1` 权重。
  - 官方 README 说明其训练数据基于 CatPred-DB split，因此本次把它作为一个可直接加载官方权重的 off-the-shelf 方法评测。
- 输入与权重整理：
  - SELFprot 官方 notebook 的 kcat 预测输入需要两列：蛋白序列 `sequence` 和底物 SMILES `smiles`。
  - 已将 benchmark 输入整理为：
    - `data/final/selfprot/selfprot_kcat_input.csv`
    - `data/final/selfprot/selfprot_kcat_input_valid_smiles.csv`
    - `data/final/selfprot/selfprot_kcat_input_truth.csv`
    - `data/final/selfprot/selfprot_kcat_input_metadata.csv`
  - 全量 978 行中，977 行 SMILES 可解析；唯一不可解析记录仍为 `ecoli / Quinate / 192.167`。
  - 官方权重已解压到 `external_methods/SELFprot/weights/models`，本次默认使用 `models_fold1`。
  - 由于官方 notebook 需要 `facebook/esm2_t12_35M_UR50D` tokenizer，本地已补齐 tokenizer 到 `external_methods/SELFprot/weights/esm2_t12_35M_UR50D_tokenizer`，避免队列任务运行时联网。
- 脚本与运行：
  - 新增 `src/35_prepare_selfprot_eval.py`，用于生成 SELFprot 输入、truth、metadata 和无效 SMILES 记录。
  - 新增 `src/36_run_selfprot_predictions.py`，按官方 notebook 的模型结构批量加载 `chem_model.pt`、`prot_model.pt`、`joint_layer3x.pt`、`kcat_head.pt`、`position_encoding.pt` 等权重，并输出预测。
  - 新增 `src/37_evaluate_selfprot_predictions.py`，将 SELFprot 预测结果与统一 truth 表合并并计算指标。
  - 新增 `run_selfprot_predictions.sbatch`，使用 `enyrnx` 环境和 GPU 队列运行全量预测。
  - 先用 20 条样本做 CPU smoke test，确认输出值范围与 log10(kcat) 尺度一致；随后提交 Slurm 作业 `1650188`，在 `gnode10` 上完成 977 条全量预测。
- 输出文件：
  - 全量预测：`data/final/selfprot/selfprot_kcat_input_output.csv`
  - 评估明细：`data/final/selfprot/selfprot_kcat_predictions_evaluated.csv`
  - 无效/未预测记录：`data/final/selfprot/selfprot_invalid_or_unpredicted_rows.csv`
  - 无效/未预测汇总：`reports/tables/selfprot_invalid_or_unpredicted_summary.csv`
  - 指标表：`reports/tables/selfprot_eval_metrics.csv`
  - 跨方法总表：`reports/tables/method_eval_summary.csv`
- SELFprot 全量评估结果：

| 分组 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---:|---:|---:|---:|---:|---:|
| all | 977 | 0.9541 | 1.2452 | 0.3694 | 0.3706 | 0.6172 |
| ecoli | 512 | 0.9749 | 1.2885 | 0.2975 | 0.2858 | 0.6191 |
| yeast | 465 | 0.9310 | 1.1956 | 0.4389 | 0.4671 | 0.6151 |
| BRENDA | 586 | 0.9225 | 1.2185 | 0.3700 | 0.3521 | 0.6485 |
| BRENDA;SABIO-RK | 81 | 0.8457 | 1.0696 | 0.4119 | 0.3567 | 0.6049 |
| SABIO-RK | 310 | 1.0419 | 1.3345 | 0.3343 | 0.3866 | 0.5613 |

- 当前跨方法 overall 汇总：
  - CatPred：913 行，MAE log10 0.8487，RMSE log10 1.1854，Pearson 0.4051，Spearman 0.4044，10 倍内比例 0.6933。
  - CataPro：977 行，MAE log10 0.7763，RMSE log10 0.9977，Pearson 0.5177，Spearman 0.5158，10 倍内比例 0.7247。
  - PMAK：780 行，MAE log10 0.7227，RMSE log10 1.0091，Pearson 0.3828，Spearman 0.4303，10 倍内比例 0.7679。
  - KinForm：563 行，MAE log10 0.7827，RMSE log10 0.9894，Pearson 0.6028，Spearman 0.6297，10 倍内比例 0.7158。
  - KcatNet：977 行，MAE log10 0.6841，RMSE log10 0.9962，Pearson 0.4738，Spearman 0.5149，10 倍内比例 0.7789。
  - PreTKcat：977 行，MAE log10 0.8623，RMSE log10 1.1105，Pearson 0.4087，Spearman 0.4581，10 倍内比例 0.6489。
  - DEKP-public-retrained：977 行，MAE log10 1.0454，RMSE log10 1.2458，Pearson 0.2112，Spearman 0.2375，10 倍内比例 0.5056。
  - SELFprot：977 行，MAE log10 0.9541，RMSE log10 1.2452，Pearson 0.3694，Spearman 0.3706，10 倍内比例 0.6172。
- 结论与注意事项：
  - SELFprot 已完成官方权重版评测，覆盖 977/978 行；缺失原因不是模型问题，而是原始 benchmark 中 1 条 SMILES 字符串无法被 RDKit 解析。
  - SELFprot 的输出头在官方 notebook 中命名为 `parameter` 和 `SD`，本次根据样本输出范围和 kcat 数据尺度，将 `parameter` 解释为 `prediction_log10`，`SD` 保留为 `selfprot_prediction_sd`。
  - 从当前统一测试集看，SELFprot 的整体 MAE/RMSE 弱于 KcatNet、CataPro、PMAK、KinForm、CatPred 和 PreTKcat，但优于当前 `DEKP-public-retrained`。
  - 由于 SELFprot 使用 CatPred-DB split 相关训练资源，正式论文中建议单独说明其训练数据来源，避免把训练集相近导致的潜在数据重叠问题掩盖掉。

## 2026-06-23 早期四类方法与当前 benchmark 口径核对

- 核对对象：
  - 早期原型方法：`DLKcat`、`MTLKP/MPEK`、`TurNuP`、`UniKP`。
  - 当前正式评测方法：`CatPred`、`CataPro`、`PMAK`、`KinForm`、`KcatNet`、`PreTKcat`、`DEKP-public-retrained`、`SELFprot`。
- 核对结论：
  - 它们都是 kcat 预测相关方法，因此“研究问题”是同一大方向。
  - 但早期四类结果和当前 `method_eval_summary.csv` 还不是完全同一个评测维度，不能直接把早期报告里的指标与当前跨方法总表并列。
- 主要差异：
  - 早期结果主要是基于 *E. coli* 代谢模型的 reaction-level 原型分析，文件包括 `reaction_kcat_MW_DLKcat.csv`、`reaction_kcat_MW_MTLKP.csv`、`reaction_kcat_MW_TurNup.csv`、`reaction_kcat_MW_UniKP.csv` 和 `kcat_comparison_enhanced.csv`。
  - 当前正式评测是以 `experimental_kcat_truth` / `benchmark_ready_truth` 为实验真值，在 *E. coli* + yeast 的 enzyme-substrate 样本上统一计算 MAE、RMSE、Pearson、Spearman 和 10 倍误差内比例。
  - 早期 `kcat_comparison_enhanced.csv` 同时包含 `database`、`fill_method`、GO_HKP 和 ecModel 填充值信息，真值口径不如当前 BRENDA/SABIO-RK curated experimental truth 干净。
  - 早期 TurNuP 更接近 reaction-aware 方法，输入是底物/产物 InChI + enzyme；而当前多数方法是 protein sequence + single substrate SMILES。它可以纳入主文，但需要专门做 reaction 输入适配。
  - 早期 DLKcat、UniKP、MTLKP/MPEK 理论上可以按当前数据重新跑一遍，得到与 CatPred 等方法同口径的结果；不能直接复用早期原型报告中的 reaction-level 统计值。
- 后续建议：
  - 将 DLKcat、UniKP、MTLKP/MPEK、TurNuP 保留为主文候选方法或历史基线。
  - 若要进入最终 `method_eval_summary.csv`，需要用当前 978 条 benchmark 入口重新生成输入、预测和评估。
  - 文章中可把早期结果描述为“方法筛选与原型分析”，正式性能比较只采用当前统一 truth、统一样本和统一指标的结果。

## 2026-06-23 DLKcat、UniKP、MTLKP/MPEK、TurNuP 当前 benchmark 评测记录

- 本次目标：
  - 将一开始原型分析中出现过的 `DLKcat`、`MTLKP/MPEK`、`TurNuP`、`UniKP` 四类方法放到当前统一 benchmark 上重新评估。
  - 当前统一 benchmark 指的是 `data/final/benchmark_ready_catpred.csv` 中的 978 条 enzyme-substrate 实验 kcat 样本，以及对应的 `benchmark_ready_truth`。
- DLKcat：
  - 已下载官方仓库到 `external_methods/DLKcat_official`。
  - 官方 README 说明 DLKcat 的预测输入是蛋白序列和底物 SMILES/名称，输出是 kcat。
  - 已解压官方 `DeeplearningApproach/Data/input.zip`，并加载官方训练好模型：
    - `external_methods/DLKcat_official/DeeplearningApproach/Results/output/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50`
  - 新增 `src/38_run_dlkcat_official.py`，按官方 `prediction_for_input.py` 逻辑批量预测当前 benchmark。
  - 官方模型内部输出为 log2(kcat)，官方脚本会转成线性 kcat；本次最终统一保存为 `prediction_kcat` 和 `prediction_log10`。
  - 预测覆盖 977/978 行；唯一缺失仍是 `ecoli / Quinate / 192.167` 的无效 SMILES。
- UniKP：
  - 官方 README 说明 UniKP 预测 kcat/Km/kcat/Km 时输出的是 log10 转换后的值，需还原为线性值。
  - 已从 HuggingFace `HanselYu/UniKP` 下载官方 kcat 模型：
    - `external_methods/UniKP_official/models/UniKP for kcat.pkl`
  - 新增 `src/39_run_unikp_official_features.py`，使用 UniKP/SMILES Transformer 生成底物特征，并复用本项目 PreTKcat 阶段已缓存的 ProtT5 平均蛋白向量，构建 2048 维 UniKP 特征。
  - 新增 `src/40_predict_unikp_official_py36.py`，用 `condaPY36lin` 环境加载官方 ExtraTrees pickle 并预测。
  - 注意：官方 pickle 是 scikit-learn 0.24.2 保存的，本机 `condaPY36lin` 为 scikit-learn 0.23.1，加载时有版本警告，但模型成功加载并完成预测；默认 Python 的 scikit-learn 1.8.0 无法直接加载该 pickle。
  - 预测覆盖 977/978 行；唯一缺失仍是无效 SMILES。
- MTLKP/MPEK：
  - 当前工作区未找到官方代码目录或可加载的官方权重，只找到早期原型输出 `reaction_kcat_MW_MTLKP.csv`。
  - 新增 `src/42_evaluate_legacy_methods_overlap.py`，用早期 MTLKP 输出在当前 benchmark 中做严格 overlap 评估。
  - 匹配键为 `reaction_id + gene_id + canonical SMILES + sequence`，因此比单纯按 reaction id 更接近 enzyme-substrate 口径。
  - 该结果命名为 `MTLKP-legacy-overlap`，只覆盖 *E. coli* 中能和早期输出精确匹配的 76 条，不是完整官方权重版全量评测。
- TurNuP：
  - 当前工作区未找到可直接加载的官方 TurNuP 成品模型，只找到早期原型输出 `reaction_kcat_MW_TurNup.csv`。
  - TurNuP 早期输出是 reaction-level，只有 `reactions/substrates/products/enzyme/kcat`，不能像 MTLKP 一样精确到当前每个 enzyme-substrate entry。
  - 本次只能按 *E. coli* `reaction_id` 做 legacy overlap 评估，命名为 `TurNuP-legacy-overlap`。
  - 该结果覆盖 338 条当前 benchmark 样本，但它是 reaction-level 映射，不能和 DLKcat-official/UniKP-official 这种 sequence+SMILES 直接预测口径完全等同。
- 新增/修改脚本：
  - `src/38_run_dlkcat_official.py`
  - `src/39_run_unikp_official_features.py`
  - `src/40_predict_unikp_official_py36.py`
  - `src/41_evaluate_method_predictions.py`
  - `src/42_evaluate_legacy_methods_overlap.py`
  - `src/17_build_method_eval_summary.py` 已加入四个新条目。
- 输出文件：
  - DLKcat 预测：`data/final/dlkcat/dlkcat_kcat_input_output.csv`
  - DLKcat 评估明细：`data/final/dlkcat/dlkcat_kcat_predictions_evaluated.csv`
  - DLKcat 指标：`reports/tables/dlkcat_official_eval_metrics.csv`
  - UniKP 预测：`data/final/unikp/unikp_kcat_input_output.csv`
  - UniKP 评估明细：`data/final/unikp/unikp_kcat_predictions_evaluated.csv`
  - UniKP 指标：`reports/tables/unikp_official_eval_metrics.csv`
  - MTLKP overlap 指标：`reports/tables/MTLKP_legacy_overlap_eval_metrics.csv`
  - TurNuP overlap 指标：`reports/tables/TurNuP_legacy_overlap_eval_metrics.csv`
  - 跨方法总表：`reports/tables/method_eval_summary.csv`
- 本次四个方法 overall 指标：

| 方法 | 口径 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---|---:|---:|---:|---:|---:|---:|
| DLKcat-official | 官方权重，sequence + substrate SMILES | 977 | 0.9610 | 1.3208 | 0.3345 | 0.3209 | 0.6070 |
| UniKP-official | 官方 kcat pickle，sequence + substrate SMILES | 977 | 0.8914 | 1.0995 | 0.4780 | 0.5106 | 0.6622 |
| MTLKP-legacy-overlap | 早期 E. coli 输出，精确 overlap | 76 | 0.6545 | 0.8634 | 0.5223 | 0.5133 | 0.8026 |
| TurNuP-legacy-overlap | 早期 E. coli 输出，reaction-level overlap | 338 | 0.7073 | 1.0213 | 0.4683 | 0.5025 | 0.7337 |

- 当前跨方法总表已更新：
  - `DLKcat-official` 和 `UniKP-official` 可以作为当前统一 benchmark 的主表结果。
  - `MTLKP-legacy-overlap` 和 `TurNuP-legacy-overlap` 已并入总表方便查看，但写文章时必须明确它们是 overlap 子集结果，尤其 TurNuP 是 reaction-level 映射，不能与全量 sequence+SMILES 方法作简单排名。
- 如果后续要把 MTLKP/MPEK 和 TurNuP 做成完全同口径评测，需要补充：
  - MTLKP/MPEK 的官方代码仓库、训练好模型权重或可复现推理脚本。
  - TurNuP 的官方代码、模型权重，以及当前 benchmark 中每条记录对应的完整 reaction 输入适配规则。

## 2026-06-23 MTLKP 和 TurNuP 官方源码/权重版统一 benchmark 评测

- 本次目标：
  - 用户补充了 MTLKP 和 TurNuP 的源码及权重，路径分别为：
    - MTLKP：`external_methods/ecm_benchmark_end/etgems_web/script/mtlkp`
    - TurNuP：`external_methods/AI_file/turnup`
  - 因此将之前只能作为历史 overlap 的 `MTLKP/MPEK` 和 `TurNuP`，改为用官方源码/权重在当前统一 benchmark 上重新跑正式评测。
- 数据口径说明：
  - MTLKP 是 `protein sequence + substrate SMILES` 口径，和 DLKcat、UniKP、CataPro、KcatNet、PreTKcat、SELFprot 等更接近。
  - TurNuP 是 `reaction + enzyme` 口径，输入不仅要酶序列，还要反应底物侧和产物侧的分子表示。通俗说，MTLKP 只看“这个酶和这个底物”，TurNuP 要看“这个酶催化的整条反应”。
  - 当前 benchmark 总共 978 条。MTLKP 可覆盖 977 条，缺 1 条 `ecoli / Quinate / 192.167` 非法 SMILES；TurNuP 可覆盖 PMAK 阶段已经整理出完整反应 SMILES 的 780 条，缺 198 条完整反应两侧分子不足的记录。
- 新增/修改脚本：
  - `src/43_prepare_mtlkp_turnup_eval.py`：生成 MTLKP 和 TurNuP 的统一输入、metadata、truth 和缺失记录。
  - `src/44_run_mtlkp_predictions.py`：封装官方 MTLKP `prediction_api.py`，修正硬编码日志路径，并把输出统一成 `prediction_log10` / `prediction_kcat`。
  - `src/45_run_turnup_predictions.py`：封装 TurNuP 官方推理流程，修正本地路径、ESM1b 加载接口差异，并缓存 ESM1b 蛋白表示。
  - `run_prepare_mtlkp_turnup_eval.sh`
  - `run_mtlkp_predict.sh`
  - `run_turnup_predict.sh`
  - `run_mtlkp_full.sbatch`
  - `run_turnup_full.sbatch`
  - `src/17_build_method_eval_summary.py` 已加入 `MTLKP-official` 和 `TurNuP-official`。
  - `src/41_evaluate_method_predictions.py` 增加 MTLKP/TurNuP 的缺失原因识别。
- 环境处理：
  - 在 `enyrnx` 环境中补齐 `xgboost==1.6.1`、`omegaconf==2.3.0`、`unimol-tools==0.1.4.post1` 和 `sentencepiece==0.2.0`。
  - `unimol-tools` 默认会去 HuggingFace 下载 UniMol 权重；为避免联网依赖，已将 MTLKP 包内的 `mol_pre_all_h_220816.pt` 和 `mol.dict.txt` 用符号链接放到 `enyrnx` 的 `unimol_tools/weights/` 目录。
  - TurNuP 的 ESM1b 权重和 XGBoost 权重均来自用户下载目录；ESM1b 蛋白向量缓存写到 `data/final/turnup/turnup_esm1b_cache.pkl`，后续复跑不需要重新算 426 条蛋白表示。
- Slurm 作业：
  - MTLKP：提交 `1650676`，运行成功。
  - TurNuP：第一次 `1650677` 因 fair-esm 加载接口版本差异失败；第二次 `1650678` 完成 ESM1b 缓存后在官方合并函数的 numpy 向量比较处失败；修正后第三次 `1650681` 运行成功。
- 输出文件：
  - MTLKP 输入：`data/final/mtlkp/mtlkp_kcat_input.csv`
  - MTLKP 预测：`data/final/mtlkp/mtlkp_kcat_input_output.csv`
  - MTLKP 评估明细：`data/final/mtlkp/mtlkp_kcat_predictions_evaluated.csv`
  - MTLKP 指标：`reports/tables/mtlkp_eval_metrics.csv`
  - MTLKP 缺失汇总：`reports/tables/mtlkp_missing_summary.csv`
  - TurNuP 输入：`data/final/turnup/turnup_kcat_input.csv`
  - TurNuP 预测：`data/final/turnup/turnup_kcat_input_output.csv`
  - TurNuP 评估明细：`data/final/turnup/turnup_kcat_predictions_evaluated.csv`
  - TurNuP 指标：`reports/tables/turnup_eval_metrics.csv`
  - TurNuP 缺失汇总：`reports/tables/turnup_missing_summary.csv`
  - 跨方法总表：`reports/tables/method_eval_summary.csv`
- 官方权重版 overall 指标：

| 方法 | 口径 | n | MAE log10 | RMSE log10 | Pearson | Spearman | 误差 10 倍内比例 |
|---|---|---:|---:|---:|---:|---:|---:|
| MTLKP-official | 官方权重，sequence + substrate SMILES | 977 | 0.8449 | 1.0941 | 0.4628 | 0.3946 | 0.6418 |
| TurNuP-official | 官方权重，reaction + enzyme | 780 | 0.7009 | 1.0061 | 0.4024 | 0.4224 | 0.7462 |

- 缺失情况：
  - MTLKP：缺 1 条，原因是 `invalid_smiles_192.167`。
  - TurNuP：缺 198 条，原因是 `missing_reaction_smiles`，其中 yeast 153 条、ecoli 45 条。
- 结论与注意事项：
  - `MTLKP-official` 和 `TurNuP-official` 已进入 `reports/tables/method_eval_summary.csv`，可以作为当前统一 benchmark 的正式评测结果。
  - 旧的 `MTLKP-legacy-overlap` 和 `TurNuP-legacy-overlap` 仍保留在总表中，作用是追溯早期原型分析；写文章排名或主表比较时，应优先使用本次官方权重版结果。
  - TurNuP 与 PMAK 一样依赖完整反应信息，因此 n=780；它的指标不能简单理解为“在全部 978 条上覆盖更少但性能更好”，更准确的表述是“在当前已补齐反应 SMILES 的 780 条 reaction-aware 子集上表现为 MAE log10 0.7009”。
  - MTLKP 的 UniMol 阶段有少量分子 3D 构象生成警告，但官方流程仍返回了 977 条预测；这些警告已保留在 `logs/mtlkp_full_1650676.err` 中，便于追溯。
  - TurNuP 的 RDKit stderr 中出现空 `ERROR:` 和氢原子警告，但最终 780 条都完成预测和评估；这些属于 RDKit 解析日志，不是本次评估失败。

## 2026-06-23 全部 kcat 方法评测结果综合分析报告

- 本次目标：
  - 将目前已完成的所有 kcat 方法评测结果汇总成一份可直接用于论文整理的分析报告。
  - 同时生成图表，帮助比较方法误差、相关性、覆盖率、偏差、分物种表现和预测-实验值关系。
- 新增脚本：
  - `src/46_generate_benchmark_report.py`
  - 该脚本读取 `reports/tables/method_eval_summary.csv`、各方法 `*_eval_metrics.csv` 和 `*_predictions_evaluated.csv`，自动生成补充统计表、图表和 Markdown 报告。
- 输出报告：
  - `reports/kcat_benchmark_analysis_report.md`
- 新增统计表：
  - `reports/tables/method_eval_summary_annotated.csv`
  - `reports/tables/method_rank_current_benchmark.csv`
  - `reports/tables/method_rank_all_with_legacy.csv`
  - `reports/tables/species_mae_matrix.csv`
  - `reports/tables/source_database_mae_matrix.csv`
- 新增图表目录：
  - `reports/figures/kcat_benchmark_summary/`
- 已生成图表：
  - `overall_error_mae_rmse.png`：当前正式方法 MAE/RMSE 总体误差。
  - `overall_correlation.png`：Pearson/Spearman 相关性。
  - `coverage_vs_mae.png`：覆盖率与 MAE 的权衡图。
  - `within10_and_bias.png`：10 倍误差内比例与系统性偏差。
  - `error_distribution_boxplot.png`：逐条记录绝对误差分布。
  - `species_mae_heatmap.png`：按物种分层的 MAE 热图。
  - `source_database_mae_heatmap.png`：按数据来源分层的 MAE 热图。
  - `predicted_vs_true_selected.png`：代表性方法的预测值 vs 实验值散点图。
- 报告主要结论：
  - 当前正式评测中，`KcatNet` 是全量/近全量 sequence+SMILES 方法中误差最低的方法，覆盖 977 条，MAE log10 为 0.6841，10 倍误差内比例为 0.7789。
  - `TurNuP-official` 和 `PMAK` 在 reaction-aware 子集上表现较好，但它们只覆盖 780 条已整理出完整 reaction SMILES 的样本。
  - `CataPro` 在全量/近全量方法中表现稳定，覆盖 977 条，MAE log10 为 0.7763。
  - `KinForm` 的 Spearman 相关性最高，为 0.6297，但覆盖 563 条，适合描述为“可覆盖子集上排序能力强”。
  - `DEKP-public-retrained` 是公开数据重训版，不应等同于原论文最优官方模型；当前 MAE log10 为 1.0454，表现较弱。
  - `MTLKP-legacy-overlap` 和 `TurNuP-legacy-overlap` 保留用于追溯早期原型分析，不建议放入主文正式性能排名。
- 核对情况：
  - 重新生成脚本运行成功，无报错。
  - 抽查图表 `coverage_vs_mae.png` 和 `predicted_vs_true_selected.png`，确认图像非空、标签可读。
  - 所有图表尺寸正常，报告中的相对链接可从 `reports/kcat_benchmark_analysis_report.md` 正常指向图表目录。

### 2026-06-23 补充分组标注

- 用户指出“全量/近全量 sequence+SMILES、reaction-aware 子集、模型特定子集、公开数据重训版”这些分组需要在报告中更明确体现。
- 已更新 `src/46_generate_benchmark_report.py`：
  - 为每个方法增加 `group_cn`、`group_note` 和 `coverage_note`。
  - 将 `DEKP-public-retrained` 从普通 sequence+SMILES 方法中单独标为“公开数据重训版”。
  - 新增 `reports/tables/method_group_annotation.csv`，逐方法列出分组、输入口径、覆盖条数和分组理由。
  - 在 `reports/kcat_benchmark_analysis_report.md` 中新增“分组定义与方法归属”小节。
  - 在总体结果表中新增中文分组列 `group_cn` 和输入口径列 `modality`。
  - 在覆盖率-误差图说明中解释颜色与分组的对应关系。
- 当前分组定义：
  - 全量/近全量 sequence+SMILES：DLKcat-official、UniKP-official、MTLKP-official、CataPro、KcatNet、PreTKcat、SELFprot。
  - reaction-aware 子集：TurNuP-official、PMAK。
  - 模型特定子集：CatPred、KinForm。
  - 公开数据重训版：DEKP-public-retrained。
  - 历史 overlap 追溯：MTLKP-legacy-overlap、TurNuP-legacy-overlap。

## 2026-06-23 阶段性总结日报

本阶段围绕 kcat 预测方法统一评测体系，完成了从 benchmark 数据整理、方法接入、模型推理、统一指标计算到综合分析报告生成的完整闭环。首先基于 BRENDA、SABIO-RK、PubChem、CKB compound 映射和已有代谢模型信息，整理出包含实验 kcat 真值、酶序列、底物 SMILES、结构文件及 reaction SMILES 的统一评测数据，并明确了 `experimental_kcat_truth`、`benchmark_ready_truth`、`benchmark_ready_catpred` 以及各方法输入、metadata、truth 文件之间的区别。随后以 CatPred 为首个样板方法跑通评测流程，并依次接入和评估 CataPro、PMAK、KinForm、KcatNet、PreTKcat、DEKP 公开数据重训版、SELFprot、DLKcat-official、UniKP-official、MTLKP-official 和 TurNuP-official，同时保留 MTLKP/TurNuP 的早期 legacy-overlap 结果用于追溯。针对不同方法的输入口径差异，进一步将结果划分为全量/近全量 sequence+SMILES、reaction-aware 子集、模型特定子集、公开数据重训版和历史 overlap 追溯五类，避免不同覆盖范围和输入信息的方法被简单混排。最终，所有方法结果已统一汇总到 `reports/tables/method_eval_summary.csv`，并生成了 `reports/kcat_benchmark_analysis_report.md` 及 8 张配套图表，包括总体误差、相关性、覆盖率-误差权衡、10 倍误差内比例与系统性偏差、误差分布、物种分层、数据来源分层和预测值-实验值散点图。当前正式评测结果显示，KcatNet 在全量/近全量 sequence+SMILES 方法中表现最佳，覆盖 977 条，MAE log10 为 0.6841；TurNuP-official 和 PMAK 在 780 条 reaction-aware 子集上表现较好；KinForm 在其可覆盖子集上 Spearman 相关性最高；DEKP-public-retrained 作为公开数据重训版单独标注，不等同于原论文官方最优模型。整体上，本阶段已形成一套可复用、可追溯、可扩展的 kcat 方法 benchmark 评测框架，为后续论文中方法比较、图表展示和结果讨论奠定了基础。

## 2026-06-24 benchmark_ready_catpred 标准集画像与方法技术比较

- 本次目标：
  - 用户希望补充 `benchmark_ready_catpred` 标准集的详细介绍，包括大肠杆菌和酿酒酵母各自分布、反应/kcat 类型、类似 KEGG 通路的分布特征。
  - 同时希望对不同 kcat 预测方法补充技术原理和比较维度，便于后续写文章时解释“为什么这些方法不能只按一个总排名简单混排”。
- 新增脚本：
  - `src/47_generate_dataset_method_context_report.py`
  - 该脚本读取 `data/final/benchmark_ready_catpred.csv`，并合并 `data/interim/enzyme_reaction_entries.csv`、`data/interim/model_reactions.csv`、`yeast-GEM.xml` 和 DLKcat 官方 `module_ec.txt`，自动生成标准集画像、通路/功能分布、方法技术比较表和图表。
- 新增报告：
  - `reports/kcat_benchmark_dataset_and_method_context.md`
- 新增统计表：
  - `reports/tables/benchmark_ready_catpred_enriched_context.csv`
  - `reports/tables/benchmark_dataset_species_summary.csv`
  - `reports/tables/benchmark_dataset_kcat_stats_by_species.csv`
  - `reports/tables/benchmark_dataset_source_by_species.csv`
  - `reports/tables/benchmark_dataset_match_level_by_species.csv`
  - `reports/tables/benchmark_dataset_enzyme_complex_type_by_species.csv`
  - `reports/tables/benchmark_dataset_substrate_role_by_species.csv`
  - `reports/tables/benchmark_dataset_ec_class_summary.csv`
  - `reports/tables/benchmark_dataset_top_reactions.csv`
  - `reports/tables/benchmark_dataset_top_substrates.csv`
  - `reports/tables/benchmark_dataset_kegg_like_primary_group.csv`
  - `reports/tables/benchmark_dataset_kegg_like_module_membership.csv`
  - `reports/tables/benchmark_dataset_direct_yeast_kegg_pathways.csv`
  - `reports/tables/method_technical_comparison.csv`
  - `reports/tables/method_comparison_dimensions.csv`
- 新增图表目录：
  - `reports/figures/kcat_dataset_context/`
- 已生成图表：
  - `species_distribution.png`：E. coli 与酿酒酵母记录数。
  - `source_by_species.png`：BRENDA/SABIO-RK 来源按物种分布。
  - `kcat_log10_distribution_by_species.png`：实验 kcat 的 log10 分布。
  - `ec_class_distribution.png`：EC 大类分布。
  - `top_reactions.png`：标准集中出现次数最多的反应。
  - `kegg_like_group_by_species.png`：按 EC-to-module 推断的 KEGG-like 功能大类分布。
  - `method_coverage_by_scope.png`：不同方法按比较口径的覆盖率。
- 标准集画像结果：
  - `benchmark_ready_catpred.csv` 当前共 978 条记录。
  - 大肠杆菌 513 条，占 52.5%；酿酒酵母 465 条，占 47.5%。
  - 覆盖 675 个模型反应、495 个基因/UniProt、347 个底物名称和 358 种 EC 注释字符串。
  - 大肠杆菌唯一反应数 451、唯一基因数 327、唯一底物数 232；酿酒酵母唯一反应数 224、唯一基因数 168、唯一底物数 116。
  - 数据来源方面，大肠杆菌主要来自 BRENDA 335 条、SABIO-RK 113 条、二者同时支持 65 条；酿酒酵母来自 BRENDA 252 条、SABIO-RK 197 条、二者同时支持 16 条。
  - 以底物名称粗略识别，343 条记录属于 ATP、NADH、H+、H2O、CoA 等“货币代谢物或辅因子类”分子。通俗说，这些分子是很多反应都会用到的通用小分子，适合统一跑模型，但在生物学解释时不能简单当作每条反应的主要特异性底物。
- kcat 与反应分布：
  - 全部记录的 `true_kcat` 中位数为 19.4 s^-1，log10(kcat) 中位数为 1.288。
  - 大肠杆菌 log10(kcat) 中位数为 1.166，酿酒酵母为 1.602。
  - EC 大类上，大肠杆菌以水解酶、转移酶、氧化还原酶和裂合酶为主；酿酒酵母以氧化还原酶和转移酶为主。
  - 出现次数较多的酵母反应包括 hexokinase、aldehyde dehydrogenase、pyruvate carboxylase、citrate synthase 等；top 底物中 H+、ATP、NADH、NADPH、H2O 等通用分子占比较高。
- 通路/功能注释口径：
  - 第一层使用 DLKcat 官方 `module_ec.txt` 做 EC-to-module 映射，得到跨物种可比较的 KEGG-like 功能大类。它不是直接 KEGG pathway ID，而是按 EC 号推断的功能模块大类。
  - 第二层使用 `yeast-GEM.xml` 中反应自带的 `kegg.pathway` ID，得到酿酒酵母的直接 KEGG pathway 分布。当前 E. coli 的 `eciML1515.json` 有 KEGG reaction ID，但没有直接 pathway 字段，因此 E. coli 暂时只能用 EC-to-module 口径做跨物种比较。
  - EC-to-module 主功能大类中，大肠杆菌约 48.3% 未能映射到精确 module，约 23.8% 属于 primary amino acids/fatty acids/nucleotides，约 11.9% 属于 intermediate；酵母中 secondary_other、intermediate、primary amino acids/fatty acids/nucleotides 和 primary carbohydrate/energy 都有明显覆盖。
  - 酵母直接 KEGG pathway 中 top 项包括 `sce01110`、`sce01130`、`sce00010`、`sce00350`、`sce00071`、`sce01200` 和 `sce01230`。报告中已提醒 `sce01110/sce01130` 这类全局通路覆盖面很广，解释时应更关注具体代谢通路。
- 方法技术比较：
  - `reports/tables/method_technical_comparison.csv` 逐方法整理了技术原理、通俗解释、输入需求、表示学习方式、模型家族、评测口径和注意事项。
  - 明确了 sequence+SMILES 方法、reaction-aware 方法、模型特定子集方法和公开数据重训版方法之间的比较边界。
  - `reports/tables/method_comparison_dimensions.csv` 总结了论文中建议使用的比较维度：输入覆盖、信息粒度、模型来源、评估指标、训练集重叠风险和生物学解释性。
- 核对情况：
  - 运行 `python src/47_generate_dataset_method_context_report.py` 成功。
  - 运行 `python -m py_compile src/47_generate_dataset_method_context_report.py` 成功。
  - 抽查 `reports/kcat_benchmark_dataset_and_method_context.md`，确认物种分布、来源分布、EC 分布、top 反应、top 底物、KEGG-like 分布和方法技术比较均已写入。
  - 图表重跑后已去掉中文字体依赖警告，适合直接纳入后续报告。

## 2026-06-24 增加 GO-HKP 功能相似性 kcat 赋值基线

- 本次目标：
  - 用户希望在已有 AI kcat 预测方法之外，增加一个基于 GO 功能相似性的直接赋值方法 GO-HKP。
  - 目的是检验“如果不训练 AI 回归模型，只按功能层级给 kcat 赋值，是否已经优于 AI 预测”。
- 本地方法目录：
  - `external_methods/GO-HKP/`
  - 其中 `script/GO_Kcat_analysis.py` 和 README 说明了 GO-HKP 的逻辑：先用 DeepGO-SE 为蛋白预测 GO term，再沿 GO 层级在 GO-kcat 统计库中寻找可参考节点，最后给反应赋一个 kcat 统计值。
- 数据核对：
  - 当前可直接复用的本地结果是 `external_methods/GO-HKP/analysis/DeepGO-SE/iML1515R/go_kcat_mean_parent_process_Total_median.json`。
  - 该结果对应 E. coli iML1515/iML1515R 反应体系；未发现可直接对应当前酿酒酵母 benchmark 的 GO-HKP/DeepGO-SE 赋值结果。
  - 本 benchmark 的 E. coli 部分有 513 条、451 个唯一反应，均能在 GO-HKP iML1515R 结果中找到反应级赋值。
  - 酿酒酵母 465 条暂时无法评测，缺失原因统一记录为 `species_without_local_deepgo_se_assignment`。如果后续要补齐 yeast，需要提供或生成 yeast-GEM 对应基因/蛋白的 DeepGO-SE GO 预测，并用 GO-HKP 流程输出反应级 kcat 赋值。
- 新增脚本：
  - `src/48_prepare_go_hkp_eval.py`
  - 该脚本读取 `data/final/benchmark_ready_catpred.csv` 和 GO-HKP 的 E. coli 反应级 JSON 赋值，生成统一评测所需的 input、output、metadata、truth 文件。
- 新增/更新输出：
  - `data/final/go_hkp/go_hkp_kcat_input.csv`
  - `data/final/go_hkp/go_hkp_kcat_input_output.csv`
  - `data/final/go_hkp/go_hkp_kcat_input_metadata.csv`
  - `data/final/go_hkp/go_hkp_kcat_all_metadata.csv`
  - `data/final/go_hkp/go_hkp_kcat_input_truth.csv`
  - `data/final/go_hkp/go_hkp_kcat_predictions_evaluated.csv`
  - `data/final/go_hkp/go_hkp_invalid_or_unpredicted_rows.csv`
  - `reports/tables/go_hkp_eval_readiness.csv`
  - `reports/tables/go_hkp_eval_metrics.csv`
  - `reports/tables/go_hkp_missing_summary.csv`
- 评测覆盖：
  - 总 benchmark：978 条。
  - GO-HKP 可评测：513 条，覆盖率 52.45%，全部来自 E. coli。
  - 缺失：465 条，全部来自 yeast，因为当前缺少本地 yeast 对应的 DeepGO-SE/GO-HKP 赋值。
- 主要指标：
  - MAE log10 = 1.2476。
  - RMSE log10 = 1.5755。
  - Pearson = 0.1045。
  - Spearman = 0.1361。
  - 10 倍误差内比例 = 0.5049。
  - bias log10 = +1.0837，表示整体偏高估。
- 结论解释：
  - GO-HKP 的预测值集中在约 100-200 s^-1 附近，对低 kcat 反应容易明显高估。
  - 在当前 E. coli 子集上，GO 功能层级直接赋值没有优于 KcatNet、CataPro、PMAK、TurNuP 等主要 AI 方法。
  - 这个结果仍然有价值，因为它提供了一个“非 AI、基于功能相似性”的朴素生物学基线。写文章时可以用它说明：仅靠 GO 功能相似性赋参考 kcat 还不足以替代模型预测。
- 已更新的汇总脚本和报告：
  - `src/17_build_method_eval_summary.py`：将 `GO-HKP` 加入统一方法总表。
  - `src/41_evaluate_method_predictions.py`：补充 GO-HKP 缺失原因字段。
  - `src/46_generate_benchmark_report.py`：将 GO-HKP 加入当前正式评测、方法分组、排序表和图表。
  - `src/47_generate_dataset_method_context_report.py`：将 GO-HKP 加入方法技术原理比较表，并新增“AI 预测 vs 直接赋值”比较维度。
- 已重新生成的表格和图表：
  - `reports/tables/method_eval_summary.csv`
  - `reports/tables/method_eval_summary_annotated.csv`
  - `reports/tables/method_group_annotation.csv`
  - `reports/tables/method_rank_current_benchmark.csv`
  - `reports/tables/method_rank_all_with_legacy.csv`
  - `reports/tables/method_technical_comparison.csv`
  - `reports/tables/method_comparison_dimensions.csv`
  - `reports/figures/kcat_benchmark_summary/`
  - `reports/figures/kcat_dataset_context/method_coverage_by_scope.png`
- 已重新生成的报告：
  - `reports/kcat_benchmark_analysis_report.md`
  - `reports/kcat_benchmark_dataset_and_method_context.md`
- 核对情况：
  - 运行 `python src/48_prepare_go_hkp_eval.py` 成功。
  - 运行统一评估脚本 `src/41_evaluate_method_predictions.py` 成功。
  - 运行 `python src/17_build_method_eval_summary.py` 成功。
  - 运行 `python src/46_generate_benchmark_report.py` 成功。
  - 运行 `python src/47_generate_dataset_method_context_report.py` 成功。
  - 运行 `python -m py_compile src/17_build_method_eval_summary.py src/41_evaluate_method_predictions.py src/46_generate_benchmark_report.py src/47_generate_dataset_method_context_report.py src/48_prepare_go_hkp_eval.py` 成功。
  - 抽查 `coverage_vs_mae.png`、`predicted_vs_true_selected.png` 和 `method_coverage_by_scope.png`，确认 GO-HKP 已出现在图中且图像非空。

## 2026-06-24 补齐 GO-HKP 的 yeast 评测

- 本次目标：
  - 用户询问是否能把 GO-HKP 的 yeast 部分也补上。
  - 目标是把 GO 功能相似性赋值基线从只覆盖 E. coli 扩展到当前 978 条统一 benchmark 全覆盖。
- 核对结果：
  - `external_methods/GO-HKP/` 只提供了 iML1515R、iBsu1147R、iCW773R、iDL1450 等本地 DeepGO-SE 预测结果，没有 yeast-GEM 对应的 DeepGO-SE 预测输出。
  - 当前 benchmark 的 yeast 部分有 465 条、168 个唯一 UniProt，`data/interim/uniprot_sequences.csv` 中这 168 个序列全部存在。
  - GO-HKP 自带的四套 DeepGO-SE 预测文件对这 168 个 yeast UniProt 的覆盖为 0，因此无法直接拼接已有 DeepGO-SE 结果。
  - `yeast-GEM.xml` 中有 UniProt 注释，但没有完整 GO 注释；旧 GO-HKP 数据目录中有 `gene2go`，但本地缺少完整 ORF/UniProt 到 NCBI GeneID 的离线映射。
- 采用方案：
  - 使用 UniProt REST 批量接口获取当前 168 个 yeast UniProt 的 GO term 和 GeneID 映射。
  - 将下载结果缓存为 `data/raw/go_hkp/yeast_uniprot_go.tsv`，以后重跑时可离线复用。
  - 由于这条路线使用的是 UniProt GO 注释，而不是 DeepGO-SE 预测，所以在所有表格和报告中明确标注为 `GO-HKP UniProt GO annotation yeast organism-filtered Total median`。
  - E. coli 仍沿用原本的 `GO-HKP DeepGO-SE iML1515R reaction Total median` 结果。
- 新增脚本：
  - `src/49_prepare_go_hkp_with_yeast_uniprot_go.py`
  - 该脚本支持两种运行方式：
    - 加 `--download-uniprot-go`：联网下载/刷新 yeast UniProt-GO TSV。
    - 不加该参数：直接使用本地缓存 TSV 离线重跑。
  - 脚本会读取 GO-HKP 的 `GO_kcat_tree_total.csv` 和 `go-basic.obo`，按 GO 层级为 yeast 基因/蛋白赋 kcat。
- 输出更新：
  - `data/raw/go_hkp/yeast_uniprot_go.tsv`
  - `data/final/go_hkp/go_hkp_kcat_input.csv`
  - `data/final/go_hkp/go_hkp_kcat_input_output.csv`
  - `data/final/go_hkp/go_hkp_kcat_input_metadata.csv`
  - `data/final/go_hkp/go_hkp_kcat_all_metadata.csv`
  - `data/final/go_hkp/go_hkp_kcat_input_truth.csv`
  - `data/final/go_hkp/go_hkp_kcat_predictions_evaluated.csv`
  - `reports/tables/go_hkp_eval_readiness.csv`
  - `reports/tables/go_hkp_eval_metrics.csv`
  - `reports/tables/go_hkp_missing_summary.csv`
- 覆盖变化：
  - 补齐前：GO-HKP 只评测 E. coli 513/978 条，yeast 465 条缺失。
  - 补齐后：GO-HKP 评测 978/978 条，缺失 0 条。
  - E. coli：513/513 ready，451 个唯一反应。
  - yeast：465/465 ready，224 个唯一反应，168 个唯一 UniProt 全部拿到 GO 映射。
- 补齐后的主要指标：
  - overall：MAE log10 = 0.9820，RMSE log10 = 1.3668，Pearson = 0.3499，Spearman = 0.4130，10 倍误差内比例 = 0.6247，bias = +0.8497。
  - E. coli：MAE log10 = 1.2476，Spearman = 0.1361，10 倍误差内比例 = 0.5049。
  - yeast：MAE log10 = 0.6891，Spearman = 0.6351，10 倍误差内比例 = 0.7570。
- 结果解释：
  - 补上 yeast 后，GO-HKP 从子集评测变成全量评测，整体 MAE 明显下降。
  - yeast 的 UniProt GO 注释赋值表现明显好于 E. coli 的本地 DeepGO-SE 反应赋值。
  - 但整体上 GO-HKP 仍弱于 KcatNet、CataPro 等更强的 AI 方法，且 bias 为正，说明整体仍偏高估。
  - 写文章时需要明确：GO-HKP 是非 AI 功能赋值基线；当前 E. coli 和 yeast 的 GO 来源不同，E. coli 为 DeepGO-SE，yeast 为 UniProt GO annotation/GOATOOLS-style 补充。
- 已更新的汇总脚本和报告：
  - `src/46_generate_benchmark_report.py`：将 GO-HKP 从“只覆盖 E. coli 的 GO 赋值子集”更新为“978/978 全覆盖的 GO 功能赋值基线”，并在正文标注 E. coli/yeast 来源差异。
  - `src/47_generate_dataset_method_context_report.py`：方法技术比较表中更新 GO-HKP 的技术原理、输入需求、注意事项和比较维度。
  - `reports/kcat_benchmark_analysis_report.md`
  - `reports/kcat_benchmark_dataset_and_method_context.md`
  - `reports/tables/method_eval_summary.csv`
  - `reports/tables/method_eval_summary_annotated.csv`
  - `reports/tables/method_group_annotation.csv`
  - `reports/tables/method_rank_current_benchmark.csv`
  - `reports/tables/method_technical_comparison.csv`
  - `reports/figures/kcat_benchmark_summary/`
  - `reports/figures/kcat_dataset_context/method_coverage_by_scope.png`
- 核对情况：
  - 运行 `python -m py_compile src/46_generate_benchmark_report.py src/47_generate_dataset_method_context_report.py src/49_prepare_go_hkp_with_yeast_uniprot_go.py` 成功。
  - 运行 `python src/49_prepare_go_hkp_with_yeast_uniprot_go.py --download-uniprot-go` 成功。
  - 运行统一评估脚本 `src/41_evaluate_method_predictions.py` 成功。
  - 运行 `python src/17_build_method_eval_summary.py` 成功。
  - 运行 `python src/46_generate_benchmark_report.py` 成功。
  - 运行 `python src/47_generate_dataset_method_context_report.py` 成功。
  - 抽查 `coverage_vs_mae.png`、`predicted_vs_true_selected.png` 和 `method_coverage_by_scope.png`，确认 GO-HKP 已显示为 978/978 全覆盖。

## 2026-06-24 导出两份报告涉及的表格

- 本次目标：
  - 用户希望将主报告 `reports/kcat_benchmark_analysis_report.md` 和方法技术比较报告 `reports/kcat_benchmark_dataset_and_method_context.md` 中涉及到的表格单独放到一个目录里，便于后续写文章和查找。
- 新增脚本：
  - `src/50_export_report_tables.py`
  - 该脚本会把两份报告用到的 CSV 表格复制到统一导出目录，并额外生成报告正文中动态拼出的精简表。
- 新增导出目录：
  - `reports/report_tables/`
- 目录结构：
  - `reports/report_tables/main_report/`：主评测分析报告相关表。
  - `reports/report_tables/main_report/method_metrics/`：主报告生成汇总和热图时用到的各方法分层指标表。
  - `reports/report_tables/dataset_method_context_report/`：标准集画像与方法技术比较报告相关表。
  - `reports/report_tables/dataset_method_context_report/core_benchmark_tables/`：报告中引用的核心 benchmark 表，包括 `experimental_kcat_truth.csv`、`benchmark_ready_truth.csv` 和 `benchmark_ready_catpred.csv`。
- 新增索引文件：
  - `reports/report_tables/manifest.csv`
  - `reports/report_tables/README.md`
- 导出内容：
  - 共导出 47 个表格文件。
  - 包括 `method_eval_summary.csv`、`method_eval_summary_annotated.csv`、`method_group_annotation.csv`、`method_rank_current_benchmark.csv`、`method_rank_all_with_legacy.csv`、`species_mae_matrix.csv`、`source_database_mae_matrix.csv`、`method_technical_comparison.csv`、`method_comparison_dimensions.csv`、标准集物种/来源/EC/通路/top reaction/top substrate 等统计表。
  - 额外生成了 `main_report_group_definitions.csv`、`main_report_overall_results.csv`、`main_report_legacy_overlap_results.csv` 和 `benchmark_file_role_definitions.csv`，这些对应报告正文中直接展示但原本没有单独 CSV 的表格。
- 核对情况：
  - 运行 `python -m py_compile src/50_export_report_tables.py` 成功。
  - 运行 `python src/50_export_report_tables.py` 成功。
  - 抽查 `reports/report_tables/manifest.csv`、`reports/report_tables/main_report/main_report_overall_results.csv`、`reports/report_tables/dataset_method_context_report/method_technical_comparison.csv` 和 `reports/report_tables/dataset_method_context_report/core_benchmark_tables/benchmark_ready_catpred.csv`，确认文件存在且行数正常。

## 2026-06-29 补全数据获取/处理方法与项目目录说明

- 本次目标：
  - 将 benchmark 的数据来源、获取、清洗、匹配和筛选方法写入 `reports/kcat_benchmark_dataset_and_method_context.md`。
  - 补充反应类型分布、GO/KEGG-like/yeast 直接 KEGG pathway 注释口径。
  - 把项目目录按初始 kcat 分析、GO、KEGG、MAE、species-level、method-level、结构补齐、方法输入输出和发布资产等用途分类。
- 报告生成脚本更新：
  - `src/47_generate_dataset_method_context_report.py`
  - 新增 `make_benchmark_build_funnel()`，导出 `reports/tables/benchmark_build_funnel.csv`。
  - 新增 `make_project_directory_map()`，导出 `reports/tables/project_directory_analysis_map.csv`。
- benchmark 构建方法在报告中明确为：
  - 从 `eciML1515.json` 和 `yeast-GEM.xml` 解析反应、GPR、EC、UniProt 和候选底物。
  - UniProt REST 补序列；模型注释、CKB、MetaNetX、PubChem 补底物 SMILES；reaction SMILES 和蛋白结构只作为方法特定输入。
  - 主实验真值只采用 BRENDA turnover number 和 SABIO-RK kcat；正值统一为 `s^-1`，BRENDA 默认排除 mutant/variant。
  - 匹配优先级：`species+EC+UniProt+substrate ID` > `species+EC+substrate ID` > `species+EC+UniProt+substrate name` > `species+EC+substrate name`。
  - 同一 entry 只保留最高匹配层级，在原始 kcat 尺度取中位数，再计算 log10。同步将 `configs/matching_rules.yaml` 从错误的 `median_log10_kcat` 改为 `median_raw_kcat_then_log10`。
  - `experimental_kcat_truth.csv` 为 1072 条实验真值全集；要求单蛋白序列和底物 SMILES 后得到 978 条统一 benchmark。
- benchmark 数量漏斗：
  - E. coli：5883 个模型反应，6115 个 enzyme-substrate entries，554 条实验真值，513 条 benchmark，451 个唯一反应，327 个基因。
  - yeast：4131 个模型反应，7816 个 entries，518 条实验真值，465 条 benchmark，224 个唯一反应，168 个基因。
- GO/KEGG 口径：
  - E. coli GO-HKP 使用 iML1515R DeepGO-SE 反应级赋值。
  - yeast GO-HKP 使用 UniProt GO annotation，再沿 GO OBO 层级匹配 GO-kcat 统计值。
  - 跨物种通路画像使用 EC-to-module 的 KEGG-like 分类；yeast 额外解析 yeast-GEM 自带 `kegg.pathway`，E. coli 因模型没有同口径 pathway 字段，不做直接 KEGG pathway 强行对比。
- 表格导出更新：
  - `src/50_export_report_tables.py` 加入上述两个新表。
  - 重新导出 `reports/report_tables/`，目前共 49 个报告表格。
- 核对：
  - `python -m py_compile src/47_generate_dataset_method_context_report.py src/50_export_report_tables.py` 成功。
  - `python src/47_generate_dataset_method_context_report.py` 成功。
  - `python src/50_export_report_tables.py` 成功。

## 2026-06-29 GitHub 与 Zenodo 公开发布准备

- 公开仓库目标：
  - `https://github.com/dengxiao01/kcat_benchmark_analysis.git`
  - 远端核对时仓库可访问但没有 `HEAD`，判断为尚无提交的空仓库。
- 新增/更新公开文件：
  - 重写根目录 `README.md`，说明 978 条双物种标准集、三类真值文件、方法分组、统一指标、完整流水线、目录结构、GitHub/Zenodo 分工和复现命令。
  - 新增 `.gitignore`，排除 token、`.env`、原始/中间数据、方法缓存、第三方嵌套仓库、模型和压缩包。
  - 新增 `requirements.txt`、`LICENSE`、`CITATION.cff`、`THIRD_PARTY_NOTICES.md`。
  - 新增 `external_methods/METHOD_SOURCES.md`，记录 13 个方法的上游地址、评测 commit、本地目标目录、模型包和观察到的上游许可。
  - BRENDA 官方许可页核对为 CC BY 4.0；Zenodo 混合许可记录采用词表中的 `other-open`，第三方模型仍服从各自上游条款。
- 新增发布工具：
  - `scripts/prepare_zenodo_release.py`：建立模型/结果归档并计算 SHA256。
  - `scripts/download_zenodo_assets.py`：按公开 manifest 下载、校验和恢复资产。
  - `scripts/publish_zenodo_release.py`：从本地 `zenodo.txt` 读取 token，创建草稿、上传、远端按文件名/字节数核对并发布；token 不写入状态文件和 Git。
- Zenodo 资产设计最终调整：
  - 初稿曾包含完整 ProtT5、基础 ESM1b 和 CatPred 全论文复现 capsule，共 36.89 GiB。
  - 核对实际推理依赖后，改为只发布任务专用权重；通用 ProtT5/ESM1b 从其官方源下载，路径写入 `METHOD_SOURCES.md`。
  - CatPred 只从 9.5 GiB capsule 中选择 `data/pretrained/production/kcat` 的 39 个条目，约 90 MiB，不再复制全部消融/复现实验。
  - 最终 7 个资产总计约 9.678 GiB：
    - `kcat_benchmark_core_data_and_results.tar.gz`：约 0.024 GiB。
    - `kcat_benchmark_catpred_kcat_assets.tar.gz`：约 0.088 GiB。
    - `catpred_db.tar.gz`：约 1.066 GiB。
    - `kinform_results.tar.gz`：约 1.546 GiB。
    - `kcat_benchmark_mtlkp_kcat_assets.tar`：约 3.421 GiB，不含 ProtT5-XL。
    - `kcat_benchmark_turnup_kcat_assets.tar`：约 2.563 GiB，不含基础 ESM1b checkpoint。
    - `kcat_benchmark_other_model_assets.tar`：约 0.970 GiB。
  - 机器可读大小、SHA256、恢复动作和目标目录写入 `zenodo_assets_manifest.csv`。
- 敏感信息与 Git 检查：
  - 确认排除 `zenodo.txt`、`external_methods/CatPred/.env.vercel`、`hf_DHRT.../`、`release/`、`data/raw/` 和第三方 `.git`。
  - 待发布 Git 清单 316 个文件，最大单文件为 `yeast-GEM.xml` 约 11.7 MiB，没有超过 GitHub 100 MiB 限制的文件。
  - 对待发布文件扫描 `ZENODO_ACCESS_TOKEN`、`MODAL_TOKEN`、`HF_TOKEN` 和常见 access-token 赋值模式，无敏感值命中。
  - 已在项目根目录执行 `git init -b main` 和 `git add .`；尚未 commit/push，等待 Zenodo DOI 和下载 URL 写回 README/manifest 后统一提交。
- Zenodo 当前远端状态：
  - 已成功创建草稿 deposition `21024684`。
  - 精简前的核心结果包曾成功上传，草稿尚未发布；精简后的核心包字节数不同，恢复上传时会自动替换。
  - 再次发起外部上传时，执行环境提示“外部操作审批额度已用尽”，并明确禁止绕过；因此当前未继续上传，也未发布记录。
- 下一步恢复命令：
  - 审批恢复并得到用户明确同意后，从项目根目录运行 `python -u scripts/publish_zenodo_release.py --publish`。
  - 脚本会读取现有 `release/zenodo/zenodo_state.json`，继续草稿 `21024684`，按文件从小到大上传，全部远端核验通过后发布并自动回写 README DOI 与 `zenodo_assets_manifest.csv` 下载链接。
  - 随后重新 `git add .`、提交并推送到目标 GitHub 仓库，最后用 GitHub/Zenodo 公共链接做远端复核。

## 2026-06-29 Zenodo 正式发布与分片上传补充

- 用户明确同意继续上传 Zenodo 并推送 GitHub 后，恢复草稿 `21024684`。
- 大文件直接 PUT 的问题：
  - Python `requests` 上传近 1 GiB 单文件时，在约 0.47 GiB 处收到 Zenodo nginx `502 Bad Gateway`。
  - 为避免整包重传，将超过 128 MiB 的逻辑包切成固定大小分片；7 个逻辑包最终对应 81 个 Zenodo 文件，低于 Zenodo 每记录 100 文件限制。
  - `zenodo_assets_manifest.csv` 新增 `bundle_name`、`part_index`、`part_count`、`bundle_size_bytes` 和 `bundle_sha256`，同时保留每片 SHA256。
  - `scripts/download_zenodo_assets.py` 更新为按逻辑包展示和下载，自动按 `part_index` 拼接，并在恢复前同时校验分片和整包哈希。
- 上传实现调整：
  - `scripts/publish_zenodo_release.py` 改用 Zenodo 官方示例推荐的 `curl --upload-file`。
  - Zenodo token 通过权限受限的临时 curl 配置传入，不出现在进程命令行、日志、manifest 或 Git 中。
  - Python `requests` 仍用于草稿元数据、远端文件清单核对和发布动作；TLS 校验始终开启。
  - 81 个远端文件均按文件名和字节数核对通过，随后发布成功。
- 正式 Zenodo 信息：
  - record：`https://zenodo.org/records/21024684`
  - DOI：`10.5281/zenodo.21024684`
  - 逻辑包：7 个。
  - 远端文件：81 个。
  - 总大小：10,392,182,635 bytes，约 9.678 GiB。
- 已自动更新：
  - `README.md` 的 DOI badge、记录链接和下载说明。
  - `zenodo_assets_manifest.csv` 的 `zenodo_record_id` 与逐文件 `download_url`。
  - `CITATION.cff` 的版本、发布日期、DOI 和 Zenodo URL。
- 上传完成后曾因执行环境外部审批额度到限而暂缓公共 API 抽查；审批恢复后访问 `https://zenodo.org/api/records/21024684` 返回 HTTP 200，确认 DOI、标题、81 个文件和 10,392,182,635 bytes 均与本地 manifest 一致。
