# benchmark_ready_catpred 标准集数据画像与方法技术比较

生成时间：2026-06-30 10:25

## 1. 先说结论

`benchmark_ready_catpred.csv` 当前有 978 条可评测记录，其中大肠杆菌 513 条（52.5%），酿酒酵母 465 条（47.5%）。它覆盖 675 个模型反应、495 个基因/UniProt、347 个底物名称和 358 种 EC 注释字符串。

这个文件名里的 `catpred` 是历史命名：CatPred 是第一个被打通的评测对象，但该文件现在不是 CatPred 专用数据，而是当前统一 benchmark 的 sequence+SMILES 标准集。

按底物名称粗略识别，343 条记录的待预测底物属于 ATP/NADH/H+/H2O/CoA 等“货币代谢物或辅因子类”分子。它们适合用于模型输入统一评测，但写生物学解释时应和真正的主底物区分开。

## 2. 文件定位与字段含义

- 标准集：`data/final/benchmark_ready_catpred.csv`
- 行粒度：一行是一个 `酶/基因组反应/候选底物/实验 kcat` 记录。
- 关键输入字段：`sequence` 是酶蛋白序列，`SMILES` 是底物结构字符串，`reaction_id` 是模型反应 ID。
- 关键真值字段：`true_kcat` 是实验 kcat 原值，`true_kcat_log10` 是 log10 变换后的真值，评估主指标都在 log10 空间计算。
- 溯源字段：`source_database`、`match_level`、`reference`、`n_measurements` 用于说明数据来自 BRENDA/SABIO-RK 以及匹配依据。

输入、metadata、truth 的关系可以这样理解：`*_input.csv` 只给模型看；`*_metadata.csv` 记录每一行的物种、反应、来源和处理状态；`*_truth.csv` 只在评估时使用，避免把答案混进模型输入。

## 3. 数据获取、清洗与 benchmark 确定方法

### 3.1 从代谢模型定义候选酶-反应-底物条目

- 大肠杆菌使用项目根目录的 `eciML1515.json`，酿酒酵母使用 `yeast-GEM.xml`。`src/01_parse_models.py` 解析反应、方向、GPR、EC、UniProt 和代谢物数据库编号。
- GPR 是 gene-protein-reaction 规则，通俗说就是一条反应由哪些基因编码的酶负责。脚本把 `or` 拆成同工酶候选，把 `and` 保留为多亚基复合物。
- 每行候选 entry 的粒度是 `物种 + 模型反应 + GPR 基因组 + 候选底物`。优先选择非辅因子反应物；如果没有，再退回全部反应物，因此 ATP/NADH/H2O 等通用分子仍可能出现在标准集。

### 3.2 补齐蛋白序列、小分子结构和反应结构

- 蛋白序列：根据模型中的 UniProt accession，通过 UniProt REST 批量获取，并缓存为 `data/raw/uniprot_sequences.fasta` 和 `data/interim/uniprot_sequences.csv`。
- 底物 SMILES：先使用模型已有注释，再通过 BiGG、KEGG、ChEBI、MetaNetX 等交叉编号在 CKB compound 数据库中映射；仍缺失时调用 PubChem PUG REST，并记录查询缓存和无法映射原因。
- 完整 reaction SMILES：只为 TurNuP/PMAK 等 reaction-aware 方法准备，不作为所有方法进入统一 benchmark 的硬条件。
- 蛋白结构：只为 DEKP 等结构感知方法收集，当前优先使用 AlphaFold/本地结构缓存，也不是统一 sequence+SMILES benchmark 的硬条件。

### 3.3 获取并整理实验 kcat 真值

- 主真值只使用 BRENDA turnover number 和 SABIO-RK kcat。早期推断或数据库填充得到的数值不进入统一 benchmark，也不作为当前公开仓库产物。
- 仅保留目标物种、正的 kcat，统一单位为 `s^-1`；BRENDA 默认排除注释为 mutant/mutation/variant 的记录。范围值取区间均值。
- 匹配先限定 `species + EC`，再比较底物数据库 ID/规范化名称以及 UniProt。优先级从高到低为：`species_ec_uniprot_substrate_id`、`species_ec_substrate_id`、`species_ec_uniprot_substrate_name`、`species_ec_substrate_name`。
- 同一 entry 只保留最高匹配层级的实验记录；多条实验值在 kcat 原始尺度取中位数，再计算 `log10(kcat)`。pH 和温度也取可用记录的中位数，并保留来源、参考文献和测量条数。

### 3.4 确定最终 benchmark

- `experimental_kcat_truth.csv` 是匹配到模型 entry 的实验真值全集，共 1072 行。
- `benchmark_ready_truth.csv` 进一步要求 entry 能进入统一模型输入，即有单蛋白序列和可用底物 SMILES，共 978 行。
- `benchmark_ready_catpred.csv` 在这 978 行上合并 sequence、SMILES、真值和溯源字段。文件名保留 `catpred` 只是因为 CatPred 是第一个打通的方法，并不表示该标准集只服务于 CatPred。
- 方法评测时从这个母表提取各自需要的输入列，真值列只在推理结束后用于评分，避免答案泄漏到模型输入。

从模型到最终 benchmark 的数量漏斗如下。注意 `enzyme_substrate_entries` 可以多于模型反应数，因为一条反应可能拆成多个基因组和多个候选底物。

| species | model_total_reactions | reactions_with_gpr | reactions_with_ec | enzyme_substrate_entries | entries_with_uniprot_sequence | entries_with_substrate_smiles | experimental_truth_rows | benchmark_ready_rows | benchmark_unique_reactions | benchmark_unique_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ecoli | 5883 | 5402 | 1479 | 6115 | 5446 | 5789 | 554 | 513 | 451 | 327 |
| yeast | 4131 | 2709 | 2436 | 7816 | 7181 | 7071 | 518 | 465 | 224 | 168 |

## 4. 物种、实验来源与匹配层级

| species | rows | percent | unique_reactions | unique_genes | unique_substrates | median_log10_kcat | pH_available_rows | temperature_available_rows | currency_or_cofactor_like_rows_by_name |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ecoli | 513 | 52.454 | 451 | 327 | 232 | 1.166 | 435 | 422 | 81 |
| yeast | 465 | 47.546 | 224 | 168 | 116 | 1.602 | 423 | 408 | 262 |

实验来源按物种分布如下：

| species | source_database | rows | percent_of_benchmark |
| --- | --- | --- | --- |
| ecoli | BRENDA | 335 | 34.254 |
| yeast | BRENDA | 252 | 25.767 |
| yeast | SABIO-RK | 197 | 20.143 |
| ecoli | SABIO-RK | 113 | 11.554 |
| ecoli | BRENDA;SABIO-RK | 65 | 6.646 |
| yeast | BRENDA;SABIO-RK | 16 | 1.636 |

匹配层级分布如下，`species_ec_uniprot_substrate_id` 通常代表物种、EC、UniProt 和底物 ID 都能对上，是最严格的一类匹配：

| species | match_level | rows | percent_of_benchmark |
| --- | --- | --- | --- |
| yeast | species_ec_substrate_id | 338 | 34.56 |
| ecoli | species_ec_substrate_id | 290 | 29.652 |
| ecoli | species_ec_uniprot_substrate_id | 220 | 22.495 |
| yeast | species_ec_uniprot_substrate_id | 127 | 12.986 |
| ecoli | species_ec_substrate_name | 3 | 0.307 |

相关图：

- `reports/figures/kcat_dataset_context/species_distribution.png`
- `reports/figures/kcat_dataset_context/source_by_species.png`
- `reports/figures/kcat_dataset_context/kcat_log10_distribution_by_species.png`

## 5. kcat 数值范围与反应分布

kcat 跨度非常大，因此评估使用 log10 空间。`log10(kcat)=0` 表示 1 s^-1，`log10(kcat)=2` 表示 100 s^-1，`log10(kcat)=-2` 表示 0.01 s^-1。

| group | n | true_kcat_min | true_kcat_q25 | true_kcat_median | true_kcat_mean | true_kcat_q75 | true_kcat_max | log10_kcat_min | log10_kcat_q25 | log10_kcat_median | log10_kcat_mean | log10_kcat_q75 | log10_kcat_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | 978 | 1.67e-04 | 3.32 | 19.4 | 1,444.14 | 102.65 | 570,000.00 | -3.777 | 0.521 | 1.288 | 1.198 | 2.011 | 5.756 |
| ecoli | 513 | 0.003 | 2.4 | 14.65 | 2,633.83 | 64 | 570,000.00 | -2.523 | 0.38 | 1.166 | 1.054 | 1.806 | 5.756 |
| yeast | 465 | 1.67e-04 | 5.6 | 40 | 131.639 | 176 | 3,830.00 | -3.777 | 0.748 | 1.602 | 1.357 | 2.246 | 3.583 |

EC 大类分布如下。这里按 EC membership 统计：一条记录如果有多个 EC，可能贡献到多个类别。

| species | ec_class | row_memberships | percent_of_benchmark |
| --- | --- | --- | --- |
| ecoli | 3 Hydrolases / 水解酶 | 169 | 17.28 |
| ecoli | 2 Transferases / 转移酶 | 165 | 16.871 |
| ecoli | 1 Oxidoreductases / 氧化还原酶 | 72 | 7.362 |
| ecoli | 4 Lyases / 裂合酶 | 71 | 7.26 |
| ecoli | 6 Ligases / 连接酶 | 25 | 2.556 |
| ecoli | 5 Isomerases / 异构酶 | 20 | 2.045 |
| yeast | 1 Oxidoreductases / 氧化还原酶 | 187 | 19.121 |
| yeast | 2 Transferases / 转移酶 | 167 | 17.076 |
| yeast | 3 Hydrolases / 水解酶 | 37 | 3.783 |
| yeast | 4 Lyases / 裂合酶 | 37 | 3.783 |
| yeast | 6 Ligases / 连接酶 | 32 | 3.272 |
| yeast | 5 Isomerases / 异构酶 | 9 | 0.92 |

出现次数最多的反应记录如下：

| species | reaction_id | reaction_name | rows | unique_substrates | unique_genes | median_log10_kcat |
| --- | --- | --- | --- | --- | --- | --- |
| yeast | r_0534 | hexokinase (D-glucose:ATP) | 10 | 2 | 5 | 1.009 |
| yeast | r_4652 | aldehyde dehydrogenase (1-propanol, NAD) | 9 | 3 | 3 | 2.246 |
| yeast | r_4654 | aldehyde dehydrogenase (1-propanol, NADP) | 9 | 3 | 3 | 1.842 |
| yeast | r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP) | 9 | 3 | 3 | 1.852 |
| yeast | r_4661 | aldehyde dehydrogenase (methionol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_0958 | pyruvate carboxylase | 6 | 3 | 2 | 1.778 |
| yeast | r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol) | 6 | 3 | 2 | 2.246 |
| yeast | r_4672 | aldehyde dehydrogenase (tyrosol, NADP) | 6 | 2 | 3 | 1.803 |
| yeast | r_0186 | aldehyde dehydrogenase (tryptophol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_0543 | homocitrate synthase | 6 | 3 | 2 | -0.347 |
| yeast | r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_4670 | aldehyde dehydrogenase (tyrosol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_0165 | mitochondrial alcohol dehydrogenase | 6 | 3 | 2 | 2.246 |
| yeast | r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD) | 6 | 2 | 3 | 2.246 |
| yeast | r_4663 | aldehyde dehydrogenase (methionol, NADP) | 6 | 2 | 3 | 1.803 |
| yeast | r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP) | 6 | 3 | 2 | 1.852 |
| yeast | r_0300 | citrate synthase | 6 | 3 | 2 | 1.045 |
| yeast | r_4653 | aldehyde dehydrogenase (1-propanol, NAD) | 6 | 3 | 2 | 2.246 |

出现次数最多的底物名称如下：

| substrate_name | rows | species_count | unique_reactions | currency_or_cofactor_like_by_name | median_kcat |
| --- | --- | --- | --- | --- | --- |
| H+ | 64 | 1 | 27 | True | 176 |
| ATP | 62 | 1 | 44 | True | 2.77 |
| NADH | 52 | 1 | 22 | True | 176 |
| NADPH | 23 | 1 | 12 | True | 56.8 |
| H2O H2O | 18 | 1 | 18 | True | 1.107 |
| CTP | 15 | 1 | 15 | True | 8.92 |
| 2-Oxoglutarate | 13 | 1 | 13 | False | 30.65 |
| L-Aspartate | 13 | 1 | 13 | False | 1.14 |
| pyruvate | 10 | 1 | 4 | False | 65 |
| H2O | 9 | 1 | 6 | True | 8.2 |
| acetyl-CoA | 9 | 1 | 6 | True | 0.45 |
| GTP C10H12N5O14P3 | 9 | 1 | 9 | True | 14.65 |
| GMP C10H12N5O8P | 9 | 1 | 9 | True | 11 |
| L-Glutamate | 8 | 1 | 8 | False | 21.925 |
| propanal | 8 | 1 | 3 | False | 0.73 |
| 2-oxoglutarate | 8 | 1 | 5 | False | 3.345 |
| polyphosphate | 8 | 1 | 4 | True | 1,350.00 |
| Malonyl-[acyl-carrier protein] | 8 | 1 | 8 | False | 0.011 |
| NAD | 8 | 1 | 5 | True | 7.4 |
| AMP C10H12N5O7P | 7 | 1 | 7 | True | 386.05 |

底物角色粗分布如下：

| species | currency_or_cofactor_like_by_name | rows | percent_of_benchmark |
| --- | --- | --- | --- |
| ecoli | False | 432 | 44.172 |
| yeast | True | 262 | 26.789 |
| yeast | False | 203 | 20.757 |
| ecoli | True | 81 | 8.282 |

相关图：

- `reports/figures/kcat_dataset_context/ec_class_distribution.png`
- `reports/figures/kcat_dataset_context/top_reactions.png`

## 6. GO、KEGG-like 与通路/功能注释

这里把 GO 和 KEGG 分开使用，避免把两类功能注释混成同一个概念：GO 更接近“蛋白/基因做什么功能”，KEGG 更接近“反应位于哪类代谢通路”。

- GO-HKP 功能赋值：E. coli 使用 GO-HKP 自带的 iML1515R DeepGO-SE 反应级结果；yeast 使用 UniProt GO 注释补齐。两者都沿 `go-basic.obo` 的 GO 层级在 `GO_kcat_tree_total.csv` 中寻找可参考 kcat，并取 Total median 作为赋值。由于两物种 GO 来源不同，报告和 metadata 中分别标注来源。
- 跨物种 KEGG-like 注释：使用 `DLKcat_official/DeeplearningApproach/Data/subsystem/module_ec.txt` 的 EC-to-module 功能大类。它不是直接 KEGG pathway ID，但可用同一口径比较 E. coli 和 yeast。
- yeast 直接 KEGG pathway：解析 `yeast-GEM.xml` 中反应自带的 `kegg.pathway` ID。E. coli 的 `eciML1515.json` 只有 KEGG reaction ID，没有系统的 pathway 字段，因此当前不能用同样方式做直接 pathway 统计。
- 这些注释用于描述标准集覆盖和构建 GO 赋值基线，不参与其他 AI 模型的真值筛选，也不会改变实验 kcat。

按 EC 推断的 KEGG-like 主功能大类如下，每行只归到一个主类，因此合计等于各物种样本数：

| species | kegg_like_primary_group_short | rows | percent_of_species | percent_of_benchmark |
| --- | --- | --- | --- | --- |
| ecoli | Unmapped | 248 | 48.343 | 25.358 |
| ecoli | Primary AA/FA/nt | 122 | 23.782 | 12.474 |
| ecoli | Intermediate | 61 | 11.891 | 6.237 |
| ecoli | Secondary other | 50 | 9.747 | 5.112 |
| ecoli | Primary carbohydrate/energy | 24 | 4.678 | 2.454 |
| ecoli | Secondary | 7 | 1.365 | 0.716 |
| ecoli | Unclassified module | 1 | 0.195 | 0.102 |
| yeast | Secondary other | 120 | 25.806 | 12.27 |
| yeast | Unmapped | 108 | 23.226 | 11.043 |
| yeast | Intermediate | 92 | 19.785 | 9.407 |
| yeast | Primary AA/FA/nt | 76 | 16.344 | 7.771 |
| yeast | Primary carbohydrate/energy | 64 | 13.763 | 6.544 |
| yeast | Secondary | 3 | 0.645 | 0.307 |
| yeast | Unclassified module | 2 | 0.43 | 0.204 |

yeast-GEM 直接 KEGG pathway ID 的 top 分布如下。`sce01100/sce01110/sce01130` 这类全局通路会覆盖很多反应，解释时应更关注具体代谢通路，例如碳代谢、氨基酸生物合成、嘌呤/嘧啶代谢等。

| kegg_pathway_id | kegg_pathway_name | rows | unique_reactions | percent_of_benchmark |
| --- | --- | --- | --- | --- |
| sce01110 | Biosynthesis of secondary metabolites | 275 | 96 | 28.119 |
| sce01130 | Biosynthesis of antibiotics | 274 | 82 | 28.016 |
| sce00010 | Glycolysis / Gluconeogenesis | 198 | 44 | 20.245 |
| sce00350 | Tyrosine metabolism | 127 | 26 | 12.986 |
| sce00071 | Fatty acid degradation | 127 | 25 | 12.986 |
| sce01200 | Carbon metabolism | 97 | 26 | 9.918 |
| sce01230 | Biosynthesis of amino acids | 87 | 37 | 8.896 |
| sce04070 |  | 46 | 38 | 4.703 |
| sce00620 | Pyruvate metabolism | 37 | 13 | 3.783 |
| sce01210 | 2-Oxocarboxylic acid metabolism | 37 | 14 | 3.783 |
| sce00562 | Inositol phosphate metabolism | 33 | 25 | 3.374 |
| sce00561 | Glycerolipid metabolism | 28 | 6 | 2.863 |
| sce00040 | Pentose and glucuronate interconversions | 26 | 6 | 2.658 |
| sce00520 | Amino sugar and nucleotide sugar metabolism | 25 | 9 | 2.556 |
| sce00300 | Lysine biosynthesis | 25 | 9 | 2.556 |
| sce00052 | Galactose metabolism | 21 | 7 | 2.147 |
| sce00270 | Cysteine and methionine metabolism | 21 | 12 | 2.147 |
| sce00970 |  | 21 | 16 | 2.147 |
| sce00260 | Glycine, serine and threonine metabolism | 19 | 11 | 1.943 |
| sce00051 | Fructose and mannose metabolism | 19 | 6 | 1.943 |

相关图：

- `reports/figures/kcat_dataset_context/kegg_like_group_by_species.png`

## 7. 项目目录结构与分析类型

下面按“目录承担什么工作”整理当前公开项目结构。统一 benchmark 位于 `data/final/`，方法级结果位于 `data/final/<method>/`，汇总表和图位于 `reports/`，可执行入口统一放在 `scripts/runners/`。

| category | path | directory_type | contents |
| --- | --- | --- | --- |
| Benchmark construction | src/01_parse_models.py to src/11_finalize_benchmark_data.py; configs/ | code and rules | GEM parsing, GPR/EC/substrate extraction, UniProt sequence retrieval, SMILES mapping, BRENDA/SABIO-RK truth matching, and final benchmark filtering. |
| Raw source data | data/raw/ | large source/cache data | BRENDA, SABIO-RK, compound/CKB, UniProt FASTA, GO mappings, and method source assets. Large files are distributed through Zenodo rather than Git. |
| Intermediate curation | data/interim/ | rebuildable intermediate data | Reaction-entry tables, sequence/SMILES queues, caches, review lists, reaction SMILES, and method input preparation tables. |
| Unified benchmark and method outputs | data/final/ | final data products | Experimental truth, benchmark-ready tables, and per-method inputs, metadata, predictions, missing rows, structures, and evaluated rows. |
| Pipeline runners | scripts/runners/ | local and Slurm entry points | Portable shell entry points for benchmark construction, method input preparation, prediction, and full cluster jobs. Run from the repository root after activating the required environment. |
| GO analysis | external_methods/GO-HKP/; data/raw/go_hkp/; data/final/go_hkp/ | functional-assignment analysis | GO hierarchy and GO-kcat resources, E. coli DeepGO-SE assignments, yeast UniProt GO mappings, GO-HKP evaluated rows, readiness, and species-level metrics. |
| KEGG/EC/pathway analysis | src/47_generate_dataset_method_context_report.py; reports/tables/benchmark_dataset_kegg* | functional distribution analysis | EC-to-module KEGG-like groups across species and direct yeast-GEM KEGG pathway annotations. |
| MAE and error analyses | reports/tables/*_eval_metrics.csv; reports/figures/kcat_benchmark_summary/ | performance analysis | Overall and grouped MAE/RMSE, correlation, bias, within-fold error, coverage-error tradeoff, error distributions, and predicted-versus-true plots. |
| Species-level analysis | reports/tables/species_mae_matrix.csv; reports/tables/benchmark_dataset_*_by_species.csv | species-stratified analysis | E. coli versus yeast counts, truth distributions, sources, matching levels, pathway groups, and method MAE. |
| Method-level analysis | reports/tables/method_*.csv; data/final/<method>/ | method comparison | Method principles, input requirements, comparison groups, coverage, rankings, evaluated predictions, and method-specific limitations. |
| Reports and publication tables | reports/; reports/report_tables/; docs/ | human-readable and manuscript material | Main analysis report, dataset/method context report, figures, standalone report tables, work log, and manuscript assets. |
| Third-party methods and model assets | external_methods/ | third-party code and large assets | Published method source code, checkpoints, model bundles, dependency snapshots, and caches. Only lightweight reproducibility code belongs in Git; large assets belong in Zenodo. |

重点分析文件可以快速定位为：

- 数据准备、方法输入、预测和 Slurm 队列入口：`scripts/runners/`。
- GO 分析：`external_methods/GO-HKP/`、`data/raw/go_hkp/`、`data/final/go_hkp/`、`reports/tables/go_hkp_*`。
- KEGG/EC/通路分析：`reports/tables/benchmark_dataset_kegg_like_*`、`benchmark_dataset_direct_yeast_kegg_pathways.csv`。
- MAE/RMSE/bias/within-fold 分析：`reports/tables/*_eval_metrics.csv` 和 `reports/figures/kcat_benchmark_summary/`。
- Species-level 分析：`species_mae_matrix.csv`、`benchmark_dataset_*_by_species.csv` 和 species heatmap。
- Method-level 分析：`method_eval_summary*.csv`、`method_rank*.csv`、`method_technical_comparison.csv`。
- 写文章用独立表格：`reports/report_tables/`，其中 `manifest.csv` 记录每张表的来源。

## 8. 不同预测方法的技术原理与比较维度

下面这张表把“模型看了什么信息”和“它适合在哪个维度比较”放在一起。通俗地说，sequence+SMILES 方法是只看酶和单个底物；reaction-aware 方法还看产物，因此信息更多但需要更完整的数据；结构感知方法还看蛋白结构，但公开复现难度也更高；GO-HKP 则不是 AI 回归模型，而是用 GO 功能相似性做直接 kcat 赋值。本项目里 GO-HKP 的 E. coli 部分来自本地 DeepGO-SE 反应赋值，yeast 部分用 UniProt GO 注释补齐。

| method | group_cn | n | coverage_percent | input_needed | model_family | plain_language_cn | main_caveat |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GO-HKP | 功能相似性 GO 赋值基线 | 978 | 100 | protein/gene GO terms or DeepGO-SE GO predictions + reaction/gene mapping | functional-similarity assignment baseline | 它不是训练一个 AI 回归模型，而是找功能相似的酶/GO 节点，把已有 kcat 作为参考值直接赋给目标反应或基因。 | 当前已覆盖 978/978，但 E. coli 与 yeast 的 GO 来源不同；yeast 是 UniProt GO 注释路线，不是 DeepGO-SE 预测路线。 |
| UniKP-official | 全量/近全量 sequence+SMILES | 977 | 99.898 | enzyme sequence + substrate SMILES | pretrained embeddings + machine learning | 先让两个预训练模型分别读懂蛋白和小分子，再把特征交给回归器。 | 依赖预计算特征和模型版本；不显式看产物侧。 |
| MTLKP-official | 全量/近全量 sequence+SMILES | 977 | 99.898 | enzyme sequence + substrate SMILES | deep multitask network, BAN/MMoE-style fusion | 把蛋白和底物都变成高维向量，并让模型学习哪些残基和哪些原子更相关。 | 官方权重直接推理；仍是单底物视角。 |
| CataPro | 全量/近全量 sequence+SMILES | 977 | 99.898 | enzyme sequence + substrate SMILES + wild/mutant type | pretrained embeddings + neural regression | 同时让模型读懂酶序列、小分子字符串和传统分子指纹。 | 需要外部预训练权重；不显式建模产物侧。 |
| DLKcat-official | 全量/近全量 sequence+SMILES | 977 | 99.898 | enzyme sequence + substrate name/SMILES | deep learning, GNN + CNN | 看一个底物长什么样、酶序列长什么样，然后学习二者组合对应的转换速度。 | 较早一代模型；不显式使用完整反应物/产物信息。 |
| DEKP-public-retrained | 公开数据重训版 | 977 | 99.898 | sequence + SMILES + protein structure/graph assets | structure-aware deep learning | 它比普通 sequence+SMILES 方法多看蛋白结构相关信息，但我们没有原论文最优私有权重，因此重新训练了公开可复现版本。 | 不是官方最优权重；结果反映当前公开复现流程，不等同于论文声称上限。 |
| SELFprot | 全量/近全量 sequence+SMILES | 977 | 99.898 | sequence + SMILES | protein-ligand deep model | 把酶和小分子看成一对相互作用对象，学习它们共同决定的动力学数值。 | README 信息较简略；需要在论文写作时继续补齐正式引用与模型细节。 |
| PreTKcat | 全量/近全量 sequence+SMILES | 977 | 99.898 | sequence + SMILES + temperature | pretrained embeddings + ExtraTrees | 除了酶和底物，还把实验温度作为影响 kcat 的因素放进去。 | 缺失温度需要填补默认值；本项目使用公开数据可复现流程。 |
| KcatNet | 全量/近全量 sequence+SMILES | 977 | 99.898 | enzyme sequence + substrate SMILES | geometric/deep learning | 把蛋白和底物都编码成结构化特征，再用深度网络学习它们的匹配关系。 | 序列会按模型规则截断；不显式使用产物侧。 |
| CatPred | 模型特定子集 | 913 | 93.354 | SMILES + sequence + unique pdbpath/protein record | deep learning ensemble | 这是一个面向多种酶动力学参数的统一深度学习框架，本项目只取其中 kcat 模型来评测。 | 官方流程有额外可处理范围限制；当前覆盖 913/978。 |
| TurNuP-official | reaction-aware 子集 | 780 | 79.755 | reactant SMILES + product SMILES + enzyme sequence | reaction-aware ML, XGBoost | 不仅看底物，还看反应前后分子怎么变，再结合酶序列预测速度。 | 必须有完整反应 SMILES；当前只覆盖 780/978。 |
| PMAK | reaction-aware 子集 | 780 | 79.755 | reaction SMILES + enzyme sequence | reaction-aware deep learning | 让模型同时看酶、反应，并尝试关注对催化更关键的残基位置。 | 依赖完整反应 SMILES；当前覆盖 780/978。 |
| KinForm | 模型特定子集 | 563 | 57.566 | sequence + SMILES, plus cached embeddings/assets | embedding-based deep/ML models | 尽量把整条蛋白和潜在结合位点的信息都编码进去，再和底物信息融合。 | 受官方 Zenodo bundle/缓存资产覆盖限制；当前覆盖 563/978。 |

比较维度建议如下：

| dimension | plain_explanation_cn | use_in_paper |
| --- | --- | --- |
| 输入覆盖 | 方法能不能吃下标准集的 978 行。缺一行通常是 SMILES 非法；缺更多则说明方法需要额外信息或官方资产不全。 | 报告覆盖率 n/978，并按全量、reaction-aware、模型特定子集、公开重训版分开解释。 |
| 信息粒度 | 只看单个底物，还是看完整反应，或者还看蛋白结构。 | sequence+SMILES 方法可互相直接比较；reaction-aware 方法需要单独说明其信息更多但覆盖更窄。 |
| 模型来源 | 是官方权重直接推理，还是我们用公开数据重训/复现。 | DEKP-public-retrained 不能和官方最优权重画等号，应标成公开可复现版本。 |
| AI 预测 vs 直接赋值 | AI 模型会从输入特征中学习连续 kcat 数值；GO-HKP 这类方法则用功能相似性把已有 kcat 统计值赋给目标反应。 | GO-HKP 可作为非 AI 生物学基线，单独回答“简单功能赋值是否已经优于 AI 预测”。 |
| 评估指标 | MAE/RMSE 看误差大小，Pearson/Spearman 看相关性，within10 看是否落在 10 倍误差内。 | 主表同时给覆盖率和误差，避免只看某一个指标。 |
| 训练集重叠风险 | 如果测试样本和方法训练集重叠，指标可能虚高。 | 后续写文章时应继续按序列、SMILES、sequence-SMILES pair 做查重标注。 |
| 生物学解释性 | 能否解释到反应、残基、底物原子或通路层面。 | PMAK/TurNuP 更适合讨论反应变化；结构方法可讨论结构资产，但需谨慎。 |

相关图：

- `reports/figures/kcat_dataset_context/method_coverage_by_scope.png`

## 9. 写文章时建议怎么表述

1. 主文可以把 `benchmark_ready_catpred.csv` 称为 unified sequence+SMILES benchmark，而不是 CatPred 专用输入。
2. 结果表必须同时给覆盖率和误差指标；否则 KinForm、CatPred、PMAK/TurNuP 这种子集方法会和全量方法混在一起，容易误导。
3. 通路分布建议分两句话写：跨物种用 EC-to-module 的 KEGG-like 功能大类；酿酒酵母另有直接 KEGG pathway ID 作为补充。
4. 如果要做更严格的生物学解释，下一步应补充 E. coli 的 KEGG reaction-to-pathway 映射，或用 BioCyc/MetaCyc 子系统统一重注释两套模型。

## 10. 输出文件清单

- `reports/tables/benchmark_dataset_species_summary.csv`
- `reports/tables/benchmark_dataset_kcat_stats_by_species.csv`
- `reports/tables/benchmark_dataset_source_by_species.csv`
- `reports/tables/benchmark_dataset_match_level_by_species.csv`
- `reports/tables/benchmark_dataset_ec_class_summary.csv`
- `reports/tables/benchmark_dataset_top_reactions.csv`
- `reports/tables/benchmark_dataset_top_substrates.csv`
- `reports/tables/benchmark_dataset_kegg_like_primary_group.csv`
- `reports/tables/benchmark_dataset_direct_yeast_kegg_pathways.csv`
- `reports/tables/method_technical_comparison.csv`
- `reports/tables/method_comparison_dimensions.csv`
- `reports/tables/benchmark_build_funnel.csv`
- `reports/tables/project_directory_analysis_map.csv`
- `reports/report_tables/`
