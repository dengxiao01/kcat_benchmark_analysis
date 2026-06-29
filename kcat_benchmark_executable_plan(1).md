# AI kcat 预测方法系统评测项目：可执行方案

## 核对后的执行调整

根据当前目录中的实际数据和脚本，本计划总体方向正确，但执行顺序需要更稳一点：

1. 先落地 **Phase 1 的本地可执行部分**：从 `eciML1515.json` 和 `yeast-GEM.xml` 中抽取 reaction、GPR、gene、UniProt、EC、底物候选信息，形成统一的 enzyme-reaction/substrate entry 表。通俗说，就是先把“要评测的对象清单”整理干净。
2. `reaction_kcat_MW_databasefill.csv` 只能作为 legacy sanity check。它里面有 Database、fill、Other_species 等混合来源，不能直接当作主文里的实验真值。
3. 主真值仍然应来自 BRENDA/SABIO-RK，并且必须保留匹配等级、实验条件和文献来源。当前工作区没有这些原始下载文件或 API 凭据，所以本轮先生成匹配所需的输入表和空 schema。
4. “热力学不对称性”建议在后续脚本和图题中改成 **directional sensitivity / reversible-pair directional asymmetry**。只看 kcat 的正反向差异，不能完整判断热力学一致性；完整判断还需要 Km 和 Keq。
5. 新方法优先级建议改为：CatPred、CataPro 先接入；PMAK 作为第三批但要先确认是否有稳定预训练推理入口；KinForm/KcatNet 放在后续扩展，因为它们依赖外部模型资产或更重的特征生成流程。
6. 第一版文章应明确定位为 **off-the-shelf predictor benchmark**，即“把已有方法当作现成预测器来测”，暂不承诺统一重训所有模型。

## 0. 项目目标

当前项目已经完成了一个以 *Escherichia coli* 代谢模型为基础的 kcat 预测方法比较原型，主要比较了 DLKcat、MTLKP/MPEK、TurNuP 和 UniKP 四类方法，并生成了 HTML 报告、PPT、统计表和图。现在的目标不是继续美化原始报告，而是把它升级成一个可以投稿的小文章：

> 构建一个跨物种、跨方法、以实验 kcat 为真值的 AI kcat 预测方法系统 benchmark。

建议文章定位为：

**Cross-species benchmark of AI-predicted enzyme turnover numbers against curated experimental kcat values**

或者：

**Benchmarking AI-predicted enzyme turnover numbers across prokaryotic and eukaryotic metabolic networks**

中文逻辑是：

> 本研究构建了一个以原核和真核模式生物代谢网络为背景的 kcat 预测评测框架，系统比较现有 AI kcat 预测方法在覆盖率、准确性、一致性、系统性偏差、跨物种泛化和适用范围方面的差异。

---

## 1. 当前项目现状

已有项目结构大致如下：

```text
kcat_benchmark/
├── kcat预测方法评估报告.pptx
└── kcat_benchmark_analysis/
    ├── eciML1515.json
    ├── GO_HKP.json
    ├── reaction_kcat_MW_DLKcat.csv
    ├── reaction_kcat_MW_MTLKP.csv
    ├── reaction_kcat_MW_TurNup.csv
    ├── reaction_kcat_MW_UniKP.csv
    ├── reaction_kcat_MW_databasefill.csv
    ├── kcat_comparison.csv
    ├── kcat_comparison_with_gpr.csv
    ├── kcat_comparison_enhanced.csv
    ├── integrated_kcat_simple.json
    ├── scripts/
    │   ├── 01_integrate_kcat.py
    │   ├── 02_create_kcat_csv.py
    │   ├── 03_filter_kcat_add_gpr.py
    │   ├── 04_add_database_fill.py
    │   ├── 05_comprehensive_analysis.py
    │   └── run_all.sh
    └── analysis_results/
        ├── 01_global_distribution.png
        ├── 02_correlation_analysis.png
        ├── 03_thermodynamic_asymmetry.png
        ├── 04_isozyme_specificity.png
        ├── 05_complex_handling.png
        ├── 06_substrate_specificity.png
        ├── 07_coverage_analysis.png
        ├── 08_ground_truth_benchmark.png
        ├── 09_bias_detection.png
        ├── 10_ensemble_comparison.png
        └── corresponding statistics csv files
```

当前项目已经能完成以下事情：

1. 合并 DLKcat、MTLKP/MPEK、TurNuP、UniKP、GO_HKP 的 kcat 预测结果；
2. 读取 `eciML1515.json`，加入 reaction-level GPR 信息；
3. 合并 `database` 和 `fill_method` 信息；
4. 生成增强版比较表 `kcat_comparison_enhanced.csv`；
5. 生成 10 类统计分析图；
6. 输出 HTML 报告和 PPT。

当前报告的分析维度包括：

1. 全局分布；
2. 方法间一致性；
3. 正反向反应方向敏感性；
4. 同工酶特异性；
5. 酶复合物处理能力；
6. 底物特异性；
7. 覆盖范围；
8. 真值基准测试；
9. 系统性偏差；
10. 简单集成建模。

这些分析方向是有价值的，可以保留一部分作为最终文章的结果框架。但当前版本还不能直接投稿，因为真值来源、物种范围和方法覆盖都需要升级。

---

## 2. 当前版本的主要问题

### 2.1 物种单一

当前只分析了 *E. coli*，即便结果完整，也容易被审稿人认为只是某个模型或某个物种上的现象。建议至少加入一个真核模式生物：

| 类型 | 物种 | 推荐模型 |
|---|---|---|
| 原核 | *Escherichia coli* | iML1515 / eciML1515 |
| 真核 | *Saccharomyces cerevisiae* | yeast-GEM / Yeast8 / Yeast9 |

这样文章可以从“单个模型分析”升级成“原核—真核跨物种 benchmark”。

### 2.2 真值来源不够干净

当前分析中使用了 `reaction_kcat_MW_databasefill.csv` 和 ecModel 相关 kcat 信息。这些值可能混合了：

1. 实验测定值；
2. BRENDA/SABIO-RK 数据库值；
3. 同源推断值；
4. 模型填充值；
5. 人工默认值；
6. 异常填充值，例如当前报告中识别出的 `112645.641475366`。

投稿时不能把这些混合来源都称为“真实 kcat”。最终主 benchmark 必须回到原始实验数据逻辑：

```text
代谢模型反应
→ GPR / gene / enzyme / EC / substrate
→ 通过 EC + organism + substrate 从 BRENDA/SABIO-RK 获取实验 kcat
→ 构建 curated experimental kcat truth
→ 与 AI 预测值比较
```

`ecYeast` 或 ecModel 里的 kcat 可以作为辅助对照，但不建议作为主真值。

### 2.3 方法偏旧或不完整

当前项目只有 4 个主要方法：

```text
DLKcat
MTLKP/MPEK
TurNuP
UniKP
```

这不够支撑“现有 AI kcat 方法系统评测”的标题。至少应增加：

```text
CatPred
CataPro
PMAK
KinForm
PreTKcat
KcatNet
DEKP
SELFprot
```

其中第一批最值得接入的是：

```text
CatPred
CataPro
PMAK
```

这三个方法能显著提升文章的新意和可信度。

### 2.4 热力学表述需要更严谨

当前报告中把正反向 kcat 是否相同解释为“是否违背热力学基本原理”。这个表述偏强。严格来说，Haldane 关系需要同时考虑：

```text
Keq
kcat_forward
kcat_reverse
Km_forward
Km_reverse
```

所以最终文章中建议把当前的“热力学不对称性”改成：

```text
directional asymmetry
reversible-pair directional sensitivity
```

更稳妥的表述是：

> 部分方法对正反向反应给出近似相同的 kcat，提示其方向敏感性不足。但单独 kcat 不能完整判断热力学一致性，完整判断需要结合 Km 和 Keq。

### 2.5 训练数据泄漏风险

很多 kcat 预测方法本身就是用 BRENDA、SABIO-RK 或其整理版本训练的。如果我们再用 BRENDA/SABIO-RK 数据评测，审稿人一定会问：

> 这些测试样本是否已经出现在模型训练集中？

因此最终文章至少需要设置：

```text
Benchmark A: all curated experimental kcat
Benchmark B: remove exact duplicated enzyme-substrate pairs from available training datasets
Benchmark C: low-sequence-identity subset, e.g. <40% or <60% to training sequences
```

如果某些方法没有公开训练集，就在文章中诚实说明：

> For methods without accessible training data, potential train-test overlap was noted as a limitation.

---

## 3. 推荐纳入的 AI kcat 预测方法

### 3.1 主评测方法

| 方法 | 是否已有 | 是否纳入主文 | 作用 |
|---|---:|---:|---|
| DLKcat | 已有 | 是 | 经典深度学习 kcat 预测方法，适合作为历史基线 |
| UniKP | 已有 | 是 | 统一预测 kcat、Km、kcat/Km 的重要方法 |
| TurNuP | 已有 | 是 | reaction-aware / natural reaction kcat 预测代表方法 |
| MPEK/MTLKP | 已有 | 是 | 多任务 kcat/Km 预测路线 |
| CatPred | 未有 | 强烈建议 | 当前最适合做主 benchmark 的新方法之一 |
| CataPro | 未有 | 强烈建议 | 高影响力方法，支持 kcat/Km/kcat/Km |
| PMAK | 未有 | 强烈建议 | 引入 reaction SMILES 和 residue-aware attention |
| KinForm | 未有 | 建议 | 强调 kcat/KM、特征优化和低相似序列泛化 |
| KcatNet | 未有 | 建议 | genome-wide kcat prediction，适合跨物种分析 |
| PreTKcat | 未有 | 建议 | 轻量 baseline，复现成本较低 |
| DEKP | 未有 | 补充 | 多模态/结构信息方法 |
| SELFprot | 未有 | 补充 | 多任务 protein parameter prediction |

### 3.2 不建议作为主评测的方法

| 方法 | 原因 |
|---|---|
| ESP | 预测 enzyme-substrate pair，不是 kcat 数值，不适合放入 kcat 回归主 benchmark |
| UniKineG | 需要 docking 和 3D complex，工程量大，适合后续扩展 |
| TCNeKP | 如果代码和数据不稳定，先只在综述表中列出 |
| GELKcat / Methods interpretability framework | 更偏解释性，若无完整代码不适合主 benchmark |
| NNKcat | 若无稳定代码，作为补充讨论即可 |

---

## 4. 文献和代码拿到后，我们到底要做什么？

这是项目最关键的地方。

拿到这些推荐文献和 GitHub 以后，**不是简单复现原文所有实验**。我们的目标不是证明 CatPred、CataPro 或 PMAK 的原文结果是否正确，而是让它们在我们的统一 benchmark 数据集上产生预测值，然后与我们整理的实验 kcat 真值进行公平比较。

换句话说，我们要做的是：

```text
不是：完整复现每篇文章的训练、消融、案例、湿实验
而是：把每个方法当成一个预测器，对我们构建的 enzyme-reaction entries 统一跑推理
```

具体分三种层级。

### 4.1 最低层级：只跑预训练模型推理

这是最优先、最现实的方式。

流程：

```text
我们的数据表
→ 转成该方法要求的输入格式
→ 调用作者提供的预训练模型或 Web/API/脚本
→ 得到 pred_kcat
→ 转回统一输出格式
→ 加入 benchmark_long.csv
```

适合：

```text
CatPred
CataPro
UniKP
DLKcat
PMAK
PreTKcat
```

优点：

1. 成本最低；
2. 最接近普通用户使用这些方法的真实场景；
3. 不需要重新训练；
4. 适合小文章；
5. 便于比较“现有工具实际可用性”。

缺点：

1. 可能受到训练数据泄漏影响；
2. 不同方法训练集不同；
3. 作者预训练模型可能不是完全同一数据体系。

### 4.2 中等层级：在统一数据集上重新训练

流程：

```text
我们的 curated experimental kcat dataset
→ 统一 train/valid/test split
→ 用各方法代码重训
→ 在同一测试集评估
```

优点：

1. 最公平；
2. 最适合做严肃 benchmark；
3. 可以控制训练数据泄漏。

缺点：

1. 工程量非常大；
2. 每个方法环境不同；
3. 有些方法不提供完整训练代码；
4. 算力消耗较高；
5. 可能偏离“用户直接使用工具”的实际场景。

### 4.3 最高层级：重新训练 + 消融 + OOD 泛化

这相当于一篇更大的方法学 benchmark 文章。

需要做：

```text
random split
enzyme-unseen split
substrate-unseen split
EC-unseen split
species-transfer split
sequence identity split
training set deduplication
uncertainty calibration
high-kcat underestimation analysis
```

这适合做大文章，但不适合当前“小文章快速发出”的目标。

### 4.4 推荐策略

建议采用两阶段策略：

```text
第一版小文章：
以“预训练模型统一推理 + curated experimental truth”为主。
重点评估各方法作为现有工具时的实际表现。

后续扩展文章：
再做统一重训、OOD split 和严格去泄漏 benchmark。
```

第一版文章可以明确写：

> We evaluated the off-the-shelf performance of existing AI kcat predictors on a curated cross-species enzyme-reaction benchmark.

这句话很重要，因为它把研究范围限定为“现有工具即插即用表现”，避免审稿人要求你完整重训所有方法。

---

## 5. 以 CatPred 为例：拿到文献和 GitHub 后如何得到我们项目需要的数据？

### 5.1 CatPred 在项目中的作用

CatPred 是当前最值得优先接入的方法之一。它的价值有三点：

1. 预测对象直接包括 kcat；
2. 有论文、GitHub、Web app 和本地运行脚本；
3. 方法定位就是 in vitro enzyme kinetic parameter prediction，和我们的项目高度匹配。

CatPred 文章标题为：

```text
CatPred: a comprehensive framework for deep learning in vitro enzyme kinetic parameters
```

论文页面：

```text
https://www.nature.com/articles/s41467-025-57215-9
```

GitHub：

```text
https://github.com/maranasgroup/CatPred
```

CatPred 官方网站：

```text
https://www.catpred.com/
```

GitHub README 显示 CatPred 支持 kcat/Km 和 Ki 预测，并支持 CSV import/export 批量工作流。仓库中也提供了本地/demo 推理脚本，例如 `demo_run.py`，其用法是：

```bash
python demo_run.py --parameter <kcat|km|ki> --input_file <path_to_input_csv>
```

注意：实际字段名、环境和模型权重路径需要以 GitHub 当前 README 和脚本说明为准。

### 5.2 我们需要 CatPred 输出什么？

我们需要的是：

```text
每一个 enzyme-substrate 或 enzyme-reaction entry 的 predicted kcat
```

统一输出格式应该是：

```csv
species,reaction_id,gene_id,ec_number,substrate_name,substrate_smiles,protein_sequence,method,pred_kcat,pred_log10_kcat,status,error_message
ecoli,R_PFK,b3916,2.7.1.11,D-fructose 6-phosphate,SMILES...,MXXXX,CatPred,123.4,2.091,success,
```

最终这些结果会并入：

```text
data/final/benchmark_long.csv
```

### 5.3 CatPred 输入数据从哪里来？

CatPred 需要的核心输入通常是：

```text
protein sequence
substrate SMILES
parameter = kcat
```

因此我们必须先从模型和数据库构建一张基础表：

```csv
species,reaction_id,gene_id,ec_number,substrate_name,substrate_smiles,protein_sequence
ecoli,R_PFK,b3916,2.7.1.11,D-fructose 6-phosphate,SMILES...,MXXXX
```

这张表不是 CatPred 给我们的，而是我们自己从代谢模型和注释数据库整理出来的。

来源如下：

| 字段 | 来源 |
|---|---|
| species | 当前模型，例如 ecoli / yeast |
| reaction_id | GEM 模型反应 ID |
| gene_id | GEM 的 GPR |
| ec_number | 模型注释 / UniProt / Rhea / KEGG / MetaNetX |
| substrate_name | 反应底物 |
| substrate_smiles | 模型代谢物注释、BiGG、MetaNetX、PubChem |
| protein_sequence | UniProt 或模型基因注释映射 |

### 5.4 CatPred 数据准备步骤

#### Step 1：生成 enzyme-reaction entries

从 E. coli 和 yeast 模型中提取：

```text
reaction_id
reaction equation
GPR
gene list
EC number
metabolites
candidate substrates
```

输出：

```text
data/interim/enzyme_reaction_entries.csv
```

示例：

```csv
species,reaction_id,gene_id,ec_number,substrate_name,metabolite_id
ecoli,R_PFK,b3916,2.7.1.11,D-fructose 6-phosphate,f6p_c
yeast,R_PFK,YGR240C,2.7.1.11,D-fructose 6-phosphate,s_0775[c]
```

#### Step 2：补 protein sequence

通过 gene → UniProt → sequence 获取蛋白序列。

输出：

```text
data/interim/enzyme_reaction_entries_with_sequence.csv
```

示例：

```csv
species,reaction_id,gene_id,uniprot_id,protein_sequence
ecoli,R_PFK,b3916,P0A796,MKRIAVLTSGGDAPGMNAAIR...
```

#### Step 3：补 substrate SMILES

通过 metabolite ID → BiGG / MetaNetX / PubChem 获取 SMILES。

输出：

```text
data/interim/enzyme_reaction_entries_with_smiles.csv
```

示例：

```csv
species,reaction_id,substrate_name,substrate_smiles
ecoli,R_PFK,D-fructose 6-phosphate,O[C@H]1...
```

#### Step 4：生成 CatPred 输入文件

写一个 adapter：

```text
src/predictors/prepare_catpred_input.py
```

输入：

```text
data/interim/enzyme_reaction_entries_with_sequence_smiles.csv
```

输出：

```text
data/interim/prediction_inputs/catpred_kcat_input.csv
```

字段需要按 CatPred 当前脚本要求确定。逻辑上至少包括：

```csv
id,sequence,smiles,parameter
ecoli|R_PFK|b3916|f6p,MKRIAVLTSGGDAPGMNAAIR...,O[C@H]1...,kcat
```

如果 CatPred 要求多底物输入或 primary substrate 标记，则需要额外字段，例如：

```csv
id,sequence,substrate_smiles,primary_substrate,parameter
```

具体字段以 CatPred GitHub 的 README、demo 文件和示例输入为准。

#### Step 5：运行 CatPred 推理

本地运行示意：

```bash
git clone https://github.com/maranasgroup/CatPred.git
cd CatPred

# 按 README 创建环境
conda env create -f environment.yml
conda activate catpred

# 运行 kcat 推理
python demo_run.py \
  --parameter kcat \
  --input_file ../data/interim/prediction_inputs/catpred_kcat_input.csv
```

如果 CatPred 提供 Web app 批量上传，也可以先用 Web app 测试小批量数据，确认输入格式没问题，再在本地或服务器上跑全量数据。

#### Step 6：整理 CatPred 输出

写一个收集脚本：

```text
src/predictors/collect_catpred_output.py
```

把 CatPred 输出转成统一格式：

```csv
species,reaction_id,gene_id,ec_number,substrate_name,method,pred_kcat,pred_log10_kcat,status
ecoli,R_PFK,b3916,2.7.1.11,D-fructose 6-phosphate,CatPred,123.4,2.091,success
```

输出到：

```text
data/raw/method_outputs/catpred_predictions.csv
```

然后并入总表：

```text
data/final/benchmark_long.csv
```

### 5.5 CatPred 需要不要重新训练？

第一版小文章不建议一上来就重训 CatPred。

推荐做法是：

```text
先跑 CatPred 的预训练模型，得到我们的 E. coli + yeast 数据上的 kcat 预测值。
```

原因：

1. 我们的研究重点是“现有方法即插即用表现”；
2. 重训会显著增加工作量；
3. 重训后很难保证每个方法都公平重训；
4. 有些方法可能不给完整训练脚本；
5. 小文章不需要从一开始做到最高规格。

但需要注意：

如果 CatPred 的训练集包含了我们的 BRENDA 真值样本，那么评测可能存在 train-test overlap。解决方式不是第一时间重训，而是：

1. 下载 CatPred-DB 或其公开训练/测试数据；
2. 把我们的 benchmark truth 与 CatPred 训练数据做 exact deduplication；
3. 标注哪些样本可能在 CatPred 训练中出现；
4. 分别报告：
   - all benchmark；
   - non-overlap benchmark。

### 5.6 CatPred 的去泄漏处理

推荐构建一个去泄漏脚本：

```text
src/quality_control/check_catpred_overlap.py
```

匹配规则：

```text
Exact overlap:
same protein sequence + same substrate SMILES + same kcat value

Relaxed overlap:
same UniProt ID + same substrate SMILES

EC-level overlap:
same EC + same substrate
```

最终在 benchmark 表中增加字段：

```csv
catpred_exact_overlap
catpred_relaxed_overlap
catpred_ec_substrate_overlap
```

主结果可以报告：

```text
CatPred on all benchmark samples
CatPred on non-overlap samples
```

这样比简单说“可能有泄漏”更有说服力。

### 5.7 CatPred 接入完成后的最小可交付结果

CatPred 接入成功后，至少应该生成：

```text
data/interim/prediction_inputs/catpred_kcat_input.csv
data/raw/method_outputs/catpred_predictions.csv
data/final/benchmark_long.csv
reports/tables/catpred_coverage.csv
reports/tables/catpred_metrics.csv
reports/figures/catpred_pred_vs_true.png
reports/figures/catpred_bias.png
```

核心指标：

```text
coverage
RMSE
MAE
R²
Pearson r
Spearman rho
median bias
high-kcat underestimation
E. coli performance
yeast performance
common-set performance
non-overlap performance
```

---

## 6. 实验 kcat 真值获取路线

### 6.1 不建议继续直接用 ecModel kcat 作为真值

`reaction_kcat_MW_databasefill.csv` 可以保留，但应降级为 legacy / supplementary source。主 benchmark 真值应来自：

```text
BRENDA
SABIO-RK
```

如果为了快速推进，也可以使用 ecYeast 或 ecModel 的 kcat 作为初步结果，但正式文章中必须说清楚：

```text
These values were used only for preliminary comparison or sanity check, not as the primary experimental benchmark.
```

### 6.2 BRENDA/SABIO-RK 获取策略

输入：

```text
EC number
organism
substrate
UniProt ID
```

输出：

```text
kcat
unit
organism
substrate
pH
temperature
reference
comment
source database
```

推荐匹配等级：

| 等级 | 匹配条件 | 是否用于主 benchmark |
|---|---|---|
| Level 1 | EC + exact organism + substrate + wild-type | 是 |
| Level 2 | EC + exact organism，但 substrate 未完全匹配 | 是，但单独标注 |
| Level 3 | EC + substrate，但不是同物种 | 补充分析 |
| Level 4 | EC only | 不进入主 benchmark |
| Level 5 | ecModel/ecYeast/inferred/fill | 不作为主真值 |

多个实验值处理方式：

```text
对 kcat 做 log10 转换
对同一 enzyme-substrate entry 取 median(log10(kcat))
同时保留 n_measurements、min、max、pH、temperature、reference_count
```

---

## 7. 数据结构重构

当前项目主要是 reaction-level 宽表：

```csv
reaction,gpr,DLKcat,MTLKP,TurNup,UniKP,database,fill_method
```

最终建议改成 long-format：

```csv
species,model_id,reaction_id,reaction_direction,gene_id,gpr_group_id,enzyme_complex_type,uniprot_id,ec_number,substrate_name,substrate_id,substrate_smiles,reaction_smiles,protein_sequence,organism,temperature,pH,true_kcat,true_kcat_log10,true_source,true_match_level,true_reference,method,pred_kcat,pred_kcat_log10,prediction_status,prediction_error_message
```

这个 long-format 表是最终项目的核心。所有图和统计都从它生成。

---

## 8. 建议的新项目结构

```text
kcat_benchmark/
├── configs/
│   ├── species.yaml
│   ├── methods.yaml
│   └── matching_rules.yaml
├── data/
│   ├── raw/
│   │   ├── models/
│   │   │   ├── ecoli_iML1515.json
│   │   │   └── yeast_gem.xml
│   │   ├── brenda/
│   │   ├── sabiork/
│   │   └── method_outputs/
│   ├── interim/
│   │   ├── enzyme_reaction_entries.csv
│   │   ├── enzyme_reaction_entries_with_sequence.csv
│   │   ├── enzyme_reaction_entries_with_smiles.csv
│   │   ├── brenda_kcat_raw.csv
│   │   ├── brenda_kcat_matched.csv
│   │   └── prediction_inputs/
│   └── final/
│       ├── benchmark_long.csv
│       ├── benchmark_common_set.csv
│       └── benchmark_summary.csv
├── src/
│   ├── 01_parse_models.py
│   ├── 02_build_enzyme_reaction_entries.py
│   ├── 03_fetch_brenda_sabiork.py
│   ├── 04_match_experimental_kcat.py
│   ├── 05_prepare_method_inputs.py
│   ├── predictors/
│   │   ├── prepare_catpred_input.py
│   │   ├── run_catpred.py
│   │   ├── collect_catpred_output.py
│   │   ├── run_dlkcat.py
│   │   ├── run_unikp.py
│   │   ├── run_turnup.py
│   │   ├── run_mpek.py
│   │   ├── run_catapro.py
│   │   ├── run_pmak.py
│   │   ├── run_kinform.py
│   │   └── run_pretkcat.py
│   ├── 06_collect_predictions.py
│   ├── 07_check_train_test_overlap.py
│   ├── 08_benchmark_metrics.py
│   ├── 09_figures.py
│   └── utils/
├── reports/
│   ├── figures/
│   ├── tables/
│   └── manuscript_draft/
└── run_all.sh
```

---

## 9. 最终分析模块

### 9.1 主分析

必须做：

```text
1. Experimental kcat coverage in E. coli and yeast
2. Prediction coverage of each AI method
3. All-available benchmark
4. Common-set benchmark
5. Species-specific performance
6. Systematic bias and high-kcat underestimation
7. Method-method agreement
8. Isozyme-level variability
9. Substrate-level variability
10. Directional sensitivity for reversible reaction pairs
```

### 9.2 指标

所有误差建议在 `log10(kcat)` 空间计算：

```text
RMSE
MAE
R²
Pearson r
Spearman rho
median bias
interquartile error
coverage
failure rate
```

统计要求：

```text
bootstrap 95% CI
paired test on common samples
FDR correction for multiple testing
separate all-available and common-set comparisons
```

---

## 10. 推荐图表

### Figure 1：研究框架

```text
E. coli + S. cerevisiae GEMs
→ EC/GPR/substrate/sequence extraction
→ BRENDA/SABIO-RK kcat curation
→ AI predictors
→ cross-species benchmark
```

### Figure 2：实验 kcat 覆盖率

```text
total reactions
EC-annotated reactions
sequence-resolved enzyme entries
substrate-resolved entries
BRENDA/SABIO-RK matched entries
final benchmark entries
```

分别展示 E. coli 和 yeast。

### Figure 3：方法覆盖率与运行成功率

```text
method × species coverage
method × failure type
method × enzyme complex coverage
method × transport reaction coverage
```

### Figure 4：主 benchmark 性能

分两个部分：

```text
all-available benchmark
common-set benchmark
```

指标：

```text
RMSE
MAE
R²
Spearman rho
```

### Figure 5：系统性偏差

```text
predicted vs experimental log10(kcat)
Bland-Altman plot
bias across true kcat quantiles
high-kcat underestimation
```

### Figure 6：方法行为差异

```text
method-method correlation
isozyme variability
substrate specificity
directional sensitivity
```

---

## 11. 文章结果结构建议

```text
Result 1:
Construction of a cross-species enzyme-reaction benchmark from E. coli and S. cerevisiae metabolic models.

Result 2:
Experimental kcat coverage remains sparse and uneven across organisms, EC classes, and reaction types.

Result 3:
AI kcat predictors differ substantially in coverage, failure modes, and prediction distributions.

Result 4:
Prediction accuracy varies across species and is strongly affected by evaluation set definition.

Result 5:
Systematic underestimation of high-turnover enzymes is a common limitation across methods.

Result 6:
Reaction-aware and sequence-aware predictors capture complementary aspects of enzyme kinetics.
```

---

## 12. 分阶段执行计划

### Phase 1：重构数据体系，不跑新方法

目标：

```text
把 E. coli + yeast 的 benchmark truth 做干净。
```

任务：

1. 下载或整理 E. coli iML1515 / eciML1515；
2. 下载或整理 yeast-GEM / Yeast8 / Yeast9；
3. 提取 reaction、GPR、gene、EC、substrate；
4. 补 UniProt sequence；
5. 补 substrate SMILES；
6. 通过 EC + organism + substrate 从 BRENDA/SABIO-RK 获取 kcat；
7. 生成 curated experimental kcat truth。

输出：

```text
data/interim/enzyme_reaction_entries.csv
data/interim/brenda_kcat_matched.csv
data/final/experimental_kcat_truth.csv
reports/tables/experimental_kcat_coverage.csv
```

### Phase 2：接入已有 4 个方法

目标：

```text
用新真值体系复现当前 4 个方法的分析。
```

方法：

```text
DLKcat
UniKP
TurNuP
MPEK/MTLKP
```

输出：

```text
data/final/benchmark_long_v1.csv
reports/tables/metrics_v1.csv
reports/figures/v1_method_comparison.png
```

### Phase 3：接入新方法

优先级：

```text
1. CatPred
2. CataPro
3. PMAK
4. KinForm
5. PreTKcat
6. KcatNet
7. DEKP
8. SELFprot
```

| 优先级 | 方法           | 论文链接                                                                                                                                                                                        | GitHub / 代码链接                                                | 备注                                                                         |
| --: | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ | -------------------------------------------------------------------------- |
|   1 | **CatPred**  | Nature Communications 2025：*CatPred: a comprehensive framework for deep learning in vitro enzyme kinetic parameters* ([Nature][1])                                                          | GitHub：`maranasgroup/CatPred` ([GitHub][2])                  | 最优先接入。支持 kcat、Km、Ki；有 Web app、预训练模型、批量 CSV 输入。                             |
|   2 | **CataPro**  | Nature Communications 2025：*Robust enzyme discovery and engineering with deep learning using CataPro* ([Nature][3])                                                                         | GitHub：`zchwang/CataPro` ([GitHub][4])                       | 很适合放主文。预测 kcat、Km、kcat/Km，影响力高，代码公开。                                       |
|   3 | **PMAK**     | Communications Biology 2026：*Enhancing kcat prediction through residue-aware attention mechanism and pre-trained representations* ([Nature][5])                                             | GitHub：`MrVincentCai/PMAK` ([GitHub][6])                     | 重点引入 **reaction SMILES + residue-aware attention**，适合补足 reaction-aware 方法。 |
|   4 | **KinForm**  | npj Systems Biology and Applications 2026：*KinForm: kinetics-informed feature optimised representation models for enzyme kcat and KM prediction* ([Nature][7])                              | GitHub：`Digital-Metabolic-Twin-Centre/KinForm` ([GitHub][8]) | 适合作为严谨 benchmark 对照，强调低相似序列泛化和特征优化。                                        |
|   5 | **KcatNet**  | Genome Biology 2026：*A geometric deep learning framework for genome-wide prediction of enzyme turnover number* ([Springer][9])                                                              | GitHub：`BioColLab/KcatNet` ([GitHub][10])                    | 和 genome-wide / metabolic model 场景最相关，但工程复杂度会比 CatPred 高。                  |
|   6 | **PreTKcat** | Computational Biology and Chemistry 2025：*PreTKcat: A pre-trained representation learning and machine learning framework for predicting enzyme turnover number* ([ACM Digital Library][11]) | GitHub：`MrVincentCai/PreTKcat` ([GitHub][12])                | 轻量 baseline，适合放补充或作为简单可复现方法。                                               |
|   7 | **DEKP**     | Briefings in Bioinformatics 2025：*DEKP: a deep learning model for enzyme kinetic parameter prediction based on pretrained models and graph neural networks* ([OUP Academic][13])            | GitHub：`wang-yi-zhen/DEKP` ([GitHub][14])                    | 用 protein sequence、substrate、protein structure 多模态特征；可作为结构信息路线代表。          |
|   8 | **SELFprot** | JCIM 2025：*SELFprot: Effective and Efficient Multitask Finetuning Methods for Protein Parameter Prediction* ([美国化学学会出版物][15])                                                               | GitHub：`marltanwilson/SELFprot` ([GitHub][16])               | 多任务 protein–ligand 参数预测，包含 kcat，但不建议作为第一批主方法。                              |

[1]: https://www.nature.com/articles/s41467-025-57215-9?utm_source=chatgpt.com "CatPred: a comprehensive framework for deep learning in ..."
[2]: https://github.com/maranasgroup/CatPred?utm_source=chatgpt.com "maranasgroup/CatPred: Machine Learning models for in ..."
[3]: https://www.nature.com/articles/s41467-025-58038-4?utm_source=chatgpt.com "Robust enzyme discovery and engineering with deep ..."
[4]: https://github.com/zchwang/CataPro?utm_source=chatgpt.com "zchwang/CataPro: A generalized enzyme kinetics ..."
[5]: https://www.nature.com/articles/s42003-026-09551-9?utm_source=chatgpt.com "Enhancing kcat prediction through residue-aware attention ..."
[6]: https://github.com/MrVincentCai/PMAK?utm_source=chatgpt.com "MrVincentCai/PMAK"
[7]: https://www.nature.com/articles/s41540-026-00692-5?utm_source=chatgpt.com "KinForm: kinetics-informed feature optimised ..."
[8]: https://github.com/Digital-Metabolic-Twin-Centre/KinForm?utm_source=chatgpt.com "Digital-Metabolic-Twin-Centre/KinForm"
[9]: https://link.springer.com/article/10.1186/s13059-026-03986-3?utm_source=chatgpt.com "A geometric deep learning framework for genome-wide ..."
[10]: https://github.com/BioColLab/KcatNet?utm_source=chatgpt.com "BioColLab/KcatNet"
[11]: https://dl.acm.org/doi/abs/10.1016/j.compbiolchem.2024.108327?utm_source=chatgpt.com "PreTKcat: : A pre-trained representation learning and ..."
[12]: https://github.com/MrVincentCai/PreTKcat?utm_source=chatgpt.com "MrVincentCai/PreTKcat"
[13]: https://academic.oup.com/bib/article/26/2/bbaf187/8119324?utm_source=chatgpt.com "DEKP: a deep learning model for enzyme kinetic parameter ..."
[14]: https://github.com/wang-yi-zhen/DEKP?utm_source=chatgpt.com "wang-yi-zhen/DEKP"
[15]: https://pubs.acs.org/doi/10.1021/acs.jcim.4c02230?utm_source=chatgpt.com "SELFprot: Effective and Efficient Multitask Finetuning Methods ..."
[16]: https://github.com/marltanwilson/SELFprot?utm_source=chatgpt.com "marltanwilson/SELFprot: finetuned protein-ligand ..."


输出：

```text
data/final/benchmark_long_v2.csv
reports/tables/method_runtime_log.csv
reports/tables/prediction_failure_log.csv
```

### Phase 4：正式统计和写文章

输出：

```text
reports/figures/Figure1_framework.png
reports/figures/Figure2_coverage.png
reports/figures/Figure3_method_coverage.png
reports/figures/Figure4_accuracy.png
reports/figures/Figure5_bias.png
reports/figures/Figure6_method_behavior.png
reports/manuscript_draft/main_text.md
```

---

## 13. 最小可发表版本

如果时间有限，最小可发表版本建议做到：

```text
Species:
E. coli + S. cerevisiae

Methods:
DLKcat
UniKP
TurNuP
MPEK/MTLKP
CatPred
CataPro
PMAK

Truth:
BRENDA/SABIO-RK curated kcat
Level 1 + Level 2 match

Main analyses:
coverage
accuracy
common-set benchmark
species-specific performance
bias
method-method agreement
```

暂时可以不做：

```text
KcatNet
KinForm
DEKP
SELFprot
full retraining
wet-lab validation
docking-based methods
```

这样文章工作量可控，也比当前 E. coli 单物种报告有明显提升。

---

## 14. 一句话总结

当前项目已经有一个不错的 E. coli 原型，但要发文章，需要从：

```text
ecModel kcat 填充值比较报告
```

升级为：

```text
E. coli + S. cerevisiae 跨物种 AI kcat 预测方法系统 benchmark
```

文献和 GitHub 拿到后，第一步不是完整复现原文所有训练实验，而是：

```text
把每个方法作为一个 off-the-shelf predictor，
对我们自己构建的 enzyme-reaction benchmark 跑统一推理，
得到 predicted kcat，
再与 BRENDA/SABIO-RK curated experimental kcat 比较。
```

以 CatPred 为例，具体就是：

```text
我们的 enzyme-reaction entries
→ 转成 CatPred 需要的 sequence + substrate SMILES 输入表
→ 调用 CatPred kcat 预训练模型或 Web/脚本
→ 得到 pred_kcat
→ 转成统一 long-format 输出
→ 合并到 benchmark_long.csv
→ 计算 coverage、RMSE、MAE、R²、Spearman、bias、cross-species performance
```

这条路线最现实、最快，也最符合“小文章”的目标。
