#!/usr/bin/env python3
"""Generate dataset context and method comparison material for the kcat benchmark."""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime
from pathlib import Path
import re
import textwrap
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries.csv"
REACTIONS = BASE / "data" / "interim" / "model_reactions.csv"
YEAST_MODEL = BASE / "yeast-GEM.xml"
MODULE_EC = BASE / "external_methods" / "DLKcat_official" / "DeeplearningApproach" / "Data" / "subsystem" / "module_ec.txt"

REPORT_DIR = BASE / "reports"
TABLE_DIR = REPORT_DIR / "tables"
FIG_DIR = REPORT_DIR / "figures" / "kcat_dataset_context"
REPORT_PATH = REPORT_DIR / "kcat_benchmark_dataset_and_method_context.md"

BENCHMARK_N = len(pd.read_csv(DATA, usecols=["entry_id"]))
SBML = "http://www.sbml.org/sbml/level3/version1/core"
RDF_RESOURCE = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"

SPECIES_LABELS = {
    "ecoli": "Escherichia coli / 大肠杆菌",
    "yeast": "Saccharomyces cerevisiae / 酿酒酵母",
}

EC_CLASS_NAMES = {
    "1": "1 Oxidoreductases / 氧化还原酶",
    "2": "2 Transferases / 转移酶",
    "3": "3 Hydrolases / 水解酶",
    "4": "4 Lyases / 裂合酶",
    "5": "5 Isomerases / 异构酶",
    "6": "6 Ligases / 连接酶",
    "7": "7 Translocases / 转位酶",
}

GROUP_SHORT_NAMES = {
    "Primary - Carbohydrate & Energy Metabolism": "Primary carbohydrate/energy",
    "Primary - amino acids, fatty acids and nucleotides": "Primary AA/FA/nt",
    "Secondary": "Secondary",
    "Secondary_other": "Secondary other",
    "Intermediate": "Intermediate",
    "x": "Unclassified module",
    "Unmapped/No exact EC-module match": "Unmapped",
}

GROUP_PRIORITY = [
    "Primary - Carbohydrate & Energy Metabolism",
    "Primary - amino acids, fatty acids and nucleotides",
    "Intermediate",
    "Secondary",
    "Secondary_other",
    "x",
]

KEGG_PATHWAY_NAMES = {
    "00010": "Glycolysis / Gluconeogenesis",
    "00020": "Citrate cycle (TCA cycle)",
    "00030": "Pentose phosphate pathway",
    "00040": "Pentose and glucuronate interconversions",
    "00051": "Fructose and mannose metabolism",
    "00052": "Galactose metabolism",
    "00061": "Fatty acid biosynthesis",
    "00062": "Fatty acid elongation",
    "00071": "Fatty acid degradation",
    "00072": "Synthesis and degradation of ketone bodies",
    "00100": "Steroid biosynthesis",
    "00130": "Ubiquinone and other terpenoid-quinone biosynthesis",
    "00220": "Arginine biosynthesis",
    "00230": "Purine metabolism",
    "00240": "Pyrimidine metabolism",
    "00250": "Alanine, aspartate and glutamate metabolism",
    "00260": "Glycine, serine and threonine metabolism",
    "00270": "Cysteine and methionine metabolism",
    "00280": "Valine, leucine and isoleucine degradation",
    "00290": "Valine, leucine and isoleucine biosynthesis",
    "00300": "Lysine biosynthesis",
    "00310": "Lysine degradation",
    "00330": "Arginine and proline metabolism",
    "00340": "Histidine metabolism",
    "00350": "Tyrosine metabolism",
    "00360": "Phenylalanine metabolism",
    "00380": "Tryptophan metabolism",
    "00400": "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "00430": "Taurine and hypotaurine metabolism",
    "00450": "Selenocompound metabolism",
    "00480": "Glutathione metabolism",
    "00500": "Starch and sucrose metabolism",
    "00520": "Amino sugar and nucleotide sugar metabolism",
    "00561": "Glycerolipid metabolism",
    "00562": "Inositol phosphate metabolism",
    "00620": "Pyruvate metabolism",
    "00630": "Glyoxylate and dicarboxylate metabolism",
    "00640": "Propanoate metabolism",
    "00650": "Butanoate metabolism",
    "00660": "C5-branched dibasic acid metabolism",
    "00730": "Thiamine metabolism",
    "00740": "Riboflavin metabolism",
    "00770": "Pantothenate and CoA biosynthesis",
    "00780": "Biotin metabolism",
    "00790": "Folate biosynthesis",
    "00900": "Terpenoid backbone biosynthesis",
    "00920": "Sulfur metabolism",
    "01040": "Biosynthesis of unsaturated fatty acids",
    "01100": "Metabolic pathways",
    "01110": "Biosynthesis of secondary metabolites",
    "01120": "Microbial metabolism in diverse environments",
    "01130": "Biosynthesis of antibiotics",
    "01200": "Carbon metabolism",
    "01210": "2-Oxocarboxylic acid metabolism",
    "01212": "Fatty acid metabolism",
    "01230": "Biosynthesis of amino acids",
}

CURRENCY_OR_COFACTOR_NAMES = {
    "h+",
    "h",
    "proton",
    "h2o",
    "water",
    "atp",
    "adp",
    "amp",
    "gtp",
    "gdp",
    "gmp",
    "ctp",
    "cdp",
    "cmp",
    "utp",
    "udp",
    "ump",
    "imp",
    "nad",
    "nadh",
    "nadp",
    "nadph",
    "fad",
    "fadh2",
    "fmn",
    "coa",
    "acetyl-coa",
    "accoa",
    "co2",
    "o2",
    "oxygen",
    "ammonia",
    "nh3",
    "nh4",
    "phosphate",
    "diphosphate",
    "pyrophosphate",
    "polyphosphate",
    "prpp",
}

METHOD_TECHNICAL_ROWS = [
    {
        "method": "DLKcat-official",
        "technical_principle_cn": "把底物 SMILES 转成分子图，把酶序列转成序列特征；用 GNN/CNN 类深度网络回归 kcat。",
        "plain_language_cn": "看一个底物长什么样、酶序列长什么样，然后学习二者组合对应的转换速度。",
        "input_needed": "enzyme sequence + substrate name/SMILES",
        "representation": "molecular graph fingerprints + protein attention-CNN",
        "model_family": "deep learning, GNN + CNN",
        "benchmark_dimension": "全量/近全量 sequence+SMILES 基线",
        "main_caveat": "较早一代模型；不显式使用完整反应物/产物信息。",
    },
    {
        "method": "UniKP-official",
        "technical_principle_cn": "用 ProtT5 表示蛋白、SMILES Transformer 表示底物，再用传统机器学习模型预测 kcat/Km/kcat/Km。",
        "plain_language_cn": "先让两个预训练模型分别读懂蛋白和小分子，再把特征交给回归器。",
        "input_needed": "enzyme sequence + substrate SMILES",
        "representation": "ProtT5 protein embedding + SMILES Transformer substrate embedding",
        "model_family": "pretrained embeddings + machine learning",
        "benchmark_dimension": "全量/近全量 sequence+SMILES",
        "main_caveat": "依赖预计算特征和模型版本；不显式看产物侧。",
    },
    {
        "method": "TurNuP-official",
        "technical_principle_cn": "用反应物和产物生成 reaction difference fingerprint，再拼接 ESM1b 蛋白向量，用 XGBoost 预测 kcat。",
        "plain_language_cn": "不仅看底物，还看反应前后分子怎么变，再结合酶序列预测速度。",
        "input_needed": "reactant SMILES + product SMILES + enzyme sequence",
        "representation": "reaction difference fingerprint + ESM1b enzyme embedding",
        "model_family": "reaction-aware ML, XGBoost",
        "benchmark_dimension": "reaction-aware 子集",
        "main_caveat": "必须有完整反应 SMILES；当前覆盖 1047/1246。",
    },
    {
        "method": "CatPred",
        "technical_principle_cn": "CatPred 官方推理流程以 SMILES、酶序列和 protein record/pdbpath 为输入，使用集成模型输出 kcat/Km/Ki 等体外动力学参数。",
        "plain_language_cn": "这是一个面向多种酶动力学参数的统一深度学习框架，本项目只取其中 kcat 模型来评测。",
        "input_needed": "SMILES + sequence + unique pdbpath/protein record",
        "representation": "molecular/protein records used by CatPred checkpoints",
        "model_family": "deep learning ensemble",
        "benchmark_dimension": "模型特定子集",
        "main_caveat": "官方流程有额外可处理范围限制；当前覆盖 1156/1246。",
    },
    {
        "method": "CataPro",
        "technical_principle_cn": "结合蛋白语言模型、小分子语言模型和分子指纹，预测 kcat/Km/kcat/Km。",
        "plain_language_cn": "同时让模型读懂酶序列、小分子字符串和传统分子指纹。",
        "input_needed": "enzyme sequence + substrate SMILES + wild/mutant type",
        "representation": "ProtT5 + MolT5/SMILES features + molecular fingerprints",
        "model_family": "pretrained embeddings + neural regression",
        "benchmark_dimension": "全量/近全量 sequence+SMILES",
        "main_caveat": "需要外部预训练权重；不显式建模产物侧。",
    },
    {
        "method": "PMAK",
        "technical_principle_cn": "引入完整反应信息和 residue-aware attention，融合蛋白预训练表示与 RXNFP 反应表示。",
        "plain_language_cn": "让模型同时看酶、反应，并尝试关注对催化更关键的残基位置。",
        "input_needed": "reaction SMILES + enzyme sequence",
        "representation": "ProtT5 enzyme embedding + RXNFP/reaction embedding",
        "model_family": "reaction-aware deep learning",
        "benchmark_dimension": "reaction-aware 子集",
        "main_caveat": "依赖完整反应 SMILES；当前覆盖 1047/1246。",
    },
    {
        "method": "KinForm",
        "technical_principle_cn": "使用多种蛋白 embedding、可选结合位点信息和底物 SMILES，训练/预测 kcat 与 KM。",
        "plain_language_cn": "尽量把整条蛋白和潜在结合位点的信息都编码进去，再和底物信息融合。",
        "input_needed": "sequence + SMILES, plus cached embeddings/assets",
        "representation": "ESM/ESMC/ProtT5 embeddings + optional binding-site features",
        "model_family": "embedding-based deep/ML models",
        "benchmark_dimension": "模型特定子集",
        "main_caveat": "受官方 Zenodo bundle/缓存资产覆盖限制；当前覆盖 729/1246。",
    },
    {
        "method": "KcatNet",
        "technical_principle_cn": "几何深度学习框架，融合 ProtT5/ESM 蛋白表示和底物 SMILES 表示来预测 turnover number。",
        "plain_language_cn": "把蛋白和底物都编码成结构化特征，再用深度网络学习它们的匹配关系。",
        "input_needed": "enzyme sequence + substrate SMILES",
        "representation": "ProtT5/ESM protein embeddings + SMILES Transformer/substrate features",
        "model_family": "geometric/deep learning",
        "benchmark_dimension": "全量/近全量 sequence+SMILES",
        "main_caveat": "序列会按模型规则截断；不显式使用产物侧。",
    },
    {
        "method": "PreTKcat",
        "technical_principle_cn": "ProtT5 编码酶序列，MolGNet 编码底物分子图，并加入温度特征，最后用 ExtraTrees 预测 kcat。",
        "plain_language_cn": "除了酶和底物，还把实验温度作为影响 kcat 的因素放进去。",
        "input_needed": "sequence + SMILES + temperature",
        "representation": "ProtT5 protein embedding + MolGNet molecular graph embedding + temperature features",
        "model_family": "pretrained embeddings + ExtraTrees",
        "benchmark_dimension": "全量/近全量 sequence+SMILES",
        "main_caveat": "缺失温度需要填补默认值；主结果使用 exact-overlap-excluded 公开重建，另报告 raw-public 和 near-excluded 敏感性结果。",
    },
    {
        "method": "DEKP-public-retrained",
        "technical_principle_cn": "DEKP 使用 ProtT5、SMILES Transformer、PST/MolFormer 及结构/图特征；本项目评测的是公开数据重训版。",
        "plain_language_cn": "它比普通 sequence+SMILES 方法多看蛋白结构相关信息，但我们没有原论文最优私有权重，因此重新训练了公开可复现版本。",
        "input_needed": "sequence + SMILES + protein structure/graph assets",
        "representation": "protein/substrate pretrained features + structure graph features",
        "model_family": "structure-aware deep learning",
        "benchmark_dimension": "公开数据重训版",
        "main_caveat": "不是官方最优权重；结果反映当前公开复现流程，不等同于论文声称上限。",
    },
    {
        "method": "SELFprot",
        "technical_principle_cn": "基于 protein-ligand interaction 预训练/微调框架，使用 CatPred-DB split 相关训练划分微调 kcat。",
        "plain_language_cn": "把酶和小分子看成一对相互作用对象，学习它们共同决定的动力学数值。",
        "input_needed": "sequence + SMILES",
        "representation": "chemical encoder + protein encoder + joint interaction layer",
        "model_family": "protein-ligand deep model",
        "benchmark_dimension": "全量/近全量 sequence+SMILES",
        "main_caveat": "README 信息较简略；需要在论文写作时继续补齐正式引用与模型细节。",
    },
    {
        "method": "GO-HKP",
        "technical_principle_cn": "按 GO 功能层级在 GO-kcat 数据库中取统计 kcat 赋值；E. coli 使用本地 DeepGO-SE 反应赋值，yeast 使用 UniProt GO 注释补齐。",
        "plain_language_cn": "它不是训练一个 AI 回归模型，而是找功能相似的酶/GO 节点，把已有 kcat 作为参考值直接赋给目标反应或基因。",
        "input_needed": "protein/gene GO terms or DeepGO-SE GO predictions + reaction/gene mapping",
        "representation": "GO hierarchy + GO-term kcat statistics",
        "model_family": "functional-similarity assignment baseline",
        "benchmark_dimension": "功能相似性 GO 赋值基线",
        "main_caveat": "当前覆盖 1236/1246；E. coli 与 yeast 的 GO 来源不同，yeast 是 UniProt GO 注释路线，不是 DeepGO-SE 预测路线。",
    },
]

METHOD_SCOPE_ORDER = [
    "全量/近全量 sequence+SMILES",
    "reaction-aware 子集",
    "模型特定子集",
    "公开数据重训版",
    "功能相似性 GO 赋值基线",
]


def ensure_dirs() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def split_ecs(value: object) -> list[str]:
    if pd.isna(value):
        return []
    out: list[str] = []
    for item in re.split(r"[;,|]", str(value)):
        ec = item.strip()
        if not ec or ec.lower() == "nan":
            continue
        ec = re.sub(r"^(ec:|EC)", "", ec, flags=re.IGNORECASE).strip()
        if re.match(r"^\d+\.\d+\.\d+\.[\d-]+$", ec) or re.match(r"^\d+\.\d+\.\d+\.-$", ec):
            out.append(ec)
    return out


def exact_ecs(value: object) -> list[str]:
    return [ec for ec in split_ecs(value) if "-" not in ec]


def ec_class(value: object) -> list[str]:
    classes = []
    for ec in split_ecs(value):
        first = ec.split(".", 1)[0]
        classes.append(EC_CLASS_NAMES.get(first, f"{first} Other EC class"))
    return sorted(set(classes))


def normalize_substrate_name(value: object) -> str:
    text = str(value).strip().lower()
    text = text.replace("−", "-")
    text = re.sub(r"\s+", " ", text)
    return text


def first_name_token(value: object) -> str:
    text = normalize_substrate_name(value)
    if not text or text == "nan":
        return ""
    text = text.replace("(", "").replace(")", "")
    text = text.replace("[", "").replace("]", "")
    return re.split(r"\s+", text)[0]


def is_currency_or_cofactor(value: object) -> bool:
    text = normalize_substrate_name(value)
    token = first_name_token(text)
    if text in CURRENCY_OR_COFACTOR_NAMES or token in CURRENCY_OR_COFACTOR_NAMES:
        return True
    if token.endswith("-coa") or token.endswith("coa"):
        return True
    return False


def load_module_ec_map() -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    ec_to_groups: dict[str, set[str]] = defaultdict(set)
    ec_to_modules: dict[str, set[str]] = defaultdict(set)
    if not MODULE_EC.exists():
        return ec_to_groups, ec_to_modules
    with MODULE_EC.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            module_id = parts[0].replace("md:", "")
            ec = parts[1].replace("EC", "").strip()
            group = parts[2].strip()
            if not ec or not group:
                continue
            ec_to_groups[ec].add(group)
            ec_to_modules[ec].add(module_id)
    return ec_to_groups, ec_to_modules


def choose_primary_group(groups: list[str]) -> str:
    if not groups:
        return "Unmapped/No exact EC-module match"
    for group in GROUP_PRIORITY:
        if group in groups:
            return group
    return sorted(groups)[0]


def collect_resources(element: ET.Element | None) -> list[str]:
    if element is None:
        return []
    return [node.attrib[RDF_RESOURCE] for node in element.iter() if RDF_RESOURCE in node.attrib]


def parse_yeast_kegg_pathways() -> dict[str, list[str]]:
    if not YEAST_MODEL.exists():
        return {}
    root = ET.parse(YEAST_MODEL).getroot()
    ns = {"sbml": SBML}
    mapping: dict[str, list[str]] = {}
    for rxn in root.findall(".//sbml:listOfReactions/sbml:reaction", ns):
        reaction_id = rxn.attrib.get("id", "")
        annotation = rxn.find("sbml:annotation", ns)
        pathways = []
        for uri in collect_resources(annotation):
            if "/kegg.pathway/" in uri:
                pathways.append(uri.rsplit("/", 1)[-1])
        if pathways:
            mapping[reaction_id] = sorted(set(pathways))
    return mapping


def pathway_number(pathway_id: str) -> str:
    match = re.search(r"(\d{5})$", str(pathway_id))
    return match.group(1) if match else str(pathway_id)


def pathway_name(pathway_id: str) -> str:
    return KEGG_PATHWAY_NAMES.get(pathway_number(pathway_id), "")


def load_enriched_benchmark() -> pd.DataFrame:
    bench = pd.read_csv(DATA)
    entries = pd.read_csv(ENTRIES)
    reactions = pd.read_csv(REACTIONS)

    entry_cols = [
        "entry_id",
        "reaction_name",
        "substrate_id",
        "substrate_selection",
        "substrate_is_cofactor_like",
        "substrate_kegg_id",
        "substrate_metanetx_id",
    ]
    reaction_cols = [
        "species",
        "reaction_id",
        "reaction_name",
        "bigg_reaction",
        "kegg_reaction",
        "rhea",
        "metanetx_reaction",
        "reactant_ids",
        "product_ids",
    ]
    entry_cols = [c for c in entry_cols if c == "entry_id" or c not in bench.columns]
    out = bench.merge(entries[[c for c in entry_cols if c in entries.columns]], on="entry_id", how="left")
    out = out.merge(
        reactions[[c for c in reaction_cols if c in reactions.columns]],
        on=["species", "reaction_id"],
        how="left",
        suffixes=("_entry", "_model"),
    )
    out["reaction_name"] = out.get("reaction_name_entry", "").fillna("")
    if "reaction_name_model" in out.columns:
        out["reaction_name"] = out["reaction_name"].where(out["reaction_name"].astype(str).str.len() > 0, out["reaction_name_model"])
    out["reaction_name"] = out["reaction_name"].fillna(out["reaction_id"])

    out["sequence_length"] = out["sequence"].fillna("").astype(str).str.len()
    out["smiles_length"] = out["SMILES"].fillna("").astype(str).str.len()
    out["species_label"] = out["species"].map(SPECIES_LABELS).fillna(out["species"])
    out["currency_or_cofactor_like_by_name"] = out["substrate_name"].map(is_currency_or_cofactor)
    out["currency_or_cofactor_like_registry"] = out["substrate_role_group"].eq("currency_or_cofactor")
    out["ec_terms"] = out["ec_number"].map(lambda x: ";".join(split_ecs(x)))
    out["ec_exact_terms"] = out["ec_number"].map(lambda x: ";".join(exact_ecs(x)))
    out["ec_classes"] = out["ec_number"].map(lambda x: ";".join(ec_class(x)))

    ec_to_groups, ec_to_modules = load_module_ec_map()

    def groups_for_row(ec_value: object) -> list[str]:
        groups: set[str] = set()
        for ec in exact_ecs(ec_value):
            groups.update(ec_to_groups.get(ec, set()))
        return sorted(groups)

    def modules_for_row(ec_value: object) -> list[str]:
        modules: set[str] = set()
        for ec in exact_ecs(ec_value):
            modules.update(ec_to_modules.get(ec, set()))
        return sorted(modules)

    out["kegg_like_module_groups"] = out["ec_number"].map(lambda x: ";".join(groups_for_row(x)))
    out["kegg_like_module_ids"] = out["ec_number"].map(lambda x: ";".join(modules_for_row(x)))
    out["kegg_like_primary_group"] = out["ec_number"].map(lambda x: choose_primary_group(groups_for_row(x)))
    out["kegg_like_primary_group_short"] = out["kegg_like_primary_group"].map(GROUP_SHORT_NAMES).fillna(out["kegg_like_primary_group"])

    yeast_pathways = parse_yeast_kegg_pathways()
    out["direct_kegg_pathways"] = out.apply(
        lambda row: ";".join(yeast_pathways.get(str(row["reaction_id"]), [])) if row["species"] == "yeast" else "",
        axis=1,
    )
    return out


def percent(n: float, denom: float | None = None) -> float:
    if denom is None:
        denom = BENCHMARK_N
    return 100.0 * n / denom if denom else np.nan


def write_csv(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


def make_benchmark_build_funnel(df: pd.DataFrame) -> pd.DataFrame:
    parse = pd.read_csv(TABLE_DIR / "model_parse_summary.csv").set_index("species")
    readiness = pd.read_csv(TABLE_DIR / "final_benchmark_readiness.csv").set_index("species")
    rows = []
    for species in sorted(df["species"].unique()):
        part = df[df["species"].eq(species)]
        rows.append(
            {
                "species": species,
                "model_total_reactions": int(parse.at[species, "total_reactions"]),
                "reactions_with_gpr": int(parse.at[species, "reactions_with_gpr"]),
                "reactions_with_ec": int(parse.at[species, "reactions_with_ec"]),
                "enzyme_substrate_entries": int(readiness.at[species, "entries"]),
                "entries_with_uniprot_sequence": int(readiness.at[species, "entries_with_sequence"]),
                "entries_with_substrate_smiles": int(readiness.at[species, "entries_with_smiles"]),
                "experimental_truth_rows": int(readiness.at[species, "truth_rows"]),
                "benchmark_ready_rows": int(len(part)),
                "benchmark_unique_reactions": int(part["reaction_id"].nunique()),
                "benchmark_unique_genes": int(part["gene_id"].nunique()),
            }
        )
    return write_csv(pd.DataFrame(rows), TABLE_DIR / "benchmark_build_funnel.csv")


def make_project_directory_map() -> pd.DataFrame:
    rows = [
        {
            "category": "Benchmark construction",
            "path": "src/01_parse_models.py through src/11_finalize_benchmark_data.py; configs/",
            "directory_type": "code and rules",
            "contents": "GEM parsing, GPR/EC/substrate extraction, UniProt sequence retrieval, SMILES mapping, BRENDA/SABIO-RK truth matching, and final benchmark filtering.",
        },
        {
            "category": "Raw source data",
            "path": "data/raw/",
            "directory_type": "large source/cache data",
            "contents": "BRENDA, SABIO-RK, compound/CKB, UniProt FASTA, GO mappings, and query caches. Source databases are downloaded separately; this directory is not published in Git.",
        },
        {
            "category": "Intermediate curation",
            "path": "data/interim/",
            "directory_type": "rebuildable intermediate data",
            "contents": "Reaction-entry tables, sequence/SMILES queues, caches, review lists, reaction SMILES, and method input preparation tables.",
        },
        {
            "category": "Unified benchmark and method outputs",
            "path": "data/final/",
            "directory_type": "final data products",
            "contents": "Git tracks the three canonical benchmark CSVs. Per-method inputs, predictions, structures, and evaluated rows are local artifacts restored or rebuilt as needed.",
        },
        {
            "category": "Pipeline runners",
            "path": "scripts/runners/",
            "directory_type": "local and Slurm entry points",
            "contents": "Portable shell entry points for benchmark construction, method input preparation, prediction, and full cluster jobs. Run from the repository root after activating the required environment.",
        },
        {
            "category": "GO analysis",
            "path": "external_methods/GO-HKP/; data/raw/go_hkp/; data/final/go_hkp/; reports/tables/go_hkp_*",
            "directory_type": "functional-assignment analysis",
            "contents": "GO hierarchy and GO-kcat resources, E. coli DeepGO-SE assignments, yeast UniProt GO mappings, GO-HKP evaluated rows, readiness, and species-level metrics.",
        },
        {
            "category": "KEGG/EC/pathway analysis",
            "path": "src/47_generate_dataset_method_context_report.py; reports/tables/benchmark_dataset_kegg*; reports/tables/benchmark_dataset_direct_yeast_kegg_pathways.csv",
            "directory_type": "functional distribution analysis",
            "contents": "EC-to-module KEGG-like groups across species and direct yeast-GEM KEGG pathway annotations.",
        },
        {
            "category": "MAE and error analyses",
            "path": "reports/tables/*_eval_metrics.csv; reports/figures/kcat_benchmark_summary/",
            "directory_type": "performance analysis",
            "contents": "Overall and grouped MAE/RMSE, correlation, bias, within-fold error, coverage-error tradeoff, error distributions, and predicted-versus-true plots.",
        },
        {
            "category": "Species-level analysis",
            "path": "reports/tables/species_mae_matrix.csv; reports/tables/benchmark_dataset_*_by_species.csv",
            "directory_type": "species-stratified analysis",
            "contents": "E. coli versus yeast counts, truth distributions, sources, matching levels, pathway groups, and method MAE.",
        },
        {
            "category": "Method-level analysis",
            "path": "reports/tables/method_*.csv; data/final/<method>/",
            "directory_type": "method comparison",
            "contents": "Method principles, input requirements, comparison groups, coverage, rankings, evaluated predictions, and method-specific limitations.",
        },
        {
            "category": "Reports and publication tables",
            "path": "reports/; reports/report_tables/; docs/",
            "directory_type": "human-readable and manuscript material",
            "contents": "Main analysis report, dataset/method context report, figures, standalone report tables, and the chronological work log.",
        },
        {
            "category": "Third-party methods and model assets",
            "path": "external_methods/",
            "directory_type": "third-party code and large assets",
            "contents": "Only METHOD_SOURCES.md is tracked in Git. Third-party source trees come from upstream repositories; selected task weights and data bundles are restored from Zenodo.",
        },
        {
            "category": "Release staging and logs",
            "path": "release/; logs/",
            "directory_type": "generated local artifacts",
            "contents": "Zenodo upload staging, publication state, scheduler stdout/stderr, and run logs. Both directories are excluded from Git.",
        },
    ]
    publication_status = {
        "Benchmark construction": "GitHub",
        "Raw source data": "Local only; gitignored and downloaded from source databases",
        "Intermediate curation": "Local only; gitignored and rebuildable",
        "Unified benchmark and method outputs": "Mixed; three core CSVs on GitHub, method artifacts local/Zenodo",
        "Pipeline runners": "GitHub",
        "GO analysis": "Mixed; summary tables on GitHub, source assets local",
        "KEGG/EC/pathway analysis": "GitHub outputs with local source assets",
        "MAE and error analyses": "GitHub summary tables and figures",
        "Species-level analysis": "GitHub summary tables and figures",
        "Method-level analysis": "Mixed; summaries on GitHub, row-level artifacts local/Zenodo",
        "Reports and publication tables": "GitHub",
        "Third-party methods and model assets": "Mixed; manifest on GitHub, assets from upstream/Zenodo",
        "Release staging and logs": "Local only; gitignored",
    }
    for row in rows:
        row["publication_status"] = publication_status[row["category"]]

    return write_csv(pd.DataFrame(rows), TABLE_DIR / "project_directory_analysis_map.csv")


def make_species_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species, part in df.groupby("species", sort=True):
        rows.append(
            {
                "species": species,
                "species_label": SPECIES_LABELS.get(species, species),
                "rows": len(part),
                "percent": percent(len(part)),
                "unique_reactions": part["reaction_id"].nunique(),
                "unique_reaction_names": part["reaction_name"].nunique(),
                "unique_genes": part["gene_id"].nunique(),
                "unique_uniprots": part["uniprot_id"].nunique(),
                "unique_substrates": part["substrate_name"].nunique(),
                "unique_smiles": part["SMILES"].nunique(),
                "unique_ec_strings": part["ec_number"].nunique(),
                "median_kcat_s-1": part["true_kcat"].median(),
                "median_log10_kcat": part["true_kcat_log10"].median(),
                "pH_available_rows": int(part["pH"].notna().sum()),
                "temperature_available_rows": int(part["temperature_c"].notna().sum()),
                "currency_or_cofactor_like_rows_by_name": int(part["currency_or_cofactor_like_by_name"].sum()),
                "currency_or_cofactor_rows_registry": int(part["currency_or_cofactor_like_registry"].sum()),
                "carrier_linked_variable_rows_registry": int(part["substrate_role_group"].eq("carrier_linked_variable").sum()),
            }
        )
    return write_csv(pd.DataFrame(rows), TABLE_DIR / "benchmark_dataset_species_summary.csv")


def make_kcat_stats(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [("all", df)] + list(df.groupby("species", sort=True))
    for name, part in groups:
        rows.append(
            {
                "group": name,
                "n": len(part),
                "true_kcat_min": part["true_kcat"].min(),
                "true_kcat_q25": part["true_kcat"].quantile(0.25),
                "true_kcat_median": part["true_kcat"].median(),
                "true_kcat_mean": part["true_kcat"].mean(),
                "true_kcat_q75": part["true_kcat"].quantile(0.75),
                "true_kcat_max": part["true_kcat"].max(),
                "log10_kcat_min": part["true_kcat_log10"].min(),
                "log10_kcat_q25": part["true_kcat_log10"].quantile(0.25),
                "log10_kcat_median": part["true_kcat_log10"].median(),
                "log10_kcat_mean": part["true_kcat_log10"].mean(),
                "log10_kcat_q75": part["true_kcat_log10"].quantile(0.75),
                "log10_kcat_max": part["true_kcat_log10"].max(),
            }
        )
    return write_csv(pd.DataFrame(rows), TABLE_DIR / "benchmark_dataset_kcat_stats_by_species.csv")


def value_count_table(df: pd.DataFrame, columns: list[str], name: str) -> pd.DataFrame:
    table = df.groupby(columns, dropna=False).size().reset_index(name="rows")
    table["percent_of_benchmark"] = table["rows"].map(percent)
    return write_csv(table.sort_values("rows", ascending=False), TABLE_DIR / f"benchmark_dataset_{name}.csv")


def make_ec_class_summary(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df[["entry_id", "species", "ec_number"]].copy()
    tmp["ec_class"] = tmp["ec_number"].map(ec_class)
    tmp = tmp.explode("ec_class")
    tmp = tmp[tmp["ec_class"].notna() & tmp["ec_class"].astype(str).ne("")]
    table = tmp.groupby(["species", "ec_class"], dropna=False).size().reset_index(name="row_memberships")
    table["percent_of_benchmark"] = table["row_memberships"].map(percent)
    return write_csv(table.sort_values(["species", "row_memberships"], ascending=[True, False]), TABLE_DIR / "benchmark_dataset_ec_class_summary.csv")


def make_top_reactions(df: pd.DataFrame) -> pd.DataFrame:
    table = (
        df.groupby(["species", "reaction_id", "reaction_name"], dropna=False)
        .agg(
            rows=("entry_id", "count"),
            unique_substrates=("substrate_name", "nunique"),
            unique_genes=("gene_id", "nunique"),
            median_kcat=("true_kcat", "median"),
            median_log10_kcat=("true_kcat_log10", "median"),
        )
        .reset_index()
    )
    table["percent_of_benchmark"] = table["rows"].map(percent)
    return write_csv(table.sort_values("rows", ascending=False), TABLE_DIR / "benchmark_dataset_top_reactions.csv")


def make_top_substrates(df: pd.DataFrame) -> pd.DataFrame:
    table = (
        df.groupby(["substrate_name"], dropna=False)
        .agg(
            rows=("entry_id", "count"),
            species_count=("species", "nunique"),
            unique_reactions=("reaction_id", "nunique"),
            median_kcat=("true_kcat", "median"),
            currency_or_cofactor_like_by_name=("currency_or_cofactor_like_by_name", "max"),
            substrate_role_group=("substrate_role_group", lambda values: ";".join(sorted(set(values)))),
        )
        .reset_index()
    )
    table["percent_of_benchmark"] = table["rows"].map(percent)
    return write_csv(table.sort_values("rows", ascending=False), TABLE_DIR / "benchmark_dataset_top_substrates.csv")


def make_pathway_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary = (
        df.groupby(["species", "kegg_like_primary_group", "kegg_like_primary_group_short"], dropna=False)
        .size()
        .reset_index(name="rows")
    )
    primary["percent_of_species"] = primary["rows"] / primary.groupby("species")["rows"].transform("sum") * 100
    primary["percent_of_benchmark"] = primary["rows"].map(percent)
    primary = write_csv(
        primary.sort_values(["species", "rows"], ascending=[True, False]),
        TABLE_DIR / "benchmark_dataset_kegg_like_primary_group.csv",
    )

    membership = df[["entry_id", "species", "kegg_like_module_groups"]].copy()
    membership["kegg_like_module_group"] = membership["kegg_like_module_groups"].str.split(";")
    membership = membership.explode("kegg_like_module_group")
    membership["kegg_like_module_group"] = membership["kegg_like_module_group"].fillna("")
    membership.loc[membership["kegg_like_module_group"].eq(""), "kegg_like_module_group"] = "Unmapped/No exact EC-module match"
    membership = (
        membership.groupby(["species", "kegg_like_module_group"], dropna=False)
        .size()
        .reset_index(name="row_memberships")
    )
    membership["percent_of_benchmark"] = membership["row_memberships"].map(percent)
    membership = write_csv(
        membership.sort_values(["species", "row_memberships"], ascending=[True, False]),
        TABLE_DIR / "benchmark_dataset_kegg_like_module_membership.csv",
    )

    direct = df[df["direct_kegg_pathways"].fillna("").astype(str).str.len() > 0][
        ["entry_id", "species", "reaction_id", "reaction_name", "direct_kegg_pathways"]
    ].copy()
    direct["kegg_pathway_id"] = direct["direct_kegg_pathways"].str.split(";")
    direct = direct.explode("kegg_pathway_id")
    direct["kegg_pathway_name"] = direct["kegg_pathway_id"].map(pathway_name)
    direct_summary = (
        direct.groupby(["species", "kegg_pathway_id", "kegg_pathway_name"], dropna=False)
        .agg(rows=("entry_id", "count"), unique_reactions=("reaction_id", "nunique"))
        .reset_index()
        .sort_values("rows", ascending=False)
    )
    direct_summary["percent_of_benchmark"] = direct_summary["rows"].map(percent)
    direct_summary = write_csv(direct_summary, TABLE_DIR / "benchmark_dataset_direct_yeast_kegg_pathways.csv")
    return primary, membership, direct_summary


def make_method_technical_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    tech = pd.DataFrame(METHOD_TECHNICAL_ROWS)
    summary = pd.read_csv(TABLE_DIR / "method_eval_summary.csv")
    summary["coverage_percent"] = summary["n"] / BENCHMARK_N * 100
    group_map = tech.set_index("method")["benchmark_dimension"].to_dict()
    summary["group_cn"] = summary["method"].map(group_map).fillna("其他")
    write_csv(summary, TABLE_DIR / "method_eval_summary_annotated.csv")
    keep = [
        "method",
        "n",
        "coverage_percent",
        "group_cn",
        "mae_log10",
        "rmse_log10",
        "pearson_log10",
        "spearman_log10",
        "within_1.0_log10_fraction",
    ]
    keep = [c for c in keep if c in summary.columns]
    tech = tech.merge(summary[keep], on="method", how="left")
    tech["group_cn"] = tech["group_cn"].fillna(tech["benchmark_dimension"])
    tech["is_current_main_benchmark"] = True
    tech = write_csv(tech, TABLE_DIR / "method_technical_comparison.csv")

    dimensions = pd.DataFrame(
        [
            {
                "dimension": "输入覆盖",
                "plain_explanation_cn": f"方法能不能处理标准集的 {BENCHMARK_N} 行。缺失行可能来自模型输入域、反应信息、结构或官方资产限制。",
                "use_in_paper": f"报告覆盖率 n/{BENCHMARK_N}，并按广覆盖、reaction-aware、模型特定子集和公开重训版分开解释。",
            },
            {
                "dimension": "信息粒度",
                "plain_explanation_cn": "只看单个底物，还是看完整反应，或者还看蛋白结构。",
                "use_in_paper": "sequence+SMILES 方法可互相直接比较；reaction-aware 方法需要单独说明其信息更多但覆盖更窄。",
            },
            {
                "dimension": "模型来源",
                "plain_explanation_cn": "是官方权重直接推理，还是我们用公开数据重训/复现。",
                "use_in_paper": "DEKP-public-retrained 不能和官方最优权重画等号，应标成公开可复现版本。",
            },
            {
                "dimension": "AI 预测 vs 直接赋值",
                "plain_explanation_cn": "AI 模型会从输入特征中学习连续 kcat 数值；GO-HKP 这类方法则用功能相似性把已有 kcat 统计值赋给目标反应。",
                "use_in_paper": "GO-HKP 可作为非 AI 生物学基线，单独回答“简单功能赋值是否已经优于 AI 预测”。",
            },
            {
                "dimension": "评估指标",
                "plain_explanation_cn": "MAE/RMSE 看误差大小，Pearson/Spearman 看相关性，within10 看是否落在 10 倍误差内。",
                "use_in_paper": "主表同时给覆盖率和误差，避免只看某一个指标。",
            },
            {
                "dimension": "训练集重叠风险",
                "plain_explanation_cn": "如果测试样本和方法训练集重叠，指标可能虚高。",
                "use_in_paper": "后续写文章时应继续按序列、SMILES、sequence-SMILES pair 做查重标注。",
            },
            {
                "dimension": "生物学解释性",
                "plain_explanation_cn": "能否解释到反应、残基、底物原子或通路层面。",
                "use_in_paper": "PMAK/TurNuP 更适合讨论反应变化；结构方法可讨论结构资产，但需谨慎。",
            },
        ]
    )
    dimensions = write_csv(dimensions, TABLE_DIR / "method_comparison_dimensions.csv")
    return tech, dimensions


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def make_figures(
    df: pd.DataFrame,
    species_summary: pd.DataFrame,
    source_table: pd.DataFrame,
    ec_table: pd.DataFrame,
    pathway_primary: pd.DataFrame,
    top_reactions: pd.DataFrame,
    method_tech: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid", context="notebook")

    plt.figure(figsize=(6.5, 4.2))
    plot_df = species_summary.copy()
    sns.barplot(data=plot_df, x="species", y="rows", hue="species", dodge=False, palette="Set2", legend=False)
    for i, row in plot_df.reset_index(drop=True).iterrows():
        plt.text(i, row["rows"] + 8, f"{int(row['rows'])}\n{row['percent']:.1f}%", ha="center", va="bottom", fontsize=10)
    plt.xlabel("Species")
    plt.ylabel("Rows")
    plt.title("Benchmark rows by species")
    savefig(FIG_DIR / "species_distribution.png")

    pivot = source_table.pivot_table(index="species", columns="source_database", values="rows", fill_value=0)
    pivot = pivot.reindex(index=sorted(pivot.index))
    pivot.plot(kind="bar", stacked=True, figsize=(7.2, 4.6), colormap="tab20c")
    plt.xlabel("Species")
    plt.ylabel("Rows")
    plt.title("Experimental source by species")
    plt.legend(title="Source", bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(FIG_DIR / "source_by_species.png")

    plt.figure(figsize=(7.2, 4.4))
    sns.boxplot(data=df, x="species", y="true_kcat_log10", hue="species", palette="Set2", legend=False)
    sns.stripplot(data=df, x="species", y="true_kcat_log10", color="0.25", alpha=0.18, size=2)
    plt.xlabel("Species")
    plt.ylabel("log10(kcat / s^-1)")
    plt.title("Experimental kcat distribution")
    savefig(FIG_DIR / "kcat_log10_distribution_by_species.png")

    ec_plot = ec_table.copy()
    ec_plot["ec_class_short"] = ec_plot["ec_class"].str.replace(r" / .*", "", regex=True)
    plt.figure(figsize=(8.6, 4.8))
    sns.barplot(data=ec_plot, y="ec_class_short", x="row_memberships", hue="species", palette="Set2")
    plt.xlabel("Row memberships")
    plt.ylabel("EC class")
    plt.title("EC class coverage")
    plt.legend(title="Species")
    savefig(FIG_DIR / "ec_class_distribution.png")

    path_plot = pathway_primary.copy()
    path_plot["group"] = path_plot["kegg_like_primary_group_short"].map(lambda x: "\n".join(textwrap.wrap(str(x), 18)))
    plt.figure(figsize=(8.8, 5.0))
    sns.barplot(data=path_plot, y="group", x="rows", hue="species", palette="Set2")
    plt.xlabel("Rows")
    plt.ylabel("KEGG-like primary group")
    plt.title("KEGG-like module group distribution")
    plt.legend(title="Species")
    savefig(FIG_DIR / "kegg_like_group_by_species.png")

    top_plot = top_reactions.head(15).copy()
    top_plot["label"] = top_plot["species"] + " | " + top_plot["reaction_id"].astype(str)
    plt.figure(figsize=(8.2, 5.4))
    sns.barplot(data=top_plot, y="label", x="rows", hue="species", dodge=False, palette="Set2", legend=False)
    plt.xlabel("Rows")
    plt.ylabel("Reaction")
    plt.title("Top reaction records")
    savefig(FIG_DIR / "top_reactions.png")

    mt = method_tech.copy()
    mt["method_display"] = mt["method"].str.replace(r"-official$", "", regex=True)
    mt["group_cn"] = mt["group_cn"].replace(
        {"全量/近全量 sequence+SMILES 基线": "全量/近全量 sequence+SMILES"}
    )
    mt["group_cn"] = pd.Categorical(mt["group_cn"], METHOD_SCOPE_ORDER, ordered=True)
    mt = mt.sort_values(["group_cn", "coverage_percent", "method"], ascending=[True, False, True])
    scope_label_en = {
        "全量/近全量 sequence+SMILES": "Broad seq+SMILES",
        "reaction-aware 子集": "Reaction-aware subset",
        "模型特定子集": "Model-specific subset",
        "公开数据重训版": "Public-data retrained",
        "功能相似性 GO 赋值基线": "GO functional assignment",
    }
    mt["scope_plot"] = mt["group_cn"].astype(str).map(scope_label_en).fillna(mt["group_cn"].astype(str))
    palette = {
        "Broad seq+SMILES": "#4C78A8",
        "Reaction-aware subset": "#F58518",
        "Model-specific subset": "#54A24B",
        "Public-data retrained": "#E45756",
        "GO functional assignment": "#72B7B2",
    }
    plt.figure(figsize=(9.2, 5.8))
    sns.barplot(data=mt, y="method_display", x="coverage_percent", hue="scope_plot", dodge=False, palette=palette)
    plt.xlabel(f"Coverage of {BENCHMARK_N}-row benchmark (%)")
    plt.ylabel("Method")
    plt.title("Method coverage by comparison scope")
    plt.xlim(0, 105)
    plt.legend(title="Scope", bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(FIG_DIR / "method_coverage_by_scope.png")


def fmt_value(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (float, np.floating)):
        if 0 < abs(value) < 0.001:
            return f"{value:.2e}"
        if abs(value) >= 1000:
            return f"{value:,.2f}"
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return str(value)


def markdown_table(df: pd.DataFrame, columns: list[str] | None = None, max_rows: int = 20) -> str:
    if columns is not None:
        show = df[[c for c in columns if c in df.columns]].copy()
    else:
        show = df.copy()
    show = show.head(max_rows)
    headers = list(show.columns)
    rows = [[fmt_value(v) for v in row] for row in show.to_numpy()]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def make_report(
    df: pd.DataFrame,
    benchmark_funnel: pd.DataFrame,
    directory_map: pd.DataFrame,
    species_summary: pd.DataFrame,
    kcat_stats: pd.DataFrame,
    source_table: pd.DataFrame,
    match_table: pd.DataFrame,
    substrate_role: pd.DataFrame,
    ec_table: pd.DataFrame,
    top_reactions: pd.DataFrame,
    top_substrates: pd.DataFrame,
    pathway_primary: pd.DataFrame,
    direct_yeast: pd.DataFrame,
    method_tech: pd.DataFrame,
    dimensions: pd.DataFrame,
) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    total = len(df)
    ecoli_n = int((df["species"] == "ecoli").sum())
    yeast_n = int((df["species"] == "yeast").sum())
    currency_n = int(df["currency_or_cofactor_like_registry"].sum())
    currency_name_only_n = int(df["currency_or_cofactor_like_by_name"].sum())
    carrier_n = int(df["substrate_role_group"].eq("carrier_linked_variable").sum())
    reaction_n = df["reaction_id"].nunique()
    gene_n = df["gene_id"].nunique()
    substrate_n = df["substrate_name"].nunique()
    ec_strings_n = df["ec_number"].nunique()

    method_show = method_tech[
        [
            "method",
            "group_cn",
            "n",
            "coverage_percent",
            "input_needed",
            "model_family",
            "plain_language_cn",
            "main_caveat",
        ]
    ].copy()
    method_show = method_show.sort_values("coverage_percent", ascending=False)

    lines = [
        "# benchmark_ready_catpred 标准集数据画像与方法技术比较",
        "",
        f"生成时间：{generated}",
        "",
        "## 1. 先说结论",
        "",
        (
            f"`benchmark_ready_catpred.csv` 当前有 {total} 条可评测记录，其中大肠杆菌 {ecoli_n} 条"
            f"（{percent(ecoli_n):.1f}%），酿酒酵母 {yeast_n} 条（{percent(yeast_n):.1f}%）。"
            f"它覆盖 {reaction_n} 个模型反应、{gene_n} 个基因/UniProt、{substrate_n} 个底物名称和"
            f" {ec_strings_n} 种 EC 注释字符串。"
        ),
        "",
        (
            "这个文件名里的 `catpred` 是历史命名：CatPred 是第一个被打通的评测对象，"
            "但该文件现在不是 CatPred 专用数据，而是当前统一 benchmark 的 sequence+SMILES 标准集。"
        ),
        "",
        (
            f"按规范化名称、BiGG/KEGG/ChEBI/MetaNetX/PubChem 标识和标准化 InChIKey 的联合注册表规则，"
            f"{currency_n} 条为 currency/cofactor，{carrier_n} 条为可能仍是可变底物的 carrier-linked metabolite；"
            f"仅名称规则得到 {currency_name_only_n} 条，只作为探索性对照。角色不会在实验底物匹配前用于删行。"
        ),
        "",
        "## 2. 文件定位与字段含义",
        "",
        "- 标准集：`data/final/benchmark_ready_catpred.csv`",
        "- 行粒度：一行是一个 `酶/基因组反应/候选底物/实验 kcat` 记录。",
        "- 关键输入字段：`sequence` 是酶蛋白序列，`SMILES` 是底物结构字符串，`reaction_id` 是模型反应 ID。",
        "- 关键真值字段：`true_kcat` 是实验 kcat 原值，`true_kcat_log10` 是 log10 变换后的真值，评估主指标都在 log10 空间计算。",
        "- 溯源字段：`source_database`、`match_level`、`reference`、`n_measurements` 用于说明数据来自 BRENDA/SABIO-RK 以及匹配依据。",
        "",
        "输入、metadata、truth 的关系可以这样理解：`*_input.csv` 只给模型看；`*_metadata.csv` 记录每一行的物种、反应、来源和处理状态；`*_truth.csv` 只在评估时使用，避免把答案混进模型输入。",
        "",
        "## 3. 数据获取、清洗与 benchmark 确定方法",
        "",
        "### 3.1 从代谢模型定义候选酶-反应-底物条目",
        "",
        "- 大肠杆菌使用项目根目录的 `eciML1515.json`，酿酒酵母使用 `yeast-GEM.xml`。`src/01_parse_models.py` 解析反应、方向、GPR、EC、UniProt 和代谢物数据库编号。",
        "- GPR 是 gene-protein-reaction 规则，通俗说就是一条反应由哪些基因编码的酶负责。脚本把 `or` 拆成同工酶候选，把 `and` 保留为多亚基复合物。",
        "- 每行候选 entry 的粒度是 `物种 + 模型反应 + GPR 基因组 + 候选底物`。解析阶段保留模型编码的全部反应物，再由实验底物名称/标识决定匹配；currency/cofactor 角色只用于事后分层，不用于提前删行。",
        "",
        "### 3.2 补齐蛋白序列、小分子结构和反应结构",
        "",
        "- 蛋白序列：根据模型中的 UniProt accession，通过 UniProt REST 批量获取，并缓存为 `data/raw/uniprot_sequences.fasta` 和 `data/interim/uniprot_sequences.csv`。",
        "- 底物 SMILES：先使用模型已有注释，再通过 BiGG、KEGG、ChEBI、MetaNetX 等交叉编号在 CKB compound 数据库中映射；仍缺失时调用 PubChem PUG REST，并记录查询缓存和无法映射原因。",
        "- 完整 reaction SMILES：只为 TurNuP/PMAK 等 reaction-aware 方法准备，不作为所有方法进入统一 benchmark 的硬条件。",
        "- 蛋白结构：只为 DEKP 等结构感知方法收集，当前优先使用 AlphaFold/本地结构缓存，也不是统一 sequence+SMILES benchmark 的硬条件。",
        "",
        "### 3.3 获取并整理实验 kcat 真值",
        "",
        "- 主真值只使用 BRENDA turnover number 和 SABIO-RK kcat。早期推断或数据库填充得到的数值不进入统一 benchmark，也不作为当前公开仓库产物。",
        "- 仅保留目标物种、正的 kcat，统一单位为 `s^-1`；BRENDA 默认排除注释为 mutant/mutation/variant 的记录。范围值取区间均值。",
        "- 匹配先限定 `species + EC`，再比较底物数据库 ID/规范化名称以及 UniProt。优先级从高到低为：`species_ec_uniprot_substrate_id`、`species_ec_substrate_id`、`species_ec_uniprot_substrate_name`、`species_ec_substrate_name`。",
        "- 同一 entry 只保留最高匹配层级的实验记录；多条实验值在 kcat 原始尺度取中位数，再计算 `log10(kcat)`。pH 和温度也取可用记录的中位数，并保留来源、参考文献和测量条数。",
        "",
        "### 3.4 确定最终 benchmark",
        "",
        f"- `experimental_kcat_truth.csv` 是匹配到模型 entry 的实验真值全集，共 {len(pd.read_csv(BASE / 'data/final/experimental_kcat_truth.csv'))} 行。",
        f"- `benchmark_ready_truth.csv` 进一步要求 entry 能进入统一模型输入，即有单蛋白序列和可用底物 SMILES，共 {total} 行。",
        f"- `benchmark_ready_catpred.csv` 在这 {total} 行上合并 sequence、SMILES、真值和逐记录溯源字段。文件名保留 `catpred` 只是因为 CatPred 是第一个打通的方法，并不表示该标准集只服务于 CatPred。",
        f"- 其中 {int(df['experimental_substrate_support'].eq('substrate_supported').sum())} 条有 BRENDA 底物特异证据；其余 {int(df['experimental_substrate_support'].eq('participant_ambiguous').sum())} 条为 SABIO-RK participant-ambiguous 完整资源记录，单独做敏感性分析。",
        "- 方法评测时从这个母表提取各自需要的输入列，真值列只在推理结束后用于评分，避免答案泄漏到模型输入。",
        "",
        "从模型到最终 benchmark 的数量漏斗如下。注意 `enzyme_substrate_entries` 可以多于模型反应数，因为一条反应可能拆成多个基因组和多个候选底物。",
        "",
        markdown_table(benchmark_funnel, max_rows=10),
        "",
        "## 4. 物种、实验来源与匹配层级",
        "",
        markdown_table(
            species_summary,
            [
                "species",
                "rows",
                "percent",
                "unique_reactions",
                "unique_genes",
                "unique_substrates",
                "median_log10_kcat",
                "pH_available_rows",
                "temperature_available_rows",
                "currency_or_cofactor_like_rows_by_name",
                "currency_or_cofactor_rows_registry",
                "carrier_linked_variable_rows_registry",
            ],
        ),
        "",
        "实验来源按物种分布如下：",
        "",
        markdown_table(source_table, ["species", "source_database", "rows", "percent_of_benchmark"], max_rows=12),
        "",
        "匹配层级分布如下，`species_ec_uniprot_substrate_id` 通常代表物种、EC、UniProt 和底物 ID 都能对上，是最严格的一类匹配：",
        "",
        markdown_table(match_table, ["species", "match_level", "rows", "percent_of_benchmark"], max_rows=20),
        "",
        "相关图：",
        "",
        "- `reports/figures/kcat_dataset_context/species_distribution.png`",
        "- `reports/figures/kcat_dataset_context/source_by_species.png`",
        "- `reports/figures/kcat_dataset_context/kcat_log10_distribution_by_species.png`",
        "",
        "## 5. kcat 数值范围与反应分布",
        "",
        "kcat 跨度非常大，因此评估使用 log10 空间。`log10(kcat)=0` 表示 1 s^-1，`log10(kcat)=2` 表示 100 s^-1，`log10(kcat)=-2` 表示 0.01 s^-1。",
        "",
        markdown_table(kcat_stats, max_rows=5),
        "",
        "EC 大类分布如下。这里按 EC membership 统计：一条记录如果有多个 EC，可能贡献到多个类别。",
        "",
        markdown_table(ec_table, ["species", "ec_class", "row_memberships", "percent_of_benchmark"], max_rows=20),
        "",
        "出现次数最多的反应记录如下：",
        "",
        markdown_table(
            top_reactions,
            ["species", "reaction_id", "reaction_name", "rows", "unique_substrates", "unique_genes", "median_log10_kcat"],
            max_rows=20,
        ),
        "",
        "出现次数最多的底物名称如下：",
        "",
        markdown_table(
            top_substrates,
            ["substrate_name", "rows", "species_count", "unique_reactions", "substrate_role_group", "currency_or_cofactor_like_by_name", "median_kcat"],
            max_rows=20,
        ),
        "",
        "底物角色联合证据分布如下；carrier-linked 单列，不与严格 currency/cofactor 合并：",
        "",
        markdown_table(substrate_role, ["species", "substrate_role_group", "rows", "percent_of_benchmark"], max_rows=12),
        "",
        "相关图：",
        "",
        "- `reports/figures/kcat_dataset_context/ec_class_distribution.png`",
        "- `reports/figures/kcat_dataset_context/top_reactions.png`",
        "",
        "## 6. GO、KEGG-like 与通路/功能注释",
        "",
        "这里把 GO 和 KEGG 分开使用，避免把两类功能注释混成同一个概念：GO 更接近“蛋白/基因做什么功能”，KEGG 更接近“反应位于哪类代谢通路”。",
        "",
        "- GO-HKP 功能赋值：E. coli 使用 GO-HKP 自带的 iML1515R DeepGO-SE 反应级结果；yeast 使用 UniProt GO 注释补齐。两者都沿 `go-basic.obo` 的 GO 层级在 `GO_kcat_tree_total.csv` 中寻找可参考 kcat，并取 Total median 作为赋值。由于两物种 GO 来源不同，报告和 metadata 中分别标注来源。",
        "- 跨物种 KEGG-like 注释：使用 `DLKcat_official/DeeplearningApproach/Data/subsystem/module_ec.txt` 的 EC-to-module 功能大类。它不是直接 KEGG pathway ID，但可用同一口径比较 E. coli 和 yeast。",
        "- yeast 直接 KEGG pathway：解析 `yeast-GEM.xml` 中反应自带的 `kegg.pathway` ID。E. coli 的 `eciML1515.json` 只有 KEGG reaction ID，没有系统的 pathway 字段，因此当前不能用同样方式做直接 pathway 统计。",
        "- 这些注释用于描述标准集覆盖和构建 GO 赋值基线，不参与其他 AI 模型的真值筛选，也不会改变实验 kcat。",
        "",
        "按 EC 推断的 KEGG-like 主功能大类如下，每行只归到一个主类，因此合计等于各物种样本数：",
        "",
        markdown_table(
            pathway_primary,
            ["species", "kegg_like_primary_group_short", "rows", "percent_of_species", "percent_of_benchmark"],
            max_rows=20,
        ),
        "",
        "yeast-GEM 直接 KEGG pathway ID 的 top 分布如下。`sce01100/sce01110/sce01130` 这类全局通路会覆盖很多反应，解释时应更关注具体代谢通路，例如碳代谢、氨基酸生物合成、嘌呤/嘧啶代谢等。",
        "",
        markdown_table(
            direct_yeast,
            ["kegg_pathway_id", "kegg_pathway_name", "rows", "unique_reactions", "percent_of_benchmark"],
            max_rows=20,
        ),
        "",
        "相关图：",
        "",
        "- `reports/figures/kcat_dataset_context/kegg_like_group_by_species.png`",
        "",
        "## 7. 项目目录结构与分析类型",
        "",
        "下面按“目录承担什么工作”整理项目结构。`publication_status` 用来区分 GitHub 已公开内容、本地 gitignored 运行目录，以及需要从 Zenodo 或上游来源恢复的资产。GitHub 中的 `data/final/` 只跟踪三张核心 benchmark CSV，`data/final/<method>/` 属于本地方法级产物。",
        "",
        markdown_table(directory_map, max_rows=30),
        "",
        "重点分析文件可以快速定位为：",
        "",
        "- 数据准备、方法输入、预测和 Slurm 队列入口：`scripts/runners/`。",
        "- GO 分析：`external_methods/GO-HKP/`、`data/raw/go_hkp/`、`data/final/go_hkp/`、`reports/tables/go_hkp_*`。",
        "- KEGG/EC/通路分析：`reports/tables/benchmark_dataset_kegg_like_*`、`reports/tables/benchmark_dataset_direct_yeast_kegg_pathways.csv`。",
        "- MAE/RMSE/bias/within-fold 分析：`reports/tables/*_eval_metrics.csv` 和 `reports/figures/kcat_benchmark_summary/`。",
        "- Species-level 分析：`reports/tables/species_mae_matrix.csv`、`reports/tables/benchmark_dataset_*_by_species.csv` 和 species heatmap。",
        "- Method-level 分析：`reports/tables/method_eval_summary*.csv`、`reports/tables/method_rank*.csv`、`reports/tables/method_technical_comparison.csv`。",
        "- 写文章用独立表格：`reports/report_tables/`，其中 `manifest.csv` 记录每张表的来源。",
        "",
        "## 8. 不同预测方法的技术原理与比较维度",
        "",
        "下面这张表把“模型看了什么信息”和“它适合在哪个维度比较”放在一起。通俗地说，sequence+SMILES 方法是只看酶和单个底物；reaction-aware 方法还看产物，因此信息更多但需要更完整的数据；结构感知方法还看蛋白结构，但公开复现难度也更高；GO-HKP 则不是 AI 回归模型，而是用 GO 功能相似性做直接 kcat 赋值。本项目里 GO-HKP 的 E. coli 部分来自本地 DeepGO-SE 反应赋值，yeast 部分用 UniProt GO 注释补齐。",
        "",
        markdown_table(method_show, max_rows=20),
        "",
        "比较维度建议如下：",
        "",
        markdown_table(dimensions, max_rows=10),
        "",
        "相关图：",
        "",
        "- `reports/figures/kcat_dataset_context/method_coverage_by_scope.png`",
        "",
        "## 9. 写文章时建议怎么表述",
        "",
        "1. 主文可以把 `benchmark_ready_catpred.csv` 称为 unified sequence+SMILES benchmark，而不是 CatPred 专用输入。",
        "2. 结果表必须同时给覆盖率和误差指标；否则 KinForm、CatPred、PMAK/TurNuP 这种子集方法会和全量方法混在一起，容易误导。",
        "3. 通路分布建议分两句话写：跨物种用 EC-to-module 的 KEGG-like 功能大类；酿酒酵母另有直接 KEGG pathway ID 作为补充。",
        "4. 如果要做更严格的生物学解释，下一步应补充 E. coli 的 KEGG reaction-to-pathway 映射，或用 BioCyc/MetaCyc 子系统统一重注释两套模型。",
        "",
        "## 10. 输出文件清单",
        "",
        "- `reports/tables/benchmark_dataset_species_summary.csv`",
        "- `reports/tables/benchmark_dataset_kcat_stats_by_species.csv`",
        "- `reports/tables/benchmark_dataset_source_by_species.csv`",
        "- `reports/tables/benchmark_dataset_match_level_by_species.csv`",
        "- `reports/tables/benchmark_dataset_ec_class_summary.csv`",
        "- `reports/tables/benchmark_dataset_top_reactions.csv`",
        "- `reports/tables/benchmark_dataset_top_substrates.csv`",
        "- `reports/tables/benchmark_dataset_kegg_like_primary_group.csv`",
        "- `reports/tables/benchmark_dataset_direct_yeast_kegg_pathways.csv`",
        "- `reports/tables/method_technical_comparison.csv`",
        "- `reports/tables/method_comparison_dimensions.csv`",
        "- `reports/tables/benchmark_build_funnel.csv`",
        "- `reports/tables/project_directory_analysis_map.csv`",
        "- `reports/report_tables/`",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    df = load_enriched_benchmark()
    write_csv(df, TABLE_DIR / "benchmark_ready_catpred_enriched_context.csv")

    benchmark_funnel = make_benchmark_build_funnel(df)
    directory_map = make_project_directory_map()
    species_summary = make_species_summary(df)
    kcat_stats = make_kcat_stats(df)
    source_table = value_count_table(df, ["species", "source_database"], "source_by_species")
    match_table = value_count_table(df, ["species", "match_level"], "match_level_by_species")
    value_count_table(df, ["species", "enzyme_complex_type"], "enzyme_complex_type_by_species")
    substrate_role = value_count_table(df, ["species", "substrate_role_group"], "substrate_role_by_species")
    ec_table = make_ec_class_summary(df)
    top_reactions = make_top_reactions(df)
    top_substrates = make_top_substrates(df)
    pathway_primary, _, direct_yeast = make_pathway_tables(df)
    method_tech, dimensions = make_method_technical_tables()

    make_figures(df, species_summary, source_table, ec_table, pathway_primary, top_reactions, method_tech)
    make_report(
        df,
        benchmark_funnel,
        directory_map,
        species_summary,
        kcat_stats,
        source_table,
        match_table,
        substrate_role,
        ec_table,
        top_reactions,
        top_substrates,
        pathway_primary,
        direct_yeast,
        method_tech,
        dimensions,
    )
    print(f"Wrote report: {REPORT_PATH}")
    print(f"Wrote tables under: {TABLE_DIR}")
    print(f"Wrote figures under: {FIG_DIR}")


if __name__ == "__main__":
    main()
