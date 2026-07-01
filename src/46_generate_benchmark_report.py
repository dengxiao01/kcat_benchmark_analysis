#!/usr/bin/env python3
"""Generate summary tables, figures, and a Chinese benchmark analysis report."""

from __future__ import annotations

from datetime import date
from pathlib import Path
import textwrap

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BASE = Path(__file__).resolve().parent.parent
REPORT_DIR = BASE / "reports"
TABLE_DIR = REPORT_DIR / "tables"
FIGURE_DIR = REPORT_DIR / "figures" / "kcat_benchmark_summary"
REPORT_PATH = REPORT_DIR / "kcat_benchmark_analysis_report.md"
BENCHMARK_N = 978


METHOD_META = {
    "DLKcat-official": {
        "scope": "Broad sequence+SMILES",
        "group_cn": "全量/近全量 sequence+SMILES",
        "group_note": "输入为酶序列和单底物 SMILES；除无效 SMILES 外基本覆盖整个 benchmark。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 977/978，缺 1 条非法 SMILES。",
        "row_file": BASE / "data/final/dlkcat/dlkcat_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "dlkcat_official_eval_metrics.csv",
    },
    "UniKP-official": {
        "scope": "Broad sequence+SMILES",
        "group_cn": "全量/近全量 sequence+SMILES",
        "group_note": "输入为酶序列和单底物 SMILES；除无效 SMILES 外基本覆盖整个 benchmark。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 977/978，缺 1 条非法 SMILES。",
        "row_file": BASE / "data/final/unikp/unikp_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "unikp_official_eval_metrics.csv",
    },
    "TurNuP-official": {
        "scope": "Reaction-aware subset",
        "group_cn": "reaction-aware 子集",
        "group_note": "输入需要完整反应信息，即底物侧 SMILES、产物侧 SMILES 和酶序列。",
        "role": "current",
        "modality": "reaction + enzyme",
        "coverage_note": "覆盖 780/978；缺失来自未补齐完整 reaction SMILES 的记录。",
        "row_file": BASE / "data/final/turnup/turnup_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "turnup_eval_metrics.csv",
    },
    "CatPred": {
        "scope": "Model-specific subset",
        "group_cn": "模型特定子集",
        "group_note": "方法输入或官方推理流程有额外限制，因此只在可被该模型有效处理的子集上评估。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 913/978；缺失来自 CatPred 官方流程可处理范围。",
        "row_file": BASE / "data/final/catpred/catpred_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "catpred_eval_metrics.csv",
    },
    "CataPro": {
        "scope": "Broad sequence+SMILES",
        "group_cn": "全量/近全量 sequence+SMILES",
        "group_note": "输入为酶序列和单底物 SMILES；除无效 SMILES 外基本覆盖整个 benchmark。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 977/978，缺 1 条非法 SMILES。",
        "row_file": BASE / "data/final/catapro/catapro_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "catapro_eval_metrics.csv",
    },
    "PMAK": {
        "scope": "Reaction-aware subset",
        "group_cn": "reaction-aware 子集",
        "group_note": "输入需要完整反应信息，即底物侧 SMILES、产物侧 SMILES 和酶序列。",
        "role": "current",
        "modality": "reaction + enzyme",
        "coverage_note": "覆盖 780/978；缺失来自未补齐完整 reaction SMILES 的记录。",
        "row_file": BASE / "data/final/pmak/pmak_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "pmak_eval_metrics.csv",
    },
    "KinForm": {
        "scope": "Model-specific subset",
        "group_cn": "模型特定子集",
        "group_note": "方法输入或官方推理流程有额外限制，因此只在可被该模型有效处理的子集上评估。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 563/978；主要受 KinForm 可处理输入范围限制。",
        "row_file": BASE / "data/final/kinform/kinform_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "kinform_eval_metrics.csv",
    },
    "KcatNet": {
        "scope": "Broad sequence+SMILES",
        "group_cn": "全量/近全量 sequence+SMILES",
        "group_note": "输入为酶序列和单底物 SMILES；除无效 SMILES 外基本覆盖整个 benchmark。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 977/978，缺 1 条非法 SMILES。",
        "row_file": BASE / "data/final/kcatnet/kcatnet_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "kcatnet_eval_metrics.csv",
    },
    "PreTKcat": {
        "scope": "Broad sequence+SMILES",
        "group_cn": "全量/近全量 sequence+SMILES",
        "group_note": "输入为酶序列和单底物 SMILES；除无效 SMILES 外基本覆盖整个 benchmark。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 977/978，缺 1 条非法 SMILES。",
        "row_file": BASE / "data/final/pretkcat/pretkcat_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "pretkcat_eval_metrics.csv",
    },
    "DEKP-public-retrained": {
        "scope": "Public-data retrained",
        "group_cn": "公开数据重训版",
        "group_note": "不是原论文官方最优权重，而是用公开数据和当前可复现流程重新训练/补齐后的版本。",
        "role": "current_retrained",
        "modality": "sequence + substrate SMILES + structure",
        "coverage_note": "覆盖 977/978；缺 1 条非法 SMILES，同时需要结构/图特征补齐。",
        "row_file": BASE / "data/final/dekp/dekp_public_retrained_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "dekp_public_retrained_eval_metrics.csv",
    },
    "SELFprot": {
        "scope": "Broad sequence+SMILES",
        "group_cn": "全量/近全量 sequence+SMILES",
        "group_note": "输入为酶序列和单底物 SMILES；除无效 SMILES 外基本覆盖整个 benchmark。",
        "role": "current",
        "modality": "sequence + substrate SMILES",
        "coverage_note": "覆盖 977/978，缺 1 条非法 SMILES。",
        "row_file": BASE / "data/final/selfprot/selfprot_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "selfprot_eval_metrics.csv",
    },
    "GO-HKP": {
        "scope": "Function-assignment GO baseline",
        "group_cn": "功能相似性 GO 赋值基线",
        "group_note": "不是深度回归模型，而是基于 GO 功能层级给反应/基因赋参考 kcat；E. coli 使用本地 DeepGO-SE 反应赋值，yeast 使用 UniProt GO 注释补齐。",
        "role": "current_baseline",
        "modality": "GO hierarchy + functional assignment",
        "coverage_note": "覆盖 978/978；E. coli 为 GO-HKP DeepGO-SE 反应级赋值，yeast 为 UniProt GO 注释的 GOATOOLS-style 补充赋值。",
        "row_file": BASE / "data/final/go_hkp/go_hkp_kcat_predictions_evaluated.csv",
        "metric_file": TABLE_DIR / "go_hkp_eval_metrics.csv",
    },
}

CURRENT_METHODS = [
    "DLKcat-official",
    "UniKP-official",
    "TurNuP-official",
    "CatPred",
    "CataPro",
    "PMAK",
    "KinForm",
    "KcatNet",
    "PreTKcat",
    "DEKP-public-retrained",
    "SELFprot",
    "GO-HKP",
]

MAIN_SCATTER_METHODS = [
    "KcatNet",
    "CataPro",
    "TurNuP-official",
    "PMAK",
    "CatPred",
    "KinForm",
    "DEKP-public-retrained",
    "GO-HKP",
]

SCOPE_COLORS = {
    "Broad sequence+SMILES": "#4C78A8",
    "Reaction-aware subset": "#F58518",
    "Model-specific subset": "#54A24B",
    "Public-data retrained": "#E45756",
    "Function-assignment GO baseline": "#72B7B2",
}


def method_color(method: str) -> str:
    return SCOPE_COLORS[METHOD_META[method]["scope"]]


def load_summary() -> pd.DataFrame:
    summary = pd.read_csv(TABLE_DIR / "method_eval_summary.csv")
    summary["coverage_fraction"] = summary["n"] / BENCHMARK_N
    summary["coverage_percent"] = summary["coverage_fraction"] * 100
    summary["scope"] = summary["method"].map(lambda m: METHOD_META.get(m, {}).get("scope", "Other"))
    summary["group_cn"] = summary["method"].map(lambda m: METHOD_META.get(m, {}).get("group_cn", "其他"))
    summary["group_note"] = summary["method"].map(lambda m: METHOD_META.get(m, {}).get("group_note", ""))
    summary["coverage_note"] = summary["method"].map(lambda m: METHOD_META.get(m, {}).get("coverage_note", ""))
    summary["role"] = summary["method"].map(lambda m: METHOD_META.get(m, {}).get("role", "other"))
    summary["modality"] = summary["method"].map(lambda m: METHOD_META.get(m, {}).get("modality", ""))
    summary["is_current"] = summary["method"].isin(CURRENT_METHODS)
    summary.to_csv(TABLE_DIR / "method_eval_summary_annotated.csv", index=False)
    return summary


def build_method_annotation(summary: pd.DataFrame) -> pd.DataFrame:
    order = CURRENT_METHODS
    ann = summary[summary["method"].isin(CURRENT_METHODS)].copy()
    ann["method"] = pd.Categorical(ann["method"], categories=order, ordered=True)
    ann = ann.sort_values("method")
    columns = [
        "method",
        "group_cn",
        "scope",
        "role",
        "modality",
        "n",
        "coverage_percent",
        "coverage_note",
        "group_note",
    ]
    ann = ann[columns].copy()
    ann.to_csv(TABLE_DIR / "method_group_annotation.csv", index=False)
    return ann


def load_rows(methods: list[str]) -> pd.DataFrame:
    frames = []
    for method in methods:
        path = METHOD_META[method]["row_file"]
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if "abs_error_log10" not in df.columns:
            df["abs_error_log10"] = (df["prediction_log10"] - df["true_kcat_log10"]).abs()
        df["method"] = method
        df["scope"] = METHOD_META[method]["scope"]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_group_metrics(methods: list[str], group_type: str, value: str = "mae_log10") -> pd.DataFrame:
    rows = []
    for method in methods:
        path = METHOD_META[method]["metric_file"]
        if not path.exists():
            continue
        df = pd.read_csv(path)
        part = df[df["group_type"].eq(group_type)].copy()
        for _, row in part.iterrows():
            rows.append({"method": method, "group": row["group"], value: row.get(value, np.nan)})
    return pd.DataFrame(rows)


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def plot_overall_error(summary: pd.DataFrame) -> Path:
    current = summary[summary["method"].isin(CURRENT_METHODS)].copy()
    current = current.sort_values("mae_log10", ascending=True)
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    y = np.arange(len(current))
    height = 0.36
    colors = [method_color(m) for m in current["method"]]
    ax.barh(y - height / 2, current["mae_log10"], height, color=colors, alpha=0.95, label="MAE")
    ax.barh(y + height / 2, current["rmse_log10"], height, color=colors, alpha=0.45, label="RMSE")
    ax.set_yticks(y)
    ax.set_yticklabels(current["method"])
    ax.invert_yaxis()
    ax.set_xlabel("Error on log10(kcat)")
    ax.set_title("Overall Error, Current Benchmark Methods")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(loc="lower right")
    return FIGURE_DIR / "overall_error_mae_rmse.png"


def plot_correlation(summary: pd.DataFrame) -> Path:
    current = summary[summary["method"].isin(CURRENT_METHODS)].copy()
    current = current.sort_values("spearman_log10", ascending=False)
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    x = np.arange(len(current))
    width = 0.38
    colors = [method_color(m) for m in current["method"]]
    ax.bar(x - width / 2, current["pearson_log10"], width, color=colors, alpha=0.95, label="Pearson")
    ax.bar(x + width / 2, current["spearman_log10"], width, color=colors, alpha=0.45, label="Spearman")
    ax.set_xticks(x)
    ax.set_xticklabels(current["method"], rotation=50, ha="right")
    ax.set_ylabel("Correlation on log10(kcat)")
    ax.set_title("Correlation, Current Benchmark Methods")
    ax.set_ylim(0, max(0.7, float(current[["pearson_log10", "spearman_log10"]].max().max()) + 0.05))
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    return FIGURE_DIR / "overall_correlation.png"


def plot_coverage_vs_mae(summary: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(12.2, 6.8))
    plot_df = summary.copy()
    label_offsets = {
        "DEKP-public-retrained": (-10, -10),
        "DLKcat-official": (8, 7),
        "SELFprot": (8, -9),
        "UniKP-official": (8, 0),
        "PreTKcat": (8, 8),
        "CataPro": (8, 5),
        "KcatNet": (8, -7),
        "CatPred": (8, 5),
        "PMAK": (8, 6),
        "TurNuP-official": (8, -9),
        "GO-HKP": (8, 6),
    }
    for scope, part in plot_df.groupby("scope"):
        ax.scatter(
            part["coverage_percent"],
            part["mae_log10"],
            s=np.sqrt(part["n"].clip(lower=50)) * 30,
            alpha=0.78,
            color=SCOPE_COLORS.get(scope, "#888888"),
            label=scope,
            edgecolor="white",
            linewidth=0.8,
        )
        for _, row in part.iterrows():
            offset = label_offsets.get(row["method"], (5, 3))
            ax.annotate(
                row["method"],
                (row["coverage_percent"], row["mae_log10"]),
                xytext=offset,
                textcoords="offset points",
                fontsize=8,
                ha="right" if offset[0] < 0 else "left",
            )
    ax.set_xlabel("Coverage of 978 benchmark rows (%)")
    ax.set_ylabel("MAE on log10(kcat), lower is better")
    ax.set_title("Accuracy-Coverage Tradeoff")
    ax.set_xlim(0, 108)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=True, borderaxespad=0)
    return FIGURE_DIR / "coverage_vs_mae.png"


def plot_within10_bias(summary: pd.DataFrame) -> Path:
    current = summary[summary["method"].isin(CURRENT_METHODS)].copy()
    current = current.sort_values("within_1.0_log10_fraction", ascending=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), gridspec_kw={"width_ratios": [1.05, 1]})
    y = np.arange(len(current))
    colors = [method_color(m) for m in current["method"]]
    axes[0].barh(y, current["within_1.0_log10_fraction"], color=colors, alpha=0.9)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(current["method"])
    axes[0].set_xlabel("Fraction within 10-fold error")
    axes[0].set_xlim(0, 1)
    axes[0].grid(axis="x", alpha=0.25)

    bias_colors = np.where(current["bias_log10"] < 0, "#D65F5F", "#4C78A8")
    axes[1].barh(y, current["bias_log10"], color=bias_colors, alpha=0.9)
    axes[1].axvline(0, color="black", linewidth=0.8)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("Bias on log10(kcat)")
    axes[1].grid(axis="x", alpha=0.25)
    fig.suptitle("Practical Accuracy and Bias")
    return FIGURE_DIR / "within10_and_bias.png"


def plot_error_distribution(rows: pd.DataFrame, summary: pd.DataFrame) -> Path:
    current_rows = rows[rows["method"].isin(CURRENT_METHODS)].copy()
    order = (
        summary[summary["method"].isin(CURRENT_METHODS)]
        .sort_values("median_abs_error_log10")
        ["method"]
        .tolist()
    )
    fig, ax = plt.subplots(figsize=(11.5, 6.4))
    palette = {method: method_color(method) for method in order}
    sns.boxplot(
        data=current_rows,
        x="method",
        y="abs_error_log10",
        hue="method",
        order=order,
        palette=palette,
        ax=ax,
        showfliers=False,
        width=0.65,
        legend=False,
    )
    sns.stripplot(
        data=current_rows.sample(min(len(current_rows), 4000), random_state=42),
        x="method",
        y="abs_error_log10",
        order=order,
        ax=ax,
        color="black",
        size=1.2,
        alpha=0.18,
        jitter=0.22,
    )
    ax.axhline(0.3, color="#777777", linestyle="--", linewidth=1, label="~2-fold")
    ax.axhline(1.0, color="#333333", linestyle=":", linewidth=1.2, label="10-fold")
    ax.set_xlabel("")
    ax.set_ylabel("Absolute error on log10(kcat)")
    ax.set_title("Per-row Error Distribution")
    ax.tick_params(axis="x", rotation=50)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    return FIGURE_DIR / "error_distribution_boxplot.png"


def plot_species_heatmap() -> Path:
    df = load_group_metrics(CURRENT_METHODS, "species")
    pivot = df.pivot(index="method", columns="group", values="mae_log10")
    method_order = [m for m in CURRENT_METHODS if m in pivot.index]
    pivot = pivot.loc[method_order]
    pivot.to_csv(TABLE_DIR / "species_mae_matrix.csv")
    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".3f",
        cmap="YlGnBu_r",
        linewidths=0.5,
        cbar_kws={"label": "MAE log10"},
        ax=ax,
    )
    ax.set_xlabel("Species")
    ax.set_ylabel("")
    ax.set_title("Species-level MAE")
    return FIGURE_DIR / "species_mae_heatmap.png"


def plot_source_heatmap() -> Path:
    df = load_group_metrics(CURRENT_METHODS, "source_database")
    pivot = df.pivot(index="method", columns="group", values="mae_log10")
    method_order = [m for m in CURRENT_METHODS if m in pivot.index]
    pivot = pivot.loc[method_order]
    pivot.to_csv(TABLE_DIR / "source_database_mae_matrix.csv")
    fig, ax = plt.subplots(figsize=(8.4, 6.7))
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        linewidths=0.5,
        cbar_kws={"label": "MAE log10"},
        ax=ax,
    )
    ax.set_xlabel("Source database")
    ax.set_ylabel("")
    ax.set_title("Source-database-level MAE")
    return FIGURE_DIR / "source_database_mae_heatmap.png"


def plot_scatter_selected(rows: pd.DataFrame, summary: pd.DataFrame) -> Path:
    selected_rows = rows[rows["method"].isin(MAIN_SCATTER_METHODS)].copy()
    vmin = np.nanpercentile(
        selected_rows[["true_kcat_log10", "prediction_log10"]].to_numpy().ravel(),
        1,
    )
    vmax = np.nanpercentile(
        selected_rows[["true_kcat_log10", "prediction_log10"]].to_numpy().ravel(),
        99,
    )
    pad = 0.25
    vmin -= pad
    vmax += pad
    n_methods = len(MAIN_SCATTER_METHODS)
    ncols = 3 if n_methods <= 9 else 4
    nrows = int(np.ceil(n_methods / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.4 * nrows), sharex=True, sharey=True)
    axes = np.array(axes).ravel()
    metric_lookup = summary.set_index("method")
    for idx, (ax, method) in enumerate(zip(axes, MAIN_SCATTER_METHODS)):
        part = selected_rows[selected_rows["method"].eq(method)]
        ax.scatter(
            part["true_kcat_log10"],
            part["prediction_log10"],
            s=9,
            alpha=0.35,
            color=method_color(method),
            edgecolors="none",
        )
        x = np.linspace(vmin, vmax, 100)
        ax.plot(x, x, color="black", linewidth=0.9)
        ax.plot(x, x + 1, color="#777777", linewidth=0.7, linestyle=":")
        ax.plot(x, x - 1, color="#777777", linewidth=0.7, linestyle=":")
        row = metric_lookup.loc[method]
        ax.set_title(f"{method}\nn={int(row['n'])}, MAE={row['mae_log10']:.3f}", fontsize=10)
        ax.grid(alpha=0.2)
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)
    for ax in axes[n_methods:]:
        ax.axis("off")
    for idx, ax in enumerate(axes[:n_methods]):
        ax.set_xlabel("True log10(kcat)" if idx >= (nrows - 1) * ncols else "")
        ax.set_ylabel("Predicted log10(kcat)" if idx % ncols == 0 else "")
    fig.suptitle("Predicted vs True kcat, Representative Methods")
    return FIGURE_DIR / "predicted_vs_true_selected.png"


def build_rank_tables(summary: pd.DataFrame) -> pd.DataFrame:
    current = summary[summary["method"].isin(CURRENT_METHODS)].copy()
    current = current.sort_values(["mae_log10", "rmse_log10"], ascending=True)
    rank_cols = [
        "method",
        "scope",
        "modality",
        "group_cn",
        "group_note",
        "n",
        "coverage_percent",
        "mae_log10",
        "rmse_log10",
        "pearson_log10",
        "spearman_log10",
        "bias_log10",
        "within_1.0_log10_fraction",
    ]
    current[rank_cols].to_csv(TABLE_DIR / "method_rank_current_benchmark.csv", index=False)
    return current


def fig_link(path: Path) -> str:
    return str(path.relative_to(REPORT_DIR))


def fmt_pct(value: float) -> str:
    return f"{value * 100:.1f}%"


def write_report(
    summary: pd.DataFrame,
    current_rank: pd.DataFrame,
    method_annotation: pd.DataFrame,
    fig_paths: dict[str, Path],
) -> None:
    top_mae = current_rank.iloc[0]
    broad = current_rank[current_rank["scope"].eq("Broad sequence+SMILES")]
    top_broad = broad.iloc[0]
    top_corr = summary[summary["method"].isin(CURRENT_METHODS)].sort_values("spearman_log10", ascending=False).iloc[0]
    go_note = ""
    go_rows = summary[summary["method"].eq("GO-HKP")]
    if not go_rows.empty:
        go_row = go_rows.iloc[0]
        go_note = (
            "\n\nGO-HKP 是“功能相似性直接赋值”非 AI 基线：当前覆盖 "
            f"{int(go_row['n'])}/978 条记录，E. coli 使用本地 GO-HKP DeepGO-SE 反应级赋值，"
            "yeast 使用 UniProt GO 注释做 GOATOOLS-style 补充赋值。整体 MAE log10 为 "
            f"{go_row['mae_log10']:.4f}，"
            f"Spearman 为 {go_row['spearman_log10']:.4f}，10 倍误差内比例为 "
            f"{go_row['within_1.0_log10_fraction']:.4f}，bias 为 {go_row['bias_log10']:.4f}。"
            "这说明 GO 功能赋值可以作为有解释性的非 AI 基线，但整体仍没有优于最强的 AI 预测方法，"
            "并且仍有偏高估趋势。"
        )

    current_table = current_rank[
        [
            "method",
            "group_cn",
            "modality",
            "n",
            "coverage_percent",
            "mae_log10",
            "rmse_log10",
            "pearson_log10",
            "spearman_log10",
            "within_1.0_log10_fraction",
            "bias_log10",
        ]
    ].copy()
    for col in ["mae_log10", "rmse_log10", "pearson_log10", "spearman_log10", "within_1.0_log10_fraction", "bias_log10"]:
        current_table[col] = current_table[col].map(lambda x: f"{x:.4f}")
    current_table["coverage_percent"] = current_table["coverage_percent"].map(lambda x: f"{x:.1f}%")
    current_table_md = current_table.to_markdown(index=False)

    group_def = pd.DataFrame(
        [
            {
                "分组": "全量/近全量 sequence+SMILES",
                "判定标准": "输入主要是酶序列和单底物 SMILES，除 1 条非法 SMILES 外基本能覆盖 978 条 benchmark。",
                "方法": "DLKcat-official, UniKP-official, CataPro, KcatNet, PreTKcat, SELFprot",
            },
            {
                "分组": "reaction-aware 子集",
                "判定标准": "模型需要完整反应信息，即底物侧、产物侧和酶序列；当前只有 780 条补齐了 reaction SMILES。",
                "方法": "TurNuP-official, PMAK",
            },
            {
                "分组": "模型特定子集",
                "判定标准": "方法官方推理流程或输入限制导致只能在该模型可处理子集上评估。",
                "方法": "CatPred, KinForm",
            },
            {
                "分组": "公开数据重训版",
                "判定标准": "不是原论文官方最优权重，而是用公开数据和当前可复现流程重新训练/补齐后的版本。",
                "方法": "DEKP-public-retrained",
            },
            {
                "分组": "功能相似性 GO 赋值基线",
                "判定标准": "不训练深度回归模型，而是用 GO 功能层级和已有 GO-kcat 统计值给反应/基因直接赋 kcat；当前 E. coli 和 yeast 的 GO 来源不同，需在正文中标注。",
                "方法": "GO-HKP",
            },
        ]
    )
    group_def_md = group_def.to_markdown(index=False)

    ann_table = method_annotation[
        [
            "method",
            "group_cn",
            "modality",
            "n",
            "coverage_percent",
            "coverage_note",
        ]
    ].copy()
    ann_table["coverage_percent"] = ann_table["coverage_percent"].map(lambda x: f"{x:.1f}%")
    ann_table_md = ann_table.to_markdown(index=False)

    report = f"""# kcat 预测方法统一评测分析报告

生成日期：{date.today().isoformat()}

## 一句话结论

在当前 978 条 *E. coli* + yeast enzyme-substrate 实验 kcat benchmark 上，若按当前可比的正式评测结果看，`{top_broad['method']}` 是“全量/近全量 sequence+SMILES”方法中误差最低的方法，MAE log10 为 {top_broad['mae_log10']:.4f}，覆盖 {int(top_broad['n'])} 条；`{top_mae['method']}` 在所有当前方法中 MAE 最低，MAE log10 为 {top_mae['mae_log10']:.4f}，它所属分组为“{top_mae['group_cn']}”；`{top_corr['method']}` 的 Spearman 相关性最高，为 {top_corr['spearman_log10']:.4f}。

这里的 MAE log10 可以通俗理解为“预测值和实验值差了几个 10 倍单位”。例如误差 1.0 表示大约差 10 倍，误差 0.3 表示大约差 2 倍。
{go_note}

## 评测口径

- 真值集合：`data/final/benchmark_ready_catpred.csv`，共 978 条实验 kcat 记录。
- 统一指标：MAE、RMSE、Pearson、Spearman、bias，以及误差在 10 倍以内的比例，全部在 log10(kcat) 尺度上计算。
- 当前正式评测方法：`DLKcat-official`、`UniKP-official`、`TurNuP-official`、`CatPred`、`CataPro`、`PMAK`、`KinForm`、`KcatNet`、`PreTKcat`、`DEKP-public-retrained`、`SELFprot`、`GO-HKP`。
- 其中 `GO-HKP` 是功能相似性直接赋值基线，不是 AI 回归模型；它用 GO 层级把功能相近的酶/反应归到可参考的 kcat 统计值上，用来回答“直接赋值是否已经足够强”这个问题。本项目中 E. coli 用 GO-HKP 已有 DeepGO-SE 结果，yeast 用 UniProt GO 注释补齐。

## 分组定义与方法归属

这些分组不是按模型名字主观划分，而是按“输入信息是否一致、覆盖范围是否一致、权重来源是否一致、是否为 AI 回归模型”来划分。通俗说，只有输入口径和覆盖范围接近的方法，才适合直接横向排名。

{group_def_md}

每个方法的详细标注如下：

{ann_table_md}

## 总体结果

{current_table_md}

## 图表解读

### 1. 整体误差：MAE/RMSE

![Overall error]({fig_link(fig_paths['overall_error'])})

MAE 更接近日常理解中的“平均偏差”，RMSE 会对特别大的错误惩罚更重。`KcatNet` 在覆盖 977 条的情况下 MAE 最低；`TurNuP-official` 和 `PMAK` 的 MAE 也低，但它们只覆盖完整 reaction SMILES 的 780 条，因此更适合作为 reaction-aware 子集比较。

### 2. 排序相关性：Pearson/Spearman

![Correlation]({fig_link(fig_paths['correlation'])})

Pearson 看线性相关，Spearman 更看排序是否一致。`KinForm` 的 Spearman 最高，说明在它能覆盖的 563 条子集上，预测排序与实验排序较一致；但它不是全量覆盖方法。全量/近全量方法中，`CataPro`、`KcatNet`、`UniKP-official` 的相关性相对靠前。

### 3. 准确率与覆盖率权衡

![Coverage vs MAE]({fig_link(fig_paths['coverage_vs_mae'])})

这张图的右下角最理想：覆盖率高、误差低。颜色对应上面的分组：蓝色是全量/近全量 sequence+SMILES，橙色是 reaction-aware 子集，绿色是模型特定子集，红色是公开数据重训版，青色是 GO 功能赋值基线。`KcatNet` 位于较理想区域；`CataPro` 覆盖完整且表现稳定；`TurNuP-official`、`PMAK` 误差较低但覆盖率受 reaction SMILES 限制；`KinForm` 相关性好但覆盖条数更少；`GO-HKP` 已覆盖全 benchmark，主要作为非 AI 赋值基线。

### 4. 10 倍误差内比例与系统性偏差

![Within10 and bias]({fig_link(fig_paths['within10_bias'])})

左图是预测落在实验值 10 倍范围内的比例，越高越好；右图是 bias，负值表示整体偏低估，正值表示整体偏高估。`KcatNet` 的 10 倍内比例最高，`PMAK` 和 `TurNuP-official` 也较高。`DLKcat-official`、`UniKP-official`、`PreTKcat`、`DEKP-public-retrained` 整体有不同程度的低估趋势；`GO-HKP` 则明显偏高估。

### 5. 单条记录误差分布

![Error distribution]({fig_link(fig_paths['error_distribution'])})

箱线图展示每个方法在逐条样本上的绝对误差。黑色虚线约等于 2 倍误差，黑色点线是 10 倍误差。`KcatNet` 的中位误差最低；`DEKP-public-retrained` 的中位误差和长尾错误都偏大，说明公开数据重训版还没有达到理想状态。

### 6. 按物种表现

![Species heatmap]({fig_link(fig_paths['species_heatmap'])})

物种分层可以帮助判断模型是否只在某个物种上表现好。总体看，多数方法在 *E. coli* 和 yeast 上有差异；写文章时建议保留 species-level 指标，避免单个 overall 数字掩盖物种偏差。

### 7. 按数据来源表现

![Source heatmap]({fig_link(fig_paths['source_heatmap'])})

BRENDA 和 SABIO-RK 的数据来源、实验条件记录方式不同，分层后能看到模型对不同来源数据的适应性。这个图适合放补充材料，主文可以简述“不同来源之间存在方法表现差异”。

### 8. 预测值 vs 实验值

![Predicted vs true]({fig_link(fig_paths['pred_true'])})

对角线是理想预测，点线是相差 10 倍的范围。这个图能直观看到哪些方法存在压缩动态范围、整体偏高/偏低或极端错误。`KcatNet` 和 `CataPro` 的点云相对更贴近对角线；`DEKP-public-retrained` 的偏离更明显。

## 推荐写作口径

1. 方法可按输入口径和模型性质分组：全量/近全量 sequence+SMILES、reaction-aware 子集、模型特定子集、公开数据重训版、功能相似性 GO 赋值基线。
2. 如果只强调全量覆盖和误差，`KcatNet` 是当前最强基线；如果强调完整 reaction 信息方法，`TurNuP-official` 和 `PMAK` 应在 780 条 reaction-aware 子集内比较。
3. `KinForm` 的相关性最好但覆盖有限，适合描述为“在可覆盖子集上排序能力强”。
4. `DEKP-public-retrained` 应明确是公开数据重训版，不等价于原论文最优官方模型。
5. `GO-HKP` 应明确是非 AI 直接赋值基线；当前结果说明“按 GO 层级直接赋值”在本 benchmark 上不能替代主流 AI 预测方法，但很适合作为一个朴素生物学基线。还需要说明 E. coli 和 yeast 的 GO 来源不同。

## 文件索引

- 总表：`reports/tables/method_eval_summary.csv`
- 注释版总表：`reports/tables/method_eval_summary_annotated.csv`
- 方法分组注释表：`reports/tables/method_group_annotation.csv`
- 当前方法排序表：`reports/tables/method_rank_current_benchmark.csv`
- 物种 MAE 矩阵：`reports/tables/species_mae_matrix.csv`
- 数据来源 MAE 矩阵：`reports/tables/source_database_mae_matrix.csv`
- 图表目录：`reports/figures/kcat_benchmark_summary/`
"""

    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    sns.set_theme(style="whitegrid", font_scale=0.95)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    summary = load_summary()
    rows = load_rows(list(METHOD_META))
    method_annotation = build_method_annotation(summary)
    current_rank = build_rank_tables(summary)

    fig_paths = {}

    path = plot_overall_error(summary)
    savefig(path)
    fig_paths["overall_error"] = path

    path = plot_correlation(summary)
    savefig(path)
    fig_paths["correlation"] = path

    path = plot_coverage_vs_mae(summary)
    savefig(path)
    fig_paths["coverage_vs_mae"] = path

    path = plot_within10_bias(summary)
    savefig(path)
    fig_paths["within10_bias"] = path

    path = plot_error_distribution(rows, summary)
    savefig(path)
    fig_paths["error_distribution"] = path

    path = plot_species_heatmap()
    savefig(path)
    fig_paths["species_heatmap"] = path

    path = plot_source_heatmap()
    savefig(path)
    fig_paths["source_heatmap"] = path

    path = plot_scatter_selected(rows, summary)
    savefig(path)
    fig_paths["pred_true"] = path

    write_report(summary, current_rank, method_annotation, fig_paths)
    print(f"Wrote report: {REPORT_PATH}")
    print(f"Wrote figures: {FIGURE_DIR}")


if __name__ == "__main__":
    main()
