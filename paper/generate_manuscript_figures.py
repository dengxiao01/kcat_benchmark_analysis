#!/usr/bin/env python3
"""Generate manuscript Figures 1-4 from versioned benchmark tables."""

from __future__ import annotations

import json
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BASE = Path(__file__).resolve().parent.parent
TABLE_DIR = BASE / "paper" / "tables_v1.2.0"
OUTPUT_DIR = BASE / "paper" / "figures"
RELEASE_FILE = BASE / "configs" / "benchmark_release.json"

METHOD_ORDER = [
    "KcatNet", "CataPro", "PreTKcat", "UniKP", "SELFprot", "DLKcat",
    "TurNuP", "PMAK", "KinForm-L", "CatPred", "DEKP-public-retrained", "GO-HKP",
]
REGIME_BY_METHOD = {
    "KcatNet": "Released sequence+substrate checkpoint",
    "CataPro": "Released sequence+substrate checkpoint",
    "PreTKcat": "Temperature-conditioned public retraining",
    "UniKP": "Released sequence+substrate checkpoint",
    "SELFprot": "Released sequence+substrate checkpoint",
    "DLKcat": "Released sequence+substrate checkpoint",
    "TurNuP": "Reaction-aware checkpoint",
    "PMAK": "Reaction-aware checkpoint",
    "KinForm-L": "Method-specific checkpoint subset",
    "CatPred": "Method-specific checkpoint subset",
    "DEKP-public-retrained": "Structure-aware public retraining",
    "GO-HKP": "Functional-assignment baseline",
}
REGIME_COLORS = {
    "Released sequence+substrate checkpoint": "#326FA8",
    "Temperature-conditioned public retraining": "#A24B52",
    "Reaction-aware checkpoint": "#20938F",
    "Method-specific checkpoint subset": "#D17A16",
    "Structure-aware public retraining": "#8157A3",
    "Functional-assignment baseline": "#747474",
}

METHOD_COLORS = {
    method: REGIME_COLORS[REGIME_BY_METHOD[method]] for method in METHOD_ORDER
}
METHOD_DISPLAY = {
    method: "DEKP" if method == "DEKP-public-retrained" else method.removesuffix("-official")
    for method in METHOD_ORDER
}


def method_display(method: str) -> str:
    return METHOD_DISPLAY.get(method, method)


def read_table(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLE_DIR / f"{name}.csv")


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 10.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def panel_label(axis: plt.Axes, label: str, x: float = -0.10, y: float = 1.05) -> None:
    axis.text(x, y, label, transform=axis.transAxes, fontsize=14,
              fontweight="bold", va="top", ha="left")


def add_release_footer(fig: plt.Figure, release: dict) -> None:
    fig.text(
        0.995, 0.004,
        f"Benchmark v{release['benchmark_version']} | data freeze {release['data_freeze_date']}",
        ha="right", va="bottom", fontsize=7.5, color="#5C6470",
    )


def save_figure(fig: plt.Figure, number: int, release: dict) -> None:
    layout_engine = fig.get_layout_engine()
    if layout_engine is not None:
        layout_engine.set(rect=(0.0, 0.025, 1.0, 1.0))
    add_release_footer(fig, release)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUTPUT_DIR / f"Figure{number}"
    description = (
        f"Benchmark version {release['benchmark_version']}; "
        f"data freeze {release['data_freeze_date']}"
    )
    fig.savefig(
        stem.with_suffix(".png"), dpi=300, bbox_inches="tight",
        metadata={
            "Title": f"kcat benchmark Figure {number}",
            "Author": "kcat_benchmark_analysis",
            "Description": description,
        },
    )
    fig.savefig(
        stem.with_suffix(".pdf"), bbox_inches="tight",
        metadata={
            "Title": f"kcat benchmark Figure {number}",
            "Author": "kcat_benchmark_analysis",
            "Subject": description,
        },
    )
    plt.close(fig)


def workflow_box(axis: plt.Axes, xy: tuple[float, float], width: float, height: float,
                 title: str, body: str, facecolor: str) -> None:
    x, y = xy
    axis.add_patch(FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.014,rounding_size=0.015",
        linewidth=1.2, edgecolor="#536173", facecolor=facecolor,
        transform=axis.transAxes, clip_on=False,
    ))
    axis.text(x + width / 2, y + height * 0.66, title,
              transform=axis.transAxes, ha="center", va="center",
              fontsize=11.5, fontweight="bold")
    axis.text(x + width / 2, y + height * 0.31, body,
              transform=axis.transAxes, ha="center", va="center",
              fontsize=9.2, linespacing=1.15)


def workflow_arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(FancyArrowPatch(
        start, end, transform=axis.transAxes, arrowstyle="-|>",
        mutation_scale=13, linewidth=1.2, color="#68778A",
        shrinkA=0, shrinkB=0,
    ))


def figure1(release: dict) -> None:
    truth = pd.read_csv(BASE / "data" / "final" / "benchmark_ready_catpred.csv")
    source = read_table("S4_Source_support")
    matching = read_table("S3_Matching_levels")
    conditions = read_table("S5_Condition_metadata")
    funnel = pd.read_csv(BASE / "reports" / "tables" / "benchmark_build_funnel.csv")

    fig = plt.figure(figsize=(14.0, 13.2), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1.25, 1.0, 0.9])
    ax_a = fig.add_subplot(grid[0, :])
    ax_a.set_axis_off()
    panel_label(ax_a, "a", x=0.00, y=1.02)
    ax_a.text(0.08, 1.0, "Benchmark construction and unified evaluation workflow",
              transform=ax_a.transAxes, va="top", fontsize=13, fontweight="bold")
    width, height = 0.255, 0.31
    x_positions = [0.02, 0.3725, 0.725]
    y_top, y_bottom = 0.58, 0.08
    candidate_ecoli = int(funnel.loc[funnel["species"].eq("ecoli"),
                                      "enzyme_substrate_entries"].iloc[0])
    candidate_yeast = int(funnel.loc[funnel["species"].eq("yeast"),
                                      "enzyme_substrate_entries"].iloc[0])
    total_n = len(truth)
    ecoli_n = int(truth["species"].eq("ecoli").sum())
    yeast_n = int(truth["species"].eq("yeast").sum())
    supported_n = int(truth["experimental_substrate_support"].eq("substrate_supported").sum())
    workflow_box(ax_a, (x_positions[0], y_top), width, height,
                 "1  Genome-scale models", "eciML1515 and yeast-GEM", "#EDF3F7")
    workflow_box(ax_a, (x_positions[1], y_top), width, height,
                 "2  Candidate generation",
                 f"{candidate_ecoli:,} E. coli pairs\n{candidate_yeast:,} yeast pairs", "#EDF3F7")
    workflow_box(ax_a, (x_positions[2], y_top), width, height,
                 "3  Molecular mapping", "UniProt sequence\nstandardized reactant SMILES", "#EDF3F7")
    workflow_box(ax_a, (x_positions[2], y_bottom), width, height,
                 "4  Experimental matching", "BRENDA and SABIO-RK\nhierarchical identifiers", "#EDF3F7")
    workflow_box(ax_a, (x_positions[1], y_bottom), width, height,
                 "5  Final benchmark",
                 f"{total_n:,} complete-resource records\n{supported_n:,} substrate-supported\n"
                 f"{ecoli_n:,} E. coli; {yeast_n:,} yeast", "#EAF5EE")
    workflow_box(ax_a, (x_positions[0], y_bottom), width, height,
                 "6  Unified evaluation", "12 methods\n6 inference regimes", "#F8F0E6")
    workflow_arrow(ax_a, (0.275, 0.735), (0.3725, 0.735))
    workflow_arrow(ax_a, (0.6275, 0.735), (0.725, 0.735))
    workflow_arrow(ax_a, (0.8525, 0.58), (0.8525, 0.39))
    workflow_arrow(ax_a, (0.725, 0.235), (0.6275, 0.235))
    workflow_arrow(ax_a, (0.3725, 0.235), (0.275, 0.235))

    species_order = ["ecoli", "yeast"]
    species_labels = ["E. coli", "S. cerevisiae"]
    ax_b = fig.add_subplot(grid[1, 0])
    panel_label(ax_b, "b")
    ax_b.set_title("Experimental source support", loc="left")
    source_order = ["BRENDA", "SABIO-RK", "BRENDA;SABIO-RK"]
    source_labels = ["BRENDA", "SABIO-RK", "Both sources"]
    source_colors = ["#4A78C2", "#73AB43", "#F07C2B"]
    bottom = np.zeros(2)
    for database, label, color in zip(source_order, source_labels, source_colors):
        values = [int(source.loc[
            source["Species"].eq(species) & source["Source database"].eq(database), "Records"
        ].sum()) for species in species_order]
        ax_b.bar(species_labels, values, bottom=bottom, color=color, width=0.56, label=label)
        bottom += np.asarray(values)
    for index, value in enumerate(bottom):
        ax_b.text(index, value + 8, f"{int(value)}", ha="center", va="bottom", fontsize=9)
    ax_b.set_ylabel("Benchmark records")
    ax_b.set_ylim(0, max(bottom) * 1.16)
    ax_b.legend(frameon=False, loc="upper center", ncol=3)
    ax_b.grid(axis="y", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_b.set_axisbelow(True)

    ax_c = fig.add_subplot(grid[1, 1])
    panel_label(ax_c, "c")
    ax_c.set_title("Hierarchical matching support", loc="left")
    match_order = [
        "Species + EC + UniProt + substrate ID",
        "Species + EC + substrate ID",
        "Species + EC + substrate name",
    ]
    match_colors = ["#3378A5", "#85B1C9", "#D5DCE1"]
    bottom = np.zeros(2)
    for level, color in zip(match_order, match_colors):
        values = [int(matching.loc[
            matching["Species"].eq(species) & matching["Matching level"].eq(level), "Records"
        ].sum()) for species in species_order]
        ax_c.bar(species_labels, values, bottom=bottom, color=color, width=0.56,
                 label=level.replace("Species + EC + ", ""))
        bottom += np.asarray(values)
    ax_c.set_ylabel("Benchmark records")
    ax_c.set_ylim(0, max(bottom) * 1.10)
    ax_c.legend(frameon=False, loc="upper right")
    ax_c.grid(axis="y", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_c.set_axisbelow(True)

    ax_d = fig.add_subplot(grid[2, 0])
    panel_label(ax_d, "d")
    ax_d.set_title("Experimental-condition metadata", loc="left")
    x = np.arange(2)
    width_bar = 0.32
    ph = [float(conditions.loc[conditions["Species"].eq(species),
                               "pH available (%)"].iloc[0]) for species in species_order]
    temperature = [float(conditions.loc[conditions["Species"].eq(species),
                                        "Temperature available (%)"].iloc[0])
                   for species in species_order]
    bars1 = ax_d.bar(x - width_bar / 2, ph, width_bar, color="#5781E6", label="pH")
    bars2 = ax_d.bar(x + width_bar / 2, temperature, width_bar,
                     color="#5FCBA0", label="Temperature")
    for bars in (bars1, bars2):
        for bar in bars:
            ax_d.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.2,
                      f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=8.8)
    ax_d.set_xticks(x, species_labels)
    ax_d.set_ylabel("Metadata available (%)")
    ax_d.set_ylim(0, 105)
    ax_d.legend(frameon=False, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.22))
    ax_d.grid(axis="y", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_d.set_axisbelow(True)

    ax_e = fig.add_subplot(grid[2, 1])
    panel_label(ax_e, "e")
    ax_e.set_title("Benchmark scale and biochemical breadth", loc="left")
    ax_e.set_axis_off()
    stats = [
        (f"{total_n:,}", "sequence-associated records"),
        (f"{truth['reaction_id'].nunique():,}", "distinct reactions"),
        (f"{truth['uniprot_id'].nunique():,}", "unique proteins"),
        (f"{truth['ec_number'].nunique():,}", "EC annotations"),
        (">9 orders", "experimental kcat range"),
        (f"{truth['true_kcat_log10'].median():.3f}", "median log10(kcat)"),
    ]
    positions = [(0.18, 0.78), (0.64, 0.78), (0.18, 0.43),
                 (0.64, 0.43), (0.18, 0.08), (0.64, 0.08)]
    for (value, label), (x_value, y_value) in zip(stats, positions):
        ax_e.text(x_value, y_value + 0.11, value, transform=ax_e.transAxes,
                  fontsize=20, fontweight="bold", color="#245A82", ha="center")
        ax_e.text(x_value, y_value, label, transform=ax_e.transAxes,
                  fontsize=9.5, color="#626B75", ha="center")
    save_figure(fig, 1, release)


def metric_ci(table: pd.DataFrame, metric: str) -> pd.DataFrame:
    return table.loc[table["metric"].eq(metric)].set_index("method").loc[METHOD_ORDER].reset_index()


def figure2(release: dict) -> None:
    full = read_table("S1_Full_metrics").set_index("Method").loc[METHOD_ORDER].reset_index()
    cluster = read_table("S18_Cluster_bootstrap")
    pair_ci = cluster.loc[
        cluster["analysis_scope"].eq("achieved_evaluation_set")
        & cluster["cluster_type"].eq("pair")
    ].set_index("method").loc[METHOD_ORDER].reset_index()

    fig, axes = plt.subplots(2, 2, figsize=(15.2, 12.2), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    panel_label(ax_a, "a")
    ax_a.set_title("Achieved coverage and descriptive MAE", loc="left")
    offsets = {
        "KcatNet": (5, -2), "CataPro": (5, -8), "PreTKcat": (5, 3),
        "UniKP": (5, -8), "SELFprot": (-42, 7), "DLKcat": (-42, -8),
        "GO-HKP": (-46, -8), "DEKP-public-retrained": (-34, -8),
        "TurNuP": (-45, 8), "PMAK": (7, -8), "KinForm-L": (7, 4), "CatPred": (7, 4),
    }
    for _, row in full.iterrows():
        method = row["Method"]
        coverage, mae = float(row["Coverage (%)"]), float(row["MAE log10"])
        ax_a.scatter(coverage, mae, s=78, color=METHOD_COLORS[method],
                     edgecolor="white", linewidth=0.8, zorder=3)
        ax_a.annotate(method_display(method), (coverage, mae), xytext=offsets[method],
                      textcoords="offset points", fontsize=8.6)
    ax_a.set_xlim(53, 103)
    ax_a.set_ylim(1.13, 0.64)
    ax_a.set_xlabel("Benchmark coverage (%)")
    ax_a.set_ylabel("MAE in log10(kcat), lower is better")
    ax_a.grid(color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_a.set_axisbelow(True)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=color,
               markeredgecolor="white", markersize=7, label=regime)
        for regime, color in REGIME_COLORS.items()
    ]
    ax_a.legend(handles=handles, frameon=False, ncol=1, loc="lower left", fontsize=7.8)

    panel_label(ax_b, "b")
    ax_b.set_title("Applicable benchmark records", loc="left")
    coverage_order = full.sort_values(["Coverage (%)", "MAE log10"], ascending=[True, False])
    y = np.arange(len(coverage_order))
    ax_b.barh(y, coverage_order["Coverage (%)"],
              color=[METHOD_COLORS[m] for m in coverage_order["Method"]], height=0.66)
    ax_b.set_yticks(y, [method_display(method) for method in coverage_order["Method"]])
    ax_b.set_xlim(0, 108)
    ax_b.set_xlabel(f"Coverage of the {int(release['canonical_rows']):,}-record benchmark (%)")
    for index, (_, row) in enumerate(coverage_order.iterrows()):
        ax_b.text(float(row["Coverage (%)"]) + 1.0, index, f"{int(row['n'])}",
                  va="center", fontsize=8.5)
    ax_b.grid(axis="x", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_b.set_axisbelow(True)

    panel_label(ax_c, "c")
    ax_c.set_title("Two-fold and ten-fold agreement", loc="left")
    accuracy_order = full.sort_values("Within 10-fold (%)", ascending=True)
    y = np.arange(len(accuracy_order))
    two = accuracy_order["Within 2-fold (%)"].to_numpy(float)
    ten = accuracy_order["Within 10-fold (%)"].to_numpy(float)
    ax_c.hlines(y, two, ten, color="#B9C2CC", linewidth=2)
    ax_c.scatter(two, y, color="#326FA8", s=45, label="Within two-fold", zorder=3)
    ax_c.scatter(ten, y, color="#20938F", s=45, label="Within ten-fold", zorder=3)
    ax_c.set_yticks(y, [method_display(method) for method in accuracy_order["Method"]])
    ax_c.set_xlim(0, 86)
    ax_c.set_xlabel("Predictions within the experimental value (%)")
    ax_c.legend(frameon=False, loc="lower right")
    ax_c.grid(axis="x", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_c.set_axisbelow(True)

    panel_label(ax_d, "d")
    ax_d.set_title("Pair-cluster uncertainty on achieved sets", loc="left")
    pair_ci = pair_ci.sort_values("row_weighted_mae_log10", ascending=False)
    y = np.arange(len(pair_ci))
    estimate = pair_ci["row_weighted_mae_log10"].to_numpy(float)
    lower = estimate - pair_ci["cluster_bootstrap_ci_low_95"].to_numpy(float)
    upper = pair_ci["cluster_bootstrap_ci_high_95"].to_numpy(float) - estimate
    ax_d.errorbar(estimate, y, xerr=np.vstack([lower, upper]), fmt="none",
                  ecolor="#8A98A8", elinewidth=1.3, capsize=2)
    ax_d.scatter(estimate, y,
                 c=[METHOD_COLORS[m] for m in pair_ci["method"]], s=42, zorder=3)
    ax_d.set_yticks(y, [method_display(method) for method in pair_ci["method"]])
    ax_d.set_xlabel("MAE in log10(kcat), 95% pair-cluster CI")
    ax_d.grid(axis="x", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_d.set_axisbelow(True)
    save_figure(fig, 2, release)


def build_scope_heatmap() -> tuple[pd.DataFrame, pd.DataFrame]:
    full = read_table("S1_Full_metrics").set_index("Method")
    reaction = read_table("S6_Reaction_subset").set_index("method")
    matched = read_table("S7_Available_case")
    catpred_scope = next(
        value for value in matched["subset"].unique() if str(value).startswith("catpred_accessible_scope_")
    )
    kinform_scope = next(
        value for value in matched["subset"].unique() if str(value).startswith("kinform_accessible_scope_")
    )
    rows, counts = [], []
    for method in METHOD_ORDER:
        catpred = matched.loc[matched["subset"].eq(catpred_scope)
                              & matched["method"].eq(method)].iloc[0]
        kinform = matched.loc[matched["subset"].eq(kinform_scope)
                              & matched["method"].eq(method)].iloc[0]
        rows.append({
            "method": method,
            "Achieved\npredictions": float(full.loc[method, "MAE log10"]),
            "Reaction scope\navailable case": float(reaction.loc[method, "mae_log10_common_subset"]),
            "CatPred scope\navailable case": float(catpred["mae_log10"]),
            "KinForm-L scope\navailable case": float(kinform["mae_log10"]),
        })
        counts.append({
            "method": method,
            "Achieved\npredictions": int(full.loc[method, "n"]),
            "Reaction scope\navailable case": int(reaction.loc[method, "n_common_subset"]),
            "CatPred scope\navailable case": int(catpred["n"]),
            "KinForm-L scope\navailable case": int(kinform["n"]),
        })
    return pd.DataFrame(rows).set_index("method"), pd.DataFrame(counts).set_index("method")


def conventional_one_decimal(value: float) -> str:
    return str(Decimal(str(value)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))


def figure3(release: dict) -> None:
    heatmap, counts = build_scope_heatmap()
    paired = read_table("S13_Paired_bootstrap")
    ranks = read_table("S12_Rank_stability")
    reaction_scope = next(
        value for value in ranks["comparison_set"].unique() if str(value).startswith("reaction_aware_common_")
    )
    reaction_n = int(read_table("S6_Reaction_subset")["common_subset_total_rows"].iloc[0])
    ranks = ranks.loc[ranks["comparison_set"].eq(reaction_scope)]
    fig = plt.figure(figsize=(14.8, 9.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0])

    ax_a = fig.add_subplot(grid[:, 0])
    panel_label(ax_a, "a")
    ax_a.set_title("Available-case MAE within defined scopes", loc="left", pad=12)
    matrix = heatmap.to_numpy(float)
    image = ax_a.imshow(matrix, cmap="YlGnBu_r", aspect="auto", vmin=0.64, vmax=1.10)
    ax_a.set_xticks(np.arange(heatmap.shape[1]), heatmap.columns)
    ax_a.set_yticks(
        np.arange(heatmap.shape[0]),
        [method_display(method) for method in heatmap.index],
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            color = "white" if value < 0.77 else "#17212B"
            ax_a.text(j, i, f"{value:.3f}\n(n={int(counts.iloc[i, j])})",
                      ha="center", va="center", fontsize=8.5, color=color)
    colorbar = fig.colorbar(image, ax=ax_a, fraction=0.045, pad=0.03)
    colorbar.set_label("MAE in log10(kcat)")
    ax_a.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax_a.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax_a.grid(which="minor", color="white", linewidth=0.8, alpha=0.6)
    ax_a.tick_params(which="minor", bottom=False, left=False)

    ax_b = fig.add_subplot(grid[0, 1])
    panel_label(ax_b, "b")
    ax_b.set_title(f"Strictly paired reaction-aware subset (n = {reaction_n:,})", loc="left")
    labels = [f"{row.method_a} - {row.method_b}" for row in paired.itertuples()]
    y = np.arange(len(paired))[::-1]
    differences = paired["mae_difference_a_minus_b"].to_numpy(float)
    low = paired["bootstrap_ci_low_95"].to_numpy(float)
    high = paired["bootstrap_ci_high_95"].to_numpy(float)
    ax_b.errorbar(differences, y, xerr=np.vstack([differences - low, high - differences]),
                  fmt="o", color="#326FA8", ecolor="#326FA8",
                  elinewidth=1.4, capsize=3, markersize=5)
    ax_b.axvline(0, color="#526171", linewidth=0.9)
    ax_b.set_yticks(y, labels)
    ax_b.set_xlabel("Paired MAE difference (method A - method B)")
    ax_b.grid(axis="x", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_b.set_axisbelow(True)

    ax_c = fig.add_subplot(grid[1, 1])
    panel_label(ax_c, "c")
    ax_c.set_title("Rank-one frequency in paired bootstrap", loc="left")
    rank_order = ["PMAK", "TurNuP", "KcatNet"]
    ranks = ranks.set_index("method").loc[rank_order].reset_index()
    values = 100.0 * ranks["rank_1_bootstrap_frequency"].to_numpy(float)
    y = np.arange(len(ranks))
    bars = ax_c.barh(y, values, color=[METHOD_COLORS[m] for m in ranks["method"]], height=0.6)
    ax_c.set_yticks(y, ranks["method"])
    ax_c.set_xlim(0, 70)
    ax_c.set_xlabel("Bootstrap replicates ranked first (%)")
    for bar, value in zip(bars, values):
        ax_c.text(value + 0.8, bar.get_y() + bar.get_height() / 2,
                  conventional_one_decimal(value), va="center", fontsize=9)
    ax_c.grid(axis="x", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_c.set_axisbelow(True)
    save_figure(fig, 3, release)


def subgroup_value(table: pd.DataFrame, method: str, feature: str,
                   group: str | bool, column: str) -> float:
    subset = table.loc[table["method"].eq(method) & table["feature"].eq(feature)].copy()
    if isinstance(group, bool):
        group_values = subset["group"].astype(str).str.lower().map({"true": True, "false": False})
        subset = subset.loc[group_values.eq(group)]
    else:
        subset = subset.loc[subset["group"].astype(str).eq(group)]
    return float(subset[column].iloc[0])


def figure4(release: dict) -> None:
    stratified = read_table("S9_Error_stratification")
    methods = ["KcatNet", "CataPro", "TurNuP", "PMAK", "KinForm-L", "CatPred"]
    categories = [
        ("species", "ecoli", "E. coli"),
        ("species", "yeast", "S. cerevisiae"),
        ("experimental_substrate_support", "substrate_supported", "Substrate\nsupported"),
        ("experimental_substrate_support", "participant_ambiguous", "Participant\nambiguous"),
        ("substrate_role_group_substrate_supported", "other_reactant", "Supported\nother reactant"),
        ("substrate_role_group_substrate_supported", "currency_or_cofactor", "Supported\ncurrency/cofactor"),
        ("substrate_role_group_substrate_supported", "carrier_linked_variable", "Supported\ncarrier-linked"),
    ]
    matrix = np.asarray([[
        subgroup_value(stratified, method, feature, group, "mean_abs_error_log10")
        for feature, group, _ in categories
    ] for method in methods])
    count_matrix = np.asarray([[
        subgroup_value(stratified, method, feature, group, "n")
        for feature, group, _ in categories
    ] for method in methods], dtype=int)
    fig = plt.figure(figsize=(14.8, 10.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.18, 1.0])

    ax_a = fig.add_subplot(grid[0, :])
    panel_label(ax_a, "a", x=-0.07)
    ax_a.set_title("Mean absolute error by species, substrate evidence and registry role", loc="left")
    finite = matrix[np.isfinite(matrix)]
    image = ax_a.imshow(
        matrix,
        cmap="YlOrRd",
        aspect="auto",
        vmin=max(0.0, float(finite.min()) - 0.03),
        vmax=float(finite.max()) + 0.03,
    )
    ax_a.set_yticks(np.arange(len(methods)), methods)
    ax_a.set_xticks(np.arange(len(categories)), [item[2] for item in categories])
    ax_a.tick_params(axis="x", rotation=27)
    for label in ax_a.get_xticklabels():
        label.set_ha("right")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            midpoint = (float(finite.min()) + float(finite.max())) / 2
            color = "white" if matrix[i, j] > midpoint else "#17212B"
            ax_a.text(
                j, i, f"{matrix[i, j]:.2f}\n(n={count_matrix[i, j]})",
                ha="center", va="center", fontsize=7.8, color=color, linespacing=1.05,
            )
    colorbar = fig.colorbar(image, ax=ax_a, fraction=0.025, pad=0.02)
    colorbar.set_label("MAE in log10(kcat)")
    ax_a.text(1.0, -0.30,
              "Role columns use only substrate-supported records and the joint name/identifier/standardized-structure registry; carrier-linked metabolites are separate.",
              transform=ax_a.transAxes, ha="right", va="top", fontsize=8.2, color="#5C6470")

    ax_b = fig.add_subplot(grid[1, 0])
    panel_label(ax_b, "b")
    ax_b.set_title("Subgroup MAE contrasts", loc="left")
    y = np.arange(len(methods))[::-1]
    species_contrast, role_contrast = [], []
    for method in methods:
        species_contrast.append(
            subgroup_value(stratified, method, "species", "ecoli", "mean_abs_error_log10")
            - subgroup_value(stratified, method, "species", "yeast", "mean_abs_error_log10")
        )
        role_contrast.append(
            subgroup_value(stratified, method, "substrate_role_group_substrate_supported", "other_reactant",
                           "mean_abs_error_log10")
            - subgroup_value(stratified, method, "substrate_role_group_substrate_supported", "currency_or_cofactor",
                             "mean_abs_error_log10")
        )
    for row_y, species_value, role_value in zip(y, species_contrast, role_contrast):
        ax_b.hlines(
            row_y,
            min(species_value, role_value),
            max(species_value, role_value),
            color="#B9C2CC",
            linewidth=1.8,
            zorder=1,
        )
    ax_b.scatter(species_contrast, y, color="#326FA8", s=48,
                 label="E. coli - yeast", zorder=3)
    ax_b.scatter(role_contrast, y, color="#E67E22", marker="D", s=42,
                 label="Supported other - currency/cofactor", zorder=3)
    ax_b.axvline(0, color="#526171", linewidth=0.9)
    ax_b.set_yticks(y, methods)
    ax_b.set_xlabel("Difference in MAE (log10 units)")
    ax_b.set_ylim(-0.5, len(methods) + 0.8)
    ax_b.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        ncol=2,
        fontsize=8.4,
        columnspacing=0.9,
        handletextpad=0.4,
    )
    ax_b.grid(axis="x", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_b.set_axisbelow(True)

    ax_c = fig.add_subplot(grid[1, 1])
    panel_label(ax_c, "c")
    ax_c.set_title("Large-error frequency by species", loc="left")
    x = np.arange(len(methods))
    width = 0.36
    ecoli = [100.0 * subgroup_value(stratified, method, "species", "ecoli", "outlier_fraction")
             for method in methods]
    yeast = [100.0 * subgroup_value(stratified, method, "species", "yeast", "outlier_fraction")
             for method in methods]
    ax_c.bar(x - width / 2, ecoli, width, color="#4C78A8", label="E. coli")
    ax_c.bar(x + width / 2, yeast, width, color="#F58518", label="S. cerevisiae")
    ax_c.set_xticks(x, methods, rotation=25, ha="right")
    ax_c.set_ylabel("Rows with |error| > 1 (%)")
    ax_c.set_ylim(0, 42)
    ax_c.legend(frameon=False, loc="upper center", ncol=2)
    ax_c.grid(axis="y", color="#D9DEE5", linewidth=0.6, alpha=0.8)
    ax_c.set_axisbelow(True)
    save_figure(fig, 4, release)


def main() -> None:
    configure_style()
    release = json.loads(RELEASE_FILE.read_text(encoding="utf-8"))
    figure1(release)
    figure2(release)
    figure3(release)
    figure4(release)
    print(f"Wrote manuscript figures to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
