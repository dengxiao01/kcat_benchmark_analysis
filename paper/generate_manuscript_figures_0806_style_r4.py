#!/usr/bin/env python3
"""Generate v1.2.0-r4 Figures 1-4 in the visual language of the 0806 draft."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle
from matplotlib.transforms import Bbox


BASE = Path(__file__).resolve().parent.parent
TABLE_DIR = BASE / "paper" / "tables_v1.2.0"
OUTPUT_DIR = BASE / "paper" / "0806" / "polished_v1.2.0-r4" / "figures"
PANEL_DIR = OUTPUT_DIR / "panels"
SOURCE_DIR = OUTPUT_DIR / "source_data"
PACKAGE_DIR = OUTPUT_DIR.parent / "figure_packages"
TRUTH_FILE = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
FUNNEL_FILE = BASE / "reports" / "tables" / "benchmark_build_funnel.csv"
TABLE0_FILE = TABLE_DIR / "Table0.csv"
AUDIT_FILE = TABLE_DIR / "Record_audit.csv"


METHOD_ORDER = [
    "KcatNet", "CataPro", "UniKP", "SELFprot", "DLKcat", "PreTKcat",
    "TurNuP", "PMAK", "KinForm-L", "CatPred", "DEKP-public-retrained", "GO-HKP",
]
REGIME_ORDER = [
    "Released sequence+substrate checkpoint",
    "Temperature-conditioned public retraining",
    "Reaction-aware checkpoint",
    "Method-specific checkpoint subset",
    "Structure-aware public retraining",
    "Functional-assignment baseline",
]
REGIME_COLORS = {
    "Released sequence+substrate checkpoint": "#1F5CC4",
    "Temperature-conditioned public retraining": "#B24C54",
    "Reaction-aware checkpoint": "#159D91",
    "Method-specific checkpoint subset": "#F57C00",
    "Structure-aware public retraining": "#7D3FA0",
    "Functional-assignment baseline": "#7A7A7A",
}
REGIME_SHORT = {
    "Released sequence+substrate checkpoint": "Released sequence+SMILES",
    "Temperature-conditioned public retraining": "Temperature-conditioned retraining",
    "Reaction-aware checkpoint": "Reaction-aware",
    "Method-specific checkpoint subset": "Method-specific subset",
    "Structure-aware public retraining": "Structure-aware retraining",
    "Functional-assignment baseline": "Functional-assignment baseline",
}
METHOD_DISPLAY = {
    method: "DEKP" if method == "DEKP-public-retrained" else method
    for method in METHOD_ORDER
}
METHOD_COLORS = {
    method: REGIME_COLORS[
        {
            "KcatNet": REGIME_ORDER[0], "CataPro": REGIME_ORDER[0],
            "UniKP": REGIME_ORDER[0], "SELFprot": REGIME_ORDER[0],
            "DLKcat": REGIME_ORDER[0], "PreTKcat": REGIME_ORDER[1],
            "TurNuP": REGIME_ORDER[2], "PMAK": REGIME_ORDER[2],
            "KinForm-L": REGIME_ORDER[3], "CatPred": REGIME_ORDER[3],
            "DEKP-public-retrained": REGIME_ORDER[4], "GO-HKP": REGIME_ORDER[5],
        }[method]
    ]
    for method in METHOD_ORDER
}


def read_table(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLE_DIR / f"{name}.csv")


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.7,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.3,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.2,
        "legend.fontsize": 7.5,
        "legend.title_fontsize": 8.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def display_method(method: str) -> str:
    return METHOD_DISPLAY.get(method, method)


def panel_label(ax: plt.Axes, label: str, x: float = -0.14, y: float = 1.08) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=15, fontweight="bold",
            ha="left", va="top")


def style_axes(ax: plt.Axes, grid_axis: str = "both") -> None:
    ax.grid(True, axis=grid_axis, linestyle=(0, (3, 3)), linewidth=0.55,
            color="#D7D7D7", zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def regime_legend(ax: plt.Axes, loc: str = "best", bbox_to_anchor=None, ncol: int = 1):
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.5,
               markerfacecolor=REGIME_COLORS[regime], markeredgecolor="white",
               label=REGIME_SHORT[regime])
        for regime in REGIME_ORDER
    ]
    return ax.legend(handles=handles, title="Inference regime", frameon=True,
                     framealpha=1, edgecolor="#C8C8C8", loc=loc,
                     bbox_to_anchor=bbox_to_anchor, ncol=ncol)


def group_separators(ax: plt.Axes, y: np.ndarray) -> None:
    for index in [4, 5, 7, 9, 10]:
        if index + 1 < len(y):
            ax.axhline((y[index] + y[index + 1]) / 2, color="#BDBDBD",
                       linewidth=0.7, linestyle=(0, (4, 3)), zorder=0)


def save_composite(fig: plt.Figure, number: int) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.08}
        if suffix == "png":
            kwargs["dpi"] = 350
        fig.savefig(OUTPUT_DIR / f"Figure{number}.{suffix}", **kwargs)
    plt.close(fig)


def save_axis_panel(fig: plt.Figure, ax: plt.Axes, stem: str) -> None:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    bbox = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    bbox = bbox.expanded(1.04, 1.08)
    fig.savefig(PANEL_DIR / f"{stem}.png", dpi=350, bbox_inches=bbox, pad_inches=0.04)
    fig.savefig(PANEL_DIR / f"{stem}.pdf", bbox_inches=bbox, pad_inches=0.04)


def save_axes_panel(
    fig: plt.Figure,
    axes: list[plt.Axes],
    stem: str,
    extra_artists: list[object] | None = None,
) -> None:
    """Save one logical panel assembled from several Matplotlib axes."""
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bboxes = [axis.get_tightbbox(renderer) for axis in axes]
    for artist in extra_artists or []:
        bboxes.append(artist.get_window_extent(renderer))
    bbox = Bbox.union(bboxes).transformed(fig.dpi_scale_trans.inverted())
    bbox = bbox.expanded(1.025, 1.08)
    fig.savefig(PANEL_DIR / f"{stem}.png", dpi=350, bbox_inches=bbox, pad_inches=0.04)
    fig.savefig(PANEL_DIR / f"{stem}.pdf", bbox_inches=bbox, pad_inches=0.04)


def draw_network_icon(ax: plt.Axes, cx: float, cy: float, color: str) -> None:
    points = np.array([[-0.035, 0], [-0.012, 0.034], [0.022, 0.025], [0.040, -0.012], [0.005, -0.035]])
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (1, 4), (2, 4)]:
        ax.plot([cx + points[i, 0], cx + points[j, 0]],
                [cy + points[i, 1], cy + points[j, 1]], color=color, lw=0.8,
                transform=ax.transAxes, clip_on=False)
    for dx, dy in points:
        ax.add_patch(Circle((cx + dx, cy + dy), 0.006, transform=ax.transAxes,
                            facecolor="white", edgecolor=color, lw=0.8))


def draw_document_icon(ax: plt.Axes, cx: float, cy: float, color: str) -> None:
    ax.add_patch(Rectangle((cx - 0.035, cy - 0.04), 0.048, 0.072,
                           transform=ax.transAxes, facecolor="white", edgecolor=color, lw=0.9))
    for offset in [0.018, 0.0, -0.018]:
        ax.plot([cx - 0.026, cx + 0.004], [cy + offset, cy + offset],
                transform=ax.transAxes, color=color, lw=0.8)
    for offset in [0.025, 0.0, -0.025]:
        ax.add_patch(Circle((cx + 0.036, cy + offset), 0.007, transform=ax.transAxes,
                            facecolor="#5FC4B8", edgecolor=color, lw=0.6))


def draw_mapping_icon(ax: plt.Axes, cx: float, cy: float, color: str) -> None:
    for offset in [-0.025, -0.010, 0.005, 0.020]:
        ax.plot([cx - 0.045, cx - 0.005], [cy + offset, cy + offset],
                transform=ax.transAxes, color=color, lw=1.0)
    theta = np.linspace(0, 2 * np.pi, 7)
    ax.plot(cx + 0.025 + 0.025 * np.cos(theta), cy + 0.022 * np.sin(theta),
            transform=ax.transAxes, color=color, lw=1.0)


def draw_database_icon(ax: plt.Axes, cx: float, cy: float, color: str) -> None:
    for offset in [-0.025, 0, 0.025]:
        ax.add_patch(Rectangle((cx - 0.035, cy + offset - 0.014), 0.07, 0.028,
                               transform=ax.transAxes, facecolor="#EEE7F8",
                               edgecolor=color, lw=0.8))


def draw_filter_icon(ax: plt.Axes, cx: float, cy: float, color: str) -> None:
    poly = Polygon([[cx - 0.045, cy + 0.035], [cx + 0.045, cy + 0.035],
                    [cx + 0.012, cy - 0.005], [cx + 0.012, cy - 0.045],
                    [cx - 0.012, cy - 0.035], [cx - 0.012, cy - 0.005]],
                   closed=True, transform=ax.transAxes, facecolor="#F4E5B8",
                   edgecolor=color, lw=0.9)
    ax.add_patch(poly)


def draw_chart_icon(ax: plt.Axes, cx: float, cy: float, color: str) -> None:
    for index, height in enumerate([0.025, 0.045, 0.065]):
        ax.add_patch(Rectangle((cx - 0.045 + index * 0.022, cy - 0.04), 0.014, height,
                               transform=ax.transAxes, facecolor=color, alpha=0.75))
    ax.add_patch(Circle((cx + 0.035, cy + 0.018), 0.024, transform=ax.transAxes,
                        fill=False, edgecolor=color, lw=1.1))
    ax.plot([cx + 0.052, cx + 0.069], [cy, cy - 0.022], transform=ax.transAxes,
            color=color, lw=1.2)


def draw_workflow_card(ax: plt.Axes, x: float, width: float, number: int,
                       title: str, body: str, color: str, icon) -> None:
    y, height = 0.08, 0.82
    ax.add_patch(FancyBboxPatch((x, y), width, height,
                               boxstyle="round,pad=0.008,rounding_size=0.022",
                               transform=ax.transAxes, facecolor="#FBFCFE",
                               edgecolor=color, linewidth=1.1))
    ax.add_patch(Circle((x + 0.026, y + height - 0.035), 0.016,
                        transform=ax.transAxes, facecolor=color, edgecolor="none"))
    ax.text(x + 0.026, y + height - 0.035, str(number), transform=ax.transAxes,
            ha="center", va="center", color="white", fontsize=7.2, fontweight="bold")
    icon(ax, x + width / 2, y + height * 0.72, color)
    ax.text(x + width / 2, y + height * 0.50, title, transform=ax.transAxes,
            ha="center", va="center", fontsize=7.7, fontweight="bold", color=color,
            linespacing=1.05)
    ax.plot([x + 0.018, x + width - 0.018], [y + height * 0.40, y + height * 0.40],
            transform=ax.transAxes, color=color, alpha=0.35, lw=0.7)
    ax.text(x + width / 2, y + height * 0.19, body, transform=ax.transAxes,
            ha="center", va="center", fontsize=6.4, color="#263238", linespacing=1.15)


def figure1() -> None:
    truth = pd.read_csv(TRUTH_FILE)
    funnel = pd.read_csv(FUNNEL_FILE).set_index("species")
    matching = read_table("S3_Matching_levels")
    source = read_table("S4_Source_support")
    condition = read_table("S5_Condition_metadata")
    total = len(truth)
    ecoli_n = int((truth["species"] == "ecoli").sum())
    yeast_n = int((truth["species"] == "yeast").sum())
    supported_n = int((truth["experimental_substrate_support"] == "substrate_supported").sum())
    ambiguous_n = total - supported_n
    candidate_n = int(funnel["enzyme_substrate_entries"].sum())
    matched_n = int(funnel["experimental_truth_rows"].sum())

    fig = plt.figure(figsize=(14.2, 9.3))
    outer = fig.add_gridspec(2, 5, height_ratios=[1.00, 1.05],
                             left=0.045, right=0.985, bottom=0.07, top=0.96,
                             hspace=0.38, wspace=0.42)
    ax_a = fig.add_subplot(outer[0, :3])
    ax_b = fig.add_subplot(outer[0, 3:])
    ax_c = fig.add_subplot(outer[1, 0:2])
    ax_d = fig.add_subplot(outer[1, 2:4])
    ax_e = fig.add_subplot(outer[1, 4])

    ax_a.set_axis_off()
    panel_label(ax_a, "a", x=-0.03, y=1.04)
    ax_a.set_title("Benchmark construction workflow", loc="left", pad=8)
    cards = [
        ("Genome-scale\nmodels", "eciML1515 (E. coli)\nYeast9 (S. cerevisiae)", "#225EA8", draw_network_icon),
        ("All-reactant\ncandidates", f"{candidate_n:,} associations\n8,797 E. coli; 7,816 yeast", "#168C7B", draw_document_icon),
        ("Protein and\nmolecular mapping", "UniProt single-protein sequence\nmodel-reactant SMILES", "#2D64B3", draw_mapping_icon),
        ("Experimental kcat\nmatching", f"BRENDA + SABIO-RK\n{matched_n:,} positive matches", "#6F49A8", draw_database_icon),
        ("Final benchmark\nfiltering", f"{matched_n - total} records removed\nsequence + valid SMILES", "#B17A16", draw_filter_icon),
        ("Final benchmark\nand evaluation", f"{total:,} records\n12 methods; 6 regimes", "#1F4E8C", draw_chart_icon),
    ]
    gap = 0.014
    width = (1.0 - gap * 5) / 6
    for index, (title, body, color, icon) in enumerate(cards):
        x = index * (width + gap)
        draw_workflow_card(ax_a, x, width, index + 1, title, body, color, icon)
        if index < 5:
            ax_a.annotate("", xy=(x + width + gap * 0.82, 0.49),
                          xytext=(x + width + gap * 0.18, 0.49),
                          xycoords=ax_a.transAxes, textcoords=ax_a.transAxes,
                          arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#37474F"))

    ax_b.set_axis_off()
    panel_label(ax_b, "b", x=-0.08, y=1.04)
    ax_b.set_title("Benchmark scale and biochemical breadth", loc="left", pad=8)
    cards_b = [
        ("Benchmark records", f"{total:,}", "model-linked records", "#225EA8"),
        ("Distinct reactions", f"{truth['reaction_id'].nunique():,}", "model reactions", "#E87512"),
        ("Unique protein sequences", f"{truth['sequence'].nunique():,}", "sequences", "#2E9B43"),
        ("EC annotations", f"{truth['ec_number'].nunique():,}", "distinct strings", "#D92F2F"),
        ("Experimental kcat range", "1.67×10⁻⁴ – 5.7×10⁵", "s⁻¹", "#8D5BC1"),
        ("Median log10(kcat)", f"{truth['true_kcat_log10'].median():.3f}", ">9 orders of magnitude", "#8C5A4A"),
    ]
    box_w, box_h = 0.30, 0.23
    for idx, (title, value, note, color) in enumerate(cards_b):
        row, col = divmod(idx, 3)
        x, y = 0.015 + col * 0.325, 0.64 - row * 0.29
        ax_b.add_patch(FancyBboxPatch((x, y), box_w, box_h,
                                     boxstyle="round,pad=0.008,rounding_size=0.018",
                                     transform=ax_b.transAxes, facecolor="#FBFCFE",
                                     edgecolor=color, lw=1.0))
        ax_b.text(x + box_w / 2, y + box_h * 0.75, title, transform=ax_b.transAxes,
                  ha="center", va="center", fontsize=7.0, fontweight="bold")
        value_size = 12.5 if idx != 4 else 8.3
        ax_b.text(x + box_w / 2, y + box_h * 0.46, value, transform=ax_b.transAxes,
                  ha="center", va="center", fontsize=value_size, fontweight="bold", color=color)
        ax_b.text(x + box_w / 2, y + box_h * 0.17, note, transform=ax_b.transAxes,
                  ha="center", va="center", fontsize=6.2, color="#37474F")
    ax_b.text(0.5, 0.135, "Species composition", transform=ax_b.transAxes,
              ha="center", fontsize=7.3, fontweight="bold")
    ax_b.plot([0.15, 0.43], [0.08, 0.08], transform=ax_b.transAxes, color="#168C7B", lw=4)
    ax_b.plot([0.57, 0.85], [0.08, 0.08], transform=ax_b.transAxes, color="#225EA8", lw=4)
    ax_b.text(0.29, 0.045, f"{ecoli_n}  E. coli", transform=ax_b.transAxes,
              ha="center", color="#168C7B", fontsize=8, fontweight="bold")
    ax_b.text(0.71, 0.045, f"{yeast_n}  S. cerevisiae", transform=ax_b.transAxes,
              ha="center", color="#225EA8", fontsize=8, fontweight="bold")
    ax_b.text(0.5, -0.025, f"Experimental substrate support: {supported_n} supported | {ambiguous_n} participant-ambiguous",
              transform=ax_b.transAxes, ha="center", fontsize=6.5, color="#4E5A63")

    species = ["ecoli", "yeast"]
    labels = [f"E. coli\n(total {ecoli_n})", f"S. cerevisiae\n(total {yeast_n})"]
    match_levels = [
        ("Species + EC + UniProt + substrate ID", "UniProt + substrate ID", "#2C7FB8"),
        ("Species + EC + substrate ID", "Substrate ID", "#74A9CF"),
        ("Species + EC + substrate name", "Substrate name", "#D0E1F2"),
    ]
    bottom = np.zeros(2)
    for level, legend, color in match_levels:
        values = np.array([int(matching.loc[(matching["Species"] == sp)
                                            & (matching["Matching level"] == level), "Records"].sum())
                           for sp in species])
        bars = ax_c.bar(labels, values, bottom=bottom, width=0.58, color=color, label=legend)
        for index, value in enumerate(values):
            if value:
                ax_c.text(index, bottom[index] + value / 2, str(value), ha="center", va="center",
                          color="white" if color != "#D0E1F2" else "#222222",
                          fontsize=7.5, fontweight="bold")
        bottom += values
    panel_label(ax_c, "c")
    ax_c.set_title("Hierarchical matching support")
    ax_c.set_ylabel("Benchmark records")
    ax_c.set_ylim(0, 1000)
    ax_c.legend(frameon=True, framealpha=0.92, edgecolor="#D0D0D0",
                loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=1, fontsize=6.4)
    style_axes(ax_c, "y")

    source_levels = [("BRENDA", "BRENDA only", "#2C7FB8"),
                     ("SABIO-RK", "SABIO-RK only", "#2CA02C"),
                     ("BRENDA;SABIO-RK", "Both databases", "#9467BD")]
    bottom = np.zeros(2)
    for level, legend, color in source_levels:
        values = np.array([int(source.loc[(source["Species"] == sp)
                                          & (source["Source database"] == level), "Records"].sum())
                           for sp in species])
        ax_d.bar(labels, values, bottom=bottom, width=0.58, color=color, label=legend)
        for index, value in enumerate(values):
            if value:
                ax_d.text(index, bottom[index] + value / 2, str(value), ha="center", va="center",
                          color="white", fontsize=7.5, fontweight="bold")
        bottom += values
    panel_label(ax_d, "d")
    ax_d.set_title("Experimental source-database support by species")
    ax_d.set_ylabel("Benchmark records")
    ax_d.set_ylim(0, 1000)
    ax_d.legend(frameon=True, framealpha=0.92, edgecolor="#D0D0D0",
                loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=1, fontsize=6.4)
    style_axes(ax_d, "y")

    x = np.arange(2)
    width_bar = 0.34
    ph = condition["pH available (%)"].to_numpy(float)
    temp = condition["Temperature available (%)"].to_numpy(float)
    bars1 = ax_e.bar(x - width_bar / 2, ph, width_bar, color="#2C7FB8", label="pH")
    bars2 = ax_e.bar(x + width_bar / 2, temp, width_bar, color="#2CA02C", label="Temperature")
    panel_label(ax_e, "e", x=-0.28)
    ax_e.set_title("Experimental-condition metadata", fontsize=9.2)
    ax_e.set_ylabel("Availability (%)")
    ax_e.set_xticks(x, ["E. coli", "S. cerevisiae"], rotation=20, ha="right")
    ax_e.set_ylim(0, 105)
    ax_e.legend(frameon=False, loc="upper center", ncol=2, fontsize=6.8)
    for bars in (bars1, bars2):
        for bar in bars:
            ax_e.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.2,
                      f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=6.5,
                      fontweight="bold")
    style_axes(ax_e, "y")

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    matching.to_csv(SOURCE_DIR / "Figure1c_matching_support.csv", index=False)
    source.to_csv(SOURCE_DIR / "Figure1d_source_support.csv", index=False)
    condition.to_csv(SOURCE_DIR / "Figure1e_condition_metadata.csv", index=False)
    pd.DataFrame({
        "metric": ["records", "reactions", "protein_sequences", "uniprot_accessions", "ec_annotations",
                   "substrate_supported", "participant_ambiguous"],
        "value": [total, truth["reaction_id"].nunique(), truth["sequence"].nunique(),
                  truth["uniprot_id"].nunique(), truth["ec_number"].nunique(), supported_n, ambiguous_n],
    }).to_csv(SOURCE_DIR / "Figure1ab_summary.csv", index=False)
    for label, axis in zip("abcde", [ax_a, ax_b, ax_c, ax_d, ax_e]):
        save_axis_panel(fig, axis, f"Figure1{label}")
    save_composite(fig, 1)


def figure2() -> None:
    full = read_table("S1_Full_metrics").set_index("Method").loc[METHOD_ORDER].reset_index()
    row_ci = read_table("S11_Bootstrap_CI")
    cluster = read_table("S18_Cluster_bootstrap")
    pair_ci = cluster.loc[(cluster["analysis_scope"] == "achieved_evaluation_set")
                          & (cluster["cluster_type"] == "pair")].set_index("method")
    full["display_method"] = full["Method"].map(METHOD_DISPLAY)
    full["color"] = full["Method"].map(METHOD_COLORS)
    full["inference_regime"] = full["Inference regime"]

    fig, axes = plt.subplots(3, 2, figsize=(13.6, 13.7))
    fig.subplots_adjust(left=0.075, right=0.93, bottom=0.055, top=0.965,
                        wspace=0.40, hspace=0.46)
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes.ravel()

    panel_label(ax_a, "a")
    style_axes(ax_a)
    for regime in REGIME_ORDER:
        part = full.loc[full["inference_regime"] == regime]
        ax_a.scatter(part["Coverage (%)"], part["MAE log10"], s=45,
                     color=REGIME_COLORS[regime], edgecolor="white", lw=0.7, zorder=3)
    offsets = {
        "KcatNet": (14, 0), "CataPro": (14, -2), "UniKP": (14, 2),
        "SELFprot": (14, -7), "DLKcat": (14, 5), "PreTKcat": (-48, -3),
        "TurNuP": (-48, 8), "PMAK": (14, -5), "KinForm-L": (-10, 0),
        "CatPred": (14, 7), "DEKP-public-retrained": (-42, 7), "GO-HKP": (-46, -8),
    }
    for _, row in full.iterrows():
        method = str(row["Method"])
        dx, dy = offsets[method]
        ax_a.annotate(display_method(method),
                      (float(row["Coverage (%)"]), float(row["MAE log10"])),
                      xytext=(dx, dy), textcoords="offset points",
                      ha="right" if dx < 0 else "left", va="center", fontsize=8.0,
                      arrowprops=dict(arrowstyle="-", color="#444444", lw=0.55,
                                      shrinkA=2, shrinkB=4), annotation_clip=False)
    ax_a.set_title("Coverage versus available-case MAE", pad=8)
    ax_a.set_xlabel("Coverage of the 1,246-record benchmark (%)")
    ax_a.set_ylabel(r"MAE of log$_{10}$($k_{cat}$)")
    ax_a.set_xlim(53, 102)
    ax_a.set_ylim(0.66, 1.10)
    regime_legend(ax_a, loc="upper left")

    panel_label(ax_b, "b")
    y = np.arange(len(full))
    ax_b.barh(y, full["n"], color=full["color"], height=0.62, zorder=2)
    ax_b.set_yticks(y, full["display_method"])
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 1430)
    ax_b.set_xlabel("Evaluated records (n of 1,246)")
    ax_b.set_title("Applicable benchmark records", pad=8)
    ax_b.axvline(1246, color="#8C8C8C", linestyle=(0, (4, 3)), lw=0.8)
    for index, row in full.iterrows():
        ax_b.text(float(row["n"]) + 15, index,
                  f"{int(row['n']):,} ({float(row['Coverage (%)']):.1f}%)",
                  va="center", fontsize=7.7)
    group_separators(ax_b, y)
    style_axes(ax_b, "x")

    panel_label(ax_c, "c")
    estimate = full["MAE log10"].to_numpy(float)
    low = np.array([pair_ci.loc[m, "cluster_bootstrap_ci_low_95"] for m in full["Method"]])
    high = np.array([pair_ci.loc[m, "cluster_bootstrap_ci_high_95"] for m in full["Method"]])
    y = np.arange(len(full))[::-1]
    for value, low_value, high_value, ypos, color in zip(
            estimate, low, high, y, full["color"]):
        ax_c.errorbar(value, ypos,
                      xerr=np.array([[value - low_value], [high_value - value]]),
                      fmt="o", color=color, ecolor=color, elinewidth=1.0,
                      capsize=2.2, markersize=4.2, zorder=3)
    ax_c.set_yticks(y, full["display_method"])
    for value, ypos in zip(estimate, y):
        ax_c.text(value + 0.008, ypos, f"{value:.3f}", va="center", fontsize=7.0)
    ax_c.set_title("Absolute prediction error", pad=8)
    ax_c.set_xlabel(r"MAE of log$_{10}$($k_{cat}$)")
    group_separators(ax_c, y[::-1])
    style_axes(ax_c, "x")

    panel_label(ax_d, "d")
    spearman = full["Spearman"].to_numpy(float)
    lows, highs = [], []
    for method in full["Method"]:
        row = row_ci.loc[(row_ci["method"] == method) & (row_ci["metric"] == "spearman_log10")].iloc[0]
        lows.append(float(row["bootstrap_ci_low_95"]))
        highs.append(float(row["bootstrap_ci_high_95"]))
    lows, highs = np.array(lows), np.array(highs)
    for value, low_value, high_value, ypos, color in zip(
            spearman, lows, highs, y, full["color"]):
        ax_d.errorbar(value, ypos,
                      xerr=np.array([[value - low_value], [high_value - value]]),
                      fmt="o", color=color, ecolor=color, elinewidth=1.0,
                      capsize=2.2, markersize=4.2, zorder=3)
    ax_d.set_yticks(y, full["display_method"])
    for value, ypos in zip(spearman, y):
        ax_d.text(value + 0.008, ypos, f"{value:.3f}", va="center", fontsize=7.0)
    ax_d.set_title("Rank correlation", pad=8)
    ax_d.set_xlabel(r"Spearman $\rho$")
    ax_d.set_xlim(0.05, 0.64)
    group_separators(ax_d, y[::-1])
    style_axes(ax_d, "x")

    panel_label(ax_e, "e")
    two = full["Within 2-fold (%)"].to_numpy(float)
    ten = full["Within 10-fold (%)"].to_numpy(float)
    for values, metric, filled, label in [
        (two, "within_0.3_fraction", False, "within two-fold"),
        (ten, "within_1.0_fraction", True, "within ten-fold"),
    ]:
        lows, highs = [], []
        for method in full["Method"]:
            row = row_ci.loc[(row_ci["method"] == method) & (row_ci["metric"] == metric)].iloc[0]
            lows.append(100 * float(row["bootstrap_ci_low_95"]))
            highs.append(100 * float(row["bootstrap_ci_high_95"]))
        lows, highs = np.array(lows), np.array(highs)
        for value, low_value, high_value, ypos, color in zip(
                values, lows, highs, y, full["color"]):
            ax_e.errorbar(value, ypos,
                          xerr=np.array([[value - low_value], [high_value - value]]),
                          fmt="o", color=color, ecolor=color, elinewidth=1.0,
                          markerfacecolor=color if filled else "white",
                          markeredgecolor=color, capsize=2.2, markersize=4.8,
                          zorder=3)
    ax_e.set_yticks(y, full["display_method"])
    ax_e.set_title("Predictions within two- and ten-fold", pad=8)
    ax_e.set_xlabel("Predictions within threshold (%)")
    group_separators(ax_e, y[::-1])
    style_axes(ax_e, "x")
    ax_e.legend(handles=[
        Line2D([0], [0], marker="o", color="#777777", markerfacecolor="white",
               linestyle="none", label="within two-fold"),
        Line2D([0], [0], marker="o", color="#777777", markerfacecolor="#777777",
               linestyle="none", label="within ten-fold"),
    ], frameon=True, loc="lower right")

    panel_label(ax_f, "f")
    bias = full["Bias log10"].to_numpy(float)
    ax_f.scatter(bias, y, c=full["color"], s=35, zorder=3)
    ax_f.axvline(0, color="#333333", lw=0.9)
    ax_f.set_yticks(y, full["display_method"])
    for value, ypos in zip(bias, y):
        ax_f.text(value + (0.025 if value >= 0 else -0.025), ypos, f"{value:.3f}",
                  va="center", ha="left" if value >= 0 else "right", fontsize=7.1)
    ax_f.set_xlim(-0.95, 1.0)
    ax_f.set_title("Directional calibration bias", pad=8)
    ax_f.set_xlabel(r"Mean signed error (predicted − observed log$_{10}$($k_{cat}$))")
    group_separators(ax_f, y[::-1])
    style_axes(ax_f, "x")
    ax_f.text(0.02, -0.14, "Underestimation", transform=ax_f.transAxes, style="italic",
              fontsize=7.2, color="#777777")
    ax_f.text(0.98, -0.14, "Overestimation", transform=ax_f.transAxes, style="italic",
              fontsize=7.2, color="#777777", ha="right")

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    full.drop(columns=["color"]).to_csv(SOURCE_DIR / "Figure2_plot_data.csv", index=False)
    pair_ci.reset_index().to_csv(SOURCE_DIR / "Figure2_pair_cluster_intervals.csv", index=False)
    row_ci.to_csv(SOURCE_DIR / "Figure2_row_bootstrap_intervals.csv", index=False)
    for label, axis in zip("abcdef", [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]):
        save_axis_panel(fig, axis, f"Figure2{label}")
    save_composite(fig, 2)


def available_scope_rows() -> dict[str, pd.DataFrame]:
    full = read_table("S1_Full_metrics").rename(columns={"Method": "method", "MAE log10": "mae_log10"})
    full["n_scope"] = full["n"]
    reaction = read_table("S6_Reaction_subset").rename(
        columns={"mae_log10_common_subset": "mae_log10", "n_common_subset": "n_scope"}
    )
    available = read_table("S7_Available_case")
    catpred_scope = next(x for x in available["subset"].unique() if str(x).startswith("catpred_accessible_scope_"))
    kinform_scope = next(x for x in available["subset"].unique() if str(x).startswith("kinform_accessible_scope_"))
    return {
        "Achieved set": full[["method", "mae_log10", "n_scope"]],
        "Reaction scope": reaction[["method", "mae_log10", "n_scope"]],
        "CatPred-accessible": available.loc[available["subset"] == catpred_scope,
                                             ["method", "mae_log10", "n"]].rename(columns={"n": "n_scope"}),
        "KinForm-L-accessible": available.loc[available["subset"] == kinform_scope,
                                               ["method", "mae_log10", "n"]].rename(columns={"n": "n_scope"}),
    }


def figure3() -> None:
    scopes = available_scope_rows()
    selected = {
        "Achieved set": METHOD_ORDER,
        "Reaction scope": ["KcatNet", "CataPro", "TurNuP", "PMAK"],
        "CatPred-accessible": ["KcatNet", "CataPro", "TurNuP", "PMAK", "CatPred"],
        "KinForm-L-accessible": ["KcatNet", "CataPro", "TurNuP", "PMAK", "KinForm-L"],
    }
    paired = read_table("S13_Paired_bootstrap")
    ranks = read_table("S12_Rank_stability")
    rank_scope = next(x for x in ranks["comparison_set"].unique()
                      if str(x).startswith("reaction_aware_common_"))
    ranks = ranks.loc[ranks["comparison_set"] == rank_scope].copy()

    fig = plt.figure(figsize=(13.7, 9.4))
    outer = fig.add_gridspec(2, 2, height_ratios=[1.42, 0.86],
                             left=0.055, right=0.985, bottom=0.075, top=0.94,
                             hspace=0.38, wspace=0.34)
    top = outer[0, :].subgridspec(1, 4, wspace=0.38)
    scope_axes = [fig.add_subplot(top[0, i]) for i in range(4)]
    ax_b = fig.add_subplot(outer[1, 0])
    ax_c = fig.add_subplot(outer[1, 1])
    panel_a_label = fig.text(0.018, 0.975, "a", fontsize=15, fontweight="bold", va="top")
    panel_a_title = fig.text(0.055, 0.972, "Available-case MAE within defined scopes",
                             fontsize=11, fontweight="bold", va="top")

    source_rows = []
    scope_subtitles = {
        "Achieved set": "Achieved predictions",
        "Reaction scope": "Reaction scope (1,047)",
        "CatPred-accessible": "CatPred scope (1,156)",
        "KinForm-L-accessible": "KinForm-L scope (729)",
    }
    for axis, scope_name in zip(scope_axes, scopes):
        frame = scopes[scope_name].set_index("method").loc[selected[scope_name]].reset_index()
        frame = frame.sort_values("mae_log10", ascending=False)
        y = np.arange(len(frame))
        colors = [METHOD_COLORS[m] for m in frame["method"]]
        bars = axis.barh(y, frame["mae_log10"], color=colors, height=0.62)
        axis.set_yticks(y, [display_method(m) for m in frame["method"]])
        axis.set_xlim(0, 1.13)
        axis.set_xlabel(r"MAE in log$_{10}$($k_{cat}$)")
        axis.set_title(scope_subtitles[scope_name], fontsize=9.2, pad=7)
        style_axes(axis, "x")
        for bar, row in zip(bars, frame.itertuples()):
            axis.text(min(row.mae_log10 - 0.015, 1.08), bar.get_y() + bar.get_height() / 2,
                      f"{row.mae_log10:.3f} (n={int(row.n_scope):,})",
                      ha="right", va="center", color="white" if row.mae_log10 > 0.72 else "#222222",
                      fontsize=6.4, fontweight="bold")
        frame = frame.assign(scope=scope_name)
        source_rows.append(frame)

    panel_label(ax_b, "b", x=-0.12)
    labels = [f"{row.method_a} − {row.method_b}" for row in paired.itertuples()]
    y = np.arange(len(paired))[::-1]
    estimates = paired["mae_difference_a_minus_b"].to_numpy(float)
    low = paired["bootstrap_ci_low_95"].to_numpy(float)
    high = paired["bootstrap_ci_high_95"].to_numpy(float)
    ax_b.errorbar(estimates, y, xerr=np.vstack([estimates - low, high - estimates]),
                  fmt="o", color="#1F5CC4", ecolor="#1F5CC4", capsize=3,
                  markersize=4.8, lw=1.0)
    ax_b.axvline(0, color="#777777", linestyle="--", lw=0.9)
    ax_b.set_yticks(y, labels)
    ax_b.set_title("Paired MAE differences on 1,047 common records", pad=8)
    ax_b.set_xlabel("MAE difference (first minus second)")
    ax_b.text(0.5, -0.17, "Negative values favor the first-named method.",
              transform=ax_b.transAxes, ha="center", fontsize=7.3, style="italic")
    style_axes(ax_b, "x")

    panel_label(ax_c, "c", x=-0.12)
    order = ["PMAK", "TurNuP", "KcatNet"]
    ranks = ranks.set_index("method").loc[order].reset_index()
    values = 100 * ranks["rank_1_bootstrap_frequency"].to_numpy(float)
    y = np.arange(len(ranks))
    for ypos, value, method in zip(y, values, ranks["method"]):
        ax_c.hlines(ypos, 0, value, color=METHOD_COLORS[method], lw=1.5)
        ax_c.scatter(value, ypos, s=38, color=METHOD_COLORS[method], zorder=3)
        ax_c.text(value + 1, ypos, f"{value:.1f}%", va="center", fontsize=8)
    ax_c.set_yticks(y, ranks["method"])
    ax_c.set_xlim(0, 72)
    ax_c.set_title("Rank-one frequency in paired bootstrap", pad=8)
    ax_c.set_xlabel("Bootstrap replicates ranked first (%)")
    ax_c.text(0.5, -0.17, "Resampling stability within this shared set; not universal-best probability.",
              transform=ax_c.transAxes, ha="center", fontsize=7.0, style="italic")
    style_axes(ax_c, "x")

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat(source_rows, ignore_index=True).to_csv(SOURCE_DIR / "Figure3a_available_case_scopes.csv", index=False)
    paired.to_csv(SOURCE_DIR / "Figure3b_paired_bootstrap.csv", index=False)
    ranks.to_csv(SOURCE_DIR / "Figure3c_rank_stability.csv", index=False)
    save_axes_panel(fig, scope_axes, "Figure3a", [panel_a_label, panel_a_title])
    save_axis_panel(fig, ax_b, "Figure3b")
    save_axis_panel(fig, ax_c, "Figure3c")
    save_composite(fig, 3)


def subgroup_value(frame: pd.DataFrame, method: str, feature: str, group: str,
                   column: str) -> float:
    row = frame.loc[(frame["method"] == method) & (frame["feature"] == feature)
                    & (frame["group"].astype(str) == group)]
    if len(row) != 1:
        raise ValueError(f"Missing subgroup: {method}, {feature}, {group}")
    return float(row.iloc[0][column])


def protein_cluster_large_error_ci(methods: list[str], replicates: int = 2000,
                                   seed: int = 20260811) -> pd.DataFrame:
    table0 = pd.read_csv(TABLE0_FILE, usecols=["method", "entry_id", "species",
                                               "prediction_status", "absolute_error_log10"])
    audit = pd.read_csv(AUDIT_FILE, usecols=["entry_id", "protein_cluster"])
    frame = table0.loc[table0["method"].isin(methods)
                       & (table0["prediction_status"] == "predicted")].merge(
                           audit, on="entry_id", validate="many_to_one")
    rng = np.random.default_rng(seed)
    rows = []
    for method in methods:
        for species in ["ecoli", "yeast"]:
            sub = frame.loc[(frame["method"] == method) & (frame["species"] == species)].copy()
            groups = [group["absolute_error_log10"].to_numpy(float)
                      for _, group in sub.groupby("protein_cluster", sort=False)]
            observed = float((sub["absolute_error_log10"] > 1).mean())
            samples = np.empty(replicates)
            for index in range(replicates):
                chosen = rng.integers(0, len(groups), len(groups))
                values = np.concatenate([groups[i] for i in chosen])
                samples[index] = np.mean(values > 1)
            rows.append({
                "method": method,
                "species": species,
                "n": len(sub),
                "n_protein_clusters": len(groups),
                "large_error_fraction": observed,
                "ci_low_95": float(np.quantile(samples, 0.025)),
                "ci_high_95": float(np.quantile(samples, 0.975)),
                "bootstrap_replicates": replicates,
                "seed": seed,
            })
    return pd.DataFrame(rows)


def figure4() -> None:
    stratified = read_table("S9_Error_stratification")
    methods = ["KcatNet", "CataPro", "TurNuP", "PMAK", "CatPred", "KinForm-L"]
    method_labels = [display_method(m) for m in methods]
    categories = [
        ("species", "ecoli", "E. coli"),
        ("species", "yeast", "S. cerevisiae"),
        ("experimental_substrate_support", "substrate_supported", "Substrate\nsupported"),
        ("experimental_substrate_support", "participant_ambiguous", "Participant\nambiguous"),
        ("substrate_role_group_substrate_supported", "other_reactant", "Supported\nother reactant"),
        ("substrate_role_group_substrate_supported", "currency_or_cofactor", "Supported\ncurrency/cofactor"),
        ("substrate_role_group_substrate_supported", "carrier_linked_variable", "Supported\ncarrier-linked"),
    ]
    matrix = np.array([[subgroup_value(stratified, method, feature, group,
                                       "mean_abs_error_log10")
                        for feature, group, _ in categories] for method in methods])
    counts = np.array([[int(subgroup_value(stratified, method, feature, group, "n"))
                        for feature, group, _ in categories] for method in methods])
    large = protein_cluster_large_error_ci(methods)

    fig = plt.figure(figsize=(13.7, 8.7))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.92],
                            left=0.075, right=0.96, bottom=0.075, top=0.94,
                            hspace=0.42, wspace=0.34)
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])

    panel_label(ax_a, "a", x=-0.055)
    im = ax_a.imshow(matrix, aspect="auto", cmap="YlOrRd",
                     vmin=np.nanmin(matrix) - 0.03, vmax=np.nanmax(matrix) + 0.03)
    ax_a.set_xticks(np.arange(len(categories)), [x[2] for x in categories])
    ax_a.set_yticks(np.arange(len(methods)), method_labels)
    ax_a.set_title("Stratified mean absolute error", loc="left", pad=30)
    midpoint = (np.nanmin(matrix) + np.nanmax(matrix)) / 2
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax_a.text(j, i, f"{matrix[i, j]:.3f}\n(n={counts[i, j]})",
                      ha="center", va="center", fontsize=7.1,
                      color="white" if matrix[i, j] > midpoint else "black")
    for xline in [1.5, 3.5]:
        ax_a.axvline(xline, color="white", lw=2.5)
    for center, title in [(0.5, "Species"), (2.5, "Experimental substrate support"),
                          (5.0, "Registry-defined role; supported records only")]:
        ax_a.text(center, -0.66, title, ha="center", va="center", fontsize=8.2,
                  fontweight="bold")
    colorbar = fig.colorbar(im, ax=ax_a, fraction=0.022, pad=0.018)
    colorbar.set_label(r"MAE in log$_{10}$($k_{cat}$)")

    panel_label(ax_b, "b", x=-0.12)
    y = np.arange(len(methods))[::-1]
    species_contrast, role_contrast = [], []
    for method in methods:
        species_contrast.append(
            subgroup_value(stratified, method, "species", "ecoli", "mean_abs_error_log10")
            - subgroup_value(stratified, method, "species", "yeast", "mean_abs_error_log10")
        )
        role_contrast.append(
            subgroup_value(stratified, method, "substrate_role_group_substrate_supported",
                           "other_reactant", "mean_abs_error_log10")
            - subgroup_value(stratified, method, "substrate_role_group_substrate_supported",
                             "currency_or_cofactor", "mean_abs_error_log10")
        )
    for ypos, value_a, value_b in zip(y, species_contrast, role_contrast):
        ax_b.hlines(ypos, min(value_a, value_b), max(value_a, value_b), color="#AEB8C2", lw=1.3)
    ax_b.scatter(species_contrast, y, color="#1F5CC4", s=34,
                 label="E. coli − S. cerevisiae", zorder=3)
    ax_b.scatter(role_contrast, y, color="#F57C00", marker="D", s=30,
                 label="Supported other − currency/cofactor", zorder=3)
    ax_b.axvline(0, color="#777777", linestyle="--", lw=0.9)
    ax_b.set_yticks(y, method_labels)
    ax_b.set_title("Descriptive subgroup MAE contrasts", pad=8)
    ax_b.set_xlabel("MAE difference (log10 units)")
    ax_b.legend(frameon=True, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    style_axes(ax_b, "x")

    panel_label(ax_c, "c", x=-0.12)
    offsets = {"ecoli": 0.14, "yeast": -0.14}
    species_colors = {"ecoli": "#1F5CC4", "yeast": "#D95F5F"}
    species_labels = {"ecoli": "E. coli", "yeast": "S. cerevisiae"}
    for method, ypos in zip(methods, y):
        for species in ["ecoli", "yeast"]:
            row = large.loc[(large["method"] == method) & (large["species"] == species)].iloc[0]
            value = float(row["large_error_fraction"])
            ax_c.errorbar(value, ypos + offsets[species],
                          xerr=np.array([[value - float(row["ci_low_95"])],
                                         [float(row["ci_high_95"]) - value]]),
                          fmt="o", color=species_colors[species], capsize=2.3, lw=1.0,
                          markersize=4.2, label=species_labels[species]
                          if method == methods[0] else None)
    ax_c.set_yticks(y, method_labels)
    ax_c.set_title("Large-error fraction by species", pad=8)
    ax_c.set_xlabel(r"Fraction with absolute log$_{10}$ error > 1")
    ax_c.legend(frameon=True, loc="lower right")
    style_axes(ax_c, "x")

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    stratified.to_csv(SOURCE_DIR / "Figure4a_stratified_mae.csv", index=False)
    pd.DataFrame({
        "method": methods,
        "ecoli_minus_yeast_mae": species_contrast,
        "supported_other_minus_currency_cofactor_mae": role_contrast,
    }).to_csv(SOURCE_DIR / "Figure4b_mae_contrasts.csv", index=False)
    large.to_csv(SOURCE_DIR / "Figure4c_large_error_protein_cluster_bootstrap.csv", index=False)
    for label, axis in zip("abc", [ax_a, ax_b, ax_c]):
        save_axis_panel(fig, axis, f"Figure4{label}")
    save_composite(fig, 4)


def make_packages() -> None:
    PACKAGE_DIR.mkdir(parents=True, exist_ok=True)
    script = Path(__file__)
    descriptions = {
        1: "Figure1ab_summary.csv Figure1c_matching_support.csv Figure1d_source_support.csv Figure1e_condition_metadata.csv",
        2: "Figure2_plot_data.csv Figure2_pair_cluster_intervals.csv Figure2_row_bootstrap_intervals.csv",
        3: "Figure3a_available_case_scopes.csv Figure3b_paired_bootstrap.csv Figure3c_rank_stability.csv",
        4: "Figure4a_stratified_mae.csv Figure4b_mae_contrasts.csv Figure4c_large_error_protein_cluster_bootstrap.csv",
    }
    for number, names in descriptions.items():
        archive = PACKAGE_DIR / f"Figure{number}_data_code_and_panels.zip"
        with ZipFile(archive, "w", compression=ZIP_DEFLATED) as zf:
            zf.write(script, arcname=script.name)
            zf.writestr(
                "README.txt",
                "Generated from audited v1.2.0 tables in the visual language of the 0806 draft.\n"
                "Run from the repository root with:\n"
                "  python paper/generate_manuscript_figures_0806_style_r4.py\n"
                f"Source tables included: {names}\n"
                "Figure 1a and 1b were redrawn using paper/0806/0806/Figure1_1a.png and "
                "Figure1_1b.png as visual references, with all counts updated.\n",
            )
            for source_name in names.split():
                path = SOURCE_DIR / source_name
                zf.write(path, arcname=f"data/{path.name}")
            zf.write(OUTPUT_DIR / f"Figure{number}.png", arcname=f"Figure{number}.png")
            zf.write(OUTPUT_DIR / f"Figure{number}.pdf", arcname=f"Figure{number}.pdf")
            for path in sorted(PANEL_DIR.glob(f"Figure{number}*")):
                zf.write(path, arcname=f"panels/{path.name}")


def main() -> None:
    configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    figure1()
    figure2()
    figure3()
    figure4()
    make_packages()
    shutil.copy2(Path(__file__), OUTPUT_DIR.parent / Path(__file__).name)
    print(f"Wrote 0806-style v1.2.0-r4 figures to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
