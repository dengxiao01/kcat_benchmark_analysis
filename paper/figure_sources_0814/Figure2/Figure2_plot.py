#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reproduce Figure 2a-f as separate publication-ready panels.
Input: panel_a.csv ... panel_f_ratio.csv
Output: PNG, PDF and SVG for every panel.

Tested with Python 3.10+, matplotlib, numpy.
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
OUT = HERE / "Figure2_panels"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.labelsize": 12,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 8.5,
    "axes.linewidth": 0.9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

COLORS = {
    "Released sequence+substrate checkpoint": "#1f77b4",
    "Temperature-conditioned public retraining": "#2ca02c",
    "Reaction-aware checkpoint": "#ff7f0e",
    "Method-specific checkpoint subset": "#9467bd",
    "Structure-aware public retraining": "#d62728",
    "Functional-assignment baseline": "#7f7f7f",
}

MARKERS = {
    "Released sequence+substrate checkpoint": "o",
    "Temperature-conditioned public retraining": "s",
    "Reaction-aware checkpoint": "D",
    "Method-specific checkpoint subset": "^",
    "Structure-aware public retraining": "P",
    "Functional-assignment baseline": "X",
}

METHOD_ORDER = [
    "KcatNet", "TurNuP", "PMAK", "KinForm-L", "CataPro", "CatPred",
    "UniKP", "PreTKcat", "SELFprot", "DLKcat", "GO-HKP",
    "DEKP-public-retrained"
]

def read_csv(name):
    with open(HERE / name, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def save_all(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

def panel_label(ax, letter):
    ax.text(-0.12, 1.03, letter, transform=ax.transAxes,
            fontsize=18, fontweight="bold", va="bottom")

# ---------------------------------------------------------
# a. Coverage–accuracy landscape
# ---------------------------------------------------------
def plot_a():
    data = read_csv("panel_a.csv")
    fig, ax = plt.subplots(figsize=(7.2, 5.7))

    offsets = {
        "KcatNet": (0.9, -0.020),
        "CataPro": (0.9, 0.012),
        "PreTKcat": (0.9, 0.012),
        "UniKP": (0.9, -0.010),
        "SELFprot": (0.9, 0.014),
        "DLKcat": (0.9, 0.016),
        "TurNuP": (0.9, -0.006),
        "PMAK": (0.9, 0.014),
        "KinForm-L": (1.0, 0.014),
        "CatPred": (1.0, 0.012),
        "GO-HKP": (-9.0, -0.010),
        "DEKP-public-retrained": (-12.0, 0.012),
    }

    for r in data:
        x = float(r["Coverage (%)"])
        y = float(r["MAE log10"])
        reg = r["Inference regime"]
        ax.scatter(x, y, s=85, marker=MARKERS[reg], color=COLORS[reg],
                   edgecolor="black", linewidth=0.6, zorder=3)

        dx, dy = offsets[r["Method"]]
        label = r["Method"].replace("DEKP-public-retrained", "DEKP-public\n-retrained")
        ax.text(x + dx, y + dy, label, fontsize=8.5, va="center")

    handles = [
        Line2D([0], [0], marker=MARKERS[k], linestyle="None",
               markerfacecolor=COLORS[k], markeredgecolor="black",
               markersize=7.5, label={
                   "Released sequence+substrate checkpoint": "Released sequence+substrate",
                   "Temperature-conditioned public retraining": "Temperature-conditioned retraining",
                   "Reaction-aware checkpoint": "Reaction-aware",
                   "Method-specific checkpoint subset": "Method-specific subset",
                   "Structure-aware public retraining": "Structure-aware retraining",
                   "Functional-assignment baseline": "Functional-assignment",
               }[k])
        for k in COLORS
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left")
    ax.set_xlabel("Benchmark coverage (%)")
    ax.set_ylabel(r"MAE on log$_{10}$($k_{cat}$)")
    ax.set_xlim(54, 105)
    ax.set_ylim(0.67, 1.12)
    ax.grid(True, linestyle="--", linewidth=0.55, alpha=0.4)
    panel_label(ax, "a")
    fig.tight_layout()
    save_all(fig, "Figure2a_coverage_MAE")

# ---------------------------------------------------------
# b. MAE with row-bootstrap CI
# ---------------------------------------------------------
def plot_b():
    data = {r["Method"]: r for r in read_csv("panel_b.csv")}
    fig, ax = plt.subplots(figsize=(7.1, 5.9))
    y = np.arange(len(METHOD_ORDER))

    for i, method in enumerate(METHOD_ORDER):
        r = data[method]
        est = float(r["MAE estimate"])
        lo = float(r["CI low 95%"])
        hi = float(r["CI high 95%"])
        reg = r["Inference regime"]
        c = COLORS[reg]
        ax.hlines(i, lo, hi, color=c, linewidth=2.2)
        ax.plot(est, i, "o", color=c, markersize=7.5,
                markeredgecolor="black", markeredgewidth=0.5)
        ax.text(hi + 0.012, i, f'n={r["n"]}', va="center", fontsize=8.3)

    ax.set_yticks(y)
    ax.set_yticklabels(METHOD_ORDER)
    ax.invert_yaxis()
    ax.set_xlabel(r"MAE on log$_{10}$($k_{cat}$) (95% bootstrap CI)")
    ax.set_xlim(0.64, 1.19)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.55, alpha=0.4)
    panel_label(ax, "b")
    fig.tight_layout()
    save_all(fig, "Figure2b_MAE_CI")

# ---------------------------------------------------------
# c. Spearman with row-bootstrap CI
# ---------------------------------------------------------
def plot_c():
    data = {r["Method"]: r for r in read_csv("panel_c.csv")}
    fig, ax = plt.subplots(figsize=(7.1, 5.9))
    y = np.arange(len(METHOD_ORDER))

    for i, method in enumerate(METHOD_ORDER):
        r = data[method]
        est = float(r["Spearman estimate"])
        lo = float(r["CI low 95%"])
        hi = float(r["CI high 95%"])
        reg = r["Inference regime"]
        c = COLORS[reg]
        ax.hlines(i, lo, hi, color=c, linewidth=2.2)
        ax.plot(est, i, "o", color=c, markersize=7.5,
                markeredgecolor="black", markeredgewidth=0.5)
        ax.text(hi + 0.008, i, f'n={r["n"]}', va="center", fontsize=8.3)

    ax.set_yticks(y)
    ax.set_yticklabels(METHOD_ORDER)
    ax.invert_yaxis()
    ax.set_xlabel(r"Spearman correlation on log$_{10}$($k_{cat}$) (95% bootstrap CI)")
    ax.set_xlim(0.04, 0.66)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.55, alpha=0.4)
    panel_label(ax, "c")
    fig.tight_layout()
    save_all(fig, "Figure2c_Spearman_CI")

# ---------------------------------------------------------
# d. Fold agreement
# ---------------------------------------------------------
def plot_d():
    data = {r["Method"]: r for r in read_csv("panel_d.csv")}
    fig, ax = plt.subplots(figsize=(7.3, 5.9))
    y = np.arange(len(METHOD_ORDER))
    offset = 0.14

    for i, method in enumerate(METHOD_ORDER):
        r = data[method]

        est2 = float(r["Within 2-fold (%)"])
        lo2 = float(r["2-fold CI low"])
        hi2 = float(r["2-fold CI high"])
        ax.hlines(i - offset, lo2, hi2, color="#4C78A8", linewidth=2.1)
        ax.plot(est2, i - offset, "o", color="#4C78A8", markersize=6.8,
                markeredgecolor="black", markeredgewidth=0.45)

        est10 = float(r["Within 10-fold (%)"])
        lo10 = float(r["10-fold CI low"])
        hi10 = float(r["10-fold CI high"])
        ax.hlines(i + offset, lo10, hi10, color="#F58518", linewidth=2.1)
        ax.plot(est10, i + offset, "o", color="#F58518", markersize=6.8,
                markeredgecolor="black", markeredgewidth=0.45)

    ax.axvline(50, color="black", linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(METHOD_ORDER)
    ax.invert_yaxis()
    ax.set_xlabel("Predictions within fold threshold (%)")
    ax.set_xlim(10, 82)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.55, alpha=0.35)
    ax.legend(
        [
            Line2D([0], [0], marker="o", color="#4C78A8", label="Within 2-fold"),
            Line2D([0], [0], marker="o", color="#F58518", label="Within 10-fold"),
        ],
        ["Within 2-fold", "Within 10-fold"],
        frameon=False, loc="lower right"
    )
    panel_label(ax, "d")
    fig.tight_layout()
    save_all(fig, "Figure2d_fold_agreement")

# ---------------------------------------------------------
# e. Directional calibration bias
# ---------------------------------------------------------
def plot_e():
    data = {r["Method"]: r for r in read_csv("panel_e.csv")}
    fig, ax = plt.subplots(figsize=(7.1, 5.9))
    y = np.arange(len(METHOD_ORDER))

    for i, method in enumerate(METHOD_ORDER):
        r = data[method]
        val = float(r["Mean signed error log10"])
        reg = r["Inference regime"]
        c = COLORS[reg]
        ax.hlines(i, 0, val, color=c, linewidth=2.3)
        ax.plot(val, i, "o", color=c, markersize=7.5,
                markeredgecolor="black", markeredgewidth=0.5)

    ax.axvline(0, color="black", linewidth=1.1)
    ax.set_yticks(y)
    ax.set_yticklabels(METHOD_ORDER)
    ax.invert_yaxis()
    ax.set_xlabel(r"Mean signed error on log$_{10}$($k_{cat}$)")
    ax.set_xlim(-0.75, 1.02)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.55, alpha=0.35)
    panel_label(ax, "e")
    fig.tight_layout()
    save_all(fig, "Figure2e_signed_error")

# ---------------------------------------------------------
# f. Dependence-aware uncertainty heatmap
# ---------------------------------------------------------
def plot_f():
    rows = read_csv("panel_f_ratio.csv")
    methods = [r["Method"] for r in rows]
    columns = ["Protein", "Enzyme–substrate pair", "Reaction", "Reference", "Label assignment"]
    matrix = np.array([[float(r[c]) for c in columns] for r in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(7.3, 5.8))
    cmap = LinearSegmentedColormap.from_list(
        "ci_ratio_blues", ["#F7FBFF", "#6BAED6", "#08306B"]
    )
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0.9, vmax=3.35)

    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(["Protein", "Enzyme–substrate\npair", "Reaction", "Reference", "Label\nassignment"])
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            color = "white" if val >= 1.55 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=8.4, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.035)
    cbar.set_label("CI width ratio\n(cluster bootstrap / row bootstrap)")
    panel_label(ax, "f")
    fig.tight_layout()
    save_all(fig, "Figure2f_cluster_bootstrap_CI_ratio")

if __name__ == "__main__":
    plot_a()
    plot_b()
    plot_c()
    plot_d()
    plot_e()
    plot_f()
    print(f"Saved Figure 2 panels to: {OUT}")
