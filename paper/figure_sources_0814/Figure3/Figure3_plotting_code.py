#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reproduce Figure 3 panels a-e.

Required packages:
    pandas
    numpy
    matplotlib

Folder structure:
    Figure3_plotting_code.py
    Figure3_plotting_data/
        Panel_a.csv
        Panel_b.csv
        Panel_c.csv
        Panel_d.csv
        Panel_e.csv

Outputs:
    figure3_panels/Figure3a.png/.pdf
    ...
    figure3_panels/Figure3e.png/.pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

DATA_DIR = Path("Figure3_plotting_data")
OUTDIR = Path("figure3_panels")
OUTDIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

METHOD_COLORS = {
    "KcatNet": "#1f77b4",
    "TurNuP": "#ff7f0e",
    "PMAK": "#2ca02c",
    "CataPro": "#9467bd",
    "CatPred": "#d62728",
    "KinForm-L": "#17becf",
    "PreTKcat": "#8c564b",
    "UniKP": "#6b7c93",
    "SELFprot": "#7f7f7f",
    "DLKcat": "#9a9567",
    "DEKP-public-retrained": "#6f6f6f",
}

def save_panel(fig, stem):
    fig.savefig(OUTDIR / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)

def clean_axes(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, linestyle=(0, (1.5, 2.5)),
                linewidth=0.65, alpha=0.35)
        ax.set_axisbelow(True)

# ============================================================
# Figure 3a
# ============================================================
a = pd.read_csv(DATA_DIR / "Panel_a.csv")
a = a.sort_values("mae_log10_common_subset").reset_index(drop=True)

fig, ax = plt.subplots(figsize=(7.2, 5.0))
y = np.arange(len(a))

for i, row in a.iterrows():
    x = row["mae_log10_common_subset"]
    method = row["method"]
    ax.hlines(i, 0.60, x, color="#CFCFCF", linewidth=1.2, zorder=1)
    ax.scatter(x, i, s=48, color=METHOD_COLORS.get(method, "#777777"),
               zorder=3)
    ax.text(x + 0.008, i, f"{x:.3f}", va="center", fontsize=8.5)

ax.set_yticks(y)
ax.set_yticklabels(a["method"])
ax.invert_yaxis()
ax.set_xlim(0.60, 1.12)
ax.set_xlabel(r"MAE in $\log_{10}(k_{\mathrm{cat}})$")
ax.set_title(
    "a   Strict paired comparison on the 1,047-record reaction-common set",
    loc="left", fontweight="bold", pad=8
)
ax.spines["left"].set_visible(False)
clean_axes(ax, "x")
fig.tight_layout()
save_panel(fig, "Figure3a")

# ============================================================
# Figure 3b
# ============================================================
b = pd.read_csv(DATA_DIR / "Panel_b.csv")

scope_defs = [
    ("CatPred-accessible scope", 1156),
    ("KinForm-L-accessible scope", 729),
]

fig, axes = plt.subplots(2, 1, figsize=(7.0, 8.0), sharex=True)

for ax, (scope, total) in zip(axes, scope_defs):
    d = b[b["subset"] == scope].sort_values("mae_log10").reset_index(drop=True)
    y = np.arange(len(d))

    for i, row in d.iterrows():
        x = row["mae_log10"]
        method = row["method"]
        full = int(row["n"]) == int(row["subset_total_rows"])

        ax.hlines(i, 0.62, x, color="#CFCFCF", linewidth=1.1, zorder=1)
        ax.scatter(
            x, i, s=48,
            facecolors=METHOD_COLORS[method] if full else "white",
            edgecolors=METHOD_COLORS[method],
            linewidths=1.6,
            zorder=3,
        )
        ax.text(x + 0.008, i, f"{x:.3f}", va="center", fontsize=8.3)
        ax.text(
            1.01, i,
            f"n={int(row['n'])}/{int(row['subset_total_rows'])}",
            transform=ax.get_yaxis_transform(),
            va="center", ha="left", fontsize=8.0
        )

    ax.set_yticks(y)
    ax.set_yticklabels(d["method"])
    ax.invert_yaxis()
    ax.set_xlim(0.62, 0.90)
    ax.set_title(f"{scope} ({total:,} rows)", loc="left", fontsize=10)
    ax.spines["left"].set_visible(False)
    clean_axes(ax, "x")

axes[0].text(-0.13, 1.16, "b", transform=axes[0].transAxes,
             fontsize=14, fontweight="bold", va="top")
axes[0].text(
    -0.03, 1.16,
    "Available-case summaries for CatPred and KinForm-L scopes",
    transform=axes[0].transAxes, fontsize=11, fontweight="bold", va="top"
)
axes[-1].set_xlabel(r"MAE in $\log_{10}(k_{\mathrm{cat}})$")
fig.subplots_adjust(left=0.20, right=0.79, top=0.90,
                    bottom=0.08, hspace=0.24)
save_panel(fig, "Figure3b")

# ============================================================
# Figure 3c
# ============================================================
c = pd.read_csv(DATA_DIR / "Panel_c.csv")
c["comparison"] = c["method_a"] + " − " + c["method_b"]
order = ["KcatNet − TurNuP", "KcatNet − PMAK", "TurNuP − PMAK"]
c = c.set_index("comparison").loc[order].reset_index()

fig, ax = plt.subplots(figsize=(7.2, 4.2))
y = np.arange(len(c))
ax.axvline(0, color="black", linestyle="--", linewidth=0.9)

for i, row in c.iterrows():
    est = row["mae_difference_a_minus_b"]
    lo = row["bootstrap_ci_low_95"]
    hi = row["bootstrap_ci_high_95"]
    ax.hlines(i, lo, hi, color="#4D4D4D", linewidth=1.5)
    ax.scatter(est, i, s=46, color=METHOD_COLORS["KcatNet"], zorder=3)
    ax.text(hi + 0.0025, i, f"{est:+.4f}", va="center", fontsize=8.5)

ax.set_yticks(y)
ax.set_yticklabels(c["comparison"])
ax.invert_yaxis()
ax.set_xlim(-0.045, 0.045)
ax.set_xlabel("Difference in MAE (a − b)")
ax.set_title(
    "c   Paired row-bootstrap intervals on the 1,047-record common set",
    loc="left", fontweight="bold", pad=8
)
ax.spines["left"].set_visible(False)
clean_axes(ax, "x")
fig.tight_layout()
save_panel(fig, "Figure3c")

# ============================================================
# Figure 3d
# ============================================================
d = pd.read_csv(DATA_DIR / "Panel_d.csv")

row_order = [
    ("KcatNet", "TurNuP"),
    ("KcatNet", "PMAK"),
    ("TurNuP", "PMAK"),
]
col_order = ["protein", "pair", "reaction", "reference", "label_assignment"]
col_labels = [
    "Protein", "Seq–substrate\npair", "Reaction",
    "Reference", "Label\nassignment"
]

matrix = np.zeros((3, 5))
q_text = [[""] * 5 for _ in range(3)]

for i, (ma, mb) in enumerate(row_order):
    for j, cluster in enumerate(col_order):
        r = d[
            (d["method_a"] == ma)
            & (d["method_b"] == mb)
            & (d["cluster_type"] == cluster)
        ].iloc[0]
        matrix[i, j] = r["error_difference_a_minus_b"]
        sig = str(r["significant_bh_fdr_0.05"]).lower() == "true"
        q_text[i][j] = f"q={r['p_value_bh_global']:.3f}" + ("*" if sig else "")

cmap = LinearSegmentedColormap.from_list(
    "blue_white_orange", ["#2C7FB8", "#F7F7F7", "#F28E2B"]
)
norm = TwoSlopeNorm(vmin=-0.06, vcenter=0.0, vmax=0.06)

fig, ax = plt.subplots(figsize=(7.5, 4.5))
im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")

ax.set_xticks(np.arange(5))
ax.set_xticklabels(col_labels)
ax.set_yticks(np.arange(3))
ax.set_yticklabels([
    "KcatNet vs TurNuP",
    "KcatNet vs PMAK",
    "TurNuP vs PMAK"
])

for i in range(3):
    for j in range(5):
        ax.text(j, i, q_text[i][j], ha="center", va="center", fontsize=8.5)

ax.set_title(
    "d   Cluster-aware Wilcoxon sensitivity across dependence definitions",
    loc="left", fontweight="bold", pad=8
)

cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
cbar.set_label("Cluster-mean error difference\n(method_a − method_b)")

fig.text(
    0.17, 0.01,
    "Blue: method_a lower error; orange: method_b lower error. "
    "Values are BH-adjusted q values; * q < 0.05.",
    fontsize=7.5
)
fig.tight_layout(rect=(0, 0.04, 1, 1))
save_panel(fig, "Figure3d")

# ============================================================
# Figure 3e
# ============================================================
e = pd.read_csv(DATA_DIR / "Panel_e.csv")

unit_order = [
    "Benchmark rows",
    "Unique seq–substrate pairs",
    "Unique label assignments",
]
methods = ["KcatNet", "TurNuP", "PMAK"]
x = np.arange(3)

fig, ax = plt.subplots(figsize=(7.2, 4.8))

for method in methods:
    dm = e[e["method"] == method].copy()
    dm["analysis_unit"] = pd.Categorical(
        dm["analysis_unit"], categories=unit_order, ordered=True
    )
    dm = dm.sort_values("analysis_unit")
    vals = dm["mae_log10"].to_numpy()

    ax.plot(x, vals, marker="o", markersize=6, linewidth=1.7,
            color=METHOD_COLORS[method], label=method)
    for xi, yi in zip(x, vals):
        ax.text(xi, yi + 0.0045, f"{yi:.3f}",
                ha="center", fontsize=8.2)

ax.set_xticks(x)
ax.set_xticklabels([
    "Benchmark rows\n(1,047)",
    "Unique seq–substrate pairs\n(796)",
    "Unique label assignments\n(761)",
])
ax.set_ylabel(r"MAE in $\log_{10}(k_{\mathrm{cat}})$")
ax.set_ylim(0.68, 0.78)
ax.set_title(
    "e   Method ordering changes with the statistical unit",
    loc="left", fontweight="bold", pad=8
)
ax.legend(frameon=False, loc="upper right")
clean_axes(ax, "y")
fig.tight_layout()
save_panel(fig, "Figure3e")

print("Finished. Output folder:", OUTDIR.resolve())
