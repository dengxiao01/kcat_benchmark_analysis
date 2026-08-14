#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

METHODS = ["KcatNet", "CataPro", "TurNuP", "PMAK", "KinForm-L", "CatPred"]
METHOD_COLORS = dict(zip(METHODS, plt.get_cmap("tab10").colors[:len(METHODS)]))
PROX_METHODS = [
    "KcatNet","CataPro","PreTKcat","UniKP","DLKcat",
    "PMAK","KinForm-L","CatPred","DEKP-public-retrained"
]
PROX_COLORS = dict(zip(PROX_METHODS, plt.get_cmap("tab10").colors[:len(PROX_METHODS)]))

def read_csv(name):
    with open(DATA / name, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def fnum(x):
    if x is None or x == "":
        return np.nan
    return float(x)

def save_all(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)

def panel_a():
    rows = read_csv("figure4a_heatmap.csv")
    row_labels = [r["method"] for r in rows]
    col_labels = [
        "E. coli", "S. cerevisiae",
        "BRENDA only", "SABIO-RK only", "BRENDA +\nSABIO-RK",
        "Other\nreactant", "Currency/\ncofactor"
    ]
    groups = [("Species", 0, 1), ("Evidence source", 2, 4), ("Substrate role", 5, 6)]
    mat = np.array([
        [
            float(r["species_ecoli"]), float(r["species_yeast"]),
            float(r["source_BRENDA_only"]), float(r["source_SABIO_RK_only"]), float(r["source_BRENDA_plus_SABIO_RK"]),
            float(r["role_other_reactant"]), float(r["role_currency_cofactor"])
        ]
        for r in rows
    ])

    fig = plt.figure(figsize=(8.2, 4.8))
    ax = fig.add_axes([0.12, 0.18, 0.76, 0.64])
    im = ax.imshow(mat, aspect="auto", cmap="coolwarm", vmin=0.55, vmax=0.95)

    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(row_labels)

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
    for label in ax.get_xticklabels():
        label.set_fontsize(8)
        label.set_fontweight("bold")
    for label in ax.get_yticklabels():
        label.set_fontsize(9)

    # gridlines between cells
    ax.set_xticks(np.arange(-.5, mat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, mat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    # numbers
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            text_color = "white" if (val < 0.60 or val > 0.87) else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center", color=text_color, fontsize=8.5)

    # group headers with bracket lines
    y_line = -1.45
    y_text = -1.80
    for title, start, end in groups:
        x0 = start - 0.42
        x1 = end + 0.42
        ax.plot([x0, x1], [y_line, y_line], color="black", lw=1.0, clip_on=False)
        ax.plot([x0, x0], [y_line, y_line+0.12], color="black", lw=1.0, clip_on=False)
        ax.plot([x1, x1], [y_line, y_line+0.12], color="black", lw=1.0, clip_on=False)
        ax.text((start+end)/2, y_text, title, ha="center", va="center", fontsize=10, fontweight="bold", clip_on=False)

    cax = fig.add_axes([0.90, 0.19, 0.02, 0.62])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"MAE in log$_{10}$(kcat)")

    fig.text(0.02, 0.965, "a", fontsize=14, fontweight="bold", va="top")
    save_all(fig, "Figure4a_heatmap")

def panel_b():
    rows = read_csv("figure4b_contrasts.csv")
    y = np.arange(len(rows))
    fig = plt.figure(figsize=(7.4, 4.5))
    ax1 = fig.add_axes([0.12, 0.18, 0.38, 0.67])
    ax2 = fig.add_axes([0.60, 0.18, 0.35, 0.67])

    for i, r in enumerate(rows):
        m = r["method"]
        c = METHOD_COLORS[m]
        v1 = float(r["source_contrast_BRENDA_minus_SABIO"])
        v2 = float(r["substrate_contrast_other_minus_currency"])
        ax1.hlines(i, min(0, v1), max(0, v1), linewidth=2, color=c)
        ax1.scatter(v1, i, s=32, color=c, zorder=3)
        ax2.hlines(i, min(0, v2), max(0, v2), linewidth=2, color=c)
        ax2.scatter(v2, i, s=32, color=c, zorder=3)
        ax1.annotate(f"{v1:+.3f}", (v1, i), xytext=(5 if v1 >= 0 else -5, 0), textcoords="offset points",
                     ha="left" if v1 >= 0 else "right", va="center", fontsize=7.5)
        ax2.annotate(f"{v2:+.3f}", (v2, i), xytext=(5 if v2 >= 0 else -5, 0), textcoords="offset points",
                     ha="left" if v2 >= 0 else "right", va="center", fontsize=7.5)

    for ax in [ax1, ax2]:
        ax.axvline(0, linewidth=0.8, color="black")
        ax.set_yticks(y)
        ax.set_yticklabels(METHODS)
        ax.invert_yaxis()
        ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax1.set_xlim(-0.11, 0.065)
    ax2.set_xlim(-0.05, 0.115)
    ax1.set_title("Source contrast")
    ax2.set_title("Substrate-role contrast")
    ax1.set_xlabel("BRENDA only − SABIO-RK only MAE")
    ax2.set_xlabel("Other reactant − currency/cofactor MAE")
    fig.text(0.53, 0.92, "Carrier-linked substrates are excluded from the main contrast because n is small.",
             ha="center", va="center", fontsize=8, style="italic")
    fig.text(0.02, 0.965, "b", fontsize=14, fontweight="bold", va="top")
    save_all(fig, "Figure4b_contrasts")

def panel_c():
    rows = read_csv("figure4c_label_provenance.csv")
    d = {(r["section"], r["metric"]): r["value"] for r in rows}
    fig = plt.figure(figsize=(7.4, 4.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.05, 0.26, 0.27, 0.60, "Within-record dispersion"),
        (0.365, 0.26, 0.27, 0.60, "Cross-database difference"),
        (0.68, 0.26, 0.27, 0.60, "Source provenance audit"),
    ]
    for x, y, w, h, title in boxes:
        patch = FancyBboxPatch((x,y), w,h, boxstyle="round,pad=0.012,rounding_size=0.015",
                               linewidth=1.0, facecolor="white")
        ax.add_patch(patch)
        ax.text(x+w/2, y+h-0.07, title, ha="center", va="top", fontsize=10, fontweight="bold")

    ax.text(0.185, 0.69,
            f'{int(float(d[("within_record_dispersion","records_with_ge_2_measurements")]))} records with ≥2\nselected measurements',
            ha="center", va="center", fontsize=9)
    ax.text(0.185, 0.52, f'Median = {float(d[("within_record_dispersion","median_dispersion_log10")]):.3f} log10',
            ha="center", fontsize=9)
    ax.text(0.185, 0.42, f'Mean = {float(d[("within_record_dispersion","mean_dispersion_log10")]):.3f} log10',
            ha="center", fontsize=9)
    ax.text(0.50, 0.69,
            f'{int(float(d[("cross_database_difference","records_supported_by_both_databases")]))} records supported by\nboth BRENDA and SABIO-RK',
            ha="center", va="center", fontsize=9)
    ax.text(0.50, 0.50,
            f'Median |difference| =\n{float(d[("cross_database_difference","median_abs_database_difference_log10")]):.3f} log10',
            ha="center", va="center", fontsize=9)

    br_raw = int(float(d[("BRENDA_provenance","raw_candidates")]))
    br_exc = int(float(d[("BRENDA_provenance","excluded_by_mutation_variant_screen")]))
    br_ret = int(float(d[("BRENDA_provenance","retained_after_screen")]))
    sa_raw = int(float(d[("SABIO_RK_provenance","raw_candidates")]))

    ax.text(0.815, 0.71, f"BRENDA candidates: {br_raw:,}", ha="center", fontsize=9, fontweight="bold")
    ax.text(0.815, 0.61, f"Excluded by mutation/variant screen:\n{br_exc:,}", ha="center", fontsize=8.5)
    ax.annotate("", xy=(0.815,0.49), xytext=(0.815,0.56), arrowprops=dict(arrowstyle="->", lw=1))
    ax.text(0.815, 0.45, f"Retained after screen: {br_ret:,}", ha="center", fontsize=9, fontweight="bold")
    ax.text(0.815, 0.34, f"SABIO-RK candidates: {sa_raw:,}\nEquivalent mutation-status field unavailable",
            ha="center", fontsize=8.5)

    ax.text(0.5, 0.10,
            "Record-level residuals reflect both model limitations and uncertainty from\n"
            "experimental aggregation and asymmetric source annotation.",
            ha="center", va="center", fontsize=8.5)
    fig.text(0.02, 0.965, "c", fontsize=14, fontweight="bold", va="top")
    save_all(fig, "Figure4c_provenance")

def panel_d():
    rows = read_csv("figure4d_training_proximity.csv")
    y = np.arange(len(rows))
    fig = plt.figure(figsize=(7.6, 5.0))
    ax = fig.add_axes([0.20, 0.15, 0.75, 0.72])
    markers = {"exact": "o", "near": "^", "none": "s"}

    for i, r in enumerate(rows):
        m = r["method"]
        c = PROX_COLORS[m]
        pts = []
        exact_n = fnum(r["exact_pair_n"])
        exact_mae = fnum(r["exact_pair_MAE_log10"])
        near_n = fnum(r["joint_near_neighbor_n"])
        near_mae = fnum(r["joint_near_neighbor_MAE_log10"])
        none_n = fnum(r["no_joint_neighbor_n"])
        none_mae = fnum(r["no_joint_neighbor_MAE_log10"])

        if not math.isnan(exact_mae) and exact_n > 0:
            pts.append(("exact", exact_mae, int(exact_n)))
        if not math.isnan(near_mae) and near_n > 0:
            pts.append(("near", near_mae, int(near_n)))
        if not math.isnan(none_mae) and none_n > 0:
            pts.append(("none", none_mae, int(none_n)))

        xs = [p[1] for p in pts]
        if len(xs) > 1:
            ax.plot(xs, [i]*len(xs), linewidth=1.5, color=c, alpha=0.8)

        for j, (kind, x, n) in enumerate(pts):
            ax.scatter(x, i, s=34, marker=markers[kind], color=c, zorder=3)
            offset_y = 8 if j % 2 == 0 else -13
            va = "bottom" if offset_y > 0 else "top"
            ax.annotate(f"{kind}\n(n={n})", (x, i), xytext=(4, offset_y), textcoords="offset points",
                        ha="left", va=va, fontsize=7)

    ax.set_yticks(y)
    ax.set_yticklabels([r["method"] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0.25, 1.18)
    ax.set_xlabel(r"MAE (log$_{10}$ $k_{\mathrm{cat}}$)")
    ax.set_title("Prediction performance across public training-corpus proximity strata")
    ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Line2D([0],[0], marker=markers["exact"], linestyle="none", label="exact"),
        Line2D([0],[0], marker=markers["near"], linestyle="none", label="joint near neighbour"),
        Line2D([0],[0], marker=markers["none"], linestyle="none", label="no joint neighbour"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right")
    fig.text(0.20, 0.055,
             "Exact = exact sequence–parent pair; near = joint sequence/chemical neighbour; none = no joint neighbour.",
             fontsize=7.5)
    fig.text(0.02, 0.965, "d", fontsize=14, fontweight="bold", va="top")
    save_all(fig, "Figure4d_training_proximity")

def panel_e():
    rows = read_csv("figure4e_pretkcat_sensitivity.csv")
    labels = [r["variant"] for r in rows]
    x = np.arange(len(rows))
    mae = [float(r["MAE_log10"]) for r in rows]
    rho = [float(r["Spearman_log10"]) for r in rows]
    train_n = [int(float(r["fitted_training_n"])) for r in rows]

    fig = plt.figure(figsize=(6.4, 4.5))
    ax1 = fig.add_axes([0.13, 0.18, 0.70, 0.68])
    ax2 = ax1.twinx()
    line1, = ax1.plot(x, mae, marker="o", linewidth=1.8, label="MAE")
    line2, = ax2.plot(x, rho, marker="s", linestyle="--", linewidth=1.6, label="Spearman")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel(r"MAE (log$_{10}$ $k_{\mathrm{cat}}$)")
    ax2.set_ylabel("Spearman correlation")
    ax1.set_ylim(0.84, 0.98)
    ax2.set_ylim(0.34, 0.42)
    ax1.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.5)
    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    for i in range(len(rows)):
        ax1.annotate(f"{mae[i]:.3f}\ntrain n={train_n[i]:,}", (x[i], mae[i]), xytext=(0, 8),
                     textcoords="offset points", ha="center", va="bottom", fontsize=7.5)
        ax2.annotate(f"ρ={rho[i]:.3f}", (x[i], rho[i]), xytext=(0, -11),
                     textcoords="offset points", ha="center", va="top", fontsize=7.5)

    ax1.axvline(1, linestyle=":", linewidth=1.0, alpha=0.6)
    ax1.set_title("PreTKcat sensitivity to progressively stricter training-corpus exclusion")
    ax1.legend([line1, line2], ["MAE", "Spearman correlation"], frameon=False, loc="lower right")
    fig.text(0.51, 0.065, "Exact-excluded is the primary PreTKcat result.", ha="center", fontsize=8)
    fig.text(0.02, 0.965, "e", fontsize=14, fontweight="bold", va="top")
    save_all(fig, "Figure4e_PreTKcat_sensitivity")

if __name__ == "__main__":
    panel_a()
    panel_b()
    panel_c()
    panel_d()
    panel_e()
    print(f"Figure 4 panels written to: {OUT}")
