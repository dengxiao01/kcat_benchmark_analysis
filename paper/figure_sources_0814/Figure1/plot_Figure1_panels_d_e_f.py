#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reproduce Figure 1 panels d, e and f.

Inputs
------
panel_d_matching_evidence.csv
panel_e_database_support.csv
panel_f_metadata_completeness.csv

Outputs
-------
Figure1d_matching_evidence.png/.pdf/.svg
Figure1e_database_support.png/.pdf/.svg
Figure1f_metadata_completeness.png/.pdf/.svg
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

BASE = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 0.9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def read_csv_dict(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def save_formats(fig, stem):
    fig.savefig(BASE / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(BASE / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(BASE / f"{stem}.svg", bbox_inches="tight")


# ------------------------------------------------------------------
# Figure 1d — Experimental matching evidence by species
# ------------------------------------------------------------------
def plot_panel_d():
    rows = read_csv_dict(BASE / "panel_d_matching_evidence.csv")
    categories = ["Accession-supported", "EC + substrate ID", "Name-supported"]
    plot_species = ["S. cerevisiae", "E. coli"]

    pct = {sp: {cat: 0.0 for cat in categories} for sp in plot_species}
    count = {sp: {cat: 0 for cat in categories} for sp in plot_species}

    for r in rows:
        sp = r["Species"]
        cat = r["Evidence"]
        pct[sp][cat] = float(r["Percent_of_species"])
        count[sp][cat] = int(r["Records"])

    fig, ax = plt.subplots(figsize=(6.4, 3.5))
    y = np.arange(len(plot_species))
    left = np.zeros(len(plot_species))

    for cat in categories:
        vals = np.array([pct[sp][cat] for sp in plot_species])
        bars = ax.barh(
            y, vals, left=left, height=0.52,
            edgecolor="white", linewidth=0.9, label=cat
        )

        for i, (bar, value) in enumerate(zip(bars, vals)):
            sp = plot_species[i]
            n = count[sp][cat]

            if value >= 8:
                ax.text(
                    left[i] + value / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{n} ({value:.1f}%)",
                    ha="center", va="center",
                    color="white", fontweight="bold", fontsize=9
                )
            elif value > 0:
                ax.annotate(
                    f"{n} ({value:.1f}%)",
                    xy=(left[i] + value,
                        bar.get_y() + bar.get_height() / 2),
                    xytext=(7, 0),
                    textcoords="offset points",
                    ha="left", va="center", fontsize=8.5,
                    arrowprops=dict(arrowstyle="-", linewidth=0.7)
                )
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels([
        r"$\it{S.\ cerevisiae}$",
        r"$\it{E.\ coli}$"
    ])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_xlabel("Proportion of records")
    ax.set_title(
        "Experimental matching evidence by species",
        fontsize=11, fontweight="bold", pad=8
    )
    ax.legend(
        frameon=False, ncol=3, loc="upper center",
        bbox_to_anchor=(0.5, -0.22), fontsize=8.6
    )
    ax.tick_params(axis="y", length=0)

    fig.subplots_adjust(
        bottom=0.28, left=0.23, right=0.97, top=0.84
    )
    save_formats(fig, "Figure1d_matching_evidence")
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 1e — Experimental source-database support by species
# ------------------------------------------------------------------
def plot_panel_e():
    rows = read_csv_dict(BASE / "panel_e_database_support.csv")
    categories = ["BRENDA only", "Both databases", "SABIO-RK only"]
    plot_species = ["S. cerevisiae", "E. coli"]

    pct = {sp: {cat: 0.0 for cat in categories} for sp in plot_species}
    count = {sp: {cat: 0 for cat in categories} for sp in plot_species}

    for r in rows:
        sp = r["Species"]
        cat = r["Database_support"]
        pct[sp][cat] = float(r["Percent_of_species"])
        count[sp][cat] = int(r["Records"])

    fig, ax = plt.subplots(figsize=(6.4, 3.5))
    y = np.arange(len(plot_species))
    left = np.zeros(len(plot_species))

    for cat in categories:
        vals = np.array([pct[sp][cat] for sp in plot_species])
        bars = ax.barh(
            y, vals, left=left, height=0.52,
            edgecolor="white", linewidth=0.9, label=cat
        )

        for i, (bar, value) in enumerate(zip(bars, vals)):
            sp = plot_species[i]
            n = count[sp][cat]
            ax.text(
                left[i] + value / 2,
                bar.get_y() + bar.get_height() / 2,
                f"{n}\n({value:.1f}%)",
                ha="center", va="center",
                color="white", fontweight="bold",
                fontsize=8.4 if value >= 7 else 7.4,
                linespacing=1.1
            )
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels([
        r"$\it{S.\ cerevisiae}$",
        r"$\it{E.\ coli}$"
    ])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_xlabel("Proportion of records")
    ax.set_title(
        "Experimental source-database support by species",
        fontsize=11, fontweight="bold", pad=8
    )
    ax.legend(
        frameon=False, ncol=3, loc="upper center",
        bbox_to_anchor=(0.5, -0.22), fontsize=8.6
    )
    ax.tick_params(axis="y", length=0)

    fig.subplots_adjust(
        bottom=0.28, left=0.23, right=0.97, top=0.84
    )
    save_formats(fig, "Figure1e_database_support")
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 1f — Experimental-condition metadata completeness
# ------------------------------------------------------------------
def plot_panel_f():
    rows = read_csv_dict(BASE / "panel_f_metadata_completeness.csv")

    species = [r["Species"] for r in rows]
    totals = np.array([int(r["Total_records"]) for r in rows])
    ph_n = np.array([int(r["pH_available_n"]) for r in rows])
    temp_n = np.array([int(r["Temperature_available_n"]) for r in rows])
    ph = np.array([float(r["pH_available_percent"]) for r in rows])
    temp = np.array([
        float(r["Temperature_available_percent"]) for r in rows
    ])

    x = np.arange(len(species))
    width = 0.30

    fig, ax = plt.subplots(figsize=(5.2, 3.8))

    bars1 = ax.bar(
        x - width / 2, ph, width,
        label="pH available",
        edgecolor="white", linewidth=0.8
    )
    bars2 = ax.bar(
        x + width / 2, temp, width,
        label="Temperature available",
        edgecolor="white", linewidth=0.8
    )

    for bars, values, counts in [
        (bars1, ph, ph_n),
        (bars2, temp, temp_n)
    ]:
        for i, (bar, value, n) in enumerate(zip(bars, values, counts)):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 1.5,
                f"{value:.1f}%\n({n}/{totals[i]})",
                ha="center", va="bottom",
                fontsize=8.5, linespacing=1.1
            )

    ax.set_xticks(x)
    ax.set_xticklabels([
        r"$\it{E.\ coli}$",
        r"$\it{S.\ cerevisiae}$"
    ])
    ax.set_ylim(0, 100)
    ax.set_ylabel("Records with metadata (%)")
    ax.set_title(
        "Experimental-condition metadata completeness",
        fontsize=11, fontweight="bold", pad=8
    )
    ax.yaxis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    ax.set_axisbelow(True)
    ax.legend(
        frameon=False, ncol=2, loc="upper center",
        bbox_to_anchor=(0.5, -0.18), fontsize=8.6
    )

    fig.subplots_adjust(
        bottom=0.25, left=0.16, right=0.97, top=0.84
    )
    save_formats(fig, "Figure1f_metadata_completeness")
    plt.close(fig)


if __name__ == "__main__":
    plot_panel_d()
    plot_panel_e()
    plot_panel_f()
    print("Figure 1d–f generated successfully.")
