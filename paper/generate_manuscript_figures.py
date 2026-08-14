#!/usr/bin/env python3
"""Rebuild the public manuscript-matched Figures 1-4 and standalone panels.

Figure 1a/b use versioned high-resolution base panels, with only benchmark
values replaced. Figures 3 and 4 follow the layouts accepted for manuscript
artifact revision 1.2.0-r5.
"""

from __future__ import annotations

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib

matplotlib.use("Agg")
import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from PIL import Image, ImageDraw, ImageFont

import generate_manuscript_figures_0806_style_r4 as r4


BASE = Path(__file__).resolve().parent.parent
TABLE_DIR = BASE / "paper" / "tables_v1.2.0"
OUTPUT_ROOT = BASE / "paper"
OUTPUT_DIR = OUTPUT_ROOT / "figures"
PANEL_DIR = OUTPUT_DIR / "panels"
SOURCE_DIR = OUTPUT_DIR / "source_data"
PACKAGE_DIR = OUTPUT_ROOT / "figure_packages"
REFERENCE_A = BASE / "paper" / "figure_references" / "Figure1a_0806_base.png"
REFERENCE_B = BASE / "paper" / "figure_references" / "Figure1b_0806_base.png"
TABLE0 = TABLE_DIR / "Table0.csv"
RECORD_AUDIT = TABLE_DIR / "Record_audit.csv"
HELPER_SCRIPT = Path(r4.__file__).resolve()


OLD_METHOD_COLORS = {
    "KcatNet": "#0072B2",
    "TurNuP": "#E69F00",
    "PMAK": "#009E73",
    "CataPro": "#CC79A7",
    "CatPred": "#8C564B",
    "KinForm-L": "#D55E00",
}
METHODS_4 = ["KcatNet", "CataPro", "TurNuP", "PMAK", "CatPred", "KinForm-L"]
CLUSTER_TYPES = [
    ("protein", "Sequence"),
    ("pair", "Seq–substrate pair"),
    ("reaction", "Reaction"),
    ("reference", "Reference"),
    ("label_assignment", "Label assignment"),
]
PAIR_ORDER = [("KcatNet", "TurNuP"), ("KcatNet", "PMAK"), ("TurNuP", "PMAK")]


def read_table(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLE_DIR / f"{name}.csv")


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_composite(fig: plt.Figure, number: int) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / f"Figure{number}.png", dpi=350, bbox_inches="tight",
                pad_inches=0.05, facecolor="white")
    fig.savefig(OUTPUT_DIR / f"Figure{number}.pdf", bbox_inches="tight",
                pad_inches=0.05, facecolor="white")
    plt.close(fig)


def save_axes_panel(
    fig: plt.Figure,
    axes: list[plt.Axes],
    stem: str,
    extra_artists: list[object] | None = None,
    hidden_artists: list[object] | None = None,
) -> None:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    hidden = hidden_artists or []
    visibility = [artist.get_visible() for artist in hidden]
    for artist in hidden:
        artist.set_visible(False)
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [axis.get_tightbbox(renderer) for axis in axes]
        boxes.extend(artist.get_window_extent(renderer) for artist in (extra_artists or []))
        bbox = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted()).expanded(1.025, 1.08)
        fig.savefig(PANEL_DIR / f"{stem}.png", dpi=350, bbox_inches=bbox,
                    pad_inches=0.04, facecolor="white")
        fig.savefig(PANEL_DIR / f"{stem}.pdf", bbox_inches=bbox,
                    pad_inches=0.04, facecolor="white")
    finally:
        for artist, visible in zip(hidden, visibility):
            artist.set_visible(visible)


def panel_label(fig: plt.Figure, x: float, y: float, label: str) -> object:
    return fig.text(x, y, label, fontsize=15, fontweight="bold", va="top", ha="left")


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "LiberationSans-Bold.ttf" if bold else "LiberationSans-Regular.ttf"
    return ImageFont.truetype(f"/usr/share/fonts/liberation-sans/{name}", size=size)


def fill_box(image: Image.Image, box: tuple[int, int, int, int], color=(250, 252, 253)) -> None:
    ImageDraw.Draw(image).rectangle(box, fill=color)


def restore_card_background(
    image: Image.Image,
    box: tuple[int, int, int, int],
    feather: int = 1,
) -> Image.Image:
    """Replace old text with a softly blended model of the card background."""
    x1, y1, x2, y2 = box
    array = np.asarray(image).astype(np.float32)
    height = y2 - y1
    width = x2 - x1
    margin = 48
    sx1, sx2 = max(0, x1 - margin), min(array.shape[1], x2 + margin)
    sy1, sy2 = max(0, y1 - margin), min(array.shape[0], y2 + margin)
    sample_pixels = array[sy1:sy2, sx1:sx2]
    sample_y, sample_x = np.mgrid[sy1:sy2, sx1:sx2]
    frame = ((sample_x < x1) | (sample_x >= x2)
             | (sample_y < y1) | (sample_y >= y2))
    gray = sample_pixels.mean(axis=2)
    chroma = sample_pixels.max(axis=2) - sample_pixels.min(axis=2)
    valid = frame & (gray > 228) & (chroma < 38)
    if int(valid.sum()) < 200:
        raise ValueError(f"Insufficient clean background pixels around {box}")

    center_x = (x1 + x2 - 1) / 2
    center_y = (y1 + y2 - 1) / 2
    scale_x = max(width, 1)
    scale_y = max(height, 1)

    def design(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
        xn = (xs - center_x) / scale_x
        yn = (ys - center_y) / scale_y
        return np.column_stack([
            np.ones_like(xn), xn, yn, xn * yn, xn ** 2, yn ** 2,
        ])

    fit_x = sample_x[valid].astype(np.float64)
    fit_y = sample_y[valid].astype(np.float64)
    fit_design = design(fit_x, fit_y)
    target_y, target_x = np.mgrid[y1:y2, x1:x2]
    target_design = design(target_x.ravel(), target_y.ravel())
    replacement_channels = []
    for channel in range(3):
        coefficients, *_ = np.linalg.lstsq(
            fit_design,
            sample_pixels[:, :, channel][valid].astype(np.float64),
            rcond=None,
        )
        replacement_channels.append(
            target_design.dot(coefficients).reshape(height, width)
        )
    replacement = np.stack(replacement_channels, axis=2).astype(np.float32)

    x_distance = np.minimum(np.arange(width) + 1, np.arange(width, 0, -1))
    y_distance = np.minimum(np.arange(height) + 1, np.arange(height, 0, -1))
    alpha = np.minimum.outer(y_distance, x_distance).astype(np.float32)
    alpha = np.clip(alpha / feather, 0.0, 1.0)[:, :, None]
    original = array[y1:y2, x1:x2]
    array[y1:y2, x1:x2] = original * (1.0 - alpha) + replacement * alpha
    return Image.fromarray(np.clip(array, 0, 255).astype(np.uint8))


def restore_complete_card_content(
    image: Image.Image,
    box: tuple[int, int, int, int],
    bottom_blend: int = 18,
) -> Image.Image:
    """Rebuild a complete card-content area without introducing an inner rectangle."""
    x1, y1, x2, y2 = box
    source = np.asarray(image).astype(np.float32)
    crop = source[y1:y2, x1:x2]
    gray = crop.mean(axis=2)
    chroma = crop.max(axis=2) - crop.min(axis=2)
    clean = (gray > 235) & (chroma < 32)
    fallback = np.median(crop[clean], axis=0)

    row_profile = np.empty((crop.shape[0], 3), dtype=np.float32)
    for index in range(crop.shape[0]):
        pixels = crop[index][clean[index]]
        row_profile[index] = np.median(pixels, axis=0) if len(pixels) >= 40 else fallback
    row_profile = cv2.GaussianBlur(row_profile[:, None, :], (1, 31), 0)[:, 0, :]

    column_profile = np.empty((crop.shape[1], 3), dtype=np.float32)
    for index in range(crop.shape[1]):
        pixels = crop[:, index][clean[:, index]]
        column_profile[index] = np.median(pixels, axis=0) if len(pixels) >= 40 else fallback
    column_profile = cv2.GaussianBlur(column_profile[None, :, :], (31, 1), 0)[0]

    replacement = row_profile[:, None, :] + column_profile[None, :, :] - fallback
    alpha = np.ones(crop.shape[:2], dtype=np.float32)
    if bottom_blend:
        alpha[-bottom_blend:] = np.linspace(1.0, 0.0, bottom_blend)[:, None]
    result = source.copy()
    result[y1:y2, x1:x2] = (
        crop * (1.0 - alpha[:, :, None])
        + replacement * alpha[:, :, None]
    )
    return Image.fromarray(np.clip(result, 0, 255).astype(np.uint8))


def centered_multiline(
    image: Image.Image,
    center: tuple[int, int],
    text: str,
    text_font: ImageFont.FreeTypeFont,
    fill: str,
    spacing: int = 5,
) -> None:
    ImageDraw.Draw(image).multiline_text(
        center,
        text,
        font=text_font,
        fill=fill,
        anchor="mm",
        align="center",
        spacing=spacing,
    )


def exact_figure1_reference_panels() -> tuple[Image.Image, Image.Image]:
    """Load and update the versioned high-resolution Figure 1 base panels."""
    if not REFERENCE_A.exists() or not REFERENCE_B.exists():
        raise FileNotFoundError("The versioned Figure 1 base panels are required")
    panel_a = Image.open(REFERENCE_A).convert("RGB")
    panel_b = Image.open(REFERENCE_B).convert("RGB")

    # Panel a: preserve cards, icons, headings, colors and spacing; replace only
    # the three data blocks whose values changed in benchmark v1.2.0.
    # Rebuild the complete lower content area of card 2. Replacing only the
    # old text bounding box produces a visible light rectangle on this card.
    panel_a = restore_complete_card_content(panel_a, (414, 538, 718, 870))
    ImageDraw.Draw(panel_a).line((405, 535, 690, 535), fill="#A8D6CC", width=2)
    centered_multiline(
        panel_a,
        (548, 655),
        "16,613 model-linked\nenzyme–reaction–\nmetabolite associations",
        font(25),
        "#222222",
        spacing=6,
    )
    centered_multiline(panel_a, (548, 810), "8,797 E. coli | 7,816 yeast",
                       font(23, bold=True), "#148A78")

    panel_a = restore_card_background(panel_a, (1185, 500, 1480, 720))
    ImageDraw.Draw(panel_a).line((1205, 560, 1460, 560), fill="#B8A8CF", width=2)
    centered_multiline(panel_a, (1333, 520), "BRENDA + SABIO-RK",
                       font(24), "#222222")
    centered_multiline(panel_a, (1333, 630), "1,354 matched\nassociations",
                       font(27, bold=True), "#6543A8", spacing=5)

    panel_a = restore_card_background(panel_a, (1950, 535, 2275, 904))
    ImageDraw.Draw(panel_a).line((1970, 535, 2255, 535), fill="#B7C9E8", width=2)
    centered_multiline(panel_a, (2112, 595), "1,246 benchmark records",
                       font(25, bold=True), "#1647B3")
    centered_multiline(panel_a, (2112, 710), "781 E. coli | 465 S. cerevisiae",
                       font(22, bold=True), "#1647B3")
    centered_multiline(panel_a, (2112, 830), "12 methods | 6\ninference regimes",
                       font(24, bold=True), "#16276B", spacing=4)

    # Panel b: preserve the original metric-card artwork and update only values.
    blue = "#174F92"
    teal = "#158E82"
    card_text = [
        ((175, 45, 410, 226), (282, 101), "1,246", 72, blue,
         (282, 181), "benchmark\nrecords"),
        ((610, 45, 820, 226), (702, 101), "773", 72, blue,
         (702, 181), "distinct\nreactions"),
        ((1020, 45, 1240, 226), (1121, 101), "514", 72, blue,
         (1121, 181), "unique protein\nsequences"),
        ((185, 282, 405, 470), (285, 369), "390", 68, blue,
         (285, 435), "EC\nannotations"),
    ]
    for box, value_center, value, size, color, note_center, note in card_text:
        panel_b = restore_card_background(panel_b, box)
        centered_multiline(panel_b, value_center, value, font(size, bold=True), color)
        centered_multiline(panel_b, note_center, note, font(25), "#222222", spacing=2)

    panel_b = restore_card_background(panel_b, (1020, 370, 1238, 480))
    centered_multiline(panel_b, (1116, 430), "1.279", font(64, bold=True), blue)

    panel_b = restore_card_background(panel_b, (105, 575, 260, 732))
    centered_multiline(panel_b, (189, 635), "781", font(59, bold=True), teal)
    centered_multiline(panel_b, (189, 694), "E. coli", font(27, bold=True), teal)

    ImageDraw.Draw(panel_b).ellipse((426, 616, 509, 720), fill=(253, 253, 252))
    centered_multiline(panel_b, (467, 668), "Total\n1246", font(27, bold=True),
                       "#222222", spacing=1)

    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    panel_a.save(PANEL_DIR / "Figure1a.png", dpi=(350, 350))
    panel_a.save(PANEL_DIR / "Figure1a.pdf", resolution=350)
    panel_b.save(PANEL_DIR / "Figure1b.png", dpi=(350, 350))
    panel_b.save(PANEL_DIR / "Figure1b.pdf", resolution=350)
    return panel_a, panel_b


def draw_matching_panel(ax: plt.Axes) -> None:
    data = read_table("S3_Matching_levels")
    species = ["ecoli", "yeast"]
    labels = ["E. coli", "S. cerevisiae"]
    levels = [
        ("Species + EC + UniProt + substrate ID", "UniProt + substrate ID", "#0B55B5", 1.0),
        ("Species + EC + substrate ID", "Substrate ID", "#4B84D5", 1.0),
        ("Species + EC + substrate name", "Substrate name", "#B8CFF1", 1.0),
    ]
    x = np.arange(2)
    bottoms = np.zeros(2)
    totals = []
    for species_name in species:
        totals.append(int(data.loc[data["Species"].eq(species_name), "Records"].sum()))
    for level, legend, color, alpha in levels:
        values = np.asarray([
            int(data.loc[data["Species"].eq(species_name) & data["Matching level"].eq(level),
                         "Records"].sum())
            for species_name in species
        ])
        bars = ax.bar(x, values, 0.55, bottom=bottoms, color=color, alpha=alpha,
                      edgecolor="#0B55B5", linewidth=0.5, label=legend)
        for bar, value, bottom in zip(bars, values, bottoms):
            if value >= 20:
                ax.text(bar.get_x() + bar.get_width() / 2, bottom + value / 2, f"{value}",
                        ha="center", va="center", color="white", fontsize=7.5,
                        fontweight="bold")
            elif value:
                ax.text(bar.get_x() + bar.get_width() / 2, bottom + value + 12, f"{value}",
                        ha="center", va="bottom", fontsize=7.3, fontweight="bold")
            elif level == "Species + EC + substrate name":
                ax.text(bar.get_x() + bar.get_width() / 2, bottom + 12, "0",
                        ha="center", va="bottom", fontsize=7.3, fontweight="bold")
        bottoms += values
    ax.set_ylim(0, 850)
    ax.set_ylabel("Number of benchmark records")
    ax.set_xticks(x, [f"{name}\n(total {total})" for name, total in zip(labels, totals)])
    ax.grid(axis="y", linestyle="--", linewidth=0.55, alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.58), fontsize=7.3)


def draw_source_panel(ax: plt.Axes) -> None:
    data = read_table("S4_Source_support")
    species = ["ecoli", "yeast"]
    labels = ["E. coli", "S. cerevisiae"]
    sources = [
        ("BRENDA", "BRENDA only", "#155CB5"),
        ("SABIO-RK", "SABIO-RK only", "#159789"),
        ("BRENDA;SABIO-RK", "Both databases", "#7650A7"),
    ]
    x = np.arange(2)
    bottoms = np.zeros(2)
    totals = []
    for species_name in species:
        totals.append(int(data.loc[data["Species"].eq(species_name), "Records"].sum()))
    for source, legend, color in sources:
        values = np.asarray([
            int(data.loc[data["Species"].eq(species_name) & data["Source database"].eq(source),
                         "Records"].sum())
            for species_name in species
        ])
        bars = ax.bar(x, values, 0.55, bottom=bottoms, color=color, edgecolor="white",
                      linewidth=0.5, label=legend)
        for bar, value, bottom in zip(bars, values, bottoms):
            ax.text(bar.get_x() + bar.get_width() / 2, bottom + value / 2, f"{value}",
                    ha="center", va="center", color="white", fontsize=7.5,
                    fontweight="bold")
        bottoms += values
    for xpos, total in zip(x, totals):
        ax.text(xpos, total + 18, f"{total}", ha="center", va="bottom",
                fontsize=7.5, fontweight="bold")
    ax.set_ylim(0, 850)
    ax.set_ylabel("Number of benchmark records")
    ax.set_xticks(x, [f"{name}\n(total {total})" for name, total in zip(labels, totals)])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=6.8, labelspacing=0.8)


def draw_metadata_panel(ax: plt.Axes) -> None:
    data = read_table("S5_Condition_metadata")
    x = np.arange(len(data))
    width = 0.31
    ph = data["pH available (%)"].to_numpy(float)
    temperature = data["Temperature available (%)"].to_numpy(float)
    bars_ph = ax.bar(x - width / 2, ph, width, color="#25599A", label="pH available")
    bars_temperature = ax.bar(x + width / 2, temperature, width, color="#158676",
                              label="Temperature available")
    for bars, values, counts, color in [
        (bars_ph, ph, data["pH available (n)"], "#25599A"),
        (bars_temperature, temperature, data["Temperature available (n)"], "#158676"),
    ]:
        for bar, value, count, total in zip(bars, values, counts, data["Records"]):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 1.6,
                    f"{value:.2f}%\n({int(count)}/{int(total)})", ha="center", va="bottom",
                    fontsize=7.0, color=color, fontweight="bold", linespacing=1.05)
    ax.set_ylim(0, 104)
    ax.set_ylabel("Availability (%)")
    ax.set_xticks(x, ["E. coli\n(n = 781)", "S. cerevisiae\n(n = 465)"])
    ax.grid(axis="y", linestyle="--", linewidth=0.55, alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=2)


def figure1() -> None:
    panel_a, panel_b = exact_figure1_reference_panels()
    fig = plt.figure(figsize=(9.5, 9.6), facecolor="white")
    grid = fig.add_gridspec(
        3, 2,
        height_ratios=[1.12, 0.90, 0.90],
        width_ratios=[1.0, 1.0],
        left=0.055, right=0.975, bottom=0.065, top=0.985,
        hspace=0.32, wspace=0.34,
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[2, 0])
    ax_e = fig.add_subplot(grid[2, 1])
    ax_a.imshow(panel_a)
    ax_b.imshow(panel_b)
    ax_a.axis("off")
    ax_b.axis("off")
    draw_matching_panel(ax_c)
    draw_source_panel(ax_d)
    draw_metadata_panel(ax_e)
    labels = [
        panel_label(fig, 0.012, 0.985, "a"),
        panel_label(fig, 0.012, 0.625, "b"),
        panel_label(fig, 0.505, 0.625, "c"),
        panel_label(fig, 0.012, 0.335, "d"),
        panel_label(fig, 0.505, 0.335, "e"),
    ]

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            ["all-reactant candidates", 16613],
            ["experimentally matched records", 1354],
            ["final benchmark records", 1246],
            ["E. coli records", 781],
            ["S. cerevisiae records", 465],
            ["distinct reactions", 773],
            ["unique protein sequences", 514],
            ["distinct EC annotations", 390],
            ["median log10(kcat)", 1.279],
            ["inference regimes", 6],
        ],
        columns=["metric", "value"],
    ).to_csv(SOURCE_DIR / "Figure1ab_updated_values.csv", index=False)
    read_table("S3_Matching_levels").to_csv(SOURCE_DIR / "Figure1c_matching_support.csv", index=False)
    read_table("S4_Source_support").to_csv(SOURCE_DIR / "Figure1d_source_support.csv", index=False)
    read_table("S5_Condition_metadata").to_csv(SOURCE_DIR / "Figure1e_condition_metadata.csv", index=False)

    save_axes_panel(fig, [ax_c], "Figure1c", hidden_artists=[labels[2]])
    save_axes_panel(fig, [ax_d], "Figure1d", hidden_artists=[labels[3]])
    save_axes_panel(fig, [ax_e], "Figure1e", hidden_artists=[labels[4]])
    save_composite(fig, 1)


def available_scope_data() -> dict[str, pd.DataFrame]:
    reaction = read_table("S6_Reaction_subset").rename(
        columns={"mae_log10_common_subset": "mae_log10", "n_common_subset": "n"}
    )
    available = read_table("S7_Available_case")
    catpred_scope = next(value for value in available["subset"].unique()
                         if str(value).startswith("catpred_accessible_scope_"))
    kinform_scope = next(value for value in available["subset"].unique()
                         if str(value).startswith("kinform_accessible_scope_"))
    return {
        "Reaction-complete subset\n(scope n=1,047; strictly common)": reaction,
        "CatPred-accessible scope\n(scope n=1,156; available case)":
            available.loc[available["subset"].eq(catpred_scope)],
        "KinForm-L-accessible scope\n(scope n=729; available case)":
            available.loc[available["subset"].eq(kinform_scope)],
    }


def cluster_bootstrap_differences(replicates: int = 2000, seed: int = 20260812) -> pd.DataFrame:
    table0 = pd.read_csv(TABLE0, usecols=["entry_id", "method", "absolute_error_log10"])
    wide = table0.pivot(index="entry_id", columns="method", values="absolute_error_log10")
    common = wide[["KcatNet", "TurNuP", "PMAK"]].dropna().reset_index()
    audit = pd.read_csv(
        RECORD_AUDIT,
        usecols=["entry_id", "protein_cluster", "pair_cluster", "reaction_cluster",
                 "reference_cluster", "label_assignment_cluster"],
    )
    common = common.merge(audit, on="entry_id", how="left", validate="one_to_one")
    if len(common) != 1047:
        raise ValueError(f"Expected 1,047 reaction-common rows, observed {len(common)}")
    rows = []
    for pair_index, (method_a, method_b) in enumerate(PAIR_ORDER):
        common["_difference"] = common[method_a] - common[method_b]
        point = float(common["_difference"].mean())
        for cluster_index, (cluster_type, cluster_label) in enumerate(CLUSTER_TYPES):
            cluster_column = f"{cluster_type}_cluster"
            grouped = common.groupby(cluster_column, sort=True)["_difference"].agg(["sum", "size"])
            sums = grouped["sum"].to_numpy(float)
            sizes = grouped["size"].to_numpy(float)
            rng = np.random.default_rng(seed + pair_index * 100 + cluster_index)
            indices = rng.integers(0, len(grouped), size=(replicates, len(grouped)))
            boot = sums[indices].sum(axis=1) / sizes[indices].sum(axis=1)
            low, high = np.percentile(boot, [2.5, 97.5])
            rows.append(
                {
                    "method_a": method_a,
                    "method_b": method_b,
                    "cluster_type": cluster_type,
                    "cluster_label": cluster_label,
                    "n_common_rows": len(common),
                    "n_clusters": len(grouped),
                    "mae_difference_a_minus_b": point,
                    "difference_ci_low_95": float(low),
                    "difference_ci_high_95": float(high),
                    "bootstrap_replicates": replicates,
                    "seed": seed + pair_index * 100 + cluster_index,
                    "estimand": "row-weighted MAE difference with intact clusters resampled",
                }
            )
    return pd.DataFrame(rows)


def weighting_sensitivity() -> pd.DataFrame:
    reaction = read_table("S6_Reaction_subset").set_index("method")
    sensitivity = read_table("S17_Sensitivity_subsets")
    rows = []
    scopes = [
        ("Row-weighted", None),
        ("Unique seq–substrate pair", "reaction_common_unique_sequence_substrate_pairs"),
        ("Label assignment", "reaction_common_unique_label_assignments"),
    ]
    for label, scope in scopes:
        for method in ["KcatNet", "TurNuP", "PMAK"]:
            if scope is None:
                value = float(reaction.loc[method, "mae_log10_common_subset"])
                units = int(reaction.loc[method, "n_common_subset"])
            else:
                row = sensitivity.loc[sensitivity["analysis_scope"].eq(scope)
                                      & sensitivity["method"].eq(method)].iloc[0]
                value = float(row["mae_log10"])
                units = int(row["scope_units"])
            rows.append({"analysis": label, "method": method, "mae_log10": value,
                         "scope_units": units})
    return pd.DataFrame(rows)


def draw_scope_bars(ax: plt.Axes, frame: pd.DataFrame, methods: list[str], title: str) -> None:
    frame = frame.set_index("method").loc[methods].reset_index()
    frame = frame.iloc[::-1]
    y = np.arange(len(frame))
    bars = ax.barh(y, frame["mae_log10"], height=0.58,
                   color=[OLD_METHOD_COLORS[method] for method in frame["method"]])
    ax.set_yticks(y, [method for method in frame["method"]])
    ax.set_xlim(0, 1.12)
    ax.set_xticks(np.arange(0, 1.2, 0.2))
    ax.set_xlabel(r"MAE in log$_{10}$($k_{cat}$)")
    ax.set_title(title, fontsize=9.3, fontweight="bold", pad=8)
    ax.grid(axis="x", color="#D9D9D9", linestyle="--", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)
    for bar, row in zip(bars, frame.itertuples()):
        ax.text(row.mae_log10 + 0.014, bar.get_y() + bar.get_height() / 2,
                f"{row.mae_log10:.3f} (n={int(row.n):,})", ha="left", va="center",
                fontsize=6.8)


def draw_weighting_panel(ax: plt.Axes, data: pd.DataFrame) -> None:
    y_positions = {
        ("Row-weighted", "KcatNet"): 8,
        ("Row-weighted", "TurNuP"): 7,
        ("Row-weighted", "PMAK"): 6,
        ("Unique seq–substrate pair", "KcatNet"): 4,
        ("Unique seq–substrate pair", "TurNuP"): 3,
        ("Unique seq–substrate pair", "PMAK"): 2,
        ("Label assignment", "KcatNet"): 0,
        ("Label assignment", "TurNuP"): -1,
        ("Label assignment", "PMAK"): -2,
    }
    xmin = min(0.66, float(data["mae_log10"].min()) - 0.015)
    xmax = max(0.79, float(data["mae_log10"].max()) + 0.025)
    for row in data.itertuples():
        y = y_positions[(row.analysis, row.method)]
        color = OLD_METHOD_COLORS[row.method]
        ax.hlines(y, xmin, row.mae_log10, color=color, linewidth=1.4)
        ax.scatter(row.mae_log10, y, color=color, s=34, zorder=3)
        ax.text(row.mae_log10 + 0.004, y, f"{row.mae_log10:.3f}", va="center", fontsize=7.2)
        ax.text(-0.035, y, row.method, transform=ax.get_yaxis_transform(), ha="right",
                va="center", color=color, fontsize=7.5)
    ax.axhline(5, color="#9E9E9E", linestyle=":", linewidth=0.9)
    ax.axhline(1, color="#9E9E9E", linestyle=":", linewidth=0.9)
    for label, ypos in [("Row-weighted", 7), ("Unique\nseq–substrate pair", 3),
                        ("Label\nassignment", -1)]:
        ax.text(-0.35, ypos, label, transform=ax.get_yaxis_transform(), ha="center",
                va="center", fontsize=8.0, fontweight="bold")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-2.8, 8.8)
    ax.set_yticks([])
    ax.set_xlabel(r"MAE in log$_{10}$($k_{cat}$)")
    ax.set_title("Sensitivity to the unit of record weighting", fontsize=9.5,
                 fontweight="bold", pad=8)
    ax.grid(axis="x", color="#D9D9D9", linestyle="--", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)


def figure3() -> None:
    scopes = available_scope_data()
    cluster = cluster_bootstrap_differences()
    weighting = weighting_sensitivity()
    fig = plt.figure(figsize=(14.2, 9.8), facecolor="white")
    outer = fig.add_gridspec(
        2, 2,
        height_ratios=[1.28, 0.92], width_ratios=[1.85, 1.0],
        left=0.065, right=0.985, bottom=0.075, top=0.955,
        hspace=0.34, wspace=0.38,
    )
    top = outer[0, :].subgridspec(1, 3, wspace=0.42)
    axes_a = [fig.add_subplot(top[0, index]) for index in range(3)]
    methods = [
        ["KcatNet", "TurNuP", "PMAK", "CataPro"],
        ["KcatNet", "TurNuP", "PMAK", "CataPro", "CatPred"],
        ["KcatNet", "TurNuP", "PMAK", "CataPro", "KinForm-L"],
    ]
    source_a = []
    for axis, (title, frame), method_order in zip(axes_a, scopes.items(), methods):
        draw_scope_bars(axis, frame, method_order, title)
        selected = frame.set_index("method").loc[method_order].reset_index().copy()
        selected["scope"] = title.replace("\n", " ")
        source_a.append(selected)

    bottom_b = outer[1, 0].subgridspec(1, 3, wspace=0.32)
    axes_b = [fig.add_subplot(bottom_b[0, index]) for index in range(3)]
    y = np.arange(len(CLUSTER_TYPES))[::-1]
    for index, (axis, pair) in enumerate(zip(axes_b, PAIR_ORDER)):
        frame = cluster.loc[cluster["method_a"].eq(pair[0]) & cluster["method_b"].eq(pair[1])]
        frame = frame.set_index("cluster_label").loc[[label for _key, label in CLUSTER_TYPES]].reset_index()
        values = frame["mae_difference_a_minus_b"].to_numpy(float)
        low = frame["difference_ci_low_95"].to_numpy(float)
        high = frame["difference_ci_high_95"].to_numpy(float)
        axis.errorbar(values, y, xerr=np.vstack([values - low, high - values]), fmt="o",
                      color="#0072B2", ecolor="#0072B2", capsize=2.3, lw=1.0,
                      markersize=4.3)
        axis.axvline(0, color="#777777", linestyle="--", linewidth=0.85)
        limit = max(abs(cluster["difference_ci_low_95"].min()),
                    abs(cluster["difference_ci_high_95"].max())) * 1.12
        axis.set_xlim(-limit, limit)
        axis.set_yticks(y)
        if index == 0:
            axis.set_yticklabels([label for _key, label in CLUSTER_TYPES])
        else:
            axis.set_yticklabels([])
        axis.set_xlabel("MAE difference\n(first minus second)")
        axis.set_title(f"{pair[0]} – {pair[1]}", fontsize=8.5, fontweight="bold")
        axis.grid(axis="x", color="#E0E0E0", linewidth=0.5, alpha=0.8)
        axis.set_axisbelow(True)
    ax_c = fig.add_subplot(outer[1, 1])
    draw_weighting_panel(ax_c, weighting)

    label_a = panel_label(fig, 0.016, 0.977, "a")
    label_b = panel_label(fig, 0.016, 0.405, "b")
    label_c = panel_label(fig, 0.705, 0.405, "c")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat(source_a, ignore_index=True).to_csv(SOURCE_DIR / "Figure3a_scope_bars.csv", index=False)
    cluster.to_csv(SOURCE_DIR / "Figure3b_cluster_bootstrap_differences.csv", index=False)
    weighting.to_csv(SOURCE_DIR / "Figure3c_weighting_sensitivity.csv", index=False)
    save_axes_panel(fig, axes_a, "Figure3a", hidden_artists=[label_a])
    save_axes_panel(fig, axes_b, "Figure3b", hidden_artists=[label_b])
    hidden_artists = [*axes_a, *axes_b, label_a, label_b]
    for artist in hidden_artists:
        artist.set_visible(False)
    save_axes_panel(fig, [ax_c], "Figure3c", hidden_artists=[label_c])
    for artist in hidden_artists:
        artist.set_visible(True)
    save_composite(fig, 3)


def subgroup_value(frame: pd.DataFrame, method: str, feature: str, group: str,
                   column: str) -> float:
    row = frame.loc[frame["method"].eq(method) & frame["feature"].eq(feature)
                    & frame["group"].astype(str).eq(group)]
    if len(row) != 1:
        raise ValueError(f"Expected one row for {method}/{feature}/{group}, observed {len(row)}")
    return float(row.iloc[0][column])


def cluster_group_bootstrap(
    group: pd.DataFrame,
    value_column: str,
    replicates: int,
    seed: int,
    statistic: str,
) -> tuple[float, float, float, int]:
    grouped = group.groupby("protein_cluster")[value_column].agg(["sum", "size"])
    sums = grouped["sum"].to_numpy(float)
    sizes = grouped["size"].to_numpy(float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(grouped), size=(replicates, len(grouped)))
    boots = sums[indices].sum(axis=1) / sizes[indices].sum(axis=1)
    point = float(group[value_column].mean())
    low, high = np.percentile(boots, [2.5, 97.5])
    return point, float(low), float(high), len(grouped)


def independent_cluster_difference_bootstrap(
    first: pd.DataFrame,
    second: pd.DataFrame,
    value_column: str,
    replicates: int,
    seed: int,
) -> tuple[np.ndarray, dict[str, int | str]]:
    """Bootstrap two disjoint subgroup means using clusters within each group."""

    def bootstrap_mean(group: pd.DataFrame, local_seed: int) -> tuple[np.ndarray, int]:
        grouped = group.groupby("protein_cluster")[value_column].agg(["sum", "size"])
        sums = grouped["sum"].to_numpy(float)
        sizes = grouped["size"].to_numpy(float)
        rng = np.random.default_rng(local_seed)
        indices = rng.integers(0, len(grouped), size=(replicates, len(grouped)))
        return sums[indices].sum(axis=1) / sizes[indices].sum(axis=1), len(grouped)

    first_boot, first_clusters = bootstrap_mean(first, seed)
    second_boot, second_clusters = bootstrap_mean(second, seed + 1)
    shared = len(set(first["protein_cluster"]) & set(second["protein_cluster"]))
    if shared:
        raise ValueError("Independent subgroup bootstrap requires disjoint protein clusters")
    return first_boot - second_boot, {
        "n_clusters_first": first_clusters,
        "n_clusters_second": second_clusters,
        "n_shared_clusters": 0,
        "n_union_clusters": first_clusters + second_clusters,
        "bootstrap_scheme": "independent protein-cluster resampling within disjoint groups",
        "seed_first": seed,
        "seed_second": seed + 1,
    }


def joint_cluster_difference_bootstrap(
    first: pd.DataFrame,
    second: pd.DataFrame,
    value_column: str,
    replicates: int,
    seed: int,
) -> tuple[np.ndarray, dict[str, int | str]]:
    """Resample the cluster union once so shared proteins retain joint dependence."""
    first_grouped = first.groupby("protein_cluster")[value_column].agg(["sum", "size"])
    second_grouped = second.groupby("protein_cluster")[value_column].agg(["sum", "size"])
    clusters = first_grouped.index.union(second_grouped.index)
    first_grouped = first_grouped.reindex(clusters, fill_value=0)
    second_grouped = second_grouped.reindex(clusters, fill_value=0)

    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(clusters), size=(replicates, len(clusters)))
    first_denominator = first_grouped["size"].to_numpy(float)[indices].sum(axis=1)
    second_denominator = second_grouped["size"].to_numpy(float)[indices].sum(axis=1)
    if np.any(first_denominator == 0) or np.any(second_denominator == 0):
        raise ValueError("A joint cluster-bootstrap replicate omitted an entire subgroup")
    first_boot = first_grouped["sum"].to_numpy(float)[indices].sum(axis=1) / first_denominator
    second_boot = second_grouped["sum"].to_numpy(float)[indices].sum(axis=1) / second_denominator
    shared = len(first_grouped.index.intersection(
        first.groupby("protein_cluster").size().index
    ).intersection(second.groupby("protein_cluster").size().index))
    return first_boot - second_boot, {
        "n_clusters_first": int(first["protein_cluster"].nunique()),
        "n_clusters_second": int(second["protein_cluster"].nunique()),
        "n_shared_clusters": shared,
        "n_union_clusters": len(clusters),
        "bootstrap_scheme": "joint protein-cluster resampling across both groups",
        "seed_first": seed,
        "seed_second": seed,
    }


def figure4_bootstrap(replicates: int = 2000, seed: int = 20260812) -> tuple[pd.DataFrame, pd.DataFrame]:
    table0 = pd.read_csv(
        TABLE0,
        usecols=["entry_id", "method", "species", "experimental_substrate_support",
                 "substrate_role_group", "absolute_error_log10"],
    ).dropna(subset=["absolute_error_log10"])
    audit = pd.read_csv(RECORD_AUDIT, usecols=["entry_id", "protein_cluster"])
    table0 = table0.merge(audit, on="entry_id", how="left", validate="many_to_one")
    table0["large_error"] = (table0["absolute_error_log10"] > 1).astype(float)
    large_rows = []
    contrast_rows = []
    for method_index, method in enumerate(METHODS_4):
        method_rows = table0.loc[table0["method"].eq(method)].copy()
        species_boots: dict[str, tuple[float, np.ndarray]] = {}
        for species_index, species in enumerate(["ecoli", "yeast"]):
            group = method_rows.loc[method_rows["species"].eq(species)]
            point, low, high, clusters = cluster_group_bootstrap(
                group, "large_error", replicates,
                seed + method_index * 1000 + species_index,
                "fraction",
            )
            large_rows.append(
                {
                    "method": method,
                    "species": species,
                    "n_rows": len(group),
                    "n_clusters": clusters,
                    "large_error_fraction": point,
                    "ci_low_95": low,
                    "ci_high_95": high,
                    "bootstrap_unit": "protein_cluster",
                    "bootstrap_replicates": replicates,
                    "seed": seed + method_index * 1000 + species_index,
                }
            )

        contrast_specs = [
            (
                "E. coli – S. cerevisiae",
                method_rows.loc[method_rows["species"].eq("ecoli")],
                method_rows.loc[method_rows["species"].eq("yeast")],
                "independent",
            ),
            (
                "Supported other reactant – supported currency/cofactor",
                method_rows.loc[
                    method_rows["experimental_substrate_support"].eq("substrate_supported")
                    & method_rows["substrate_role_group"].eq("other_reactant")
                ],
                method_rows.loc[
                    method_rows["experimental_substrate_support"].eq("substrate_supported")
                    & method_rows["substrate_role_group"].eq("currency_or_cofactor")
                ],
                "joint",
            ),
        ]
        for contrast_index, (contrast, first, second, scheme) in enumerate(contrast_specs):
            local_seed = seed + 10000 + method_index * 1000 + contrast_index * 10
            if scheme == "independent":
                difference, bootstrap_metadata = independent_cluster_difference_bootstrap(
                    first, second, "absolute_error_log10", replicates, local_seed
                )
            else:
                difference, bootstrap_metadata = joint_cluster_difference_bootstrap(
                    first, second, "absolute_error_log10", replicates, local_seed
                )
            low, high = np.percentile(difference, [2.5, 97.5])
            row = {
                "method": method,
                "contrast": contrast,
                "n_first": len(first),
                "n_second": len(second),
                "mae_difference": float(first["absolute_error_log10"].mean()
                                        - second["absolute_error_log10"].mean()),
                "ci_low_95": float(low),
                "ci_high_95": float(high),
                "bootstrap_unit": "protein_cluster",
                "bootstrap_replicates": replicates,
            }
            row.update(bootstrap_metadata)
            contrast_rows.append(row)
    return pd.DataFrame(large_rows), pd.DataFrame(contrast_rows)


def draw_figure4_heatmap(ax: plt.Axes, stratified: pd.DataFrame) -> plt.Axes:
    categories = [
        ("species", "ecoli", "E. coli"),
        ("species", "yeast", "S. cerevisiae"),
        ("source_database", "BRENDA", "BRENDA only"),
        ("source_database", "SABIO-RK", "SABIO-RK only"),
        ("source_database", "BRENDA;SABIO-RK", "Both"),
        ("substrate_role_group_substrate_supported", "other_reactant", "Supported\nother reactant"),
        ("substrate_role_group_substrate_supported", "currency_or_cofactor", "Supported\ncurrency/cofactor"),
        ("substrate_role_group_substrate_supported", "carrier_linked_variable", "Supported\ncarrier-linked"),
    ]
    matrix = np.asarray([
        [subgroup_value(stratified, method, feature, group, "mean_abs_error_log10")
         for feature, group, _label in categories]
        for method in METHODS_4
    ])
    counts = np.asarray([
        [int(subgroup_value(stratified, method, feature, group, "n"))
         for feature, group, _label in categories]
        for method in METHODS_4
    ])
    finite = matrix[np.isfinite(matrix)]
    vmin = max(0.45, float(finite.min()) - 0.03)
    vmax = min(1.15, float(finite.max()) + 0.03)
    image = ax.imshow(matrix, cmap="YlOrRd", aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(METHODS_4)), METHODS_4)
    ax.set_xticks(np.arange(len(categories)), [label for _f, _g, label in categories])
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=5)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            threshold = vmin + 0.66 * (vmax - vmin)
            color = "white" if matrix[i, j] > threshold else "black"
            ax.text(j, i, f"{matrix[i, j]:.3f}\n(n={counts[i, j]})",
                    ha="center", va="center", fontsize=7.1, color=color,
                    linespacing=1.1)
    for boundary in [1.5, 4.5]:
        ax.axvline(boundary, color="white", linewidth=2.2)
    group_specs = [(-0.45, 1.45, 0.5, "Species"),
                   (1.55, 4.45, 3.0, "Experimental source"),
                   (4.55, 7.45, 6.0, "Substrate role")]
    for start, end, center, label in group_specs:
        ax.plot([start, end], [1.105, 1.105], transform=ax.get_xaxis_transform(),
                color="#4F4F4F", linewidth=0.8, clip_on=False)
        ax.text(center, 1.145, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    colorbar = ax.figure.colorbar(image, ax=ax, fraction=0.018, pad=0.018)
    colorbar.set_label(r"MAE in log$_{10}$($k_{cat}$)")
    return colorbar.ax


def figure4() -> None:
    stratified = read_table("S9_Error_stratification")
    large, contrasts = figure4_bootstrap()
    fig = plt.figure(figsize=(14.2, 9.8), facecolor="white")
    outer = fig.add_gridspec(
        2, 2,
        height_ratios=[1.18, 0.92], width_ratios=[0.92, 1.48],
        left=0.075, right=0.985, bottom=0.075, top=0.89,
        hspace=0.34, wspace=0.32,
    )
    ax_a = fig.add_subplot(outer[0, :])
    ax_b = fig.add_subplot(outer[1, 0])
    bottom_c = outer[1, 1].subgridspec(1, 2, wspace=0.30)
    axes_c = [fig.add_subplot(bottom_c[0, index]) for index in range(2)]
    colorbar_ax = draw_figure4_heatmap(ax_a, stratified)

    y = np.arange(len(METHODS_4))[::-1]
    offsets = {"ecoli": 0.10, "yeast": -0.10}
    colors = {"ecoli": "#0072B2", "yeast": "#D62728"}
    labels = {"ecoli": "E. coli", "yeast": "S. cerevisiae"}
    for ypos, method in zip(y, METHODS_4):
        method_rows = large.loc[large["method"].eq(method)].set_index("species")
        values = [float(method_rows.loc[species, "large_error_fraction"])
                  for species in ["ecoli", "yeast"]]
        ax_b.hlines(ypos, min(values), max(values), color="#AEB8C2", linewidth=1.2,
                    zorder=1)
        for species in ["ecoli", "yeast"]:
            row = method_rows.loc[species]
            value = float(row["large_error_fraction"])
            ax_b.errorbar(
                value,
                ypos + offsets[species],
                xerr=np.asarray([[value - float(row["ci_low_95"])],
                                 [float(row["ci_high_95"]) - value]]),
                fmt="o", color=colors[species], ecolor=colors[species], capsize=2.2,
                lw=1.0, markersize=4.2, label=labels[species]
                if method == METHODS_4[0] else None,
            )
    ax_b.set_yticks(y, METHODS_4)
    ax_b.set_xlabel(r"Fraction with absolute log$_{10}$ error > 1")
    ax_b.set_title("Large-error fraction by species", fontsize=9.5,
                   fontweight="bold", pad=8)
    ax_b.grid(axis="x", color="#D9D9D9", linestyle="--", linewidth=0.55, alpha=0.8)
    ax_b.set_axisbelow(True)
    ax_b.legend(frameon=True, loc="upper right", bbox_to_anchor=(0.985, 0.985), ncol=1)

    contrast_order = [
        "E. coli – S. cerevisiae",
        "Supported other reactant – supported currency/cofactor",
    ]
    for axis_index, (axis, contrast) in enumerate(zip(axes_c, contrast_order)):
        frame = contrasts.loc[contrasts["contrast"].eq(contrast)].set_index("method").loc[METHODS_4]
        values = frame["mae_difference"].to_numpy(float)
        low = frame["ci_low_95"].to_numpy(float)
        high = frame["ci_high_95"].to_numpy(float)
        axis.errorbar(values, y, xerr=np.vstack([values - low, high - values]), fmt="o",
                      color="#0072B2", ecolor="#0072B2", capsize=2.2, lw=1.0,
                      markersize=4.2)
        axis.axvline(0, color="#777777", linestyle="--", linewidth=0.85)
        axis.set_yticks(y)
        if axis_index == 0:
            axis.set_yticklabels(METHODS_4)
        else:
            axis.set_yticklabels([])
        axis.set_xlabel("MAE difference")
        axis.set_title(contrast, fontsize=8.5, pad=8)
        axis.grid(axis="x", color="#E0E0E0", linewidth=0.5, alpha=0.8)
        axis.set_axisbelow(True)
    all_bounds = np.r_[contrasts["ci_low_95"].to_numpy(float),
                       contrasts["ci_high_95"].to_numpy(float)]
    limit = max(abs(all_bounds.min()), abs(all_bounds.max())) * 1.08
    for axis in axes_c:
        axis.set_xlim(-limit, limit)

    label_a = panel_label(fig, 0.016, 0.965, "a")
    label_b = panel_label(fig, 0.016, 0.405, "b")
    label_c = panel_label(fig, 0.505, 0.405, "c")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    stratified.to_csv(SOURCE_DIR / "Figure4a_stratified_mae.csv", index=False)
    large.to_csv(SOURCE_DIR / "Figure4b_large_error_protein_cluster_bootstrap.csv", index=False)
    contrasts.to_csv(SOURCE_DIR / "Figure4c_subgroup_mae_difference_bootstrap.csv", index=False)
    save_axes_panel(fig, [ax_a, colorbar_ax], "Figure4a", hidden_artists=[label_a])
    save_axes_panel(fig, [ax_b], "Figure4b", hidden_artists=[label_b])
    save_axes_panel(fig, axes_c, "Figure4c", hidden_artists=[label_c])
    save_composite(fig, 4)


def figure2() -> None:
    """Retain the accepted Figure 2 and make each exported panel self-contained."""
    r4.OUTPUT_DIR = OUTPUT_DIR
    r4.PANEL_DIR = PANEL_DIR
    r4.SOURCE_DIR = SOURCE_DIR
    r4.PACKAGE_DIR = PACKAGE_DIR
    r4.configure_style()
    original_save_axis_panel = r4.save_axis_panel

    def save_axis_panel_with_legend(fig: plt.Figure, ax: plt.Axes, stem: str) -> None:
        panel_identifier = [
            artist for artist in ax.texts
            if artist.get_text().strip() == stem[-1]
            and artist.get_fontweight() == "bold"
        ]
        regime = None
        extra_artists = [*ax.get_yticklabels()]
        if stem in {"Figure2b", "Figure2c", "Figure2d", "Figure2e", "Figure2f"}:
            existing_legend = ax.get_legend()
            if existing_legend is not None:
                ax.add_artist(existing_legend)
            anchor_y = -0.25 if stem == "Figure2f" else -0.20
            regime = r4.regime_legend(
                ax,
                loc="upper center",
                bbox_to_anchor=(0.5, anchor_y),
                ncol=3,
            )
            extra_artists.insert(0, regime)
        other_axes = [other for other in fig.axes if other is not ax]
        for other in other_axes:
            other.set_visible(False)
        try:
            save_axes_panel(
                fig,
                [ax],
                stem,
                extra_artists=extra_artists,
                hidden_artists=panel_identifier,
            )
        finally:
            for other in other_axes:
                other.set_visible(True)
            if regime is not None:
                regime.remove()

    r4.save_axis_panel = save_axis_panel_with_legend
    try:
        r4.figure2()
    finally:
        r4.save_axis_panel = original_save_axis_panel


def make_packages() -> None:
    PACKAGE_DIR.mkdir(parents=True, exist_ok=True)
    source_names = {
        1: ["Figure1ab_updated_values.csv", "Figure1c_matching_support.csv",
            "Figure1d_source_support.csv", "Figure1e_condition_metadata.csv"],
        2: ["Figure2_plot_data.csv", "Figure2_pair_cluster_intervals.csv",
            "Figure2_row_bootstrap_intervals.csv"],
        3: ["Figure3a_scope_bars.csv", "Figure3b_cluster_bootstrap_differences.csv",
            "Figure3c_weighting_sensitivity.csv"],
        4: ["Figure4a_stratified_mae.csv",
            "Figure4b_large_error_protein_cluster_bootstrap.csv",
            "Figure4c_subgroup_mae_difference_bootstrap.csv"],
    }
    panel_letters = {1: "abcde", 2: "abcdef", 3: "abc", 4: "abc"}
    script = Path(__file__)
    for number in range(1, 5):
        archive_path = PACKAGE_DIR / f"Figure{number}_data_code_and_panels.zip"
        with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
            archive.write(script, arcname=script.name)
            archive.write(HELPER_SCRIPT, arcname=HELPER_SCRIPT.name)
            note = (
                "Generated from benchmark v1.2.0 data using manuscript figure revision 1.2.0-r5.\n"
                "Run from the repository root with:\n"
                "  python paper/generate_manuscript_figures.py\n"
            )
            if number == 1:
                note += (
                    "Figure 1a/b use versioned high-resolution base panels; only benchmark values are replaced.\n"
                )
                archive.write(REFERENCE_A, arcname=f"references/{REFERENCE_A.name}")
                archive.write(REFERENCE_B, arcname=f"references/{REFERENCE_B.name}")
            archive.writestr("README.txt", note)
            for name in source_names[number]:
                archive.write(SOURCE_DIR / name, arcname=f"data/{name}")
            archive.write(OUTPUT_DIR / f"Figure{number}.png", arcname=f"Figure{number}.png")
            archive.write(OUTPUT_DIR / f"Figure{number}.pdf", arcname=f"Figure{number}.pdf")
            for letter in panel_letters[number]:
                for suffix in ["png", "pdf"]:
                    path = PANEL_DIR / f"Figure{number}{letter}.{suffix}"
                    archive.write(path, arcname=f"panels/{path.name}")


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
    print(f"Wrote manuscript figure revision 1.2.0-r5 to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
