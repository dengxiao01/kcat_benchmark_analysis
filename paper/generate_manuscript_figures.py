#!/usr/bin/env python3
"""Validate the 0814/V4 figure snapshot and rerun its plotting packages.

The composite PNG files are frozen author-approved raster exports. The
standalone panel set follows the final V4 order; Figure 2c-f are regenerated
from the supplied plotting package to normalize its pre-assembly ordering. The
four data-backed packages are retained so their numerical panels can be rerun.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image


PAPER = Path(__file__).resolve().parent
FIGURES = PAPER / "figures"
SOURCES = PAPER / "figure_sources_0814"
MANIFEST = SOURCES / "figure_asset_manifest.csv"

PACKAGE_COMMANDS = (
    ("Figure1", "plot_Figure1_panels_d_e_f.py"),
    ("Figure2", "Figure2_plot.py"),
    ("Figure3", "Figure3_plotting_code.py"),
    ("Figure4", "plot_figure4_panels_v2.py"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"Empty manifest: {MANIFEST}")
    return rows


def validate_snapshot() -> None:
    failures: list[str] = []
    rows = read_manifest()

    for row in rows:
        path = PAPER / row["path"]
        if not path.is_file():
            failures.append(f"missing: {path.relative_to(PAPER)}")
            continue

        observed_hash = sha256(path)
        if observed_hash != row["sha256"]:
            failures.append(
                f"hash mismatch: {path.relative_to(PAPER)} "
                f"({observed_hash} != {row['sha256']})"
            )

        if path.suffix.lower() == ".png":
            with Image.open(path) as image:
                observed_size = f"{image.width}x{image.height}"
                image.verify()
            if observed_size != row["dimensions"]:
                failures.append(
                    f"dimension mismatch: {path.relative_to(PAPER)} "
                    f"({observed_size} != {row['dimensions']})"
                )

    composites = sorted(FIGURES.glob("Figure[1-4].png"))
    panels = sorted((FIGURES / "panels").glob("Figure*.png"))
    if len(composites) != 4:
        failures.append(f"expected 4 composite PNGs, found {len(composites)}")
    if len(panels) != 20:
        failures.append(f"expected 20 panel PNGs, found {len(panels)}")

    if failures:
        raise RuntimeError("0814/V4 snapshot validation failed:\n- " + "\n- ".join(failures))

    print(
        "0814/V4 snapshot OK: "
        f"{len(composites)} composites, {len(panels)} panels, "
        f"{len(rows)} checks"
    )


def rebuild_code_backed_panels(destination: Path) -> None:
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)

    for package_name, script_name in PACKAGE_COMMANDS:
        source = SOURCES / package_name
        target = destination / package_name
        if target.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing rebuild directory: {target}"
            )
        shutil.copytree(source, target)
        subprocess.run([sys.executable, script_name], cwd=target, check=True)

    generated_pngs = sorted(destination.rglob("*.png"))
    generated_pdfs = sorted(destination.rglob("*.pdf"))
    print(
        f"Rebuilt {len(generated_pngs)} PNG and {len(generated_pdfs)} PDF "
        f"code-backed panels under {destination}"
    )
    print(
        "Figure 1a/1b are frozen artwork without a plotting script. "
        "The Figure 4 package also regenerates its provenance panel, which "
        "was not retained in the final four-panel V4 composition."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rebuild-code-panels",
        type=Path,
        metavar="OUTPUT_DIR",
        help="rerun the four 0814 plotting packages in a new output directory",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validate_snapshot()
    if args.rebuild_code_panels:
        rebuild_code_backed_panels(args.rebuild_code_panels)


if __name__ == "__main__":
    main()
