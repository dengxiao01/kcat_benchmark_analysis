#!/usr/bin/env python3
"""Collect DEKP structure files into stable benchmark/training directories."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = BASE / "data" / "final" / "dekp" / "dekp_kcat_input_metadata.csv"
DEFAULT_DEKP_DATA = BASE / "external_methods" / "DEKP" / "datasets" / "kcat_dataset.csv"
DEFAULT_AUTHOR_ROOT = BASE / "external_methods" / "DEKP" / "structures" / "author_archive"
DEFAULT_BENCHMARK_DIR = BASE / "external_methods" / "DEKP" / "structures" / "benchmark" / "AlphaFold"
DEFAULT_TRAIN_DIR = BASE / "external_methods" / "DEKP" / "structures" / "public_kcat" / "AlphaFold"
DEFAULT_REPORT = BASE / "reports" / "tables" / "dekp_structure_collection_report.csv"
DEFAULT_MISSING = BASE / "data" / "final" / "dekp" / "dekp_missing_structure_after_collection.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Link DEKP PDB files from the author archive.")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--dekp-data", type=Path, default=DEFAULT_DEKP_DATA)
    parser.add_argument("--author-root", type=Path, default=DEFAULT_AUTHOR_ROOT)
    parser.add_argument("--benchmark-dir", type=Path, default=DEFAULT_BENCHMARK_DIR)
    parser.add_argument("--train-dir", type=Path, default=DEFAULT_TRAIN_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--missing", type=Path, default=DEFAULT_MISSING)
    return parser.parse_args()


def clean_ids(series: pd.Series) -> list[str]:
    values = series.dropna().astype(str).str.strip()
    return sorted({value for value in values if value and value.lower() != "nan"})


def structure_priority(path: Path) -> tuple[int, int, str]:
    text = str(path)
    alpha_rank = 0 if "AlphaFold" in text else 1
    kcat_rank = 0 if "/kcat/" in text else 1
    return alpha_rank, kcat_rank, text


def index_author_structures(root: Path) -> dict[str, Path]:
    candidates = sorted(root.rglob("*.pdb"), key=structure_priority)
    indexed: dict[str, Path] = {}
    for path in candidates:
        indexed.setdefault(path.stem, path)
    return indexed


def link_one(source: Path, target: Path) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        if target.is_symlink() and Path(os.readlink(target)) == source:
            return "exists"
        if target.stat().st_size > 0:
            return "exists"
        target.unlink()
    target.symlink_to(source)
    return "linked"


def link_ids(ids: list[str], indexed: dict[str, Path], out_dir: Path, scope: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for uniprot in ids:
        source = indexed.get(uniprot)
        if source is None:
            rows.append(
                {
                    "scope": scope,
                    "uniprot_id": uniprot,
                    "status": "missing",
                    "source": "",
                    "target": str(out_dir / f"{uniprot}.pdb"),
                }
            )
            continue
        status = link_one(source.resolve(), out_dir / f"{uniprot}.pdb")
        rows.append(
            {
                "scope": scope,
                "uniprot_id": uniprot,
                "status": status,
                "source": str(source),
                "target": str(out_dir / f"{uniprot}.pdb"),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    metadata = pd.read_csv(args.metadata)
    dekp_train = pd.read_csv(args.dekp_data, sep="\t")
    benchmark_ids = clean_ids(metadata["uniprot_id"])
    train_ids = clean_ids(dekp_train["UniprotID"])
    indexed = index_author_structures(args.author_root)

    rows = []
    rows.extend(link_ids(benchmark_ids, indexed, args.benchmark_dir, "benchmark"))
    rows.extend(link_ids(train_ids, indexed, args.train_dir, "public_kcat_train"))
    report = pd.DataFrame(rows)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.report, index=False)

    missing_ids = sorted(
        report.loc[(report["scope"] == "benchmark") & (report["status"] == "missing"), "uniprot_id"]
        .dropna()
        .astype(str)
        .unique()
    )
    missing = metadata[metadata["uniprot_id"].astype(str).isin(missing_ids)].copy()
    args.missing.parent.mkdir(parents=True, exist_ok=True)
    missing.to_csv(args.missing, index=False)

    summary = report.groupby(["scope", "status"], dropna=False).size().reset_index(name="count")
    print(summary.to_string(index=False))
    print(f"Indexed author structures: {len(indexed)}")
    print(f"Benchmark missing unique UniProt after collection: {len(missing_ids)}")
    print(f"Wrote report: {args.report}")
    print(f"Wrote missing rows: {args.missing}")


if __name__ == "__main__":
    main()
