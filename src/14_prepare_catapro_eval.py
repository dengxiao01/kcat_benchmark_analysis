#!/usr/bin/env python3
"""Prepare CataPro-ready kcat benchmark inputs from the finalized truth table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "catapro"
DEFAULT_REPORT = BASE / "reports" / "tables" / "catapro_eval_readiness.csv"

REQUIRED_COLUMNS = {
    "entry_id",
    "species",
    "reaction_id",
    "gene_id",
    "uniprot_id",
    "ec_number",
    "substrate_name",
    "SMILES",
    "sequence",
    "true_kcat",
    "true_kcat_log10",
    "unit",
    "source_database",
    "match_level",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create CataPro input, metadata, and truth files for the kcat benchmark."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--sample-size",
        type=int,
        default=20,
        help="Also write a small balanced sample bundle for smoke tests. Use 0 to disable.",
    )
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    empty = {
        column: int(clean_text(df[column]).eq("").sum())
        for column in ["entry_id", "SMILES", "sequence"]
    }
    bad = {column: count for column, count in empty.items() if count}
    if bad:
        raise ValueError(f"{path} contains empty required values: {bad}")
    if df["entry_id"].duplicated().any():
        examples = ", ".join(df.loc[df["entry_id"].duplicated(), "entry_id"].head(5))
        raise ValueError(f"entry_id must be unique. Examples: {examples}")


def reset_row_ids(df: pd.DataFrame) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["catapro_row_id"] = range(len(out))
    return out


def catapro_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "catapro_row_id": df["catapro_row_id"],
            "Enzyme_id": df["entry_id"],
            "type": "wild",
            "sequence": clean_text(df["sequence"]),
            "smiles": clean_text(df["SMILES"]),
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "catapro_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "sequence",
        "true_kcat",
        "true_kcat_log10",
        "unit",
        "pH",
        "temperature_c",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
        "enzyme_complex_type",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "catapro_row_id",
        "entry_id",
        "species",
        "true_kcat",
        "true_kcat_log10",
        "unit",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def balanced_sample(df: pd.DataFrame, size: int) -> pd.DataFrame:
    if size <= 0 or len(df) <= size:
        return df.copy()
    groups = [(name, part.copy()) for name, part in df.groupby("species", sort=True)]
    per_group = max(1, size // max(len(groups), 1))
    selected_indexes: list[int] = []
    for _, part in groups:
        selected_indexes.extend(part.head(per_group).index.tolist())
    selected_set = set(selected_indexes)
    remaining = [idx for idx in df.index.tolist() if idx not in selected_set]
    selected_indexes.extend(remaining[: max(0, size - len(selected_indexes))])
    selected_indexes = selected_indexes[:size]
    return df.loc[selected_indexes].sort_index().copy()


def write_bundle(df: pd.DataFrame, prefix: Path) -> None:
    df = reset_row_ids(df)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    catapro_input(df).to_csv(prefix.with_suffix(".csv"), index=False)
    metadata(df).to_csv(prefix.parent / f"{prefix.name}_metadata.csv", index=False)
    truth(df).to_csv(prefix.parent / f"{prefix.name}_truth.csv", index=False)


def write_report(df: pd.DataFrame, report: Path) -> None:
    rows = []
    for group_name, part in [("all", df)] + list(df.groupby("species", sort=True)):
        source_counts = part["source_database"].value_counts(dropna=False).to_dict()
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "unique_sequences": part["sequence"].nunique(),
                "unique_smiles": part["SMILES"].nunique(),
                "unique_substrates": part["substrate_name"].nunique(),
                "brenda_rows": int(source_counts.get("BRENDA", 0)),
                "sabiork_rows": int(source_counts.get("SABIO-RK", 0)),
                "brenda_sabiork_rows": int(source_counts.get("BRENDA;SABIO-RK", 0)),
            }
        )
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report, index=False)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)

    df = df.copy()
    df["entry_id"] = clean_text(df["entry_id"])
    df["sequence"] = clean_text(df["sequence"])
    df["SMILES"] = clean_text(df["SMILES"])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_bundle(df, args.out_dir / "catapro_kcat_input")
    if args.sample_size > 0:
        sample = balanced_sample(df, args.sample_size)
        write_bundle(sample, args.out_dir / f"catapro_kcat_input_sample{len(sample)}")
    write_report(df, args.report)

    print(f"Wrote CataPro input bundle to {args.out_dir}")
    print(f"Rows: {len(df)}; unique sequences: {df['sequence'].nunique()}; unique SMILES: {df['SMILES'].nunique()}")
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
