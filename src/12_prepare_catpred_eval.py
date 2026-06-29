#!/usr/bin/env python3
"""Prepare CatPred-ready kcat benchmark inputs from the finalized truth table."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "catpred"
DEFAULT_REPORT = BASE / "reports" / "tables" / "catpred_eval_readiness.csv"

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
        description="Create CatPred input, metadata, and truth files for the kcat benchmark."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--sample-size",
        type=int,
        default=100,
        help="Also write a small balanced sample bundle for smoke tests. Use 0 to disable.",
    )
    return parser.parse_args()


def sequence_id(sequence: str) -> str:
    digest = hashlib.sha1(sequence.encode("utf-8")).hexdigest()[:16]
    return f"seq_{digest}"


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    empty = {
        column: int(clean_text(df[column]).eq("").sum())
        for column in ["entry_id", "substrate_name", "SMILES", "sequence"]
    }
    bad = {column: count for column, count in empty.items() if count}
    if bad:
        raise ValueError(f"{path} contains empty required values: {bad}")
    if df["entry_id"].duplicated().any():
        examples = ", ".join(df.loc[df["entry_id"].duplicated(), "entry_id"].head(5))
        raise ValueError(f"entry_id must be unique. Examples: {examples}")


def catpred_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "catpred_row_id": df["catpred_row_id"],
            "entry_id": clean_text(df["entry_id"]),
            "Substrate": clean_text(df["substrate_name"]),
            "SMILES": clean_text(df["SMILES"]),
            "sequence": clean_text(df["sequence"]),
            "pdbpath": clean_text(df["sequence_id"]),
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "catpred_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "sequence_id",
        "SMILES",
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
        "catpred_row_id",
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
    remaining = [idx for idx in df.index.tolist() if idx not in set(selected_indexes)]
    selected_indexes.extend(remaining[: max(0, size - len(selected_indexes))])
    selected_indexes = selected_indexes[:size]
    return df.loc[selected_indexes].sort_index().copy()


def reset_row_ids(df: pd.DataFrame) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["catpred_row_id"] = range(len(out))
    return out


def write_bundle(df: pd.DataFrame, prefix: Path) -> None:
    df = reset_row_ids(df)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    catpred_input(df).to_csv(prefix.with_suffix(".csv"), index=False)
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
    df["sequence"] = clean_text(df["sequence"])
    df["SMILES"] = clean_text(df["SMILES"])
    df["substrate_name"] = clean_text(df["substrate_name"])
    df["sequence_id"] = df["sequence"].map(sequence_id)

    # CatPred uses basename(pdbpath) as a protein-record key, so make the key
    # deterministic from the sequence rather than from a mutable source ID.
    conflicts = df.groupby("sequence_id")["sequence"].nunique()
    if conflicts.gt(1).any():
        raise ValueError("Internal error: sequence hash collision detected.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_bundle(df, args.out_dir / "catpred_kcat_input")
    if args.sample_size > 0:
        sample = balanced_sample(df, args.sample_size)
        write_bundle(sample, args.out_dir / f"catpred_kcat_input_sample{len(sample)}")
    write_report(df, args.report)

    print(f"Wrote CatPred input bundle to {args.out_dir}")
    print(f"Rows: {len(df)}; unique sequences: {df['sequence'].nunique()}; unique SMILES: {df['SMILES'].nunique()}")
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
