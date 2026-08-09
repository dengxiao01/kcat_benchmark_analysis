#!/usr/bin/env python3
"""Prepare KcatNet-ready kcat benchmark inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "kcatnet"
DEFAULT_REPORT = BASE / "reports" / "tables" / "kcatnet_eval_readiness.csv"

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
        description="Create KcatNet input, metadata, and truth files for the kcat benchmark."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--sample-size", type=int, default=20)
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


def canonical_smiles(smiles: object) -> tuple[bool, str]:
    text = str(smiles).strip()
    if not text or text == "nan":
        return False, ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return False, ""
    return True, Chem.MolToSmiles(mol, canonical=True)


def kcatnet_model_sequence(sequence: str) -> str:
    return sequence[:500] + sequence[-500:] if len(sequence) > 1000 else sequence


def enrich(df: pd.DataFrame) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["kcatnet_row_id"] = range(len(out))
    out["entry_id"] = clean_text(out["entry_id"])
    out["sequence"] = clean_text(out["sequence"])
    out["SMILES"] = clean_text(out["SMILES"])
    out["kcatnet_model_sequence"] = out["sequence"].map(kcatnet_model_sequence)
    out["kcatnet_sequence_truncated"] = out["sequence"] != out["kcatnet_model_sequence"]

    smiles_status = out["SMILES"].map(canonical_smiles)
    out["kcatnet_smiles_valid"] = [valid for valid, _ in smiles_status]
    out["kcatnet_canonical_smiles"] = [canonical for _, canonical in smiles_status]
    return out


def kcatnet_input(df: pd.DataFrame) -> pd.DataFrame:
    smiles = df["kcatnet_canonical_smiles"].where(
        df["kcatnet_canonical_smiles"].astype(str).str.len() > 0, df["SMILES"]
    )
    return pd.DataFrame(
        {
            "kcatnet_row_id": df["kcatnet_row_id"],
            "entry_id": df["entry_id"],
            "Pro_seq": df["sequence"],
            "Smile": smiles,
            "true_kcat": df["true_kcat"],
            "true_kcat_log10": df["true_kcat_log10"],
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "kcatnet_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "kcatnet_canonical_smiles",
        "kcatnet_smiles_valid",
        "sequence",
        "kcatnet_model_sequence",
        "kcatnet_sequence_truncated",
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
        "kcatnet_row_id",
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
    return df.loc[selected_indexes[:size]].sort_index().copy()


def write_bundle(df: pd.DataFrame, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    kcatnet_input(df).to_csv(prefix.with_suffix(".csv"), index=False)
    kcatnet_input(df).to_excel(prefix.with_suffix(".xlsx"), index=False)
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
                "valid_smiles_rows": int(part["kcatnet_smiles_valid"].sum()),
                "invalid_smiles_rows": int((~part["kcatnet_smiles_valid"]).sum()),
                "truncated_sequence_rows": int(part["kcatnet_sequence_truncated"].sum()),
                "unique_sequences": part["sequence"].nunique(),
                "unique_model_sequences": part["kcatnet_model_sequence"].nunique(),
                "unique_smiles": part["SMILES"].nunique(),
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
    df = enrich(df)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    valid = df[df["kcatnet_smiles_valid"]].copy()
    invalid = df[~df["kcatnet_smiles_valid"]].copy()

    write_bundle(df, args.out_dir / "kcatnet_kcat_input")
    write_bundle(valid, args.out_dir / "kcatnet_kcat_input_valid_smiles")
    metadata(invalid).to_csv(args.out_dir / "kcatnet_invalid_smiles_rows.csv", index=False)

    if args.sample_size > 0:
        sample = balanced_sample(valid, args.sample_size)
        write_bundle(sample, args.out_dir / f"kcatnet_kcat_input_sample{len(sample)}")

    write_report(df, args.report)
    print(f"Wrote KcatNet input bundle to {args.out_dir}")
    print(f"Rows: {len(df)} total; {len(valid)} valid SMILES; {len(invalid)} invalid SMILES")
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
