#!/usr/bin/env python3
"""Prepare KinForm-ready kcat benchmark inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

try:
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # Keep the script usable in lightweight Python envs.
    Chem = None


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "kinform"
DEFAULT_REPORT = BASE / "reports" / "tables" / "kinform_eval_readiness.csv"

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
        description="Create KinForm JSON input, metadata, and truth files for the kcat benchmark."
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


def kinform_model_sequence(sequence: str) -> str:
    """Match KinForm's kcat sequence truncation rule."""
    return sequence[:749] + sequence[-749:] if len(sequence) > 1499 else sequence


def canonical_smiles(smiles: object) -> tuple[bool, str]:
    text = str(smiles).strip()
    if not text or text == "nan":
        return False, ""
    if Chem is None:
        return True, text
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return False, ""
    return True, Chem.MolToSmiles(mol, canonical=True)


def enrich(df: pd.DataFrame) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["kinform_row_id"] = range(len(out))
    out["entry_id"] = clean_text(out["entry_id"])
    out["sequence"] = clean_text(out["sequence"])
    out["SMILES"] = clean_text(out["SMILES"])
    out["kinform_model_sequence"] = out["sequence"].map(kinform_model_sequence)
    out["kinform_sequence_truncated"] = out["sequence"] != out["kinform_model_sequence"]

    smiles_status = out["SMILES"].map(canonical_smiles)
    out["kinform_smiles_valid"] = [valid for valid, _ in smiles_status]
    out["kinform_canonical_smiles"] = [canonical for _, canonical in smiles_status]
    out["kinform_smiles_validation"] = "rdkit" if Chem is not None else "not_checked_rdkit_unavailable"
    return out


def input_rows(df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        rows.append(
            {
                "kinform_row_id": int(row["kinform_row_id"]),
                "entry_id": row["entry_id"],
                "sequence": row["sequence"],
                "smiles": row["kinform_canonical_smiles"] or row["SMILES"],
                "value": float(row["true_kcat"]),
            }
        )
    return rows


def write_json(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, ensure_ascii=False, indent=2)


def input_csv(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinform_row_id": df["kinform_row_id"],
            "entry_id": df["entry_id"],
            "sequence": df["sequence"],
            "smiles": df["kinform_canonical_smiles"].where(
                df["kinform_canonical_smiles"].astype(str).str.len() > 0, df["SMILES"]
            ),
            "value": df["true_kcat"],
            "true_kcat_log10": df["true_kcat_log10"],
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "kinform_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "kinform_canonical_smiles",
        "kinform_smiles_valid",
        "kinform_smiles_validation",
        "sequence",
        "kinform_model_sequence",
        "kinform_sequence_truncated",
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
        "kinform_row_id",
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


def write_report(df: pd.DataFrame, report: Path) -> None:
    rows = []
    for group_name, part in [("all", df)] + list(df.groupby("species", sort=True)):
        source_counts = part["source_database"].value_counts(dropna=False).to_dict()
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "valid_smiles_rows": int(part["kinform_smiles_valid"].sum()),
                "invalid_smiles_rows": int((~part["kinform_smiles_valid"]).sum()),
                "truncated_sequence_rows": int(part["kinform_sequence_truncated"].sum()),
                "unique_sequences": part["sequence"].nunique(),
                "unique_model_sequences": part["kinform_model_sequence"].nunique(),
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
    valid = df[df["kinform_smiles_valid"]].copy()
    invalid = df[~df["kinform_smiles_valid"]].copy()

    input_csv(df).to_csv(args.out_dir / "kinform_kcat_input.csv", index=False)
    input_csv(valid).to_csv(args.out_dir / "kinform_kcat_input_valid_smiles.csv", index=False)
    metadata(df).to_csv(args.out_dir / "kinform_kcat_input_metadata.csv", index=False)
    truth(df).to_csv(args.out_dir / "kinform_kcat_input_truth.csv", index=False)
    metadata(invalid).to_csv(args.out_dir / "kinform_invalid_smiles_rows.csv", index=False)
    write_json(input_rows(df), args.out_dir / "kinform_kcat_input.json")
    write_json(input_rows(valid), args.out_dir / "kinform_kcat_input_valid_smiles.json")

    if args.sample_size > 0:
        sample = balanced_sample(valid, args.sample_size)
        sample_prefix = args.out_dir / f"kinform_kcat_input_sample{len(sample)}"
        input_csv(sample).to_csv(sample_prefix.with_suffix(".csv"), index=False)
        metadata(sample).to_csv(args.out_dir / f"{sample_prefix.name}_metadata.csv", index=False)
        truth(sample).to_csv(args.out_dir / f"{sample_prefix.name}_truth.csv", index=False)
        write_json(input_rows(sample), sample_prefix.with_suffix(".json"))

    write_report(df, args.report)
    print(f"Wrote KinForm input bundle to {args.out_dir}")
    print(f"Rows: {len(df)} total; {len(valid)} valid SMILES; {len(invalid)} invalid SMILES")
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
