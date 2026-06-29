#!/usr/bin/env python3
"""Prepare SELFprot-ready benchmark inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "selfprot"
DEFAULT_REPORT = BASE / "reports" / "tables" / "selfprot_eval_readiness.csv"


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
    parser = argparse.ArgumentParser(description="Create SELFprot input, metadata, truth, and readiness report.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--sample-size", type=int, default=20)
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def canonical_smiles(smiles: object) -> tuple[bool, str]:
    text = str(smiles).strip()
    if not text or text.lower() == "nan":
        return False, ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return False, ""
    return True, Chem.MolToSmiles(mol, canonical=True)


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    if df["entry_id"].duplicated().any():
        examples = ", ".join(df.loc[df["entry_id"].duplicated(), "entry_id"].head(5))
        raise ValueError(f"entry_id must be unique. Examples: {examples}")


def enrich(df: pd.DataFrame) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["selfprot_row_id"] = range(len(out))
    for column in ["entry_id", "sequence", "SMILES", "uniprot_id", "species"]:
        out[column] = clean_text(out[column])
    statuses = out["SMILES"].map(canonical_smiles)
    out["selfprot_smiles_valid"] = [valid for valid, _ in statuses]
    out["selfprot_canonical_smiles"] = [canonical for _, canonical in statuses]
    out["selfprot_sequence_length"] = out["sequence"].str.len()
    out["selfprot_smiles_length"] = out["selfprot_canonical_smiles"].str.len()
    return out


def selfprot_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "selfprot_row_id": df["selfprot_row_id"],
            "entry_id": df["entry_id"],
            "sequence": df["sequence"],
            "smiles": df["selfprot_canonical_smiles"].where(
                df["selfprot_canonical_smiles"].astype(str).str.len() > 0, df["SMILES"]
            ),
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "selfprot_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "selfprot_canonical_smiles",
        "selfprot_smiles_valid",
        "sequence",
        "selfprot_sequence_length",
        "selfprot_smiles_length",
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
        "selfprot_row_id",
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
    if len(selected_indexes) < size:
        remainder = df.loc[~df.index.isin(selected_set)].head(size - len(selected_indexes))
        selected_indexes.extend(remainder.index.tolist())
    return df.loc[selected_indexes[:size]].copy()


def readiness(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, part in [("all", df), *[(str(k), v) for k, v in df.groupby("species", sort=True)]]:
        valid = part[part["selfprot_smiles_valid"]].copy()
        rows.append(
            {
                "group": group,
                "rows": len(part),
                "valid_smiles_rows": int(part["selfprot_smiles_valid"].sum()),
                "invalid_smiles_rows": int((~part["selfprot_smiles_valid"]).sum()),
                "predictable_rows": len(valid),
                "unique_sequences": valid["sequence"].nunique(),
                "unique_smiles": valid["selfprot_canonical_smiles"].nunique(),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)
    enriched = enrich(df)
    valid = enriched[enriched["selfprot_smiles_valid"]].copy()
    invalid = enriched[~enriched["selfprot_smiles_valid"]].copy()
    sample = balanced_sample(valid, args.sample_size)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    selfprot_input(enriched).to_csv(args.out_dir / "selfprot_kcat_input.csv", index=False)
    selfprot_input(valid).to_csv(args.out_dir / "selfprot_kcat_input_valid_smiles.csv", index=False)
    selfprot_input(sample).to_csv(args.out_dir / "selfprot_kcat_input_sample20.csv", index=False)
    metadata(enriched).to_csv(args.out_dir / "selfprot_kcat_input_metadata.csv", index=False)
    metadata(valid).to_csv(args.out_dir / "selfprot_kcat_input_valid_smiles_metadata.csv", index=False)
    metadata(sample).to_csv(args.out_dir / "selfprot_kcat_input_sample20_metadata.csv", index=False)
    truth(enriched).to_csv(args.out_dir / "selfprot_kcat_input_truth.csv", index=False)
    truth(valid).to_csv(args.out_dir / "selfprot_kcat_input_valid_smiles_truth.csv", index=False)
    truth(sample).to_csv(args.out_dir / "selfprot_kcat_input_sample20_truth.csv", index=False)
    metadata(invalid).to_csv(args.out_dir / "selfprot_invalid_smiles_rows.csv", index=False)

    report = readiness(enriched)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.report, index=False)
    print(f"Wrote SELFprot input bundle to {args.out_dir}")
    print(
        f"Rows: {len(enriched)} total; {len(valid)} valid SMILES; "
        f"{len(invalid)} invalid SMILES"
    )
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
