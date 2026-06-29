#!/usr/bin/env python3
"""Prepare PMAK-ready kcat benchmark inputs with complete reaction SMILES."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries_with_sequence_smiles.csv"
DEFAULT_REACTIONS = BASE / "data" / "interim" / "model_reactions.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "pmak"
DEFAULT_REPORT = BASE / "reports" / "tables" / "pmak_eval_readiness.csv"

REQUIRED_COLUMNS = {
    "entry_id",
    "species",
    "reaction_id",
    "substrate_name",
    "SMILES",
    "sequence",
    "true_kcat",
    "true_kcat_log10",
    "source_database",
    "match_level",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create PMAK input, metadata, and truth files.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--entries", type=Path, default=DEFAULT_ENTRIES)
    parser.add_argument("--reactions", type=Path, default=DEFAULT_REACTIONS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--sample-size", type=int, default=20)
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def split_ids(value: object) -> list[str]:
    return [item.strip() for item in str(value).split(";") if item.strip() and item.strip() != "nan"]


def canonical_smiles(smiles: object) -> str:
    text = str(smiles).strip()
    if not text or text == "nan":
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def build_metabolite_smiles(entries: pd.DataFrame) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for _, row in entries.iterrows():
        met_id = str(row.get("substrate_id", "")).strip()
        if not met_id or met_id in lookup:
            continue
        smiles = canonical_smiles(row.get("substrate_smiles", ""))
        if smiles:
            lookup[met_id] = smiles
    return lookup


def reaction_parts(row: pd.Series, metabolite_smiles: dict[str, str]) -> dict[str, object]:
    missing: list[str] = []
    invalid: list[str] = []
    sides: dict[str, list[str]] = {}
    for side, column in [("reactant", "reactant_ids"), ("product", "product_ids")]:
        smiles_parts: list[str] = []
        for met_id in split_ids(row.get(column, "")):
            raw_smiles = metabolite_smiles.get(met_id, "")
            if not raw_smiles:
                missing.append(met_id)
                continue
            smiles = canonical_smiles(raw_smiles)
            if not smiles:
                invalid.append(met_id)
                continue
            smiles_parts.append(smiles)
        sides[side] = smiles_parts

    complete = not missing and not invalid and bool(sides["reactant"]) and bool(sides["product"])
    reactant_smiles = ".".join(sides["reactant"])
    product_smiles = ".".join(sides["product"])
    reaction_smiles = f"{reactant_smiles}>>{product_smiles}" if complete else ""
    return {
        "pmak_reaction_complete": complete,
        "pmak_reactant_smiles": reactant_smiles,
        "pmak_product_smiles": product_smiles,
        "pmak_reaction_smiles": reaction_smiles,
        "pmak_missing_metabolite_ids": ";".join(missing),
        "pmak_invalid_metabolite_ids": ";".join(invalid),
    }


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


def pmak_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pmak_row_id": df["pmak_row_id"],
            "entry_id": df["entry_id"],
            "sequence": clean_text(df["sequence"]),
            "reaction_smiles": clean_text(df["pmak_reaction_smiles"]),
            "reactant_smiles": clean_text(df["pmak_reactant_smiles"]),
            "product_smiles": clean_text(df["pmak_product_smiles"]),
            "substrate_name": clean_text(df["substrate_name"]),
            "substrate_smiles": clean_text(df["SMILES"]),
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "pmak_row_id",
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
        "pmak_reaction_complete",
        "pmak_reactant_smiles",
        "pmak_product_smiles",
        "pmak_reaction_smiles",
        "pmak_missing_metabolite_ids",
        "pmak_invalid_metabolite_ids",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "pmak_row_id",
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
    selected: list[int] = []
    for _, part in groups:
        selected.extend(part.head(per_group).index.tolist())
    selected_set = set(selected)
    remaining = [idx for idx in df.index.tolist() if idx not in selected_set]
    selected.extend(remaining[: max(0, size - len(selected))])
    return df.loc[selected[:size]].sort_index().copy()


def write_report(df: pd.DataFrame, report: Path) -> None:
    rows = []
    for group_name, part in [("all", df)] + list(df.groupby("species", sort=True)):
        complete = part[part["pmak_reaction_complete"]]
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "reaction_complete_rows": len(complete),
                "reaction_incomplete_rows": len(part) - len(complete),
                "unique_complete_reactions": complete[["species", "reaction_id"]].drop_duplicates().shape[0],
                "unique_complete_reaction_smiles": complete["pmak_reaction_smiles"].nunique(),
                "unique_sequences_complete": complete["sequence"].nunique(),
            }
        )
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report, index=False)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)
    entries = pd.read_csv(args.entries)
    reactions = pd.read_csv(args.reactions)
    metabolite_smiles = build_metabolite_smiles(entries)

    df = df.copy().reset_index(drop=True)
    df["pmak_row_id"] = range(len(df))
    df["entry_id"] = clean_text(df["entry_id"])
    df["sequence"] = clean_text(df["sequence"])
    df["SMILES"] = clean_text(df["SMILES"])

    reaction_lookup = reactions.set_index(["species", "reaction_id"])
    reaction_rows = []
    for _, row in df.iterrows():
        key = (row["species"], row["reaction_id"])
        if key not in reaction_lookup.index:
            reaction_rows.append(
                {
                    "pmak_reaction_complete": False,
                    "pmak_reactant_smiles": "",
                    "pmak_product_smiles": "",
                    "pmak_reaction_smiles": "",
                    "pmak_missing_metabolite_ids": "reaction_not_found",
                    "pmak_invalid_metabolite_ids": "",
                }
            )
            continue
        rxn_row = reaction_lookup.loc[key]
        if isinstance(rxn_row, pd.DataFrame):
            rxn_row = rxn_row.iloc[0]
        reaction_rows.append(reaction_parts(rxn_row, metabolite_smiles))
    df = pd.concat([df, pd.DataFrame(reaction_rows)], axis=1)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    valid = df[df["pmak_reaction_complete"]].copy()
    missing = df[~df["pmak_reaction_complete"]].copy()

    pmak_input(valid).to_csv(args.out_dir / "pmak_kcat_input.csv", index=False)
    metadata(df).to_csv(args.out_dir / "pmak_kcat_input_metadata.csv", index=False)
    truth(df).to_csv(args.out_dir / "pmak_kcat_input_truth.csv", index=False)
    metadata(missing).to_csv(args.out_dir / "pmak_missing_reaction_smiles_rows.csv", index=False)

    if args.sample_size > 0:
        sample = balanced_sample(valid, args.sample_size)
        pmak_input(sample).to_csv(args.out_dir / f"pmak_kcat_input_sample{len(sample)}.csv", index=False)
        metadata(sample).to_csv(args.out_dir / f"pmak_kcat_input_sample{len(sample)}_metadata.csv", index=False)
        truth(sample).to_csv(args.out_dir / f"pmak_kcat_input_sample{len(sample)}_truth.csv", index=False)

    write_report(df, args.report)
    print(f"Wrote PMAK input bundle to {args.out_dir}")
    print(f"Rows: {len(df)} total; {len(valid)} with complete reaction SMILES; {len(missing)} incomplete")
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
