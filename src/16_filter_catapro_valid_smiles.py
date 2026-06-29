#!/usr/bin/env python3
"""Filter CataPro input rows to SMILES strings RDKit can parse."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "catapro" / "catapro_kcat_input.csv"
DEFAULT_METADATA = BASE / "data" / "final" / "catapro" / "catapro_kcat_input_metadata.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "catapro" / "catapro_kcat_input_valid_smiles.csv"
DEFAULT_INVALID = BASE / "data" / "final" / "catapro" / "catapro_invalid_smiles_rows.csv"
DEFAULT_REPORT = BASE / "reports" / "tables" / "catapro_valid_smiles_summary.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write a CataPro input file containing only RDKit-valid SMILES.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--invalid", type=Path, default=DEFAULT_INVALID)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def smiles_is_valid(smiles: str) -> bool:
    if not isinstance(smiles, str) or not smiles.strip():
        return False
    return Chem.MolFromSmiles(smiles.strip()) is not None


def main() -> None:
    args = parse_args()
    inp = pd.read_csv(args.input)
    meta = pd.read_csv(args.metadata)
    required = {"catapro_row_id", "Enzyme_id", "type", "sequence", "smiles"}
    missing = sorted(required.difference(inp.columns))
    if missing:
        raise ValueError(f"{args.input} is missing required columns: {', '.join(missing)}")

    inp = inp.copy()
    inp["catapro_smiles_valid"] = inp["smiles"].map(smiles_is_valid)
    valid = inp[inp["catapro_smiles_valid"]].drop(columns=["catapro_smiles_valid"])
    invalid_ids = inp.loc[~inp["catapro_smiles_valid"], "catapro_row_id"]
    invalid = meta[meta["catapro_row_id"].isin(invalid_ids)].copy()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.invalid.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    valid.to_csv(args.output, index=False)
    invalid.to_csv(args.invalid, index=False)
    pd.DataFrame(
        [
            {"group": "all", "rows": len(inp), "valid_smiles_rows": len(valid), "invalid_smiles_rows": len(invalid)}
        ]
    ).to_csv(args.report, index=False)

    print(f"Wrote valid CataPro input: {args.output}")
    print(f"Rows: {len(inp)} total; {len(valid)} RDKit-valid; {len(invalid)} invalid")
    if len(invalid):
        print(f"Invalid rows: {args.invalid}")


if __name__ == "__main__":
    main()
