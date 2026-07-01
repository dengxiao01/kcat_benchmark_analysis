#!/usr/bin/env python3
"""Prepare unified benchmark inputs for TurNuP."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
BENCHMARK = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
PMAK_INPUT = BASE / "data" / "final" / "pmak" / "pmak_kcat_input.csv"
PMAK_META = BASE / "data" / "final" / "pmak" / "pmak_kcat_input_metadata.csv"
TURNUP_DIR = BASE / "data" / "final" / "turnup"



def smiles_side_to_turnup(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().replace(".", ";")


def prepare_turnup() -> None:
    TURNUP_DIR.mkdir(parents=True, exist_ok=True)
    bench = pd.read_csv(BENCHMARK)
    pmak_input = pd.read_csv(PMAK_INPUT)
    pmak_meta = pd.read_csv(PMAK_META)

    ready = pmak_input.copy().reset_index(drop=True)
    ready.insert(0, "turnup_row_id", range(len(ready)))
    ready["substrates"] = ready["reactant_smiles"].map(smiles_side_to_turnup)
    ready["products"] = ready["product_smiles"].map(smiles_side_to_turnup)
    ready["enzyme"] = ready["sequence"].fillna("").astype(str).str.upper()

    input_df = ready[
        [
            "turnup_row_id",
            "entry_id",
            "substrates",
            "products",
            "enzyme",
            "reactant_smiles",
            "product_smiles",
            "reaction_smiles",
            "sequence",
            "substrate_name",
            "substrate_smiles",
        ]
    ].copy()

    valid_meta = pmak_meta.merge(
        ready[["entry_id", "turnup_row_id", "substrates", "products"]],
        on="entry_id",
        how="inner",
        validate="one_to_one",
    )
    valid_meta["turnup_prediction_status"] = "ready"
    valid_meta["turnup_missing_reason"] = ""

    all_meta = bench.merge(
        valid_meta[
            [
                "entry_id",
                "turnup_row_id",
                "turnup_prediction_status",
                "turnup_missing_reason",
                "substrates",
                "products",
                "pmak_reaction_complete",
                "pmak_missing_metabolite_ids",
                "pmak_invalid_metabolite_ids",
            ]
        ],
        on="entry_id",
        how="left",
    )
    all_meta["turnup_prediction_status"] = all_meta["turnup_prediction_status"].fillna("missing")
    all_meta["turnup_missing_reason"] = all_meta["turnup_missing_reason"].fillna("missing_reaction_smiles")
    all_meta.loc[
        all_meta["turnup_prediction_status"].eq("ready"),
        "turnup_missing_reason",
    ] = ""

    truth = valid_meta[["turnup_row_id", "entry_id", "true_kcat", "true_kcat_log10"]].copy()

    input_df.to_csv(TURNUP_DIR / "turnup_kcat_input.csv", index=False)
    valid_meta.to_csv(TURNUP_DIR / "turnup_kcat_input_metadata.csv", index=False)
    all_meta.to_csv(TURNUP_DIR / "turnup_kcat_all_metadata.csv", index=False)
    truth.to_csv(TURNUP_DIR / "turnup_kcat_input_truth.csv", index=False)

    missing = all_meta[~all_meta["turnup_prediction_status"].eq("ready")].copy()
    missing.to_csv(TURNUP_DIR / "turnup_missing_reaction_smiles_rows.csv", index=False)
    print(f"TurNuP ready rows: {len(input_df)}")
    print(f"TurNuP missing rows: {len(missing)}")


def main() -> None:
    prepare_turnup()


if __name__ == "__main__":
    main()
