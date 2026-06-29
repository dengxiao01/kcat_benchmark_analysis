#!/usr/bin/env python3
"""Prepare GO-HKP functional-assignment kcat predictions for the benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_GO_ASSIGNMENT = (
    BASE
    / "external_methods"
    / "GO-HKP"
    / "analysis"
    / "DeepGO-SE"
    / "iML1515R"
    / "go_kcat_mean_parent_process_Total_median.json"
)
DEFAULT_OUT_DIR = BASE / "data" / "final" / "go_hkp"
DEFAULT_REPORT = BASE / "reports" / "tables" / "go_hkp_eval_readiness.csv"

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
        description="Create GO-HKP input, metadata, truth, and prediction files."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--go-assignment", type=Path, default=DEFAULT_GO_ASSIGNMENT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    if df["entry_id"].duplicated().any():
        examples = ", ".join(df.loc[df["entry_id"].duplicated(), "entry_id"].head(5))
        raise ValueError(f"entry_id must be unique. Examples: {examples}")


def load_assignments(path: Path) -> dict[str, float]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    assignments: dict[str, float] = {}
    for reaction_id, value in raw.items():
        try:
            kcat = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(kcat) and kcat > 0:
            assignments[str(reaction_id)] = kcat
    return assignments


def enrich(df: pd.DataFrame, assignments: dict[str, float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_meta = df.reset_index(drop=True).copy()
    all_meta["entry_id"] = clean_text(all_meta["entry_id"])
    all_meta["reaction_id"] = clean_text(all_meta["reaction_id"])
    all_meta["go_hkp_prediction_status"] = "missing"
    all_meta["go_hkp_missing_reason"] = "species_without_local_deepgo_se_assignment"
    all_meta["go_hkp_assignment_source"] = ""
    all_meta["go_hkp_assignment_kcat"] = np.nan

    ecoli_mask = all_meta["species"].eq("ecoli")
    matched = ecoli_mask & all_meta["reaction_id"].isin(assignments)
    all_meta.loc[ecoli_mask, "go_hkp_missing_reason"] = "missing_reaction_assignment"
    all_meta.loc[matched, "go_hkp_prediction_status"] = "ready"
    all_meta.loc[matched, "go_hkp_missing_reason"] = ""
    all_meta.loc[matched, "go_hkp_assignment_source"] = "GO-HKP DeepGO-SE iML1515R Total median"
    all_meta.loc[matched, "go_hkp_assignment_kcat"] = all_meta.loc[matched, "reaction_id"].map(assignments)
    all_meta["go_hkp_row_id"] = pd.NA

    ready = all_meta[all_meta["go_hkp_prediction_status"].eq("ready")].copy().reset_index(drop=True)
    ready.insert(0, "go_hkp_row_id_ready", range(len(ready)))
    ready["go_hkp_row_id"] = ready["go_hkp_row_id_ready"]
    ready = ready.drop(columns=["go_hkp_row_id_ready"])

    all_meta = all_meta.drop(columns=["go_hkp_row_id"]).merge(
        ready[["entry_id", "go_hkp_row_id"]],
        on="entry_id",
        how="left",
        validate="one_to_one",
    )
    return all_meta, ready


def go_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "go_hkp_row_id": df["go_hkp_row_id"],
            "entry_id": df["entry_id"],
            "species": df["species"],
            "reaction_id": df["reaction_id"],
            "gene_id": df["gene_id"],
            "uniprot_id": df["uniprot_id"],
            "assignment_rule": "reaction_id -> GO-HKP DeepGO-SE iML1515R Total median kcat",
        }
    )


def go_output(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "go_hkp_row_id": df["go_hkp_row_id"],
            "entry_id": df["entry_id"],
            "reaction_id": df["reaction_id"],
            "prediction_kcat": df["go_hkp_assignment_kcat"].astype(float),
            "prediction_log10": np.log10(df["go_hkp_assignment_kcat"].astype(float)),
            "go_hkp_assignment_source": df["go_hkp_assignment_source"],
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "go_hkp_row_id",
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
        "go_hkp_prediction_status",
        "go_hkp_missing_reason",
        "go_hkp_assignment_source",
        "go_hkp_assignment_kcat",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "go_hkp_row_id",
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


def readiness(all_meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group_name, part in [("all", all_meta)] + list(all_meta.groupby("species", sort=True)):
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "ready_rows": int(part["go_hkp_prediction_status"].eq("ready").sum()),
                "missing_rows": int((~part["go_hkp_prediction_status"].eq("ready")).sum()),
                "unique_reactions": part["reaction_id"].nunique(),
                "ready_unique_reactions": part.loc[
                    part["go_hkp_prediction_status"].eq("ready"), "reaction_id"
                ].nunique(),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)
    assignments = load_assignments(args.go_assignment)
    all_meta, ready = enrich(df, assignments)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    go_input(ready).to_csv(args.out_dir / "go_hkp_kcat_input.csv", index=False)
    go_output(ready).to_csv(args.out_dir / "go_hkp_kcat_input_output.csv", index=False)
    metadata(ready).to_csv(args.out_dir / "go_hkp_kcat_input_metadata.csv", index=False)
    metadata(all_meta).to_csv(args.out_dir / "go_hkp_kcat_all_metadata.csv", index=False)
    truth(ready).to_csv(args.out_dir / "go_hkp_kcat_input_truth.csv", index=False)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    readiness(all_meta).to_csv(args.report, index=False)

    print(f"GO-HKP assignments loaded: {len(assignments)} reactions")
    print(f"Ready rows: {len(ready)}")
    print(f"Missing rows: {len(all_meta) - len(ready)}")
    print(f"Wrote GO-HKP files to: {args.out_dir}")


if __name__ == "__main__":
    main()
