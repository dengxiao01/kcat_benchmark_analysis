#!/usr/bin/env python3
"""Evaluate early E. coli legacy method outputs on the current benchmark where they overlap."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except ImportError:
    Chem = None


BASE = Path(__file__).resolve().parent.parent
DEFAULT_BENCHMARK = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "legacy_four_methods"
DEFAULT_REPORT_DIR = BASE / "reports" / "tables"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate early four-method outputs on current benchmark overlaps.")
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    return parser.parse_args()


def canonical_smiles(value: object) -> str:
    text = str(value).strip()
    if not text or text == "nan":
        return ""
    if Chem is None:
        return text
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def score_rows(rows: pd.DataFrame) -> pd.DataFrame:
    metrics = []
    groups = [("all", "all", rows)] + [("species", str(k), v) for k, v in rows.groupby("species", sort=True)]
    if "source_database" in rows.columns:
        groups.extend(("source_database", str(k), v) for k, v in rows.groupby("source_database", sort=True))
    for group_type, group, part in groups:
        valid = part[["true_kcat_log10", "prediction_log10"]].replace([np.inf, -np.inf], np.nan).dropna()
        n = len(valid)
        item: dict[str, object] = {"group_type": group_type, "group": group, "n": n}
        if n:
            y_true = valid["true_kcat_log10"]
            y_pred = valid["prediction_log10"]
            err = y_pred - y_true
            ss_res = float(np.square(err).sum())
            ss_tot = float(np.square(y_true - y_true.mean()).sum())
            item.update(
                {
                    "mae_log10": float(err.abs().mean()),
                    "rmse_log10": float(np.sqrt(np.square(err).mean())),
                    "bias_log10": float(err.mean()),
                    "median_abs_error_log10": float(err.abs().median()),
                    "r2_log10": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
                    "pearson_log10": float(y_true.corr(y_pred, method="pearson")) if n > 1 else np.nan,
                    "spearman_log10": float(y_true.corr(y_pred, method="spearman")) if n > 1 else np.nan,
                    "within_0.3_log10_fraction": float((err.abs() <= 0.3).mean()),
                    "within_1.0_log10_fraction": float((err.abs() <= 1.0).mean()),
                }
            )
        metrics.append(item)
    return pd.DataFrame(metrics)


def write_method_outputs(method: str, rows: pd.DataFrame, meta: pd.DataFrame, out_dir: Path, report_dir: Path) -> None:
    method_dir = out_dir / method.lower().replace("/", "_").replace("-", "_")
    method_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(method_dir / f"{method}_legacy_overlap_predictions_evaluated.csv", index=False)
    missing = meta.loc[~meta["legacy_row_id"].isin(rows["legacy_row_id"])].copy()
    missing["legacy_missing_reason"] = np.where(
        missing["species"].astype(str).ne("ecoli"),
        "legacy_ecoli_only",
        "no_legacy_overlap",
    )
    missing.to_csv(method_dir / f"{method}_legacy_overlap_missing_rows.csv", index=False)
    summary = (
        missing.groupby(["species", "legacy_missing_reason"], dropna=False)
        .size()
        .reset_index(name="rows")
        .sort_values(["rows", "species"], ascending=[False, True])
    )
    summary.to_csv(report_dir / f"{method}_legacy_overlap_missing_summary.csv", index=False)
    score_rows(rows).to_csv(report_dir / f"{method}_legacy_overlap_eval_metrics.csv", index=False)


def reaction_level_predictions(path: Path, method: str) -> pd.DataFrame:
    pred = pd.read_csv(path)
    pred["prediction_kcat"] = pd.to_numeric(pred["kcat"], errors="coerce")
    pred = pred[(pred["prediction_kcat"] > 0) & pred["reactions"].notna()].copy()
    pred["prediction_log10"] = np.log10(pred["prediction_kcat"])
    agg = (
        pred.groupby("reactions", as_index=False)["prediction_log10"]
        .median()
        .rename(columns={"reactions": "reaction_id"})
    )
    agg["method"] = method
    return agg


def mtlkp_predictions(path: Path) -> pd.DataFrame:
    pred = pd.read_csv(path)
    pred["canonical_smiles"] = pred["Substrate_SMILES"].map(canonical_smiles)
    pred["prediction_kcat"] = pd.to_numeric(pred["kcat"], errors="coerce")
    pred = pred[(pred["prediction_kcat"] > 0) & (pred["canonical_smiles"] != "")].copy()
    pred["prediction_log10"] = np.log10(pred["prediction_kcat"])
    keys = ["reactions", "genes", "canonical_smiles", "Sequence"]
    agg = pred.groupby(keys, as_index=False)["prediction_log10"].median()
    agg = agg.rename(
        columns={
            "reactions": "reaction_id",
            "genes": "gene_id",
            "Sequence": "sequence",
        }
    )
    agg["method"] = "MTLKP_legacy_overlap"
    return agg


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    bench = pd.read_csv(args.benchmark).reset_index(drop=True)
    bench["legacy_row_id"] = range(len(bench))
    bench["canonical_smiles"] = bench["SMILES"].map(canonical_smiles)

    meta_cols = [
        "legacy_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "canonical_smiles",
        "sequence",
        "true_kcat",
        "true_kcat_log10",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
    ]
    meta = bench[[col for col in meta_cols if col in bench.columns]].copy()
    meta.to_csv(args.out_dir / "legacy_four_methods_benchmark_metadata.csv", index=False)

    # MTLKP has sequence and substrate-level rows, so use the stricter exact-style key.
    mtlkp = mtlkp_predictions(BASE / "reaction_kcat_MW_MTLKP.csv")
    mtlkp_rows = meta[meta["species"].eq("ecoli")].merge(
        mtlkp,
        on=["reaction_id", "gene_id", "canonical_smiles", "sequence"],
        how="inner",
        validate="many_to_one",
    )
    mtlkp_rows["prediction_kcat"] = np.power(10.0, mtlkp_rows["prediction_log10"])
    mtlkp_rows["prediction_column"] = "legacy_exact_sequence_smiles_log10"
    mtlkp_rows["error_log10"] = mtlkp_rows["prediction_log10"] - mtlkp_rows["true_kcat_log10"]
    mtlkp_rows["abs_error_log10"] = mtlkp_rows["error_log10"].abs()
    write_method_outputs("MTLKP", mtlkp_rows, meta, args.out_dir, args.report_dir)

    # TurNuP only has reaction-level early outputs; match by E. coli reaction id.
    turnup = reaction_level_predictions(BASE / "reaction_kcat_MW_TurNup.csv", "TurNuP_legacy_overlap")
    turnup_rows = meta[meta["species"].eq("ecoli")].merge(turnup, on="reaction_id", how="inner", validate="many_to_one")
    turnup_rows["prediction_kcat"] = np.power(10.0, turnup_rows["prediction_log10"])
    turnup_rows["prediction_column"] = "legacy_reaction_level_log10"
    turnup_rows["error_log10"] = turnup_rows["prediction_log10"] - turnup_rows["true_kcat_log10"]
    turnup_rows["abs_error_log10"] = turnup_rows["error_log10"].abs()
    write_method_outputs("TurNuP", turnup_rows, meta, args.out_dir, args.report_dir)

    print("Legacy overlap evaluation complete.")
    print(f"MTLKP overlap rows: {len(mtlkp_rows)}")
    print(f"TurNuP overlap rows: {len(turnup_rows)}")


if __name__ == "__main__":
    main()
