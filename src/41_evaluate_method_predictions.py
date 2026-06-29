#!/usr/bin/env python3
"""Generic evaluator for kcat prediction tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
LOG_COLUMNS = ["prediction_log10", "pred_log10_kcat", "predicted_log10_kcat", "log10_kcat"]
LINEAR_COLUMNS = ["prediction_kcat", "pred_kcat", "predicted_kcat", "kcat", "y_pred", "prediction"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score method predictions against benchmark truth.")
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--all-metadata", type=Path, default=None)
    parser.add_argument("--out-rows", type=Path, required=True)
    parser.add_argument("--out-metrics", type=Path, required=True)
    parser.add_argument("--out-missing", type=Path, required=True)
    parser.add_argument("--out-missing-summary", type=Path, required=True)
    parser.add_argument("--prediction-column", default="")
    return parser.parse_args()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def detect_prediction(pred: pd.DataFrame, explicit: str = "") -> tuple[pd.Series, str]:
    if explicit:
        if explicit not in pred.columns:
            raise ValueError(f"Prediction column {explicit!r} not found.")
        values = numeric(pred[explicit])
        if explicit in LINEAR_COLUMNS or ("kcat" in explicit.lower() and "log" not in explicit.lower()):
            return np.log10(values.where(values > 0)), explicit
        return values, explicit
    for column in LOG_COLUMNS:
        if column in pred.columns:
            return numeric(pred[column]), column
    for column in LINEAR_COLUMNS:
        if column in pred.columns:
            values = numeric(pred[column])
            return np.log10(values.where(values > 0)), column
    raise ValueError(f"No recognized prediction column found. Columns: {list(pred.columns)}")


def merge_predictions(pred: pd.DataFrame, meta: pd.DataFrame, pred_log10: pd.Series, pred_col: str) -> pd.DataFrame:
    pred = pred.copy()
    meta = meta.copy()
    pred["prediction_log10"] = pred_log10
    pred["prediction_column"] = pred_col

    optional = [col for col in pred.columns if col.startswith("prediction_") and col not in {"prediction_log10"}]
    keep_base = ["prediction_log10", "prediction_column"] + optional

    row_ids = [col for col in meta.columns if col.endswith("_row_id") and col in pred.columns]
    if row_ids:
        row_id = row_ids[0]
        pred[row_id] = pd.to_numeric(pred[row_id], errors="coerce").astype("Int64")
        meta[row_id] = pd.to_numeric(meta[row_id], errors="coerce").astype("Int64")
        rows = meta.merge(pred[[row_id] + keep_base], on=row_id, how="inner", validate="one_to_one")
    elif "entry_id" in pred.columns and "entry_id" in meta.columns:
        rows = meta.merge(pred[["entry_id"] + keep_base], on="entry_id", how="inner", validate="one_to_one")
    else:
        if len(pred) != len(meta):
            raise ValueError("Prediction rows do not match metadata rows, and no row id or entry_id can be used.")
        rows = pd.concat([meta.reset_index(drop=True), pred[keep_base].reset_index(drop=True)], axis=1)

    rows["true_kcat_log10"] = numeric(rows["true_kcat_log10"])
    rows["true_kcat"] = numeric(rows["true_kcat"])
    rows["prediction_kcat"] = np.power(10.0, rows["prediction_log10"])
    rows["error_log10"] = rows["prediction_log10"] - rows["true_kcat_log10"]
    rows["abs_error_log10"] = rows["error_log10"].abs()
    rows = rows.replace([np.inf, -np.inf], np.nan).dropna(subset=["true_kcat_log10", "prediction_log10"]).copy()
    return rows


def missing_rows(all_meta: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    row_ids = [col for col in all_meta.columns if col.endswith("_row_id") and col in rows.columns]
    if row_ids:
        row_id = row_ids[0]
        return all_meta.loc[~all_meta[row_id].isin(rows[row_id])].copy()
    if "entry_id" in all_meta.columns and "entry_id" in rows.columns:
        return all_meta.loc[~all_meta["entry_id"].isin(rows["entry_id"])].copy()
    return all_meta.iloc[0:0].copy()


def score_group(group_type: str, group_name: str, df: pd.DataFrame) -> dict[str, float | int | str]:
    valid = df[["true_kcat_log10", "prediction_log10"]].replace([np.inf, -np.inf], np.nan).dropna()
    n = len(valid)
    row: dict[str, float | int | str] = {"group_type": group_type, "group": group_name, "n": n}
    if n == 0:
        return row
    y_true = valid["true_kcat_log10"]
    y_pred = valid["prediction_log10"]
    err = y_pred - y_true
    ss_res = float(np.square(err).sum())
    ss_tot = float(np.square(y_true - y_true.mean()).sum())
    row.update(
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
    return row


def build_metrics(rows: pd.DataFrame) -> pd.DataFrame:
    metrics = [score_group("all", "all", rows)]
    for species, part in rows.groupby("species", sort=True):
        metrics.append(score_group("species", str(species), part))
    if "source_database" in rows.columns:
        for source, part in rows.groupby("source_database", sort=True):
            metrics.append(score_group("source_database", str(source), part))
    return pd.DataFrame(metrics)


def reason_for_missing(row: pd.Series) -> str:
    smiles = str(row.get("SMILES", row.get("smiles", ""))).strip()
    if not smiles or smiles == "nan":
        return "empty_smiles"
    for col in [
        "prediction_status",
        "legacy_missing_reason",
        "unikp_sequence_feature_status",
        "kinform_unavailable_reason",
        "mtlkp_missing_reason",
        "turnup_missing_reason",
        "go_hkp_missing_reason",
    ]:
        value = str(row.get(col, "")).strip()
        if value and value.lower() not in {"nan", "success", "true"}:
            return value
    if "." in smiles:
        return "invalid_or_multicomponent_smiles"
    return "unpredicted"


def write_missing_summary(missing: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if missing.empty:
        pd.DataFrame(columns=["species", "reason", "rows"]).to_csv(out_path, index=False)
        return
    missing = missing.copy()
    missing["reason"] = missing.apply(reason_for_missing, axis=1)
    summary = (
        missing.groupby(["species", "reason"], dropna=False)
        .size()
        .reset_index(name="rows")
        .sort_values(["rows", "species"], ascending=[False, True])
    )
    summary.to_csv(out_path, index=False)


def main() -> None:
    args = parse_args()
    pred = pd.read_csv(args.predictions)
    meta = pd.read_csv(args.metadata)
    all_meta = pd.read_csv(args.all_metadata) if args.all_metadata else meta.copy()
    pred_log10, pred_col = detect_prediction(pred, args.prediction_column)
    rows = merge_predictions(pred, meta, pred_log10, pred_col)
    missing = missing_rows(all_meta, rows)
    metrics = build_metrics(rows)

    args.out_rows.parent.mkdir(parents=True, exist_ok=True)
    args.out_metrics.parent.mkdir(parents=True, exist_ok=True)
    args.out_missing.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.out_rows, index=False)
    metrics.to_csv(args.out_metrics, index=False)
    missing.to_csv(args.out_missing, index=False)
    write_missing_summary(missing, args.out_missing_summary)

    overall = metrics[(metrics["group_type"] == "all") & (metrics["group"] == "all")].iloc[0]
    print(f"Prediction column: {pred_col}")
    print(f"Evaluated rows: {int(overall['n'])}")
    print(f"Missing rows: {len(missing)}")
    print(f"MAE log10: {overall.get('mae_log10', np.nan):.4g}")
    print(f"RMSE log10: {overall.get('rmse_log10', np.nan):.4g}")
    print(f"Metrics: {args.out_metrics}")


if __name__ == "__main__":
    main()
