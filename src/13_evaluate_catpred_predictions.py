#!/usr/bin/env python3
"""Evaluate CatPred kcat predictions against the finalized benchmark truth."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_PREDICTIONS = BASE / "data" / "final" / "catpred" / "catpred_kcat_input_output.csv"
DEFAULT_METADATA = BASE / "data" / "final" / "catpred" / "catpred_kcat_input_metadata.csv"
DEFAULT_OUT_ROWS = BASE / "data" / "final" / "catpred" / "catpred_kcat_predictions_evaluated.csv"
DEFAULT_OUT_METRICS = BASE / "reports" / "tables" / "catpred_eval_metrics.csv"

LOG_COLUMNS = [
    "Prediction_log10",
    "prediction_log10",
    "log10kcat_max",
    "log10_kcat",
    "predicted_log10_kcat",
]
LINEAR_COLUMNS = [
    "Prediction_(s^(-1))",
    "Prediction",
    "prediction",
    "predicted_kcat",
    "kcat",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score CatPred output against benchmark truth.")
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--out-rows", type=Path, default=DEFAULT_OUT_ROWS)
    parser.add_argument("--out-metrics", type=Path, default=DEFAULT_OUT_METRICS)
    parser.add_argument(
        "--prediction-column",
        default="",
        help="Optional explicit prediction column. If omitted, common CatPred column names are detected.",
    )
    return parser.parse_args()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def detect_prediction(pred: pd.DataFrame, explicit: str = "") -> tuple[pd.Series, str]:
    if explicit:
        if explicit not in pred.columns:
            raise ValueError(f"Prediction column {explicit!r} not found.")
        values = numeric(pred[explicit])
        if explicit in LINEAR_COLUMNS or "s^(-1)" in explicit or "kcat" in explicit.lower() and "log" not in explicit.lower():
            return np.log10(values.where(values > 0)), explicit
        return values, explicit

    for column in LOG_COLUMNS:
        if column in pred.columns:
            return numeric(pred[column]), column
    for column in LINEAR_COLUMNS:
        if column in pred.columns:
            values = numeric(pred[column])
            return np.log10(values.where(values > 0)), column
    raise ValueError(
        "No recognized CatPred prediction column found. "
        f"Looked for log columns {LOG_COLUMNS} and linear columns {LINEAR_COLUMNS}."
    )


def combine_predictions(pred: pd.DataFrame, meta: pd.DataFrame, pred_log10: pd.Series, pred_col: str) -> pd.DataFrame:
    meta = meta.copy()
    pred = pred.copy()
    pred["prediction_log10"] = pred_log10
    pred["prediction_column"] = pred_col

    if "catpred_row_id" in pred.columns and "catpred_row_id" in meta.columns:
        pred["catpred_row_id"] = pd.to_numeric(pred["catpred_row_id"], errors="coerce").astype("Int64")
        meta["catpred_row_id"] = pd.to_numeric(meta["catpred_row_id"], errors="coerce").astype("Int64")
        combined = meta.merge(
            pred[["catpred_row_id", "prediction_log10", "prediction_column"]],
            on="catpred_row_id",
            how="inner",
            validate="one_to_one",
        )
    elif "entry_id" in pred.columns and "entry_id" in meta.columns:
        combined = meta.merge(
            pred[["entry_id", "prediction_log10", "prediction_column"]],
            on="entry_id",
            how="inner",
            validate="one_to_one",
        )
    else:
        if len(pred) != len(meta):
            raise ValueError(
                f"Prediction rows ({len(pred)}) do not match metadata rows ({len(meta)}), "
                "and no entry_id column is available for merging."
            )
        combined = pd.concat(
            [
                meta.reset_index(drop=True),
                pred[["prediction_log10", "prediction_column"]].reset_index(drop=True),
            ],
            axis=1,
        )

    combined["true_kcat_log10"] = numeric(combined["true_kcat_log10"])
    combined["true_kcat"] = numeric(combined["true_kcat"])
    combined["prediction_kcat"] = np.power(10.0, combined["prediction_log10"])
    combined["error_log10"] = combined["prediction_log10"] - combined["true_kcat_log10"]
    combined["abs_error_log10"] = combined["error_log10"].abs()
    return combined


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


def main() -> None:
    args = parse_args()
    pred = pd.read_csv(args.predictions)
    meta = pd.read_csv(args.metadata)
    pred_log10, pred_col = detect_prediction(pred, args.prediction_column)
    rows = combine_predictions(pred, meta, pred_log10, pred_col)
    metrics = build_metrics(rows)

    args.out_rows.parent.mkdir(parents=True, exist_ok=True)
    args.out_metrics.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.out_rows, index=False)
    metrics.to_csv(args.out_metrics, index=False)

    overall = metrics[(metrics["group_type"] == "all") & (metrics["group"] == "all")].iloc[0]
    print(f"Prediction column: {pred_col}")
    print(f"Evaluated rows: {int(overall['n'])}")
    print(f"MAE log10: {overall.get('mae_log10', np.nan):.4g}")
    print(f"RMSE log10: {overall.get('rmse_log10', np.nan):.4g}")
    print(f"Metrics: {args.out_metrics}")


if __name__ == "__main__":
    main()
