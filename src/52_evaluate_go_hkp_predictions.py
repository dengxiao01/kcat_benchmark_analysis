#!/usr/bin/env python3
"""Evaluate GO-HKP assignments against the finalized benchmark truth."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_PREDICTIONS = BASE / "data" / "final" / "go_hkp" / "go_hkp_kcat_input_output.csv"
DEFAULT_METADATA = BASE / "data" / "final" / "go_hkp" / "go_hkp_kcat_input_metadata.csv"
DEFAULT_ALL_METADATA = BASE / "data" / "final" / "go_hkp" / "go_hkp_kcat_all_metadata.csv"
DEFAULT_OUT_ROWS = BASE / "data" / "final" / "go_hkp" / "go_hkp_kcat_predictions_evaluated.csv"
DEFAULT_OUT_METRICS = BASE / "reports" / "tables" / "go_hkp_eval_metrics.csv"
DEFAULT_OUT_MISSING = BASE / "reports" / "tables" / "go_hkp_missing_summary.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--all-metadata", type=Path, default=DEFAULT_ALL_METADATA)
    parser.add_argument("--out-rows", type=Path, default=DEFAULT_OUT_ROWS)
    parser.add_argument("--out-metrics", type=Path, default=DEFAULT_OUT_METRICS)
    parser.add_argument("--out-missing", type=Path, default=DEFAULT_OUT_MISSING)
    return parser.parse_args()


def combine_predictions(predictions: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    required_prediction = {"entry_id", "prediction_kcat", "prediction_log10"}
    required_metadata = {"entry_id", "true_kcat", "true_kcat_log10"}
    missing_prediction = required_prediction.difference(predictions.columns)
    missing_metadata = required_metadata.difference(metadata.columns)
    if missing_prediction:
        raise ValueError(f"Prediction file is missing columns: {sorted(missing_prediction)}")
    if missing_metadata:
        raise ValueError(f"Metadata file is missing columns: {sorted(missing_metadata)}")
    if predictions["entry_id"].duplicated().any() or metadata["entry_id"].duplicated().any():
        raise ValueError("GO-HKP entry_id values must be unique")

    prediction_columns = [
        column
        for column in [
            "entry_id",
            "prediction_kcat",
            "prediction_log10",
            "go_hkp_assignment_source",
        ]
        if column in predictions.columns
    ]
    rows = metadata.merge(
        predictions[prediction_columns],
        on="entry_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_output"),
    )
    for column in ["true_kcat", "true_kcat_log10", "prediction_kcat", "prediction_log10"]:
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
    finite = np.isfinite(rows[["true_kcat_log10", "prediction_log10"]]).all(axis=1)
    positive = rows[["true_kcat", "prediction_kcat"]].gt(0).all(axis=1)
    rows = rows.loc[finite & positive].copy()
    rows["prediction_column"] = "prediction_log10"
    rows["error_log10"] = rows["prediction_log10"] - rows["true_kcat_log10"]
    rows["abs_error_log10"] = rows["error_log10"].abs()
    return rows


def score_group(group_type: str, group_name: str, rows: pd.DataFrame) -> dict[str, object]:
    truth = rows["true_kcat_log10"].to_numpy(float)
    prediction = rows["prediction_log10"].to_numpy(float)
    error = prediction - truth
    n = len(rows)
    result: dict[str, object] = {"group_type": group_type, "group": group_name, "n": n}
    if n == 0:
        return result
    ss_total = float(np.square(truth - truth.mean()).sum())
    result.update(
        {
            "mae_log10": float(np.abs(error).mean()),
            "rmse_log10": float(np.sqrt(np.square(error).mean())),
            "bias_log10": float(error.mean()),
            "median_abs_error_log10": float(np.median(np.abs(error))),
            "r2_log10": float(1.0 - np.square(error).sum() / ss_total) if ss_total else np.nan,
            "pearson_log10": float(pd.Series(truth).corr(pd.Series(prediction), method="pearson")),
            "spearman_log10": float(pd.Series(truth).corr(pd.Series(prediction), method="spearman")),
            "within_0.3_log10_fraction": float((np.abs(error) <= 0.3).mean()),
            "within_1.0_log10_fraction": float((np.abs(error) <= 1.0).mean()),
        }
    )
    return result


def build_metrics(rows: pd.DataFrame) -> pd.DataFrame:
    metrics = [score_group("all", "all", rows)]
    for species, part in rows.groupby("species", sort=True):
        metrics.append(score_group("species", str(species), part))
    for source, part in rows.groupby("source_database", sort=True):
        metrics.append(score_group("source_database", str(source), part))
    return pd.DataFrame(metrics)


def build_missing_summary(all_metadata: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    predicted_ids = set(rows["entry_id"].astype(str))
    missing = all_metadata.loc[
        ~all_metadata["entry_id"].astype(str).isin(predicted_ids)
    ].copy()
    if missing.empty:
        return pd.DataFrame(columns=["species", "reason", "rows"])
    missing["reason"] = missing["go_hkp_missing_reason"].fillna("unpredicted")
    return missing.groupby(["species", "reason"], dropna=False).size().reset_index(name="rows")


def main() -> None:
    args = parse_args()
    predictions = pd.read_csv(args.predictions)
    metadata = pd.read_csv(args.metadata)
    all_metadata = pd.read_csv(args.all_metadata)
    rows = combine_predictions(predictions, metadata)
    metrics = build_metrics(rows)
    missing = build_missing_summary(all_metadata, rows)

    for path in [args.out_rows, args.out_metrics, args.out_missing]:
        path.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.out_rows, index=False)
    metrics.to_csv(args.out_metrics, index=False)
    missing.to_csv(args.out_missing, index=False)

    overall = metrics.loc[metrics["group_type"].eq("all")].iloc[0]
    print(f"Evaluated rows: {int(overall['n'])}")
    print(f"MAE log10: {float(overall['mae_log10']):.6f}")
    print(f"Unpredicted rows: {int(missing['rows'].sum()) if len(missing) else 0}")
    print(f"Rows: {args.out_rows}")
    print(f"Metrics: {args.out_metrics}")


if __name__ == "__main__":
    main()
