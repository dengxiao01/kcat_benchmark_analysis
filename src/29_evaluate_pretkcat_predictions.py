#!/usr/bin/env python3
"""Evaluate PreTKcat predictions against the finalized benchmark truth."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_PREDICTIONS = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_input_output.csv"
DEFAULT_METADATA = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_input_metadata.csv"
DEFAULT_OUT_ROWS = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_predictions_evaluated.csv"
DEFAULT_OUT_METRICS = BASE / "reports" / "tables" / "pretkcat_eval_metrics.csv"
DEFAULT_OUT_MISSING = BASE / "data" / "final" / "pretkcat" / "pretkcat_invalid_or_unpredicted_rows.csv"
DEFAULT_OUT_MISSING_SUMMARY = BASE / "reports" / "tables" / "pretkcat_invalid_or_unpredicted_summary.csv"

LOG_COLUMNS = ["prediction_log10", "pred_log10_kcat", "predicted_log10_kcat", "log10_kcat"]
LINEAR_COLUMNS = ["prediction_kcat", "pred_kcat", "predicted_kcat", "kcat", "prediction"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score PreTKcat output against benchmark truth.")
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--out-rows", type=Path, default=DEFAULT_OUT_ROWS)
    parser.add_argument("--out-metrics", type=Path, default=DEFAULT_OUT_METRICS)
    parser.add_argument("--out-missing", type=Path, default=DEFAULT_OUT_MISSING)
    parser.add_argument("--out-missing-summary", type=Path, default=DEFAULT_OUT_MISSING_SUMMARY)
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
    raise ValueError(f"No recognized PreTKcat prediction column found. Columns: {list(pred.columns)}")


def combine_predictions(pred: pd.DataFrame, meta: pd.DataFrame, pred_log10: pd.Series, pred_col: str) -> pd.DataFrame:
    pred = pred.copy()
    meta = meta.copy()
    pred["prediction_log10"] = pred_log10
    pred["prediction_column"] = pred_col

    optional = []
    rename_optional = {
        "sequence": "pretkcat_output_sequence",
        "smiles": "pretkcat_output_smiles",
        "prediction_kcat": "pretkcat_output_predicted_kcat",
        "pretkcat_training_mode": "pretkcat_training_mode",
        "pretkcat_wrapper_version": "pretkcat_wrapper_version",
        "pretkcat_training_overlap_policy": "pretkcat_training_overlap_policy",
        "pretkcat_pair_identity_definition": "pretkcat_pair_identity_definition",
        "pretkcat_train_rows": "pretkcat_train_rows",
        "pretkcat_train_rows_raw_usable": "pretkcat_train_rows_raw_usable",
        "pretkcat_train_rows_source_exact_pair_overlap": "pretkcat_train_rows_source_exact_pair_overlap",
        "pretkcat_train_rows_removed_exact_pair": "pretkcat_train_rows_removed_exact_pair",
        "pretkcat_train_rows_removed_near_only": "pretkcat_train_rows_removed_near_only",
        "pretkcat_train_rows_removed_total": "pretkcat_train_rows_removed_total",
        "pretkcat_train_rows_after_pair_exclusion": "pretkcat_train_rows_after_pair_exclusion",
        "pretkcat_near_sequence_identity_threshold_percent": "pretkcat_near_sequence_identity_threshold_percent",
        "pretkcat_near_chemical_tanimoto_threshold": "pretkcat_near_chemical_tanimoto_threshold",
        "pretkcat_near_alignment_coverage_threshold_percent": "pretkcat_near_alignment_coverage_threshold_percent",
        "pretkcat_benchmark_rows_with_joint_neighbor": "pretkcat_benchmark_rows_with_joint_neighbor",
        "pretkcat_feature_dim": "pretkcat_feature_dim",
        "pretkcat_n_estimators": "pretkcat_n_estimators",
    }
    for column, renamed in rename_optional.items():
        if column in pred.columns:
            pred[renamed] = pred[column]
            optional.append(renamed)

    if "pretkcat_row_id" in pred.columns and "pretkcat_row_id" in meta.columns:
        pred["pretkcat_row_id"] = pd.to_numeric(pred["pretkcat_row_id"], errors="coerce").astype("Int64")
        meta["pretkcat_row_id"] = pd.to_numeric(meta["pretkcat_row_id"], errors="coerce").astype("Int64")
        keep = ["pretkcat_row_id", "prediction_log10", "prediction_column"] + optional
        combined = meta.merge(pred[keep], on="pretkcat_row_id", how="inner", validate="one_to_one")
    elif "entry_id" in pred.columns and "entry_id" in meta.columns:
        keep = ["entry_id", "prediction_log10", "prediction_column"] + optional
        combined = meta.merge(pred[keep], on="entry_id", how="inner", validate="one_to_one")
    else:
        if len(pred) != len(meta):
            raise ValueError(
                f"Prediction rows ({len(pred)}) do not match metadata rows ({len(meta)}), "
                "and no pretkcat_row_id or entry_id column is available for merging."
            )
        combined = pd.concat(
            [
                meta.reset_index(drop=True),
                pred[["prediction_log10", "prediction_column"] + optional].reset_index(drop=True),
            ],
            axis=1,
        )

    combined["true_kcat_log10"] = numeric(combined["true_kcat_log10"])
    combined["true_kcat"] = numeric(combined["true_kcat"])
    combined["prediction_kcat"] = np.power(10.0, combined["prediction_log10"])
    combined["error_log10"] = combined["prediction_log10"] - combined["true_kcat_log10"]
    combined["abs_error_log10"] = combined["error_log10"].abs()
    if "pretkcat_source_train_exact_pair_overlap" in combined.columns:
        if "pretkcat_training_overlap_policy" in combined.columns:
            strict = combined["pretkcat_training_overlap_policy"].isin(
                {
                    "strict_sequence_canonical_smiles_pair_disjoint",
                    "sequence_standardized_parent_identity_pair_disjoint",
                    "standardized_sequence_parent_identity_pair_disjoint",
                    "joint_sequence_chemical_near_neighbor_disjoint",
                }
            )
        else:
            strict = pd.Series(False, index=combined.index)
        combined["pretkcat_fit_pair_was_excluded"] = (
            combined["pretkcat_source_train_exact_pair_overlap"].astype(bool) & strict
        )
    return combined


def missing_rows(meta: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    if "pretkcat_row_id" in meta.columns and "pretkcat_row_id" in rows.columns:
        return meta.loc[~meta["pretkcat_row_id"].isin(rows["pretkcat_row_id"])].copy()
    if "entry_id" in meta.columns and "entry_id" in rows.columns:
        return meta.loc[~meta["entry_id"].isin(rows["entry_id"])].copy()
    return meta.iloc[0:0].copy()


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
    overlap_column = "pretkcat_source_train_exact_pair_overlap"
    if overlap_column not in rows.columns and "pretkcat_train_exact_pair_overlap" in rows.columns:
        overlap_column = "pretkcat_train_exact_pair_overlap"
    if overlap_column in rows.columns:
        for overlap, part in rows.groupby(overlap_column, sort=True):
            metrics.append(score_group("source_corpus_exact_pair_overlap", str(bool(overlap)), part))
    return pd.DataFrame(metrics)


def missing_reason(row: pd.Series) -> str:
    if not bool(row.get("pretkcat_smiles_valid", True)):
        return "invalid_smiles"
    return "unpredicted"


def write_missing_summary(missing: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if missing.empty:
        pd.DataFrame(columns=["species", "reason", "rows"]).to_csv(out_path, index=False)
        return
    missing = missing.copy()
    missing["reason"] = missing.apply(missing_reason, axis=1)
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
    pred_log10, pred_col = detect_prediction(pred, args.prediction_column)
    rows = combine_predictions(pred, meta, pred_log10, pred_col)
    metrics = build_metrics(rows)
    missing = missing_rows(meta, rows)

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
