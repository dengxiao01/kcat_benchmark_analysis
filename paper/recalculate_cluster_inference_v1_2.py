#!/usr/bin/env python3
"""Independently reconstruct V4 Table S14 from the released workbook."""

from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import wilcoxon


BASE = Path(__file__).resolve().parents[1]
WORKBOOK = BASE / "paper" / "Supplementary_tables.xlsx"
OUTPUT = BASE / "paper" / "independent_cluster_inference_v1.2.0-r3.csv"
METHOD_FILES = {
    "KcatNet": BASE / "data" / "final" / "kcatnet" / "kcatnet_kcat_predictions_evaluated.csv",
    "CataPro": BASE / "data" / "final" / "catapro" / "catapro_kcat_predictions_evaluated.csv",
    "TurNuP": BASE / "data" / "final" / "turnup" / "turnup_kcat_predictions_evaluated.csv",
    "PMAK": BASE / "data" / "final" / "pmak" / "pmak_kcat_predictions_evaluated.csv",
}


def bh_adjust(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def main() -> None:
    record = pd.read_excel(WORKBOOK, sheet_name="Table S22")
    method_errors = []
    for method, path in METHOD_FILES.items():
        frame = pd.read_csv(path, usecols=["entry_id", "abs_error_log10"])
        frame["method"] = method
        method_errors.append(frame)
    errors = pd.concat(method_errors, ignore_index=True).pivot(
        index="entry_id", columns="method", values="abs_error_log10"
    )
    cluster_types = ["protein", "pair", "reaction", "reference", "label_assignment"]
    data = errors.merge(
        record.set_index("entry_id")[[f"{name}_cluster" for name in cluster_types]],
        left_index=True,
        right_index=True,
        validate="one_to_one",
    )
    comparisons = [
        ("broad_common", "KcatNet", "CataPro"),
        ("reaction_common", "KcatNet", "TurNuP"),
        ("reaction_common", "KcatNet", "PMAK"),
        ("reaction_common", "TurNuP", "PMAK"),
    ]
    rows = []
    for scope, method_a, method_b in comparisons:
        common = data.loc[data[method_a].notna() & data[method_b].notna()]
        for cluster_type in cluster_types:
            grouped = common.groupby(f"{cluster_type}_cluster")[[method_a, method_b]].mean()
            test = wilcoxon(grouped[method_a], grouped[method_b], alternative="two-sided")
            rows.append(
                {
                    "comparison_scope": scope,
                    "method_a": method_a,
                    "method_b": method_b,
                    "cluster_type": cluster_type,
                    "recalculated_n_common_rows": len(common),
                    "recalculated_n_paired_clusters": len(grouped),
                    "recalculated_cluster_mean_error_a": grouped[method_a].mean(),
                    "recalculated_cluster_mean_error_b": grouped[method_b].mean(),
                    "recalculated_wilcoxon_statistic": float(test.statistic),
                    "recalculated_p_value_raw": float(test.pvalue),
                }
            )
    calculated = pd.DataFrame(rows)
    calculated["recalculated_p_value_bh_global"] = bh_adjust(
        calculated["recalculated_p_value_raw"].to_numpy()
    )
    expected = pd.read_excel(WORKBOOK, sheet_name="Table S14")
    keys = ["comparison_scope", "method_a", "method_b", "cluster_type"]
    merged = calculated.merge(expected, on=keys, validate="one_to_one")
    pairs = {
        "n_common_rows": "recalculated_n_common_rows",
        "n_paired_clusters": "recalculated_n_paired_clusters",
        "cluster_mean_error_a": "recalculated_cluster_mean_error_a",
        "cluster_mean_error_b": "recalculated_cluster_mean_error_b",
        "wilcoxon_statistic": "recalculated_wilcoxon_statistic",
        "p_value_raw": "recalculated_p_value_raw",
        "p_value_bh_global": "recalculated_p_value_bh_global",
    }
    for expected_column, calculated_column in pairs.items():
        merged[f"match_{expected_column}"] = np.isclose(
            merged[calculated_column].astype(float),
            merged[expected_column].astype(float),
            rtol=1e-10,
            atol=1e-12,
        )
    merged["match_significance"] = (
        merged["recalculated_p_value_bh_global"].lt(0.05)
        == merged["significant_bh_fdr_0.05"].astype(bool)
    )
    merged["scipy_version"] = scipy.__version__
    merged.to_csv(OUTPUT, index=False)
    match_columns = [column for column in merged if column.startswith("match_")]
    failures = int((~merged[match_columns]).sum().sum())
    print(f"Independent Table S14 rows: {len(merged)}; failed comparisons: {failures}")
    print(f"Output: {OUTPUT}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
