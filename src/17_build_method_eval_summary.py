#!/usr/bin/env python3
"""Build a compact cross-method kcat benchmark summary table."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
METHOD_METRICS = {
    "DLKcat-official": BASE / "reports" / "tables" / "dlkcat_official_eval_metrics.csv",
    "UniKP-official": BASE / "reports" / "tables" / "unikp_official_eval_metrics.csv",
    "TurNuP-official": BASE / "reports" / "tables" / "turnup_eval_metrics.csv",
    "CatPred": BASE / "reports" / "tables" / "catpred_eval_metrics.csv",
    "CataPro": BASE / "reports" / "tables" / "catapro_eval_metrics.csv",
    "PMAK": BASE / "reports" / "tables" / "pmak_eval_metrics.csv",
    "KinForm": BASE / "reports" / "tables" / "kinform_eval_metrics.csv",
    "KcatNet": BASE / "reports" / "tables" / "kcatnet_eval_metrics.csv",
    "PreTKcat": BASE / "reports" / "tables" / "pretkcat_eval_metrics.csv",
    "DEKP-public-retrained": BASE / "reports" / "tables" / "dekp_public_retrained_eval_metrics.csv",
    "SELFprot": BASE / "reports" / "tables" / "selfprot_eval_metrics.csv",
    "GO-HKP": BASE / "reports" / "tables" / "go_hkp_eval_metrics.csv",
}
DEFAULT_OUTPUT = BASE / "reports" / "tables" / "method_eval_summary.csv"


def load_overall(method: str, path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    df = pd.read_csv(path)
    overall = df[(df["group_type"] == "all") & (df["group"] == "all")]
    if overall.empty:
        return None
    row = overall.iloc[0].to_dict()
    row["method"] = method
    return row


def main() -> None:
    rows = [load_overall(method, path) for method, path in METHOD_METRICS.items()]
    rows = [row for row in rows if row is not None]
    if not rows:
        raise FileNotFoundError("No method metrics were found.")
    out = pd.DataFrame(rows)
    front = ["method", "n", "mae_log10", "rmse_log10", "pearson_log10", "spearman_log10"]
    columns = [column for column in front if column in out.columns] + [
        column for column in out.columns if column not in front
    ]
    DEFAULT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out[columns].to_csv(DEFAULT_OUTPUT, index=False)
    print(f"Wrote method summary: {DEFAULT_OUTPUT}")


if __name__ == "__main__":
    main()
