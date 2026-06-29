#!/usr/bin/env python3
"""Check which benchmark rows can be predicted by local KinForm assets."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_KINFORM_ROOT = BASE / "external_methods" / "KinForm"
DEFAULT_METADATA = BASE / "data" / "final" / "kinform" / "kinform_kcat_input_metadata.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "kinform"
DEFAULT_REPORT = BASE / "reports" / "tables" / "kinform_feature_coverage_summary.csv"

EMBEDDING_LAYERS_FOR_NO_COMPUTE = [
    "prot_t5_last",
    "prot_t5_layer_19",
    "esmc_layer_24",
    "esmc_layer_32",
    "esm2_layer_26",
    "esm2_layer_29",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the KinForm-predictable benchmark subset.")
    parser.add_argument("--kinform-root", type=Path, default=DEFAULT_KINFORM_ROOT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def load_sequence_lookup(path: Path) -> dict[str, str]:
    with path.open("rb") as handle:
        lookup = pickle.load(handle)
    if not isinstance(lookup, dict):
        raise TypeError(f"{path} did not contain a dictionary.")
    return {str(seq_id): str(seq) for seq_id, seq in lookup.items()}


def missing_embedding_files(root: Path, seq_id: str) -> list[str]:
    missing: list[str] = []
    for layer in EMBEDDING_LAYERS_FOR_NO_COMPUTE:
        for vec_type in ["mean_vecs", "weighted_vecs"]:
            rel = Path("results") / "protein_embeddings" / layer / vec_type / f"{seq_id}.npy"
            if not (root / rel).exists():
                missing.append(str(rel))
    return missing


def input_rows(df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        rows.append(
            {
                "kinform_row_id": int(row["kinform_row_id"]),
                "entry_id": row["entry_id"],
                "sequence": row["sequence"],
                "smiles": row["kinform_canonical_smiles"] or row["SMILES"],
                "value": float(row["true_kcat"]),
            }
        )
    return rows


def write_json(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, ensure_ascii=False, indent=2)


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "kinform_row_id",
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


def input_csv(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinform_row_id": df["kinform_row_id"],
            "entry_id": df["entry_id"],
            "sequence": df["sequence"],
            "smiles": df["kinform_canonical_smiles"].where(
                df["kinform_canonical_smiles"].astype(str).str.len() > 0, df["SMILES"]
            ),
            "value": df["true_kcat"],
            "true_kcat_log10": df["true_kcat_log10"],
            "kinform_sequence_id": df["kinform_sequence_id"],
        }
    )


def reason(row: pd.Series) -> str:
    reasons = []
    if not bool(row.get("kinform_smiles_valid", False)):
        reasons.append("invalid_smiles")
    if not bool(row.get("kinform_sequence_in_lookup", False)):
        reasons.append("missing_sequence_lookup")
    if not bool(row.get("kinform_precomputed_features_ready", False)):
        reasons.append("missing_precomputed_features")
    if not bool(row.get("kinform_assets_available", False)):
        reasons.append("missing_kinform_assets")
    return ";".join(reasons) if reasons else "ready"


def write_report(df: pd.DataFrame, report: Path) -> None:
    rows = []
    for group_name, part in [("all", df)] + list(df.groupby("species", sort=True)):
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "valid_smiles_rows": int(part["kinform_smiles_valid"].sum()),
                "sequence_lookup_rows": int(part["kinform_sequence_in_lookup"].sum()),
                "precomputed_feature_rows": int(part["kinform_precomputed_features_ready"].sum()),
                "predictable_rows": int(part["kinform_predictable"].sum()),
                "missing_sequence_lookup_rows": int((~part["kinform_sequence_in_lookup"]).sum()),
                "missing_precomputed_feature_rows": int(
                    part["kinform_sequence_in_lookup"].sum()
                    - part["kinform_precomputed_features_ready"].sum()
                ),
                "unique_predictable_sequences": part.loc[
                    part["kinform_predictable"], "kinform_model_sequence"
                ].nunique(),
            }
        )
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report, index=False)


def main() -> None:
    args = parse_args()
    meta = pd.read_csv(args.metadata)
    seq_lookup_path = args.kinform_root / "results" / "sequence_id_to_sequence.pkl"
    trained_model_path = args.kinform_root / "results" / "trained_models" / "kcat_KinForm-L" / "model.joblib"
    transformer_path = args.kinform_root / "results" / "trained_models" / "kcat_KinForm-L" / "transformers.joblib"
    binding_site_path = args.kinform_root / "results" / "binding_sites" / "binding_sites_all.tsv"

    assets_available = all(
        path.exists()
        for path in [seq_lookup_path, trained_model_path, transformer_path, binding_site_path]
    )
    if seq_lookup_path.exists():
        seq_id_to_seq = load_sequence_lookup(seq_lookup_path)
    else:
        seq_id_to_seq = {}
    seq_to_id = {seq: seq_id for seq_id, seq in seq_id_to_seq.items()}

    rows = []
    for _, row in meta.iterrows():
        model_seq = str(row["kinform_model_sequence"])
        seq_id = seq_to_id.get(model_seq, "")
        missing_files = missing_embedding_files(args.kinform_root, seq_id) if seq_id else []
        rows.append(
            {
                "kinform_sequence_id": seq_id,
                "kinform_assets_available": assets_available,
                "kinform_sequence_in_lookup": bool(seq_id),
                "kinform_precomputed_features_ready": bool(seq_id) and not missing_files,
                "kinform_missing_precomputed_files": ";".join(missing_files),
            }
        )

    out = pd.concat([meta.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    out["kinform_predictable"] = (
        out["kinform_assets_available"]
        & out["kinform_smiles_valid"].astype(bool)
        & out["kinform_sequence_in_lookup"]
        & out["kinform_precomputed_features_ready"]
    )
    out["kinform_unavailable_reason"] = out.apply(reason, axis=1)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    predictable = out[out["kinform_predictable"]].copy()
    unavailable = out[~out["kinform_predictable"]].copy()

    out.to_csv(args.out_dir / "kinform_feature_coverage.csv", index=False)
    predictable.to_csv(args.out_dir / "kinform_kcat_input_predictable_metadata.csv", index=False)
    input_csv(predictable).to_csv(args.out_dir / "kinform_kcat_input_predictable.csv", index=False)
    truth(predictable).to_csv(args.out_dir / "kinform_kcat_input_predictable_truth.csv", index=False)
    write_json(input_rows(predictable), args.out_dir / "kinform_kcat_input_predictable.json")
    unavailable.to_csv(args.out_dir / "kinform_unavailable_rows.csv", index=False)
    write_report(out, args.report)

    print(f"KinForm assets available: {assets_available}")
    print(f"Rows: {len(out)} total; {len(predictable)} predictable; {len(unavailable)} unavailable")
    print(f"Coverage report: {args.report}")


if __name__ == "__main__":
    main()
