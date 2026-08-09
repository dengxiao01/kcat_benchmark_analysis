#!/usr/bin/env python3
"""Prepare PreTKcat-ready kcat benchmark inputs.

PreTKcat uses enzyme sequence, substrate SMILES, and temperature features.
The public repository provides training data and feature/model code, but not a
standalone fitted kcat regressor, so this preparation step keeps enough
metadata to make the later retraining-based evaluation traceable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize


RDLogger.DisableLog("rdApp.*")
UNCHARGER = rdMolStandardize.Uncharger()

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_TRAIN = BASE / "external_methods" / "PreTKcat" / "datasets" / "DLTKcat_data" / "kcat_merge_DLTKcat.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "pretkcat"
DEFAULT_REPORT = BASE / "reports" / "tables" / "pretkcat_eval_readiness.csv"

DEFAULT_TEMP_C = 30.0
TEMP_C_MIN = 0.0
TEMP_C_MAX = 100.0
INV_TEMP_MIN = 1.0 / (TEMP_C_MAX + 273.15)
INV_TEMP_MAX = 1.0 / (TEMP_C_MIN + 273.15)

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
        description="Create PreTKcat input, metadata, and truth files for the kcat benchmark."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--train-data", type=Path, default=DEFAULT_TRAIN)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--sample-size", type=int, default=20)
    parser.add_argument("--default-temp-c", type=float, default=DEFAULT_TEMP_C)
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    empty = {
        column: int(clean_text(df[column]).eq("").sum())
        for column in ["entry_id", "SMILES", "sequence"]
    }
    bad = {column: count for column, count in empty.items() if count}
    if bad:
        raise ValueError(f"{path} contains empty required values: {bad}")
    if df["entry_id"].duplicated().any():
        examples = ", ".join(df.loc[df["entry_id"].duplicated(), "entry_id"].head(5))
        raise ValueError(f"entry_id must be unique. Examples: {examples}")


def canonical_smiles(smiles: object) -> tuple[bool, str]:
    text = str(smiles).strip()
    if not text or text == "nan":
        return False, ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return False, ""
    return True, Chem.MolToSmiles(mol, canonical=True)


def chemical_parent_key(smiles: object) -> str:
    text = str(smiles).strip()
    if not text or text == "nan":
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    try:
        parent = rdMolStandardize.FragmentParent(mol)
        parent = UNCHARGER.uncharge(parent)
        Chem.SanitizeMol(parent)
    except Exception:
        parent = Chem.Mol(mol)
        Chem.SanitizeMol(parent)
    parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    try:
        key = Chem.MolToInchiKey(parent).split("-", 1)[0]
    except Exception:
        key = ""
    return key or parent_smiles


def pretkcat_model_sequence(sequence: object) -> str:
    text = str(sequence).strip()
    if len(text) > 1000:
        return text[:500] + text[-500:]
    return text


def add_temperature_features(df: pd.DataFrame, default_temp_c: float) -> pd.DataFrame:
    out = df.copy()
    if "temperature_c" in out.columns:
        temperature = pd.to_numeric(out["temperature_c"], errors="coerce")
    else:
        temperature = pd.Series([pd.NA] * len(out), index=out.index, dtype="Float64")
    out["pretkcat_temperature_imputed"] = temperature.isna()
    out["pretkcat_temperature_c"] = temperature.fillna(default_temp_c).astype(float)
    out["pretkcat_temperature_outside_training_range"] = (
        (out["pretkcat_temperature_c"] < TEMP_C_MIN) | (out["pretkcat_temperature_c"] > TEMP_C_MAX)
    )
    out["pretkcat_temp_k"] = out["pretkcat_temperature_c"] + 273.15
    out["pretkcat_inv_temp"] = 1.0 / out["pretkcat_temp_k"]
    out["pretkcat_temp_k_norm"] = (out["pretkcat_temperature_c"] - TEMP_C_MIN) / (TEMP_C_MAX - TEMP_C_MIN)
    out["pretkcat_inv_temp_norm"] = (out["pretkcat_inv_temp"] - INV_TEMP_MIN) / (INV_TEMP_MAX - INV_TEMP_MIN)
    return out


def training_overlap_sets(train_path: Path) -> dict[str, set[str]]:
    empty = {"sequence": set(), "smiles": set(), "pair": set()}
    if not train_path.exists():
        return empty
    train = pd.read_csv(train_path)
    required = {"seq", "smiles", "kcat"}
    if not required.issubset(train.columns):
        return empty
    train = train.copy()
    train["kcat"] = pd.to_numeric(train["kcat"], errors="coerce")
    train = train[(train["kcat"] > 0) & train["seq"].notna() & train["smiles"].notna()].copy()
    train = train[~train["smiles"].astype(str).str.contains(".", regex=False)].copy()
    train["pretkcat_model_sequence"] = train["seq"].map(pretkcat_model_sequence)
    smiles_status = train["smiles"].map(canonical_smiles)
    train["smiles_valid"] = [valid for valid, _ in smiles_status]
    train["canonical_smiles"] = [canonical for _, canonical in smiles_status]
    train["chemical_parent_key"] = train["smiles"].map(chemical_parent_key)
    train = train[train["smiles_valid"] & train["chemical_parent_key"].ne("")].copy()
    sequence_set = set(train["pretkcat_model_sequence"].astype(str))
    smiles_set = set(train["chemical_parent_key"].astype(str))
    pair_set = set((train["pretkcat_model_sequence"] + "||" + train["chemical_parent_key"]).astype(str))
    return {"sequence": sequence_set, "smiles": smiles_set, "pair": pair_set}


def enrich(df: pd.DataFrame, train_path: Path, default_temp_c: float) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["pretkcat_row_id"] = range(len(out))
    out["entry_id"] = clean_text(out["entry_id"])
    out["sequence"] = clean_text(out["sequence"])
    out["SMILES"] = clean_text(out["SMILES"])
    out["pretkcat_model_sequence"] = out["sequence"].map(pretkcat_model_sequence)
    out["pretkcat_sequence_truncated"] = out["sequence"] != out["pretkcat_model_sequence"]

    smiles_status = out["SMILES"].map(canonical_smiles)
    out["pretkcat_smiles_valid"] = [valid for valid, _ in smiles_status]
    out["pretkcat_canonical_smiles"] = [canonical for _, canonical in smiles_status]
    out["pretkcat_standardized_parent_key"] = out["SMILES"].map(chemical_parent_key)
    out = add_temperature_features(out, default_temp_c)

    overlap = training_overlap_sets(train_path)
    out["pretkcat_train_exact_sequence_overlap"] = out["pretkcat_model_sequence"].isin(overlap["sequence"])
    out["pretkcat_train_exact_smiles_overlap"] = out["pretkcat_standardized_parent_key"].isin(overlap["smiles"])
    keys = out["pretkcat_model_sequence"].astype(str) + "||" + out["pretkcat_standardized_parent_key"].astype(str)
    out["pretkcat_train_exact_pair_overlap"] = keys.isin(overlap["pair"])
    # Explicit aliases clarify that overlap is measured against the original
    # public source corpus. The strict wrapper removes exact pairs before fit.
    out["pretkcat_source_train_exact_sequence_overlap"] = out["pretkcat_train_exact_sequence_overlap"]
    out["pretkcat_source_train_exact_smiles_overlap"] = out["pretkcat_train_exact_smiles_overlap"]
    out["pretkcat_source_train_exact_pair_overlap"] = out["pretkcat_train_exact_pair_overlap"]
    return out


def pretkcat_input(df: pd.DataFrame) -> pd.DataFrame:
    smiles = df["pretkcat_canonical_smiles"].where(
        df["pretkcat_canonical_smiles"].astype(str).str.len() > 0, df["SMILES"]
    )
    return pd.DataFrame(
        {
            "pretkcat_row_id": df["pretkcat_row_id"],
            "entry_id": df["entry_id"],
            "sequence": clean_text(df["sequence"]),
            "smiles": smiles,
            "temp_k_norm": df["pretkcat_temp_k_norm"],
            "inv_temp_norm": df["pretkcat_inv_temp_norm"],
            "temperature_c": df["pretkcat_temperature_c"],
            "true_kcat": df["true_kcat"],
            "true_kcat_log10": df["true_kcat_log10"],
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "pretkcat_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "pretkcat_canonical_smiles",
        "pretkcat_standardized_parent_key",
        "pretkcat_smiles_valid",
        "sequence",
        "pretkcat_model_sequence",
        "pretkcat_sequence_truncated",
        "true_kcat",
        "true_kcat_log10",
        "unit",
        "pH",
        "temperature_c",
        "pretkcat_temperature_c",
        "pretkcat_temperature_imputed",
        "pretkcat_temperature_outside_training_range",
        "pretkcat_temp_k",
        "pretkcat_inv_temp",
        "pretkcat_temp_k_norm",
        "pretkcat_inv_temp_norm",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
        "enzyme_complex_type",
        "pretkcat_train_exact_sequence_overlap",
        "pretkcat_train_exact_smiles_overlap",
        "pretkcat_train_exact_pair_overlap",
        "pretkcat_source_train_exact_sequence_overlap",
        "pretkcat_source_train_exact_smiles_overlap",
        "pretkcat_source_train_exact_pair_overlap",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "pretkcat_row_id",
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


def balanced_sample(df: pd.DataFrame, size: int) -> pd.DataFrame:
    if size <= 0 or len(df) <= size:
        return df.copy()
    groups = [(name, part.copy()) for name, part in df.groupby("species", sort=True)]
    per_group = max(1, size // max(len(groups), 1))
    selected_indexes: list[int] = []
    for _, part in groups:
        selected_indexes.extend(part.head(per_group).index.tolist())
    selected_set = set(selected_indexes)
    remaining = [idx for idx in df.index.tolist() if idx not in selected_set]
    selected_indexes.extend(remaining[: max(0, size - len(selected_indexes))])
    return df.loc[selected_indexes[:size]].sort_index().copy()


def write_bundle(df: pd.DataFrame, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    pretkcat_input(df).to_csv(prefix.with_suffix(".csv"), index=False)
    metadata(df).to_csv(prefix.parent / f"{prefix.name}_metadata.csv", index=False)
    truth(df).to_csv(prefix.parent / f"{prefix.name}_truth.csv", index=False)


def write_report(df: pd.DataFrame, report: Path) -> None:
    rows = []
    for group_name, part in [("all", df)] + list(df.groupby("species", sort=True)):
        source_counts = part["source_database"].value_counts(dropna=False).to_dict()
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "valid_smiles_rows": int(part["pretkcat_smiles_valid"].sum()),
                "invalid_smiles_rows": int((~part["pretkcat_smiles_valid"]).sum()),
                "temperature_imputed_rows": int(part["pretkcat_temperature_imputed"].sum()),
                "temperature_outside_training_range_rows": int(part["pretkcat_temperature_outside_training_range"].sum()),
                "truncated_sequence_rows": int(part["pretkcat_sequence_truncated"].sum()),
                "unique_sequences": part["sequence"].nunique(),
                "unique_model_sequences": part["pretkcat_model_sequence"].nunique(),
                "unique_smiles": part["SMILES"].nunique(),
                "source_corpus_exact_pair_overlap_rows": int(part["pretkcat_train_exact_pair_overlap"].sum()),
                "source_corpus_exact_sequence_overlap_rows": int(part["pretkcat_train_exact_sequence_overlap"].sum()),
                "source_corpus_exact_smiles_overlap_rows": int(part["pretkcat_train_exact_smiles_overlap"].sum()),
                "brenda_rows": int(source_counts.get("BRENDA", 0)),
                "sabiork_rows": int(source_counts.get("SABIO-RK", 0)),
                "brenda_sabiork_rows": int(source_counts.get("BRENDA;SABIO-RK", 0)),
            }
        )
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report, index=False)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)
    df = enrich(df, args.train_data, args.default_temp_c)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    valid = df[df["pretkcat_smiles_valid"]].copy()
    invalid = df[~df["pretkcat_smiles_valid"]].copy()

    write_bundle(df, args.out_dir / "pretkcat_kcat_input")
    write_bundle(valid, args.out_dir / "pretkcat_kcat_input_valid_smiles")
    metadata(invalid).to_csv(args.out_dir / "pretkcat_invalid_smiles_rows.csv", index=False)

    if args.sample_size > 0:
        sample = balanced_sample(valid, args.sample_size)
        write_bundle(sample, args.out_dir / f"pretkcat_kcat_input_sample{len(sample)}")

    write_report(df, args.report)
    print(f"Wrote PreTKcat input bundle to {args.out_dir}")
    print(f"Rows: {len(df)} total; {len(valid)} valid SMILES; {len(invalid)} invalid SMILES")
    print(f"Temperature imputed rows: {int(df['pretkcat_temperature_imputed'].sum())}")
    print(f"Exact pair overlaps in original public training corpus: {int(df['pretkcat_train_exact_pair_overlap'].sum())}")
    print(f"Readiness report: {args.report}")


if __name__ == "__main__":
    main()
