#!/usr/bin/env python3
"""Prepare DEKP-ready benchmark inputs and asset coverage reports.

DEKP is structure-aware: besides sequence and substrate SMILES it needs protein
structure-derived graph features. The public repository does not ship a fitted
kcat model checkpoint, so this preparation step focuses on traceable inputs and
on clearly reporting which rows have the structural assets required by DEKP.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import rarfile
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_DEKP_DATA = BASE / "external_methods" / "DEKP" / "datasets" / "kcat_dataset.csv"
DEFAULT_STRUCTURE_ARCHIVE = BASE / "external_methods" / "DEKP" / "datasets" / "protein_structure_datasets.rar"
DEFAULT_LOCAL_STRUCTURES = BASE / "external_methods" / "DEKP" / "structures" / "benchmark" / "AlphaFold"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "dekp"
DEFAULT_REPORT = BASE / "reports" / "tables" / "dekp_eval_readiness.csv"
DEFAULT_ASSET_REPORT = BASE / "reports" / "tables" / "dekp_asset_coverage.csv"

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
    parser = argparse.ArgumentParser(description="Create DEKP input, metadata, truth, and coverage files.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--dekp-data", type=Path, default=DEFAULT_DEKP_DATA)
    parser.add_argument("--structure-archive", type=Path, default=DEFAULT_STRUCTURE_ARCHIVE)
    parser.add_argument("--local-structures", type=Path, default=DEFAULT_LOCAL_STRUCTURES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--asset-report", type=Path, default=DEFAULT_ASSET_REPORT)
    parser.add_argument("--sample-size", type=int, default=20)
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    empty = {
        column: int(clean_text(df[column]).eq("").sum())
        for column in ["entry_id", "SMILES", "sequence", "uniprot_id"]
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


def archive_structure_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    ids: set[str] = set()
    with rarfile.RarFile(path) as handle:
        for info in handle.infolist():
            name = Path(info.filename).name
            if name.endswith(".pdb"):
                ids.add(name[:-4])
    return ids


def local_structure_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {item.stem for item in path.glob("*.pdb")}


def load_dekp_training(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    required = {"UniprotID", "Sequence", "Smiles", "Label"}
    if not required.issubset(df.columns):
        return pd.DataFrame()
    df = df.copy()
    df["Label"] = pd.to_numeric(df["Label"], errors="coerce")
    df = df.dropna(subset=["UniprotID", "Sequence", "Smiles", "Label"]).copy()
    df["dekp_train_smiles_valid"] = False
    df["dekp_train_canonical_smiles"] = ""
    statuses = df["Smiles"].map(canonical_smiles)
    df["dekp_train_smiles_valid"] = [valid for valid, _ in statuses]
    df["dekp_train_canonical_smiles"] = [canonical for _, canonical in statuses]
    df = df[df["dekp_train_smiles_valid"]].copy()
    return df


def enrich(df: pd.DataFrame, dekp_train: pd.DataFrame, archive_ids: set[str], local_ids: set[str]) -> pd.DataFrame:
    out = df.reset_index(drop=True).copy()
    out["dekp_row_id"] = range(len(out))
    out["entry_id"] = clean_text(out["entry_id"])
    out["sequence"] = clean_text(out["sequence"])
    out["SMILES"] = clean_text(out["SMILES"])
    out["uniprot_id"] = clean_text(out["uniprot_id"])
    statuses = out["SMILES"].map(canonical_smiles)
    out["dekp_smiles_valid"] = [valid for valid, _ in statuses]
    out["dekp_canonical_smiles"] = [canonical for _, canonical in statuses]
    out["dekp_sequence_length"] = out["sequence"].str.len()
    out["dekp_smiles_length"] = out["dekp_canonical_smiles"].str.len()

    out["dekp_author_structure_available"] = out["uniprot_id"].isin(archive_ids)
    out["dekp_local_structure_available"] = out["uniprot_id"].isin(local_ids)
    out["dekp_structure_available"] = out["dekp_author_structure_available"] | out["dekp_local_structure_available"]
    out["dekp_structure_source"] = "missing"
    out.loc[out["dekp_author_structure_available"], "dekp_structure_source"] = "dekp_author_structure_archive"
    out.loc[
        (~out["dekp_author_structure_available"]) & out["dekp_local_structure_available"],
        "dekp_structure_source",
    ] = "local_benchmark_structure"

    out["dekp_train_uniprot_overlap"] = False
    out["dekp_train_exact_pair_overlap"] = False
    out["dekp_train_exact_sequence_overlap"] = False
    out["dekp_train_exact_smiles_overlap"] = False
    if not dekp_train.empty:
        train_uniprot = set(dekp_train["UniprotID"].astype(str))
        train_sequence = set(dekp_train["Sequence"].astype(str))
        train_smiles = set(dekp_train["dekp_train_canonical_smiles"].astype(str))
        train_pair = set((dekp_train["Sequence"].astype(str) + "||" + dekp_train["dekp_train_canonical_smiles"].astype(str)))
        out["dekp_train_uniprot_overlap"] = out["uniprot_id"].isin(train_uniprot)
        out["dekp_train_exact_sequence_overlap"] = out["sequence"].isin(train_sequence)
        out["dekp_train_exact_smiles_overlap"] = out["dekp_canonical_smiles"].isin(train_smiles)
        keys = out["sequence"].astype(str) + "||" + out["dekp_canonical_smiles"].astype(str)
        out["dekp_train_exact_pair_overlap"] = keys.isin(train_pair)

    smiles_values = out["dekp_canonical_smiles"].where(
        out["dekp_canonical_smiles"].astype(str).str.len() > 0, out["SMILES"]
    )
    out["dekp_cid"] = pd.factorize(smiles_values, sort=True)[0].astype(int)
    return out


def dekp_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ECNumber": clean_text(df["ec_number"]),
            "Organism": clean_text(df["species"]),
            "Smiles": df["dekp_canonical_smiles"].where(
                df["dekp_canonical_smiles"].astype(str).str.len() > 0, df["SMILES"]
            ),
            "Substrate": clean_text(df["substrate_name"]),
            "Sequence": clean_text(df["sequence"]),
            "Type": "wildtype",
            "Label": pd.to_numeric(df["true_kcat_log10"], errors="coerce"),
            "Unit": clean_text(df["unit"]),
            "UniprotID": clean_text(df["uniprot_id"]),
            "CID": df["dekp_cid"],
            "Set": "benchmark",
            "dekp_row_id": df["dekp_row_id"],
            "entry_id": df["entry_id"],
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "dekp_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "dekp_canonical_smiles",
        "dekp_smiles_valid",
        "sequence",
        "dekp_sequence_length",
        "dekp_smiles_length",
        "dekp_cid",
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
        "dekp_author_structure_available",
        "dekp_local_structure_available",
        "dekp_structure_available",
        "dekp_structure_source",
        "dekp_train_uniprot_overlap",
        "dekp_train_exact_sequence_overlap",
        "dekp_train_exact_smiles_overlap",
        "dekp_train_exact_pair_overlap",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "dekp_row_id",
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
    dekp_input(df).to_csv(prefix.with_suffix(".csv"), index=False)
    metadata(df).to_csv(prefix.parent / f"{prefix.name}_metadata.csv", index=False)
    truth(df).to_csv(prefix.parent / f"{prefix.name}_truth.csv", index=False)


def write_readiness_report(df: pd.DataFrame, report: Path) -> None:
    rows = []
    for group_name, part in [("all", df)] + list(df.groupby("species", sort=True)):
        valid = part[part["dekp_smiles_valid"] & part["dekp_structure_available"]]
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "valid_smiles_rows": int(part["dekp_smiles_valid"].sum()),
                "invalid_smiles_rows": int((~part["dekp_smiles_valid"]).sum()),
                "structure_available_rows": int(part["dekp_structure_available"].sum()),
                "structure_missing_rows": int((~part["dekp_structure_available"]).sum()),
                "predictable_rows": len(valid),
                "unique_uniprot": part["uniprot_id"].nunique(),
                "unique_uniprot_with_structure": part.loc[part["dekp_structure_available"], "uniprot_id"].nunique(),
                "unique_sequences": part["sequence"].nunique(),
                "unique_smiles": part["dekp_canonical_smiles"].nunique(),
                "train_uniprot_overlap_rows": int(part["dekp_train_uniprot_overlap"].sum()),
                "train_exact_pair_overlap_rows": int(part["dekp_train_exact_pair_overlap"].sum()),
            }
        )
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report, index=False)


def write_asset_report(df: pd.DataFrame, archive_ids: set[str], local_ids: set[str], dekp_train: pd.DataFrame, report: Path) -> None:
    benchmark_uniprot = set(df["uniprot_id"].astype(str))
    train_uniprot = set(dekp_train["UniprotID"].astype(str)) if not dekp_train.empty else set()
    rows = [
        {"asset": "dekp_public_kcat_dataset_rows", "count": len(dekp_train)},
        {"asset": "dekp_public_kcat_dataset_unique_uniprot", "count": len(train_uniprot)},
        {"asset": "author_structure_archive_unique_pdb", "count": len(archive_ids)},
        {"asset": "local_benchmark_structure_unique_pdb", "count": len(local_ids)},
        {"asset": "benchmark_unique_uniprot", "count": len(benchmark_uniprot)},
        {"asset": "benchmark_unique_uniprot_with_any_structure", "count": len(benchmark_uniprot & (archive_ids | local_ids))},
        {"asset": "benchmark_unique_uniprot_missing_structure", "count": len(benchmark_uniprot - (archive_ids | local_ids))},
        {"asset": "benchmark_rows_with_any_structure", "count": int(df["dekp_structure_available"].sum())},
        {"asset": "benchmark_rows_missing_structure", "count": int((~df["dekp_structure_available"]).sum())},
        {"asset": "benchmark_rows_with_valid_smiles_and_structure", "count": int((df["dekp_smiles_valid"] & df["dekp_structure_available"]).sum())},
        {"asset": "benchmark_rows_exact_train_pair_overlap", "count": int(df["dekp_train_exact_pair_overlap"].sum())},
    ]
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report, index=False)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)
    dekp_train = load_dekp_training(args.dekp_data)
    archive_ids = archive_structure_ids(args.structure_archive)
    local_ids = local_structure_ids(args.local_structures)
    df = enrich(df, dekp_train, archive_ids, local_ids)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    valid_smiles = df[df["dekp_smiles_valid"]].copy()
    predictable = df[df["dekp_smiles_valid"] & df["dekp_structure_available"]].copy()
    invalid = df[~df["dekp_smiles_valid"]].copy()
    missing_structure = df[~df["dekp_structure_available"]].copy()

    write_bundle(df, args.out_dir / "dekp_kcat_input")
    write_bundle(valid_smiles, args.out_dir / "dekp_kcat_input_valid_smiles")
    write_bundle(predictable, args.out_dir / "dekp_kcat_input_structure_available")
    metadata(invalid).to_csv(args.out_dir / "dekp_invalid_smiles_rows.csv", index=False)
    metadata(missing_structure).to_csv(args.out_dir / "dekp_missing_structure_rows.csv", index=False)

    if args.sample_size > 0:
        sample = balanced_sample(predictable, args.sample_size)
        write_bundle(sample, args.out_dir / f"dekp_kcat_input_sample{len(sample)}")

    write_readiness_report(df, args.report)
    write_asset_report(df, archive_ids, local_ids, dekp_train, args.asset_report)
    print(f"Wrote DEKP input bundle to {args.out_dir}")
    print(f"Rows: {len(df)} total; {len(valid_smiles)} valid SMILES; {len(predictable)} with valid SMILES and structure")
    print(f"Missing structures: {len(missing_structure)} rows")
    print(f"Readiness report: {args.report}")
    print(f"Asset report: {args.asset_report}")


if __name__ == "__main__":
    main()
