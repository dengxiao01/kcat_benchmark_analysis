#!/usr/bin/env python3
"""Rebuild CatPred reference-corpus overlap flags for the canonical benchmark."""

from __future__ import annotations

import argparse
import io
import tarfile
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_BENCHMARK = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_ARCHIVE = BASE / "external_methods" / "CatPred_datas" / "catpred-db.tar.gz"
DEFAULT_OUTPUT = BASE / "reports" / "tables" / "catpred_db_vs_our_benchmark_overlap.csv"
CATPRED_TABLES = [
    "CatPred-DB/data/kcat/kcat-random_trainval.csv",
    "CatPred-DB/data/kcat/kcat-random_val.csv",
    "CatPred-DB/data/kcat/kcat-random_test.csv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def load_reference(archive: Path) -> pd.DataFrame:
    if not archive.exists():
        raise FileNotFoundError(f"CatPred-DB archive not found: {archive}")
    frames = []
    with tarfile.open(archive, "r:gz") as handle:
        for member_name in CATPRED_TABLES:
            member = handle.getmember(member_name)
            extracted = handle.extractfile(member)
            if extracted is None:
                raise FileNotFoundError(f"Could not read {member_name} from {archive}")
            frames.append(pd.read_csv(io.BytesIO(extracted.read())))
    return pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)


def main() -> None:
    args = parse_args()
    benchmark = pd.read_csv(args.benchmark)
    reference = load_reference(args.archive)

    sequence_full = set(zip(reference["sequence"].astype(str), reference["reactant_smiles"].astype(str)))
    uniprot_full = set(zip(reference["uniprot"].astype(str), reference["reactant_smiles"].astype(str)))
    uniprots = set(reference["uniprot"].astype(str))
    ecs = set(reference["ec"].dropna().astype(str))
    uniprot_components = {
        (str(row.uniprot), component)
        for row in reference.itertuples(index=False)
        for component in str(row.reactant_smiles).split(".")
    }

    output = benchmark[
        ["entry_id", "species", "uniprot_id", "ec_number", "substrate_name", "SMILES", "sequence"]
    ].copy()
    output["exact_sequence_reactant_smiles_overlap"] = [
        (str(row.sequence), str(row.SMILES)) in sequence_full for row in output.itertuples(index=False)
    ]
    output["exact_uniprot_reactant_smiles_overlap"] = [
        (str(row.uniprot_id), str(row.SMILES)) in uniprot_full for row in output.itertuples(index=False)
    ]
    output["uniprot_overlap"] = output["uniprot_id"].astype(str).isin(uniprots)
    output["ec_overlap"] = output["ec_number"].astype(str).map(
        lambda value: any(ec in ecs for ec in value.split(";"))
    )
    output["same_uniprot_substrate_component_overlap"] = [
        (str(row.uniprot_id), str(row.SMILES)) in uniprot_components for row in output.itertuples(index=False)
    ]
    output = output.drop(columns=["sequence"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(f"Wrote CatPred overlap audit: {args.output}")
    for column in output.columns[6:]:
        print(f"{column}: {int(output[column].sum())}")


if __name__ == "__main__":
    main()
