#!/usr/bin/env python3
"""Create the experimental-kcat truth schema from BRENDA and SABIO-RK."""

from __future__ import annotations

import csv
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
FINAL = BASE / "data" / "final"
TABLES = BASE / "reports" / "tables"
BRENDA_RAW = BASE / "data" / "raw" / "brenda" / "brenda_kcat_raw.csv"
SABIO_RAW = BASE / "data" / "raw" / "sabiork" / "sabiork_kcat_raw.csv"

TRUTH_FIELDS = [
    "entry_id",
    "species",
    "reaction_id",
    "gene_id",
    "uniprot_id",
    "ec_number",
    "substrate_name",
    "substrate_smiles",
    "true_kcat",
    "true_kcat_log10",
    "unit",
    "pH",
    "temperature_c",
    "source_database",
    "match_level",
    "reference",
    "n_measurements",
]


def write_empty_truth() -> None:
    FINAL.mkdir(parents=True, exist_ok=True)
    path = FINAL / "experimental_kcat_truth.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRUTH_FIELDS, lineterminator="\n")
        writer.writeheader()
    print(f"Wrote empty experimental truth schema to {path}")


def raw_source_status(raw_path: Path) -> tuple[bool, str]:
    if raw_path.exists():
        with raw_path.open(newline="", encoding="utf-8") as handle:
            record_count = sum(1 for _ in csv.DictReader(handle))
        return True, f"raw_records={record_count}"

    source_dir = raw_path.parent
    if source_dir.exists() and any(source_dir.iterdir()):
        return True, f"expected_file_missing={raw_path.name}"
    return False, "raw_source_not_found"


def write_truth_status() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / "experimental_truth_status.csv"
    brenda_detected, brenda_status = raw_source_status(BRENDA_RAW)
    sabiork_detected, sabiork_status = raw_source_status(SABIO_RAW)
    rows = [
        {
            "source": "BRENDA",
            "raw_file_detected": str(brenda_detected),
            "status": brenda_status,
        },
        {
            "source": "SABIO-RK",
            "raw_file_detected": str(sabiork_detected),
            "status": sabiork_status,
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "raw_file_detected", "status"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote truth-source status to {path}")


def main() -> None:
    truth_path = FINAL / "experimental_kcat_truth.csv"
    if truth_path.exists() and (BRENDA_RAW.exists() or SABIO_RAW.exists()):
        print(f"Keeping existing experimental truth table at {truth_path}")
        write_truth_status()
        return
    write_empty_truth()
    write_truth_status()


if __name__ == "__main__":
    main()
