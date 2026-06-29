#!/usr/bin/env python3
"""Create the experimental-kcat truth schema.

This script does not treat ecModel/database-fill values as primary truth.
It waits for curated BRENDA/SABIO-RK exports under data/raw/brenda/ and
data/raw/sabiork/. Until then, it writes an empty schema plus a legacy
reference table for sanity checks.
"""

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
        writer = csv.DictWriter(handle, fieldnames=TRUTH_FIELDS)
        writer.writeheader()
    print(f"Wrote empty experimental truth schema to {path}")


def write_legacy_reference() -> None:
    source = BASE / "reaction_kcat_MW_databasefill.csv"
    target = FINAL / "legacy_ecoli_kcat_reference.csv"
    if not source.exists():
        return
    with source.open("r", newline="", encoding="utf-8") as src, target.open("w", newline="", encoding="utf-8") as dst:
        reader = csv.DictReader(src)
        fieldnames = ["reaction_id", "legacy_kcat", "legacy_data_type", "legacy_use"]
        writer = csv.DictWriter(dst, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            writer.writerow(
                {
                    "reaction_id": row.get("", ""),
                    "legacy_kcat": row.get("kcat", ""),
                    "legacy_data_type": row.get("data_type", ""),
                    "legacy_use": "sanity_check_only",
                }
            )
    print(f"Wrote legacy sanity-check reference to {target}")


def write_truth_status() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / "experimental_truth_status.csv"
    rows = [
        {
            "source": "BRENDA",
            "raw_file_detected": str(any((BASE / "data" / "raw" / "brenda").glob("*"))),
            "status": "waiting_for_curated_export_or_api_credentials",
        },
        {
            "source": "SABIO-RK",
            "raw_file_detected": str(any((BASE / "data" / "raw" / "sabiork").glob("*"))),
            "status": "waiting_for_curated_export_or_api_credentials",
        },
        {
            "source": "reaction_kcat_MW_databasefill.csv",
            "raw_file_detected": str((BASE / "reaction_kcat_MW_databasefill.csv").exists()),
            "status": "legacy_sanity_check_only_not_primary_truth",
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "raw_file_detected", "status"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote truth-source status to {path}")


def main() -> None:
    truth_path = FINAL / "experimental_kcat_truth.csv"
    if truth_path.exists() and (BRENDA_RAW.exists() or SABIO_RAW.exists()):
        print(f"Keeping existing experimental truth table at {truth_path}")
        write_legacy_reference()
        write_truth_status()
        return
    write_empty_truth()
    write_legacy_reference()
    write_truth_status()


if __name__ == "__main__":
    main()
