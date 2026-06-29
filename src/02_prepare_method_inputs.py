#!/usr/bin/env python3
"""Prepare off-the-shelf predictor input templates from Phase 1 entries."""

from __future__ import annotations

import csv
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
BASE_ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries.csv"
SEQUENCE_ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries_with_sequence.csv"
SMILES_ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries_with_sequence_smiles.csv"
OUT_DIR = BASE / "data" / "interim" / "prediction_inputs"
LOG = BASE / "reports" / "tables" / "prediction_input_readiness.csv"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    source_entries = BASE_ENTRIES
    if SMILES_ENTRIES.exists() and SMILES_ENTRIES.stat().st_mtime >= BASE_ENTRIES.stat().st_mtime:
        source_entries = SMILES_ENTRIES
    elif SEQUENCE_ENTRIES.exists() and SEQUENCE_ENTRIES.stat().st_mtime >= BASE_ENTRIES.stat().st_mtime:
        source_entries = SEQUENCE_ENTRIES
    rows = read_rows(source_entries)
    catpred_rows = []
    readiness = []

    for row in rows:
        has_sequence = bool(row.get("protein_sequence"))
        has_smiles = bool(row.get("substrate_smiles"))
        is_single = row.get("enzyme_complex_type") == "single_gene"
        ready = has_sequence and has_smiles and is_single
        reason = []
        if not has_sequence:
            reason.append("missing_sequence")
        if not has_smiles:
            reason.append("missing_smiles")
        if not is_single:
            reason.append("complex_not_supported_by_default")
        readiness.append(
            {
                "entry_id": row["entry_id"],
                "method": "CatPred",
                "is_ready": str(ready),
                "reason": ";".join(reason),
            }
        )
        if ready:
            catpred_rows.append(
                {
                    "entry_id": row["entry_id"],
                    "SMILES": row["substrate_smiles"],
                    "sequence": row["protein_sequence"],
                    "pdbpath": row["uniprot_id"] or row["entry_id"],
                }
            )

    write_rows(OUT_DIR / "catpred_kcat_input.csv", catpred_rows, ["entry_id", "SMILES", "sequence", "pdbpath"])
    write_rows(LOG, readiness, ["entry_id", "method", "is_ready", "reason"])
    print(f"Wrote {len(catpred_rows)} CatPred-ready rows to {OUT_DIR / 'catpred_kcat_input.csv'}")
    print(f"Wrote readiness log for {len(readiness)} rows to {LOG}")


if __name__ == "__main__":
    main()
