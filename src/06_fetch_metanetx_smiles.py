#!/usr/bin/env python3
"""Fill substrate SMILES using MetaNetX MNXM identifiers.

The MetaNetX chemistry property table is large, so this script streams it and
keeps only rows needed by data/interim/substrate_smiles_queue.csv.
"""

from __future__ import annotations

import csv
import time
import urllib.request
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
INTERIM = BASE / "data" / "interim"
TABLES = BASE / "reports" / "tables"
SUBSTRATE_QUEUE = INTERIM / "substrate_smiles_queue.csv"
BASE_ENTRIES = INTERIM / "enzyme_reaction_entries.csv"
SEQUENCE_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence.csv"
SMILES_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence_smiles.csv"
MNX_SUBSET = INTERIM / "metanetx_smiles_subset.csv"
MNX_URL = "https://www.metanetx.org/ftp/latest/chem_prop.tsv"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def split_values(value: str) -> list[str]:
    return [item for item in str(value or "").split(";") if item and item != "nan"]


def stream_metanetx(targets: set[str], timeout: int) -> dict[str, dict[str, str]]:
    found: dict[str, dict[str, str]] = {}
    request = urllib.request.Request(MNX_URL, headers={"User-Agent": "kcat-benchmark-analysis/0.1"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        for raw in response:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            mnxm_id = parts[0]
            if mnxm_id not in targets:
                continue
            smiles = parts[8].strip()
            if not smiles:
                continue
            found[mnxm_id] = {
                "mnxm_id": mnxm_id,
                "mnx_name": parts[1].strip() if len(parts) > 1 else "",
                "mnx_reference": parts[2].strip() if len(parts) > 2 else "",
                "formula": parts[3].strip() if len(parts) > 3 else "",
                "charge": parts[4].strip() if len(parts) > 4 else "",
                "mass": parts[5].strip() if len(parts) > 5 else "",
                "inchi": parts[6].strip() if len(parts) > 6 else "",
                "inchikey": parts[7].strip() if len(parts) > 7 else "",
                "smiles": smiles,
            }
            if len(found) == len(targets):
                break
    return found


def fill_substrate_queue(queue_rows: list[dict[str, str]], mnx_rows: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    output = []
    for row in queue_rows:
        chosen_id = ""
        smiles = ""
        for mnxm_id in split_values(row.get("substrate_metanetx_id", "")):
            if mnxm_id in mnx_rows:
                chosen_id = mnxm_id
                smiles = mnx_rows[mnxm_id]["smiles"]
                break
        out = dict(row)
        out["substrate_smiles"] = smiles
        out["smiles_source"] = "MetaNetX_chem_prop" if smiles else ""
        out["smiles_source_id"] = chosen_id
        out["smiles_status"] = "smiles_mapped" if smiles else "needs_mapping"
        output.append(out)
    return output


def fill_entries(entries: list[dict[str, str]], substrate_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    smiles_by_key = {
        (row["species"], row["substrate_id"]): row.get("substrate_smiles", "")
        for row in substrate_rows
    }
    source_by_key = {
        (row["species"], row["substrate_id"]): row.get("smiles_source", "")
        for row in substrate_rows
    }
    output = []
    for row in entries:
        key = (row["species"], row["substrate_id"])
        smiles = smiles_by_key.get(key, "")
        row["substrate_smiles"] = smiles
        row["smiles_status"] = "smiles_mapped" if smiles else "needs_smiles_mapping"
        row["smiles_source"] = source_by_key.get(key, "")
        output.append(row)
    return output


def refresh_smiles_coverage(entries: list[dict[str, str]]) -> None:
    by_species: dict[str, dict[str, int]] = {}
    for row in entries:
        stats = by_species.setdefault(
            row["species"],
            {
                "enzyme_substrate_entries": 0,
                "entries_with_smiles": 0,
                "single_gene_entries_with_sequence_and_smiles": 0,
                "complex_entries_with_smiles": 0,
            },
        )
        stats["enzyme_substrate_entries"] += 1
        if row.get("substrate_smiles"):
            stats["entries_with_smiles"] += 1
            if row.get("enzyme_complex_type") == "complex":
                stats["complex_entries_with_smiles"] += 1
            if row.get("enzyme_complex_type") == "single_gene" and row.get("protein_sequence"):
                stats["single_gene_entries_with_sequence_and_smiles"] += 1
    rows = [{"species": species, **{k: str(v) for k, v in stats.items()}} for species, stats in sorted(by_species.items())]
    write_rows(
        TABLES / "smiles_coverage_by_species.csv",
        rows,
        [
            "species",
            "enzyme_substrate_entries",
            "entries_with_smiles",
            "single_gene_entries_with_sequence_and_smiles",
            "complex_entries_with_smiles",
        ],
    )

    coverage_path = TABLES / "experimental_kcat_coverage.csv"
    if coverage_path.exists():
        coverage_rows = read_rows(coverage_path)
        smiles_counts = {row["species"]: row["entries_with_smiles"] for row in rows}
        for row in coverage_rows:
            if row["stage"] == "entries_with_smiles":
                row["n"] = smiles_counts.get(row["species"], "0")
        write_rows(coverage_path, coverage_rows, ["species", "stage", "stage_label", "n"])


def main() -> None:
    start = time.time()
    queue_rows = read_rows(SUBSTRATE_QUEUE)
    targets = {
        mnxm_id
        for row in queue_rows
        for mnxm_id in split_values(row.get("substrate_metanetx_id", ""))
        if mnxm_id.startswith("MNXM")
    }
    print(f"Need MetaNetX SMILES for {len(targets)} unique MNXM IDs")

    mnx_rows = stream_metanetx(targets, timeout=600)
    write_rows(
        MNX_SUBSET,
        list(mnx_rows.values()),
        ["mnxm_id", "mnx_name", "mnx_reference", "formula", "charge", "mass", "inchi", "inchikey", "smiles"],
    )

    filled_queue = fill_substrate_queue(queue_rows, mnx_rows)
    write_rows(SUBSTRATE_QUEUE, filled_queue, list(filled_queue[0].keys()))

    entry_source = SEQUENCE_ENTRIES if SEQUENCE_ENTRIES.exists() else BASE_ENTRIES
    filled_entries = fill_entries(read_rows(entry_source), filled_queue)
    write_rows(SMILES_ENTRIES, filled_entries, list(filled_entries[0].keys()))
    refresh_smiles_coverage(filled_entries)

    summary = [
        {
            "mnxm_requested": str(len(targets)),
            "mnxm_with_smiles": str(len(mnx_rows)),
            "substrates_total": str(len(filled_queue)),
            "substrates_with_smiles": str(sum(bool(row.get("substrate_smiles")) for row in filled_queue)),
            "entries_with_smiles": str(sum(bool(row.get("substrate_smiles")) for row in filled_entries)),
            "elapsed_seconds": f"{time.time() - start:.1f}",
        }
    ]
    write_rows(
        TABLES / "metanetx_smiles_summary.csv",
        summary,
        [
            "mnxm_requested",
            "mnxm_with_smiles",
            "substrates_total",
            "substrates_with_smiles",
            "entries_with_smiles",
            "elapsed_seconds",
        ],
    )
    print(summary[0])
    print(f"Wrote updated entries to {SMILES_ENTRIES}")


if __name__ == "__main__":
    main()
