#!/usr/bin/env python3
"""Build unique curation queues for sequences, SMILES, and kcat matching."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries.csv"
INTERIM = BASE / "data" / "interim"
TABLES = BASE / "reports" / "tables"


def split_values(value: str) -> list[str]:
    return [item for item in (value or "").split(";") if item]


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_entries() -> list[dict[str, str]]:
    return read_rows(ENTRIES)


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_uniprot_queue(entries: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in entries:
        for uniprot_id in split_values(row.get("uniprot_id", "")):
            grouped[uniprot_id]["species"].add(row["species"])
            grouped[uniprot_id]["entry_id"].add(row["entry_id"])
            grouped[uniprot_id]["gpr_group_id"].add(row["gpr_group_id"])
    rows = []
    for uniprot_id, values in sorted(grouped.items()):
        rows.append(
            {
                "uniprot_id": uniprot_id,
                "species": ";".join(sorted(values["species"])),
                "n_entries": str(len(values["entry_id"])),
                "n_gpr_groups": str(len(values["gpr_group_id"])),
                "sequence_status": "needs_fetch",
            }
        )
    return rows


def build_substrate_queue(entries: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    scalar = {}
    for row in entries:
        key = (row["species"], row["substrate_id"])
        scalar[key] = {
            "substrate_name": row.get("substrate_name", ""),
            "substrate_bigg_id": row.get("substrate_bigg_id", ""),
            "substrate_kegg_id": row.get("substrate_kegg_id", ""),
            "substrate_chebi_id": row.get("substrate_chebi_id", ""),
            "substrate_metanetx_id": row.get("substrate_metanetx_id", ""),
            "substrate_is_cofactor_like": row.get("substrate_is_cofactor_like", ""),
            "substrate_role_class": row.get("substrate_role_class", ""),
            "substrate_role_evidence": row.get("substrate_role_evidence", ""),
            "substrate_role_registry_name": row.get("substrate_role_registry_name", ""),
        }
        grouped[key]["entry_id"].add(row["entry_id"])
        grouped[key]["reaction_id"].add(row["reaction_id"])
    rows = []
    for key, values in sorted(grouped.items()):
        species, substrate_id = key
        info = scalar[key]
        rows.append(
            {
                "species": species,
                "substrate_id": substrate_id,
                "substrate_name": info["substrate_name"],
                "substrate_bigg_id": info["substrate_bigg_id"],
                "substrate_kegg_id": info["substrate_kegg_id"],
                "substrate_chebi_id": info["substrate_chebi_id"],
                "substrate_metanetx_id": info["substrate_metanetx_id"],
                "substrate_is_cofactor_like": info["substrate_is_cofactor_like"],
                "substrate_role_class": info["substrate_role_class"],
                "substrate_role_evidence": info["substrate_role_evidence"],
                "substrate_role_registry_name": info["substrate_role_registry_name"],
                "n_reactions": str(len(values["reaction_id"])),
                "n_entries": str(len(values["entry_id"])),
                "smiles_status": "needs_mapping",
            }
        )
    return rows


def restore_existing_smiles(
    rows: list[dict[str, str]], existing: list[dict[str, str]]
) -> list[dict[str, str]]:
    by_key = {(row.get("species", ""), row.get("substrate_id", "")): row for row in existing}
    carry = ["substrate_smiles", "smiles_source", "smiles_source_id", "smiles_status"]
    for row in rows:
        previous = by_key.get((row.get("species", ""), row.get("substrate_id", "")))
        if not previous:
            continue
        for column in carry:
            if previous.get(column, ""):
                row[column] = previous[column]
    return rows


def build_kcat_query_queue(entries: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    names = {}
    for row in entries:
        for ec_number in split_values(row.get("ec_number", "")):
            if "-" in ec_number:
                continue
            key = (row["species"], ec_number, row.get("substrate_name", ""))
            names[key] = {
                "organism": row.get("organism", ""),
                "substrate_id": row.get("substrate_id", ""),
                "substrate_bigg_id": row.get("substrate_bigg_id", ""),
                "substrate_kegg_id": row.get("substrate_kegg_id", ""),
                "substrate_chebi_id": row.get("substrate_chebi_id", ""),
            }
            grouped[key]["entry_id"].add(row["entry_id"])
            grouped[key]["reaction_id"].add(row["reaction_id"])
            grouped[key]["uniprot_id"].update(split_values(row.get("uniprot_id", "")))
    rows = []
    for key, values in sorted(grouped.items()):
        species, ec_number, substrate_name = key
        info = names[key]
        rows.append(
            {
                "species": species,
                "organism": info["organism"],
                "ec_number": ec_number,
                "substrate_name": substrate_name,
                "substrate_id": info["substrate_id"],
                "substrate_bigg_id": info["substrate_bigg_id"],
                "substrate_kegg_id": info["substrate_kegg_id"],
                "substrate_chebi_id": info["substrate_chebi_id"],
                "uniprot_id": ";".join(sorted(values["uniprot_id"])),
                "n_reactions": str(len(values["reaction_id"])),
                "n_entries": str(len(values["entry_id"])),
                "match_status": "needs_brenda_sabiork_query",
            }
        )
    return rows


def main() -> None:
    entries = read_entries()
    previous_substrate_rows = read_rows(INTERIM / "substrate_smiles_queue.csv")
    uniprot_rows = build_uniprot_queue(entries)
    substrate_rows = restore_existing_smiles(
        build_substrate_queue(entries), previous_substrate_rows
    )
    kcat_rows = build_kcat_query_queue(entries)

    write_rows(
        INTERIM / "uniprot_sequence_queue.csv",
        uniprot_rows,
        ["uniprot_id", "species", "n_entries", "n_gpr_groups", "sequence_status"],
    )
    write_rows(
        INTERIM / "substrate_smiles_queue.csv",
        substrate_rows,
        [
            "species",
            "substrate_id",
            "substrate_name",
            "substrate_bigg_id",
            "substrate_kegg_id",
            "substrate_chebi_id",
            "substrate_metanetx_id",
            "substrate_is_cofactor_like",
            "substrate_role_class",
            "substrate_role_evidence",
            "substrate_role_registry_name",
            "n_reactions",
            "n_entries",
            "smiles_status",
            "substrate_smiles",
            "smiles_source",
            "smiles_source_id",
        ],
    )
    write_rows(
        INTERIM / "brenda_sabiork_query_queue.csv",
        kcat_rows,
        [
            "species",
            "organism",
            "ec_number",
            "substrate_name",
            "substrate_id",
            "substrate_bigg_id",
            "substrate_kegg_id",
            "substrate_chebi_id",
            "uniprot_id",
            "n_reactions",
            "n_entries",
            "match_status",
        ],
    )

    status_rows = [
        {"queue": "uniprot_sequence_queue", "n_rows": str(len(uniprot_rows))},
        {"queue": "substrate_smiles_queue", "n_rows": str(len(substrate_rows))},
        {"queue": "brenda_sabiork_query_queue", "n_rows": str(len(kcat_rows))},
    ]
    write_rows(TABLES / "curation_queue_summary.csv", status_rows, ["queue", "n_rows"])

    print(f"Wrote {len(uniprot_rows)} UniProt IDs to {INTERIM / 'uniprot_sequence_queue.csv'}")
    print(f"Wrote {len(substrate_rows)} substrates to {INTERIM / 'substrate_smiles_queue.csv'}")
    print(f"Wrote {len(kcat_rows)} EC-substrate queries to {INTERIM / 'brenda_sabiork_query_queue.csv'}")


if __name__ == "__main__":
    main()
