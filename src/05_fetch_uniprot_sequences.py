#!/usr/bin/env python3
"""Fetch protein sequences from UniProt and update Phase 1 entry tables."""

from __future__ import annotations

import argparse
import csv
import time
import urllib.parse
import urllib.request
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
QUEUE = BASE / "data" / "interim" / "uniprot_sequence_queue.csv"
ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries.csv"
RAW_FASTA = BASE / "data" / "raw" / "uniprot_sequences.fasta"
SEQUENCE_CSV = BASE / "data" / "interim" / "uniprot_sequences.csv"
UPDATED_ENTRIES = BASE / "data" / "interim" / "enzyme_reaction_entries_with_sequence.csv"
TABLES = BASE / "reports" / "tables"
API = "https://rest.uniprot.org/uniprotkb/accessions"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def chunks(items: list[str], size: int) -> list[list[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def parse_fasta(text: str) -> dict[str, tuple[str, str]]:
    result = {}
    header = ""
    accession = ""
    seq_parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if accession:
                result[accession] = (header, "".join(seq_parts))
            header = line[1:]
            parts = header.split("|")
            accession = parts[1] if len(parts) >= 2 else header.split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if accession:
        result[accession] = (header, "".join(seq_parts))
    return result


def fetch_chunk(accessions: list[str], timeout: int) -> str:
    query = urllib.parse.urlencode({"accessions": ",".join(accessions), "format": "fasta"})
    url = f"{API}?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "kcat-benchmark-analysis/0.1"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def split_values(value: str) -> list[str]:
    return [item for item in (value or "").split(";") if item]


def update_entries(entries: list[dict[str, str]], sequences: dict[str, tuple[str, str]]) -> list[dict[str, str]]:
    updated = []
    for row in entries:
        uniprots = split_values(row.get("uniprot_id", ""))
        if len(uniprots) == 1 and uniprots[0] in sequences:
            row["protein_sequence"] = sequences[uniprots[0]][1]
            row["sequence_status"] = "sequence_fetched"
        elif len(uniprots) > 1 and all(uid in sequences for uid in uniprots):
            row["sequence_status"] = "complex_sequences_available"
        elif uniprots:
            row["sequence_status"] = "sequence_missing_from_uniprot_fetch"
        else:
            row["sequence_status"] = "missing_uniprot"
        updated.append(row)
    return updated


def refresh_coverage_tables(updated_entries: list[dict[str, str]]) -> None:
    by_species: dict[str, dict[str, int]] = {}
    for row in updated_entries:
        species = row["species"]
        stats = by_species.setdefault(
            species,
            {
                "enzyme_substrate_entries": 0,
                "entries_with_single_sequence": 0,
                "complex_entries_with_sequences_available": 0,
                "missing_uniprot": 0,
                "sequence_missing_from_uniprot_fetch": 0,
                "entries_with_smiles": 0,
            },
        )
        stats["enzyme_substrate_entries"] += 1
        status = row.get("sequence_status", "")
        if status == "sequence_fetched":
            stats["entries_with_single_sequence"] += 1
        elif status == "complex_sequences_available":
            stats["complex_entries_with_sequences_available"] += 1
        elif status == "missing_uniprot":
            stats["missing_uniprot"] += 1
        elif status == "sequence_missing_from_uniprot_fetch":
            stats["sequence_missing_from_uniprot_fetch"] += 1
        if row.get("substrate_smiles"):
            stats["entries_with_smiles"] += 1

    sequence_rows = []
    for species, stats in sorted(by_species.items()):
        sequence_rows.append({"species": species, **{k: str(v) for k, v in stats.items()}})
    write_rows(
        TABLES / "sequence_coverage_by_species.csv",
        sequence_rows,
        [
            "species",
            "enzyme_substrate_entries",
            "entries_with_single_sequence",
            "complex_entries_with_sequences_available",
            "missing_uniprot",
            "sequence_missing_from_uniprot_fetch",
            "entries_with_smiles",
        ],
    )

    parse_summary_path = TABLES / "model_parse_summary.csv"
    if parse_summary_path.exists():
        parse_summary = read_rows(parse_summary_path)
        stages = [
            ("total_reactions", "All model reactions"),
            ("reactions_with_gpr", "Reactions with GPR"),
            ("reactions_with_ec", "Reactions with EC number"),
            ("enzyme_substrate_entries", "Enzyme-substrate entries"),
            ("entries_with_uniprot", "Entries with UniProt"),
            ("entries_with_sequence", "Entries with single protein sequence"),
            ("entries_with_smiles", "Entries with substrate SMILES"),
            ("experimental_kcat_matched", "Entries matched to experimental kcat"),
        ]
        coverage_rows = []
        for row in parse_summary:
            species = row["species"]
            row = dict(row)
            row["entries_with_sequence"] = str(by_species.get(species, {}).get("entries_with_single_sequence", 0))
            row["entries_with_smiles"] = str(by_species.get(species, {}).get("entries_with_smiles", 0))
            for key, label in stages:
                coverage_rows.append({"species": species, "stage": key, "stage_label": label, "n": row[key]})
        write_rows(TABLES / "experimental_kcat_coverage.csv", coverage_rows, ["species", "stage", "stage_label", "n"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--sleep", type=float, default=0.2)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--limit", type=int, default=0, help="Fetch only the first N IDs; 0 means all.")
    parser.add_argument("--use-cache", action="store_true", help="Use existing uniprot_sequences.csv instead of downloading.")
    args = parser.parse_args()

    queue_rows = read_rows(QUEUE)
    accessions = [row["uniprot_id"] for row in queue_rows if row.get("uniprot_id")]
    if args.limit:
        accessions = accessions[: args.limit]

    sequences: dict[str, tuple[str, str]] = {}
    fasta_blocks = []
    failures = []
    if args.use_cache:
        for row in read_rows(SEQUENCE_CSV):
            sequences[row["uniprot_id"]] = (row.get("fasta_header", ""), row.get("protein_sequence", ""))
        print(f"Loaded {len(sequences)} cached sequences from {SEQUENCE_CSV}")
    else:
        for index, chunk in enumerate(chunks(accessions, args.chunk_size), start=1):
            try:
                text = fetch_chunk(chunk, args.timeout)
                fasta_blocks.append(text.strip())
                sequences.update(parse_fasta(text))
                print(f"Fetched chunk {index}: requested={len(chunk)}, cumulative_sequences={len(sequences)}")
            except Exception as exc:  # noqa: BLE001 - keep fetch robust and report all failures.
                failures.append({"chunk_index": str(index), "accessions": ";".join(chunk), "error": str(exc)})
                print(f"Failed chunk {index}: {exc}")
            time.sleep(args.sleep)

        RAW_FASTA.parent.mkdir(parents=True, exist_ok=True)
        RAW_FASTA.write_text("\n".join(block for block in fasta_blocks if block) + "\n", encoding="utf-8")

    sequence_rows = [
        {
            "uniprot_id": accession,
            "fasta_header": sequences[accession][0],
            "protein_sequence": sequences[accession][1],
            "sequence_length": str(len(sequences[accession][1])),
        }
        for accession in sorted(sequences)
    ]
    write_rows(SEQUENCE_CSV, sequence_rows, ["uniprot_id", "fasta_header", "protein_sequence", "sequence_length"])

    updated_entries = update_entries(read_rows(ENTRIES), sequences)
    write_rows(UPDATED_ENTRIES, updated_entries, list(updated_entries[0].keys()))
    refresh_coverage_tables(updated_entries)

    queue_out = []
    for row in queue_rows:
        uid = row["uniprot_id"]
        out = dict(row)
        out["sequence_status"] = "fetched" if uid in sequences else "missing_or_not_requested"
        queue_out.append(out)
    write_rows(QUEUE, queue_out, list(queue_out[0].keys()))

    write_rows(
        TABLES / "uniprot_fetch_failures.csv",
        failures,
        ["chunk_index", "accessions", "error"],
    )
    summary = [
        {
            "requested_uniprot_ids": str(len(accessions)),
            "fetched_sequences": str(len(sequences)),
            "failed_chunks": str(len(failures)),
            "updated_entries_with_single_sequence": str(
                sum(row["sequence_status"] == "sequence_fetched" for row in updated_entries)
            ),
            "complex_entries_with_sequences_available": str(
                sum(row["sequence_status"] == "complex_sequences_available" for row in updated_entries)
            ),
        }
    ]
    write_rows(
        TABLES / "uniprot_fetch_summary.csv",
        summary,
        [
            "requested_uniprot_ids",
            "fetched_sequences",
            "failed_chunks",
            "updated_entries_with_single_sequence",
            "complex_entries_with_sequences_available",
        ],
    )
    print(f"Wrote FASTA to {RAW_FASTA}")
    print(f"Wrote sequence table to {SEQUENCE_CSV}")
    print(f"Wrote updated entries to {UPDATED_ENTRIES}")


if __name__ == "__main__":
    main()
