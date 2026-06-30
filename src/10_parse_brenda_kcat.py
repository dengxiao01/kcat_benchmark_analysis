#!/usr/bin/env python3
"""Parse BRENDA turnover numbers and rebuild experimental kcat truth.

The BRENDA JSON field ``turnover_number`` is treated as experimental kcat.
Rows are filtered to the model organisms and matched to model entries by
EC number, substrate identifiers/name, and UniProt where available.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sqlite3
import statistics
import tarfile
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
INTERIM = BASE / "data" / "interim"
RAW = BASE / "data" / "raw"
FINAL = BASE / "data" / "final"
TABLES = BASE / "reports" / "tables"

BRENDA_JSON_TAR = RAW / "brenda" / "brenda_2026_1.json.tar.gz"
BRENDA_RAW = RAW / "brenda" / "brenda_kcat_raw.csv"
SABIO_RAW = RAW / "sabiork" / "sabiork_kcat_raw.csv"
CKB_DB = RAW / "compounds" / "ckb" / "compounds.sqlite"
SMILES_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence_smiles.csv"
SEQUENCE_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence.csv"
BASE_ENTRIES = INTERIM / "enzyme_reaction_entries.csv"
QUERY_QUEUE = INTERIM / "brenda_sabiork_query_queue.csv"
TRUTH = FINAL / "experimental_kcat_truth.csv"

MODEL_ORGANISMS = {
    "ecoli": "Escherichia coli",
    "yeast": "Saccharomyces cerevisiae",
}

RAW_FIELDS = [
    "source_database",
    "source_record_id",
    "species",
    "organism",
    "ec_number",
    "substrate_name",
    "substrate_bigg_ids",
    "substrate_metanetx_ids",
    "substrate_kegg_ids",
    "substrate_chebi_ids",
    "enzyme_uniprot_ids",
    "kcat",
    "kcat_log10",
    "unit",
    "pH",
    "temperature_c",
    "comment",
    "reference",
]

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

MATCH_STRENGTH = {
    "species_ec_uniprot_substrate_id": 5,
    "species_ec_substrate_id": 4,
    "species_ec_uniprot_substrate_name": 3,
    "species_ec_substrate_name": 2,
}


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def split_values(value: str) -> list[str]:
    values = []
    for item in str(value or "").split(";"):
        item = item.strip()
        if item and item.lower() != "nan":
            values.append(item)
    return values


def fmt_float(value: float | None, digits: int = 6) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def parse_float(value: str) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def median_or_none(values: list[float | None]) -> float | None:
    clean = [value for value in values if value is not None and math.isfinite(value)]
    if not clean:
        return None
    return statistics.median(clean)


def normalize_name(value: str) -> str:
    text = str(value or "").lower().strip()
    text = re.sub(r"\s+c\d+h[\da-z]*n?\d*o?\d*s?\d*p?\d*$", "", text)
    text = re.sub(r"\s+\[[^\]]+\]$", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def organism_matches(brenda_organism: str, target_organism: str) -> bool:
    value = brenda_organism.strip()
    return value == target_organism or value.startswith(target_organism + " ")


def entry_source_path() -> Path:
    if SMILES_ENTRIES.exists():
        return SMILES_ENTRIES
    if SEQUENCE_ENTRIES.exists():
        return SEQUENCE_ENTRIES
    return BASE_ENTRIES


def relevant_ecs(entries: list[dict[str, str]]) -> set[str]:
    ecs = set()
    for row in entries:
        for ec_number in split_values(row.get("ec_number", "")):
            if "-" not in ec_number:
                ecs.add(ec_number)
    return ecs


def load_brenda_json(path: Path) -> dict:
    with tarfile.open(path, "r:gz") as tar:
        member = tar.getmembers()[0]
        handle = tar.extractfile(member)
        if handle is None:
            raise RuntimeError(f"Cannot read {member.name} from {path}")
        return json.load(handle)["data"]


def parse_turnover_value(value: str) -> tuple[float | None, str]:
    substrate = ""
    substrate_match = re.search(r"\{([^{}]+)\}", value)
    if substrate_match:
        substrate = substrate_match.group(1).strip()
    numeric_part = value.split("{", 1)[0]
    numbers = [float(item) for item in re.findall(r"(?<![A-Za-z])-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", numeric_part)]
    if not numbers:
        return None, substrate
    if len(numbers) >= 2 and re.search(r"\d\s*[-–]\s*\d", numeric_part):
        return statistics.mean(numbers[:2]), substrate
    return numbers[0], substrate


def parse_condition(comment: str, kind: str) -> str:
    if kind == "pH":
        match = re.search(r"\bpH\s*([0-9]+(?:\.[0-9]+)?)", comment, flags=re.IGNORECASE)
    else:
        match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:°\s*)?C\b", comment)
    return match.group(1) if match else ""


def is_mutant_comment(comment: str) -> bool:
    return bool(re.search(r"\b(mutant|mutation|variant)\b", comment, flags=re.IGNORECASE))


def ckb_name_candidates(name: str) -> list[str]:
    candidates = [name.strip()]
    if name:
        candidates.extend([name.lower(), name.upper(), name.title()])
    seen = set()
    output = []
    for item in candidates:
        if item and item not in seen:
            output.append(item)
            seen.add(item)
    return output


def map_substrate_names_with_ckb(substrate_names: set[str], db_path: Path) -> dict[str, dict[str, str]]:
    if not db_path.exists() or not substrate_names:
        return {}
    lookup_values = sorted({candidate for name in substrate_names for candidate in ckb_name_candidates(name)})
    conn = sqlite3.connect(db_path)
    try:
        name_to_compounds: dict[str, set[int]] = {}
        for start in range(0, len(lookup_values), 500):
            chunk = lookup_values[start : start + 500]
            placeholders = ",".join(["?"] * len(chunk))
            sql = f"""
                select accession, compound_id
                from compound_identifiers
                where registry_id = 3
                  and accession in ({placeholders})
            """
            for accession, compound_id in conn.execute(sql, chunk):
                name_to_compounds.setdefault(normalize_name(accession), set()).add(compound_id)

        compound_ids = sorted({cid for ids in name_to_compounds.values() for cid in ids})
        compound_to_ids: dict[int, dict[str, set[str]]] = {
            cid: {"bigg": set(), "metanetx": set(), "kegg": set(), "chebi": set()} for cid in compound_ids
        }
        registry_to_field = {
            13: "bigg",
            4: "metanetx",
            8: "kegg",
            7: "chebi",
        }
        for start in range(0, len(compound_ids), 500):
            chunk = compound_ids[start : start + 500]
            placeholders = ",".join(["?"] * len(chunk))
            sql = f"""
                select compound_id, registry_id, accession
                from compound_identifiers
                where registry_id in (13, 4, 8, 7)
                  and compound_id in ({placeholders})
            """
            for compound_id, registry_id, accession in conn.execute(sql, chunk):
                field = registry_to_field.get(registry_id)
                if field:
                    compound_to_ids.setdefault(compound_id, {"bigg": set(), "metanetx": set(), "kegg": set(), "chebi": set()})[field].add(accession)

        output = {}
        for name in substrate_names:
            ids = name_to_compounds.get(normalize_name(name), set())
            merged = {"bigg": set(), "metanetx": set(), "kegg": set(), "chebi": set()}
            for compound_id in ids:
                for field, values in compound_to_ids.get(compound_id, {}).items():
                    merged[field].update(values)
            output[name] = {field: ";".join(sorted(values)) for field, values in merged.items()}
        return output
    finally:
        conn.close()


def reference_text(ec_number: str, index: int, refs: list[str], reference_table: dict) -> str:
    parts = []
    for ref_id in refs:
        ref = reference_table.get(ref_id, {})
        token = f"BRENDA:{ec_number}:TN:{index}:ref:{ref_id}"
        pmid = ref.get("pmid")
        if pmid:
            token += f"|PubMed:{pmid}"
        parts.append(token)
    return ";".join(parts) if parts else f"BRENDA:{ec_number}:TN:{index}"


def parse_brenda_records(entries: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    data = load_brenda_json(Path(args.brenda_json))
    ecs = relevant_ecs(entries)
    substrate_names: set[str] = set()
    intermediate = []

    for ec_number in sorted(ecs):
        ec_data = data.get(ec_number)
        if not ec_data:
            continue
        proteins = ec_data.get("protein", {})
        references = ec_data.get("reference", {})
        for index, row in enumerate(ec_data.get("turnover_number", []), start=1):
            kcat, substrate = parse_turnover_value(str(row.get("value", "")))
            if kcat is None or kcat <= 0 or not substrate:
                continue
            comment = row.get("comment", "")
            if not args.include_mutants and is_mutant_comment(comment):
                continue
            for protein_id in row.get("proteins", []):
                protein = proteins.get(str(protein_id), {})
                organism = protein.get("organism", "")
                species = next((sid for sid, target in MODEL_ORGANISMS.items() if organism_matches(organism, target)), "")
                if not species:
                    continue
                accessions = protein.get("accessions", []) or []
                intermediate.append(
                    {
                        "source_database": "BRENDA",
                        "source_record_id": f"{ec_number}:TN:{index}:protein:{protein_id}",
                        "species": species,
                        "organism": organism,
                        "ec_number": ec_number,
                        "substrate_name": substrate,
                        "enzyme_uniprot_ids": ";".join(sorted(set(accessions))),
                        "kcat": fmt_float(kcat),
                        "kcat_log10": fmt_float(math.log10(kcat)),
                        "unit": "s^-1",
                        "pH": parse_condition(comment, "pH"),
                        "temperature_c": parse_condition(comment, "temperature"),
                        "comment": comment,
                        "reference": reference_text(ec_number, index, row.get("references", []), references),
                    }
                )
                substrate_names.add(substrate)

    name_to_ids = map_substrate_names_with_ckb(substrate_names, Path(args.ckb_db))
    records = []
    for row in intermediate:
        ids = name_to_ids.get(row["substrate_name"], {})
        row["substrate_bigg_ids"] = ids.get("bigg", "")
        row["substrate_metanetx_ids"] = ids.get("metanetx", "")
        row["substrate_kegg_ids"] = ids.get("kegg", "")
        row["substrate_chebi_ids"] = ids.get("chebi", "")
        records.append(row)
    return records


def read_sabio_as_generic() -> list[dict[str, str]]:
    rows = []
    for row in read_rows(SABIO_RAW):
        rows.append(
            {
                "source_database": "SABIO-RK",
                "source_record_id": row.get("sabiork_entry_id", ""),
                "species": row.get("query_species", ""),
                "organism": row.get("query_organism", ""),
                "ec_number": row.get("query_ec_number", ""),
                "substrate_name": row.get("participant_names", ""),
                "substrate_bigg_ids": "",
                "substrate_metanetx_ids": "",
                "substrate_kegg_ids": row.get("participant_kegg_ids", ""),
                "substrate_chebi_ids": row.get("participant_chebi_ids", ""),
                "enzyme_uniprot_ids": row.get("enzyme_uniprot_ids", ""),
                "kcat": row.get("kcat", ""),
                "kcat_log10": row.get("kcat_log10", ""),
                "unit": row.get("unit", "s^-1"),
                "pH": row.get("pH", ""),
                "temperature_c": row.get("temperature_c", ""),
                "comment": row.get("buffer", ""),
                "reference": row.get("source_url", ""),
            }
        )
    return rows


def substrate_match(entry: dict[str, str], record: dict[str, str]) -> str:
    entry_ids = (
        set(split_values(entry.get("substrate_kegg_id", "")))
        | set(split_values(entry.get("substrate_chebi_id", "")))
        | set(split_values(entry.get("substrate_bigg_id", "")))
        | set(split_values(entry.get("substrate_metanetx_id", "")))
    )
    record_ids = (
        set(split_values(record.get("substrate_kegg_ids", "")))
        | set(split_values(record.get("substrate_chebi_ids", "")))
        | set(split_values(record.get("substrate_bigg_ids", "")))
        | set(split_values(record.get("substrate_metanetx_ids", "")))
    )
    id_match = bool(entry_ids & record_ids)

    entry_uniprot = set(split_values(entry.get("uniprot_id", "")))
    record_uniprot = set(split_values(record.get("enzyme_uniprot_ids", "")))
    uniprot_match = bool(entry_uniprot & record_uniprot)

    entry_name = normalize_name(entry.get("substrate_name", ""))
    record_names = [normalize_name(name) for name in split_values(record.get("substrate_name", ""))]
    name_match = bool(entry_name) and entry_name in record_names

    if id_match and uniprot_match:
        return "species_ec_uniprot_substrate_id"
    if id_match:
        return "species_ec_substrate_id"
    if name_match and uniprot_match:
        return "species_ec_uniprot_substrate_name"
    if name_match:
        return "species_ec_substrate_name"
    return ""


def build_combined_truth(entries: list[dict[str, str]], records: list[dict[str, str]]) -> list[dict[str, str]]:
    records_by_species_ec: dict[tuple[str, str], list[dict[str, str]]] = {}
    for record in records:
        records_by_species_ec.setdefault((record["species"], record["ec_number"]), []).append(record)

    truth_rows = []
    for entry in entries:
        matches = []
        for ec_number in split_values(entry.get("ec_number", "")):
            for record in records_by_species_ec.get((entry["species"], ec_number), []):
                level = substrate_match(entry, record)
                if not level:
                    continue
                matches.append((MATCH_STRENGTH[level], level, record))
        if not matches:
            continue
        max_strength = max(item[0] for item in matches)
        selected = [item for item in matches if item[0] == max_strength]
        kcats = [parse_float(item[2]["kcat"]) for item in selected]
        kcats = [value for value in kcats if value is not None and value > 0]
        if not kcats:
            continue
        median_kcat = statistics.median(kcats)
        p_h = median_or_none([parse_float(item[2].get("pH", "")) for item in selected])
        temperature = median_or_none([parse_float(item[2].get("temperature_c", "")) for item in selected])
        sources = sorted({item[2]["source_database"] for item in selected})
        references = sorted({ref for item in selected for ref in split_values(item[2].get("reference", ""))})
        levels = sorted({item[1] for item in selected}, key=lambda level: -MATCH_STRENGTH[level])
        truth_rows.append(
            {
                "entry_id": entry["entry_id"],
                "species": entry["species"],
                "reaction_id": entry["reaction_id"],
                "gene_id": entry.get("gene_id", ""),
                "uniprot_id": entry.get("uniprot_id", ""),
                "ec_number": entry.get("ec_number", ""),
                "substrate_name": entry.get("substrate_name", ""),
                "substrate_smiles": entry.get("substrate_smiles", ""),
                "true_kcat": fmt_float(median_kcat),
                "true_kcat_log10": fmt_float(math.log10(median_kcat)),
                "unit": "s^-1",
                "pH": fmt_float(p_h),
                "temperature_c": fmt_float(temperature),
                "source_database": ";".join(sources),
                "match_level": ";".join(levels),
                "reference": ";".join(references),
                "n_measurements": str(len(selected)),
            }
        )
    return truth_rows


def write_reports(brenda_records: list[dict[str, str]], combined_records: list[dict[str, str]], truth_rows: list[dict[str, str]]) -> None:
    by_species: dict[str, dict[str, int]] = {}
    for row in brenda_records:
        stats = by_species.setdefault(row["species"], {"brenda_raw_records": 0, "combined_raw_records": 0, "truth_entries": 0})
        stats["brenda_raw_records"] += 1
    for row in combined_records:
        stats = by_species.setdefault(row["species"], {"brenda_raw_records": 0, "combined_raw_records": 0, "truth_entries": 0})
        stats["combined_raw_records"] += 1
    for row in truth_rows:
        stats = by_species.setdefault(row["species"], {"brenda_raw_records": 0, "combined_raw_records": 0, "truth_entries": 0})
        stats["truth_entries"] += 1
    write_rows(
        TABLES / "brenda_kcat_summary.csv",
        [{"species": species, **{key: str(value) for key, value in stats.items()}} for species, stats in sorted(by_species.items())],
        ["species", "brenda_raw_records", "combined_raw_records", "truth_entries"],
    )

    status_rows = [
        {
            "source": "BRENDA",
            "raw_file_detected": str(BRENDA_RAW.exists()),
            "status": f"parsed_raw_records={len(brenda_records)}",
        },
        {
            "source": "SABIO-RK",
            "raw_file_detected": str(SABIO_RAW.exists()),
            "status": f"raw_records={sum(1 for row in combined_records if row['source_database'] == 'SABIO-RK')}",
        },
    ]
    write_rows(TABLES / "experimental_truth_status.csv", status_rows, ["source", "raw_file_detected", "status"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--brenda-json", default=str(BRENDA_JSON_TAR))
    parser.add_argument("--ckb-db", default=str(CKB_DB))
    parser.add_argument("--include-mutants", action="store_true")
    args = parser.parse_args()

    entries = read_rows(entry_source_path())
    brenda_records = parse_brenda_records(entries, args)
    write_rows(BRENDA_RAW, brenda_records, RAW_FIELDS)
    combined_records = brenda_records + read_sabio_as_generic()
    truth_rows = build_combined_truth(entries, combined_records)
    write_rows(TRUTH, truth_rows, TRUTH_FIELDS)
    write_reports(brenda_records, combined_records, truth_rows)

    print(f"Wrote {len(brenda_records)} BRENDA raw kcat records to {BRENDA_RAW}")
    print(f"Wrote {len(truth_rows)} combined experimental truth rows to {TRUTH}")


if __name__ == "__main__":
    main()
