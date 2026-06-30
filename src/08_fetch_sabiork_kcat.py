#!/usr/bin/env python3
"""Fetch and match experimental kcat values from SABIO-RK.

SABIO-RK returns SBML/XML, not a flat table. This script queries by
EC number + organism + kcat, parses the kinetic-law XML, and matches
records back to model enzyme-substrate entries.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
INTERIM = BASE / "data" / "interim"
RAW = BASE / "data" / "raw" / "sabiork"
FINAL = BASE / "data" / "final"
TABLES = BASE / "reports" / "tables"

QUERY_QUEUE = INTERIM / "brenda_sabiork_query_queue.csv"
BASE_ENTRIES = INTERIM / "enzyme_reaction_entries.csv"
SEQUENCE_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence.csv"
SMILES_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence_smiles.csv"
QUERY_CACHE = RAW / "sabiork_query_entry_ids.csv"
RAW_RECORDS = RAW / "sabiork_kcat_raw.csv"
TRUTH = FINAL / "experimental_kcat_truth.csv"

SABIO_BASE = "https://sabiork.h-its.org/sabioRestWebServices"
USER_AGENT = "kcat-benchmark-analysis/0.1"

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

RAW_FIELDS = [
    "query_species",
    "query_organism",
    "query_ec_number",
    "sabiork_entry_id",
    "reaction_ec_numbers",
    "kcat",
    "kcat_log10",
    "unit",
    "pH",
    "temperature_c",
    "buffer",
    "participant_names",
    "participant_kegg_ids",
    "participant_chebi_ids",
    "enzyme_uniprot_ids",
    "pubmed_ids",
    "source_url",
]

QUERY_CACHE_FIELDS = [
    "query_key",
    "species",
    "organism",
    "ec_number",
    "status",
    "n_entry_ids",
    "entry_ids",
    "error",
]

MATCH_STRENGTH = {
    "species_ec_uniprot_substrate_id": 5,
    "species_ec_substrate_id": 4,
    "species_ec_uniprot_substrate_name": 3,
    "species_ec_substrate_name": 2,
    "species_ec_uniprot_only": 1,
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
    output = []
    for item in str(value or "").split(";"):
        item = item.strip()
        if item and item.lower() != "nan":
            output.append(item)
    return output


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def normalize_name(value: str) -> str:
    text = str(value or "").lower().strip()
    text = re.sub(r"\s+c\d+h[\da-z]*n?\d*o?\d*s?\d*p?\d*$", "", text)
    text = re.sub(r"\s+\[[^\]]+\]$", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_float(value: str) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def fmt_float(value: float | None, digits: int = 6) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def median_or_none(values: list[float]) -> float | None:
    clean = [v for v in values if v is not None and math.isfinite(v)]
    if not clean:
        return None
    return statistics.median(clean)


def extract_resources(element: ET.Element) -> list[str]:
    resources = []
    for child in element.iter():
        for key, value in child.attrib.items():
            if local_name(key) == "resource":
                resources.append(value)
    return resources


def resource_tail(resource: str) -> str:
    return resource.rstrip("/").rsplit("/", 1)[-1]


def resource_ids(resources: list[str], token: str) -> list[str]:
    values = []
    for resource in resources:
        if token not in resource:
            continue
        tail = urllib.parse.unquote(resource_tail(resource))
        if token == "ec-code":
            tail = tail.replace("ec-code:", "")
        values.append(tail)
    return sorted(set(values))


def canonical_unit(unit: str) -> str:
    if unit == "swedgeone":
        return "s^-1"
    return unit or ""


def find_sabiork_text(element: ET.Element, tag_name: str) -> str:
    for child in element.iter():
        if local_name(child.tag) == tag_name:
            return " ".join("".join(child.itertext()).split())
    return ""


def build_species_map(root: ET.Element) -> dict[str, dict[str, object]]:
    species_map: dict[str, dict[str, object]] = {}
    for species in root.iter():
        if local_name(species.tag) != "species":
            continue
        species_id = species.attrib.get("id", "")
        resources = extract_resources(species)
        species_map[species_id] = {
            "name": species.attrib.get("name", ""),
            "resources": resources,
            "kegg": resource_ids(resources, "kegg.compound"),
            "chebi": resource_ids(resources, "chebi"),
            "uniprot": resource_ids(resources, "uniprot"),
        }
    return species_map


def reaction_species_ids(reaction: ET.Element, tag_name: str) -> list[str]:
    ids = []
    for child in reaction.iter():
        if local_name(child.tag) == tag_name:
            species_id = child.attrib.get("species", "")
            if species_id:
                ids.append(species_id)
    return ids


def parse_sabiork_sbml(xml_bytes: bytes, query: dict[str, str]) -> list[dict[str, str]]:
    root = ET.fromstring(xml_bytes)
    species_map = build_species_map(root)
    records: list[dict[str, str]] = []

    for reaction in root.iter():
        if local_name(reaction.tag) != "reaction":
            continue
        reaction_resources = extract_resources(reaction)
        reaction_ecs = resource_ids(reaction_resources, "ec-code")
        participant_ids = reaction_species_ids(reaction, "speciesReference")
        modifier_ids = reaction_species_ids(reaction, "modifierSpeciesReference")
        enzyme_ids = [sid for sid in modifier_ids if sid.startswith("ENZ_")]

        participant_names = sorted({str(species_map.get(sid, {}).get("name", "")) for sid in participant_ids})
        participant_names = [name for name in participant_names if name]
        participant_kegg = sorted(
            {
                item
                for sid in participant_ids
                for item in species_map.get(sid, {}).get("kegg", [])
            }
        )
        participant_chebi = sorted(
            {
                item
                for sid in participant_ids
                for item in species_map.get(sid, {}).get("chebi", [])
            }
        )
        enzyme_uniprots = sorted(
            {
                item
                for sid in enzyme_ids
                for item in species_map.get(sid, {}).get("uniprot", [])
            }
        )

        kinetic_law = next((child for child in reaction if local_name(child.tag) == "kineticLaw"), None)
        if kinetic_law is None:
            continue
        entry_id = find_sabiork_text(kinetic_law, "kineticLawID")
        if not entry_id:
            entry_id = kinetic_law.attrib.get("metaid", "").replace("META_KL_", "")
        p_h = find_sabiork_text(kinetic_law, "startValuepH")
        temperature = find_sabiork_text(kinetic_law, "startValueTemperature")
        buffer_text = find_sabiork_text(kinetic_law, "buffer")
        pubmed_ids = resource_ids(extract_resources(kinetic_law), "pubmed")

        for parameter in kinetic_law.iter():
            if local_name(parameter.tag) != "localParameter":
                continue
            name = parameter.attrib.get("name", "") or parameter.attrib.get("id", "")
            if name.lower() != "kcat":
                continue
            kcat = parse_float(parameter.attrib.get("value", ""))
            if kcat is None or kcat <= 0:
                continue
            records.append(
                {
                    "query_species": query["species"],
                    "query_organism": query["organism"],
                    "query_ec_number": query["ec_number"],
                    "sabiork_entry_id": entry_id,
                    "reaction_ec_numbers": ";".join(reaction_ecs),
                    "kcat": fmt_float(kcat),
                    "kcat_log10": fmt_float(math.log10(kcat)),
                    "unit": canonical_unit(parameter.attrib.get("units", "")),
                    "pH": p_h,
                    "temperature_c": temperature,
                    "buffer": buffer_text,
                    "participant_names": ";".join(participant_names),
                    "participant_kegg_ids": ";".join(participant_kegg),
                    "participant_chebi_ids": ";".join(participant_chebi),
                    "enzyme_uniprot_ids": ";".join(enzyme_uniprots),
                    "pubmed_ids": ";".join(pubmed_ids),
                    "source_url": f"{SABIO_BASE}/kineticLaws/{entry_id}",
                }
            )
    return records


def request_text(url: str, timeout: int) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", "replace")


def request_bytes(url: str, timeout: int) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def query_entry_ids(query: dict[str, str], timeout: int) -> tuple[str, list[str], str]:
    query_string = f'ECNumber:{query["ec_number"]} AND Organism:"{query["organism"]}" AND Parametertype:"kcat"'
    url = f"{SABIO_BASE}/searchKineticLaws/entryIDs?format=txt&q={urllib.parse.quote(query_string)}"
    try:
        text = request_text(url, timeout)
    except urllib.error.HTTPError as exc:
        return "error", [], f"HTTP {exc.code}"
    except Exception as exc:  # noqa: BLE001 - batch fetching should continue.
        return "error", [], str(exc)
    if "no data found" in text.lower():
        return "no_data", [], ""
    ids = re.findall(r"\b\d+\b", text)
    return "found", sorted(set(ids), key=int), ""


def load_query_cache() -> dict[str, dict[str, str]]:
    return {row["query_key"]: row for row in read_rows(QUERY_CACHE)}


def write_query_cache(cache: dict[str, dict[str, str]]) -> None:
    rows = [cache[key] for key in sorted(cache)]
    write_rows(QUERY_CACHE, rows, QUERY_CACHE_FIELDS)


def query_key(query: dict[str, str]) -> str:
    return f'{query["species"]}|{query["organism"]}|{query["ec_number"]}'


def entry_source_path() -> Path:
    if SMILES_ENTRIES.exists():
        return SMILES_ENTRIES
    if SEQUENCE_ENTRIES.exists():
        return SEQUENCE_ENTRIES
    return BASE_ENTRIES


def build_queries(
    queue_rows: list[dict[str, str]],
    entries: list[dict[str, str]],
    only_ready: bool,
    species_filter: str,
) -> list[dict[str, str]]:
    ready_ecs: set[tuple[str, str]] = set()
    if only_ready:
        for row in entries:
            if not row.get("protein_sequence"):
                continue
            if not row.get("substrate_smiles"):
                continue
            if row.get("enzyme_complex_type") != "single_gene":
                continue
            for ec_number in split_values(row.get("ec_number", "")):
                if "-" not in ec_number:
                    ready_ecs.add((row["species"], ec_number))

    queries = {}
    for row in queue_rows:
        if species_filter and row.get("species") != species_filter:
            continue
        ec_number = row.get("ec_number", "")
        if not ec_number or "-" in ec_number:
            continue
        if only_ready and (row["species"], ec_number) not in ready_ecs:
            continue
        key = (row["species"], row.get("organism", ""), ec_number)
        queries[key] = {"species": key[0], "organism": key[1], "ec_number": key[2]}
    return [queries[key] for key in sorted(queries)]


def fetch_raw_records(queries: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    RAW.mkdir(parents=True, exist_ok=True)
    cache = load_query_cache()
    all_records: list[dict[str, str]] = read_rows(RAW_RECORDS)
    processed_uncached = 0

    for index, query in enumerate(queries, start=1):
        key = query_key(query)
        cached = cache.get(key)
        if cached:
            status = cached["status"]
            ids = split_values(cached.get("entry_ids", ""))
        else:
            if args.max_ec_queries and processed_uncached >= args.max_ec_queries:
                break
            status, ids, error = query_entry_ids(query, args.timeout)
            processed_uncached += 1
            cache[key] = {
                "query_key": key,
                "species": query["species"],
                "organism": query["organism"],
                "ec_number": query["ec_number"],
                "status": status,
                "n_entry_ids": str(len(ids)),
                "entry_ids": ";".join(ids),
                "error": error,
            }
            if processed_uncached % args.save_every == 0:
                write_query_cache(cache)
            time.sleep(args.sleep)

        if status != "found" or not ids:
            continue
        if args.max_ids_per_query:
            ids = ids[: args.max_ids_per_query]
        for start in range(0, len(ids), args.batch_size):
            batch = ids[start : start + args.batch_size]
            ids_param = urllib.parse.quote(",".join(batch), safe="")
            url = f"{SABIO_BASE}/kineticLaws?kinlawids={ids_param}"
            try:
                xml_bytes = request_bytes(url, args.timeout)
                all_records.extend(parse_sabiork_sbml(xml_bytes, query))
            except Exception as exc:  # noqa: BLE001 - record the query and keep going.
                print(f"Batch fetch failed for {key} ids {batch[:3]}...: {exc}; retrying one by one.")
                for entry_id in batch:
                    single_url = f"{SABIO_BASE}/kineticLaws/{entry_id}"
                    try:
                        xml_bytes = request_bytes(single_url, args.timeout)
                        all_records.extend(parse_sabiork_sbml(xml_bytes, query))
                    except Exception as single_exc:  # noqa: BLE001 - keep the long batch robust.
                        print(f"Failed to fetch/parse {key} id {entry_id}: {single_exc}")
                    time.sleep(args.sleep)
            time.sleep(args.sleep)

        if index % args.save_every == 0:
            write_rows(RAW_RECORDS, dedupe_raw_records(all_records), RAW_FIELDS)
            print(f"Processed {index}/{len(queries)} query groups; raw records={len(all_records)}")

    write_query_cache(cache)
    records = dedupe_raw_records(all_records)
    write_rows(RAW_RECORDS, records, RAW_FIELDS)
    return records


def dedupe_raw_records(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen = set()
    output = []
    for row in rows:
        key = (
            row.get("query_species", ""),
            row.get("query_ec_number", ""),
            row.get("sabiork_entry_id", ""),
            row.get("kcat", ""),
            row.get("participant_kegg_ids", ""),
            row.get("participant_chebi_ids", ""),
        )
        if key in seen:
            continue
        seen.add(key)
        output.append(row)
    return output


def substrate_match(entry: dict[str, str], record: dict[str, str], allow_uniprot_only: bool) -> str:
    entry_kegg = set(split_values(entry.get("substrate_kegg_id", "")))
    entry_chebi = set(split_values(entry.get("substrate_chebi_id", "")))
    record_kegg = set(split_values(record.get("participant_kegg_ids", "")))
    record_chebi = set(split_values(record.get("participant_chebi_ids", "")))
    id_match = bool(entry_kegg & record_kegg) or bool(entry_chebi & record_chebi)

    entry_uniprot = set(split_values(entry.get("uniprot_id", "")))
    record_uniprot = set(split_values(record.get("enzyme_uniprot_ids", "")))
    uniprot_match = bool(entry_uniprot & record_uniprot)

    query_name = normalize_name(entry.get("substrate_name", ""))
    participant_names = [normalize_name(name) for name in split_values(record.get("participant_names", ""))]
    name_match = bool(query_name) and query_name in participant_names

    if id_match and uniprot_match:
        return "species_ec_uniprot_substrate_id"
    if id_match:
        return "species_ec_substrate_id"
    if name_match and uniprot_match:
        return "species_ec_uniprot_substrate_name"
    if name_match:
        return "species_ec_substrate_name"
    if allow_uniprot_only and uniprot_match:
        return "species_ec_uniprot_only"
    return ""


def build_truth(entries: list[dict[str, str]], records: list[dict[str, str]], allow_uniprot_only: bool) -> list[dict[str, str]]:
    records_by_species_ec: dict[tuple[str, str], list[dict[str, str]]] = {}
    for record in records:
        records_by_species_ec.setdefault((record["query_species"], record["query_ec_number"]), []).append(record)

    truth_rows = []
    for entry in entries:
        matches = []
        for ec_number in split_values(entry.get("ec_number", "")):
            for record in records_by_species_ec.get((entry["species"], ec_number), []):
                level = substrate_match(entry, record, allow_uniprot_only)
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
        references = []
        for _, _, record in selected:
            source_id = f'SABIO-RK:{record["sabiork_entry_id"]}'
            pubmed = record.get("pubmed_ids", "")
            if pubmed:
                source_id += f'|PubMed:{pubmed}'
            references.append(source_id)
        references = sorted(set(references))
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
                "source_database": "SABIO-RK",
                "match_level": ";".join(levels),
                "reference": ";".join(references),
                "n_measurements": str(len(selected)),
            }
        )
    return truth_rows


def write_reports(queries: list[dict[str, str]], raw_records: list[dict[str, str]], truth_rows: list[dict[str, str]]) -> None:
    by_species: dict[str, dict[str, int]] = {}
    query_cache_rows = read_rows(QUERY_CACHE)
    report_queries = query_cache_rows if query_cache_rows else queries
    for query in report_queries:
        stats = by_species.setdefault(query["species"], {"query_ecs": 0, "raw_records": 0, "truth_entries": 0})
        stats["query_ecs"] += 1
    for record in raw_records:
        by_species.setdefault(record["query_species"], {"query_ecs": 0, "raw_records": 0, "truth_entries": 0})
        by_species[record["query_species"]]["raw_records"] += 1
    for row in truth_rows:
        by_species.setdefault(row["species"], {"query_ecs": 0, "raw_records": 0, "truth_entries": 0})
        by_species[row["species"]]["truth_entries"] += 1
    write_rows(
        TABLES / "sabiork_fetch_summary.csv",
        [{"species": species, **{key: str(value) for key, value in stats.items()}} for species, stats in sorted(by_species.items())],
        ["species", "query_ecs", "raw_records", "truth_entries"],
    )

    status_rows = [
        {
            "source": "BRENDA",
            "raw_file_detected": str(any((BASE / "data" / "raw" / "brenda").glob("*"))),
            "status": "waiting_for_curated_export_or_api_credentials",
        },
        {
            "source": "SABIO-RK",
            "raw_file_detected": str(RAW_RECORDS.exists()),
            "status": f"fetched_raw_records={len(raw_records)};matched_truth_entries={len(truth_rows)}",
        },
    ]
    write_rows(TABLES / "experimental_truth_status.csv", status_rows, ["source", "raw_file_detected", "status"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--only-ready", action="store_true", help="Query only ECs currently ready for predictor input.")
    parser.add_argument("--species", choices=["ecoli", "yeast"], help="Only query and match one species.")
    parser.add_argument("--max-ec-queries", type=int, default=0, help="Maximum uncached EC/organism queries to send.")
    parser.add_argument("--max-ids-per-query", type=int, default=0, help="Cap SABIO entry IDs per EC/organism query.")
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--sleep", type=float, default=0.1)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--save-every", type=int, default=25)
    parser.add_argument("--match-only", action="store_true", help="Do not query SABIO-RK; rebuild truth/reports from cached raw records.")
    parser.add_argument("--allow-uniprot-only", action="store_true", help="Allow EC+organism+UniProt matches without substrate match.")
    args = parser.parse_args()

    queue_rows = read_rows(QUERY_QUEUE)
    entries = read_rows(entry_source_path())
    queries = build_queries(queue_rows, entries, only_ready=args.only_ready, species_filter=args.species or "")
    raw_records = read_rows(RAW_RECORDS) if args.match_only else fetch_raw_records(queries, args)
    truth_rows = build_truth(entries, raw_records, args.allow_uniprot_only)
    write_rows(TRUTH, truth_rows, TRUTH_FIELDS)
    write_reports(queries, raw_records, truth_rows)

    print(f"Wrote {len(raw_records)} raw SABIO-RK kcat records to {RAW_RECORDS}")
    print(f"Wrote {len(truth_rows)} matched experimental truth rows to {TRUTH}")


if __name__ == "__main__":
    main()
