#!/usr/bin/env python3
"""Fill missing substrate SMILES using PubChem PUG-REST.

Lookup priority:
  1. ChEBI IDs as PubChem RegistryID
  2. KEGG compound IDs as PubChem RegistryID
  3. Cleaned substrate name
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
INTERIM = BASE / "data" / "interim"
TABLES = BASE / "reports" / "tables"
SUBSTRATE_QUEUE = INTERIM / "substrate_smiles_queue.csv"
BASE_ENTRIES = INTERIM / "enzyme_reaction_entries.csv"
SEQUENCE_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence.csv"
SMILES_ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence_smiles.csv"
PUBCHEM_CACHE = INTERIM / "pubchem_smiles_cache.csv"
CKB_DB = BASE / "data" / "raw" / "compounds" / "ckb" / "compounds.sqlite"

CKB_REGISTRIES = [
    ("bigg.metabolite", "substrate_bigg_id", 13),
    ("metanetx.chemical", "substrate_metanetx_id", 4),
    ("kegg", "substrate_kegg_id", 8),
    ("chebi", "substrate_chebi_id", 7),
    ("synonyms", "substrate_name", 3),
]


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


def clean_name(name: str) -> str:
    cleaned = str(name or "").strip()
    cleaned = re.sub(r"\s+C\d+H[\dA-Za-z]*N?\d*O?\d*S?\d*P?\d*$", "", cleaned)
    cleaned = re.sub(r"\s+\[[^\]]+\]$", "", cleaned)
    cleaned = re.sub(r"\s+\(n-C\d+:\d+ACP\)$", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+\(n-C\d+:\d+\)$", "", cleaned, flags=re.IGNORECASE)
    return cleaned.strip()


def cache_key(kind: str, value: str) -> str:
    return f"{kind}:{value}"


def read_cache() -> dict[str, dict[str, str]]:
    if not PUBCHEM_CACHE.exists():
        return {}
    rows = read_rows(PUBCHEM_CACHE)
    return {cache_key(row["query_type"], row["query_value"]): row for row in rows}


def write_cache(cache: dict[str, dict[str, str]]) -> None:
    rows = list(cache.values())
    fieldnames = ["query_type", "query_value", "status", "cid", "smiles", "connectivity_smiles", "error"]
    write_rows(PUBCHEM_CACHE, rows, fieldnames)


def pubchem_url(query_type: str, query_value: str) -> str:
    value = urllib.parse.quote(query_value, safe="")
    if query_type == "name":
        return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{value}/property/CanonicalSMILES,IsomericSMILES/JSON"
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{value}/property/CanonicalSMILES,IsomericSMILES/JSON"


def fetch_pubchem(query_type: str, query_value: str, timeout: int) -> dict[str, str]:
    url = pubchem_url(query_type, query_value)
    request = urllib.request.Request(url, headers={"User-Agent": "kcat-benchmark-analysis/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        return {
            "query_type": query_type,
            "query_value": query_value,
            "status": "not_found" if exc.code == 404 else "error",
            "cid": "",
            "smiles": "",
            "connectivity_smiles": "",
            "error": f"HTTP {exc.code}",
        }
    except Exception as exc:  # noqa: BLE001 - keep long batch robust.
        return {
            "query_type": query_type,
            "query_value": query_value,
            "status": "error",
            "cid": "",
            "smiles": "",
            "connectivity_smiles": "",
            "error": str(exc),
        }

    props = payload.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return {
            "query_type": query_type,
            "query_value": query_value,
            "status": "not_found",
            "cid": "",
            "smiles": "",
            "connectivity_smiles": "",
            "error": "",
        }
    first = props[0]
    smiles = first.get("SMILES") or first.get("IsomericSMILES") or first.get("CanonicalSMILES") or ""
    return {
        "query_type": query_type,
        "query_value": query_value,
        "status": "found" if smiles else "not_found",
        "cid": str(first.get("CID", "")),
        "smiles": smiles,
        "connectivity_smiles": first.get("ConnectivitySMILES") or first.get("CanonicalSMILES") or "",
        "error": "",
    }


def query_candidates(row: dict[str, str], include_name: bool = True, max_chebi: int = 5) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []
    for kegg_id in split_values(row.get("substrate_kegg_id", "")):
        candidates.append(("registry", kegg_id))
    for chebi_id in split_values(row.get("substrate_chebi_id", ""))[:max_chebi]:
        candidates.append(("registry", chebi_id))
    name = clean_name(row.get("substrate_name", ""))
    if include_name and name:
        candidates.append(("name", name))
    seen = set()
    unique = []
    for item in candidates:
        if item not in seen:
            unique.append(item)
            seen.add(item)
    return unique


def mapping_keys(row: dict[str, str]) -> list[tuple[str, str]]:
    keys: list[tuple[str, str]] = []
    for field, label in [
        ("substrate_kegg_id", "kegg"),
        ("substrate_chebi_id", "chebi"),
        ("substrate_metanetx_id", "metanetx"),
        ("substrate_bigg_id", "bigg"),
    ]:
        for value in split_values(row.get(field, "")):
            keys.append((label, value))
    name = clean_name(row.get("substrate_name", "")).lower()
    if name:
        keys.append(("name", name))
    seen = set()
    unique = []
    for item in keys:
        if item not in seen:
            unique.append(item)
            seen.add(item)
    return unique


def ckb_accessions(row: dict[str, str], field: str) -> list[str]:
    if field != "substrate_name":
        return split_values(row.get(field, ""))
    name = str(row.get(field, "") or "").strip()
    cleaned = clean_name(name)
    no_state = re.sub(r"\s+\((?:reduced|oxidized|oxidised|no Fe\\(III\\))\)$", "", cleaned, flags=re.IGNORECASE)
    variants = [
        name,
        cleaned,
        no_state,
        cleaned.lower(),
        cleaned.upper(),
        cleaned.title(),
        no_state.lower(),
        no_state.upper(),
        no_state.title(),
    ]
    seen = set()
    output = []
    for item in variants:
        item = item.strip()
        if item and item not in seen:
            output.append(item)
            seen.add(item)
    return output


def hard_smiles_reason(row: dict[str, str]) -> str:
    name = str(row.get("substrate_name", "") or "")
    bigg_id = str(row.get("substrate_bigg_id", "") or "")
    text = f"{name} {bigg_id}"
    checks = [
        ("protein_or_redox_carrier", r"\bACP\b|acyl-carrier|\[acyl-carrier protein\]|apoprotein|flavodoxin|glutaredoxin|ferredoxin|cytochrome|thioredoxin|IscS|SufBCD|SufSE|disulfide isomerase|scaffold complex"),
        ("nucleic_acid_polymer", r"\bDNA\b|\bRNA\b|tRNA|rRNA|mRNA|oligoribonucleotide|polynucleotide"),
        ("polymer_or_glycan", r"glucan|cellulose|chitin|starch|dextrin|mannan|glycogen|poly|oligosaccharide|peptidoglycan|lipopolysaccharide|murein|^G\d{5}$"),
        ("metal_cluster_or_cofactor_complex", r"\[\dFe-|iron-sulfur|molybdenum cofactor|molybdopterin|tungsten bispterin|damaged iron"),
        ("ambiguous_or_variable_structure", r"\?|unknown|unspecified|electron acceptor|electron donor|radical"),
    ]
    for reason, pattern in checks:
        if re.search(pattern, text, flags=re.IGNORECASE):
            return reason
    return ""


def propagate_existing_mappings(queue_rows: list[dict[str, str]]) -> int:
    mapped_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in queue_rows:
        if not row.get("substrate_smiles"):
            continue
        for key in mapping_keys(row):
            mapped_by_key.setdefault(key, row)

    changed = 0
    for row in queue_rows:
        if row.get("substrate_smiles"):
            continue
        for key in mapping_keys(row):
            source = mapped_by_key.get(key)
            if not source:
                continue
            row["substrate_smiles"] = source.get("substrate_smiles", "")
            row["smiles_source"] = source.get("smiles_source", "")
            row["smiles_source_id"] = source.get("smiles_source_id", "")
            row["smiles_status"] = "smiles_mapped"
            changed += 1
            break
    return changed


def query_ckb(accessions_by_registry: dict[int, set[str]], db_path: Path) -> dict[tuple[int, str], dict[str, str]]:
    if not db_path.exists():
        return {}
    conn = sqlite3.connect(db_path)
    try:
        hits: dict[tuple[int, str], dict[str, str]] = {}
        for registry_id, accessions in accessions_by_registry.items():
            values = sorted(accessions)
            for start in range(0, len(values), 500):
                chunk = values[start : start + 500]
                if not chunk:
                    continue
                placeholders = ",".join(["?"] * len(chunk))
                sql = f"""
                    select ci.accession, ci.compound_id, c.isomeric_smiles, c.canonical_smiles
                    from compound_identifiers ci
                    join compounds c on c.id = ci.compound_id
                    where ci.registry_id = ?
                      and ci.accession in ({placeholders})
                      and (c.isomeric_smiles is not null or c.canonical_smiles is not null)
                """
                for accession, compound_id, isomeric_smiles, canonical_smiles in conn.execute(sql, [registry_id] + chunk):
                    key = (registry_id, accession)
                    if key in hits:
                        continue
                    smiles = isomeric_smiles or canonical_smiles or ""
                    if smiles:
                        hits[key] = {
                            "compound_id": str(compound_id),
                            "smiles": smiles,
                        }
        return hits
    finally:
        conn.close()


def fill_queue_from_ckb(queue_rows: list[dict[str, str]], args: argparse.Namespace) -> int:
    accessions_by_registry: dict[int, set[str]] = {}
    for row in queue_rows:
        if args.species and row.get("species") != args.species:
            continue
        if row.get("substrate_smiles"):
            continue
        if args.only_non_cofactor and row.get("substrate_is_cofactor_like") == "True":
            continue
        for _, field, registry_id in CKB_REGISTRIES:
            accessions_by_registry.setdefault(registry_id, set()).update(ckb_accessions(row, field))

    hits = query_ckb(accessions_by_registry, Path(args.ckb_db))
    changed = 0
    for row in queue_rows:
        if args.species and row.get("species") != args.species:
            continue
        if row.get("substrate_smiles"):
            continue
        if args.only_non_cofactor and row.get("substrate_is_cofactor_like") == "True":
            continue
        for namespace, field, registry_id in CKB_REGISTRIES:
            chosen = None
            chosen_accession = ""
            for accession in ckb_accessions(row, field):
                chosen = hits.get((registry_id, accession))
                if chosen:
                    chosen_accession = accession
                    break
            if not chosen:
                continue
            row["substrate_smiles"] = chosen["smiles"]
            row["smiles_source"] = "CKB_compounds_sqlite"
            row["smiles_source_id"] = f"{namespace}:{chosen_accession}|compound_id:{chosen['compound_id']}"
            row["smiles_status"] = "smiles_mapped"
            changed += 1
            break
    return changed


def fill_queue(queue_rows: list[dict[str, str]], cache: dict[str, dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    processed = 0
    changed = 0
    for row in queue_rows:
        if args.only_pubchem_candidate and hard_smiles_reason(row):
            continue
        if args.species and row.get("species") != args.species:
            continue
        if row.get("substrate_smiles"):
            continue
        if args.skip_name and row.get("smiles_status") == "pubchem_registry_not_found":
            continue
        if not args.skip_name and row.get("smiles_status") == "pubchem_not_found":
            continue
        if args.only_non_cofactor and row.get("substrate_is_cofactor_like") == "True":
            continue
        if args.limit and processed >= args.limit:
            break
        processed += 1
        chosen = None
        for query_type, query_value in query_candidates(row, include_name=not args.skip_name, max_chebi=args.max_chebi):
            key = cache_key(query_type, query_value)
            if key not in cache:
                cache[key] = fetch_pubchem(query_type, query_value, args.timeout)
                time.sleep(args.sleep)
            if cache[key].get("status") == "found" and cache[key].get("smiles"):
                chosen = cache[key]
                break
        if chosen:
            row["substrate_smiles"] = chosen["smiles"]
            row["smiles_source"] = "PubChem_PUG_REST"
            row["smiles_source_id"] = f"CID:{chosen['cid']}"
            row["smiles_status"] = "smiles_mapped"
            changed += 1
        else:
            row["substrate_smiles"] = row.get("substrate_smiles", "")
            row["smiles_source"] = row.get("smiles_source", "")
            row["smiles_source_id"] = row.get("smiles_source_id", "")
            row["smiles_status"] = "pubchem_registry_not_found" if args.skip_name else "pubchem_not_found"
        if processed % args.save_every == 0:
            write_cache(cache)
            write_rows(SUBSTRATE_QUEUE, queue_rows, list(queue_rows[0].keys()))
            print(f"Processed {processed}, newly mapped {changed}, cache={len(cache)}")
    write_cache(cache)
    print(f"Processed {processed}, newly mapped {changed}, cache={len(cache)}")
    return queue_rows


def fill_entries(entries: list[dict[str, str]], substrate_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_key = {(row["species"], row["substrate_id"]): row for row in substrate_rows}
    output = []
    for row in entries:
        substrate = by_key.get((row["species"], row["substrate_id"]), {})
        smiles = substrate.get("substrate_smiles", "")
        row["substrate_smiles"] = smiles
        row["smiles_status"] = "smiles_mapped" if smiles else "needs_smiles_mapping"
        row["smiles_source"] = substrate.get("smiles_source", "")
        row["smiles_source_id"] = substrate.get("smiles_source_id", "")
        output.append(row)
    return output


def refresh_reports(queue_rows: list[dict[str, str]], entries: list[dict[str, str]]) -> None:
    summary = [
        {
            "substrates_total": str(len(queue_rows)),
            "substrates_with_smiles": str(sum(bool(row.get("substrate_smiles")) for row in queue_rows)),
            "entries_total": str(len(entries)),
            "entries_with_smiles": str(sum(bool(row.get("substrate_smiles")) for row in entries)),
            "catpred_ready_rows": str(
                sum(
                    bool(row.get("substrate_smiles"))
                    and bool(row.get("protein_sequence"))
                    and row.get("enzyme_complex_type") == "single_gene"
                    for row in entries
                )
            ),
        }
    ]
    write_rows(
        TABLES / "pubchem_smiles_summary.csv",
        summary,
        ["substrates_total", "substrates_with_smiles", "entries_total", "entries_with_smiles", "catpred_ready_rows"],
    )

    by_species: dict[str, dict[str, int]] = {}
    for row in entries:
        stats = by_species.setdefault(row["species"], {"entries": 0, "entries_with_smiles": 0, "catpred_ready_rows": 0})
        stats["entries"] += 1
        if row.get("substrate_smiles"):
            stats["entries_with_smiles"] += 1
            if row.get("protein_sequence") and row.get("enzyme_complex_type") == "single_gene":
                stats["catpred_ready_rows"] += 1
    write_rows(
        TABLES / "smiles_coverage_by_species.csv",
        [{"species": s, **{k: str(v) for k, v in stats.items()}} for s, stats in sorted(by_species.items())],
        ["species", "entries", "entries_with_smiles", "catpred_ready_rows"],
    )

    coverage_path = TABLES / "experimental_kcat_coverage.csv"
    if coverage_path.exists():
        coverage_rows = read_rows(coverage_path)
        counts = {s: str(stats["entries_with_smiles"]) for s, stats in by_species.items()}
        for row in coverage_rows:
            if row["stage"] == "entries_with_smiles":
                row["n"] = counts.get(row["species"], "0")
        write_rows(coverage_path, coverage_rows, ["species", "stage", "stage_label", "n"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="Only process first N missing substrates.")
    parser.add_argument("--sleep", type=float, default=0.15)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--save-every", type=int, default=50)
    parser.add_argument("--only-non-cofactor", action="store_true")
    parser.add_argument("--skip-name", action="store_true", help="Use ChEBI/KEGG RegistryID only; faster and more conservative.")
    parser.add_argument("--max-chebi", type=int, default=5, help="Maximum ChEBI registry IDs to try per substrate.")
    parser.add_argument("--species", choices=["ecoli", "yeast"], help="Only query missing substrates for one species.")
    parser.add_argument("--ckb-db", default=str(CKB_DB), help="Local CKB compounds.sqlite path for offline SMILES lookup.")
    parser.add_argument("--ckb-only", action="store_true", help="Use local CKB only; do not query PubChem.")
    parser.add_argument("--skip-ckb", action="store_true", help="Skip local CKB lookup before PubChem.")
    parser.add_argument(
        "--only-pubchem-candidate",
        action="store_true",
        help="Skip proteins, nucleic acids, polymers, metal clusters, and clearly ambiguous substrates during PubChem queries.",
    )
    parser.add_argument("--update-only", action="store_true", help="Do not query PubChem; only synchronize entries/reports from substrate queue.")
    args = parser.parse_args()

    queue_rows = read_rows(SUBSTRATE_QUEUE)
    for row in queue_rows:
        row.setdefault("substrate_smiles", "")
        row.setdefault("smiles_source", "")
        row.setdefault("smiles_source_id", "")
    cache = read_cache()
    propagated = propagate_existing_mappings(queue_rows)
    if propagated:
        print(f"Propagated {propagated} existing SMILES mappings across matching substrate identifiers.")
    if not args.update_only:
        if not args.skip_ckb:
            ckb_changed = fill_queue_from_ckb(queue_rows, args)
            print(f"Filled {ckb_changed} missing SMILES from local CKB.")
            if ckb_changed:
                propagated = propagate_existing_mappings(queue_rows)
                if propagated:
                    print(f"Propagated {propagated} CKB SMILES mappings across matching substrate identifiers.")
        if not args.ckb_only:
            queue_rows = fill_queue(queue_rows, cache, args)
        write_rows(SUBSTRATE_QUEUE, queue_rows, list(queue_rows[0].keys()))
    elif propagated:
        write_rows(SUBSTRATE_QUEUE, queue_rows, list(queue_rows[0].keys()))

    entry_source = SMILES_ENTRIES if SMILES_ENTRIES.exists() else (SEQUENCE_ENTRIES if SEQUENCE_ENTRIES.exists() else BASE_ENTRIES)
    entries = fill_entries(read_rows(entry_source), queue_rows)
    write_rows(SMILES_ENTRIES, entries, list(entries[0].keys()))
    refresh_reports(queue_rows, entries)
    print(f"Wrote updated entries to {SMILES_ENTRIES}")


if __name__ == "__main__":
    main()
