#!/usr/bin/env python3
"""Create final benchmark-ready tables and unresolved SMILES review lists."""

from __future__ import annotations

import csv
import re
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
INTERIM = BASE / "data" / "interim"
FINAL = BASE / "data" / "final"
TABLES = BASE / "reports" / "tables"

SUBSTRATE_QUEUE = INTERIM / "substrate_smiles_queue.csv"
ENTRIES = INTERIM / "enzyme_reaction_entries_with_sequence_smiles.csv"
CATPRED_INPUT = INTERIM / "prediction_inputs" / "catpred_kcat_input.csv"
TRUTH = FINAL / "experimental_kcat_truth.csv"

UNRESOLVED = INTERIM / "unresolved_smiles_review.csv"
READY_TRUTH = FINAL / "benchmark_ready_truth.csv"
READY_CATPRED = FINAL / "benchmark_ready_catpred.csv"


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


def hard_smiles_reason(name: str, bigg_id: str = "") -> str:
    text = f"{name} {bigg_id}"
    checks = [
        ("protein_or_redox_carrier", r"\bACP\b|acyl-carrier|\[acyl-carrier protein\]|apoprotein|flavodoxin|glutaredoxin|ferredoxin|cytochrome|thioredoxin|IscS|SufBCD|SufSE|disulfide isomerase|scaffold complex"),
        ("nucleic_acid_polymer", r"\bDNA\b|\bRNA\b|tRNA|rRNA|mRNA|oligoribonucleotide|polynucleotide"),
        ("polymer_or_glycan", r"glucan|cellulose|chitin|starch|dextrin|mannan|glycogen|poly|oligosaccharide|peptidoglycan|lipopolysaccharide|murein|^G\d{5}$"),
        ("metal_cluster_or_cofactor_complex", r"\[\dFe-|iron-sulfur|molybdenum cofactor|molybdopterin|tungsten bispterin|damaged iron"),
        ("ambiguous_or_variable_structure", r"\?|unknown|unspecified|electron acceptor|electron donor|radical"),
        ("curated_lipid_needed", r"phosphatidyl|cardiolipin|diacylglycerol|triglyceride|ceramide|sphing|CDP-diacylglycerol|lysocardiolipin"),
    ]
    for reason, pattern in checks:
        if re.search(pattern, text, flags=re.IGNORECASE):
            return reason
    return "no_database_mapping"


def action_required(reason: str) -> str:
    return {
        "protein_or_redox_carrier": "Do not use as a small-molecule SMILES benchmark substrate unless a method supports protein/redox-carrier substrates.",
        "nucleic_acid_polymer": "Provide a curated RNA/DNA structure representation or exclude from SMILES-based predictors.",
        "polymer_or_glycan": "Provide exact oligomer/glycan SMILES/InChI or exclude from small-molecule SMILES predictors.",
        "metal_cluster_or_cofactor_complex": "Provide a curated structure for the exact metal/cofactor complex or exclude.",
        "ambiguous_or_variable_structure": "Provide exact stereochemistry/structure and database cross-reference.",
        "curated_lipid_needed": "Provide a lipid structure library mapping this model lipid ID to exact SMILES/InChI.",
        "no_database_mapping": "Provide SMILES/InChI or a trusted KEGG/ChEBI/MetaNetX/BiGG cross-reference.",
    }[reason]


def bool_field(value: str) -> bool:
    return bool(str(value or "").strip()) and str(value).strip().lower() != "nan"


def main() -> None:
    substrates = read_rows(SUBSTRATE_QUEUE)
    entries = read_rows(ENTRIES)
    truth = read_rows(TRUTH)
    catpred = read_rows(CATPRED_INPUT)

    entries_by_id = {row["entry_id"]: row for row in entries}
    truth_by_entry = {}
    for row in truth:
        truth_by_entry.setdefault(row["entry_id"], []).append(row)
    catpred_by_id = {row["entry_id"]: row for row in catpred}
    catpred_ready_ids = set(catpred_by_id)

    unresolved_rows = []
    for substrate in substrates:
        if bool_field(substrate.get("substrate_smiles", "")):
            continue
        reason = hard_smiles_reason(substrate.get("substrate_name", ""), substrate.get("substrate_bigg_id", ""))
        related_entries = [
            row
            for row in entries
            if row.get("species") == substrate.get("species") and row.get("substrate_id") == substrate.get("substrate_id")
        ]
        blocked_truth = sum(len(truth_by_entry.get(row["entry_id"], [])) for row in related_entries)
        blocked_catpred_candidates = sum(
            bool_field(row.get("protein_sequence", ""))
            and row.get("enzyme_complex_type") == "single_gene"
            for row in related_entries
        )
        unresolved_rows.append(
            {
                **substrate,
                "unresolved_reason": reason,
                "action_required": action_required(reason),
                "blocked_truth_rows": str(blocked_truth),
                "blocked_catpred_candidate_entries": str(blocked_catpred_candidates),
            }
        )

    ready_truth = [row for row in truth if row["entry_id"] in catpred_ready_ids]
    ready_catpred = []
    for row in ready_truth:
        cat = catpred_by_id[row["entry_id"]]
        entry = entries_by_id.get(row["entry_id"], {})
        output_row = dict(row)
        output_row["SMILES"] = cat["SMILES"]
        output_row["sequence"] = cat["sequence"]
        output_row["pdbpath"] = cat["pdbpath"]
        output_row["enzyme_complex_type"] = row.get(
            "enzyme_complex_type", entry.get("enzyme_complex_type", "")
        )
        ready_catpred.append(output_row)

    write_rows(UNRESOLVED, unresolved_rows, list(unresolved_rows[0].keys()) if unresolved_rows else [])
    write_rows(READY_TRUTH, ready_truth, list(truth[0].keys()) if truth else [])
    core_fields = [
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "sequence",
        "pdbpath",
        "true_kcat",
        "true_kcat_log10",
        "unit",
        "pH",
        "temperature_c",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
        "enzyme_complex_type",
    ]
    provenance_fields = [field for field in truth[0].keys() if field not in core_fields] if truth else []
    write_rows(READY_CATPRED, ready_catpred, core_fields + provenance_fields)

    species_values = sorted({row["species"] for row in entries} | {row["species"] for row in truth})
    readiness_rows = []
    for species in species_values:
        species_entries = [row for row in entries if row["species"] == species]
        species_truth = [row for row in truth if row["species"] == species]
        species_ready_truth = [row for row in ready_truth if row["species"] == species]
        readiness_rows.append(
            {
                "species": species,
                "entries": str(len(species_entries)),
                "entries_with_sequence": str(sum(bool_field(row.get("protein_sequence", "")) for row in species_entries)),
                "entries_with_smiles": str(sum(bool_field(row.get("substrate_smiles", "")) for row in species_entries)),
                "catpred_ready_rows": str(sum(row["entry_id"] in catpred_ready_ids for row in species_entries)),
                "truth_rows": str(len(species_truth)),
                "truth_with_smiles": str(sum(bool_field(row.get("substrate_smiles", "")) for row in species_truth)),
                "truth_catpred_ready_rows": str(len(species_ready_truth)),
                "unresolved_substrates": str(sum(row["species"] == species for row in unresolved_rows)),
            }
        )
    write_rows(
        TABLES / "final_benchmark_readiness.csv",
        readiness_rows,
        [
            "species",
            "entries",
            "entries_with_sequence",
            "entries_with_smiles",
            "catpred_ready_rows",
            "truth_rows",
            "truth_with_smiles",
            "truth_catpred_ready_rows",
            "unresolved_substrates",
        ],
    )

    summary = {}
    for row in unresolved_rows:
        key = (row["species"], row["unresolved_reason"])
        stats = summary.setdefault(key, {"substrates": 0, "entries": 0, "blocked_truth_rows": 0, "blocked_catpred_candidate_entries": 0})
        stats["substrates"] += 1
        stats["entries"] += int(row.get("n_entries", "0") or 0)
        stats["blocked_truth_rows"] += int(row.get("blocked_truth_rows", "0") or 0)
        stats["blocked_catpred_candidate_entries"] += int(row.get("blocked_catpred_candidate_entries", "0") or 0)
    summary_rows = [
        {"species": species, "unresolved_reason": reason, **{key: str(value) for key, value in stats.items()}}
        for (species, reason), stats in sorted(summary.items())
    ]
    write_rows(
        TABLES / "unresolved_smiles_summary.csv",
        summary_rows,
        ["species", "unresolved_reason", "substrates", "entries", "blocked_truth_rows", "blocked_catpred_candidate_entries"],
    )

    print(f"Wrote {len(ready_catpred)} CatPred/CataPro-ready truth rows to {READY_CATPRED}")
    print(f"Wrote {len(unresolved_rows)} unresolved substrate rows to {UNRESOLVED}")


if __name__ == "__main__":
    main()
