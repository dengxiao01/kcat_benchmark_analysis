#!/usr/bin/env python3
"""Finalize substrate-role annotations after chemical structures are available.

Model identifiers and names are useful early in the pipeline, but yeast-GEM
uses opaque ``s_XXXX`` identifiers. This pass therefore recomputes each role
from normalized names, external identifiers, PubChem CID, and a standardized
largest-fragment InChIKey. Stoichiometry remains record-level audit context; it
is not treated as sufficient evidence that a participant is a currency
metabolite.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from substrate_roles import DEFAULT_REGISTRY, classify_participant, load_registry, split_values


BASE = Path(__file__).resolve().parent.parent
QUEUE = BASE / "data" / "interim" / "substrate_smiles_queue.csv"
AUDIT = BASE / "reports" / "tables" / "substrate_role_registry_audit.csv"
SUMMARY = BASE / "reports" / "tables" / "substrate_role_by_species.csv"

RDLogger.DisableLog("rdApp.*")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Recompute substrate roles using identifiers and structures.")
    parser.add_argument("--queue", type=Path, default=QUEUE)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_pubchem_cids(source_id: object) -> str:
    values = re.findall(r"(?:^|[|;])CID:(\d+)", str(source_id or ""), flags=re.IGNORECASE)
    return ";".join(dict.fromkeys(values))


def standardized_parent_keys(smiles: object) -> tuple[str, str, str]:
    text = str(smiles or "").strip()
    if not text or text.lower() == "nan":
        return "", "", "missing_smiles"
    try:
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            return "", "", "invalid_smiles"
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        key = Chem.MolToInchiKey(mol)
        if not key:
            return "", "", "inchikey_unavailable"
        return key, key.split("-", 1)[0], "ok"
    except Exception as exc:  # noqa: BLE001 - retain row-level audit instead of aborting the batch.
        return "", "", f"error:{type(exc).__name__}"


def evidence_types(evidence: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(item.split(":", 1)[0] for item in evidence))


def confidence_label(evidence: tuple[str, ...]) -> str:
    kinds = set(evidence_types(evidence))
    if not kinds:
        return "not_registry_matched"
    if len(kinds) >= 2:
        return "multiple_evidence_types"
    if kinds == {"normalized_name"}:
        return "name_only"
    return "identifier_or_structure_only"


def broad_role_group(role_class: str, is_role_like: bool) -> str:
    if role_class == "carrier_linked_metabolite":
        return "carrier_linked_variable"
    if is_role_like:
        return "currency_or_cofactor"
    return "other_reactant"


def registry_structure_consistency(
    registry_by_name: dict[str, dict[str, str]], registry_name: str, connectivity_key: str
) -> str:
    registered = {
        value.split("-", 1)[0].lower()
        for value in split_values(registry_by_name.get(registry_name, {}).get("inchikey_connectivity", ""))
    }
    if not registered or not connectivity_key:
        return "not_assessable"
    return "consistent" if connectivity_key.lower() in registered else "mismatch"


def main() -> None:
    args = parse_args()
    rows = read_rows(args.queue)
    if not rows:
        raise RuntimeError(f"No substrate rows found in {args.queue}")

    registry_rows = list(load_registry(str(args.registry)))
    if any(None in row for row in registry_rows):
        raise ValueError("Currency/cofactor registry contains an unquoted delimiter or extra column.")
    registry_by_name = {row.get("canonical_name", ""): row for row in registry_rows}

    audit_rows: list[dict[str, object]] = []
    changed = 0
    for row in rows:
        prior_class = row.get("substrate_role_class", "")
        pubchem_cids = parse_pubchem_cids(row.get("smiles_source_id", ""))
        parent_key, connectivity_key, structure_status = standardized_parent_keys(row.get("substrate_smiles", ""))
        role = classify_participant(
            metabolite_id=row.get("substrate_id", ""),
            name=row.get("substrate_name", ""),
            bigg_ids=row.get("substrate_bigg_id", ""),
            kegg_ids=row.get("substrate_kegg_id", ""),
            chebi_ids=row.get("substrate_chebi_id", ""),
            metanetx_ids=row.get("substrate_metanetx_id", ""),
            pubchem_cids=pubchem_cids,
            inchikey_connectivity=connectivity_key,
            registry_path=args.registry,
        )
        role_group = broad_role_group(role.role_class, role.is_currency_or_cofactor_like)
        confidence = confidence_label(role.evidence)
        kinds = evidence_types(role.evidence)
        consistency = registry_structure_consistency(
            registry_by_name, role.registry_name, connectivity_key
        )
        if consistency == "mismatch":
            confidence = f"{confidence}_structure_conflict"

        row["substrate_pubchem_cid"] = pubchem_cids
        row["substrate_parent_inchikey"] = parent_key
        row["substrate_parent_inchikey_connectivity"] = connectivity_key
        row["substrate_structure_standardization_status"] = structure_status
        row["substrate_is_cofactor_like"] = str(role.is_currency_or_cofactor_like)
        row["substrate_role_class"] = role.role_class
        row["substrate_role_group"] = role_group
        row["substrate_role_evidence"] = ";".join(role.evidence)
        row["substrate_role_evidence_types"] = ";".join(kinds)
        row["substrate_role_evidence_count"] = str(len(kinds))
        row["substrate_role_confidence"] = confidence
        row["substrate_role_registry_name"] = role.registry_name
        row["substrate_role_registry_structure_consistency"] = consistency
        changed += int(prior_class != role.role_class)

        audit_rows.append(
            {
                "species": row.get("species", ""),
                "substrate_id": row.get("substrate_id", ""),
                "substrate_name": row.get("substrate_name", ""),
                "substrate_bigg_id": row.get("substrate_bigg_id", ""),
                "substrate_kegg_id": row.get("substrate_kegg_id", ""),
                "substrate_chebi_id": row.get("substrate_chebi_id", ""),
                "substrate_metanetx_id": row.get("substrate_metanetx_id", ""),
                "substrate_pubchem_cid": pubchem_cids,
                "substrate_parent_inchikey_connectivity": connectivity_key,
                "structure_standardization_status": structure_status,
                "prior_role_class": prior_class,
                "substrate_role_class": role.role_class,
                "substrate_role_group": role_group,
                "substrate_role_registry_name": role.registry_name,
                "substrate_role_evidence": ";".join(role.evidence),
                "substrate_role_evidence_types": ";".join(kinds),
                "substrate_role_confidence": confidence,
                "registry_structure_consistency": consistency,
                "role_class_changed": str(prior_class != role.role_class),
                "n_reactions": row.get("n_reactions", ""),
                "n_entries": row.get("n_entries", ""),
            }
        )

    queue_fields = list(rows[0].keys())
    write_rows(args.queue, rows, queue_fields)
    write_rows(args.audit, audit_rows, list(audit_rows[0].keys()))

    grouped: dict[tuple[str, str, str, str], dict[str, int]] = {}
    for row in rows:
        key = (
            row.get("species", ""),
            row.get("substrate_role_group", ""),
            row.get("substrate_role_class", ""),
            row.get("substrate_role_confidence", ""),
        )
        stats = grouped.setdefault(key, {"substrates": 0, "substrates_with_smiles": 0, "entries": 0})
        stats["substrates"] += 1
        stats["substrates_with_smiles"] += int(bool(row.get("substrate_smiles", "")))
        stats["entries"] += int(float(row.get("n_entries", "0") or 0))
    summary_rows = [
        {
            "species": key[0],
            "substrate_role_group": key[1],
            "substrate_role_class": key[2],
            "substrate_role_confidence": key[3],
            **stats,
        }
        for key, stats in sorted(grouped.items())
    ]
    write_rows(args.summary, summary_rows, list(summary_rows[0].keys()))

    mismatches = sum(row["registry_structure_consistency"] == "mismatch" for row in audit_rows)
    print(f"Finalized roles for {len(rows)} substrate rows; changed class for {changed} rows.")
    print(f"Registry/structure mismatches requiring review: {mismatches}")
    print(f"Wrote queue: {args.queue}")
    print(f"Wrote audit: {args.audit}")


if __name__ == "__main__":
    main()
