#!/usr/bin/env python3
"""Classify model reaction participants using auditable chemical evidence.

The classifier deliberately separates participant-role annotation from
experimental-substrate matching. A currency/cofactor-like flag can support
stratified sensitivity analyses, but it must not override a substrate-specific
experimental record.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable


BASE = Path(__file__).resolve().parent.parent
DEFAULT_REGISTRY = BASE / "configs" / "currency_cofactor_registry.csv"


def split_values(value: object) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def normalize_name(value: object) -> str:
    text = str(value or "").strip().lower()
    text = text.replace("−", "-").replace("_", " ")
    text = re.sub(r"[\[\](){}]", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def normalize_identifier(value: object) -> str:
    text = str(value or "").strip().lower()
    text = text.replace("https://identifiers.org/", "")
    return text


def base_metabolite_id(value: object) -> str:
    return re.sub(r"_[a-z][a-z0-9]*$", "", str(value or "").strip().lower())


@dataclass(frozen=True)
class RoleEvidence:
    is_currency_or_cofactor_like: bool
    role_class: str
    evidence: tuple[str, ...]
    registry_name: str = ""


@lru_cache(maxsize=4)
def load_registry(path: str = str(DEFAULT_REGISTRY)) -> tuple[dict[str, str], ...]:
    registry_path = Path(path)
    if not registry_path.exists():
        raise FileNotFoundError(f"Currency/cofactor registry not found: {registry_path}")
    with registry_path.open(newline="", encoding="utf-8") as handle:
        return tuple(dict(row) for row in csv.DictReader(handle))


def _field_values(row: dict[str, str], column: str, *, names: bool = False) -> set[str]:
    normalizer = normalize_name if names else normalize_identifier
    return {normalizer(value) for value in split_values(row.get(column, "")) if normalizer(value)}


def _name_matches(observed: str, aliases: set[str]) -> bool:
    if not observed:
        return False
    for alias in aliases:
        if observed == alias:
            return True
        # BiGG JSON names sometimes append a molecular formula, for example
        # "ATP C10H12N5O13P3" and "H2O H2O".
        if observed.startswith(alias + " ") and re.match(r"^[a-z0-9+\-]+$", observed[len(alias) + 1 :]):
            return True
    return False


def classify_participant(
    *,
    metabolite_id: object = "",
    name: object = "",
    bigg_ids: object = "",
    kegg_ids: object = "",
    chebi_ids: object = "",
    metanetx_ids: object = "",
    pubchem_cids: object = "",
    inchikey_connectivity: object = "",
    registry_path: Path = DEFAULT_REGISTRY,
) -> RoleEvidence:
    """Return a role flag plus every registry field that supported it."""

    observed_name = normalize_name(name)
    observed = {
        "bigg_ids": {base_metabolite_id(metabolite_id)}
        | {normalize_identifier(value) for value in split_values(bigg_ids)},
        "kegg_ids": {normalize_identifier(value) for value in split_values(kegg_ids)},
        "chebi_ids": {normalize_identifier(value) for value in split_values(chebi_ids)},
        "metanetx_ids": {normalize_identifier(value) for value in split_values(metanetx_ids)},
        "pubchem_cids": {normalize_identifier(value).removeprefix("cid:") for value in split_values(pubchem_cids)},
        "inchikey_connectivity": {
            normalize_identifier(value).split("-", 1)[0] for value in split_values(inchikey_connectivity)
        },
    }
    observed = {key: {value for value in values if value} for key, values in observed.items()}

    candidates: list[tuple[int, int, dict[str, str], tuple[str, ...]]] = []
    for registry_index, row in enumerate(load_registry(str(registry_path))):
        evidence: list[str] = []
        aliases = _field_values(row, "aliases", names=True) | {normalize_name(row.get("canonical_name", ""))}
        if _name_matches(observed_name, aliases):
            evidence.append("normalized_name")
        for column in observed:
            registered = _field_values(row, column)
            if column == "pubchem_cids":
                registered = {value.removeprefix("cid:") for value in registered}
            if column == "inchikey_connectivity":
                registered = {value.split("-", 1)[0] for value in registered}
            overlap = observed[column] & registered
            if overlap:
                evidence.append(f"{column}:{'|'.join(sorted(overlap))}")
        if evidence:
            evidence_types = {item.split(":", 1)[0] for item in evidence}
            identifier_types = evidence_types - {"normalized_name", "inchikey_connectivity"}
            has_name = "normalized_name" in evidence_types
            has_structure = "inchikey_connectivity" in evidence_types
            # A single stale cross-reference must not turn an unrelated
            # metabolite into a cofactor. Without a name match, require either
            # two identifier namespaces or an identifier plus structure.
            if not has_name and len(identifier_types) < 2 and not (identifier_types and has_structure):
                continue
            score = sum(
                100 if item == "normalized_name" else 5 if item.startswith("inchikey_connectivity:") else 50
                for item in evidence
            )
            candidates.append((score, -registry_index, row, tuple(evidence)))

    if candidates:
        _, _, row, evidence = max(candidates, key=lambda item: (item[0], item[1]))
        return RoleEvidence(
            is_currency_or_cofactor_like=True,
            role_class=row.get("role_class", "currency_or_cofactor_like"),
            evidence=evidence,
            registry_name=row.get("canonical_name", ""),
        )

    # CoA-linked metabolites are kept as their own broad class. They can be
    # primary variable substrates, so this flag is descriptive rather than an
    # instruction to remove them from experimental matching.
    if re.search(r"(?:^|[-\s])coa(?:$|\s)", observed_name) or observed_name.endswith("coa"):
        return RoleEvidence(True, "carrier_linked_metabolite", ("normalized_name:coa_linked",))

    return RoleEvidence(False, "other_reactant", tuple())


def choose_reactants_for_matching(participant_ids: Iterable[str]) -> list[str]:
    """Keep every model reactant; experimental evidence selects the substrate."""

    return [str(identifier) for identifier in participant_ids if str(identifier)]
