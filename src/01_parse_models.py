#!/usr/bin/env python3
"""Build cross-species enzyme-reaction entry tables from local GEM files.

The output is intentionally conservative: sequence, SMILES, and experimental
kcat values are left blank until they are fetched from primary sources.
"""

from __future__ import annotations

import csv
import json
import re
import xml.etree.ElementTree as ET
from itertools import product
from pathlib import Path
from typing import Any, Iterable

from substrate_roles import RoleEvidence, choose_reactants_for_matching, classify_participant


BASE = Path(__file__).resolve().parent.parent
INTERIM = BASE / "data" / "interim"
TABLES = BASE / "reports" / "tables"

RDF_RESOURCE = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"
FBC = "http://www.sbml.org/sbml/level3/version1/fbc/version2"
SBML = "http://www.sbml.org/sbml/level3/version1/core"


def as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value if v not in (None, "")]
    return [str(value)] if str(value) else []


def clean_identifier_from_uri(uri: str, namespace: str) -> str:
    marker = f"/{namespace}/"
    if marker in uri:
        return uri.rsplit(marker, 1)[-1]
    return uri.rsplit("/", 1)[-1]


def annotation_value(annotation: dict[str, Any], key: str) -> str:
    return ";".join(as_list(annotation.get(key)))


def annotation_list(annotation: dict[str, Any], key: str) -> list[str]:
    return as_list(annotation.get(key))


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


class GPRParser:
    """Parse a simple GPR boolean expression into OR-of-AND gene groups."""

    TOKEN_RE = re.compile(r"\(|\)|\band\b|\bor\b|[A-Za-z0-9_.:-]+", re.IGNORECASE)

    def __init__(self, text: str):
        self.tokens = [tok for tok in self.TOKEN_RE.findall(text or "") if tok.strip()]
        self.pos = 0

    def parse(self) -> list[tuple[str, ...]]:
        if not self.tokens:
            return []
        groups = self._parse_or()
        return sorted({tuple(sorted(set(group))) for group in groups if group})

    def _peek(self) -> str | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def _consume(self) -> str | None:
        token = self._peek()
        if token is not None:
            self.pos += 1
        return token

    def _parse_or(self) -> list[tuple[str, ...]]:
        groups = self._parse_and()
        while (self._peek() or "").lower() == "or":
            self._consume()
            groups.extend(self._parse_and())
        return groups

    def _parse_and(self) -> list[tuple[str, ...]]:
        groups = self._parse_factor()
        while (self._peek() or "").lower() == "and":
            self._consume()
            right = self._parse_factor()
            groups = [tuple(list(a) + list(b)) for a, b in product(groups, right)]
        return groups

    def _parse_factor(self) -> list[tuple[str, ...]]:
        token = self._consume()
        if token is None:
            return []
        if token == "(":
            groups = self._parse_or()
            if self._peek() == ")":
                self._consume()
            return groups
        if token == ")":
            return []
        if token.lower() in {"and", "or"}:
            return []
        return [(token,)]


def parse_gpr_groups(gpr: str) -> list[tuple[str, ...]]:
    try:
        return GPRParser(gpr).parse()
    except Exception:
        genes = re.findall(r"[A-Za-z][A-Za-z0-9_.:-]*", gpr or "")
        return [tuple(sorted(set(g for g in genes if g.lower() not in {"and", "or"})))]


def read_ecoli_model(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        model = json.load(handle)

    genes = {g["id"]: g for g in model.get("genes", [])}
    metabolites = {m["id"]: m for m in model.get("metabolites", [])}

    reaction_rows: list[dict[str, str]] = []
    entry_rows: list[dict[str, str]] = []

    for rxn in model.get("reactions", []):
        reaction_id = rxn.get("id", "")
        annotation = rxn.get("annotation", {}) or {}
        gpr = (rxn.get("gene_reaction_rule") or "").strip()
        ec_numbers = annotation_list(annotation, "ec-code")
        reactant_stoichiometry = {
            mid: abs(float(coeff)) for mid, coeff in (rxn.get("metabolites") or {}).items() if coeff < 0
        }
        reactants = list(reactant_stoichiometry)
        products_ = [mid for mid, coeff in (rxn.get("metabolites") or {}).items() if coeff > 0]
        substrates = choose_reactants_for_matching(reactants)
        substrate_selection = "all_model_reactants_for_experimental_matching"
        groups = parse_gpr_groups(gpr)

        reaction_rows.append(
            {
                "species": "ecoli",
                "model_id": "eciML1515",
                "reaction_id": reaction_id,
                "reaction_name": rxn.get("name", ""),
                "reaction_direction": "reverse" if "_reverse" in reaction_id else "forward",
                "is_reversible": str(float(rxn.get("lower_bound", 0)) < 0),
                "gpr": gpr,
                "n_gpr_groups": str(len(groups)),
                "ec_number": ";".join(ec_numbers),
                "bigg_reaction": annotation_value(annotation, "bigg.reaction"),
                "kegg_reaction": annotation_value(annotation, "kegg.reaction"),
                "rhea": annotation_value(annotation, "rhea"),
                "metanetx_reaction": annotation_value(annotation, "metanetx.reaction"),
                "reactant_ids": ";".join(reactants),
                "reactant_stoichiometry": ";".join(
                    f"{reactant_id}:{reactant_stoichiometry[reactant_id]:g}" for reactant_id in reactants
                ),
                "product_ids": ";".join(products_),
            }
        )

        for group_index, gene_group in enumerate(groups, start=1):
            uniprots = []
            gene_names = []
            for gene_id in gene_group:
                gene = genes.get(gene_id, {})
                gene_names.append(gene.get("name", ""))
                uniprots.extend(annotation_list(gene.get("annotation", {}) or {}, "uniprot"))

            for substrate_id in substrates:
                met = metabolites.get(substrate_id, {})
                met_ann = met.get("annotation", {}) or {}
                role = classify_participant(
                    metabolite_id=substrate_id,
                    name=met.get("name", ""),
                    bigg_ids=annotation_value(met_ann, "bigg.metabolite"),
                    kegg_ids=annotation_value(met_ann, "kegg.compound"),
                    chebi_ids=annotation_value(met_ann, "chebi"),
                    metanetx_ids=annotation_value(met_ann, "metanetx.chemical"),
                )
                entry_rows.append(
                    make_entry_row(
                        species="ecoli",
                        organism="Escherichia coli",
                        model_id="eciML1515",
                        reaction_id=reaction_id,
                        reaction_name=rxn.get("name", ""),
                        reaction_direction="reverse" if "_reverse" in reaction_id else "forward",
                        gpr=gpr,
                        group_index=group_index,
                        gene_group=gene_group,
                        gene_names=gene_names,
                        uniprots=uniprots,
                        ec_numbers=ec_numbers,
                        substrate_id=substrate_id,
                        substrate_name=met.get("name", ""),
                        substrate_compartment=met.get("compartment", ""),
                        substrate_stoichiometry=reactant_stoichiometry.get(substrate_id, 1.0),
                        substrate_selection=substrate_selection,
                        substrate_role=role,
                        substrate_annotation=met_ann,
                        bigg_key="bigg.metabolite",
                        kegg_key="kegg.compound",
                        chebi_key="chebi",
                        mnx_key="metanetx.chemical",
                    )
                )

    return reaction_rows, entry_rows


def collect_resources(element: ET.Element | None) -> list[str]:
    if element is None:
        return []
    return [node.attrib[RDF_RESOURCE] for node in element.iter() if RDF_RESOURCE in node.attrib]


def resources_by_namespace(element: ET.Element | None, namespace: str) -> list[str]:
    values = []
    for uri in collect_resources(element):
        if f"/{namespace}/" in uri:
            values.append(clean_identifier_from_uri(uri, namespace))
    return values


def first_child(element: ET.Element | None) -> ET.Element | None:
    if element is None:
        return None
    for child in list(element):
        return child
    return None


def gpa_to_expression(element: ET.Element | None) -> str:
    node = first_child(element)
    return _gpa_node_to_expression(node)


def _gpa_node_to_expression(node: ET.Element | None) -> str:
    if node is None:
        return ""
    name = local_name(node.tag)
    if name == "geneProductRef":
        return node.attrib.get(f"{{{FBC}}}geneProduct", "")
    if name in {"and", "or"}:
        parts = [_gpa_node_to_expression(child) for child in list(node)]
        parts = [part for part in parts if part]
        joiner = f" {name} "
        return "(" + joiner.join(parts) + ")" if len(parts) > 1 else (parts[0] if parts else "")
    return ""


def read_yeast_model(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    root = ET.parse(path).getroot()
    ns = {"sbml": SBML, "fbc": FBC}

    species_map = {}
    for species in root.findall(".//sbml:listOfSpecies/sbml:species", ns):
        species_map[species.attrib.get("id", "")] = species

    gene_uniprot: dict[str, list[str]] = {}
    gene_label: dict[str, str] = {}
    for gene in root.findall(".//fbc:listOfGeneProducts/fbc:geneProduct", ns):
        gene_id = gene.attrib.get(f"{{{FBC}}}id", "")
        gene_label[gene_id] = gene.attrib.get(f"{{{FBC}}}label", "")
        gene_uniprot[gene_id] = resources_by_namespace(gene.find("sbml:annotation", ns), "uniprot")

    reaction_rows: list[dict[str, str]] = []
    entry_rows: list[dict[str, str]] = []

    for rxn in root.findall(".//sbml:listOfReactions/sbml:reaction", ns):
        reaction_id = rxn.attrib.get("id", "")
        reaction_name = rxn.attrib.get("name", "")
        annotation = rxn.find("sbml:annotation", ns)
        gpr = gpa_to_expression(rxn.find("fbc:geneProductAssociation", ns))
        groups = parse_gpr_groups(gpr)
        ec_numbers = resources_by_namespace(annotation, "ec-code")
        reactant_refs = rxn.findall("sbml:listOfReactants/sbml:speciesReference", ns)
        reactants = [ref.attrib.get("species", "") for ref in reactant_refs]
        reactant_stoichiometry = {
            ref.attrib.get("species", ""): abs(float(ref.attrib.get("stoichiometry", "1") or 1))
            for ref in reactant_refs
        }
        products_ = [
            ref.attrib.get("species", "")
            for ref in rxn.findall("sbml:listOfProducts/sbml:speciesReference", ns)
        ]
        substrates = choose_reactants_for_matching(reactants)
        substrate_selection = "all_model_reactants_for_experimental_matching"

        reaction_rows.append(
            {
                "species": "yeast",
                "model_id": "yeast-GEM_v9.0.2",
                "reaction_id": reaction_id,
                "reaction_name": reaction_name,
                "reaction_direction": "forward",
                "is_reversible": rxn.attrib.get("reversible", ""),
                "gpr": gpr,
                "n_gpr_groups": str(len(groups)),
                "ec_number": ";".join(ec_numbers),
                "bigg_reaction": ";".join(resources_by_namespace(annotation, "bigg.reaction")),
                "kegg_reaction": ";".join(resources_by_namespace(annotation, "kegg.reaction")),
                "rhea": ";".join(resources_by_namespace(annotation, "rhea")),
                "metanetx_reaction": ";".join(resources_by_namespace(annotation, "metanetx.reaction")),
                "reactant_ids": ";".join(reactants),
                "reactant_stoichiometry": ";".join(
                    f"{reactant_id}:{reactant_stoichiometry[reactant_id]:g}" for reactant_id in reactants
                ),
                "product_ids": ";".join(products_),
            }
        )

        for group_index, gene_group in enumerate(groups, start=1):
            uniprots = []
            gene_names = []
            for gene_id in gene_group:
                gene_names.append(gene_label.get(gene_id, ""))
                uniprots.extend(gene_uniprot.get(gene_id, []))

            for substrate_id in substrates:
                species_node = species_map.get(substrate_id)
                ann_node = species_node.find("sbml:annotation", ns) if species_node is not None else None
                substrate_annotation = {
                    "bigg.metabolite": resources_by_namespace(ann_node, "bigg.metabolite"),
                    "kegg.compound": resources_by_namespace(ann_node, "kegg.compound"),
                    "chebi": resources_by_namespace(ann_node, "chebi"),
                    "metanetx.chemical": resources_by_namespace(ann_node, "metanetx.chemical"),
                }
                role = classify_participant(
                    metabolite_id=substrate_id,
                    name=species_node.attrib.get("name", "") if species_node is not None else "",
                    bigg_ids=";".join(substrate_annotation["bigg.metabolite"]),
                    kegg_ids=";".join(substrate_annotation["kegg.compound"]),
                    chebi_ids=";".join(substrate_annotation["chebi"]),
                    metanetx_ids=";".join(substrate_annotation["metanetx.chemical"]),
                )
                entry_rows.append(
                    make_entry_row(
                        species="yeast",
                        organism="Saccharomyces cerevisiae",
                        model_id="yeast-GEM_v9.0.2",
                        reaction_id=reaction_id,
                        reaction_name=reaction_name,
                        reaction_direction="forward",
                        gpr=gpr,
                        group_index=group_index,
                        gene_group=gene_group,
                        gene_names=gene_names,
                        uniprots=uniprots,
                        ec_numbers=ec_numbers,
                        substrate_id=substrate_id,
                        substrate_name=species_node.attrib.get("name", "") if species_node is not None else "",
                        substrate_compartment=species_node.attrib.get("compartment", "") if species_node is not None else "",
                        substrate_stoichiometry=reactant_stoichiometry.get(substrate_id, 1.0),
                        substrate_selection=substrate_selection,
                        substrate_role=role,
                        substrate_annotation=substrate_annotation,
                        bigg_key="bigg.metabolite",
                        kegg_key="kegg.compound",
                        chebi_key="chebi",
                        mnx_key="metanetx.chemical",
                    )
                )

    return reaction_rows, entry_rows


def make_entry_row(
    *,
    species: str,
    organism: str,
    model_id: str,
    reaction_id: str,
    reaction_name: str,
    reaction_direction: str,
    gpr: str,
    group_index: int,
    gene_group: Iterable[str],
    gene_names: Iterable[str],
    uniprots: Iterable[str],
    ec_numbers: Iterable[str],
    substrate_id: str,
    substrate_name: str,
    substrate_compartment: str,
    substrate_stoichiometry: float,
    substrate_selection: str,
    substrate_role: RoleEvidence,
    substrate_annotation: dict[str, Any],
    bigg_key: str,
    kegg_key: str,
    chebi_key: str,
    mnx_key: str,
) -> dict[str, str]:
    gene_group = tuple(gene_group)
    uniprots = tuple(sorted(set(u for u in uniprots if u)))
    enzyme_complex_type = "single_gene" if len(gene_group) == 1 else "complex"
    entry_id = f"{species}|{reaction_id}|g{group_index}|{substrate_id}"
    return {
        "entry_id": entry_id,
        "species": species,
        "model_id": model_id,
        "organism": organism,
        "reaction_id": reaction_id,
        "reaction_name": reaction_name,
        "reaction_direction": reaction_direction,
        "gpr": gpr,
        "gpr_group_id": f"{species}|{reaction_id}|g{group_index}",
        "gene_id": ";".join(gene_group),
        "gene_name": ";".join(gene_names),
        "enzyme_complex_type": enzyme_complex_type,
        "uniprot_id": ";".join(uniprots),
        "ec_number": ";".join(sorted(set(ec_numbers))),
        "substrate_id": substrate_id,
        "substrate_name": substrate_name,
        "substrate_compartment": substrate_compartment,
        "substrate_stoichiometry": f"{substrate_stoichiometry:g}",
        "substrate_selection": substrate_selection,
        "candidate_selection_policy": "all_reactants_then_experimental_substrate_matching",
        "substrate_is_cofactor_like": str(substrate_role.is_currency_or_cofactor_like),
        "substrate_role_class": substrate_role.role_class,
        "substrate_role_evidence": ";".join(substrate_role.evidence),
        "substrate_role_registry_name": substrate_role.registry_name,
        "substrate_bigg_id": ";".join(as_list(substrate_annotation.get(bigg_key))),
        "substrate_kegg_id": ";".join(as_list(substrate_annotation.get(kegg_key))),
        "substrate_chebi_id": ";".join(as_list(substrate_annotation.get(chebi_key))),
        "substrate_metanetx_id": ";".join(as_list(substrate_annotation.get(mnx_key))),
        "substrate_smiles": "",
        "reaction_smiles": "",
        "protein_sequence": "",
        "sequence_status": "needs_sequence_fetch" if uniprots else "missing_uniprot",
        "smiles_status": "needs_smiles_mapping",
        "truth_status": "needs_brenda_sabiork_match",
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarize(reactions: list[dict[str, str]], entries: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for species in sorted({r["species"] for r in reactions}):
        rxn_rows = [r for r in reactions if r["species"] == species]
        entry_rows = [e for e in entries if e["species"] == species]
        rows.append(
            {
                "species": species,
                "total_reactions": str(len(rxn_rows)),
                "reactions_with_gpr": str(sum(bool(r["gpr"]) for r in rxn_rows)),
                "reactions_with_ec": str(sum(bool(r["ec_number"]) for r in rxn_rows)),
                "enzyme_substrate_entries": str(len(entry_rows)),
                "entries_with_ec": str(sum(bool(e["ec_number"]) for e in entry_rows)),
                "entries_with_uniprot": str(sum(bool(e["uniprot_id"]) for e in entry_rows)),
                "single_gene_entries": str(sum(e["enzyme_complex_type"] == "single_gene" for e in entry_rows)),
                "complex_entries": str(sum(e["enzyme_complex_type"] == "complex" for e in entry_rows)),
                "entries_with_sequence": "0",
                "entries_with_smiles": "0",
                "experimental_kcat_matched": "0",
            }
        )
    return rows


def make_stage_coverage(summary_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    stages = [
        ("total_reactions", "All model reactions"),
        ("reactions_with_gpr", "Reactions with GPR"),
        ("reactions_with_ec", "Reactions with EC number"),
        ("enzyme_substrate_entries", "Enzyme-substrate entries"),
        ("entries_with_uniprot", "Entries with UniProt"),
        ("entries_with_sequence", "Entries with protein sequence"),
        ("entries_with_smiles", "Entries with substrate SMILES"),
        ("experimental_kcat_matched", "Entries matched to experimental kcat"),
    ]
    rows = []
    for row in summary_rows:
        for key, label in stages:
            rows.append({"species": row["species"], "stage": key, "stage_label": label, "n": row[key]})
    return rows


def main() -> None:
    print("Phase 1: parsing local GEM models")
    ecoli_reactions, ecoli_entries = read_ecoli_model(BASE / "eciML1515.json")
    yeast_reactions, yeast_entries = read_yeast_model(BASE / "yeast-GEM.xml")

    reaction_rows = ecoli_reactions + yeast_reactions
    entry_rows = ecoli_entries + yeast_entries
    summary_rows = summarize(reaction_rows, entry_rows)

    write_csv(INTERIM / "model_reactions.csv", reaction_rows)
    write_csv(INTERIM / "enzyme_reaction_entries.csv", entry_rows)
    write_csv(TABLES / "model_parse_summary.csv", summary_rows)
    write_csv(TABLES / "experimental_kcat_coverage.csv", make_stage_coverage(summary_rows))

    print(f"Wrote {len(reaction_rows)} reactions to {INTERIM / 'model_reactions.csv'}")
    print(f"Wrote {len(entry_rows)} entries to {INTERIM / 'enzyme_reaction_entries.csv'}")
    print(f"Wrote summary to {TABLES / 'model_parse_summary.csv'}")
    for row in summary_rows:
        print(
            "{species}: reactions={total_reactions}, with_gpr={reactions_with_gpr}, "
            "with_ec={reactions_with_ec}, entries={enzyme_substrate_entries}, "
            "entries_with_uniprot={entries_with_uniprot}".format(**row)
        )


if __name__ == "__main__":
    main()
