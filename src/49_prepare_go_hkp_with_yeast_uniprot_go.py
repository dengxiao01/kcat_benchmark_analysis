#!/usr/bin/env python3
"""Prepare GO-HKP predictions, adding yeast via UniProt GO annotations."""

from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_ECOLI_ASSIGNMENT = (
    BASE
    / "external_methods"
    / "GO-HKP"
    / "analysis"
    / "DeepGO-SE"
    / "iML1515R"
    / "go_kcat_mean_parent_process_Total_median.json"
)
DEFAULT_GO_KCAT = BASE / "external_methods" / "GO-HKP" / "data" / "GO" / "GO_kcat_tree_total.csv"
DEFAULT_OBO = BASE / "external_methods" / "GO-HKP" / "data" / "GO" / "go-basic.obo"
DEFAULT_YEAST_GO_TSV = BASE / "data" / "raw" / "go_hkp" / "yeast_uniprot_go.tsv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "go_hkp"
DEFAULT_REPORT = BASE / "reports" / "tables" / "go_hkp_eval_readiness.csv"

REQUIRED_COLUMNS = {
    "entry_id",
    "species",
    "reaction_id",
    "gene_id",
    "uniprot_id",
    "ec_number",
    "substrate_name",
    "SMILES",
    "sequence",
    "true_kcat",
    "true_kcat_log10",
    "unit",
    "source_database",
    "match_level",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create GO-HKP benchmark files. E. coli uses the local GO-HKP DeepGO-SE "
            "reaction assignment; yeast uses UniProt GO annotations as a GOATOOLS-style supplement."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--ecoli-assignment", type=Path, default=DEFAULT_ECOLI_ASSIGNMENT)
    parser.add_argument("--go-kcat", type=Path, default=DEFAULT_GO_KCAT)
    parser.add_argument("--obo", type=Path, default=DEFAULT_OBO)
    parser.add_argument("--yeast-go-tsv", type=Path, default=DEFAULT_YEAST_GO_TSV)
    parser.add_argument("--download-uniprot-go", action="store_true")
    parser.add_argument("--yeast-organism-filter", default="Saccharomyces cerevisiae")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_input(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    if df["entry_id"].duplicated().any():
        examples = ", ".join(df.loc[df["entry_id"].duplicated(), "entry_id"].head(5))
        raise ValueError(f"entry_id must be unique. Examples: {examples}")


def load_ecoli_assignments(path: Path) -> dict[str, float]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    assignments: dict[str, float] = {}
    for reaction_id, value in raw.items():
        try:
            kcat = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(kcat) and kcat > 0:
            assignments[str(reaction_id)] = kcat
    return assignments


def download_uniprot_go(accessions: list[str], out_path: Path, chunk_size: int = 40) -> None:
    if not accessions:
        raise ValueError("No UniProt accessions were provided.")
    lines: list[str] = []
    header = ""
    for start in range(0, len(accessions), chunk_size):
        chunk = accessions[start : start + chunk_size]
        query = "(" + " OR ".join(f"accession:{accession}" for accession in chunk) + ")"
        params = {
            "compressed": "false",
            "format": "tsv",
            "fields": "accession,go_id,xref_geneid",
            "query": query,
        }
        url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(params)
        request = urllib.request.Request(url, headers={"User-Agent": "kcat-benchmark-go-hkp/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            text = response.read().decode("utf-8").strip()
        if not text:
            continue
        parts = text.splitlines()
        if not header:
            header = parts[0]
        lines.extend(parts[1:])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(header + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def split_semicolon_terms(value: object) -> list[str]:
    if pd.isna(value):
        return []
    terms = []
    for item in str(value).replace(",", ";").split(";"):
        item = item.strip()
        if item and item.lower() != "nan":
            terms.append(item)
    return terms


def load_uniprot_go(path: Path) -> dict[str, dict[str, list[str]]]:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} does not exist. Re-run with --download-uniprot-go or provide this TSV."
        )
    df = pd.read_csv(path, sep="\t")
    expected = {"Entry", "Gene Ontology IDs", "GeneID"}
    missing = expected.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    mapping: dict[str, dict[str, list[str]]] = {}
    for _, row in df.iterrows():
        entry = str(row["Entry"]).strip()
        if not entry:
            continue
        mapping[entry] = {
            "go_terms": split_semicolon_terms(row["Gene Ontology IDs"]),
            "geneids": split_semicolon_terms(row["GeneID"]),
        }
    return mapping


def parse_float_or_range(value: str) -> float | None:
    token = str(value).strip()
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        pass
    if "-" not in token:
        return None
    parts = [part.strip() for part in token.split("-", 1)]
    if len(parts) != 2:
        return None
    try:
        return (float(parts[0]) + float(parts[1])) / 2.0
    except ValueError:
        return None


def load_go_kcat(path: Path, organism_filter: str) -> dict[str, list[float]]:
    df = pd.read_csv(path)
    if organism_filter:
        df = df[df["Organism"].fillna("").astype(str).str.contains(organism_filter, regex=False)]
    df = df[df["Kcat"].notna()]
    direct: dict[str, list[float]] = defaultdict(list)
    for _, row in df.iterrows():
        go_term = str(row["GO_term"]).strip()
        if not go_term:
            continue
        for token in str(row["Kcat"]).split(";"):
            value = parse_float_or_range(token)
            if value is not None and np.isfinite(value) and value > 0:
                direct[go_term].append(float(value))
    return dict(direct)


def parse_go_children(obo_path: Path) -> dict[str, set[str]]:
    children: dict[str, set[str]] = defaultdict(set)
    current: str | None = None
    with obo_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "[Term]":
                current = None
                continue
            if line.startswith("id: GO:"):
                current = line.split("id: ", 1)[1].strip()
                children.setdefault(current, set())
                continue
            if current and line.startswith("is_a: GO:"):
                parent = line.split("is_a: ", 1)[1].split()[0].strip()
                children[parent].add(current)
            elif current and line.startswith("relationship: part_of GO:"):
                parent = line.split("relationship: part_of ", 1)[1].split()[0].strip()
                children[parent].add(current)
    return dict(children)


def descendant_terms(go_term: str, children: dict[str, set[str]]) -> set[str]:
    seen: set[str] = set()
    stack = list(children.get(go_term, set()))
    while stack:
        child = stack.pop()
        if child in seen:
            continue
        seen.add(child)
        stack.extend(children.get(child, set()))
    return seen


def build_go_kcat_estimator(
    go_kcat: dict[str, list[float]],
    children: dict[str, set[str]],
) -> tuple[callable, callable]:
    value_cache: dict[str, float | None] = {}
    support_cache: dict[str, int] = {}

    def estimate(go_term: str) -> float | None:
        if go_term in value_cache:
            return value_cache[go_term]
        values = list(go_kcat.get(go_term, []))
        if not values:
            for child in descendant_terms(go_term, children):
                values.extend(go_kcat.get(child, []))
        support_cache[go_term] = len(values)
        if values:
            value = float(np.mean(values))
            value_cache[go_term] = value
            return value
        value_cache[go_term] = None
        return None

    def support(go_term: str) -> int:
        if go_term not in support_cache:
            estimate(go_term)
        return support_cache.get(go_term, 0)

    return estimate, support


def assign_rows(
    df: pd.DataFrame,
    ecoli_assignments: dict[str, float],
    uniprot_go: dict[str, dict[str, list[str]]],
    estimate_go_kcat,
    support_go_kcat,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_meta = df.reset_index(drop=True).copy()
    all_meta["entry_id"] = clean_text(all_meta["entry_id"])
    all_meta["reaction_id"] = clean_text(all_meta["reaction_id"])
    all_meta["uniprot_id"] = clean_text(all_meta["uniprot_id"])
    all_meta["go_hkp_prediction_status"] = "missing"
    all_meta["go_hkp_missing_reason"] = "unsupported_species"
    all_meta["go_hkp_assignment_source"] = ""
    all_meta["go_hkp_assignment_kcat"] = np.nan
    all_meta["go_hkp_go_terms"] = ""
    all_meta["go_hkp_go_terms_with_kcat"] = ""
    all_meta["go_hkp_go_terms_with_kcat_count"] = 0
    all_meta["go_hkp_geneid"] = ""
    all_meta["go_hkp_row_id"] = pd.NA

    ecoli_mask = all_meta["species"].eq("ecoli")
    ecoli_matched = ecoli_mask & all_meta["reaction_id"].isin(ecoli_assignments)
    all_meta.loc[ecoli_mask, "go_hkp_missing_reason"] = "missing_reaction_assignment"
    all_meta.loc[ecoli_matched, "go_hkp_prediction_status"] = "ready"
    all_meta.loc[ecoli_matched, "go_hkp_missing_reason"] = ""
    all_meta.loc[ecoli_matched, "go_hkp_assignment_source"] = (
        "GO-HKP DeepGO-SE iML1515R reaction Total median"
    )
    all_meta.loc[ecoli_matched, "go_hkp_assignment_kcat"] = all_meta.loc[
        ecoli_matched, "reaction_id"
    ].map(ecoli_assignments)

    yeast_mask = all_meta["species"].eq("yeast")
    all_meta.loc[yeast_mask, "go_hkp_missing_reason"] = "missing_uniprot_go_terms"
    for idx, row in all_meta[yeast_mask].iterrows():
        accession = str(row["uniprot_id"]).strip()
        info = uniprot_go.get(accession, {})
        go_terms = sorted(set(info.get("go_terms", [])))
        geneids = sorted(set(info.get("geneids", [])))
        all_meta.at[idx, "go_hkp_go_terms"] = ";".join(go_terms)
        all_meta.at[idx, "go_hkp_geneid"] = ";".join(geneids)
        if not go_terms:
            continue
        kcat_by_go = []
        supported_terms = []
        for go_term in go_terms:
            value = estimate_go_kcat(go_term)
            if value is not None and np.isfinite(value) and value > 0:
                kcat_by_go.append(value)
                supported_terms.append(f"{go_term}:{support_go_kcat(go_term)}")
        all_meta.at[idx, "go_hkp_go_terms_with_kcat"] = ";".join(supported_terms)
        all_meta.at[idx, "go_hkp_go_terms_with_kcat_count"] = len(supported_terms)
        if not kcat_by_go:
            all_meta.at[idx, "go_hkp_missing_reason"] = "go_terms_without_kcat_support"
            continue
        all_meta.at[idx, "go_hkp_prediction_status"] = "ready"
        all_meta.at[idx, "go_hkp_missing_reason"] = ""
        all_meta.at[idx, "go_hkp_assignment_source"] = (
            "GO-HKP UniProt GO annotation yeast organism-filtered Total median"
        )
        all_meta.at[idx, "go_hkp_assignment_kcat"] = float(np.median(kcat_by_go))

    ready = all_meta[all_meta["go_hkp_prediction_status"].eq("ready")].copy().reset_index(drop=True)
    ready.insert(0, "go_hkp_row_id_ready", range(len(ready)))
    ready["go_hkp_row_id"] = ready["go_hkp_row_id_ready"]
    ready = ready.drop(columns=["go_hkp_row_id_ready"])

    all_meta = all_meta.drop(columns=["go_hkp_row_id"]).merge(
        ready[["entry_id", "go_hkp_row_id"]],
        on="entry_id",
        how="left",
        validate="one_to_one",
    )
    return all_meta, ready


def go_input(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "go_hkp_row_id": df["go_hkp_row_id"],
            "entry_id": df["entry_id"],
            "species": df["species"],
            "reaction_id": df["reaction_id"],
            "gene_id": df["gene_id"],
            "uniprot_id": df["uniprot_id"],
            "assignment_rule": df["go_hkp_assignment_source"],
        }
    )


def go_output(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "go_hkp_row_id": df["go_hkp_row_id"],
            "entry_id": df["entry_id"],
            "reaction_id": df["reaction_id"],
            "prediction_kcat": df["go_hkp_assignment_kcat"].astype(float),
            "prediction_log10": np.log10(df["go_hkp_assignment_kcat"].astype(float)),
            "go_hkp_assignment_source": df["go_hkp_assignment_source"],
        }
    )


def metadata(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "go_hkp_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "sequence",
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
        "go_hkp_prediction_status",
        "go_hkp_missing_reason",
        "go_hkp_assignment_source",
        "go_hkp_assignment_kcat",
        "go_hkp_go_terms_with_kcat_count",
        "go_hkp_go_terms",
        "go_hkp_go_terms_with_kcat",
        "go_hkp_geneid",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def truth(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "go_hkp_row_id",
        "entry_id",
        "species",
        "true_kcat",
        "true_kcat_log10",
        "unit",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
    ]
    return df[[column for column in keep if column in df.columns]].copy()


def readiness(all_meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [("all", all_meta)] + list(all_meta.groupby("species", sort=True))
    for group_name, part in groups:
        rows.append(
            {
                "group": group_name,
                "rows": len(part),
                "ready_rows": int(part["go_hkp_prediction_status"].eq("ready").sum()),
                "missing_rows": int((~part["go_hkp_prediction_status"].eq("ready")).sum()),
                "unique_reactions": part["reaction_id"].nunique(),
                "ready_unique_reactions": part.loc[
                    part["go_hkp_prediction_status"].eq("ready"), "reaction_id"
                ].nunique(),
                "median_go_terms_with_kcat": float(
                    part.loc[
                        part["go_hkp_prediction_status"].eq("ready"),
                        "go_hkp_go_terms_with_kcat_count",
                    ].median()
                )
                if part["go_hkp_prediction_status"].eq("ready").any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_input(df, args.input)
    df["species"] = clean_text(df["species"])
    df["uniprot_id"] = clean_text(df["uniprot_id"])

    yeast_accessions = sorted(set(df.loc[df["species"].eq("yeast"), "uniprot_id"]) - {""})
    if args.download_uniprot_go:
        download_uniprot_go(yeast_accessions, args.yeast_go_tsv)

    ecoli_assignments = load_ecoli_assignments(args.ecoli_assignment)
    uniprot_go = load_uniprot_go(args.yeast_go_tsv)
    go_kcat = load_go_kcat(args.go_kcat, args.yeast_organism_filter)
    children = parse_go_children(args.obo)
    estimate_go_kcat, support_go_kcat = build_go_kcat_estimator(go_kcat, children)
    all_meta, ready = assign_rows(
        df,
        ecoli_assignments,
        uniprot_go,
        estimate_go_kcat,
        support_go_kcat,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    go_input(ready).to_csv(args.out_dir / "go_hkp_kcat_input.csv", index=False)
    go_output(ready).to_csv(args.out_dir / "go_hkp_kcat_input_output.csv", index=False)
    metadata(ready).to_csv(args.out_dir / "go_hkp_kcat_input_metadata.csv", index=False)
    metadata(all_meta).to_csv(args.out_dir / "go_hkp_kcat_all_metadata.csv", index=False)
    truth(ready).to_csv(args.out_dir / "go_hkp_kcat_input_truth.csv", index=False)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    readiness(all_meta).to_csv(args.report, index=False)

    print(f"E. coli GO-HKP assignments loaded: {len(ecoli_assignments)} reactions")
    print(f"Yeast UniProt GO accessions loaded: {len(uniprot_go)}")
    print(f"Yeast organism-filtered GO-kcat terms loaded: {len(go_kcat)}")
    print(f"Ready rows: {len(ready)}")
    print(f"Missing rows: {len(all_meta) - len(ready)}")
    print(f"Wrote GO-HKP files to: {args.out_dir}")


if __name__ == "__main__":
    main()
