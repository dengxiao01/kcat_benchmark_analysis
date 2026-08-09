#!/usr/bin/env python3
"""Build submission-review sensitivity, dependence, and overlap audits.

Run this script from an activated environment containing RDKit, SciPy, and pandas:

    python paper/build_submission_audits.py

The analysis deliberately distinguishes a negative result from unavailable
evidence.  Methods without released record-level training corpora are marked
``unknown`` rather than being treated as overlap-free.
"""

from __future__ import annotations

import argparse
import importlib.util
import io
import json
import math
import subprocess
import tarfile
import tempfile
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdChemReactions
from scipy.stats import wilcoxon


BASE = Path(__file__).resolve().parent.parent
TRUTH_PATH = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
CONTEXT_PATH = BASE / "reports" / "tables" / "benchmark_ready_catpred_enriched_context.csv"
ENTRY_PATH = BASE / "data" / "interim" / "enzyme_reaction_entries_with_sequence_smiles.csv"
MODEL_REACTION_PATH = BASE / "data" / "interim" / "model_reactions.csv"
TABLE_DIR = BASE / "paper" / "tables_v1.2.0"
REPORT_TABLE_DIR = BASE / "reports" / "tables"
DETAIL_DIR = BASE / "paper" / "submission_audit_details_v1.2.0"
MATCHER_PATH = BASE / "src" / "10_parse_brenda_kcat.py"

BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_SEED = 20260724
NEAR_SEQUENCE_IDENTITY = 80.0
NEAR_CHEMICAL_TANIMOTO = 0.80

METHOD_FILES = {
    "KcatNet": BASE / "data" / "final" / "kcatnet" / "kcatnet_kcat_predictions_evaluated.csv",
    "CataPro": BASE / "data" / "final" / "catapro" / "catapro_kcat_predictions_evaluated.csv",
    "PreTKcat": BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_predictions_evaluated.csv",
    "UniKP": BASE / "data" / "final" / "unikp" / "unikp_kcat_predictions_evaluated.csv",
    "SELFprot": BASE / "data" / "final" / "selfprot" / "selfprot_kcat_predictions_evaluated.csv",
    "DLKcat": BASE / "data" / "final" / "dlkcat" / "dlkcat_kcat_predictions_evaluated.csv",
    "TurNuP": BASE / "data" / "final" / "turnup" / "turnup_kcat_predictions_evaluated.csv",
    "PMAK": BASE / "data" / "final" / "pmak" / "pmak_kcat_predictions_evaluated.csv",
    "KinForm-L": BASE / "data" / "final" / "kinform" / "kinform_kcat_predictions_evaluated.csv",
    "CatPred": BASE / "data" / "final" / "catpred" / "catpred_kcat_predictions_evaluated.csv",
    "DEKP-public-retrained": (
        BASE
        / "data"
        / "final"
        / "dekp"
        / "dekp_public_retrained_kcat_predictions_evaluated.csv"
    ),
    "GO-HKP": BASE / "data" / "final" / "go_hkp" / "go_hkp_kcat_predictions_evaluated.csv",
}
METHOD_ORDER = list(METHOD_FILES)
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
UNCHARGER = rdMolStandardize.Uncharger()
RDLogger.DisableLog("rdApp.*")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAP_REPLICATES)
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    parser.add_argument(
        "--skip-neighbor-search",
        action="store_true",
        help="Build all audits except DIAMOND-coupled sequence/chemical neighbor classes.",
    )
    parser.add_argument(
        "--skip-brenda-reparse",
        action="store_true",
        help="Do not reparse BRENDA with mutant records included.",
    )
    return parser.parse_args()


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna("").astype(str).str.lower().isin({"true", "1", "yes"})


def load_methods() -> dict[str, pd.DataFrame]:
    methods = {}
    for method, path in METHOD_FILES.items():
        frame = pd.read_csv(path)
        frame["entry_id"] = frame["entry_id"].astype(str)
        methods[method] = frame
    return methods


def clean_sequence(value: object) -> str:
    return "".join(str(value or "").split()).upper()


def truncate_1000(sequence: str) -> str:
    return sequence if len(sequence) <= 1000 else sequence[:500] + sequence[-500:]


def mol_identity(value: object) -> dict[str, object] | None:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return None
    molecule = Chem.MolFromSmiles(text)
    if molecule is None:
        return None
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    try:
        parent = rdMolStandardize.FragmentParent(molecule)
        parent = UNCHARGER.uncharge(parent)
        Chem.SanitizeMol(parent)
    except Exception:
        parent = Chem.Mol(molecule)
        Chem.SanitizeMol(parent)
    parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    try:
        connectivity_key = Chem.MolToInchiKey(parent).split("-", 1)[0]
    except Exception:
        connectivity_key = ""
    if not connectivity_key:
        connectivity_key = parent_smiles
    return {
        "canonical": canonical,
        "parent_smiles": parent_smiles,
        "connectivity_key": connectivity_key,
        "fingerprint": MORGAN.GetFingerprint(parent),
    }


def component_identities(value: object) -> list[dict[str, object]]:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return []
    components = []
    for token in text.split("."):
        identity = mol_identity(token)
        if identity is not None:
            components.append(identity)
    return components


def canonical_reaction(value: object) -> tuple[str, object] | None:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return None
    if ">>" not in text:
        return None
    left, right = text.split(">>", 1)
    left_keys = sorted(
        identity["connectivity_key"] for identity in component_identities(left)
    )
    right_keys = sorted(
        identity["connectivity_key"] for identity in component_identities(right)
    )
    if not left_keys or not right_keys:
        return None
    canonical = ".".join(left_keys) + ">>" + ".".join(right_keys)
    try:
        reaction = rdChemReactions.ReactionFromSmarts(text, useSmiles=True)
        fingerprint = rdChemReactions.CreateDifferenceFingerprintForReaction(reaction)
    except Exception:
        return None
    return canonical, fingerprint


def metrics_from_arrays(truth: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    error = prediction - truth
    abs_error = np.abs(error)
    if len(truth) > 1:
        spearman = float(pd.Series(truth).corr(pd.Series(prediction), method="spearman"))
    else:
        spearman = np.nan
    return {
        "n": int(len(truth)),
        "mae_log10": float(abs_error.mean()),
        "rmse_log10": float(np.sqrt(np.square(error).mean())),
        "spearman_log10": spearman,
        "bias_log10": float(error.mean()),
        "within_0.3_fraction": float((abs_error <= 0.3).mean()),
        "within_1.0_fraction": float((abs_error <= 1.0).mean()),
    }


def frame_metrics(frame: pd.DataFrame) -> dict[str, float]:
    return metrics_from_arrays(
        frame["true_kcat_log10"].to_numpy(float),
        frame["prediction_log10"].to_numpy(float),
    )


def add_pair_keys(truth: pd.DataFrame) -> pd.DataFrame:
    truth = truth.copy()
    identities = [mol_identity(value) for value in truth["SMILES"]]
    if any(identity is None for identity in identities):
        raise ValueError("Final benchmark contains an invalid substrate structure")
    truth["chemical_parent_key"] = [
        str(identity["connectivity_key"]) for identity in identities if identity is not None
    ]
    truth["pair_cluster"] = truth["sequence"].map(clean_sequence) + "|" + truth["chemical_parent_key"]
    truth["protein_cluster"] = truth["sequence"].map(clean_sequence)
    truth["reaction_cluster"] = truth["species"].astype(str) + "|" + truth["reaction_id"].astype(str)
    truth["reference_cluster"] = truth["reference"].fillna("").astype(str)
    empty_ref = truth["reference_cluster"].eq("")
    truth.loc[empty_ref, "reference_cluster"] = "missing_reference|" + truth.loc[
        empty_ref, "entry_id"
    ].astype(str)
    accession_matched = truth["match_level"].eq("species_ec_uniprot_substrate_id")
    truth["label_assignment_cluster"] = "direct|" + truth["entry_id"].astype(str)
    truth.loc[~accession_matched, "label_assignment_cluster"] = (
        "weak|"
        + truth.loc[~accession_matched, "species"].astype(str)
        + "|"
        + truth.loc[~accession_matched, "ec_number"].astype(str)
        + "|"
        + truth.loc[~accession_matched, "chemical_parent_key"].astype(str)
        + "|"
        + truth.loc[~accession_matched, "true_kcat_log10"].astype(str)
        + "|"
        + truth.loc[~accession_matched, "reference"].fillna("").astype(str)
    )
    return truth


def build_label_audit(truth: pd.DataFrame, context: pd.DataFrame) -> pd.DataFrame:
    def row(
        item: str,
        count: int,
        denominator: int,
        interpretation: str,
        status: str = "observed",
    ) -> dict[str, object]:
        return {
            "audit_item": item,
            "status": status,
            "count": int(count),
            "denominator": int(denominator),
            "percent": 100.0 * count / denominator if denominator else np.nan,
            "interpretation": interpretation,
        }

    total = len(truth)
    pair_sizes = truth.groupby("pair_cluster").size()
    reference_sizes = truth.groupby("reference_cluster").size()
    weaker = truth[~truth["match_level"].eq("species_ec_uniprot_substrate_id")].copy()
    assignment = weaker.groupby("label_assignment_cluster").agg(
        rows=("entry_id", "size"),
        sequences=("sequence", "nunique"),
        proteins=("uniprot_id", "nunique"),
    )
    multi = assignment[assignment["sequences"] > 1]
    affected_keys = set(multi.index)
    affected_rows = int(weaker["label_assignment_cluster"].isin(affected_keys).sum())

    merged = truth[
        [
            "entry_id",
            "substrate_role_group",
            "substrate_role_class",
            "substrate_role_confidence",
        ]
    ].merge(
        context[
            [
                "entry_id",
                "currency_or_cofactor_like_by_name",
            ]
        ],
        on="entry_id",
        how="left",
        validate="one_to_one",
    )
    hplus = truth["substrate_name"].fillna("").str.lower().isin({"h+", "proton"})
    water = truth["substrate_name"].fillna("").str.lower().isin({"h2o", "water"})

    rows = [
        row("benchmark_rows", total, total, "Table rows; these are not all independent experiments."),
        row(
            "accession_matched_rows",
            int(truth["match_level"].eq("species_ec_uniprot_substrate_id").sum()),
            total,
            "The experimental record and evaluated sequence share a UniProt accession.",
        ),
        row(
            "species_ec_substrate_id_rows",
            int(truth["match_level"].eq("species_ec_substrate_id").sum()),
            total,
            "The value is linked at species + EC + substrate-ID level, not accession level.",
        ),
        row(
            "species_ec_substrate_name_rows",
            int(truth["match_level"].eq("species_ec_substrate_name").sum()),
            total,
            "The weakest retained matching level uses normalized substrate name.",
        ),
        row(
            "unique_sequence_substrate_pairs",
            int(truth["pair_cluster"].nunique()),
            total,
            "Unique full sequence + standardized substrate-parent identities.",
        ),
        row(
            "rows_in_duplicated_pair_clusters",
            int(truth["pair_cluster"].isin(pair_sizes[pair_sizes > 1].index).sum()),
            total,
            f"Maximum rows sharing one pair: {int(pair_sizes.max())}.",
        ),
        row(
            "unique_reaction_clusters",
            int(truth["reaction_cluster"].nunique()),
            total,
            "Species-specific metabolic-model reaction identifiers.",
        ),
        row(
            "unique_reference_strings",
            int(truth["reference_cluster"].nunique()),
            total,
            "A semicolon-joined reference field is an approximate literature cluster.",
        ),
        row(
            "rows_in_duplicated_reference_clusters",
            int(truth["reference_cluster"].isin(reference_sizes[reference_sizes > 1].index).sum()),
            total,
            f"Maximum rows sharing one reference string: {int(reference_sizes.max())}.",
        ),
        row(
            "unique_label_assignment_clusters",
            int(truth["label_assignment_cluster"].nunique()),
            total,
            (
                "Accession-matched rows are individual units; weaker matches sharing species, EC, "
                "substrate parent, selected label, and reference are one assignment cluster."
            ),
        ),
        row(
            "weaker_evidence_assignment_groups_with_multiple_sequences",
            len(multi),
            len(assignment),
            (
                f"{affected_rows} benchmark rows are affected; the largest group maps one aggregated "
                f"label to {int(assignment['sequences'].max())} sequences."
            ),
        ),
        row(
            "rows_in_shared_multisequence_label_assignments",
            affected_rows,
            total,
            "Rows whose weaker-evidence label-assignment cluster contains multiple sequences.",
        ),
        row(
            "heteromeric_multigene_complex_rows",
            int((truth["enzyme_complex_type"] != "single_gene").sum()),
            total,
            "The final benchmark excludes heteromeric multi-gene entries; oligomeric state is not annotated.",
        ),
        row("monatomic_proton_rows", int(hplus.sum()), total, "Selected model reactant is H+."),
        row("water_rows", int(water.sum()), total, "Selected model reactant is water."),
        row(
            "name_heuristic_currency_or_cofactor_rows",
            int(bool_series(merged["currency_or_cofactor_like_by_name"]).sum()),
            total,
            "Exploratory name-only classification; it is not used for the primary role analysis.",
        ),
        row(
            "registry_currency_or_cofactor_rows",
            int(merged["substrate_role_group"].eq("currency_or_cofactor").sum()),
            total,
            "Joint normalized-name, database-identifier, and standardized-structure registry classification.",
        ),
        row(
            "registry_carrier_linked_variable_rows",
            int(merged["substrate_role_group"].eq("carrier_linked_variable").sum()),
            total,
            "Carrier-linked metabolites are retained as a separate role because they may be the variable substrate.",
        ),
    ]
    return pd.DataFrame(rows)


def build_sensitivity(
    truth: pd.DataFrame,
    context: pd.DataFrame,
    methods: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    context_small = context[["entry_id", "currency_or_cofactor_like_by_name"]].copy()
    joined = truth.merge(context_small, on="entry_id", how="left", validate="one_to_one")
    names = joined["substrate_name"].fillna("").str.lower()
    no_hplus_water = ~names.isin({"h+", "proton", "h2o", "water"})
    primary = joined["substrate_role_group"].eq("other_reactant")
    accession = joined["match_level"].eq("species_ec_uniprot_substrate_id")
    sabiork_participant_ambiguous = joined["experimental_substrate_support"].eq(
        "participant_ambiguous"
    )
    subsets = {
        "all_rows": np.ones(len(joined), dtype=bool),
        "substrate_supported_rows": (~sabiork_participant_ambiguous).to_numpy(),
        "participant_ambiguous_rows": sabiork_participant_ambiguous.to_numpy(),
        "accession_matched_rows": accession.to_numpy(),
        "exclude_hplus_and_water_rows": no_hplus_water.to_numpy(),
        "registry_other_reactant_rows": primary.to_numpy(),
        "registry_currency_or_cofactor_rows": joined["substrate_role_group"].eq("currency_or_cofactor").to_numpy(),
        "registry_carrier_linked_variable_rows": joined["substrate_role_group"].eq("carrier_linked_variable").to_numpy(),
        "accession_registry_other_reactant_rows": (accession & primary).to_numpy(),
        "exclude_sabio_participant_ambiguous": (~sabiork_participant_ambiguous).to_numpy(),
    }
    truth_indexed = joined.set_index("entry_id")
    rows = []
    row_aggregation_rule = "one benchmark row per unit; no within-unit aggregation"
    pair_aggregation_rule = (
        "median observed and predicted log10(kcat) separately within each pair cluster; "
        "MAE averages cluster absolute errors with equal cluster weight"
    )
    label_aggregation_rule = (
        "median observed and predicted log10(kcat) separately within each label-assignment cluster; "
        "MAE averages cluster absolute errors with equal cluster weight"
    )

    for subset_name, mask in subsets.items():
        identifiers = set(joined.loc[mask, "entry_id"].astype(str))
        for method in METHOD_ORDER:
            frame = methods[method]
            part = frame[frame["entry_id"].isin(identifiers)].copy()
            result = frame_metrics(part)
            rows.append(
                {
                    "analysis_scope": subset_name,
                    "statistical_unit": "benchmark_row",
                    "aggregation_rule": row_aggregation_rule,
                    "scope_units": len(identifiers),
                    "method": method,
                    "n": result["n"],
                    "coverage_within_scope_percent": 100.0 * result["n"] / len(identifiers),
                    **{key: value for key, value in result.items() if key != "n"},
                }
            )

    grouped_scopes = {
        "unique_sequence_substrate_pairs": set(joined["entry_id"]),
        "accession_unique_sequence_substrate_pairs": set(
            joined.loc[accession, "entry_id"].astype(str)
        ),
    }
    for subset_name, identifiers in grouped_scopes.items():
        scope_truth = truth_indexed.loc[sorted(identifiers)].reset_index()
        total_groups = scope_truth["pair_cluster"].nunique()
        for method in METHOD_ORDER:
            merged = scope_truth[
                ["entry_id", "pair_cluster", "true_kcat_log10"]
            ].merge(
                methods[method][["entry_id", "prediction_log10"]],
                on="entry_id",
                how="inner",
                validate="one_to_one",
            )
            grouped = (
                merged.groupby("pair_cluster", as_index=False)
                .agg(
                    true_kcat_log10=("true_kcat_log10", "median"),
                    prediction_log10=("prediction_log10", "median"),
                )
            )
            result = frame_metrics(grouped)
            rows.append(
                {
                    "analysis_scope": subset_name,
                    "statistical_unit": "unique_sequence_substrate_pair",
                    "aggregation_rule": pair_aggregation_rule,
                    "scope_units": total_groups,
                    "method": method,
                    "n": result["n"],
                    "coverage_within_scope_percent": 100.0 * result["n"] / total_groups,
                    **{key: value for key, value in result.items() if key != "n"},
                }
            )

    total_label_groups = joined["label_assignment_cluster"].nunique()
    for method in METHOD_ORDER:
        merged = joined[
            ["entry_id", "label_assignment_cluster", "true_kcat_log10"]
        ].merge(
            methods[method][["entry_id", "prediction_log10"]],
            on="entry_id",
            how="inner",
            validate="one_to_one",
        )
        grouped = (
            merged.groupby("label_assignment_cluster", as_index=False)
            .agg(
                true_kcat_log10=("true_kcat_log10", "median"),
                prediction_log10=("prediction_log10", "median"),
            )
        )
        result = frame_metrics(grouped)
        rows.append(
            {
                "analysis_scope": "unique_label_assignments",
                "statistical_unit": "label_assignment_cluster",
                "aggregation_rule": label_aggregation_rule,
                "scope_units": total_label_groups,
                "method": method,
                "n": result["n"],
                "coverage_within_scope_percent": 100.0 * result["n"] / total_label_groups,
                **{key: value for key, value in result.items() if key != "n"},
            }
        )
    reaction_methods = ("KcatNet", "TurNuP", "PMAK")
    reaction_ids = set.intersection(
        *(set(methods[method]["entry_id"].astype(str)) for method in reaction_methods)
    )
    reaction_truth = truth_indexed.loc[sorted(reaction_ids)].reset_index()
    reaction_group_scopes = (
        (
            "reaction_common_unique_sequence_substrate_pairs",
            "pair_cluster",
            "unique_sequence_substrate_pair",
            pair_aggregation_rule,
        ),
        (
            "reaction_common_unique_label_assignments",
            "label_assignment_cluster",
            "label_assignment_cluster",
            label_aggregation_rule,
        ),
    )
    for analysis_scope, cluster_column, statistical_unit, aggregation_rule in reaction_group_scopes:
        total_groups = reaction_truth[cluster_column].nunique()
        for method in reaction_methods:
            merged = reaction_truth[
                ["entry_id", cluster_column, "true_kcat_log10"]
            ].merge(
                methods[method][["entry_id", "prediction_log10"]],
                on="entry_id",
                how="inner",
                validate="one_to_one",
            )
            grouped = (
                merged.groupby(cluster_column, as_index=False)
                .agg(
                    true_kcat_log10=("true_kcat_log10", "median"),
                    prediction_log10=("prediction_log10", "median"),
                )
            )
            result = frame_metrics(grouped)
            rows.append(
                {
                    "analysis_scope": analysis_scope,
                    "statistical_unit": statistical_unit,
                    "aggregation_rule": aggregation_rule,
                    "scope_units": total_groups,
                    "method": method,
                    "n": result["n"],
                    "coverage_within_scope_percent": 100.0 * result["n"] / total_groups,
                    **{key: value for key, value in result.items() if key != "n"},
                }
            )
    return pd.DataFrame(rows)


def cluster_bootstrap(
    truth: pd.DataFrame,
    methods: dict[str, pd.DataFrame],
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    cluster_columns = [
        "protein_cluster",
        "pair_cluster",
        "reaction_cluster",
        "reference_cluster",
        "label_assignment_cluster",
    ]
    rows = []

    def append_scope(
        method: str,
        frame: pd.DataFrame,
        analysis_scope: str,
        seed_offset: int,
    ) -> None:
        for cluster_index, cluster_column in enumerate(cluster_columns):
            grouped = frame.groupby(cluster_column)["abs_error_log10"].agg(["sum", "count"])
            cluster_ids = grouped.index.to_numpy()
            sums = grouped["sum"].to_numpy(float)
            counts = grouped["count"].to_numpy(float)
            local_seed = seed + seed_offset + cluster_index
            rng = np.random.default_rng(local_seed)
            estimates = np.empty(replicates)
            for replicate in range(replicates):
                selected = rng.integers(0, len(cluster_ids), size=len(cluster_ids))
                estimates[replicate] = sums[selected].sum() / counts[selected].sum()
            low, high = np.percentile(estimates, [2.5, 97.5])
            rows.append(
                {
                    "analysis_scope": analysis_scope,
                    "method": method,
                    "cluster_type": cluster_column.replace("_cluster", ""),
                    "n_rows": len(frame),
                    "n_clusters": len(cluster_ids),
                    "row_weighted_mae_log10": float(frame["abs_error_log10"].mean()),
                    "cluster_bootstrap_ci_low_95": float(low),
                    "cluster_bootstrap_ci_high_95": float(high),
                    "bootstrap_replicates": replicates,
                    "seed": local_seed,
                    "estimand": "row-weighted MAE with clusters resampled as intact units",
                }
            )

    for method_index, method in enumerate(METHOD_ORDER):
        frame = methods[method][["entry_id", "abs_error_log10"]].merge(
            truth[["entry_id", *cluster_columns]],
            on="entry_id",
            how="left",
            validate="one_to_one",
        )
        append_scope(method, frame, "achieved_evaluation_set", method_index * 100)

    reaction_ids = set.intersection(
        *(set(methods[method]["entry_id"].astype(str)) for method in ["KcatNet", "TurNuP", "PMAK"])
    )
    kcatnet_reaction = methods["KcatNet"].loc[
        methods["KcatNet"]["entry_id"].isin(reaction_ids),
        ["entry_id", "abs_error_log10"],
    ].merge(
        truth[["entry_id", *cluster_columns]],
        on="entry_id",
        how="left",
        validate="one_to_one",
    )
    append_scope(
        "KcatNet",
        kcatnet_reaction,
        f"reaction_aware_common_{len(reaction_ids)}",
        5000,
    )

    sabiork_ambiguous_ids = set(
        truth.loc[
            truth["experimental_substrate_support"].eq("participant_ambiguous"),
            "entry_id",
        ].astype(str)
    )
    for method_index, method in enumerate(METHOD_ORDER):
        frame = methods[method].loc[
            ~methods[method]["entry_id"].isin(sabiork_ambiguous_ids),
            ["entry_id", "abs_error_log10"],
        ].merge(
            truth[["entry_id", *cluster_columns]],
            on="entry_id",
            how="left",
            validate="one_to_one",
        )
        append_scope(
            method,
            frame,
            "exclude_sabio_participant_ambiguous",
            10000 + method_index * 100,
        )
    return pd.DataFrame(rows)


def bh_adjust(values: list[float]) -> list[float]:
    count = len(values)
    order = np.argsort(values)
    adjusted = np.ones(count)
    running = 1.0
    for rank in range(count, 0, -1):
        index = int(order[rank - 1])
        running = min(running, values[index] * count / rank)
        adjusted[index] = running
    return adjusted.tolist()


def cluster_wilcoxon(truth: pd.DataFrame, methods: dict[str, pd.DataFrame]) -> pd.DataFrame:
    comparisons = [
        ("broad_common", "KcatNet", "CataPro"),
        ("reaction_common", "KcatNet", "TurNuP"),
        ("reaction_common", "KcatNet", "PMAK"),
        ("reaction_common", "TurNuP", "PMAK"),
    ]
    cluster_columns = [
        "protein_cluster",
        "pair_cluster",
        "reaction_cluster",
        "reference_cluster",
        "label_assignment_cluster",
    ]
    rows = []
    p_values = []
    for scope, method_a, method_b in comparisons:
        left = methods[method_a][["entry_id", "abs_error_log10"]].rename(
            columns={"abs_error_log10": "error_a"}
        )
        right = methods[method_b][["entry_id", "abs_error_log10"]].rename(
            columns={"abs_error_log10": "error_b"}
        )
        common = left.merge(right, on="entry_id", validate="one_to_one").merge(
            truth[["entry_id", *cluster_columns]],
            on="entry_id",
            validate="one_to_one",
        )
        for cluster_column in cluster_columns:
            grouped = common.groupby(cluster_column)[["error_a", "error_b"]].mean()
            try:
                test = wilcoxon(grouped["error_a"], grouped["error_b"], alternative="two-sided")
                statistic = float(test.statistic)
                p_value = float(test.pvalue)
            except ValueError:
                statistic = np.nan
                p_value = 1.0
            p_values.append(p_value)
            rows.append(
                {
                    "comparison_scope": scope,
                    "method_a": method_a,
                    "method_b": method_b,
                    "cluster_type": cluster_column.replace("_cluster", ""),
                    "n_common_rows": len(common),
                    "n_paired_clusters": len(grouped),
                    "cluster_mean_error_a": float(grouped["error_a"].mean()),
                    "cluster_mean_error_b": float(grouped["error_b"].mean()),
                    "wilcoxon_statistic": statistic,
                    "p_value_raw": p_value,
                    "estimand": "paired distribution of cluster-mean absolute errors",
                }
            )
    adjusted = bh_adjust(p_values)
    for row, q_value in zip(rows, adjusted):
        row["p_value_bh_global"] = q_value
        row["significant_bh_fdr_0.05"] = bool(q_value < 0.05)
    return pd.DataFrame(rows)


def load_matcher_module():
    spec = importlib.util.spec_from_file_location("kcat_matcher", MATCHER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {MATCHER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def selected_measurements(
    matcher,
    entries: list[dict[str, str]],
    records: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    records_by_species_ec: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for record in records:
        records_by_species_ec[(record["species"], record["ec_number"])].append(record)
    selected_by_entry = {}
    for entry in entries:
        matches = []
        for ec_number in matcher.split_values(entry.get("ec_number", "")):
            for record in records_by_species_ec.get((entry["species"], ec_number), []):
                level = matcher.substrate_match(entry, record)
                if level:
                    matches.append((matcher.MATCH_STRENGTH[level], record))
        if matches:
            strength = max(item[0] for item in matches)
            selected_by_entry[entry["entry_id"]] = [
                item[1] for item in matches if item[0] == strength
            ]
    return selected_by_entry


def build_measurement_audit(
    truth: pd.DataFrame,
    skip_brenda_reparse: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matcher = load_matcher_module()
    entries = matcher.read_rows(ENTRY_PATH)
    brenda = matcher.read_rows(matcher.BRENDA_RAW)
    sabio = matcher.read_sabio_as_generic()
    selected = selected_measurements(matcher, entries, brenda + sabio)
    detail_rows = []
    for _, benchmark_row in truth.iterrows():
        records = selected.get(str(benchmark_row["entry_id"]), [])
        values = np.asarray(
            [
                float(record["kcat"])
                for record in records
                if float(record.get("kcat", 0) or 0) > 0
            ],
            dtype=float,
        )
        logs = np.log10(values) if len(values) else np.asarray([], dtype=float)
        median = float(np.median(logs)) if len(logs) else np.nan
        sources = sorted({record["source_database"] for record in records})
        source_medians = {}
        for source in sources:
            source_values = [
                math.log10(float(record["kcat"]))
                for record in records
                if record["source_database"] == source and float(record["kcat"]) > 0
            ]
            source_medians[source] = float(np.median(source_values))
        p_h = [
            matcher.parse_float(record.get("pH", ""))
            for record in records
        ]
        temperature = [
            matcher.parse_float(record.get("temperature_c", ""))
            for record in records
        ]
        p_h = [value for value in p_h if value is not None]
        temperature = [value for value in temperature if value is not None]
        detail_rows.append(
            {
                "entry_id": benchmark_row["entry_id"],
                "selected_measurements": len(values),
                "source_databases": ";".join(sources),
                "log10_range": float(logs.max() - logs.min()) if len(logs) > 1 else np.nan,
                "log10_median_absolute_deviation": (
                    float(np.median(np.abs(logs - median))) if len(logs) > 1 else np.nan
                ),
                "log10_mean_absolute_deviation_from_entry_median": (
                    float(np.mean(np.abs(logs - median))) if len(logs) > 1 else np.nan
                ),
                "brenda_sabiork_median_difference_abs_log10": (
                    abs(source_medians["BRENDA"] - source_medians["SABIO-RK"])
                    if {"BRENDA", "SABIO-RK"}.issubset(source_medians)
                    else np.nan
                ),
                "pH_range": max(p_h) - min(p_h) if len(p_h) > 1 else np.nan,
                "temperature_range_c": (
                    max(temperature) - min(temperature) if len(temperature) > 1 else np.nan
                ),
            }
        )
    detail = pd.DataFrame(detail_rows)
    repeated = detail[detail["selected_measurements"] >= 2]
    cross = detail[detail["brenda_sabiork_median_difference_abs_log10"].notna()]
    summary = pd.DataFrame(
        [
            {
                "analysis": "records_with_multiple_selected_measurements",
                "n": len(repeated),
                "denominator": len(detail),
                "median_abs_dispersion_log10": repeated[
                    "log10_median_absolute_deviation"
                ].median(),
                "mean_abs_dispersion_log10": repeated[
                    "log10_mean_absolute_deviation_from_entry_median"
                ].mean(),
                "median_cross_database_difference_log10": np.nan,
                "dispersion_formula": (
                    "For record r: d_r = median_i(|log10(kcat_ri) - "
                    "median_j(log10(kcat_rj))|); reported value = median_r(d_r)."
                ),
                "interpretation": (
                    "Dispersion among measurements selected at the strongest available matching level; "
                    "only records with at least two selected positive measurements enter this summary."
                ),
            },
            {
                "analysis": "records_with_brenda_and_sabiork_selected_measurements",
                "n": len(cross),
                "denominator": len(detail),
                "median_abs_dispersion_log10": np.nan,
                "mean_abs_dispersion_log10": np.nan,
                "median_cross_database_difference_log10": cross[
                    "brenda_sabiork_median_difference_abs_log10"
                ].median(),
                "dispersion_formula": (
                    "For record r: |median_i_in_BRENDA(log10(kcat_ri)) - "
                    "median_j_in_SABIO-RK(log10(kcat_rj))|; reported value = median across records."
                ),
                "interpretation": (
                    "Absolute difference between source-specific medians; source support is not "
                    "treated as independent replication."
                ),
            },
        ]
    )

    current_mutant_comments = sum(
        matcher.is_mutant_comment(record.get("comment", "")) for record in brenda
    )
    mutation_rows = [
        {
            "source": "BRENDA",
            "audit_status": "verified_from_parser_and_cached_records",
            "counting_stage": (
                "Parsed BRENDA turnover-number candidate records after EC/substrate parsing; "
                "the before/after counts differ only by the comment-based mutation screen."
            ),
            "raw_before_filter": np.nan,
            "excluded_by_mutation_screen": np.nan,
            "retained_after_screen": len(brenda),
            "mutation_flagged_in_retained_cache": current_mutant_comments,
            "benchmark_policy": (
                "Comments matching mutant, mutation, or variant are excluded by default."
            ),
        },
        {
            "source": "SABIO-RK",
            "audit_status": "unknown_mutation_status",
            "counting_stage": (
                "Cached parsed SABIO-RK rows; no mutation-status field is available, so no "
                "mutation screen can be applied or quantified."
            ),
            "raw_before_filter": len(sabio),
            "excluded_by_mutation_screen": np.nan,
            "retained_after_screen": len(sabio),
            "mutation_flagged_in_retained_cache": np.nan,
            "benchmark_policy": (
                "Cached export has no mutation/variant field; wild-type status cannot be verified."
            ),
        },
    ]
    if not skip_brenda_reparse:
        args = argparse.Namespace(
            brenda_json=str(matcher.BRENDA_JSON_TAR),
            ckb_db=str(matcher.CKB_DB),
            include_mutants=True,
        )
        brenda_with_mutants = matcher.parse_brenda_records(entries, args)
        current_ids = {record["source_record_id"] for record in brenda}
        excluded = [
            record for record in brenda_with_mutants
            if record["source_record_id"] not in current_ids
            and matcher.is_mutant_comment(record.get("comment", ""))
        ]
        mutation_rows[0]["raw_before_filter"] = len(brenda_with_mutants)
        mutation_rows[0]["excluded_by_mutation_screen"] = len(excluded)
        reconciled = len(brenda_with_mutants) - len(excluded)
        if reconciled != len(brenda):
            raise ValueError(
                "BRENDA mutation-stage counts do not reconcile: "
                f"before={len(brenda_with_mutants)}, excluded={len(excluded)}, retained={len(brenda)}"
            )
    return summary, detail, pd.DataFrame(mutation_rows)


def build_substrate_direction_audit(truth: pd.DataFrame, context: pd.DataFrame) -> pd.DataFrame:
    model_reactions = pd.read_csv(MODEL_REACTION_PATH)
    frame = truth.merge(
        context[["entry_id", "currency_or_cofactor_like_by_name"]],
        on="entry_id",
        how="left",
        validate="one_to_one",
    ).merge(
        model_reactions[["species", "reaction_id", "is_reversible"]],
        on=["species", "reaction_id"],
        how="left",
        validate="many_to_one",
    )
    if frame["is_reversible"].isna().any():
        raise ValueError("Missing model reversibility annotation for final benchmark rows")
    names = frame["substrate_name"].fillna("").str.lower()
    reaction_ids = frame["reaction_id"].fillna("").astype(str)
    reaction_names = frame["reaction_name"].fillna("").str.lower()
    reversible = bool_series(frame["is_reversible"])
    rows = []

    def add(item: str, mask: pd.Series, interpretation: str) -> None:
        rows.append(
            {
                "audit_item": item,
                "count": int(mask.sum()),
                "denominator": len(frame),
                "percent": 100.0 * mask.mean(),
                "interpretation": interpretation,
            }
        )

    for direction, count in frame["reaction_direction"].fillna("unknown").value_counts().items():
        add(
            f"reaction_direction_{direction}",
            frame["reaction_direction"].fillna("unknown").eq(direction),
            (
                "Model-side encoding used to construct the candidate reactant input; this is not "
                "evidence that the experimental assay measured the same direction."
            ),
        )
    add(
        "model_reversible_rows",
        reversible,
        "Rows linked to a model reaction whose lower bound permits reverse flux.",
    )
    reversible_reactions = frame.loc[
        reversible, ["species", "reaction_id"]
    ].drop_duplicates()
    all_reactions = frame[["species", "reaction_id"]].drop_duplicates()
    rows.append(
        {
            "audit_item": "model_reversible_unique_reactions",
            "count": len(reversible_reactions),
            "denominator": len(all_reactions),
            "percent": 100.0 * len(reversible_reactions) / len(all_reactions),
            "interpretation": (
                "Unique model reactions represented by at least one benchmark row and marked reversible "
                "by model bounds. Experimental assay direction remains unavailable."
            ),
        }
    )
    add(
        "experimental_assay_direction_unverified",
        pd.Series(True, index=frame.index),
        "The cached experimental records do not provide a direction that can be reconciled to model encoding.",
    )
    for selection, count in frame["substrate_selection"].fillna("unknown").value_counts().items():
        add(
            f"substrate_selection_{selection}",
            frame["substrate_selection"].fillna("unknown").eq(selection),
            (
                "All model-encoded reactants are retained before experimental substrate matching; "
                "this field records that candidate policy and is not a biochemical role label."
            ),
        )
    add(
        "monatomic_proton_selected",
        names.isin({"h+", "proton"}),
        "Chemically valid model reactant but generally not a defensible primary assay substrate.",
    )
    add(
        "water_selected",
        names.isin({"h2o", "water"}),
        "Chemically valid model reactant but often a currency participant.",
    )
    add(
        "registry_currency_or_cofactor",
        frame["substrate_role_group"].eq("currency_or_cofactor"),
        "Joint registry classification using normalized names, external identifiers, and standardized structures.",
    )
    add(
        "registry_carrier_linked_variable",
        frame["substrate_role_group"].eq("carrier_linked_variable"),
        "Carrier-linked metabolites are reported separately because they may be the variable substrate.",
    )
    add(
        "name_only_currency_or_cofactor_exploratory",
        bool_series(frame["currency_or_cofactor_like_by_name"]),
        "Exploratory name-only classification; not used for the primary role stratification.",
    )
    add(
        "exchange_demand_or_sink_reaction",
        reaction_ids.str.match(r"^(EX_|DM_|SK_)"),
        "Boundary reactions were not expected to survive the gene/sequence requirements.",
    )
    add(
        "transport_named_reaction",
        reaction_names.str.contains("transport|transporter", regex=True),
        "No reaction with transport or transporter in its name survived the final filters.",
    )
    add(
        "spontaneous_named_reaction",
        reaction_names.str.contains("spontaneous|spont", regex=True),
        "No reaction marked spontaneous by its name survived the final filters.",
    )
    add(
        "sabiork_only_participant_level_label",
        frame["experimental_substrate_support"].eq("participant_ambiguous"),
        (
            "SABIO-RK export lists all reaction participants; the cached data do not identify which "
            "participant was the assayed variable substrate."
        ),
    )
    add(
        "brenda_specific_substrate_support_present",
        frame["experimental_substrate_support"].eq("substrate_supported"),
        "BRENDA turnover-number records provide a substrate-specific field.",
    )
    return pd.DataFrame(rows)


def load_json_frame(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as handle:
        return pd.DataFrame(json.load(handle))


def load_catpred_corpus() -> pd.DataFrame:
    archive_path = BASE / "external_methods" / "CatPred_datas" / "catpred-db.tar.gz"
    member = "CatPred-DB/data/kcat/kcat-random_trainvaltest.csv"
    with tarfile.open(archive_path, "r:gz") as archive:
        handle = archive.extractfile(member)
        if handle is None:
            raise FileNotFoundError(member)
        return pd.read_csv(io.BytesIO(handle.read()), low_memory=False)


def load_training_corpora(truth: pd.DataFrame) -> dict[str, dict[str, object]]:
    corpora: dict[str, dict[str, object]] = {}

    kcatnet = pd.read_pickle(BASE / "external_methods" / "KcatNet" / "Dataset" / "KcatNet_traindf.pkl")
    corpora["KcatNet"] = {
        "frame": kcatnet.rename(columns={"Pro_seq": "sequence", "Smile": "chemical"}),
        "mode": "substrate",
        "source": "KcatNet_traindf.pkl (fit split)",
        "status": "available_fit_split",
    }

    catapro = pd.read_csv(
        BASE / "external_methods" / "CataPro" / "datasets" / "kcat-data_0.4simi-10fold.csv"
    )
    corpora["CataPro"] = {
        "frame": catapro.rename(columns={"Sequence": "sequence", "Smiles": "chemical"}),
        "mode": "substrate",
        "source": "kcat-data_0.4simi-10fold.csv (union across released folds)",
        "status": "available_source_corpus_union",
    }

    pretkcat = pd.read_csv(
        BASE
        / "external_methods"
        / "PreTKcat"
        / "datasets"
        / "DLTKcat_data"
        / "kcat_merge_DLTKcat.csv"
    ).rename(columns={"seq": "sequence", "smiles": "chemical"})
    benchmark_exclusion = {
        truncate_1000(clean_sequence(row.sequence))
        + "|"
        + str(mol_identity(row.SMILES)["connectivity_key"])
        for row in truth.itertuples()
    }
    pretkcat["_fit_key"] = [
        truncate_1000(clean_sequence(sequence))
        + "|"
        + (
            str(identity["connectivity_key"]) if identity is not None else ""
        )
        for sequence, identity in zip(
            pretkcat["sequence"],
            [mol_identity(value) for value in pretkcat["chemical"]],
        )
    ]
    pretkcat_effective = pretkcat[~pretkcat["_fit_key"].isin(benchmark_exclusion)].copy()
    corpora["PreTKcat"] = {
        "frame": pretkcat_effective,
        "mode": "substrate",
        "source": "kcat_merge_DLTKcat.csv after benchmark exact-pair exclusion",
        "status": "available_effective_fit_corpus",
        "source_rows_before_exclusion": len(pretkcat),
    }

    dlkcat_path = (
        BASE
        / "external_methods"
        / "DLKcat_official"
        / "DeeplearningApproach"
        / "Data"
        / "database"
        / "Kcat_combination_0918.json"
    )
    dlkcat = load_json_frame(dlkcat_path).rename(
        columns={"Sequence": "sequence", "Smiles": "chemical"}
    )
    corpora["DLKcat"] = {
        "frame": dlkcat,
        "mode": "substrate",
        "source": "Kcat_combination_0918.json (source corpus; effective checkpoint split unavailable)",
        "status": "available_source_corpus_only",
    }

    unikp_path = (
        BASE
        / "external_methods"
        / "AI_file"
        / "UniKP"
        / "datasets"
        / "Kcat_combination_0918_wildtype_mutant.json"
    )
    unikp = load_json_frame(unikp_path).rename(
        columns={"Sequence": "sequence", "Smiles": "chemical"}
    )
    corpora["UniKP"] = {
        "frame": unikp,
        "mode": "substrate",
        "source": "Kcat_combination_0918_wildtype_mutant.json",
        "status": "available_source_corpus",
    }

    kinform = load_json_frame(
        BASE
        / "external_methods"
        / "KinForm"
        / "data"
        / "EITLEM_data"
        / "KCAT"
        / "kcat_data.json"
    ).rename(columns={"sequence": "sequence", "smiles": "chemical"})
    corpora["KinForm-L"] = {
        "frame": kinform,
        "mode": "substrate",
        "source": "EITLEM_data/KCAT/kcat_data.json (training on all data in released code)",
        "status": "available_fit_corpus",
    }

    catpred = load_catpred_corpus().rename(
        columns={"sequence": "sequence", "reactant_smiles": "chemical"}
    )
    corpora["CatPred"] = {
        "frame": catpred,
        "mode": "reactant_components",
        "source": "CatPred-DB kcat-random_trainvaltest.csv (production all-data corpus)",
        "status": "available_production_corpus",
    }

    pmak_paths = sorted(
        (
            BASE
            / "external_methods"
            / "PMAK"
            / "supplement"
            / "data"
            / "turnup"
            / "cold_reaction"
        ).glob("kcat_train_fold_*_with_sub.csv")
    )
    pmak_frames = []
    for path in pmak_paths:
        frame = pd.read_csv(path, low_memory=False)
        frame["chemical"] = frame["reactant_smiles"].astype(str) + ">>" + frame[
            "product_smiles"
        ].astype(str)
        pmak_frames.append(frame[["sequence", "chemical"]])
    pmak = pd.concat(pmak_frames, ignore_index=True).drop_duplicates()
    corpora["PMAK"] = {
        "frame": pmak,
        "mode": "reaction",
        "source": "union of five released reaction-cold training folds",
        "status": "available_checkpoint_training_union",
    }

    dekp = pd.read_csv(
        BASE / "external_methods" / "DEKP" / "datasets" / "kcat_dataset.csv",
        sep="	",
    )
    sequence_column = next(
        column for column in ["sequence", "Sequence", "seq"] if column in dekp.columns
    )
    chemical_column = next(
        column for column in ["smiles", "Smiles", "substrate_smiles"] if column in dekp.columns
    )
    dekp = dekp.rename(columns={sequence_column: "sequence", chemical_column: "chemical"})
    benchmark_pairs = set(truth["pair_cluster"])
    dekp["_pair_key"] = [
        clean_sequence(sequence)
        + "|"
        + (
            str(identity["connectivity_key"]) if identity is not None else ""
        )
        for sequence, identity in zip(
            dekp["sequence"],
            [mol_identity(value) for value in dekp["chemical"]],
        )
    ]
    dekp_effective = dekp[~dekp["_pair_key"].isin(benchmark_pairs)].copy()
    corpora["DEKP-public-retrained"] = {
        "frame": dekp_effective,
        "mode": "substrate",
        "source": "DEKP kcat_dataset.csv after benchmark exact-pair exclusion",
        "status": "available_effective_fit_corpus",
        "source_rows_before_exclusion": len(dekp),
    }
    return corpora


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="ascii") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n{sequence}\n")


def diamond_hits(
    query_sequences: dict[str, str],
    subject_sequences: dict[str, str],
) -> dict[str, list[tuple[str, float]]]:
    with tempfile.TemporaryDirectory(prefix="kcat_overlap_") as tmp:
        tmp_path = Path(tmp)
        query_fasta = tmp_path / "query.fasta"
        subject_fasta = tmp_path / "subject.fasta"
        database = tmp_path / "training"
        output = tmp_path / "hits.tsv"
        write_fasta(query_fasta, query_sequences)
        write_fasta(subject_fasta, subject_sequences)
        subprocess.run(
            ["diamond", "makedb", "--in", str(subject_fasta), "-d", str(database), "--quiet"],
            check=True,
        )
        subprocess.run(
            [
                "diamond",
                "blastp",
                "-d",
                str(database),
                "-q",
                str(query_fasta),
                "-o",
                str(output),
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "--id",
                "30",
                "--query-cover",
                "50",
                "--subject-cover",
                "50",
                "--max-target-seqs",
                "250",
                "--sensitive",
                "--threads",
                "8",
                "--quiet",
            ],
            check=True,
        )
        hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
        if output.exists() and output.stat().st_size:
            frame = pd.read_csv(output, sep="\t", names=["query", "subject", "identity"])
            for row in frame.itertuples(index=False):
                hits[str(row.query)].append((str(row.subject), float(row.identity)))
        return hits


def prepare_training_rows(
    frame: pd.DataFrame,
    mode: str,
) -> tuple[pd.DataFrame, dict[str, list[int]], list[object], list[str]]:
    rows = []
    fingerprints: list[object] = []
    keys: list[str] = []
    for sequence, chemical in zip(frame["sequence"], frame["chemical"]):
        sequence = clean_sequence(sequence)
        if not sequence:
            continue
        if mode == "reaction":
            identity = canonical_reaction(chemical)
            if identity is None:
                continue
            key, fingerprint = identity
            component_fps = [fingerprint]
            component_keys = [key]
        elif mode == "reactant_components":
            components = component_identities(chemical)
            if not components:
                continue
            component_fps = [component["fingerprint"] for component in components]
            component_keys = [str(component["connectivity_key"]) for component in components]
        else:
            identity = mol_identity(chemical)
            if identity is None:
                continue
            component_fps = [identity["fingerprint"]]
            component_keys = [str(identity["connectivity_key"])]
        row_index = len(rows)
        rows.append(
            {
                "sequence": sequence,
                "fingerprints": component_fps,
                "keys": component_keys,
            }
        )
        fingerprints.extend(component_fps)
        keys.extend(component_keys)
    prepared = pd.DataFrame(rows)
    sequence_to_rows: dict[str, list[int]] = defaultdict(list)
    for index, sequence in enumerate(prepared["sequence"]):
        sequence_to_rows[sequence].append(index)
    return prepared, sequence_to_rows, fingerprints, keys


def query_chemical_identity(
    row: pd.Series,
    method: str,
    mode: str,
) -> tuple[str, object] | None:
    if mode == "reaction":
        pmak_input = pd.read_csv(
            BASE / "data" / "final" / "pmak" / "pmak_kcat_input.csv"
        ).set_index("entry_id")
        if row["entry_id"] not in pmak_input.index:
            return None
        return canonical_reaction(pmak_input.loc[row["entry_id"], "reaction_smiles"])
    identity = mol_identity(row["SMILES"])
    if identity is None:
        return None
    return str(identity["connectivity_key"]), identity["fingerprint"]


def overlap_for_method(
    method: str,
    evaluation: pd.DataFrame,
    truth: pd.DataFrame,
    corpus: dict[str, object],
    skip_neighbor_search: bool,
) -> tuple[pd.DataFrame, dict[str, object]]:
    mode = str(corpus["mode"])
    training, sequence_to_rows, all_fingerprints, all_keys = prepare_training_rows(
        corpus["frame"], mode
    )
    evaluation_truth = truth[truth["entry_id"].isin(evaluation["entry_id"])].copy()
    query_sequence_ids = {
        f"q{index}": sequence
        for index, sequence in enumerate(
            sorted(set(evaluation_truth["sequence"].map(clean_sequence)))
        )
    }
    sequence_to_query_id = {sequence: identifier for identifier, sequence in query_sequence_ids.items()}
    subject_sequences = {
        f"s{index}": sequence for index, sequence in enumerate(sorted(sequence_to_rows))
    }
    subject_id_to_sequence = {
        identifier: sequence for identifier, sequence in subject_sequences.items()
    }
    if skip_neighbor_search:
        hits = {}
    else:
        hits = diamond_hits(query_sequence_ids, subject_sequences)
    all_key_set = set(all_keys)
    detail_rows = []

    pmak_input = None
    if mode == "reaction":
        pmak_input = pd.read_csv(
            BASE / "data" / "final" / "pmak" / "pmak_kcat_input.csv"
        ).set_index("entry_id")

    for _, row in evaluation_truth.iterrows():
        sequence = clean_sequence(row["sequence"])
        if mode == "reaction":
            if row["entry_id"] not in pmak_input.index:
                continue
            query_identity = canonical_reaction(
                pmak_input.loc[row["entry_id"], "reaction_smiles"]
            )
        else:
            molecule = mol_identity(row["SMILES"])
            query_identity = (
                (str(molecule["connectivity_key"]), molecule["fingerprint"])
                if molecule is not None
                else None
            )
        if query_identity is None:
            continue
        query_key, query_fp = query_identity
        exact_sequence = sequence in sequence_to_rows
        exact_chemical_any = query_key in all_key_set
        exact_pair = False
        if exact_sequence:
            exact_pair = any(
                query_key in training.iloc[index]["keys"]
                for index in sequence_to_rows[sequence]
            )
        # The global maximum against every training molecule is not used by
        # the joint-neighbor classification and dominates runtime. Preserve
        # exact chemical identity as 1.0; nonidentity maxima are intentionally
        # left uncomputed while sequence-conditioned similarities below remain
        # fully evaluated.
        chemical_any = 1.0 if exact_chemical_any else np.nan

        max_sequence_identity = 100.0 if exact_sequence else np.nan
        max_joint_chemical_at_seq80 = 0.0
        max_joint_sequence_at_chem80 = 0.0
        if exact_sequence:
            for index in sequence_to_rows[sequence]:
                similarities = DataStructs.BulkTanimotoSimilarity(
                    query_fp, training.iloc[index]["fingerprints"]
                )
                if similarities:
                    maximum = float(max(similarities))
                    max_joint_chemical_at_seq80 = max(max_joint_chemical_at_seq80, maximum)
                    if maximum >= NEAR_CHEMICAL_TANIMOTO:
                        max_joint_sequence_at_chem80 = 100.0

        query_id = sequence_to_query_id[sequence]
        for subject_id, identity in hits.get(query_id, []):
            max_sequence_identity = (
                identity if not np.isfinite(max_sequence_identity) else max(max_sequence_identity, identity)
            )
            subject_sequence = subject_id_to_sequence[subject_id]
            for index in sequence_to_rows[subject_sequence]:
                similarities = DataStructs.BulkTanimotoSimilarity(
                    query_fp, training.iloc[index]["fingerprints"]
                )
                if not similarities:
                    continue
                maximum = float(max(similarities))
                if identity >= NEAR_SEQUENCE_IDENTITY:
                    max_joint_chemical_at_seq80 = max(max_joint_chemical_at_seq80, maximum)
                if maximum >= NEAR_CHEMICAL_TANIMOTO:
                    max_joint_sequence_at_chem80 = max(max_joint_sequence_at_chem80, identity)

        near_pair = (
            not exact_pair
            and max_joint_chemical_at_seq80 >= NEAR_CHEMICAL_TANIMOTO
            and max_joint_sequence_at_chem80 >= NEAR_SEQUENCE_IDENTITY
        )
        if exact_pair:
            category = "exact_overlap"
        elif near_pair:
            category = "near_neighbor"
        elif skip_neighbor_search:
            category = "neighbor_not_computed"
        else:
            category = "no_joint_neighbor_under_thresholds"
        detail_rows.append(
            {
                "method": method,
                "entry_id": row["entry_id"],
                "training_proximity_class": category,
                "exact_sequence": exact_sequence,
                "exact_chemical_identity_anywhere": exact_chemical_any,
                "exact_sequence_chemical_pair": exact_pair,
                "max_sequence_identity_percent": max_sequence_identity,
                "max_chemical_similarity_anywhere": chemical_any,
                "max_chemical_similarity_with_sequence_identity_ge_80": (
                    max_joint_chemical_at_seq80
                ),
                "max_sequence_identity_with_chemical_similarity_ge_0.80": (
                    max_joint_sequence_at_chem80
                ),
            }
        )
    detail = pd.DataFrame(detail_rows)
    error_map = evaluation.set_index("entry_id")["abs_error_log10"]
    detail["abs_error_log10"] = detail["entry_id"].map(error_map)
    if detail["abs_error_log10"].isna().any():
        raise ValueError(f"Missing evaluated absolute errors in {method} overlap audit")
    counts = detail["training_proximity_class"].value_counts()
    class_metrics = {}
    for category, prefix in [
        ("exact_overlap", "exact_overlap"),
        ("near_neighbor", "near_neighbor"),
        ("no_joint_neighbor_under_thresholds", "no_joint_neighbor"),
    ]:
        values = detail.loc[
            detail["training_proximity_class"].eq(category), "abs_error_log10"
        ]
        class_metrics[f"{prefix}_mae_log10"] = float(values.mean()) if len(values) else np.nan
        class_metrics[f"{prefix}_within_1.0_fraction"] = (
            float((values <= 1.0).mean()) if len(values) else np.nan
        )
    summary = {
        "method": method,
        "training_corpus_status": corpus["status"],
        "training_corpus_source": corpus["source"],
        "training_rows_loaded": len(corpus["frame"]),
        "training_rows_standardized": len(training),
        "evaluated_rows": len(evaluation_truth),
        "audited_rows": len(detail),
        "exact_sequence_rows": int(detail["exact_sequence"].sum()),
        "exact_chemical_identity_rows": int(detail["exact_chemical_identity_anywhere"].sum()),
        "exact_pair_overlap_rows": int(detail["exact_sequence_chemical_pair"].sum()),
        "near_neighbor_rows": int(counts.get("near_neighbor", 0)),
        "no_joint_neighbor_rows_under_thresholds": int(
            counts.get("no_joint_neighbor_under_thresholds", 0)
        ),
        "neighbor_not_computed_rows": int(counts.get("neighbor_not_computed", 0)),
        "near_sequence_identity_threshold_percent": NEAR_SEQUENCE_IDENTITY,
        "near_chemical_tanimoto_threshold": NEAR_CHEMICAL_TANIMOTO,
        "chemical_identity": (
            "uncharged largest-fragment connectivity InChIKey; reaction mode uses direction-aware "
            "standardized reaction sides"
        ),
        "global_nonidentity_chemical_similarity_status": (
            "not_computed; exact identity retained and joint sequence-conditioned similarities computed"
        ),
        "publication_overlap_status": "unavailable_in_released_record-level corpus",
        "interpretation": (
            "Exact and thresholded-neighbor proximity, not proof that a specific record was used "
            "for gradient fitting unless the corpus status says fit split. The no-joint-neighbor "
            "class means no single training record met both thresholds; exact sequence-only or "
            "chemical-identity-only matches may still be present."
        ),
        **class_metrics,
    }
    if "source_rows_before_exclusion" in corpus:
        summary["source_rows_before_exclusion"] = corpus["source_rows_before_exclusion"]
    return detail, summary


def build_training_overlap(
    truth: pd.DataFrame,
    methods: dict[str, pd.DataFrame],
    skip_neighbor_search: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    corpora = load_training_corpora(truth)
    summaries = []
    details = []
    for method in METHOD_ORDER:
        if method in corpora:
            detail, summary = overlap_for_method(
                method,
                methods[method],
                truth,
                corpora[method],
                skip_neighbor_search,
            )
            details.append(detail)
            summaries.append(summary)
            continue
        if method == "GO-HKP":
            status = "not_applicable_functional_assignment_baseline"
            interpretation = "No sequence-substrate regression training corpus."
        else:
            status = "unknown_training_corpus_not_released_or_not_identified"
            interpretation = (
                "No auditable record-level corpus was found in the local release package; "
                "this must not be interpreted as zero overlap."
            )
        summaries.append(
            {
                "method": method,
                "training_corpus_status": status,
                "training_corpus_source": "",
                "training_rows_loaded": np.nan,
                "training_rows_standardized": np.nan,
                "evaluated_rows": len(methods[method]),
                "audited_rows": 0,
                "exact_sequence_rows": np.nan,
                "exact_chemical_identity_rows": np.nan,
                "exact_pair_overlap_rows": np.nan,
                "near_neighbor_rows": np.nan,
                "no_joint_neighbor_rows_under_thresholds": np.nan,
                "neighbor_not_computed_rows": np.nan,
                "near_sequence_identity_threshold_percent": NEAR_SEQUENCE_IDENTITY,
                "near_chemical_tanimoto_threshold": NEAR_CHEMICAL_TANIMOTO,
                "chemical_identity": "",
                "publication_overlap_status": "unavailable",
                "interpretation": interpretation,
                "exact_overlap_mae_log10": np.nan,
                "exact_overlap_within_1.0_fraction": np.nan,
                "near_neighbor_mae_log10": np.nan,
                "near_neighbor_within_1.0_fraction": np.nan,
                "no_joint_neighbor_mae_log10": np.nan,
                "no_joint_neighbor_within_1.0_fraction": np.nan,
            }
        )
    detail_frame = pd.concat(details, ignore_index=True) if details else pd.DataFrame()
    return pd.DataFrame(summaries), detail_frame


def build_pretkcat_variant_audit() -> pd.DataFrame:
    variants = [
        (
            "raw-public",
            BASE / "data/final/pretkcat/variants/raw-public/training_overlap_audit.json",
            BASE / "reports/tables/pretkcat_variants/raw_public_eval_metrics.csv",
            False,
        ),
        (
            "exact-excluded",
            BASE / "data/final/pretkcat/variants/exact-excluded/training_overlap_audit.json",
            BASE / "reports/tables/pretkcat_eval_metrics.csv",
            True,
        ),
        (
            "near-excluded",
            BASE / "data/final/pretkcat/variants/near-excluded/training_overlap_audit.json",
            BASE / "reports/tables/pretkcat_variants/near_excluded_eval_metrics.csv",
            False,
        ),
    ]
    rows = []
    for variant, audit_path, metric_path, is_primary in variants:
        with audit_path.open("r", encoding="utf-8") as handle:
            audit = json.load(handle)
        metrics = pd.read_csv(metric_path)
        overall = metrics.loc[
            metrics["group_type"].eq("all") & metrics["group"].eq("all")
        ].iloc[0]
        rows.append(
            {
                "variant": variant,
                "manuscript_role": "primary PreTKcat result" if is_primary else "sensitivity analysis",
                "training_overlap_policy": audit["training_overlap_policy"],
                "pair_identity_definition": audit["pair_identity_definition"],
                "near_neighbor_definition": audit["near_neighbor_definition"],
                "training_rows_raw_usable": audit["train_rows_raw_usable"],
                "training_rows_removed_exact_pair": audit["train_rows_removed_exact_pair"],
                "training_rows_removed_near_only": audit["train_rows_removed_near_only"],
                "training_rows_removed_total": audit["train_rows_removed_total"],
                "training_rows_fitted": audit["train_rows_fitted"],
                "benchmark_rows_overlapping_raw_source_exact_pair": audit[
                    "benchmark_rows_source_exact_pair_overlap"
                ],
                "benchmark_rows_with_raw_source_joint_neighbor": audit[
                    "benchmark_rows_with_joint_neighbor"
                ],
                "fitted_corpus_exact_pair_overlap_status": (
                    "present" if variant == "raw-public" else "excluded"
                ),
                "fitted_corpus_joint_near_neighbor_status": (
                    "excluded" if variant == "near-excluded" else "not_excluded"
                ),
                "n": int(overall["n"]),
                "mae_log10": float(overall["mae_log10"]),
                "rmse_log10": float(overall["rmse_log10"]),
                "pearson_log10": float(overall["pearson_log10"]),
                "spearman_log10": float(overall["spearman_log10"]),
                "bias_log10": float(overall["bias_log10"]),
                "within_1.0_log10_fraction": float(overall["within_1.0_log10_fraction"]),
                "audit_file": str(audit_path.relative_to(BASE)),
                "metrics_file": str(metric_path.relative_to(BASE)),
            }
        )
    frame = pd.DataFrame(rows)
    primary_mae = float(frame.loc[frame["variant"].eq("exact-excluded"), "mae_log10"].iloc[0])
    frame["mae_difference_vs_exact_excluded"] = frame["mae_log10"] - primary_mae
    return frame



def build_record_audit(
    truth: pd.DataFrame,
    context: pd.DataFrame,
    measurement_detail: pd.DataFrame,
    overlap_detail: pd.DataFrame,
) -> pd.DataFrame:
    """Materialize every row-level field needed to reproduce S16-S24."""
    model_reactions = pd.read_csv(MODEL_REACTION_PATH)
    audit = truth.copy()
    audit["label_assignment_group_size"] = audit.groupby(
        "label_assignment_cluster"
    )["entry_id"].transform("size")
    audit["label_assignment_sequence_count"] = audit.groupby(
        "label_assignment_cluster"
    )["sequence"].transform("nunique")

    audit = audit.merge(
        model_reactions[["species", "reaction_id", "is_reversible"]],
        on=["species", "reaction_id"],
        how="left",
        validate="many_to_one",
    )
    context_columns = [
        "entry_id",
        "currency_or_cofactor_like_by_name",
        "kegg_like_primary_group_short",
        "direct_kegg_pathways",
    ]
    audit = audit.merge(
        context[context_columns], on="entry_id", how="left", validate="one_to_one"
    )
    audit["currency_or_cofactor_like_registry"] = audit["substrate_role_group"].eq(
        "currency_or_cofactor"
    )
    audit["carrier_linked_variable_registry"] = audit["substrate_role_group"].eq(
        "carrier_linked_variable"
    )
    audit["sabiork_participant_ambiguous"] = audit[
        "experimental_substrate_support"
    ].eq("participant_ambiguous")
    audit["model_forward_encoding"] = audit["reaction_direction"].eq("forward")
    audit["experimental_assay_direction_status"] = "not_verifiable_from_cached_records"
    audit = audit.merge(
        measurement_detail, on="entry_id", how="left", validate="one_to_one"
    )

    overlap_fields = [
        "training_proximity_class",
        "exact_sequence",
        "exact_chemical_identity_anywhere",
        "exact_sequence_chemical_pair",
    ]
    if not overlap_detail.empty and overlap_detail.duplicated(["method", "entry_id"]).any():
        raise ValueError("Duplicate method-entry rows in training-overlap detail")
    for method in METHOD_ORDER:
        slug = method.lower().replace("-", "_")
        method_detail = overlap_detail.loc[
            overlap_detail["method"].eq(method), ["entry_id", *overlap_fields]
        ]
        if method_detail.empty:
            default_class = (
                "not_applicable_functional_assignment_baseline"
                if method == "GO-HKP"
                else "unknown_training_corpus"
            )
            audit[f"training_proximity_class__{slug}"] = default_class
            for field in overlap_fields[1:]:
                audit[f"training_{field}__{slug}"] = np.nan
            continue
        indexed = method_detail.set_index("entry_id")
        audit[f"training_proximity_class__{slug}"] = audit["entry_id"].map(
            indexed["training_proximity_class"]
        ).fillna("not_evaluated")
        for field in overlap_fields[1:]:
            audit[f"training_{field}__{slug}"] = audit["entry_id"].map(indexed[field])

    preferred = [
        "entry_id",
        "species",
        "reaction_id",
        "reaction_name",
        "reaction_direction",
        "model_forward_encoding",
        "is_reversible",
        "experimental_assay_direction_status",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_id",
        "substrate_name",
        "substrate_compartment",
        "substrate_stoichiometry",
        "substrate_reaction_side",
        "SMILES",
        "chemical_parent_key",
        "smiles_source",
        "smiles_source_id",
        "substrate_bigg_id",
        "substrate_kegg_id",
        "substrate_chebi_id",
        "substrate_metanetx_id",
        "substrate_pubchem_cid",
        "substrate_parent_inchikey",
        "substrate_parent_inchikey_connectivity",
        "substrate_structure_standardization_status",
        "substrate_selection",
        "candidate_selection_policy",
        "substrate_is_cofactor_like",
        "substrate_role_class",
        "substrate_role_group",
        "substrate_role_evidence",
        "substrate_role_evidence_types",
        "substrate_role_evidence_count",
        "substrate_role_confidence",
        "substrate_role_registry_name",
        "substrate_role_registry_structure_consistency",
        "currency_or_cofactor_like_by_name",
        "currency_or_cofactor_like_registry",
        "carrier_linked_variable_registry",
        "sabiork_participant_ambiguous",
        "experimental_substrate_support",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
        "source_record_ids",
        "selected_measurement_kcat_values",
        "selected_measurements_json",
        "measurement_log10_median_abs_deviation",
        "measurement_log10_mean_abs_deviation",
        "measurement_log10_range",
        "selected_measurements",
        "source_databases",
        "log10_range",
        "log10_median_absolute_deviation",
        "log10_mean_absolute_deviation_from_entry_median",
        "brenda_sabiork_median_difference_abs_log10",
        "pH_range",
        "temperature_range_c",
        "true_kcat",
        "true_kcat_log10",
        "sequence",
        "enzyme_complex_type",
        "protein_cluster",
        "pair_cluster",
        "reaction_cluster",
        "reference_cluster",
        "label_assignment_cluster",
        "label_assignment_group_size",
        "label_assignment_sequence_count",
        "kegg_like_primary_group_short",
        "direct_kegg_pathways",
    ]
    remaining = [column for column in audit.columns if column not in preferred]
    audit = audit[[*preferred, *remaining]]
    if len(audit) != len(truth) or audit["entry_id"].nunique() != len(truth):
        raise ValueError(
            f"Record audit does not preserve the {len(truth)} unique benchmark rows"
        )
    return audit

def write_outputs(
    frames: dict[str, pd.DataFrame],
    measurement_detail: pd.DataFrame,
    overlap_detail: pd.DataFrame,
    record_audit: pd.DataFrame,
) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    DETAIL_DIR.mkdir(parents=True, exist_ok=True)
    for name, frame in frames.items():
        frame.to_csv(TABLE_DIR / f"{name}.csv", index=False)
    measurement_detail.to_csv(
        DETAIL_DIR / "measurement_dispersion_per_entry.csv", index=False
    )
    overlap_detail.to_csv(
        DETAIL_DIR / "training_overlap_per_method_record.csv", index=False
    )
    record_audit.to_csv(TABLE_DIR / "Record_audit.csv", index=False)
    record_audit.to_csv(DETAIL_DIR / "record_audit.csv", index=False)
    frames["S16_Label_audit"].to_csv(
        REPORT_TABLE_DIR / "submission_label_and_independence_audit.csv", index=False
    )
    frames["S20_Training_overlap"].to_csv(
        REPORT_TABLE_DIR / "submission_training_overlap_summary.csv", index=False
    )


def main() -> None:
    args = parse_args()
    truth = add_pair_keys(pd.read_csv(TRUTH_PATH))
    context = pd.read_csv(CONTEXT_PATH)
    methods = load_methods()

    measurement_summary, measurement_detail, mutation_audit = build_measurement_audit(
        truth, args.skip_brenda_reparse
    )
    overlap_summary, overlap_detail = build_training_overlap(
        truth, methods, args.skip_neighbor_search
    )
    record_audit = build_record_audit(
        truth, context, measurement_detail, overlap_detail
    )
    frames = {
        "S16_Label_audit": build_label_audit(truth, context),
        "S17_Sensitivity_subsets": build_sensitivity(truth, context, methods),
        "S18_Cluster_bootstrap": cluster_bootstrap(
            truth, methods, args.bootstrap_replicates, args.seed
        ),
        "S19_Cluster_wilcoxon": cluster_wilcoxon(truth, methods),
        "S20_Training_overlap": overlap_summary,
        "S21_Measurement_dispersion": measurement_summary,
        "S22_Mutation_status": mutation_audit,
        "S23_Substrate_direction": build_substrate_direction_audit(truth, context),
        "S24_PreTKcat_variants": build_pretkcat_variant_audit(),
    }
    write_outputs(frames, measurement_detail, overlap_detail, record_audit)
    print(f"Wrote submission audit tables: {TABLE_DIR}")
    print(f"Wrote detailed audit records: {DETAIL_DIR}")
    print(
        overlap_summary[
            [
                "method",
                "training_corpus_status",
                "exact_pair_overlap_rows",
                "near_neighbor_rows",
                "no_joint_neighbor_rows_under_thresholds",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
