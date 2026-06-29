#!/usr/bin/env python3
"""Export tables used by the benchmark reports into one organized directory."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
REPORTS = BASE / "reports"
TABLES = REPORTS / "tables"
EXPORT_ROOT = REPORTS / "report_tables"


MAIN_REPORT_COPIES = [
    ("method_eval_summary.csv", "Overall cross-method summary."),
    ("method_eval_summary_annotated.csv", "Cross-method summary with report grouping annotations."),
    ("method_group_annotation.csv", "Method grouping, modality, coverage, and notes."),
    ("method_rank_current_benchmark.csv", "Current-method ranking table used in the main report."),
    ("method_rank_all_with_legacy.csv", "All-method ranking table including legacy overlap rows."),
    ("species_mae_matrix.csv", "Species-level MAE matrix used for the species heatmap."),
    ("source_database_mae_matrix.csv", "Source-database MAE matrix used for the source heatmap."),
    ("go_hkp_eval_readiness.csv", "GO-HKP readiness and coverage summary."),
    ("go_hkp_eval_metrics.csv", "GO-HKP overall and stratified metrics."),
]

METHOD_METRIC_COPIES = [
    "dlkcat_official_eval_metrics.csv",
    "unikp_official_eval_metrics.csv",
    "mtlkp_eval_metrics.csv",
    "turnup_eval_metrics.csv",
    "catpred_eval_metrics.csv",
    "catapro_eval_metrics.csv",
    "pmak_eval_metrics.csv",
    "kinform_eval_metrics.csv",
    "kcatnet_eval_metrics.csv",
    "pretkcat_eval_metrics.csv",
    "dekp_public_retrained_eval_metrics.csv",
    "selfprot_eval_metrics.csv",
    "go_hkp_eval_metrics.csv",
    "MTLKP_legacy_overlap_eval_metrics.csv",
    "TurNuP_legacy_overlap_eval_metrics.csv",
]

CONTEXT_REPORT_COPIES = [
    ("benchmark_build_funnel.csv", "Model-to-benchmark construction funnel by species."),
    ("project_directory_analysis_map.csv", "Project directory map organized by analysis type."),
    ("benchmark_ready_catpred_enriched_context.csv", "Enriched benchmark context table."),
    ("benchmark_dataset_species_summary.csv", "Species distribution and core dataset statistics."),
    ("benchmark_dataset_kcat_stats_by_species.csv", "Experimental kcat distribution by species."),
    ("benchmark_dataset_source_by_species.csv", "Experimental source database distribution by species."),
    ("benchmark_dataset_match_level_by_species.csv", "Truth matching level distribution by species."),
    ("benchmark_dataset_enzyme_complex_type_by_species.csv", "Enzyme complex type distribution by species."),
    ("benchmark_dataset_substrate_role_by_species.csv", "Currency/cofactor-like substrate summary by species."),
    ("benchmark_dataset_ec_class_summary.csv", "EC class coverage summary."),
    ("benchmark_dataset_top_reactions.csv", "Top reaction records in the benchmark."),
    ("benchmark_dataset_top_substrates.csv", "Top substrate records in the benchmark."),
    ("benchmark_dataset_kegg_like_primary_group.csv", "KEGG-like primary functional group distribution."),
    ("benchmark_dataset_kegg_like_module_membership.csv", "KEGG-like module membership distribution."),
    ("benchmark_dataset_direct_yeast_kegg_pathways.csv", "Direct yeast-GEM KEGG pathway summary."),
    ("method_technical_comparison.csv", "Technical comparison of prediction methods."),
    ("method_comparison_dimensions.csv", "Recommended method comparison dimensions."),
    ("method_eval_summary_annotated.csv", "Performance summary joined into the technical method table."),
]

CORE_BENCHMARK_COPIES = [
    (BASE / "data" / "final" / "experimental_kcat_truth.csv", "Experimental kcat truth table."),
    (BASE / "data" / "final" / "benchmark_ready_truth.csv", "Truth-only benchmark-ready table."),
    (BASE / "data" / "final" / "benchmark_ready_catpred.csv", "Unified sequence+SMILES benchmark master table."),
]


def copy_table(src: Path, dst: Path, manifest: list[dict[str, str]], report: str, description: str) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    manifest.append(
        {
            "report": report,
            "table_file": dst.relative_to(EXPORT_ROOT).as_posix(),
            "source_path": src.relative_to(BASE).as_posix(),
            "description": description,
            "kind": "copied",
        }
    )


def write_table(
    df: pd.DataFrame,
    dst: Path,
    manifest: list[dict[str, str]],
    report: str,
    description: str,
    source: str,
) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dst, index=False)
    manifest.append(
        {
            "report": report,
            "table_file": dst.relative_to(EXPORT_ROOT).as_posix(),
            "source_path": source,
            "description": description,
            "kind": "generated",
        }
    )


def main_group_definitions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "group": "Broad sequence+SMILES",
                "criterion": "Enzyme sequence plus single-substrate SMILES; nearly full 978-row coverage except invalid SMILES.",
                "methods": "DLKcat-official, UniKP-official, MTLKP-official, CataPro, KcatNet, PreTKcat, SELFprot",
            },
            {
                "group": "Reaction-aware subset",
                "criterion": "Requires complete reaction information: reactant SMILES, product SMILES, and enzyme sequence.",
                "methods": "TurNuP-official, PMAK",
            },
            {
                "group": "Model-specific subset",
                "criterion": "Official inference flow or model assets restrict the processable benchmark subset.",
                "methods": "CatPred, KinForm",
            },
            {
                "group": "Public-data retrained",
                "criterion": "Publicly reproducible retraining or completion rather than original private/official best weights.",
                "methods": "DEKP-public-retrained",
            },
            {
                "group": "Function-assignment GO baseline",
                "criterion": "Non-regression GO hierarchy assignment baseline; E. coli and yeast use different GO sources.",
                "methods": "GO-HKP",
            },
            {
                "group": "Legacy overlap",
                "criterion": "Early E. coli overlap outputs retained for traceability, not formal ranking.",
                "methods": "MTLKP-legacy-overlap, TurNuP-legacy-overlap",
            },
        ]
    )


def main_overall_results() -> pd.DataFrame:
    rank = pd.read_csv(TABLES / "method_rank_current_benchmark.csv")
    keep = [
        "method",
        "group_cn",
        "modality",
        "n",
        "coverage_percent",
        "mae_log10",
        "rmse_log10",
        "pearson_log10",
        "spearman_log10",
        "within_1.0_log10_fraction",
        "bias_log10",
    ]
    return rank[[column for column in keep if column in rank.columns]].copy()


def main_legacy_results() -> pd.DataFrame:
    summary = pd.read_csv(TABLES / "method_eval_summary_annotated.csv")
    legacy = summary[summary["role"].eq("legacy")].copy()
    keep = [
        "method",
        "n",
        "coverage_percent",
        "mae_log10",
        "rmse_log10",
        "pearson_log10",
        "spearman_log10",
        "within_1.0_log10_fraction",
    ]
    return legacy[[column for column in keep if column in legacy.columns]].copy()


def benchmark_file_roles() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "file": "experimental_kcat_truth.csv",
                "purpose": "Experimental kcat truth table covering all matched BRENDA/SABIO-RK model records.",
                "model_input_recommendation": "Do not use directly as model input; sequence/SMILES readiness is not guaranteed.",
            },
            {
                "file": "benchmark_ready_truth.csv",
                "purpose": "Truth-only 978-row benchmark-ready subset used for evaluation.",
                "model_input_recommendation": "Do not feed to prediction models because it intentionally omits model input fields.",
            },
            {
                "file": "benchmark_ready_catpred.csv",
                "purpose": "Unified 978-row sequence+SMILES benchmark master table with truth and metadata.",
                "model_input_recommendation": "Use as the benchmark mother table, then extract method-specific input columns without truth leakage.",
            },
        ]
    )


def write_readme(manifest: pd.DataFrame) -> None:
    lines = [
        "# Report Table Exports",
        "",
        "This directory collects the tables used by:",
        "",
        "- `reports/kcat_benchmark_analysis_report.md`",
        "- `reports/kcat_benchmark_dataset_and_method_context.md`",
        "",
        "Subdirectories:",
        "",
        "- `main_report/`: tables for the unified benchmark performance report.",
        "- `main_report/method_metrics/`: method-level metric tables used to build summary matrices and figures.",
        "- `dataset_method_context_report/`: tables for dataset context and method technical comparison.",
        "- `dataset_method_context_report/core_benchmark_tables/`: core benchmark CSVs referenced by the context report.",
        "",
        "Use `manifest.csv` for the source path and a short description of each exported table.",
        "",
        f"Exported table count: {len(manifest)}",
        "",
    ]
    (EXPORT_ROOT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    manifest: list[dict[str, str]] = []

    for filename, description in MAIN_REPORT_COPIES:
        copy_table(
            TABLES / filename,
            EXPORT_ROOT / "main_report" / filename,
            manifest,
            "kcat_benchmark_analysis_report",
            description,
        )
    for filename in METHOD_METRIC_COPIES:
        copy_table(
            TABLES / filename,
            EXPORT_ROOT / "main_report" / "method_metrics" / filename,
            manifest,
            "kcat_benchmark_analysis_report",
            "Method-specific metric table used by the main report.",
        )
    write_table(
        main_group_definitions(),
        EXPORT_ROOT / "main_report" / "main_report_group_definitions.csv",
        manifest,
        "kcat_benchmark_analysis_report",
        "Group definition table embedded in the main report.",
        "generated from src/46_generate_benchmark_report.py report logic",
    )
    write_table(
        main_overall_results(),
        EXPORT_ROOT / "main_report" / "main_report_overall_results.csv",
        manifest,
        "kcat_benchmark_analysis_report",
        "Overall result table embedded in the main report.",
        "reports/tables/method_rank_current_benchmark.csv",
    )
    write_table(
        main_legacy_results(),
        EXPORT_ROOT / "main_report" / "main_report_legacy_overlap_results.csv",
        manifest,
        "kcat_benchmark_analysis_report",
        "Legacy overlap table embedded in the main report.",
        "reports/tables/method_eval_summary_annotated.csv",
    )

    for filename, description in CONTEXT_REPORT_COPIES:
        copy_table(
            TABLES / filename,
            EXPORT_ROOT / "dataset_method_context_report" / filename,
            manifest,
            "kcat_benchmark_dataset_and_method_context",
            description,
        )
    for src, description in CORE_BENCHMARK_COPIES:
        copy_table(
            src,
            EXPORT_ROOT / "dataset_method_context_report" / "core_benchmark_tables" / src.name,
            manifest,
            "kcat_benchmark_dataset_and_method_context",
            description,
        )
    write_table(
        benchmark_file_roles(),
        EXPORT_ROOT / "dataset_method_context_report" / "benchmark_file_role_definitions.csv",
        manifest,
        "kcat_benchmark_dataset_and_method_context",
        "Benchmark file role table embedded in the context report.",
        "generated from reports/kcat_benchmark_dataset_and_method_context.md",
    )

    manifest_df = pd.DataFrame(manifest).sort_values(["report", "table_file"])
    manifest_df.to_csv(EXPORT_ROOT / "manifest.csv", index=False)
    write_readme(manifest_df)
    print(f"Wrote report table export directory: {EXPORT_ROOT}")
    print(f"Exported tables: {len(manifest_df)}")


if __name__ == "__main__":
    main()
