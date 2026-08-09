# Publication Data and Reproducibility Artifacts

This directory contains the public, manuscript-matched data products for
benchmark version 1.2.0, artifact revision 1.2.0-r3. Author drafts, archived
manuscripts, internal review notes, and editor-generated files are intentionally
not part of the public repository.

## Public Contents

- `tables_v1.2.0/`: Supplementary Data Table 0, Table 0W, Record_audit,
  Tables 1-2, and Supplementary Tables S1-S24 as CSV files.
- `figures/`: Figures 1-4 in PNG and PDF formats.
- `submission_audit_details_v1.2.0/`: record-level measurement-dispersion,
  mapping/dependence, and training-proximity audit tables.
- `kcat_benchmark_reorganized_tables_reviewed_v1.2.0.xlsx`: the consolidated
  supplementary workbook.
- `kcat_benchmark_audit_checks_v1.2.0.csv`: machine-readable results from the
  299-check artifact audit.
- `paper_statistics_v1.2.0.json`: compact facts used by manuscript and figure
  builders.
- `independent_cluster_inference_v1.2.0-r3.csv`: an independent reconstruction
  of all 20 S19 dependence-aware comparisons.

## Rebuilding Tables and Figures

Run from the repository root after restoring the method-level result assets
listed in `zenodo_assets_manifest.csv`:

```bash
python paper/build_submission_audits.py
python paper/build_table0.py
python paper/rebuild_paper_tables.py
python paper/generate_manuscript_figures.py
python paper/recalculate_cluster_inference_v1_2.py
```

`build_submission_audits.py` additionally requires RDKit, SciPy, DIAMOND, and
the publicly available method training corpora documented in
`external_methods/METHOD_SOURCES.md`. Missing upstream training corpora are
reported as unknown rather than interpreted as zero overlap.

The reviewed manuscript itself is managed through the journal submission
workflow and is not distributed as repository source. All numerical claims
needed to audit the benchmark are materialized in the CSV and XLSX files above.
