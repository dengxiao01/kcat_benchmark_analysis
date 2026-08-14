# Publication Data and Reproducibility Artifacts

This directory contains the public, manuscript-matched products for benchmark
version 1.2.0. The numerical data artifact remains revision 1.2.0-r3. The
current manuscript presentation is the author-approved 0814/V4 snapshot dated
2026-08-14.

Author drafts, cover letters, archived manuscripts, internal review notes, and
editor-generated metadata are intentionally not part of the public repository.

## 0814/V4 Submission Artifacts

- [`Supplementary_tables.xlsx`](Supplementary_tables.xlsx): exact formatted
  0814 workbook with `Index` followed by `Table S1` through `Table S23`.
- [`supplementary_tables_0814/`](supplementary_tables_0814/): deterministic CSV
  mirror of every workbook sheet for browser preview and programmatic reuse.
- [`figures/`](figures/): four final high-resolution composite PNGs and 20
  author-supplied standalone panel PNGs.
- [`FIGURE_CAPTIONS.md`](FIGURE_CAPTIONS.md): complete Figure 1-4 captions and
  interpretation boundaries from the V4 manuscript.
- [`figure_sources_0814/`](figure_sources_0814/): expanded plotting code and
  compact input data from the four author-supplied 0814 figure packages.
- [`generate_manuscript_figures.py`](generate_manuscript_figures.py): validates
  exact hashes/dimensions and can rerun every data-backed plotting package.
- [`export_supplementary_tables_0814.py`](export_supplementary_tables_0814.py):
  recreates the CSV mirror directly from the formatted workbook.

The final composite images are authoritative for panel lettering and layout.
The source packages retained two pre-assembly letter sequences after panels
were removed; the final V4 mapping is documented in
[`figure_sources_0814/README.md`](figure_sources_0814/README.md).

## Machine Audit Layer

- `tables_v1.2.0/`: Table 0, Table 0W, Record_audit, Tables 1-2, and the
  historical S1-S24 machine-audit exports.
- `submission_audit_details_v1.2.0/`: record-level measurement-dispersion,
  mapping/dependence, and training-proximity audit tables.
- `kcat_benchmark_reorganized_tables_reviewed_v1.2.0.xlsx`: the earlier
  consolidated audit workbook retained for the existing 299-check pipeline;
  it is not the current submission-formatted supplementary workbook.
- `kcat_benchmark_audit_checks_v1.2.0.csv`: machine-readable results from the
  data-artifact audit.
- `paper_statistics_v1.2.0.json`: compact facts used by benchmark builders.
- `independent_cluster_inference_v1.2.0-r3.csv`: independent reconstruction of
  the 20 dependence-aware comparisons in the audit layer.

## Validation and Rebuilding

Run from the repository root:

```bash
# Validate the exact 0814/V4 workbook and figure snapshot.
python paper/generate_manuscript_figures.py

# Rerun the four data-backed plotting packages outside the repository.
python paper/generate_manuscript_figures.py \
  --rebuild-code-panels /tmp/kcat_0814_rebuild

# Recreate GitHub-previewable CSV files from the formatted workbook.
python paper/export_supplementary_tables_0814.py
```

The figure validator requires Pillow. Recalculation requires Matplotlib, NumPy,
and pandas; local typography can vary when Arial is unavailable. The workbook
export requires openpyxl.

The deeper benchmark build remains available through
`build_submission_audits.py`, `build_table0.py`, and `rebuild_paper_tables.py`.
Those scripts require the method-level result assets listed in
`zenodo_assets_manifest.csv`; the audit builder additionally requires RDKit,
SciPy, DIAMOND, and the public method training corpora documented in
`external_methods/METHOD_SOURCES.md`.
