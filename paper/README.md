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
- [`figures/`](figures/): four final high-resolution composite PNGs and 20
  matching standalone panel PNGs. Figure 2c-f were regenerated from the
  supplied data and plotting code so their filenames, visible letters, and
  content follow the final V4 order.
- [`FIGURE_CAPTIONS.md`](FIGURE_CAPTIONS.md): complete Figure 1-4 captions and
  interpretation boundaries from the V4 manuscript.
- [`figure_sources_0814/`](figure_sources_0814/): expanded plotting code and
  compact input data from the four author-supplied 0814 figure packages.
- [`generate_manuscript_figures.py`](generate_manuscript_figures.py): validates
  exact hashes/dimensions and can rerun every data-backed plotting package.

The final composite images are authoritative for panel lettering and layout.
Several source packages retained pre-assembly letter sequences or were
reordered during assembly; the final V4 mapping is documented in
[`figure_sources_0814/README.md`](figure_sources_0814/README.md).

## Record-Level Audit Coverage

- `Supplementary_tables.xlsx`, Table S22: the 1,246-row, 132-field source,
  structure, dependence, direction, role, measurement, and training-proximity
  audit.
- `Supplementary_tables.xlsx`, Table S23: the 14,952-row method-record input,
  prediction-status, output, and error table.
- `figure_sources_0814/Figure4/data/training_proximity_record_audit.csv`: the
  nonredundant record-level neighbor audit, including four continuous sequence
  and chemical similarity fields used to verify the proximity classes.
- `paper_statistics_v1.2.0.json`: compact facts used by benchmark builders.
- `independent_cluster_inference_v1.2.0-r3.csv`: independent reconstruction of
  the 20 dependence-aware comparisons in the audit layer.

The former `tables_v1.2.0/`, `supplementary_tables_0814/`,
`submission_audit_details_v1.2.0/`, earlier consolidated workbook, and stale
path-based audit report are intentionally not tracked. Their substantive
submission content is either present in the V4 workbook or, for the continuous
neighbor metrics, retained in the single Figure 4 audit file above.

## Validation and Rebuilding

Run from the repository root:

```bash
# Validate the exact 0814/V4 workbook and figure snapshot.
python paper/generate_manuscript_figures.py

# Rerun the four data-backed plotting packages outside the repository.
python paper/generate_manuscript_figures.py \
  --rebuild-code-panels /tmp/kcat_0814_rebuild
```

The figure validator requires Pillow. Recalculation requires Matplotlib, NumPy,
and pandas; local typography can vary when Arial is unavailable. The workbook
inspection and independent cluster check require openpyxl.

The deeper benchmark build remains available through
`build_submission_audits.py`, `build_table0.py`, and `rebuild_paper_tables.py`;
their intermediate table exports are generated locally and are not public
submission artifacts. These scripts require the method-level result assets
listed in `zenodo_assets_manifest.csv`; the audit builder additionally requires
RDKit, SciPy, DIAMOND, and the public method training corpora documented in
`external_methods/METHOD_SOURCES.md`.
