# CSV Mirror of the 0814/V4 Supplementary Tables

These CSV files are deterministic exports of the 24 worksheets in
[`../Supplementary_tables.xlsx`](../Supplementary_tables.xlsx). The XLSX file
is the authoritative submission artifact because it preserves the 0814 table
formatting, widths, wrapping, frozen panes, and number formats. The CSV mirror
is provided for GitHub preview, diffing, and programmatic reuse.

Regenerate the mirror with:

```bash
python paper/export_supplementary_tables_0814.py
```

## Table Order

| Table | Content |
| --- | --- |
| [S1](Table_S1.csv) | Hierarchical experimental matching levels by species |
| [S2](Table_S2.csv) | Experimental source-database support by species |
| [S3](Table_S3.csv) | Availability of experimental pH and temperature metadata |
| [S4](Table_S4.csv) | Benchmark-wide record and prediction matrix |
| [S5](Table_S5.csv) | Method inputs, implementation status, preprocessing, and caveats |
| [S6](Table_S6.csv) | Complete performance metrics on achieved evaluation sets |
| [S7](Table_S7.csv) | Conditional row-resampling intervals |
| [S8](Table_S8.csv) | Row-resampling rank stability on common evaluation sets |
| [S9](Table_S9.csv) | Cluster-bootstrap intervals under five dependence definitions |
| [S10](Table_S10.csv) | Strict 1,047-record reaction-common comparison |
| [S11](Table_S11.csv) | CatPred- and KinForm-L-scope available-case summaries |
| [S12](Table_S12.csv) | Secondary row-level Wilcoxon comparisons |
| [S13](Table_S13.csv) | Paired row-resampling MAE differences |
| [S14](Table_S14.csv) | Cluster-aggregated paired comparisons |
| [S15](Table_S15.csv) | Substrate, accession, role, pair, and label sensitivity analyses |
| [S16](Table_S16.csv) | Error stratification by species, source, and substrate role |
| [S17](Table_S17.csv) | Recurrent large-error benchmark records |
| [S18](Table_S18.csv) | Experimental-label dispersion |
| [S19](Table_S19.csv) | Mutation-status provenance audit |
| [S20](Table_S20.csv) | Performance by public training-corpus proximity |
| [S21](Table_S21.csv) | PreTKcat overlap-exclusion sensitivity |
| [S22](Table_S22.csv) | Record-level provenance and audit fields |
| [S23](Table_S23.csv) | Long-format method input, prediction, output, and error table |

The exact titles and worksheet names are also retained in [Index.csv](Index.csv).

The older `tables_v1.2.0/` directory remains available as the machine audit
layer used to build the benchmark. Its historical S1-S24 filenames are not the
numbering of the 0814/V4 submission workbook.
