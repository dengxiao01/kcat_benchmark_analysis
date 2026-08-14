# Figure 4 reproducible plotting package (v2)

This package contains the updated plotting data and reproducible Matplotlib code for the revised Figure 4.

## Update
Figure 4a has been redesigned as a heatmap.

## Files
- `Figure4_plotting_data_v2.xlsx`
- `data/figure4a_heatmap.csv`
- `data/figure4b_contrasts.csv`
- `data/figure4c_label_provenance.csv`
- `data/figure4d_training_proximity.csv`
- `data/figure4e_pretkcat_sensitivity.csv`
- `data/training_proximity_record_audit.csv`
- `plot_figure4_panels_v2.py`

The record-level training-proximity audit adds the maximum sequence-identity
and chemical-similarity values underlying the Figure 4c proximity classes. The
class labels and exact-overlap flags are also present in Supplementary Table
S22, but these four continuous neighbor metrics are not duplicated in the
formatted workbook.

## Run
```bash
python plot_figure4_panels_v2.py
```

## Output
The script writes each panel separately to `output/` in PNG (600 dpi), PDF and SVG.

## Note
Figure 4a uses grouped columns:
- Species
- Evidence source
- Substrate role

All values are MAE on the log10(kcat) scale.
