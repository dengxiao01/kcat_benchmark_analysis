# Figure Sources for the 0814/V4 Manuscript Snapshot

This directory is the expanded, Git-friendly form of the four plotting
packages supplied with the author-approved `paper/0814` manuscript snapshot.
The original ZIP containers are intentionally not tracked; their code, input
tables, and documentation are tracked as ordinary files instead.

## Final Images and Recalculation

The authoritative submitted composite raster images are in `../figures/`:

- four final V4 composite figures in `../figures/Figure1.png` through
  `Figure4.png`;
- 20 corresponding panel exports in `../figures/panels/`; Figure 2c-f are
  regenerated from the supplied package after normalizing its panel order to
  the final V4 manuscript;
- SHA-256 hashes and dimensions in `figure_asset_manifest.csv`.

The scripts under `Figure1/` through `Figure4/` recalculate the data-backed
panels from the small CSV/XLSX inputs supplied in the 0814 packages. They were
executed successfully during publication preparation. Their regenerated raster
appearance can vary slightly with the installed Matplotlib and Arial fonts, so
the frozen composite PNG files define the exact 0814 manuscript appearance.
The standalone Figure 2c-f files are the normalized reference rerender from the
same supplied data and code.

Run every package in a new directory with:

```bash
python paper/generate_manuscript_figures.py \
  --rebuild-code-panels /tmp/kcat_0814_rebuild
```

Running the command without `--rebuild-code-panels` validates every frozen
figure and the supplementary workbook against the manifest.

## Panel Mapping

The V4 manuscript removed intermediate package panels and also reordered four
Figure 2 panels during final assembly:

| Public V4 panel | Package-internal panel | Content |
| --- | --- | --- |
| Figure 1c | Figure 1d | Experimental matching evidence by species |
| Figure 1d | Figure 1e | Experimental source-database support |
| Figure 1e | Figure 1f | Experimental-condition metadata completeness |
| Figure 2c | Figure 2f | Dependence-aware CI-width-ratio heatmap |
| Figure 2d | Figure 2c | Spearman correlation with row-bootstrap CI |
| Figure 2e | Figure 2d | Within-fold agreement |
| Figure 2f | Figure 2e | Mean signed error |
| Figure 4c | Figure 4d | Public training-corpus proximity |
| Figure 4d | Figure 4e | PreTKcat exclusion sensitivity |

The public Figure 2 script, input filenames, and standalone outputs have been
normalized to the final V4 order. The table above records why their ordering
differs from the originally supplied package.

The Figure 4 package can also regenerate its package-internal panel c on label
provenance. That panel is retained as reproducibility material but is not part
of the final four-panel V4 Figure 4 or its caption. Figure 1a and Figure 1b are
author-designed raster artwork and do not have a plotting script in the 0814
package.

The complete final captions and statistical interpretation boundaries are in
[`../FIGURE_CAPTIONS.md`](../FIGURE_CAPTIONS.md).
