# Figures for the 0814/V4 Manuscript Snapshot

`Figure1.png` through `Figure4.png` are the four final high-resolution
composite figures from the author-approved 0814/V4 manuscript. The `panels/`
directory contains the 20 corresponding standalone raster panels under stable
Figure/panel filenames. Figure 2c-f are code-regenerated from the supplied 0814
data package so the standalone content and visible letters follow the final V4
order: c=CI-width ratio, d=Spearman correlation, e=fold agreement, and f=signed
error.

The final composite figures are authoritative for panel lettering and layout.
Several source plotting packages retained pre-assembly letters or were
reordered during final assembly; the mapping to final V4 letters is documented
in [`../figure_sources_0814/README.md`](../figure_sources_0814/README.md).

Full captions are stored in
[`../FIGURE_CAPTIONS.md`](../FIGURE_CAPTIONS.md), and plotting code plus compact
source data are under `../figure_sources_0814/`.

Validate the exact PNG snapshot with:

```bash
python paper/generate_manuscript_figures.py
```
