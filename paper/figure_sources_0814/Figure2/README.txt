Figure 2 reproducibility package

Data files
----------
panel_a.csv : benchmark coverage and MAE
panel_b.csv : MAE estimates and row-bootstrap 95% CIs
panel_c.csv : Spearman estimates and row-bootstrap 95% CIs
panel_d.csv : within 2-fold / 10-fold fractions and 95% CIs
panel_e.csv : mean signed error
panel_f_ratio.csv : cluster-bootstrap CI-width / row-bootstrap CI-width ratios

Code
----
Figure2_plot.py

Run
---
python Figure2_plot.py

Dependencies
------------
numpy
matplotlib

Output
------
The script saves Figure 2a-f separately as PNG, PDF and SVG.

Source
------
Supplementary Tables S5-S9.

Panel f definition:
CI width ratio = cluster-bootstrap MAE CI width / row-bootstrap MAE CI width.
