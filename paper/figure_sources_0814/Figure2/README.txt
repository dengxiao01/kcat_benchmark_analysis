Figure 2 reproducibility package

Data files
----------
panel_a.csv : benchmark coverage and MAE
panel_b.csv : MAE estimates and row-bootstrap 95% CIs
panel_c_ci_width_ratio.csv : cluster-bootstrap CI-width / row-bootstrap CI-width ratios
panel_d_spearman.csv : Spearman estimates and row-bootstrap 95% CIs
panel_e_fold_agreement.csv : within 2-fold / 10-fold fractions and 95% CIs
panel_f_signed_error.csv : mean signed error

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
The script saves Figure 2a-f separately as PNG, PDF and SVG in the final V4
manuscript order.

Source
------
Supplementary Tables S5-S9.

Panel c definition:
CI width ratio = cluster-bootstrap MAE CI width / row-bootstrap MAE CI width.

Version note
------------
The supplied 0814 package used the pre-assembly order c=Spearman,
d=fold agreement, e=signed error and f=CI-width ratio. The final V4 manuscript
reordered these panels to c=CI-width ratio, d=Spearman, e=fold agreement and
f=signed error. This public script and the filenames above follow the final V4
manuscript order.
