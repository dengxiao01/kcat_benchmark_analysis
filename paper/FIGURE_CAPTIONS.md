# Figure Captions for the 0814/V4 Manuscript Snapshot

These captions reproduce the final Figure 1-4 descriptions in the 0814/V4
manuscript and correspond to the composite PNG files under `paper/figures/`.

## Figure 1

**Construction, composition and experimental provenance of the kcat
benchmark.**

**(a)** Benchmark construction workflow. Enzyme-reaction associations were
parsed from the genome-scale metabolic models eciML1515 (*E. coli*) and Yeast9
(*S. cerevisiae*). All model-encoded reactants were retained during candidate
generation, yielding 16,613 enzyme-reactant candidates (8,797 for *E. coli* and
7,816 for *S. cerevisiae*). Protein sequences and reactant structures were
mapped before hierarchical matching to experimental kcat records from BRENDA
and SABIO-RK. Positive experimental matches were subsequently filtered by kcat
validity and unit convertibility, availability of a single protein sequence,
and an RDKit-parseable reactant SMILES, resulting in the final benchmark used
for unified predictor evaluation. **(b)** Benchmark scale and biochemical
breadth. The final dataset comprises 1,246 records, including 781 from *E.
coli* and 465 from *S. cerevisiae*, spanning 773 model reactions, 518 UniProt
accessions, 390 EC annotation strings and 288 unique reactant structures.
Experimental kcat values range from 1.67 x 10^-4 to 5.7 x 10^5 s^-1, with a
median log10(kcat) of 1.279 and a dynamic range exceeding nine orders of
magnitude. **(c)** Experimental matching evidence stratified by species.
Records were classified according to the most specific evidence supporting the
experimental assignment: species + EC + UniProt + substrate identifier,
species + EC + substrate identifier, or species + EC + normalized substrate
name. **(d)** Experimental source-database support by species, showing records
supported exclusively by BRENDA, exclusively by SABIO-RK, or by both databases.
**(e)** Completeness of experimental-condition metadata. Bars indicate the
proportions of benchmark records with reported pH or temperature values for
each species; counts are shown together with the corresponding percentages.

## Figure 2

**Coverage-aware evaluation of predictive accuracy, ranking performance,
calibration and statistical uncertainty across current kcat prediction
methods.**

**(a)** Relationship between benchmark coverage and mean absolute error (MAE)
on the log10(kcat) scale for the 12 evaluated methods. Marker colors and shapes
indicate inference regimes, including released sequence-substrate checkpoint
models, temperature-conditioned public retraining, reaction-aware models,
method-specific checkpoint subsets, structure-aware public retraining and the
functional-assignment baseline. **(b)** MAE estimates with 95% confidence
intervals obtained by row-level nonparametric bootstrap resampling (2,000
replicates). The number of evaluated benchmark records is indicated for each
method. **(c)** Sensitivity of MAE uncertainty to dependence among benchmark
records. Heatmap values show the ratio of the cluster-bootstrap 95% confidence-
interval width to the corresponding row-bootstrap interval width for five
clustering schemes: protein, enzyme-substrate pair, reaction, experimental
reference and experimental-label assignment. Values greater than 1 indicate
broader uncertainty after resampling related records as intact clusters.
Cluster-bootstrap estimates were calculated on each method's achieved
evaluation set using 2,000 replicates. Because methods differ in inference
requirements and applicability, performance estimates obtained on restricted
evaluation subsets should be interpreted within their corresponding coverage
domains rather than as direct full-benchmark rankings. **(d)** Spearman rank
correlations between predicted and experimental log10(kcat), with 95% row-
bootstrap confidence intervals. **(e)** Fractions of predictions falling within
two-fold and ten-fold of the corresponding experimental kcat values, with 95%
row-bootstrap confidence intervals; the dashed vertical line marks 50%
agreement. **(f)** Mean signed error on the log10(kcat) scale, showing
directional calibration bias. Positive values indicate systematic
overestimation and negative values indicate underestimation.

## Figure 3

**Matched-scope and dependence-aware comparisons reveal instability in the
ordering of leading kcat predictors.**

**(a)** Mean absolute error (MAE) for methods evaluated on the strict
1,047-record reaction-common set, allowing direct comparison on identical
benchmark records. Methods are ordered by MAE. **(b)** Available-case
performance within the CatPred-accessible (1,156 records) and
KinForm-L-accessible (729 records) scopes. Filled symbols denote methods
evaluated on all records in the indicated scope, whereas open symbols denote
methods with partial coverage; the number of evaluated records is shown for
each method. These panels therefore summarize performance within each
accessible scope but should not be interpreted as strict rankings between
methods evaluated on different subsets. **(c)** Paired row-level bootstrap
estimates of MAE differences among KcatNet, TurNuP and PMAK on the 1,047-record
common set. Points indicate observed MAE differences (method a - method b) and
horizontal lines indicate percentile 95% confidence intervals from 2,000
paired bootstrap resamples; the dashed vertical line marks zero difference.
**(d)** Sensitivity of pairwise comparisons to the definition of the independent
biological or experimental unit. For each method pair, absolute errors were
aggregated within protein, standardized sequence-substrate pair, reaction,
literature-reference or experimental-label-assignment clusters and compared
using paired Wilcoxon tests. Cell color represents the difference in cluster-
mean absolute error (method a - method b), with blue indicating lower error for
method a and orange indicating lower error for method b. Values are globally
Benjamini-Hochberg-adjusted q values; asterisks indicate q < 0.05. **(e)** MAE
for KcatNet, TurNuP and PMAK after changing the unit assigned equal weight: the
original 1,047 benchmark rows, 796 unique standardized sequence-substrate
pairs, or 761 unique experimental-label assignments. For the latter two
analyses, observed and predicted log10(kcat) values were separately aggregated
by the median within each cluster before MAE calculation. Together, the panels
show that the small differences among the leading predictors depend on the
shared evaluation scope, dependence structure and weighting unit.

## Figure 4

**Biological context, experimental provenance and training-corpus proximity
shape prediction error.**

**(a)** Mean absolute error (MAE, log10(kcat)) for six representative predictors
stratified by species, experimental source and substrate role. Species
comparisons include *E. coli* and *S. cerevisiae*. Experimental-source groups
comprise BRENDA-only, SABIO-RK-only and records supported by both databases.
Substrate-role comparisons are restricted to records for which the
experimental substrate was supported by the matching evidence and include
other reactants and currency/cofactor-like substrates. Values are calculated
within each method's native applicability domain and therefore represent
within-method error patterns rather than coverage-matched rankings. **(b)**
Direction and magnitude of source- and substrate-role-associated MAE contrasts.
Source contrast is defined as MAE (BRENDA only) - MAE (SABIO-RK only), and
substrate-role contrast as MAE (other reactant) - MAE (currency/cofactor).
Positive and negative values indicate the direction of the corresponding error
difference. Carrier-linked substrates are omitted from the main contrast
because of their small sample size. **(c)** Prediction MAE stratified by
proximity to auditable public training corpora. Exact denotes an exact
standardized sequence-parent pair; near denotes a joint sequence/chemical
neighbor meeting the predefined sequence-identity and chemical-similarity
thresholds; none denotes the absence of such a joint neighbor. The auditable
corpus does not have the same relationship to the evaluated model for every
method - for example, it may represent a released fit split, a union of
training folds, a source corpus or a production corpus - so these comparisons
quantify public-corpus proximity and should not be interpreted uniformly as
evidence of training-set leakage. **(d)** Sensitivity of PreTKcat to
progressively stricter exclusion of benchmark-proximal training records.
Performance is shown for the unfiltered public reconstruction, the exact-pair-
excluded training set used for the primary benchmark result, and a more
stringent joint-near-neighbor-excluded training set. Numbers denote fitted
training-set sizes; MAE and Spearman correlation are evaluated on the same
1,246-record benchmark.
