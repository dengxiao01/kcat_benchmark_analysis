# Changelog

## [1.2.0] - 2026-08-09

### Manuscript figure revision 1.2.0-r5 - 2026-08-14

- Published the manuscript-accepted Figure 1-4 layouts together with 17
  standalone panels and panel-level source-data CSV files.
- Made the public figure builder independent of the private manuscript DOCX by
  versioning the two high-resolution Figure 1 base panels.
- Removed panel-letter identifiers from standalone exports while retaining them
  in composite figures, restored standalone legends, added the zero name-only
  match label in Figure 1c, and retained the Figure 4a heatmap colorbar during
  standalone cropping.
- Preserved all benchmark rows, predictions, statistics, and conclusions; this
  revision changes figure presentation and reproducibility assets only.

### Artifact revision 1.2.0-r3 - 2026-08-09

- Independently recalculated all method-level metrics from Table 0W and all 20
  dependence-aware S19 comparisons from record-level errors; no numerical or
  significance-decision mismatches were found.
- Standardized 28 cited references and added direct citations for DeepGO-SE,
  the eQuilibrator compound database, and RDKit.
- Clarified evaluated-reactant structures, actual method sequence inputs,
  reaction-equation provenance, and long-sequence truncation policies.
- Expanded the manuscript-specific audit to 119 checks and the complete
  artifact audit to 299 checks; both pass in full.
- Published the manuscript-matched code, canonical data, Table 0/Table 0W,
  Record_audit, S1-S24, figures, workbook, and reproducibility assets as the
  v1.2.0 GitHub and Zenodo release.

### Artifact revision 1.2.0-r2 - 2026-08-07

- Rebased the reviewed manuscript on the author-supplied 0806 V3 section logic
  while refreshing every quantitative statement, manuscript table, figure, and
  method applicability description from the v1.2.0 artifacts.
- Added strictly common 1,047-record reaction-set equal-weight analyses to S17:
  796 sequence-substrate pair clusters and 761 label-assignment clusters for
  KcatNet, TurNuP, and PMAK.
- Added a versioned 0806 submission bundle and manuscript-specific audit; the
  final V3 audit passes 84/84 checks and the expanded project audit passes
  294/294 checks.

### Corrected benchmark construction and substrate roles

- Retained every model-encoded reactant before experimental substrate matching
  instead of using model metabolite IDs to prefilter cofactors.
- Expanded the complete-resource benchmark from 978 to 1,246 rows while
  preserving every previous entry ID. The 268 added rows are all E. coli; the
  465-row yeast entry set is unchanged.
- Corrected Quinate to the neutral PubChem CID 6508 structure. All 1,246 final
  SMILES are RDKit-parseable.
- Replaced the yeast-GEM `s_XXXX`-dependent role flag with a versioned registry
  using normalized names, external identifiers, and standardized parent
  connectivity. The final roles are 606 currency/cofactor, 34 carrier-linked
  variable, and 606 other-reactant rows.
- Defined the 1,246-row complete resource as the primary coverage analysis and
  the 778-row BRENDA substrate-supported set as the formal sensitivity
  analysis; 468 SABIO-RK-only participant-ambiguous rows are explicit.

### Reconstructed methods and statistical analyses

- Rebuilt PreTKcat under raw-public, exact-pair-excluded, and joint
  near-neighbor-excluded policies. The exact-excluded model is primary and
  removes 248 training rows before fitting 16,001 rows.
- Completed the DEKP public-data reconstruction with structures and predictions
  for all 1,246 rows after removing 230 exact-overlap training rows.
- Recomputed all 12 method outputs, coverage, metrics, paired tests, five
  cluster definitions, sensitivity analyses, and training-proximity audits.
- Updated the strict reaction-aware intersection to 1,047 rows, CatPred to
  1,156 rows, KinForm-L to 729 rows, and GO-HKP to 1,236 rows.

### Submission artifacts

- Rebuilt Table 0 as 14,952 method-record rows with actual sequence input,
  model-forward reaction equations, evaluated metabolite, experimental kcat,
  prediction status, and predicted kcat; Table 0W remains the 1,246-row matrix.
- Expanded Record_audit to 1,246 rows and 132 fields and added Supplementary
  Table S24 for the three PreTKcat reconstruction policies.
- Rebuilt Figures 1-4 and the reviewed XLSX and DOCX as versioned v1.2.0
  outputs without overwriting the archived author draft.
- Rewrote benchmark filtering, substrate-role, PreTKcat, DEKP, and Data and code
  availability sections; dated public-link verification was removed.
- Added `paper/audit_paper_artifacts_v1_2.py`; all 280 encoded checks pass.


## [1.1.0] - 2026-07-17

### Artifact revision 1.1.0-r7 - 2026-08-03

- Removed manuscript-build history from visible article prose, including the
  manuscript revision date and time-stamped public-link verification.
- Reframed Quinate and CKB text as stable mapping methods and an explicit
  reproducibility boundary rather than a correction log.
- Separated article-supplied benchmark data, GitHub-hosted project code, and
  Zenodo-hosted large inference assets in Data and code availability.
- Updated the canonical GitHub address to dxg-9527/kcat_benchmark_analysis
  after the public account redirect.
- Replaced remaining internal audit/local-run wording with scientific
  assessment and study-specific reconstruction terminology.
- Clarified the CatPred checkpoint ensemble/cache-key description and the
  provenance of DEKP protein structures.
- Preserved all 978 benchmark rows, predictions, metrics, tables, and figures;
  expanded the artifact audit to 525 scoped checks.

### Artifact revision 1.1.0-r6 - 2026-08-02

- Reorganized Supplementary Data Table 0 into 11,736 method-record rows with
  method, canonical sequence, materialized method sequence, sequence policy,
  model-forward reaction equations, evaluated metabolite, prediction status,
  experimental kcat, and predicted kcat.
- Preserved the former 978-row, 12-prediction-column matrix as Supplementary
  Data Table 0W (`Table0_wide`) for backward-compatible inspection.
- Retained explicit unscored rows for method-specific applicability gaps rather
  than dropping them from the long table.
- Rewrote `Final benchmark filtering and provenance fields` as separate
  filtering, reaction/metabolite-provenance, and Table 0 interpretation
  paragraphs.
- Corrected KcatNet long-sequence metadata to the actual first-500 plus
  last-500 rule used by its embedding code; predictions and metrics are
  unchanged.
- Increased the table schema version to 1.1 and expanded the artifact audit to
  490 scoped checks.

### Artifact revision 1.1.0-r5 - 2026-07-24

- Corrected the remaining Figure 3 section heading and short caption so the
  available-case panel is not described as a matched subset.
- Reported KcatNet-TurNuP and KcatNet-PMAK inference across protein, pair,
  reaction, reference, and label-assignment clustering rather than selecting
  only pair-cluster results.
- Added an `aggregation_rule` field to S17 and independently recomputed all 36
  method results for unique-pair, accession-unique-pair, and unique-label
  analyses using separate within-cluster medians and equal cluster weight.
- Retained the 978-row complete-resource analysis as primary and named the
  668-row SABIO ambiguity exclusion the substrate-supported sensitivity set.
- Repaired the author-name audit, added exact Figure 3 title and footer checks,
  and expanded the artifact audit to 402 scoped checks.
- Reduced displayed metric precision and wrapped long fields in S17 and
  S20-S23; full numeric values remain available in CSV and cell storage.
- Confirmed unauthenticated HTTP 200 access to the public GitHub repository and
  Zenodo record 21024684.

### Artifact revision 1.1.0-r4 - 2026-07-24

- Added a formal 668-row sensitivity scope excluding the 310 SABIO-RK-only
  participant-ambiguous labels, with method metrics and five cluster interval
  types for every predictor.
- Added label-assignment clusters to bootstrap intervals and paired tests so a
  weaker-evidence label shared across several sequence-substrate pairs is
  resampled as one unit.
- Renamed the mixed-coverage Figure 3 and S7 analyses as available-case scopes,
  displayed `n` in every heatmap cell, and reserved strict paired inference for
  shared identifiers.
- Unified S9, Figure 4, S17, and the manuscript on the 359-row union of model-ID
  and name-based currency/cofactor heuristics.
- Replaced the `remote` shorthand with
  `no_joint_neighbor_under_thresholds`, which does not exclude separate
  sequence-only or chemical-only training overlap.
- Added the 978-row `Record_audit` worksheet/CSV with source, matching,
  dependence, measurement, mapping, model-direction, substrate-role, and
  method-specific training-proximity fields.
- Added explicit S21 dispersion formulas, S22 before/excluded/retained mutation
  stages, and S23 model-reversibility counts while marking experimental assay
  direction as unverified.

### Artifact revision 1.1.0-r3 - 2026-07-24

- Replaced canonical-SMILES-only pair exclusion with sequence plus uncharged
  largest-fragment connectivity identity, then retrained PreTKcat and DEKP.
- Added matching-strength, duplicated-label, sensitivity, cluster-bootstrap,
  cluster-Wilcoxon, training-proximity, measurement-
  dispersion, mutation-status, and substrate-direction audits (S16-S23).
- Split Table 2 and Figure 2 into explicit checkpoint, retraining,
  reaction-aware, method-specific, and functional-assignment scopes.
- Revised the manuscript to distinguish sequence-associated labels from
  accession-matched measurements and to use dependence-aware statistical
  language.
- Added database snapshot hashes, updated resource references, corrected the
  author spelling to Qinghan Meng, and reduced the active footer to one PAGE
  field.

### Artifact revision 1.1.0-r2 - 2026-07-24

- Rebuilt the reviewed XLSX as plain OOXML worksheets without structured table
  objects or VML comments to improve Microsoft Excel compatibility.
- Shortened the visible Figure 2/3 label from DEKP-public-retrained to DEKP;
  internal method IDs, tables, and statistics are unchanged.
- Connected each pair of Figure 4b subgroup contrasts on one horizontal line
  and reserved non-data space for the panel b/c legends.

### Corrected

- Replaced the malformed Quinate value 192.167 with the neutral PubChem CID
  6508 SMILES and rebuilt all three canonical truth/benchmark layers.
- Added validation that rejects pure numeric values in SMILES fields and falls
  back from invalid CKB isomeric SMILES to a usable canonical SMILES.
- Refreshed all structure-dependent and reaction-aware method inputs,
  predictions, metrics, manuscript tables, and figures.

### Changed

- PreTKcat public-data retraining now excludes every public training row whose
  truncated model sequence and canonical substrate SMILES exactly match a
  benchmark pair. The v1.1.0 fit removes 25 source rows representing 22 unique
  pairs; these pairs correspond to 26 benchmark records.
- The primary sequence-plus-substrate benchmark is now 978/978 RDKit-valid
  structures. PMAK and TurNuP cover 782 rows, CatPred 914 rows, and KinForm-L
  564 rows.
- Paper artifacts now report benchmark version, data-freeze date, source-model
  versions and checksums, upstream method revisions, runtime versions, and
  deterministic statistical settings.
- Added a reproducible GO-HKP evaluator and a 158-check v1.1.0 paper/data audit.

## [1.0.0] - 2026-07-01

- Initial public benchmark release and unified evaluation of 12 methods.
