# Benchmark Context Tables

These English-language tables add dataset-characterization fields and summaries
to the three canonical benchmark files in `data/final/`.

| File | Content |
| --- | --- |
| `benchmark_ready_catpred_enriched_context.csv` | 1,246 benchmark records with EC classes, reaction identifiers, reaction names, pathway groups, direct Yeast9 KEGG pathways, and sequence/SMILES lengths |
| `benchmark_build_funnel.csv` | Species-level candidate and filtering counts |
| `benchmark_dataset_ec_class_summary.csv` | EC-class counts by species |
| `benchmark_dataset_top_reactions.csv` | Reaction-level record counts and identifiers |
| `benchmark_dataset_kegg_like_primary_group.csv` | KEGG-like primary-group counts by species |
| `benchmark_dataset_direct_yeast_kegg_pathways.csv` | Direct Yeast9 KEGG pathway counts and reaction coverage |

Rebuild the internal report and these public context tables after restoring the
released method and context assets:

```bash
python src/47_generate_dataset_method_context_report.py
```
