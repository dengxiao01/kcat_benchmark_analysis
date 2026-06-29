# Third-party data and model notices

This project integrates data, software, and trained model assets created by other groups. The repository MIT license applies only to original benchmark code and documentation. It does not replace any upstream license.

## Experimental and annotation data

- **BRENDA:** BRENDA states that its copyrightable data are available under CC BY 4.0. This benchmark uses turnover-number records after organism/EC/substrate matching and aggregation. Source: <https://www.brenda-enzymes.org/license.php>.
- **SABIO-RK:** kcat records were queried through the SABIO-RK web service. Cite SABIO-RK and follow the current database terms: <https://sabio.h-its.org/>.
- **UniProt:** protein sequences and the yeast GO mapping were obtained from UniProt REST endpoints. Cite the UniProt release used for a final publication: <https://www.uniprot.org/>.
- **PubChem:** missing compound structures were queried through PUG REST: <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>.
- **Metabolic models:** `eciML1515.json` and `yeast-GEM.xml` retain the provenance and license of their respective model projects.

The benchmark CSV files are derived research outputs. Reusers should cite this benchmark and the source databases represented by the `source_database` field.

## Prediction methods and weights

The exact upstream repositories, evaluated revisions, expected local paths, and observed licenses are listed in [`external_methods/METHOD_SOURCES.md`](external_methods/METHOD_SOURCES.md). Model weights and bundled foundation models remain subject to their original licenses and acceptable-use terms. In particular, bundles may contain components from ProtT5, Uni-Mol, ESM, RXNFP, or other upstream projects.

Redistribution in the Zenodo companion record is for research reproducibility and does not assert ownership of third-party assets. If an upstream license or author request conflicts with the companion bundle, the upstream terms control and the affected file should be obtained directly from its official source.

