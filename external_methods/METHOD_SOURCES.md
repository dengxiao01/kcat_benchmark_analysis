# Third-party method sources

`external_methods/` is intentionally excluded from Git because it contains nested repositories, method-specific environments, large weights, and caches. Clone each upstream project into the exact local path below, then restore the corresponding Zenodo asset when required.

| method | upstream source | evaluated revision | expected local path | large asset | upstream license observed locally |
| --- | --- | --- | --- | --- | --- |
| CatPred | <https://github.com/maranasgroup/CatPred.git> | `8e72d324e9e6f7a9a24c3f8a720884c7c1740a9b` | `external_methods/CatPred` | `kcat_benchmark_catpred_kcat_assets.tar.gz`, `catpred_db.tar.gz` | MIT |
| CataPro | <https://github.com/zchwang/CataPro.git> | `cc89b2c81768665cf6fd76dfda607ce88691f601` | `external_methods/CataPro` | `kcat_benchmark_other_model_assets.tar` | MIT |
| DEKP | <https://github.com/wang-yi-zhen/DEKP.git> | `d2b8c1372b5c1855fd2de9aaadde19cf8cc7fa8d` | `external_methods/DEKP` | public-retrained model in `kcat_benchmark_other_model_assets.tar` | check upstream repository |
| DLKcat | <https://github.com/SysBioChalmers/DLKcat.git> | `7c15d0d4a7ac029f9d75564d9f2a93874aeaaec7` | `external_methods/DLKcat_official` | upstream model files | GPL-3.0 |
| GO-HKP | <https://github.com/tibbdc/GO-HKP.git> | `5d086a4ded543295250eb7db2a1ea4b1336e7f48` | `external_methods/GO-HKP` | upstream GO tree/data | MIT + Commons Clause 1.0 |
| KcatNet | <https://github.com/BioColLab/KcatNet.git> | `7d370f69f9d1bbed517655d23d4d80bd76594321` | `external_methods/KcatNet` | `kcat_benchmark_other_model_assets.tar` | GPL-3.0 |
| KinForm | <https://github.com/Digital-Metabolic-Twin-Centre/KinForm.git> | `f7a70eb1cd6723ba3a8d606432e522ea2b0fa9fd` | `external_methods/KinForm` | `kinform_results.tar.gz` | MIT |
| PMAK | <https://github.com/MrVincentCai/PMAK.git> | `1b1bea4580ef7bb908f893d3b13213a1486bbb98` | `external_methods/PMAK` | `kcat_benchmark_other_model_assets.tar` | check upstream repository |
| PreTKcat | <https://github.com/MrVincentCai/PreTKcat.git> | `b7bc0562a9b8555a201c5f6c72fc2a660dcdb76d` | `external_methods/PreTKcat` | `kcat_benchmark_other_model_assets.tar` | check upstream repository |
| SELFprot | <https://github.com/marltanwilson/SELFprot.git> | `880c2e8fd685ed0e8d574382439f7d7ca75cc9d0` | `external_methods/SELFprot` | `kcat_benchmark_other_model_assets.tar` | MIT |
| UniKP | <https://github.com/Luo-SynBioLab/UniKP> | `5cee5c4a64ba2daf59c63a5b5cbaa0cadf97ef26` | `external_methods/AI_file/UniKP` | `kcat_benchmark_other_model_assets.tar` | check upstream repository |
| TurNuP | <https://github.com/AlexanderKroll/kcat_prediction> | source archive supplied locally | `external_methods/AI_file/turnup/kcat_prediction_function-main/kcat_prediction_function-main` | `kcat_benchmark_turnup_kcat_assets.tar` | check upstream repository and ESM licenses |

The Zenodo bundles intentionally exclude generic foundation checkpoints that are already distributed by their maintainers:

| method | foundation model | official source | expected path |
| --- | --- | --- | --- |
| TurNuP | ESM-1b `esm1b_t33_650M_UR50S.pt` | <https://github.com/facebookresearch/esm> | `external_methods/AI_file/turnup/kcat_prediction_function-main/kcat_prediction_function-main/code/data/saved_models/ESM1b/esm1b_t33_650M_UR50S.pt` |

Example:

```bash
git clone https://github.com/maranasgroup/CatPred.git external_methods/CatPred
git -C external_methods/CatPred checkout 8e72d324e9e6f7a9a24c3f8a720884c7c1740a9b
python scripts/download_zenodo_assets.py --asset kcat_benchmark_catpred_kcat_assets.tar.gz --restore
```

The repository revision records what was evaluated; it does not imply that every upstream repository remains installable unchanged on newer Python/PyTorch versions. Always cite the original paper and comply with the license shipped by each upstream project/model.
