#!/usr/bin/env python3
"""Run TurNuP kcat inference on the unified benchmark reaction subset."""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch


BASE = Path(__file__).resolve().parent.parent
TURNUP_CODE = (
    BASE
    / "external_methods"
    / "AI_file"
    / "turnup"
    / "kcat_prediction_function-main"
    / "kcat_prediction_function-main"
    / "code"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run TurNuP predictions.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--esm-cache", type=Path, required=True)
    parser.add_argument("--device", default="cuda:0")
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, np.ndarray]:
    if path.exists():
        with path.open("rb") as handle:
            return pickle.load(handle)
    return {}


def save_cache(path: Path, cache: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(cache, handle, protocol=4)


def calculate_esm1b_vectors(enzyme_list: list[str], cache_path: Path, device_name: str) -> pd.DataFrame:
    import esm
    import enzyme_representations as er

    device = torch.device(device_name if torch.cuda.is_available() and device_name.startswith("cuda") else "cpu")
    df_enzyme = er.preprocess_enzymes([seq.upper() for seq in enzyme_list])
    cache = load_cache(cache_path)
    missing = [
        seq
        for seq in df_enzyme["amino acid sequence"].tolist()
        if seq not in cache and er.validate_enzyme(seq[:1022])
    ]

    if missing:
        model_location = TURNUP_CODE / "data" / "saved_models" / "ESM1b" / "esm1b_t33_650M_UR50S.pt"
        regression_location = TURNUP_CODE / "data" / "saved_models" / "ESM1b" / "esm1b_t33_650M_UR50S-contact-regression.pt"
        task_model = TURNUP_CODE / "data" / "saved_models" / "ESM1b" / "model_ESM_binary_A100_epoch_1_new_split.pkl"

        print(f"Loading ESM1b model on {device}. Missing sequences: {len(missing)}")
        model_data = torch.load(model_location, map_location="cpu")
        regression_data = torch.load(regression_location, map_location="cpu")
        model, alphabet = esm.pretrained.load_model_and_alphabet_core(
            "esm1b_t33_650M_UR50S",
            model_data,
            regression_data,
        )
        model_dict = torch.load(task_model, map_location="cpu")
        model_dict_v2 = {k.split("model.")[-1]: v for k, v in model_dict.items()}
        for key in [
            "module.fc1.weight",
            "module.fc1.bias",
            "module.fc2.weight",
            "module.fc2.bias",
            "module.fc3.weight",
            "module.fc3.bias",
        ]:
            model_dict_v2.pop(key, None)
        model.load_state_dict(model_dict_v2)
        model.eval().to(device)
        batch_converter = alphabet.get_batch_converter()

        for index, seq in enumerate(missing, start=1):
            model_input = seq[:1022]
            _, _, batch_tokens = batch_converter([(f"protein_{index}", model_input)])
            batch_tokens = batch_tokens.to(device)
            with torch.no_grad():
                result = model(batch_tokens, repr_layers=[33])
            cache[seq] = result["representations"][33][0][0].detach().cpu().numpy()
            if index % 25 == 0:
                save_cache(cache_path, cache)
                print(f"Cached {index}/{len(missing)} new ESM1b vectors")
        save_cache(cache_path, cache)
        del model
        if device.type == "cuda":
            torch.cuda.empty_cache()

    df_enzyme["enzyme rep"] = ""
    for ind in df_enzyme.index:
        seq = df_enzyme.loc[ind, "amino acid sequence"]
        if seq in cache:
            df_enzyme.at[ind, "enzyme rep"] = cache[seq]
    return df_enzyme


def is_missing_feature(value: object) -> bool:
    return isinstance(value, str) and value == ""


def robust_merging_reaction_and_enzyme_df(
    df_reaction: pd.DataFrame,
    df_enzyme: pd.DataFrame,
    df_kcat: pd.DataFrame,
) -> pd.DataFrame:
    df_kcat = df_kcat.copy()
    df_kcat["difference_fp"] = ""
    df_kcat["enzyme rep"] = ""
    df_kcat["complete"] = True

    for ind in df_kcat.index:
        reaction_match = df_reaction.loc[
            df_reaction["substrates"].eq(df_kcat.at[ind, "substrates"])
            & df_reaction["products"].eq(df_kcat.at[ind, "products"])
        ]
        enzyme_match = df_enzyme.loc[
            df_enzyme["amino acid sequence"].eq(df_kcat.at[ind, "enzyme"])
        ]
        diff_fp = reaction_match.iloc[0]["difference_fp"] if not reaction_match.empty else ""
        esm1b_rep = enzyme_match.iloc[0]["enzyme rep"] if not enzyme_match.empty else ""

        if is_missing_feature(diff_fp) or is_missing_feature(esm1b_rep):
            df_kcat.at[ind, "complete"] = False
        else:
            df_kcat.at[ind, "difference_fp"] = diff_fp
            df_kcat.at[ind, "enzyme rep"] = esm1b_rep
    return df_kcat


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.esm_cache.parent.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(TURNUP_CODE))
    import enzyme_representations as er
    import kcat_prediction_new as kp
    import metabolite_preprocessing as mp

    code_dir = str(TURNUP_CODE) + "/"
    er.CURRENT_DIR = code_dir
    mp.CURRENT_DIR = code_dir
    kp.CURRENT_DIR = code_dir
    kp.calcualte_esm1b_ts_vectors = lambda enzyme_list: calculate_esm1b_vectors(
        enzyme_list,
        args.esm_cache,
        args.device,
    )
    kp.merging_reaction_and_enzyme_df = robust_merging_reaction_and_enzyme_df

    inp = pd.read_csv(args.input)
    substrates = inp["substrates"].fillna("").astype(str).tolist()
    products = inp["products"].fillna("").astype(str).tolist()
    enzymes = inp["enzyme"].fillna("").astype(str).tolist()

    pred = kp.kcat_predicton_new(code_dir, substrates, products, enzymes)
    out = pd.concat([inp.reset_index(drop=True), pred.reset_index(drop=True)], axis=1)
    if "kcat [s^(-1)]" in out.columns:
        out["prediction_kcat"] = pd.to_numeric(out["kcat [s^(-1)]"], errors="coerce")
        out["prediction_log10"] = np.log10(out["prediction_kcat"].where(out["prediction_kcat"] > 0))
    out.to_csv(args.output, index=False)
    print(f"Wrote TurNuP predictions: {args.output}")


if __name__ == "__main__":
    main()
