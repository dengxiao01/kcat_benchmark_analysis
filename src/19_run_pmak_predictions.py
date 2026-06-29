#!/usr/bin/env python3
"""Run PMAK kcat prediction on prepared benchmark inputs."""

from __future__ import annotations

import argparse
import gc
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rxnfp.transformer_fingerprints import RXNBERTFingerprintGenerator, get_default_model_and_tokenizer
from transformers import T5EncoderModel, T5Tokenizer


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "pmak" / "pmak_kcat_input.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "pmak" / "pmak_kcat_input_output.csv"
DEFAULT_FEATURE_CACHE = BASE / "data" / "final" / "pmak" / "pmak_feature_cache.pkl"
DEFAULT_PMAK_CODE = BASE / "external_methods" / "PMAK" / "code"
DEFAULT_MODEL_DIR = DEFAULT_PMAK_CODE / "save_model" / "CV"
DEFAULT_PROTT5 = BASE / "external_methods" / "CataPro" / "models" / "prot_t5_xl_uniref50"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict log10(kcat) with PMAK reaction-cold checkpoints.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--feature-cache", type=Path, default=DEFAULT_FEATURE_CACHE)
    parser.add_argument("--pmak-code", type=Path, default=DEFAULT_PMAK_CODE)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL_DIR)
    parser.add_argument("--prott5-model", type=Path, default=DEFAULT_PROTT5)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--max-checkpoints", type=int, default=0, help="Use only the first N checkpoints. 0 means all.")
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, dict[str, np.ndarray]]:
    if not path.exists():
        return {"sequence": {}, "reaction": {}}
    with path.open("rb") as handle:
        cache = pickle.load(handle)
    cache.setdefault("sequence", {})
    cache.setdefault("reaction", {})
    return cache


def save_cache(path: Path, cache: dict[str, dict[str, np.ndarray]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(cache, handle, protocol=pickle.HIGHEST_PROTOCOL)


def truncate_sequence(sequence: str) -> str:
    if len(sequence) > 1000:
        return sequence[:500] + sequence[-500:]
    return sequence


def spaced_sequence(sequence: str) -> str:
    sequence = truncate_sequence(str(sequence).strip())
    return " ".join(sequence)


def compute_sequence_features(sequences: list[str], model_path: Path, device: torch.device) -> dict[str, np.ndarray]:
    if not sequences:
        return {}
    tokenizer = T5Tokenizer.from_pretrained(str(model_path), do_lower_case=False)
    model = T5EncoderModel.from_pretrained(str(model_path))
    gc.collect()
    model = model.to(device)
    model.eval()

    features: dict[str, np.ndarray] = {}
    for index, sequence in enumerate(sequences, start=1):
        print(f"For PMAK sequence {index}/{len(sequences)}", flush=True)
        sequence_text = re.sub(r"[UZOB]", "X", spaced_sequence(sequence))
        encoded = tokenizer.batch_encode_plus([sequence_text], add_special_tokens=True, padding=True)
        input_ids = torch.tensor(encoded["input_ids"]).to(device)
        attention_mask = torch.tensor(encoded["attention_mask"]).to(device)
        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        seq_len = int((attention_mask[0] == 1).sum().item())
        features[sequence] = embedding[0, : seq_len - 1].detach().cpu().numpy().astype(np.float32)

    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return features


def compute_reaction_features(reactions: list[str]) -> dict[str, np.ndarray]:
    if not reactions:
        return {}
    model, tokenizer = get_default_model_and_tokenizer()
    generator = RXNBERTFingerprintGenerator(model, tokenizer)
    features: dict[str, np.ndarray] = {}
    for index, reaction_smiles in enumerate(reactions, start=1):
        print(f"For PMAK reaction {index}/{len(reactions)}", flush=True)
        features[reaction_smiles] = np.asarray(generator.convert(reaction_smiles), dtype=np.float32)
    return features


def ensure_features(
    df: pd.DataFrame,
    cache_path: Path,
    prott5_model: Path,
    device: torch.device,
) -> dict[str, dict[str, np.ndarray]]:
    cache = load_cache(cache_path)
    sequences = list(dict.fromkeys(df["sequence"].astype(str)))
    reactions = list(dict.fromkeys(df["reaction_smiles"].astype(str)))

    missing_sequences = [sequence for sequence in sequences if sequence not in cache["sequence"]]
    missing_reactions = [reaction for reaction in reactions if reaction not in cache["reaction"]]

    print(
        f"Feature cache: {len(sequences) - len(missing_sequences)}/{len(sequences)} sequences, "
        f"{len(reactions) - len(missing_reactions)}/{len(reactions)} reactions already cached",
        flush=True,
    )
    if missing_sequences:
        cache["sequence"].update(compute_sequence_features(missing_sequences, prott5_model, device))
        save_cache(cache_path, cache)
    if missing_reactions:
        cache["reaction"].update(compute_reaction_features(missing_reactions))
        save_cache(cache_path, cache)
    return cache


def load_pmak_model(code_dir: Path, checkpoint: Path, device: torch.device):
    sys.path.insert(0, str(code_dir))
    from Kcat_model import InteractPre  # noqa: PLC0415

    model = InteractPre().to(device)
    state = torch.load(checkpoint, map_location=device)
    model.load_state_dict(state, strict=False)
    model.eval()
    return model


def predict_with_checkpoint(model, df: pd.DataFrame, cache: dict[str, dict[str, np.ndarray]], device: torch.device) -> np.ndarray:
    preds: list[float] = []
    with torch.no_grad():
        for _, row in df.iterrows():
            reaction = torch.tensor(cache["reaction"][str(row["reaction_smiles"])], dtype=torch.float32).unsqueeze(0).unsqueeze(0)
            protein = torch.tensor(cache["sequence"][str(row["sequence"])], dtype=torch.float32).unsqueeze(0)
            reaction = reaction.to(device)
            protein = protein.to(device)
            pred = model(reaction, protein).detach().cpu().numpy().reshape(-1)[0]
            preds.append(float(pred))
    return np.asarray(preds, dtype=np.float32)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    required = {"pmak_row_id", "entry_id", "sequence", "reaction_smiles"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{args.input} is missing required columns: {', '.join(missing)}")
    df = df.copy()

    device = torch.device(args.device if torch.cuda.is_available() or not args.device.startswith("cuda") else "cpu")
    print(f"PMAK device: {device}", flush=True)
    cache = ensure_features(df, args.feature_cache, args.prott5_model, device)

    checkpoints = sorted(args.model_dir.glob("Fold_*_reaction_cold.pth"))
    if args.max_checkpoints > 0:
        checkpoints = checkpoints[: args.max_checkpoints]
    if not checkpoints:
        raise FileNotFoundError(f"No PMAK reaction-cold checkpoints found in {args.model_dir}")

    fold_preds = []
    for checkpoint in checkpoints:
        print(f"Predicting with {checkpoint.name}", flush=True)
        model = load_pmak_model(args.pmak_code, checkpoint, device)
        fold_preds.append(predict_with_checkpoint(model, df, cache, device))
        del model
        if device.type == "cuda":
            torch.cuda.empty_cache()

    pred_matrix = np.vstack(fold_preds).T
    df["prediction_log10"] = pred_matrix.mean(axis=1)
    df["prediction_log10_std"] = pred_matrix.std(axis=1)
    df["prediction_column"] = "prediction_log10"
    for idx, checkpoint in enumerate(checkpoints):
        df[f"prediction_log10_{checkpoint.stem}"] = pred_matrix[:, idx]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Wrote PMAK predictions: {args.output}")
    print(f"Rows: {len(df)}; checkpoints: {len(checkpoints)}")


if __name__ == "__main__":
    main()
