#!/usr/bin/env python3
"""Run SELFprot kcat prediction on prepared benchmark inputs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
from transformers import AutoTokenizer


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "selfprot" / "selfprot_kcat_input_valid_smiles.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "selfprot" / "selfprot_kcat_input_output.csv"
DEFAULT_WEIGHTS = BASE / "external_methods" / "SELFprot" / "weights" / "models"
DEFAULT_PROT_TOKENIZER = BASE / "external_methods" / "SELFprot" / "weights" / "esm2_t12_35M_UR50D_tokenizer"

MODEL_FILES = [
    "chem_model.pt",
    "prot_model.pt",
    "joint_layer3x.pt",
    "kcat_head.pt",
    "position_encoding.pt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict log10(kcat) with SELFprot.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--weights", type=Path, default=DEFAULT_WEIGHTS)
    parser.add_argument("--protein-tokenizer", type=Path, default=DEFAULT_PROT_TOKENIZER)
    parser.add_argument("--folder-name", default="models", help="SELFprot folder label for output metadata.")
    parser.add_argument("--fold", type=int, default=1)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--limit", type=int, default=0)
    return parser.parse_args()


def resolve_device(requested: str) -> torch.device:
    if requested.startswith("cuda") and not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device(requested)


def load_module(path: Path, device: torch.device) -> Any:
    try:
        return torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=device)


def load_models(weights: Path, fold: int, device: torch.device) -> dict[str, Any]:
    fold_dir = weights / f"models_fold{fold}"
    if not fold_dir.exists():
        raise FileNotFoundError(f"SELFprot fold directory not found: {fold_dir}")
    missing = [name for name in MODEL_FILES if not (fold_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"{fold_dir} is missing model files: {', '.join(missing)}")
    models = {name: load_module(fold_dir / name, device) for name in MODEL_FILES}
    for module in models.values():
        if hasattr(module, "parameters"):
            for param in module.parameters():
                param.requires_grad = False
        if hasattr(module, "eval"):
            module.eval()
        if hasattr(module, "to"):
            module.to(device)
    return models


def load_input(path: Path, limit: int = 0) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"selfprot_row_id", "entry_id", "sequence", "smiles"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    df = df.copy()
    df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()
    df["smiles"] = df["smiles"].fillna("").astype(str).str.strip()
    df = df[(df["sequence"] != "") & (df["smiles"] != "")].reset_index(drop=True)
    if limit > 0:
        df = df.head(limit).copy()
    return df.reset_index(drop=True)


def tokenize_inputs(df: pd.DataFrame, sf_tokenizer_path: Path, protein_tokenizer_path: Path) -> TensorDataset:
    tokenizer_sf = AutoTokenizer.from_pretrained(str(sf_tokenizer_path))
    tokenizer_prot = AutoTokenizer.from_pretrained(str(protein_tokenizer_path))
    tokenizer_sf.model_max_length = 1024
    tokenizer_prot.model_max_length = 1024

    encoded_sf = tokenizer_sf(
        df["smiles"].astype(str).tolist(),
        padding="max_length",
        truncation=True,
        max_length=1024,
        return_tensors="pt",
    )
    encoded_prot = tokenizer_prot(
        df["sequence"].astype(str).tolist(),
        padding="max_length",
        truncation=True,
        max_length=1024,
        return_tensors="pt",
    )
    return TensorDataset(
        encoded_sf["input_ids"],
        encoded_sf["attention_mask"],
        encoded_prot["input_ids"],
        encoded_prot["attention_mask"],
    )


def predict(df: pd.DataFrame, models: dict[str, Any], weights: Path, protein_tokenizer: Path, device: torch.device, batch_size: int) -> tuple[np.ndarray, np.ndarray]:
    dataset = tokenize_inputs(df, weights / "sf_tokenizer", protein_tokenizer)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    chem_model = models["chem_model.pt"]
    prot_model = models["prot_model.pt"]
    joint_layer = models["joint_layer3x.pt"]
    kcat_head = models["kcat_head.pt"]
    position_encoding = models["position_encoding.pt"]

    positions = position_encoding(
        torch.arange(1, 2049, dtype=torch.long, device=device).view(1, 2048)
    ).to(device)
    pred_values: list[np.ndarray] = []
    pred_sd: list[np.ndarray] = []
    with torch.no_grad():
        for input_ids_sf, attention_mask_sf, input_ids_prot, attention_mask_prot in tqdm(loader, desc="SELFprot inference"):
            input_ids_sf = input_ids_sf.to(device=device, dtype=torch.long)
            attention_mask_sf = attention_mask_sf.to(device=device, dtype=torch.long)
            input_ids_prot = input_ids_prot.to(device=device, dtype=torch.long)
            attention_mask_prot = attention_mask_prot.to(device=device, dtype=torch.long)

            sf_predictions = chem_model(input_ids_sf, attention_mask=attention_mask_sf).last_hidden_state
            prot_predictions = prot_model(input_ids_prot, attention_mask=attention_mask_prot).last_hidden_state
            combined = torch.cat([prot_predictions, sf_predictions], dim=1) + positions
            joint_mask = torch.cat([attention_mask_prot, attention_mask_sf], dim=1).view(-1, 1, 1, 2048)
            mixture = joint_layer(combined, attention_mask=joint_mask)
            mean_mixture = torch.mean(mixture.last_hidden_state, dim=1)
            out = kcat_head(mean_mixture).detach().cpu().numpy()
            pred_values.append(out[:, 0])
            pred_sd.append(out[:, 1] if out.shape[1] > 1 else np.full(out.shape[0], np.nan))
    return np.concatenate(pred_values), np.concatenate(pred_sd)


def main() -> None:
    args = parse_args()
    device = resolve_device(args.device)
    print(f"SELFprot device: {device}", flush=True)
    df = load_input(args.input, args.limit)
    if df.empty:
        raise ValueError(f"No usable rows found in {args.input}")
    models = load_models(args.weights, args.fold, device)
    prediction_log10, prediction_sd = predict(df, models, args.weights, args.protein_tokenizer, device, args.batch_size)
    if len(prediction_log10) != len(df):
        raise RuntimeError(f"SELFprot returned {len(prediction_log10)} predictions for {len(df)} rows")

    out = df.copy()
    out["prediction_log10"] = prediction_log10.astype(float)
    out["prediction_kcat"] = np.power(10.0, out["prediction_log10"])
    out["selfprot_prediction_sd"] = prediction_sd.astype(float)
    out["prediction_column"] = "prediction_log10"
    out["selfprot_folder"] = args.folder_name
    out["selfprot_fold"] = args.fold

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote SELFprot predictions: {args.output}", flush=True)
    print(f"Rows: {len(out)}", flush=True)


if __name__ == "__main__":
    main()
