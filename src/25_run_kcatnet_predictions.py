#!/usr/bin/env python3
"""Run KcatNet kcat prediction on prepared benchmark inputs."""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "kcatnet" / "kcatnet_kcat_input_valid_smiles.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "kcatnet" / "kcatnet_kcat_input_output.csv"
DEFAULT_PROTEIN_CACHE = BASE / "data" / "final" / "kcatnet" / "kcatnet_protein_cache.pkl"
DEFAULT_LIGAND_CACHE = BASE / "data" / "final" / "kcatnet" / "kcatnet_ligand_cache.pkl"
DEFAULT_KCATNET_ROOT = BASE / "external_methods" / "KcatNet"
DEFAULT_PROTT5 = BASE / "external_methods" / "CataPro" / "models" / "prot_t5_xl_uniref50"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict log10(kcat) with KcatNet.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--protein-cache", type=Path, default=DEFAULT_PROTEIN_CACHE)
    parser.add_argument("--ligand-cache", type=Path, default=DEFAULT_LIGAND_CACHE)
    parser.add_argument("--kcatnet-root", type=Path, default=DEFAULT_KCATNET_ROOT)
    parser.add_argument("--prott5-model", type=Path, default=DEFAULT_PROTT5)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0, help="Predict only the first N rows. 0 means all rows.")
    return parser.parse_args()


def resolve_device(requested: str) -> torch.device:
    if requested.startswith("cuda") and not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device(requested)


def load_cache(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("rb") as handle:
        return pickle.load(handle)


def to_cpu(value: Any) -> Any:
    if torch.is_tensor(value):
        return value.detach().cpu()
    if isinstance(value, dict):
        return {key: to_cpu(item) for key, item in value.items()}
    if isinstance(value, list):
        return [to_cpu(item) for item in value]
    if isinstance(value, tuple):
        return tuple(to_cpu(item) for item in value)
    return value


def save_cache(path: Path, cache: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(to_cpu(cache), handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_input(path: Path, limit: int = 0) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"kcatnet_row_id", "entry_id", "Pro_seq", "Smile"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    df = df.copy()
    df["Pro_seq"] = df["Pro_seq"].fillna("").astype(str).str.strip()
    df["Smile"] = df["Smile"].fillna("").astype(str).str.strip()
    df = df[(df["Pro_seq"] != "") & (df["Smile"] != "")].reset_index(drop=True)
    if limit > 0:
        df = df.head(limit).copy()
    return df.reset_index(drop=True)


def prepare_imports(root: Path, prott5_model: Path, device: torch.device) -> None:
    os.environ.setdefault("KCATNET_PROTT5_DIR", str(prott5_model))
    os.environ["KCATNET_DEVICE"] = str(device)
    os.environ.setdefault("TORCH_HOME", str(BASE / "external_methods" / "torch_cache"))
    for local_pydeps in [
        BASE / "external_methods" / "kcatnet_scatter_src",
        BASE / "external_methods" / "catapro_pydeps",
    ]:
        if local_pydeps.exists():
            sys.path.insert(0, str(local_pydeps))
    sys.path.insert(0, str(root))
    os.chdir(root)


def ensure_features(
    df: pd.DataFrame,
    protein_cache_path: Path,
    ligand_cache_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    import __main__  # noqa: PLC0415
    from utils.build_vocab import Vocab, WordVocab  # noqa: PLC0415
    from utils.ligand_init import ligand_init  # noqa: PLC0415
    from utils.protein_init import protein_init  # noqa: PLC0415

    __main__.Vocab = Vocab
    __main__.WordVocab = WordVocab

    protein_cache = load_cache(protein_cache_path)
    ligand_cache = load_cache(ligand_cache_path)

    sequences = list(dict.fromkeys(df["Pro_seq"].astype(str)))
    smiles = list(dict.fromkeys(df["Smile"].astype(str)))
    missing_sequences = [sequence for sequence in sequences if sequence not in protein_cache]
    missing_smiles = [smile for smile in smiles if smile not in ligand_cache]

    print(
        f"KcatNet feature cache: {len(sequences) - len(missing_sequences)}/{len(sequences)} sequences, "
        f"{len(smiles) - len(missing_smiles)}/{len(smiles)} SMILES already cached",
        flush=True,
    )
    if missing_sequences:
        protein_cache.update(protein_init(missing_sequences))
        save_cache(protein_cache_path, protein_cache)
        print(f"Updated protein cache: {protein_cache_path}", flush=True)
    if missing_smiles:
        ligand_cache.update(ligand_init(missing_smiles))
        save_cache(ligand_cache_path, ligand_cache)
        print(f"Updated ligand cache: {ligand_cache_path}", flush=True)

    return (
        {sequence: protein_cache[sequence] for sequence in sequences},
        {smile: ligand_cache[smile] for smile in smiles},
    )


def build_model(root: Path, device: torch.device):
    from models.model_kcat import KcatNet  # noqa: PLC0415

    with (root / "config_KcatNet.json").open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    degree_dict = torch.load(root / "Dataset" / "degree.pt", map_location="cpu")
    prot_deg = degree_dict["protein_deg"]
    params = config["params"]
    model = KcatNet(
        prot_deg,
        mol_in_channels=params["mol_in_channels"],
        prot_in_channels=params["prot_in_channels"],
        prot_evo_channels=params["prot_evo_channels"],
        hidden_channels=params["hidden_channels"],
        pre_layers=params["pre_layers"],
        post_layers=params["post_layers"],
        aggregators=params["aggregators"],
        scalers=params["scalers"],
        total_layer=params["total_layer"],
        K=params["K"],
        heads=params["heads"],
        dropout=params["dropout"],
        dropout_attn_score=params["dropout_attn_score"],
        device=device,
    ).to(device)
    state = torch.load(root / "RESULT" / "model_KcatNet.pt", map_location=device)
    model.load_state_dict(state)
    model.eval()
    return model


def predict(df: pd.DataFrame, protein_dict: dict[str, Any], ligand_dict: dict[str, Any], root: Path, device: torch.device, batch_size: int) -> np.ndarray:
    from torch_geometric.loader import DataLoader  # noqa: PLC0415
    from utils.Kcat_Dataset import EnzMolDataset  # noqa: PLC0415
    from utils.trainer import pred  # noqa: PLC0415

    dataset = EnzMolDataset(df.reset_index(drop=True), ligand_dict, protein_dict)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, follow_batch=["mol_x", "prot_node_esm"])
    model = build_model(root, device)
    values = pred(model, loader, device=device)
    return np.asarray(values, dtype=np.float64)


def main() -> None:
    args = parse_args()
    device = resolve_device(args.device)
    print(f"KcatNet device: {device}", flush=True)
    prepare_imports(args.kcatnet_root, args.prott5_model, device)
    df = load_input(args.input, args.limit)
    if df.empty:
        raise ValueError(f"No usable rows found in {args.input}")

    protein_dict, ligand_dict = ensure_features(df, args.protein_cache, args.ligand_cache)
    pred_log10 = predict(df, protein_dict, ligand_dict, args.kcatnet_root, device, args.batch_size)
    if len(pred_log10) != len(df):
        raise RuntimeError(f"KcatNet returned {len(pred_log10)} predictions for {len(df)} rows")

    out = df.copy()
    out["prediction_log10"] = pred_log10
    out["prediction_kcat"] = np.power(10.0, pred_log10)
    out["Predicted Kcats"] = out["prediction_kcat"]
    out["prediction_column"] = "prediction_log10"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote KcatNet predictions: {args.output}", flush=True)
    print(f"Rows: {len(out)}", flush=True)


if __name__ == "__main__":
    main()
