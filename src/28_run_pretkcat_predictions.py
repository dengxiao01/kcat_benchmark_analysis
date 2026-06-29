#!/usr/bin/env python3
"""Run PreTKcat-style kcat prediction on prepared benchmark inputs.

The public PreTKcat repository contains feature extraction code, the public
training dataset, and instructions to download MolGNet/ProtT5, but it does not
publish a fitted ExtraTrees kcat regressor. This wrapper therefore trains an
ExtraTreesRegressor on the public PreTKcat kcat training table and then predicts
the benchmark rows with the same feature recipe:

    MolGNet substrate vector + ProtT5 sequence vector + normalized temperature.
"""

from __future__ import annotations

import argparse
import gc
import os
import pickle
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit import RDLogger
from sklearn.ensemble import ExtraTreesRegressor
from transformers import T5EncoderModel, T5Tokenizer


RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_input_valid_smiles.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_input_output.csv"
DEFAULT_TRAIN = BASE / "external_methods" / "PreTKcat" / "datasets" / "DLTKcat_data" / "kcat_merge_DLTKcat.csv"
DEFAULT_PRETKCAT_ROOT = BASE / "external_methods" / "PreTKcat"
DEFAULT_MOLGNET = DEFAULT_PRETKCAT_ROOT / "MolGNet.pt"
DEFAULT_PROTT5 = BASE / "external_methods" / "CataPro" / "models" / "prot_t5_xl_uniref50"
DEFAULT_FEATURE_CACHE = BASE / "data" / "final" / "pretkcat" / "pretkcat_feature_cache.pkl"
DEFAULT_MODEL_CACHE = BASE / "data" / "final" / "pretkcat" / "pretkcat_extratrees_model.pkl"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict log10(kcat) with a PreTKcat-style public-data fit.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--train-data", type=Path, default=DEFAULT_TRAIN)
    parser.add_argument("--pretkcat-root", type=Path, default=DEFAULT_PRETKCAT_ROOT)
    parser.add_argument("--molgnet-model", type=Path, default=DEFAULT_MOLGNET)
    parser.add_argument("--prott5-model", type=Path, default=DEFAULT_PROTT5)
    parser.add_argument("--feature-cache", type=Path, default=DEFAULT_FEATURE_CACHE)
    parser.add_argument("--model-cache", type=Path, default=DEFAULT_MODEL_CACHE)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--limit", type=int, default=0, help="Predict only the first N benchmark rows. 0 means all.")
    parser.add_argument("--sample-train-size", type=int, default=0, help="Use only first N training rows for smoke tests.")
    parser.add_argument("--n-estimators", type=int, default=100)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--force-retrain", action="store_true")
    return parser.parse_args()


def resolve_device(requested: str) -> torch.device:
    if requested.startswith("cuda") and not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device(requested)


def prepare_imports(pretkcat_root: Path) -> None:
    for path in [
        BASE / "external_methods" / "kcatnet_scatter_src",
        BASE / "external_methods" / "catapro_pydeps",
        pretkcat_root / "Pretrained_Model",
    ]:
        if path.exists():
            sys.path.insert(0, str(path))
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")


def load_cache(path: Path) -> dict[str, dict[str, np.ndarray]]:
    if not path.exists():
        return {"sequence": {}, "smiles": {}}
    with path.open("rb") as handle:
        cache = pickle.load(handle)
    cache.setdefault("sequence", {})
    cache.setdefault("smiles", {})
    return cache


def save_cache(path: Path, cache: dict[str, dict[str, np.ndarray]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(cache, handle, protocol=pickle.HIGHEST_PROTOCOL)


def canonical_smiles(smiles: object) -> str:
    text = str(smiles).strip()
    if not text or text == "nan":
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def pretkcat_model_sequence(sequence: object) -> str:
    text = str(sequence).strip()
    if len(text) > 1000:
        return text[:500] + text[-500:]
    return text


def spaced_sequence(sequence: str) -> str:
    return " ".join(pretkcat_model_sequence(sequence))


def load_input(path: Path, limit: int = 0) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"pretkcat_row_id", "entry_id", "sequence", "smiles", "temp_k_norm", "inv_temp_norm"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    df = df.copy()
    df["sequence"] = df["sequence"].fillna("").astype(str).str.strip()
    df["smiles"] = df["smiles"].fillna("").astype(str).str.strip()
    df["temp_k_norm"] = pd.to_numeric(df["temp_k_norm"], errors="coerce")
    df["inv_temp_norm"] = pd.to_numeric(df["inv_temp_norm"], errors="coerce")
    df = df[(df["sequence"] != "") & (df["smiles"] != "")].copy()
    df = df.dropna(subset=["temp_k_norm", "inv_temp_norm"]).reset_index(drop=True)
    df["pretkcat_model_sequence"] = df["sequence"].map(pretkcat_model_sequence)
    if limit > 0:
        df = df.head(limit).copy()
    return df.reset_index(drop=True)


def load_training(path: Path, sample_size: int = 0) -> pd.DataFrame:
    train = pd.read_csv(path)
    required = {"kcat", "smiles", "seq", "Temp_K_norm", "Inv_Temp_norm"}
    missing = sorted(required.difference(train.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    train = train.copy()
    train["kcat"] = pd.to_numeric(train["kcat"], errors="coerce")
    train["Temp_K_norm"] = pd.to_numeric(train["Temp_K_norm"], errors="coerce")
    train["Inv_Temp_norm"] = pd.to_numeric(train["Inv_Temp_norm"], errors="coerce")
    train["seq"] = train["seq"].fillna("").astype(str).str.strip()
    train["smiles"] = train["smiles"].fillna("").astype(str).str.strip()
    train = train[
        (train["kcat"] > 0)
        & (train["seq"] != "")
        & (train["smiles"] != "")
        & (~train["smiles"].str.contains(".", regex=False))
    ].copy()
    train = train.dropna(subset=["Temp_K_norm", "Inv_Temp_norm"]).copy()
    train["canonical_smiles"] = train["smiles"].map(canonical_smiles)
    train = train[train["canonical_smiles"] != ""].copy()
    train["pretkcat_model_sequence"] = train["seq"].map(pretkcat_model_sequence)
    train["label_log10"] = np.log10(train["kcat"].astype(float))
    if sample_size > 0:
        train = train.head(sample_size).copy()
    return train.reset_index(drop=True)


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
        if index == 1 or index % 25 == 0 or index == len(sequences):
            print(f"For PreTKcat sequence {index}/{len(sequences)}", flush=True)
        sequence_text = re.sub(r"[UZOB]", "X", spaced_sequence(sequence))
        encoded = tokenizer.batch_encode_plus([sequence_text], add_special_tokens=True, padding=True)
        input_ids = torch.tensor(encoded["input_ids"], device=device)
        attention_mask = torch.tensor(encoded["attention_mask"], device=device)
        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        seq_len = int((attention_mask[0] == 1).sum().item())
        pooled = embedding[0, : seq_len - 1].mean(dim=0)
        features[sequence] = pooled.detach().cpu().numpy().astype(np.float32)

    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return features


def load_molgnet(model_path: Path, device: torch.device):
    from MPG_util.graph_bert import MolGT  # noqa: PLC0415

    if not model_path.exists():
        raise FileNotFoundError(
            f"MolGNet checkpoint not found: {model_path}. "
            "Download MPG MolGNet.pt and place it at this path."
        )
    model = MolGT(num_layer=5, emb_dim=768, heads=12, num_message_passing=3, drop_ratio=0.5)
    state = torch.load(model_path, map_location=device)
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    model.load_state_dict(state)
    model = model.to(device)
    model.eval()
    return model


def compute_smiles_features(smiles_list: list[str], model_path: Path, device: torch.device) -> dict[str, np.ndarray]:
    if not smiles_list:
        return {}
    from MPG_util.mol2graph import mol_to_graph_data_dic  # noqa: PLC0415
    from torch_geometric.nn import global_add_pool  # noqa: PLC0415

    model = load_molgnet(model_path, device)
    features: dict[str, np.ndarray] = {}
    with torch.no_grad():
        for index, smiles in enumerate(smiles_list, start=1):
            if index == 1 or index % 50 == 0 or index == len(smiles_list):
                print(f"For PreTKcat SMILES {index}/{len(smiles_list)}", flush=True)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES reached MolGNet feature extraction: {smiles}")
            data = mol_to_graph_data_dic(mol).to(device)
            node_features = model(data)
            batch = torch.zeros(data.x.size(0), dtype=torch.long, device=device)
            pooled = global_add_pool(node_features, batch).reshape(-1)
            features[smiles] = pooled.detach().cpu().numpy().astype(np.float32)
    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return features


def ensure_features(
    train: pd.DataFrame,
    benchmark: pd.DataFrame,
    cache_path: Path,
    molgnet_model: Path,
    prott5_model: Path,
    device: torch.device,
) -> dict[str, dict[str, np.ndarray]]:
    cache = load_cache(cache_path)
    sequences = list(
        dict.fromkeys(
            train["pretkcat_model_sequence"].astype(str).tolist()
            + benchmark["pretkcat_model_sequence"].astype(str).tolist()
        )
    )
    smiles = list(
        dict.fromkeys(
            train["canonical_smiles"].astype(str).tolist()
            + benchmark["smiles"].astype(str).tolist()
        )
    )
    missing_sequences = [sequence for sequence in sequences if sequence not in cache["sequence"]]
    missing_smiles = [smile for smile in smiles if smile not in cache["smiles"]]
    print(
        f"PreTKcat feature cache: {len(sequences) - len(missing_sequences)}/{len(sequences)} sequences, "
        f"{len(smiles) - len(missing_smiles)}/{len(smiles)} SMILES already cached",
        flush=True,
    )
    if missing_sequences:
        cache["sequence"].update(compute_sequence_features(missing_sequences, prott5_model, device))
        save_cache(cache_path, cache)
        print(f"Updated sequence feature cache: {cache_path}", flush=True)
    if missing_smiles:
        cache["smiles"].update(compute_smiles_features(missing_smiles, molgnet_model, device))
        save_cache(cache_path, cache)
        print(f"Updated SMILES feature cache: {cache_path}", flush=True)
    return cache


def build_feature_matrix(df: pd.DataFrame, cache: dict[str, dict[str, np.ndarray]], *, training: bool) -> np.ndarray:
    sequence_column = "pretkcat_model_sequence"
    smiles_column = "canonical_smiles" if training else "smiles"
    temp_columns = ["Temp_K_norm", "Inv_Temp_norm"] if training else ["temp_k_norm", "inv_temp_norm"]
    smiles = np.vstack([cache["smiles"][str(value)] for value in df[smiles_column]]).astype(np.float32)
    sequence = np.vstack([cache["sequence"][str(value)] for value in df[sequence_column]]).astype(np.float32)
    temperature = df[temp_columns].to_numpy(dtype=np.float32)
    return np.concatenate([smiles, sequence, temperature], axis=1).astype(np.float32)


def model_metadata(args: argparse.Namespace, train_rows: int, feature_dim: int) -> dict[str, Any]:
    return {
        "train_data": str(args.train_data.resolve()),
        "sample_train_size": args.sample_train_size,
        "train_rows": train_rows,
        "feature_dim": feature_dim,
        "n_estimators": args.n_estimators,
        "random_state": args.random_state,
    }


def load_cached_model(path: Path, expected_meta: dict[str, Any], force_retrain: bool) -> ExtraTreesRegressor | None:
    if force_retrain or not path.exists():
        return None
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict) or "model" not in payload or "metadata" not in payload:
        return None
    cached_meta = payload["metadata"]
    keys = ["train_data", "sample_train_size", "train_rows", "feature_dim", "n_estimators", "random_state"]
    if all(cached_meta.get(key) == expected_meta.get(key) for key in keys):
        print(f"Reusing cached PreTKcat ExtraTrees model: {path}", flush=True)
        return payload["model"]
    return None


def save_model(path: Path, model: ExtraTreesRegressor, metadata: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump({"model": model, "metadata": metadata}, handle, protocol=pickle.HIGHEST_PROTOCOL)


def fit_or_load_model(args: argparse.Namespace, x_train: np.ndarray, y_train: np.ndarray) -> ExtraTreesRegressor:
    metadata = model_metadata(args, len(y_train), x_train.shape[1])
    model = load_cached_model(args.model_cache, metadata, args.force_retrain)
    if model is not None:
        return model
    print(
        f"Training PreTKcat ExtraTrees on {len(y_train)} rows, feature_dim={x_train.shape[1]}, "
        f"n_estimators={args.n_estimators}",
        flush=True,
    )
    model = ExtraTreesRegressor(
        n_estimators=args.n_estimators,
        random_state=args.random_state,
        n_jobs=args.n_jobs,
    )
    model.fit(x_train, y_train)
    save_model(args.model_cache, model, metadata)
    print(f"Saved PreTKcat ExtraTrees model: {args.model_cache}", flush=True)
    return model


def main() -> None:
    args = parse_args()
    device = resolve_device(args.device)
    print(f"PreTKcat device: {device}", flush=True)
    prepare_imports(args.pretkcat_root)

    benchmark = load_input(args.input, args.limit)
    if benchmark.empty:
        raise ValueError(f"No usable benchmark rows found in {args.input}")
    train = load_training(args.train_data, args.sample_train_size)
    if train.empty:
        raise ValueError(f"No usable training rows found in {args.train_data}")
    print(f"Benchmark rows: {len(benchmark)}; training rows: {len(train)}", flush=True)

    cache = ensure_features(train, benchmark, args.feature_cache, args.molgnet_model, args.prott5_model, device)
    x_train = build_feature_matrix(train, cache, training=True)
    y_train = train["label_log10"].to_numpy(dtype=np.float32)
    x_benchmark = build_feature_matrix(benchmark, cache, training=False)
    model = fit_or_load_model(args, x_train, y_train)

    pred_log10 = model.predict(x_benchmark).astype(np.float64)
    out = benchmark.copy()
    out["prediction_log10"] = pred_log10
    out["prediction_kcat"] = np.power(10.0, pred_log10)
    out["prediction_column"] = "prediction_log10"
    out["pretkcat_training_mode"] = "public_pretkcat_training_data_extratrees_fit"
    out["pretkcat_train_rows"] = len(train)
    out["pretkcat_feature_dim"] = x_train.shape[1]
    out["pretkcat_n_estimators"] = args.n_estimators

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote PreTKcat predictions: {args.output}", flush=True)
    print(f"Rows: {len(out)}", flush=True)


if __name__ == "__main__":
    main()
