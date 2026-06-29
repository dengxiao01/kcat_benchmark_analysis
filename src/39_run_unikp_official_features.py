#!/usr/bin/env python3
"""Build official UniKP feature matrix for the finalized kcat benchmark."""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch


BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "unikp"
DEFAULT_UNIKP_CODE = BASE / "external_methods" / "CatPred" / "external" / "UniKP"
DEFAULT_PRETKCAT_CACHE = BASE / "data" / "final" / "pretkcat" / "pretkcat_feature_cache.pkl"
DEFAULT_KCATNET_CACHE = BASE / "data" / "final" / "kcatnet" / "kcatnet_protein_cache.pkl"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create official UniKP input features.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--unikp-code", type=Path, default=DEFAULT_UNIKP_CODE)
    parser.add_argument("--pretkcat-cache", type=Path, default=DEFAULT_PRETKCAT_CACHE)
    parser.add_argument("--kcatnet-cache", type=Path, default=DEFAULT_KCATNET_CACHE)
    return parser.parse_args()


def model_sequence(sequence: str) -> str:
    sequence = str(sequence).strip()
    return sequence[:500] + sequence[-500:] if len(sequence) > 1000 else sequence


def load_pretkcat_vectors(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with path.open("rb") as handle:
        cache = pickle.load(handle)
    return {str(k): np.asarray(v, dtype=np.float32) for k, v in cache.get("sequence", {}).items()}


def load_kcatnet_vectors(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with path.open("rb") as handle:
        cache = pickle.load(handle)
    vectors = {}
    for sequence, value in cache.items():
        if not isinstance(value, dict) or "token_representation" not in value:
            continue
        tensor = value["token_representation"]
        array = tensor.detach().cpu().numpy() if torch.is_tensor(tensor) else np.asarray(tensor)
        vectors[str(sequence)] = array.mean(axis=0).astype(np.float32)
    return vectors


def load_sequence_vectors(df: pd.DataFrame, pretkcat_cache: Path, kcatnet_cache: Path) -> tuple[np.ndarray, list[str]]:
    pretkcat = load_pretkcat_vectors(pretkcat_cache)
    kcatnet = load_kcatnet_vectors(kcatnet_cache)
    features = []
    status = []
    for _, row in df.iterrows():
        sequence = str(row["sequence"]).strip()
        seq_model = model_sequence(sequence)
        if seq_model in pretkcat:
            features.append(pretkcat[seq_model])
            status.append("pretkcat_prott5_mean_cache")
        elif sequence in pretkcat:
            features.append(pretkcat[sequence])
            status.append("pretkcat_prott5_mean_cache_original")
        elif sequence in kcatnet:
            features.append(kcatnet[sequence])
            status.append("kcatnet_prott5_mean_cache")
        else:
            features.append(np.full(1024, np.nan, dtype=np.float32))
            status.append("missing_sequence_feature")
    return np.vstack(features), status


def smiles_to_vec(smiles: list[str], code_dir: Path) -> np.ndarray:
    sys.path.insert(0, str(code_dir))
    import __main__  # noqa: PLC0415
    from build_vocab import WordVocab  # noqa: PLC0415
    from pretrain_trfm import TrfmSeq2seq  # noqa: PLC0415
    from utils import split  # noqa: PLC0415

    __main__.WordVocab = WordVocab
    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    vocab = WordVocab.load_vocab(str(code_dir / "vocab.pkl"))

    def get_inputs(sm: str):
        seq_len = 220
        tokens = split(str(sm).strip())
        if len(tokens) > 218:
            tokens = tokens[:109] + tokens[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in tokens]
        ids = [sos_index] + ids + [eos_index]
        seg = [1] * len(ids)
        padding = [pad_index] * (seq_len - len(ids))
        ids.extend(padding)
        seg.extend(padding)
        return ids, seg

    x_id, x_seg = [], []
    for sm in smiles:
        ids, seg = get_inputs(sm)
        x_id.append(ids)
        x_seg.append(seg)

    trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
    state = torch.load(code_dir / "trfm_12_23000.pkl", map_location="cpu")
    trfm.load_state_dict(state)
    trfm.eval()
    with torch.no_grad():
        return trfm.encode(torch.tensor(x_id, dtype=torch.long).T)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input).reset_index(drop=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df["unikp_row_id"] = range(len(df))
    df["unikp_model_sequence"] = df["sequence"].map(model_sequence)
    df["unikp_smiles_usable"] = ~df["SMILES"].fillna("").astype(str).str.contains(".", regex=False)

    metadata_cols = [
        "unikp_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "sequence",
        "unikp_model_sequence",
        "true_kcat",
        "true_kcat_log10",
        "unit",
        "pH",
        "temperature_c",
        "source_database",
        "match_level",
        "reference",
        "n_measurements",
        "enzyme_complex_type",
    ]

    seq_vec, seq_status = load_sequence_vectors(df, args.pretkcat_cache, args.kcatnet_cache)
    df["unikp_sequence_feature_status"] = seq_status
    valid = df[df["unikp_smiles_usable"] & np.isfinite(seq_vec).all(axis=1)].copy()
    valid_indexes = valid.index.to_numpy()

    print(f"UniKP total rows: {len(df)}")
    print(f"UniKP feature-ready rows: {len(valid)}")
    print("Computing SMILES Transformer features...", flush=True)
    smiles_vec = smiles_to_vec(valid["SMILES"].astype(str).tolist(), args.unikp_code)
    fused = np.concatenate([smiles_vec, seq_vec[valid_indexes]], axis=1).astype(np.float32)

    np.save(args.out_dir / "unikp_official_features.npy", fused)
    valid[["unikp_row_id", "entry_id", "SMILES", "sequence"]].to_csv(
        args.out_dir / "unikp_official_feature_rows.csv", index=False
    )
    df[[col for col in metadata_cols + ["unikp_smiles_usable", "unikp_sequence_feature_status"] if col in df.columns]].to_csv(
        args.out_dir / "unikp_kcat_input_metadata.csv", index=False
    )
    valid[[col for col in metadata_cols + ["unikp_smiles_usable", "unikp_sequence_feature_status"] if col in valid.columns]].to_csv(
        args.out_dir / "unikp_kcat_input_predictable_metadata.csv", index=False
    )
    df[[col for col in ["unikp_row_id", "entry_id", "species", "true_kcat", "true_kcat_log10", "source_database", "match_level"] if col in df.columns]].to_csv(
        args.out_dir / "unikp_kcat_input_truth.csv", index=False
    )
    print(f"Wrote UniKP feature matrix: {args.out_dir / 'unikp_official_features.npy'} {fused.shape}")


if __name__ == "__main__":
    main()
