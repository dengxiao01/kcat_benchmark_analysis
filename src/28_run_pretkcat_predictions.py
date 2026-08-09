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
import hashlib
import json
import os
import pickle
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit import RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.ensemble import ExtraTreesRegressor
from transformers import T5EncoderModel, T5Tokenizer


RDLogger.DisableLog("rdApp.*")
UNCHARGER = rdMolStandardize.Uncharger()
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
NEAR_SEQUENCE_IDENTITY = 80.0
NEAR_CHEMICAL_TANIMOTO = 0.80
NEAR_ALIGNMENT_COVERAGE = 50.0

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_input_valid_smiles.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "pretkcat" / "pretkcat_kcat_input_output.csv"
DEFAULT_TRAIN = BASE / "external_methods" / "PreTKcat" / "datasets" / "DLTKcat_data" / "kcat_merge_DLTKcat.csv"
DEFAULT_PRETKCAT_ROOT = BASE / "external_methods" / "PreTKcat"
DEFAULT_MOLGNET = DEFAULT_PRETKCAT_ROOT / "MolGNet.pt"
DEFAULT_PROTT5 = BASE / "external_methods" / "CataPro" / "models" / "prot_t5_xl_uniref50"
DEFAULT_FEATURE_CACHE = BASE / "data" / "final" / "pretkcat" / "pretkcat_feature_cache.pkl"
DEFAULT_MODEL_CACHE = BASE / "data" / "final" / "pretkcat" / "pretkcat_extratrees_model.pkl"
PRETKCAT_WRAPPER_VERSION = "1.3.0"


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
    parser.add_argument(
        "--overlap-policy",
        choices=["raw-public", "exact-excluded", "near-excluded"],
        default="exact-excluded",
        help=(
            "Training-corpus policy: keep the raw public rows, remove exact standardized "
            "sequence-parent pairs, or remove joint sequence/chemical near neighbors."
        ),
    )
    parser.add_argument(
        "--keep-exact-benchmark-pairs",
        action="store_true",
        help="Deprecated alias for --overlap-policy raw-public.",
    )
    parser.add_argument("--diamond", default="diamond")
    parser.add_argument("--near-sequence-identity", type=float, default=NEAR_SEQUENCE_IDENTITY)
    parser.add_argument("--near-chemical-tanimoto", type=float, default=NEAR_CHEMICAL_TANIMOTO)
    parser.add_argument("--near-alignment-coverage", type=float, default=NEAR_ALIGNMENT_COVERAGE)
    parser.add_argument(
        "--overlap-audit-output",
        type=Path,
        default=None,
        help="Optional JSON path for the training-overlap audit.",
    )
    parser.add_argument(
        "--audit-overlap-only",
        action="store_true",
        help="Stop after overlap filtering and audit generation; do not build features or fit the model.",
    )
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


def chemical_parent_key(smiles: object) -> str:
    text = str(smiles).strip()
    if not text or text == "nan":
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    try:
        parent = rdMolStandardize.FragmentParent(mol)
        parent = UNCHARGER.uncharge(parent)
        Chem.SanitizeMol(parent)
    except Exception:
        parent = Chem.Mol(mol)
        Chem.SanitizeMol(parent)
    parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    try:
        key = Chem.MolToInchiKey(parent).split("-", 1)[0]
    except Exception:
        key = ""
    return key or parent_smiles


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
    df["canonical_smiles"] = df["smiles"].map(canonical_smiles)
    df["chemical_parent_key"] = df["smiles"].map(chemical_parent_key)
    df = df[(df["canonical_smiles"] != "") & (df["chemical_parent_key"] != "")].copy()
    df["smiles"] = df["canonical_smiles"]
    if limit > 0:
        df = df.head(limit).copy()
    return df.reset_index(drop=True)


def load_training(path: Path) -> pd.DataFrame:
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
    train["chemical_parent_key"] = train["smiles"].map(chemical_parent_key)
    train = train[(train["canonical_smiles"] != "") & (train["chemical_parent_key"] != "")].copy()
    train["pretkcat_model_sequence"] = train["seq"].map(pretkcat_model_sequence)
    train["label_log10"] = np.log10(train["kcat"].astype(float))
    return train.reset_index(drop=True)


def pair_keys(df: pd.DataFrame) -> pd.Series:
    return df["pretkcat_model_sequence"].astype(str) + "||" + df["chemical_parent_key"].astype(str)


def parent_fingerprint(smiles: object):
    text = str(smiles).strip()
    molecule = Chem.MolFromSmiles(text)
    if molecule is None:
        return None
    try:
        parent = rdMolStandardize.FragmentParent(molecule)
        parent = UNCHARGER.uncharge(parent)
        Chem.SanitizeMol(parent)
    except Exception:
        parent = Chem.Mol(molecule)
        Chem.SanitizeMol(parent)
    return MORGAN.GetFingerprint(parent)


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="ascii") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n{sequence}\n")


def near_neighbor_training_rows(
    train: pd.DataFrame,
    benchmark: pd.DataFrame,
    *,
    diamond: str,
    sequence_identity: float,
    chemical_tanimoto: float,
    alignment_coverage: float,
) -> tuple[set[int], set[int], dict[str, Any]]:
    sequence_to_train_rows: dict[str, list[int]] = defaultdict(list)
    for index, sequence in enumerate(train["pretkcat_model_sequence"].astype(str)):
        sequence_to_train_rows[sequence].append(index)

    query_sequences = {
        f"q{index}": sequence
        for index, sequence in enumerate(sorted(set(benchmark["pretkcat_model_sequence"].astype(str))))
    }
    subject_sequences = {
        f"s{index}": sequence for index, sequence in enumerate(sorted(sequence_to_train_rows))
    }
    sequence_to_query = {sequence: identifier for identifier, sequence in query_sequences.items()}
    subject_to_sequence = {identifier: sequence for identifier, sequence in subject_sequences.items()}

    with tempfile.TemporaryDirectory(prefix="pretkcat_near_exclusion_") as tmp:
        tmp_path = Path(tmp)
        query_fasta = tmp_path / "benchmark.fasta"
        subject_fasta = tmp_path / "training.fasta"
        database = tmp_path / "training"
        output = tmp_path / "hits.tsv"
        write_fasta(query_fasta, query_sequences)
        write_fasta(subject_fasta, subject_sequences)
        subprocess.run(
            [diamond, "makedb", "--in", str(subject_fasta), "-d", str(database), "--quiet"],
            check=True,
        )
        subprocess.run(
            [
                diamond,
                "blastp",
                "-d",
                str(database),
                "-q",
                str(query_fasta),
                "-o",
                str(output),
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "qcovhsp",
                "scovhsp",
                "--id",
                str(sequence_identity),
                "--query-cover",
                str(alignment_coverage),
                "--subject-cover",
                str(alignment_coverage),
                "--max-target-seqs",
                str(max(1, len(subject_sequences))),
                "--sensitive",
                "--threads",
                "8",
                "--quiet",
            ],
            check=True,
        )
        hits: dict[str, list[str]] = defaultdict(list)
        hit_count = 0
        if output.exists() and output.stat().st_size:
            frame = pd.read_csv(
                output,
                sep="\t",
                names=["query", "subject", "identity", "query_coverage", "subject_coverage"],
            )
            hit_count = len(frame)
            for row in frame.itertuples(index=False):
                hits[str(row.query)].append(str(row.subject))

    train_fingerprints = [parent_fingerprint(value) for value in train["canonical_smiles"]]
    removed: set[int] = set()
    benchmark_rows_with_joint_neighbor: set[int] = set()
    for benchmark_index, row in benchmark.iterrows():
        query_fp = parent_fingerprint(row["smiles"])
        if query_fp is None:
            continue
        query_id = sequence_to_query[str(row["pretkcat_model_sequence"])]
        candidate_indices: list[int] = []
        for subject_id in hits.get(query_id, []):
            candidate_indices.extend(sequence_to_train_rows[subject_to_sequence[subject_id]])
        if not candidate_indices:
            continue
        similarities = DataStructs.BulkTanimotoSimilarity(
            query_fp, [train_fingerprints[index] for index in candidate_indices]
        )
        matched = [
            index for index, similarity in zip(candidate_indices, similarities)
            if float(similarity) >= chemical_tanimoto
        ]
        if matched:
            removed.update(matched)
            benchmark_rows_with_joint_neighbor.add(int(benchmark_index))

    return removed, benchmark_rows_with_joint_neighbor, {
        "near_sequence_identity_threshold_percent": sequence_identity,
        "near_chemical_tanimoto_threshold": chemical_tanimoto,
        "near_alignment_coverage_threshold_percent": alignment_coverage,
        "diamond_sequence_hits": hit_count,
        "benchmark_rows_with_diamond_joint_neighbor": len(benchmark_rows_with_joint_neighbor),
        "benchmark_rows_with_joint_neighbor": len(benchmark_rows_with_joint_neighbor),
        "benchmark_rows_exact_added_to_near_exclusion": 0,
        "train_rows_joint_neighbor": len(removed),
        "near_neighbor_definition": (
            "same training row with DIAMOND identity and bidirectional alignment coverage at or "
            "above thresholds plus Morgan-radius-2 parent-structure Tanimoto at or above threshold"
        ),
    }


def filter_training_pairs(
    train: pd.DataFrame,
    benchmark: pd.DataFrame,
    *,
    overlap_policy: str,
    sample_size: int,
    diamond: str,
    near_sequence_identity: float,
    near_chemical_tanimoto: float,
    near_alignment_coverage: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    benchmark_pair_series = pair_keys(benchmark)
    train_pair_series = pair_keys(train)
    benchmark_keys = set(benchmark_pair_series)
    train_keys = set(train_pair_series)
    exact_overlap_mask = train_pair_series.isin(benchmark_keys)
    benchmark_exact_overlap_mask = benchmark_pair_series.isin(train_keys)
    raw_usable_rows = len(train)
    source_overlap_rows = int(exact_overlap_mask.sum())
    source_overlap_unique_pairs = int(train_pair_series.loc[exact_overlap_mask].nunique())
    benchmark_overlap_rows = int(benchmark_exact_overlap_mask.sum())
    near_audit: dict[str, Any] = {
        "near_sequence_identity_threshold_percent": near_sequence_identity,
        "near_chemical_tanimoto_threshold": near_chemical_tanimoto,
        "near_alignment_coverage_threshold_percent": near_alignment_coverage,
        "diamond_sequence_hits": 0,
        "benchmark_rows_with_diamond_joint_neighbor": 0,
        "benchmark_rows_with_joint_neighbor": 0,
        "benchmark_rows_exact_added_to_near_exclusion": 0,
        "train_rows_joint_neighbor": 0,
        "train_rows_exact_added_to_near_exclusion": 0,
        "near_neighbor_definition": "not_computed_for_this_policy",
    }

    if overlap_policy == "raw-public":
        removal_mask = pd.Series(False, index=train.index)
        policy = "raw_public_corpus_no_benchmark_exclusion"
    elif overlap_policy == "exact-excluded":
        removal_mask = exact_overlap_mask
        policy = "standardized_sequence_parent_identity_pair_disjoint"
    elif overlap_policy == "near-excluded":
        near_rows, near_benchmark_rows, near_audit = near_neighbor_training_rows(
            train,
            benchmark,
            diamond=diamond,
            sequence_identity=near_sequence_identity,
            chemical_tanimoto=near_chemical_tanimoto,
            alignment_coverage=near_alignment_coverage,
        )
        exact_indices = set(train.index[exact_overlap_mask])
        exact_benchmark_indices = set(benchmark.index[benchmark_exact_overlap_mask])
        exact_rows_not_recovered = exact_indices.difference(near_rows)
        exact_benchmark_rows_not_recovered = exact_benchmark_indices.difference(near_benchmark_rows)
        near_rows.update(exact_indices)
        near_benchmark_rows.update(exact_benchmark_indices)
        near_audit["train_rows_exact_added_to_near_exclusion"] = len(exact_rows_not_recovered)
        near_audit["benchmark_rows_exact_added_to_near_exclusion"] = len(
            exact_benchmark_rows_not_recovered
        )
        near_audit["train_rows_joint_neighbor"] = len(near_rows)
        near_audit["benchmark_rows_with_joint_neighbor"] = len(near_benchmark_rows)
        removal_mask = pd.Series(train.index.isin(near_rows), index=train.index)
        policy = "joint_sequence_chemical_near_neighbor_disjoint"
    else:
        raise ValueError(f"Unknown overlap policy: {overlap_policy}")

    filtered = train.loc[~removal_mask].copy()
    removed_rows = int(removal_mask.sum())
    removed_exact_rows = int((removal_mask & exact_overlap_mask).sum())
    removed_near_only_rows = removed_rows - removed_exact_rows
    rows_after_exclusion = len(filtered)
    if sample_size > 0:
        filtered = filtered.head(sample_size).copy()
    audit = {
        "training_overlap_policy": policy,
        "training_overlap_variant": overlap_policy,
        "train_rows_raw_usable": raw_usable_rows,
        "train_rows_source_exact_pair_overlap": source_overlap_rows,
        "train_unique_source_exact_pair_overlap": source_overlap_unique_pairs,
        "benchmark_rows_source_exact_pair_overlap": benchmark_overlap_rows,
        "train_rows_removed_exact_pair": removed_exact_rows,
        "train_rows_removed_near_only": removed_near_only_rows,
        "train_rows_removed_total": removed_rows,
        "train_rows_after_pair_exclusion": rows_after_exclusion,
        "train_rows_fitted": len(filtered),
        "benchmark_rows": len(benchmark),
        "benchmark_exact_pair_keys": len(benchmark_keys),
        "benchmark_pair_key_sha256": hashlib.sha256(
            "\n".join(sorted(benchmark_keys)).encode("utf-8")
        ).hexdigest(),
        "pair_identity_definition": "model_sequence_plus_uncharged_largest_fragment_connectivity_identity",
        **near_audit,
    }
    return filtered.reset_index(drop=True), audit


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


def model_metadata(
    args: argparse.Namespace,
    train_rows: int,
    feature_dim: int,
    audit: dict[str, Any],
) -> dict[str, Any]:
    return {
        "wrapper_version": PRETKCAT_WRAPPER_VERSION,
        "train_data": str(args.train_data.resolve()),
        "sample_train_size": args.sample_train_size,
        "train_rows": train_rows,
        "feature_dim": feature_dim,
        "n_estimators": args.n_estimators,
        "random_state": args.random_state,
        **audit,
    }


def load_cached_model(path: Path, expected_meta: dict[str, Any], force_retrain: bool) -> ExtraTreesRegressor | None:
    if force_retrain or not path.exists():
        return None
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict) or "model" not in payload or "metadata" not in payload:
        return None
    cached_meta = payload["metadata"]
    keys = [
        "wrapper_version",
        "train_data",
        "sample_train_size",
        "train_rows",
        "feature_dim",
        "n_estimators",
        "random_state",
        "training_overlap_policy",
        "train_rows_removed_exact_pair",
        "train_rows_removed_near_only",
        "train_rows_removed_total",
        "benchmark_pair_key_sha256",
        "near_sequence_identity_threshold_percent",
        "near_chemical_tanimoto_threshold",
        "near_alignment_coverage_threshold_percent",
    ]
    if all(cached_meta.get(key) == expected_meta.get(key) for key in keys):
        print(f"Reusing cached PreTKcat ExtraTrees model: {path}", flush=True)
        return payload["model"]
    return None


def save_model(path: Path, model: ExtraTreesRegressor, metadata: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump({"model": model, "metadata": metadata}, handle, protocol=pickle.HIGHEST_PROTOCOL)


def fit_or_load_model(
    args: argparse.Namespace,
    x_train: np.ndarray,
    y_train: np.ndarray,
    audit: dict[str, Any],
) -> ExtraTreesRegressor:
    metadata = model_metadata(args, len(y_train), x_train.shape[1], audit)
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
    train = load_training(args.train_data)
    if train.empty:
        raise ValueError(f"No usable training rows found in {args.train_data}")
    overlap_policy = "raw-public" if args.keep_exact_benchmark_pairs else args.overlap_policy
    train, training_audit = filter_training_pairs(
        train,
        benchmark,
        overlap_policy=overlap_policy,
        sample_size=args.sample_train_size,
        diamond=args.diamond,
        near_sequence_identity=args.near_sequence_identity,
        near_chemical_tanimoto=args.near_chemical_tanimoto,
        near_alignment_coverage=args.near_alignment_coverage,
    )
    if train.empty:
        raise ValueError("No training rows remain after benchmark-overlap exclusion.")
    print(
        f"Benchmark rows: {len(benchmark)}; raw usable training rows: "
        f"{training_audit['train_rows_raw_usable']}; exact-pair rows found: "
        f"{training_audit['train_rows_source_exact_pair_overlap']}; removed total: "
        f"{training_audit['train_rows_removed_total']}; fitted rows: {len(train)}",
        flush=True,
    )
    if args.overlap_audit_output is not None:
        args.overlap_audit_output.parent.mkdir(parents=True, exist_ok=True)
        with args.overlap_audit_output.open("w", encoding="utf-8") as handle:
            json.dump(training_audit, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(f"Wrote PreTKcat overlap audit: {args.overlap_audit_output}", flush=True)
    if args.audit_overlap_only:
        return

    cache = ensure_features(train, benchmark, args.feature_cache, args.molgnet_model, args.prott5_model, device)
    x_train = build_feature_matrix(train, cache, training=True)
    y_train = train["label_log10"].to_numpy(dtype=np.float32)
    x_benchmark = build_feature_matrix(benchmark, cache, training=False)
    model = fit_or_load_model(args, x_train, y_train, training_audit)

    pred_log10 = model.predict(x_benchmark).astype(np.float64)
    out = benchmark.copy()
    out["prediction_log10"] = pred_log10
    out["prediction_kcat"] = np.power(10.0, pred_log10)
    out["prediction_column"] = "prediction_log10"
    out["pretkcat_training_mode"] = (
        "PreTKcat-public-reconstructed_" + training_audit["training_overlap_variant"]
    )
    out["pretkcat_wrapper_version"] = PRETKCAT_WRAPPER_VERSION
    out["pretkcat_training_overlap_policy"] = training_audit["training_overlap_policy"]
    out["pretkcat_pair_identity_definition"] = training_audit["pair_identity_definition"]
    out["pretkcat_train_rows"] = len(train)
    out["pretkcat_train_rows_raw_usable"] = training_audit["train_rows_raw_usable"]
    out["pretkcat_train_rows_source_exact_pair_overlap"] = training_audit[
        "train_rows_source_exact_pair_overlap"
    ]
    out["pretkcat_train_rows_removed_exact_pair"] = training_audit["train_rows_removed_exact_pair"]
    out["pretkcat_train_rows_removed_near_only"] = training_audit["train_rows_removed_near_only"]
    out["pretkcat_train_rows_removed_total"] = training_audit["train_rows_removed_total"]
    out["pretkcat_train_rows_after_pair_exclusion"] = training_audit["train_rows_after_pair_exclusion"]
    out["pretkcat_near_sequence_identity_threshold_percent"] = training_audit[
        "near_sequence_identity_threshold_percent"
    ]
    out["pretkcat_near_chemical_tanimoto_threshold"] = training_audit[
        "near_chemical_tanimoto_threshold"
    ]
    out["pretkcat_near_alignment_coverage_threshold_percent"] = training_audit[
        "near_alignment_coverage_threshold_percent"
    ]
    out["pretkcat_benchmark_rows_with_joint_neighbor"] = training_audit[
        "benchmark_rows_with_joint_neighbor"
    ]
    out["pretkcat_feature_dim"] = x_train.shape[1]
    out["pretkcat_n_estimators"] = args.n_estimators

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote PreTKcat predictions: {args.output}", flush=True)
    print(f"Rows: {len(out)}", flush=True)


if __name__ == "__main__":
    main()
