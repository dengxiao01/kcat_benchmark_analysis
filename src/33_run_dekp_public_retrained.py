#!/usr/bin/env python3
"""Train and evaluate a DEKP public-data retrained kcat model.

This runner reuses DEKP's MetaDecoder architecture and graph feature recipe, but
keeps the workflow self-contained for the benchmark project. The official DEKP
repository does not ship final kcat weights, so this is intentionally reported as
DEKP-public-retrained.
"""

from __future__ import annotations

import argparse
import math
import os
import pickle
import random
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch_geometric
from rdkit import Chem
from rdkit import RDLogger
from scipy.spatial import cKDTree
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset
from torch_geometric.data import Data
from tqdm import tqdm


RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_TRAIN = BASE / "external_methods" / "DEKP" / "datasets" / "kcat_dataset.csv"
DEFAULT_BENCHMARK = BASE / "data" / "final" / "dekp" / "dekp_kcat_input_valid_smiles.csv"
DEFAULT_METADATA = BASE / "data" / "final" / "dekp" / "dekp_kcat_input_valid_smiles_metadata.csv"
DEFAULT_OUTPUT = BASE / "data" / "final" / "dekp" / "dekp_public_retrained_kcat_input_output.csv"
DEFAULT_MODEL = BASE / "data" / "final" / "dekp" / "dekp_public_retrained_model.pt"
DEFAULT_RUN_REPORT = BASE / "reports" / "tables" / "dekp_public_retrained_run_report.csv"
DEFAULT_FEATURE_CACHE = BASE / "data" / "final" / "dekp" / "dekp_public_retrained_feature_cache.pkl"
DEFAULT_GRAPH_CACHE = BASE / "data" / "final" / "dekp" / "dekp_public_retrained_graph_cache.pkl"
DEFAULT_TRAIN_STRUCTURES = BASE / "external_methods" / "DEKP" / "structures" / "public_kcat" / "AlphaFold"
DEFAULT_BENCHMARK_STRUCTURES = BASE / "external_methods" / "DEKP" / "structures" / "benchmark" / "AlphaFold"
DEFAULT_DEKP_ROOT = BASE / "external_methods" / "DEKP" / "DEKP"
DEFAULT_UNIKP_ROOT = BASE / "external_methods" / "CatPred" / "external" / "UniKP"

SMILES_PATTERN = re.compile(
    r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
)
AA_VOCAB = {
    "<pad>": 0,
    "<unk>": 1,
    "<bos>": 2,
    "<eos>": 3,
    "A": 4,
    "C": 5,
    "D": 6,
    "E": 7,
    "F": 8,
    "G": 9,
    "H": 10,
    "I": 11,
    "K": 12,
    "L": 13,
    "M": 14,
    "N": 15,
    "P": 16,
    "Q": 17,
    "R": 18,
    "S": 19,
    "T": 20,
    "V": 21,
    "W": 22,
    "Y": 23,
    "X": 24,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train DEKP-public-retrained and predict benchmark kcat.")
    parser.add_argument("--train-data", type=Path, default=DEFAULT_TRAIN)
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--model-out", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--run-report", type=Path, default=DEFAULT_RUN_REPORT)
    parser.add_argument("--feature-cache", type=Path, default=DEFAULT_FEATURE_CACHE)
    parser.add_argument("--graph-cache", type=Path, default=DEFAULT_GRAPH_CACHE)
    parser.add_argument("--train-structures", type=Path, default=DEFAULT_TRAIN_STRUCTURES)
    parser.add_argument("--benchmark-structures", type=Path, default=DEFAULT_BENCHMARK_STRUCTURES)
    parser.add_argument("--dekp-root", type=Path, default=DEFAULT_DEKP_ROOT)
    parser.add_argument("--unikp-root", type=Path, default=DEFAULT_UNIKP_ROOT)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--hidden", type=int, default=32)
    parser.add_argument("--n-layer", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.3)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--seed", type=int, default=3407)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--sample-train-size", type=int, default=0, help="Use first N public rows for smoke tests.")
    parser.add_argument("--limit", type=int, default=0, help="Predict only first N benchmark rows. 0 means all.")
    parser.add_argument("--force-features", action="store_true")
    parser.add_argument("--exclude-exact-benchmark-pairs", action="store_true", default=True)
    return parser.parse_args()


def prepare_imports(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(args.dekp_root))
    sys.path.insert(0, str(args.dekp_root / "Encode"))
    sys.path.insert(0, str(args.unikp_root))
    scatter = BASE / "external_methods" / "kcatnet_scatter_src"
    if scatter.exists():
        sys.path.insert(0, str(scatter))
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")


def resolve_device(requested: str) -> torch.device:
    if requested.startswith("cuda") and not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device(requested)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def canonical_smiles(smiles: object) -> str:
    text = str(smiles).strip()
    if not text or text.lower() == "nan":
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def load_public_training(args: argparse.Namespace, benchmark_pairs: set[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    train = pd.read_csv(args.train_data, sep="\t")
    required = {"ECNumber", "Organism", "Smiles", "Substrate", "Sequence", "Type", "Label", "Unit", "UniprotID"}
    missing = sorted(required.difference(train.columns))
    if missing:
        raise ValueError(f"{args.train_data} is missing required columns: {', '.join(missing)}")
    train = train.copy()
    train["Label"] = pd.to_numeric(train["Label"], errors="coerce")
    for column in ["ECNumber", "Organism", "Smiles", "Substrate", "Sequence", "Type", "Unit", "UniprotID"]:
        train[column] = train[column].fillna("").astype(str).str.strip()
    train["canonical_smiles"] = train["Smiles"].map(canonical_smiles)
    before = len(train)
    train = train[
        train["Label"].notna()
        & (train["Sequence"] != "")
        & (train["UniprotID"] != "")
        & (train["canonical_smiles"] != "")
        & (~train["canonical_smiles"].str.contains(".", regex=False))
    ].copy()
    after_clean = len(train)
    train["pair_key"] = train["Sequence"].astype(str) + "||" + train["canonical_smiles"].astype(str)
    overlap_rows = int(train["pair_key"].isin(benchmark_pairs).sum())
    if args.exclude_exact_benchmark_pairs:
        train = train[~train["pair_key"].isin(benchmark_pairs)].copy()
    if args.sample_train_size > 0:
        train = train.head(args.sample_train_size).copy()
    stats = {
        "public_rows_raw": before,
        "public_rows_after_clean": after_clean,
        "public_exact_pair_overlap_rows": overlap_rows,
        "public_rows_used": len(train),
    }
    return train.reset_index(drop=True), stats


def load_benchmark(args: argparse.Namespace) -> pd.DataFrame:
    benchmark = pd.read_csv(args.benchmark)
    required = {"ECNumber", "Organism", "Smiles", "Substrate", "Sequence", "Type", "Label", "Unit", "UniprotID", "entry_id", "dekp_row_id"}
    missing = sorted(required.difference(benchmark.columns))
    if missing:
        raise ValueError(f"{args.benchmark} is missing required columns: {', '.join(missing)}")
    benchmark = benchmark.copy()
    benchmark["Label"] = pd.to_numeric(benchmark["Label"], errors="coerce")
    for column in ["ECNumber", "Organism", "Smiles", "Substrate", "Sequence", "Type", "Unit", "UniprotID", "entry_id"]:
        benchmark[column] = benchmark[column].fillna("").astype(str).str.strip()
    benchmark["canonical_smiles"] = benchmark["Smiles"].map(canonical_smiles)
    benchmark = benchmark[(benchmark["canonical_smiles"] != "") & benchmark["Label"].notna()].copy()
    if args.limit > 0:
        benchmark = benchmark.head(args.limit).copy()
    return benchmark.reset_index(drop=True)


def assign_cids(train: pd.DataFrame, benchmark: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_smiles = pd.concat([train["canonical_smiles"], benchmark["canonical_smiles"]], ignore_index=True)
    smiles_to_cid = {smiles: index for index, smiles in enumerate(sorted(all_smiles.unique()))}
    train = train.copy()
    benchmark = benchmark.copy()
    train["CID"] = train["canonical_smiles"].map(smiles_to_cid).astype(int)
    benchmark["CID"] = benchmark["canonical_smiles"].map(smiles_to_cid).astype(int)
    train["Smiles"] = train["canonical_smiles"]
    benchmark["Smiles"] = benchmark["canonical_smiles"]
    return train, benchmark


def load_pickle(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("rb") as handle:
        return pickle.load(handle)


def save_pickle(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)


def smiles_to_trfm(smiles: list[str], unikp_root: Path) -> np.ndarray:
    import __main__  # noqa: PLC0415

    from build_vocab import Vocab, WordVocab  # noqa: PLC0415
    from pretrain_trfm import TrfmSeq2seq  # noqa: PLC0415
    from utils import split  # noqa: PLC0415

    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    __main__.Vocab = Vocab
    __main__.WordVocab = WordVocab
    vocab = WordVocab.load_vocab(str(unikp_root / "vocab.pkl"))

    def get_inputs(smiles_tokens: str) -> tuple[list[int], list[int]]:
        seq_len = 220
        tokens = smiles_tokens.split()
        if len(tokens) > 218:
            tokens = tokens[:109] + tokens[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in tokens]
        ids = [sos_index] + ids + [eos_index]
        seg = [1] * len(ids)
        padding = [pad_index] * (seq_len - len(ids))
        ids.extend(padding)
        seg.extend(padding)
        return ids, seg

    x_id = []
    for smi in smiles:
        ids, _ = get_inputs(split(smi))
        x_id.append(ids)
    xid = torch.tensor(x_id)
    model = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
    state = torch.load(unikp_root / "trfm_12_23000.pkl", map_location="cpu")
    model.load_state_dict(state)
    model.eval()
    with torch.no_grad():
        return model.encode(torch.t(xid)).astype(np.float32)


def ensure_trfm_features(df: pd.DataFrame, args: argparse.Namespace) -> dict[int, np.ndarray]:
    cache = {} if args.force_features else load_pickle(args.feature_cache)
    trfm_cache: dict[str, np.ndarray] = cache.get("trfm_by_smiles", {}) if isinstance(cache, dict) else {}
    cid_by_smiles = (
        df[["canonical_smiles", "CID"]]
        .drop_duplicates()
        .sort_values("CID")
        .set_index("canonical_smiles")["CID"]
        .to_dict()
    )
    missing = [smiles for smiles in cid_by_smiles if smiles not in trfm_cache]
    if missing:
        print(f"Generating UniKP trfm features for {len(missing)} SMILES", flush=True)
        features = smiles_to_trfm(missing, args.unikp_root)
        for smiles, feature in zip(missing, features, strict=True):
            trfm_cache[smiles] = feature
        save_pickle(args.feature_cache, {"trfm_by_smiles": trfm_cache})
    return {int(cid): np.asarray(trfm_cache[smiles], dtype=np.float32) for smiles, cid in cid_by_smiles.items()}


def structure_paths(args: argparse.Namespace) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for directory in [args.benchmark_structures, args.train_structures]:
        if not directory.exists():
            continue
        for path in directory.glob("*.pdb"):
            paths.setdefault(path.stem, path)
    return paths


def graph_from_pdb(pdb_path: Path, radius: float = 10.0, nneighbor: int = 20) -> Data:
    from extract_pdb_feature import atom_idx, get_cb, get_geo_feat  # noqa: PLC0415

    coords = torch.tensor(parse_pdb_robust(pdb_path, get_cb), dtype=torch.float32)
    if coords.ndim != 3 or coords.shape[0] < 2:
        raise ValueError(f"Not enough residues in {pdb_path}")
    query = coords[:, atom_idx["CA"], :].numpy()
    tree = cKDTree(query)
    src: list[int] = []
    dst: list[int] = []
    for index, point in enumerate(query):
        neighbors = tree.query_ball_point(point, r=radius)
        neighbors = [neighbor for neighbor in neighbors if neighbor != index]
        neighbors = sorted(neighbors, key=lambda neighbor: float(np.linalg.norm(query[neighbor] - point)))[:nneighbor]
        for neighbor in neighbors:
            src.append(index)
            dst.append(neighbor)
    if not src:
        for index in range(coords.shape[0] - 1):
            src.extend([index, index + 1])
            dst.extend([index + 1, index])
    edge_index = torch.tensor([src, dst], dtype=torch.long)
    node, edge = get_geo_feat(coords, edge_index, D_count=1)
    node = torch.nan_to_num(node.float())
    edge = torch.nan_to_num(edge.float())
    return Data(x=node, edge_index=edge_index, edge_attr=edge, name=pdb_path.stem)


def parse_pdb_robust(pdb_path: Path, get_cb_func) -> np.ndarray:
    fillna = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    residues: list[np.ndarray] = []
    current_pos: str | None = None
    current_atoms: dict[str, np.ndarray] = {}

    def flush_current() -> None:
        nonlocal current_atoms
        if not current_atoms:
            return
        side_chain = [coord for atom, coord in current_atoms.items() if atom not in {"N", "CA", "C", "O"}]
        if side_chain:
            r_group = np.stack(side_chain, axis=0).mean(axis=0)
        else:
            r_group = current_atoms.get("CA", fillna)
        residue = np.stack(
            [
                current_atoms.get("N", fillna),
                current_atoms.get("CA", fillna),
                current_atoms.get("C", fillna),
                current_atoms.get("O", fillna),
                r_group,
            ],
            axis=0,
        ).astype(np.float32)
        residues.append(residue)
        current_atoms = {}

    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[0:6].strip()
            if record == "TER":
                flush_current()
                current_pos = None
                continue
            if record != "ATOM":
                continue
            pos = line[21:26].strip()
            if current_pos is not None and pos != current_pos:
                flush_current()
            current_pos = pos
            atom = line[12:16].strip()
            if atom.startswith("H"):
                continue
            try:
                current_atoms[atom] = np.array(
                    [line[30:38].strip(), line[38:46].strip(), line[46:54].strip()],
                    dtype=np.float32,
                )
            except ValueError:
                current_atoms[atom] = fillna
    flush_current()
    if not residues:
        raise ValueError(f"No ATOM residues parsed from {pdb_path}")
    coords = np.stack(residues, axis=0)
    cb = get_cb_func(coords[:, 0], coords[:, 1], coords[:, 2])[:, None]
    return np.concatenate([coords, cb.astype(np.float32)], axis=1)


def ensure_graphs(df: pd.DataFrame, args: argparse.Namespace) -> tuple[dict[str, Data], list[str]]:
    cache = {} if args.force_features else load_pickle(args.graph_cache)
    graph_cache: dict[str, Data] = cache.get("graphs", {}) if isinstance(cache, dict) else {}
    failed: dict[str, str] = cache.get("failed", {}) if isinstance(cache, dict) else {}
    paths = structure_paths(args)
    needed = sorted(df["UniprotID"].dropna().astype(str).unique())
    missing = [uniprot for uniprot in needed if uniprot not in graph_cache]
    for uniprot in tqdm(missing, desc="Building DEKP PDB graphs"):
        pdb_path = paths.get(uniprot)
        if pdb_path is None:
            failed[uniprot] = "missing_pdb"
            continue
        try:
            graph_cache[uniprot] = graph_from_pdb(pdb_path)
            failed.pop(uniprot, None)
        except Exception as exc:  # noqa: BLE001
            failed[uniprot] = f"{type(exc).__name__}: {exc}"
        if (len(graph_cache) + len(failed)) % 100 == 0:
            save_pickle(args.graph_cache, {"graphs": graph_cache, "failed": failed})
    save_pickle(args.graph_cache, {"graphs": graph_cache, "failed": failed})
    return graph_cache, sorted(failed)


def build_smiles_vocab(smiles_values: list[str]) -> dict[str, int]:
    tokens = {"<pad>": 0, "<unk>": 1, "<bos>": 2, "<eos>": 3, "<mask>": 4}
    for smiles in smiles_values:
        for token in SMILES_PATTERN.findall(smiles):
            if token not in tokens:
                tokens[token] = len(tokens)
    return tokens


def encode_sequence(sequence: str, max_len: int) -> torch.Tensor:
    sequence = re.sub(r"[UZOB*]", "X", sequence.strip().rstrip("*"))
    ids = [AA_VOCAB["<bos>"]]
    ids.extend(AA_VOCAB.get(token, AA_VOCAB["<unk>"]) for token in sequence[: max_len - 2])
    ids.append(AA_VOCAB["<eos>"])
    if len(ids) < max_len:
        ids.extend([AA_VOCAB["<pad>"]] * (max_len - len(ids)))
    return torch.tensor(ids[:max_len], dtype=torch.long)


def encode_smiles(smiles: str, vocab: dict[str, int], max_len: int) -> torch.Tensor:
    ids = [vocab["<bos>"]]
    ids.extend(vocab.get(token, vocab["<unk>"]) for token in SMILES_PATTERN.findall(smiles)[: max_len - 2])
    ids.append(vocab["<eos>"])
    if len(ids) < max_len:
        ids.extend([vocab["<pad>"]] * (max_len - len(ids)))
    return torch.tensor(ids[:max_len], dtype=torch.long)


class DEKPDataset(Dataset):
    def __init__(
        self,
        df: pd.DataFrame,
        trfm: dict[int, np.ndarray],
        graphs: dict[str, Data],
        smiles_vocab: dict[str, int],
        max_protein_len: int,
        max_smiles_len: int,
    ) -> None:
        self.df = df.reset_index(drop=True).copy()
        self.trfm = trfm
        self.graphs = graphs
        self.smiles_vocab = smiles_vocab
        self.max_protein_len = max_protein_len
        self.max_smiles_len = max_smiles_len

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, index: int):
        row = self.df.iloc[index]
        uniprot = str(row["UniprotID"])
        cid = int(row["CID"])
        sequence = str(row["Sequence"])
        smiles = str(row["Smiles"])
        return (
            self.graphs[uniprot],
            encode_sequence(sequence, self.max_protein_len),
            encode_smiles(smiles, self.smiles_vocab, self.max_smiles_len),
            torch.tensor(self.trfm[cid], dtype=torch.float32),
            torch.tensor(float(row["Label"]), dtype=torch.float32),
            str(row.get("entry_id", "")),
            int(row.get("dekp_row_id", -1)),
            uniprot,
        )


def graph_collate_fn(batch):
    graph_batch = torch_geometric.data.Batch.from_data_list([item[0] for item in batch])
    protein = torch.stack([item[1] for item in batch], dim=0)
    smiles = torch.stack([item[2] for item in batch], dim=0)
    feature = torch.stack([item[3] for item in batch], dim=0)
    label = torch.stack([item[4] for item in batch], dim=0)
    entry_ids = [item[5] for item in batch]
    row_ids = [item[6] for item in batch]
    uniprots = [item[7] for item in batch]
    return graph_batch, protein, smiles, feature, label, entry_ids, row_ids, uniprots


def evaluate(model: nn.Module, loader: DataLoader, device: torch.device, criterion: nn.Module) -> dict[str, float]:
    model.eval()
    losses: list[float] = []
    preds: list[np.ndarray] = []
    labels: list[np.ndarray] = []
    with torch.no_grad():
        for graph, protein, smiles, feature, label, *_ in loader:
            graph = graph.to(device)
            protein = protein.to(device)
            smiles = smiles.to(device)
            feature = feature.to(device)
            label = label.to(device)
            pred = model(graph, protein, smiles, feature)
            loss = criterion(pred, label)
            losses.append(float(loss.item()))
            preds.append(pred.detach().cpu().numpy())
            labels.append(label.detach().cpu().numpy())
    y_pred = np.concatenate(preds)
    y_true = np.concatenate(labels)
    err = y_pred - y_true
    return {
        "loss": float(np.mean(losses)),
        "rmse": float(np.sqrt(np.mean(np.square(err)))),
        "mae": float(np.mean(np.abs(err))),
        "pearson": float(pd.Series(y_true).corr(pd.Series(y_pred), method="pearson")) if len(y_true) > 1 else math.nan,
    }


def train_model(
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    trfm: dict[int, np.ndarray],
    graphs: dict[str, Data],
    smiles_vocab: dict[str, int],
    args: argparse.Namespace,
    device: torch.device,
) -> tuple[nn.Module, dict[str, float]]:
    from fine_tune import MetaDecoder  # noqa: PLC0415

    max_protein_len = int(max(train_df["Sequence"].str.len().max(), valid_df["Sequence"].str.len().max()) + 2)
    max_smiles_len = int(max(train_df["Smiles"].str.len().max(), valid_df["Smiles"].str.len().max()) + 2)
    max_protein_len = min(max(max_protein_len, 16), 2500)
    max_smiles_len = min(max(max_smiles_len, 16), 500)
    train_ds = DEKPDataset(train_df, trfm, graphs, smiles_vocab, max_protein_len, max_smiles_len)
    valid_ds = DEKPDataset(valid_df, trfm, graphs, smiles_vocab, max_protein_len, max_smiles_len)
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=0, collate_fn=graph_collate_fn)
    valid_loader = DataLoader(valid_ds, batch_size=args.batch_size, shuffle=False, num_workers=0, collate_fn=graph_collate_fn)

    model = MetaDecoder(
        seq_vocab_size=len(AA_VOCAB),
        smi_vocab_size=len(smiles_vocab),
        feature_dim_list=[1024],
        hidden=args.hidden,
        num_layers=args.n_layer,
        protein_len=max_protein_len,
        smi_len=max_smiles_len,
        dropout=args.dropout,
        kernel_size=9,
    ).to(device)
    criterion = nn.MSELoss(reduction="mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=3e-4)

    best_state = None
    best = {"valid_rmse": math.inf, "epoch": 0}
    stale = 0
    for epoch in range(1, args.epochs + 1):
        model.train()
        epoch_losses: list[float] = []
        for graph, protein, smiles, feature, label, *_ in train_loader:
            graph = graph.to(device)
            protein = protein.to(device)
            smiles = smiles.to(device)
            feature = feature.to(device)
            label = label.to(device)
            optimizer.zero_grad(set_to_none=True)
            pred = model(graph, protein, smiles, feature)
            loss = criterion(pred, label)
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.item()))
        valid_metrics = evaluate(model, valid_loader, device, criterion)
        print(
            f"Epoch {epoch}/{args.epochs} train_loss={np.mean(epoch_losses):.4f} "
            f"valid_rmse={valid_metrics['rmse']:.4f} valid_mae={valid_metrics['mae']:.4f} "
            f"valid_pearson={valid_metrics['pearson']:.4f}",
            flush=True,
        )
        if valid_metrics["rmse"] < best["valid_rmse"]:
            best = {
                "valid_rmse": valid_metrics["rmse"],
                "valid_mae": valid_metrics["mae"],
                "valid_pearson": valid_metrics["pearson"],
                "train_loss": float(np.mean(epoch_losses)),
                "epoch": epoch,
                "max_protein_len": max_protein_len,
                "max_smiles_len": max_smiles_len,
                "smiles_vocab_size": len(smiles_vocab),
            }
            best_state = {key: value.detach().cpu() for key, value in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= args.patience:
                print(f"Early stopping after {epoch} epochs", flush=True)
                break
    if best_state is None:
        raise RuntimeError("Training did not produce a model state.")
    model.load_state_dict(best_state)
    return model, best


def predict_benchmark(
    model: nn.Module,
    benchmark: pd.DataFrame,
    trfm: dict[int, np.ndarray],
    graphs: dict[str, Data],
    smiles_vocab: dict[str, int],
    best: dict[str, float],
    args: argparse.Namespace,
    device: torch.device,
) -> pd.DataFrame:
    max_protein_len = int(best["max_protein_len"])
    max_smiles_len = int(best["max_smiles_len"])
    dataset = DEKPDataset(benchmark, trfm, graphs, smiles_vocab, max_protein_len, max_smiles_len)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=0, collate_fn=graph_collate_fn)
    model.eval()
    rows: list[dict[str, object]] = []
    with torch.no_grad():
        for graph, protein, smiles, feature, label, entry_ids, row_ids, uniprots in loader:
            graph = graph.to(device)
            protein = protein.to(device)
            smiles = smiles.to(device)
            feature = feature.to(device)
            pred = model(graph, protein, smiles, feature).detach().cpu().numpy()
            for value, true_value, entry_id, row_id, uniprot in zip(pred, label.numpy(), entry_ids, row_ids, uniprots, strict=True):
                rows.append(
                    {
                        "dekp_row_id": row_id,
                        "entry_id": entry_id,
                        "uniprot_id": uniprot,
                        "prediction_log10": float(value),
                        "prediction_kcat": float(np.power(10.0, value)),
                        "true_kcat_log10": float(true_value),
                        "prediction_column": "prediction_log10",
                        "dekp_training_mode": "DEKP-public-retrained",
                        "dekp_feature_mode": "trfm+sequence_cnn+structure_graph",
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    prepare_imports(args)
    from fine_tune import MetaDecoder  # noqa: F401, PLC0415

    set_seed(args.seed)
    device = resolve_device(args.device)
    print(f"DEKP-public-retrained device: {device}", flush=True)
    benchmark = load_benchmark(args)
    benchmark_pairs = set(benchmark["Sequence"].astype(str) + "||" + benchmark["canonical_smiles"].astype(str))
    train, stats = load_public_training(args, benchmark_pairs)
    train, benchmark = assign_cids(train, benchmark)
    all_rows = pd.concat([train, benchmark], ignore_index=True, sort=False)

    trfm = ensure_trfm_features(all_rows, args)
    graphs, failed_graphs = ensure_graphs(all_rows, args)
    if failed_graphs:
        failed_set = set(failed_graphs)
        before_train = len(train)
        before_benchmark = len(benchmark)
        train = train[~train["UniprotID"].isin(failed_set)].copy()
        benchmark = benchmark[~benchmark["UniprotID"].isin(failed_set)].copy()
        print(
            f"Dropped rows with failed graphs: train {before_train - len(train)}, "
            f"benchmark {before_benchmark - len(benchmark)}",
            flush=True,
        )
    if train.empty or benchmark.empty:
        raise ValueError("No usable train or benchmark rows after graph filtering.")

    train_df, valid_df = train_test_split(train, test_size=0.1, random_state=args.seed, shuffle=True)
    train_df = train_df.reset_index(drop=True)
    valid_df = valid_df.reset_index(drop=True)
    smiles_vocab = build_smiles_vocab(all_rows["Smiles"].dropna().astype(str).tolist())

    print(
        f"Training rows: {len(train_df)}; validation rows: {len(valid_df)}; "
        f"benchmark rows: {len(benchmark)}; trfm CIDs: {len(trfm)}; graphs: {len(graphs)}",
        flush=True,
    )
    model, best = train_model(train_df, valid_df, trfm, graphs, smiles_vocab, args, device)
    pred = predict_benchmark(model, benchmark, trfm, graphs, smiles_vocab, best, args, device)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pred.to_csv(args.output, index=False)
    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "best": best,
            "args": vars(args),
            "aa_vocab": AA_VOCAB,
            "smiles_vocab": smiles_vocab,
            "stats": stats,
        },
        args.model_out,
    )
    report = pd.DataFrame(
        [
            {
                **stats,
                "train_rows_after_split": len(train_df),
                "valid_rows": len(valid_df),
                "benchmark_rows_predicted": len(pred),
                "failed_graph_uniprots": len(failed_graphs),
                "best_epoch": best["epoch"],
                "best_valid_rmse": best["valid_rmse"],
                "best_valid_mae": best["valid_mae"],
                "best_valid_pearson": best["valid_pearson"],
                "feature_mode": "trfm+sequence_cnn+structure_graph",
                "excluded_exact_benchmark_pairs": bool(args.exclude_exact_benchmark_pairs),
            }
        ]
    )
    args.run_report.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.run_report, index=False)
    print(f"Wrote DEKP predictions: {args.output}", flush=True)
    print(f"Wrote DEKP model: {args.model_out}", flush=True)
    print(f"Wrote run report: {args.run_report}", flush=True)


if __name__ == "__main__":
    main()
