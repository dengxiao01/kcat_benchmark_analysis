#!/usr/bin/env python3
"""Run official DLKcat predictions on the finalized kcat benchmark."""

from __future__ import annotations

import argparse
import math
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE / "data" / "final" / "benchmark_ready_catpred.csv"
DEFAULT_OUT_DIR = BASE / "data" / "final" / "dlkcat"
DEFAULT_DLKCAT = BASE / "external_methods" / "DLKcat_official" / "DeeplearningApproach"

REQUIRED_COLUMNS = {
    "entry_id",
    "species",
    "reaction_id",
    "gene_id",
    "substrate_name",
    "SMILES",
    "sequence",
    "true_kcat",
    "true_kcat_log10",
    "source_database",
    "match_level",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict kcat with the official DLKcat trained model.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--dlkcat-root", type=Path, default=DEFAULT_DLKCAT)
    parser.add_argument("--device", default="cpu")
    return parser.parse_args()


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def load_dlkcat_model(dlkcat_root: Path, device: torch.device):
    code_dir = dlkcat_root / "Code" / "example"
    input_dir = dlkcat_root / "Data" / "input"
    model_path = (
        dlkcat_root
        / "Results"
        / "output"
        / "all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50"
    )
    if not input_dir.exists():
        raise FileNotFoundError(f"DLKcat input directory not found: {input_dir}. Unzip Data/input.zip first.")
    if not model_path.exists():
        raise FileNotFoundError(f"DLKcat trained model not found: {model_path}")

    sys.path.insert(0, str(code_dir))
    import model as dlkcat_model  # noqa: PLC0415

    dictionaries = {
        "fingerprint": load_pickle(input_dir / "fingerprint_dict.pickle"),
        "atom": load_pickle(input_dir / "atom_dict.pickle"),
        "bond": load_pickle(input_dir / "bond_dict.pickle"),
        "edge": load_pickle(input_dir / "edge_dict.pickle"),
        "word": load_pickle(input_dir / "sequence_dict.pickle"),
    }
    network = dlkcat_model.KcatPrediction(
        device,
        len(dictionaries["fingerprint"]),
        len(dictionaries["word"]),
        20,
        3,
        11,
        3,
        3,
    ).to(device)
    state = torch.load(model_path, map_location=device)
    network.load_state_dict(state)
    network.eval()
    return network, dictionaries


def split_sequence(sequence: str, word_dict: dict, ngram: int = 3) -> np.ndarray:
    sequence = "-" + str(sequence).strip() + "="
    words = []
    for index in range(len(sequence) - ngram + 1):
        token = sequence[index : index + ngram]
        words.append(word_dict.get(token, 0))
    return np.asarray(words, dtype=np.int64)


def create_atoms(mol, atom_dict: dict) -> np.ndarray:
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    for atom in mol.GetAromaticAtoms():
        atoms[atom.GetIdx()] = (atoms[atom.GetIdx()], "aromatic")
    return np.asarray([atom_dict.get(atom, 0) for atom in atoms], dtype=np.int64)


def create_ijbonddict(mol, bond_dict: dict):
    i_jbond_dict = defaultdict(list)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_id = bond_dict.get(str(bond.GetBondType()), 0)
        i_jbond_dict[i].append((j, bond_id))
        i_jbond_dict[j].append((i, bond_id))
    return i_jbond_dict


def extract_fingerprints(atoms: np.ndarray, i_jbond_dict, fingerprint_dict: dict, edge_dict: dict, radius: int = 2):
    if len(atoms) == 1 or radius == 0:
        return np.asarray([fingerprint_dict.get(int(atom), 0) for atom in atoms], dtype=np.int64)

    nodes = atoms
    i_jedge_dict = i_jbond_dict
    fingerprints = []
    for _ in range(radius):
        fingerprints = []
        for i, j_edge in i_jedge_dict.items():
            neighbors = [(nodes[j], edge) for j, edge in j_edge]
            fingerprint = (nodes[i], tuple(sorted(neighbors)))
            fingerprints.append(fingerprint_dict.get(fingerprint, 0))
        nodes = fingerprints

        next_edges = defaultdict(list)
        for i, j_edge in i_jedge_dict.items():
            for j, edge in j_edge:
                both_side = tuple(sorted((nodes[i], nodes[j])))
                next_edges[i].append((j, edge_dict.get((both_side, edge), 0)))
        i_jedge_dict = next_edges
    return np.asarray(fingerprints, dtype=np.int64)


def predict_one(model, dictionaries: dict, sequence: str, smiles: str, device: torch.device) -> tuple[float | None, str]:
    text = str(smiles).strip()
    if not text or text == "nan":
        return None, "empty_smiles"
    if "." in text:
        return None, "multicomponent_or_invalid_smiles"
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None, "invalid_smiles"

    try:
        mol_h = Chem.AddHs(mol)
        atoms = create_atoms(mol_h, dictionaries["atom"])
        bonds = create_ijbonddict(mol_h, dictionaries["bond"])
        fingerprints = extract_fingerprints(atoms, bonds, dictionaries["fingerprint"], dictionaries["edge"])
        adjacency = Chem.GetAdjacencyMatrix(mol_h)
        words = split_sequence(sequence, dictionaries["word"])

        inputs = [
            torch.LongTensor(fingerprints).to(device),
            torch.FloatTensor(adjacency).to(device),
            torch.LongTensor(words).to(device),
        ]
        with torch.no_grad():
            pred_log2 = float(model.forward(inputs).detach().cpu().reshape(-1)[0])
        pred_kcat = math.pow(2.0, pred_log2)
        return pred_kcat, "success"
    except Exception as exc:  # Keep one bad molecule from aborting the whole benchmark.
        return None, f"prediction_error:{type(exc).__name__}"


def validate(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    validate(df, args.input)
    df = df.reset_index(drop=True).copy()
    df["dlkcat_row_id"] = range(len(df))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    metadata_cols = [
        "dlkcat_row_id",
        "entry_id",
        "species",
        "reaction_id",
        "gene_id",
        "uniprot_id",
        "ec_number",
        "substrate_name",
        "SMILES",
        "sequence",
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
    df[[col for col in metadata_cols if col in df.columns]].to_csv(
        args.out_dir / "dlkcat_kcat_input_metadata.csv", index=False
    )
    df[[col for col in ["dlkcat_row_id", "entry_id", "species", "true_kcat", "true_kcat_log10", "source_database", "match_level"] if col in df.columns]].to_csv(
        args.out_dir / "dlkcat_kcat_input_truth.csv", index=False
    )
    df[["substrate_name", "SMILES", "sequence"]].rename(
        columns={"substrate_name": "Substrate Name", "SMILES": "Substrate SMILES", "sequence": "Protein Sequence"}
    ).to_csv(args.out_dir / "dlkcat_official_input.tsv", sep="\t", index=False)

    requested = torch.device(args.device)
    device = requested if requested.type != "cuda" or torch.cuda.is_available() else torch.device("cpu")
    model, dictionaries = load_dlkcat_model(args.dlkcat_root, device)

    rows = []
    for index, row in df.iterrows():
        if index == 0 or (index + 1) % 100 == 0 or index + 1 == len(df):
            print(f"DLKcat prediction {index + 1}/{len(df)}", flush=True)
        pred_kcat, status = predict_one(model, dictionaries, row["sequence"], row["SMILES"], device)
        pred_log10 = math.log10(pred_kcat) if pred_kcat and pred_kcat > 0 else np.nan
        rows.append(
            {
                "dlkcat_row_id": int(row["dlkcat_row_id"]),
                "entry_id": row["entry_id"],
                "prediction_kcat": pred_kcat,
                "prediction_log10": pred_log10,
                "prediction_status": status,
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(args.out_dir / "dlkcat_kcat_input_output.csv", index=False)
    print(f"Wrote DLKcat predictions: {args.out_dir / 'dlkcat_kcat_input_output.csv'}")
    print(f"Rows: {len(out)}; success: {(out['prediction_status'] == 'success').sum()}")


if __name__ == "__main__":
    main()
