#!/usr/bin/env python
"""Predict with the official UniKP kcat ExtraTrees pickle.

Run this script with an old scikit-learn environment, e.g. condaPY36lin.
"""

from __future__ import print_function

import argparse
import math
import pickle
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_FEATURES = BASE / "data" / "final" / "unikp" / "unikp_official_features.npy"
DEFAULT_ROWS = BASE / "data" / "final" / "unikp" / "unikp_official_feature_rows.csv"
DEFAULT_MODEL = BASE / "external_methods" / "UniKP_official" / "models" / "UniKP for kcat.pkl"
DEFAULT_OUTPUT = BASE / "data" / "final" / "unikp" / "unikp_kcat_input_output.csv"


def parse_args():
    parser = argparse.ArgumentParser(description="Predict log10(kcat) with official UniKP model.")
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main():
    args = parse_args()
    features = np.load(str(args.features))
    rows = pd.read_csv(str(args.rows))
    with args.model.open("rb") as handle:
        model = pickle.load(handle)
    prediction_log10 = model.predict(features)
    out = rows.copy()
    out["prediction_log10"] = prediction_log10
    out["prediction_kcat"] = [math.pow(10.0, value) for value in prediction_log10]
    out["prediction_status"] = "success"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(str(args.output), index=False)
    print("Wrote UniKP predictions:", args.output)
    print("Rows:", len(out))


if __name__ == "__main__":
    main()
