#!/usr/bin/env python3
"""Run official MTLKP kcat inference and normalize its output table."""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
MTLKP_ROOT = BASE / "external_methods" / "ecm_benchmark_end" / "etgems_web" / "script" / "mtlkp"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MTLKP kcat predictions.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--task", default="Kcat", choices=["Kcat", "Km"])
    return parser.parse_args()


def local_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(MTLKP_ROOT))
    os.chdir(MTLKP_ROOT)

    from omegaconf import OmegaConf
    import prediction_api

    prediction_api.setup_logger = local_logger
    config = OmegaConf.load(MTLKP_ROOT / "config.yaml")
    predictor = prediction_api.Predictor(args.task, config)
    predictor.main(str(input_path))

    official_output = Path(str(input_path).replace(".csv", f"_{args.task}_output.csv"))
    if not official_output.exists():
        raise FileNotFoundError(f"MTLKP did not create expected output: {official_output}")
    if official_output.resolve() != output_path:
        shutil.copy2(official_output, output_path)

    out = pd.read_csv(output_path)
    if f"{args.task}(log10)" in out.columns:
        out["prediction_log10"] = pd.to_numeric(out[f"{args.task}(log10)"], errors="coerce")
    if args.task in out.columns:
        out["prediction_kcat"] = pd.to_numeric(out[args.task], errors="coerce")
    out.to_csv(output_path, index=False)
    print(f"Wrote normalized MTLKP predictions: {output_path}")


if __name__ == "__main__":
    main()
