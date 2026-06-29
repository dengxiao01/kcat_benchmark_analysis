#!/usr/bin/env python3
"""Download missing AlphaFold structures for DEKP benchmark rows."""

from __future__ import annotations

import argparse
import json
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent.parent
DEFAULT_MISSING = BASE / "data" / "final" / "dekp" / "dekp_missing_structure_rows.csv"
DEFAULT_OUT_DIR = BASE / "external_methods" / "DEKP" / "structures" / "benchmark" / "AlphaFold"
DEFAULT_REPORT = BASE / "reports" / "tables" / "dekp_alphafold_download_report.csv"
ALPHAFOLD_URLS = [
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v5.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v3.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v2.pdb",
]
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download missing AlphaFold PDB files for DEKP evaluation.")
    parser.add_argument("--missing", type=Path, default=DEFAULT_MISSING)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--sleep", type=float, default=0.05)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--no-check-certificate", action="store_true")
    parser.add_argument("--limit", type=int, default=0, help="Download only the first N unique UniProt IDs. 0 means all.")
    return parser.parse_args()


def read_url(url: str, timeout: int, ssl_context: ssl.SSLContext | None = None) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "kcat-benchmark-dekp/1.0"})
    with urllib.request.urlopen(request, timeout=timeout, context=ssl_context) as response:
        return response.read()


def alphafold_pdb_urls(
    uniprot: str,
    timeout: int,
    ssl_context: ssl.SSLContext | None = None,
) -> tuple[list[str], str]:
    urls: list[str] = []
    api_error = ""
    try:
        payload = read_url(ALPHAFOLD_API.format(uniprot=uniprot), timeout, ssl_context)
        records = json.loads(payload.decode("utf-8"))
        for record in records:
            pdb_url = record.get("pdbUrl")
            if pdb_url:
                urls.append(str(pdb_url))
    except Exception as exc:  # noqa: BLE001
        api_error = f"{type(exc).__name__}: {exc}"
    for template in ALPHAFOLD_URLS:
        urls.append(template.format(uniprot=uniprot))
    return list(dict.fromkeys(urls)), api_error


def download_one(
    uniprot: str,
    out_path: Path,
    timeout: int,
    ssl_context: ssl.SSLContext | None = None,
) -> dict[str, object]:
    if out_path.exists() and out_path.stat().st_size > 0:
        return {"uniprot_id": uniprot, "status": "exists", "url": "", "bytes": out_path.stat().st_size, "error": ""}
    last_error = ""
    urls, api_error = alphafold_pdb_urls(uniprot, timeout, ssl_context)
    if api_error:
        last_error = f"API {api_error}"
    for url in urls:
        try:
            payload = read_url(url, timeout, ssl_context)
            if b"\nATOM" not in payload and not payload.startswith(b"ATOM"):
                last_error = "downloaded file has no ATOM records"
                continue
            out_path.write_bytes(payload)
            return {"uniprot_id": uniprot, "status": "downloaded", "url": url, "bytes": len(payload), "error": ""}
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code}"
        except Exception as exc:  # noqa: BLE001
            last_error = f"{type(exc).__name__}: {exc}"
    return {"uniprot_id": uniprot, "status": "failed", "url": "", "bytes": 0, "error": last_error}


def main() -> None:
    args = parse_args()
    missing = pd.read_csv(args.missing)
    if "uniprot_id" not in missing.columns:
        raise ValueError(f"{args.missing} is missing uniprot_id column")
    uniprots = sorted(set(missing["uniprot_id"].dropna().astype(str)))
    if args.limit > 0:
        uniprots = uniprots[: args.limit]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ssl_context = ssl._create_unverified_context() if args.no_check_certificate else None
    rows = []
    for index, uniprot in enumerate(uniprots, start=1):
        out_path = args.out_dir / f"{uniprot}.pdb"
        row = download_one(uniprot, out_path, args.timeout, ssl_context)
        rows.append(row)
        print(f"{index}/{len(uniprots)} {uniprot}: {row['status']} {row.get('bytes', 0)}", flush=True)
        if args.sleep > 0:
            time.sleep(args.sleep)
    report = pd.DataFrame(rows)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.report, index=False)
    print(f"Wrote AlphaFold download report: {args.report}")
    print(report["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
