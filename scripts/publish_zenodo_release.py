#!/usr/bin/env python3
"""Create, upload, verify, and optionally publish the benchmark Zenodo record."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests


BASE = Path(__file__).resolve().parent.parent
RELEASE_DIR = BASE / "release" / "zenodo"
MANIFEST_PATH = RELEASE_DIR / "assets_manifest.csv"
PUBLIC_MANIFEST = BASE / "zenodo_assets_manifest.csv"
TOKEN_PATH = BASE / "zenodo.txt"
STATE_PATH = RELEASE_DIR / "zenodo_state.json"
API_ROOT = "https://zenodo.org/api/deposit/depositions"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true", help="Publish after every asset passes remote verification.")
    parser.add_argument("--skip-upload", action="store_true", help="Only create/read the draft and update public links.")
    parser.add_argument("--force-upload", action="store_true", help="Upload files even when name and size already match.")
    parser.add_argument("--workers", type=int, default=1, help="Number of files to upload concurrently (default: 1).")
    return parser.parse_args()


def token() -> str:
    value = TOKEN_PATH.read_text(encoding="utf-8").strip()
    if not value:
        raise RuntimeError(f"Zenodo token file is empty: {TOKEN_PATH}")
    return value


def load_rows() -> list[dict[str, str]]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save_state(data: dict[str, Any]) -> None:
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def auth_headers(access_token: str) -> dict[str, str]:
    return {"Authorization": f"Bearer {access_token}"}


def check_response(response: requests.Response) -> requests.Response:
    if not response.ok:
        body = response.text[:2000]
        raise RuntimeError(f"Zenodo API {response.status_code} for {response.url}: {body}")
    return response


def deposition_metadata() -> dict[str, Any]:
    description = (
        "Unified benchmark of published kcat prediction methods on 978 experimentally supported "
        "enzyme-substrate records from Escherichia coli and Saccharomyces cerevisiae. The deposit "
        "contains canonical benchmark tables, evaluated predictions, reports, figures, and large "
        "model assets needed to reproduce method-specific inference. Third-party model files retain "
        "their upstream licenses and citation requirements; see the included asset manifest and "
        "the GitHub METHOD_SOURCES.md file."
    )
    return {
        "metadata": {
            "title": "kcat Benchmark Analysis: Unified Benchmark, Results, and Model Assets",
            "upload_type": "dataset",
            "description": description,
            "creators": [{"name": "dengxiao01"}],
            "access_right": "open",
            "license": "other-open",
            "keywords": [
                "kcat",
                "enzyme kinetics",
                "machine learning",
                "metabolic models",
                "Escherichia coli",
                "Saccharomyces cerevisiae",
            ],
            "notes": (
                "BRENDA-derived data are used under CC BY 4.0. Other database, software, and "
                "model components retain their respective upstream terms."
            ),
        }
    }


def get_or_create_deposition(access_token: str) -> dict[str, Any]:
    if STATE_PATH.exists():
        state = json.loads(STATE_PATH.read_text(encoding="utf-8"))
        deposition_id = state["id"]
        response = requests.get(f"{API_ROOT}/{deposition_id}", headers=auth_headers(access_token), timeout=60)
        return check_response(response).json()

    response = requests.post(
        API_ROOT,
        headers={**auth_headers(access_token), "Content-Type": "application/json"},
        json=deposition_metadata(),
        timeout=60,
    )
    deposition = check_response(response).json()
    save_state({"id": deposition["id"]})
    print(f"Created Zenodo draft deposition {deposition['id']}")
    return deposition


def local_path(row: dict[str, str]) -> Path:
    path = BASE / row["local_path"]
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def upload_asset(bucket_url: str, row: dict[str, str], access_token: str) -> None:
    path = local_path(row)
    destination = f"{bucket_url}/{urllib.parse.quote(row['asset_name'])}"
    with tempfile.NamedTemporaryFile("w", prefix="zenodo-curl-", suffix=".conf") as config:
        config.write(f'header = "Authorization: Bearer {access_token}"\n')
        config.write("fail-with-body\n")
        config.write("silent\n")
        config.write("show-error\n")
        config.write("connect-timeout = 60\n")
        config.write("max-time = 43200\n")
        config.flush()

        for attempt in range(1, 4):
            print(f"Uploading {path.name}, attempt {attempt}/3")
            with tempfile.NamedTemporaryFile("w+", prefix="zenodo-response-", suffix=".txt") as response_body:
                result = subprocess.run(
                    [
                        "curl",
                        "--config",
                        config.name,
                        "--upload-file",
                        str(path),
                        "--header",
                        "Content-Type: application/octet-stream",
                        "--output",
                        response_body.name,
                        "--write-out",
                        "%{http_code}\t%{size_upload}\t%{speed_upload}\t%{time_total}",
                        destination,
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                response_body.seek(0)
                body = response_body.read(2000)
            if result.returncode == 0:
                print(f"Uploaded: {row['asset_name']} (http/bytes/speed/time: {result.stdout})")
                return
            message = result.stderr.strip() or body.strip() or result.stdout.strip()
            print(f"Upload failed for {row['asset_name']} on attempt {attempt}: {message}", flush=True)
            if attempt == 3:
                raise RuntimeError(f"curl failed after 3 attempts for {row['asset_name']}")
            time.sleep(10 * attempt)


def remote_files(deposition: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {item["filename"]: item for item in deposition.get("files", [])}


def verify_remote(rows: list[dict[str, str]], deposition: dict[str, Any]) -> None:
    files = remote_files(deposition)
    errors = []
    for row in rows:
        item = files.get(row["asset_name"])
        if item is None:
            errors.append(f"missing {row['asset_name']}")
        elif int(item["filesize"]) != int(row["size_bytes"]):
            errors.append(f"size mismatch {row['asset_name']}: {item['filesize']} != {row['size_bytes']}")
    if errors:
        raise RuntimeError("Remote verification failed: " + "; ".join(errors))


def update_public_references(rows: list[dict[str, str]], deposition: dict[str, Any]) -> None:
    record_id = str(deposition["id"])
    doi = deposition.get("metadata", {}).get("prereserve_doi", {}).get("doi", "")
    for row in rows:
        row["zenodo_record_id"] = record_id
        quoted = urllib.parse.quote(row["asset_name"])
        row["download_url"] = f"https://zenodo.org/records/{record_id}/files/{quoted}?download=1"

    with PUBLIC_MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    readme = BASE / "README.md"
    content = readme.read_text(encoding="utf-8")
    start = "<!-- ZENODO_START -->"
    end = "<!-- ZENODO_END -->"
    if doi:
        block = (
            f"{start}\n"
            f"[![DOI](https://zenodo.org/badge/DOI/{doi}.svg)](https://doi.org/{doi})\n\n"
            f"**Large assets:** [Zenodo record {record_id}](https://zenodo.org/records/{record_id}) "
            f"(DOI: `{doi}`). File checksums and restore paths are in `zenodo_assets_manifest.csv`.\n"
            f"{end}"
        )
    else:
        block = (
            f"{start}\n**Large assets:** [Zenodo draft/record {record_id}]"
            f"(https://zenodo.org/records/{record_id}).\n{end}"
        )
    before, marker, remainder = content.partition(start)
    if not marker or end not in remainder:
        raise RuntimeError("README Zenodo markers are missing")
    _, _, after = remainder.partition(end)
    readme.write_text(before + block + after, encoding="utf-8")


def main() -> None:
    args = parse_args()
    access_token = token()
    rows = load_rows()
    deposition = get_or_create_deposition(access_token)
    bucket_url = deposition["links"]["bucket"]
    existing = remote_files(deposition)

    if not args.skip_upload:
        pending = []
        for row in rows:
            item = existing.get(row["asset_name"])
            same_size = item and int(item["filesize"]) == int(row["size_bytes"])
            if same_size and not args.force_upload:
                print(f"Remote file already matches by name and size: {row['asset_name']}")
                continue
            pending.append(row)
        if pending:
            pending.sort(key=lambda row: int(row["size_bytes"]))
            workers = max(1, min(args.workers, len(pending)))
            print(f"Uploading {len(pending)} files with {workers} concurrent workers")
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(upload_asset, bucket_url, row, access_token): row for row in pending
                }
                for future in as_completed(futures):
                    row = futures[future]
                    future.result()
                    print(f"Completed worker task: {row['asset_name']}")

    response = requests.get(f"{API_ROOT}/{deposition['id']}", headers=auth_headers(access_token), timeout=60)
    deposition = check_response(response).json()
    verify_remote(rows, deposition)
    update_public_references(rows, deposition)

    if args.publish:
        if deposition.get("submitted"):
            print(f"Zenodo deposition {deposition['id']} is already published")
        else:
            response = requests.post(deposition["links"]["publish"], headers=auth_headers(access_token), timeout=120)
            deposition = check_response(response).json()
            update_public_references(rows, deposition)
            print(f"Published Zenodo record {deposition['id']}")
    else:
        print(f"Zenodo draft {deposition['id']} verified; rerun with --publish to make it public")


if __name__ == "__main__":
    main()
