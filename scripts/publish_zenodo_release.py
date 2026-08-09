#!/usr/bin/env python3
"""Create, upload, verify, and optionally publish the benchmark Zenodo record."""

from __future__ import annotations

import argparse
import csv
import hashlib
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
MD5_CACHE_PATH = RELEASE_DIR / "md5_cache.json"
API_ROOT = "https://zenodo.org/api/deposit/depositions"
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true", help="Publish after every asset passes remote verification.")
    parser.add_argument(
        "--new-version",
        action="store_true",
        help="Create a new version when the saved deposition is already published.",
    )
    parser.add_argument(
        "--initialize",
        action="store_true",
        help="Create/read the draft, update metadata and README, then exit before file checks.",
    )
    parser.add_argument(
        "--prune",
        action="store_true",
        help="Delete draft files that are not present in the prepared manifest.",
    )
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


def request_with_retry(
    method: str,
    url: str,
    *,
    attempts: int = 4,
    **kwargs: Any,
) -> requests.Response:
    """Retry transient Zenodo gateway and connection failures."""
    for attempt in range(1, attempts + 1):
        try:
            response = requests.request(method, url, **kwargs)
        except (requests.ConnectionError, requests.Timeout) as exc:
            if attempt == attempts:
                raise
            wait_seconds = min(60, 5 * 2 ** (attempt - 1))
            print(
                f"Zenodo {method.upper()} connection failure on attempt {attempt}/{attempts}: "
                f"{exc}; retrying in {wait_seconds}s",
                flush=True,
            )
            time.sleep(wait_seconds)
            continue

        if response.status_code not in RETRYABLE_STATUS_CODES or attempt == attempts:
            return response
        wait_seconds = min(60, 5 * 2 ** (attempt - 1))
        print(
            f"Zenodo {method.upper()} returned HTTP {response.status_code} on attempt "
            f"{attempt}/{attempts}; retrying in {wait_seconds}s",
            flush=True,
        )
        time.sleep(wait_seconds)

    raise AssertionError("unreachable")


def deposition_metadata() -> dict[str, Any]:
    description = (
        "Version 1.2.0 of a unified benchmark of 12 published kcat prediction methods on 1,246 "
        "model-linked enzyme-reaction-metabolite records from Escherichia coli and Saccharomyces "
        "cerevisiae. The complete resource contains 778 records with substrate-specific BRENDA "
        "support and 468 SABIO-RK-only participant-ambiguous records. The deposit contains canonical "
        "benchmark tables, method-level predictions, provenance and dependence audits, reports, "
        "figures, reproducibility code, and large model assets used for method-specific inference. "
        "Third-party files retain their upstream licenses and citation requirements."
    )
    return {
        "metadata": {
            "title": "kcat Benchmark Analysis: Unified Benchmark, Results, and Model Assets",
            "upload_type": "dataset",
            "description": description,
            "creators": [{"name": "dengxiao01"}],
            "access_right": "open",
            "license": "other-open",
            "version": "1.2.0",
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
                "model components retain their respective upstream terms. Project repository: "
                "https://github.com/dxg-9527/kcat_benchmark_analysis"
            ),
            "related_identifiers": [
                {
                    "identifier": "https://github.com/dxg-9527/kcat_benchmark_analysis",
                    "relation": "isSupplementTo",
                    "scheme": "url",
                }
            ],
        }
    }


def get_or_create_deposition(access_token: str, new_version: bool) -> dict[str, Any]:
    if STATE_PATH.exists():
        state = json.loads(STATE_PATH.read_text(encoding="utf-8"))
        deposition_id = state["id"]
        response = request_with_retry(
            "get", f"{API_ROOT}/{deposition_id}", headers=auth_headers(access_token), timeout=60
        )
        deposition = check_response(response).json()
        if not deposition.get("submitted"):
            return deposition
        if not new_version:
            raise RuntimeError(
                f"Zenodo deposition {deposition_id} is already published; rerun with --new-version."
            )

        new_version_url = deposition.get("links", {}).get("newversion")
        if not new_version_url:
            new_version_url = f"{API_ROOT}/{deposition_id}/actions/newversion"
        response = request_with_retry(
            "post", new_version_url, headers=auth_headers(access_token), timeout=120
        )
        payload = check_response(response).json()
        latest_draft = payload.get("links", {}).get("latest_draft")
        if latest_draft:
            response = request_with_retry(
                "get", latest_draft, headers=auth_headers(access_token), timeout=60
            )
            deposition = check_response(response).json()
        else:
            deposition = payload
        save_state({"id": deposition["id"], "source_version_id": deposition_id})
        print(f"Created Zenodo version draft {deposition['id']} from published record {deposition_id}")
        return deposition

    response = request_with_retry(
        "post",
        API_ROOT,
        headers={**auth_headers(access_token), "Content-Type": "application/json"},
        json=deposition_metadata(),
        timeout=60,
    )
    deposition = check_response(response).json()
    save_state({"id": deposition["id"]})
    print(f"Created Zenodo draft deposition {deposition['id']}")
    return deposition


def update_deposition_metadata(deposition: dict[str, Any], access_token: str) -> dict[str, Any]:
    def without_none(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: without_none(item) for key, item in value.items() if item is not None}
        if isinstance(value, list):
            return [without_none(item) for item in value]
        return value

    target = deposition_metadata()
    current_metadata = deposition.get("metadata", {})
    if all(
        without_none(current_metadata.get(key)) == without_none(value)
        for key, value in target["metadata"].items()
    ):
        print(f"Zenodo draft {deposition['id']} metadata already matches the release target")
        return deposition

    response = request_with_retry(
        "put",
        f"{API_ROOT}/{deposition['id']}",
        headers={**auth_headers(access_token), "Content-Type": "application/json"},
        json=target,
        timeout=120,
    )
    return check_response(response).json()


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
        config.write("fail\n")
        config.write("silent\n")
        config.write("show-error\n")
        config.write("connect-timeout = 60\n")
        config.write("max-time = 1800\n")
        config.write("speed-limit = 1024\n")
        config.write("speed-time = 120\n")
        config.flush()

        max_attempts = 6
        for attempt in range(1, max_attempts + 1):
            print(f"Uploading {path.name}, attempt {attempt}/{max_attempts}")
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
            metrics = result.stdout.strip().split("\t")
            http_code = int(metrics[0]) if metrics and metrics[0].isdigit() else 0
            uploaded_bytes = int(float(metrics[1])) if len(metrics) > 1 and metrics[1] else -1
            complete = uploaded_bytes == path.stat().st_size
            if result.returncode == 0 and 200 <= http_code < 300 and complete:
                print(f"Uploaded: {row['asset_name']} (http/bytes/speed/time: {result.stdout})")
                return
            message = result.stderr.strip() or body.strip() or result.stdout.strip()
            print(
                f"Upload failed for {row['asset_name']} on attempt {attempt} "
                f"(http={http_code}, bytes={uploaded_bytes}/{path.stat().st_size}): {message}",
                flush=True,
            )
            if attempt == max_attempts:
                raise RuntimeError(
                    f"curl failed after {max_attempts} attempts for {row['asset_name']}"
                )
            time.sleep(min(120, 15 * attempt))


def remote_files(deposition: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {item["filename"]: item for item in deposition.get("files", [])}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def local_md5(row: dict[str, str], cache: dict[str, str]) -> str:
    name = row["asset_name"]
    if name not in cache:
        cache[name] = md5(local_path(row))
    return cache[name]


def load_md5_cache(rows: list[dict[str, str]]) -> dict[str, str]:
    if not MD5_CACHE_PATH.exists():
        return {}
    try:
        saved = json.loads(MD5_CACHE_PATH.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    cache: dict[str, str] = {}
    for row in rows:
        path = local_path(row)
        entry = saved.get(row["asset_name"], {})
        if (
            int(entry.get("size", -1)) == path.stat().st_size
            and int(entry.get("mtime_ns", -1)) == path.stat().st_mtime_ns
            and entry.get("md5")
        ):
            cache[row["asset_name"]] = str(entry["md5"])
    return cache


def save_md5_cache(rows: list[dict[str, str]], cache: dict[str, str]) -> None:
    saved = {}
    for row in rows:
        name = row["asset_name"]
        if name not in cache:
            continue
        path = local_path(row)
        saved[name] = {
            "size": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
            "md5": cache[name],
        }
    MD5_CACHE_PATH.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")


def remote_matches(
    item: dict[str, Any] | None,
    row: dict[str, str],
    checksum_cache: dict[str, str],
) -> bool:
    if not item or int(item["filesize"]) != int(row["size_bytes"]):
        return False
    remote_checksum = str(item.get("checksum", "")).removeprefix("md5:")
    return bool(remote_checksum) and remote_checksum == local_md5(row, checksum_cache)


def prune_remote_files(
    rows: list[dict[str, str]],
    deposition: dict[str, Any],
    access_token: str,
) -> None:
    expected = {row["asset_name"] for row in rows}
    for name, item in remote_files(deposition).items():
        if name in expected:
            continue
        delete_url = item.get("links", {}).get("self")
        if not delete_url:
            raise RuntimeError(f"No delete URL for unexpected draft file: {name}")
        response = request_with_retry(
            "delete", delete_url, headers=auth_headers(access_token), timeout=120
        )
        check_response(response)
        print(f"Removed unlisted draft file: {name}")


def verify_remote(
    rows: list[dict[str, str]],
    deposition: dict[str, Any],
    checksum_cache: dict[str, str],
) -> None:
    files = remote_files(deposition)
    errors = []
    for row in rows:
        item = files.get(row["asset_name"])
        if item is None:
            errors.append(f"missing {row['asset_name']}")
        elif int(item["filesize"]) != int(row["size_bytes"]):
            errors.append(f"size mismatch {row['asset_name']}: {item['filesize']} != {row['size_bytes']}")
        elif not remote_matches(item, row, checksum_cache):
            errors.append(f"MD5 mismatch {row['asset_name']}")
    if errors:
        raise RuntimeError("Remote verification failed: " + "; ".join(errors))


def update_readme_reference(deposition: dict[str, Any]) -> None:
    record_id = str(deposition["id"])
    doi = (
        deposition.get("doi")
        or deposition.get("metadata", {}).get("doi")
        or deposition.get("metadata", {}).get("prereserve_doi", {}).get("doi", "")
    )
    concept_doi = deposition.get("conceptdoi", "")
    if not concept_doi and deposition.get("conceptrecid"):
        concept_doi = f"10.5281/zenodo.{deposition['conceptrecid']}"
    readme = BASE / "README.md"
    content = readme.read_text(encoding="utf-8")
    start = "<!-- ZENODO_START -->"
    end = "<!-- ZENODO_END -->"
    if doi:
        concept_line = (
            f"All releases are linked by the concept DOI [`{concept_doi}`](https://doi.org/{concept_doi}).\n"
            if concept_doi
            else ""
        )
        block = (
            f"{start}\n"
            f"[![DOI](https://zenodo.org/badge/DOI/{doi}.svg)](https://doi.org/{doi})\n\n"
            f"**Large assets:** [Zenodo record {record_id}](https://zenodo.org/records/{record_id}) "
            f"(version DOI: `{doi}`). File checksums and restore paths are in `zenodo_assets_manifest.csv`.\n"
            f"{concept_line}{end}"
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


def update_public_references(rows: list[dict[str, str]], deposition: dict[str, Any]) -> None:
    record_id = str(deposition["id"])
    for row in rows:
        row["zenodo_record_id"] = record_id
        quoted = urllib.parse.quote(row["asset_name"])
        row["download_url"] = f"https://zenodo.org/records/{record_id}/files/{quoted}?download=1"

    with PUBLIC_MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    update_readme_reference(deposition)


def main() -> None:
    args = parse_args()
    access_token = token()
    deposition = get_or_create_deposition(access_token, args.new_version)
    deposition = update_deposition_metadata(deposition, access_token)
    if args.initialize:
        update_readme_reference(deposition)
        doi = deposition.get("metadata", {}).get("prereserve_doi", {}).get("doi", "")
        print(f"Initialized Zenodo draft {deposition['id']} with reserved DOI {doi}")
        return

    rows = load_rows()
    bucket_url = deposition["links"]["bucket"]
    if args.prune:
        prune_remote_files(rows, deposition, access_token)
        response = request_with_retry(
            "get", f"{API_ROOT}/{deposition['id']}", headers=auth_headers(access_token), timeout=60
        )
        deposition = check_response(response).json()
    existing = remote_files(deposition)
    checksum_cache = load_md5_cache(rows)

    if not args.skip_upload:
        pending = []
        for row in rows:
            item = existing.get(row["asset_name"])
            if remote_matches(item, row, checksum_cache) and not args.force_upload:
                print(f"Remote file already matches by name, size, and MD5: {row['asset_name']}")
                continue
            pending.append(row)
        save_md5_cache(rows, checksum_cache)
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

    response = request_with_retry(
        "get", f"{API_ROOT}/{deposition['id']}", headers=auth_headers(access_token), timeout=60
    )
    deposition = check_response(response).json()
    verify_remote(rows, deposition, checksum_cache)
    save_md5_cache(rows, checksum_cache)
    update_public_references(rows, deposition)

    if args.publish:
        if deposition.get("submitted"):
            print(f"Zenodo deposition {deposition['id']} is already published")
        else:
            draft_response = request_with_retry(
                "get",
                f"https://zenodo.org/api/records/{deposition['id']}/draft",
                headers=auth_headers(access_token),
                timeout=120,
            )
            draft = check_response(draft_response).json()
            publish_url = draft.get("links", {}).get("publish")
            if not publish_url:
                raise RuntimeError(f"Zenodo draft {deposition['id']} has no Records API publish link")
            response = request_with_retry(
                "post", publish_url, headers=auth_headers(access_token), timeout=600
            )
            check_response(response)
            response = request_with_retry(
                "get",
                f"{API_ROOT}/{deposition['id']}",
                headers=auth_headers(access_token),
                timeout=120,
            )
            deposition = check_response(response).json()
            update_public_references(rows, deposition)
            print(f"Published Zenodo record {deposition['id']}")
    else:
        print(f"Zenodo draft {deposition['id']} verified; rerun with --publish to make it public")


if __name__ == "__main__":
    main()
