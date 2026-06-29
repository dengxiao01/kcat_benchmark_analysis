#!/usr/bin/env python3
"""Download, verify, and optionally restore large benchmark assets from Zenodo."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import tarfile
import urllib.request
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
MANIFEST = BASE / "zenodo_assets_manifest.csv"
DOWNLOAD_DIR = BASE / "release" / "downloads"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List available assets and exit.")
    parser.add_argument(
        "--asset",
        action="append",
        help="Bundle or part filename to download. Repeat for multiple bundles; use 'all' for every bundle.",
    )
    parser.add_argument("--restore", action="store_true", help="Restore downloaded files to manifest targets.")
    parser.add_argument("--force", action="store_true", help="Replace existing downloads and restored files.")
    return parser.parse_args()


def load_manifest() -> list[dict[str, str]]:
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Asset manifest not found: {MANIFEST}")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path: Path, row: dict[str, str]) -> None:
    expected_size = int(row["size_bytes"])
    if path.stat().st_size != expected_size:
        raise RuntimeError(f"Size mismatch for {path.name}: {path.stat().st_size} != {expected_size}")
    observed = sha256(path)
    if observed != row["sha256"]:
        raise RuntimeError(f"SHA256 mismatch for {path.name}: {observed} != {row['sha256']}")


def verify_bundle(path: Path, row: dict[str, str]) -> None:
    expected_size = int(row["bundle_size_bytes"])
    if path.stat().st_size != expected_size:
        raise RuntimeError(f"Bundle size mismatch for {path.name}: {path.stat().st_size} != {expected_size}")
    observed = sha256(path)
    if observed != row["bundle_sha256"]:
        raise RuntimeError(f"Bundle SHA256 mismatch for {path.name}: {observed} != {row['bundle_sha256']}")


def download(row: dict[str, str], force: bool) -> Path:
    url = row.get("download_url", "").strip()
    if not url:
        raise RuntimeError(f"No Zenodo download URL is recorded for {row['asset_name']}")
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    target = DOWNLOAD_DIR / row["asset_name"]
    if target.exists() and not force:
        print(f"Verifying existing download: {target}")
        verify(target, row)
        return target

    partial = target.with_suffix(target.suffix + ".part")
    if partial.exists():
        partial.unlink()
    print(f"Downloading {row['asset_name']} ({int(row['size_bytes']) / 1024**3:.2f} GiB)")
    with urllib.request.urlopen(url) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=16 * 1024 * 1024)
    verify(partial, row)
    partial.replace(target)
    return target


def assemble_bundle(rows: list[dict[str, str]], parts: list[Path], force: bool) -> Path:
    first = rows[0]
    if len(parts) == 1 and parts[0].name == first["bundle_name"]:
        verify_bundle(parts[0], first)
        return parts[0]

    output = DOWNLOAD_DIR / first["bundle_name"]
    if output.exists() and not force:
        print(f"Verifying existing assembled bundle: {output}")
        verify_bundle(output, first)
        return output

    partial = output.with_suffix(output.suffix + ".assembling")
    if partial.exists():
        partial.unlink()
    print(f"Assembling {len(parts)} parts into {output.name}")
    with partial.open("wb") as dst:
        for part in parts:
            with part.open("rb") as src:
                shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
    verify_bundle(partial, first)
    partial.replace(output)
    return output


def validated_members(tar: tarfile.TarFile, destination: Path) -> list[tarfile.TarInfo]:
    destination = destination.resolve()
    members = tar.getmembers()
    for member in members:
        candidate = (destination / member.name).resolve()
        if destination not in candidate.parents and candidate != destination:
            raise RuntimeError(f"Unsafe archive path: {member.name}")
        if member.issym() or member.islnk():
            raise RuntimeError(f"Archive links are not restored automatically: {member.name}")
    return members


def extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as tar:
        tar.extractall(destination, members=validated_members(tar, destination))


def restore(downloaded: Path, row: dict[str, str], force: bool) -> None:
    action = row["restore_action"]
    target = (BASE / row["restore_target"]).resolve()
    if BASE.resolve() not in target.parents and target != BASE.resolve():
        raise RuntimeError(f"Restore target escapes repository: {target}")

    if action == "extract_to_repo_root":
        print(f"Extracting {downloaded.name} into {BASE}")
        extract(downloaded, BASE)
        return

    if target.exists() and not force:
        if target.is_file() and sha256(target) == sha256(downloaded):
            print(f"Already restored: {target}")
        else:
            raise FileExistsError(f"Restore target exists; use --force to replace it: {target}")
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(downloaded, target)
        print(f"Placed {downloaded.name} at {target.relative_to(BASE)}")

    if action == "place_file_then_extract":
        print(f"Extracting {target.name} into {target.parent}")
        extract(target, target.parent)
    elif action != "place_file":
        raise RuntimeError(f"Unknown restore action: {action}")


def main() -> None:
    args = parse_args()
    rows = load_manifest()
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["bundle_name"], []).append(row)
    for group in groups.values():
        group.sort(key=lambda row: int(row["part_index"]))

    if args.list:
        for bundle_name, group in groups.items():
            row = group[0]
            print(
                f"{bundle_name}\t{int(row['bundle_size_bytes']) / 1024**3:.2f} GiB\t"
                f"{len(group)} parts\t{row['required_for']}\t{row['restore_target']}"
            )
        return

    requested = args.asset or []
    if not requested:
        raise SystemExit("Choose --list or provide --asset NAME (use --asset all for every asset).")
    known_names = set(groups) | {row["asset_name"] for row in rows}
    missing = set(requested) - known_names - {"all"}
    if missing:
        raise SystemExit(f"Unknown assets: {', '.join(sorted(missing))}")
    if "all" in requested:
        selected_bundles = list(groups)
    else:
        selected_bundles = []
        for name in requested:
            bundle = name if name in groups else next(row["bundle_name"] for row in rows if row["asset_name"] == name)
            if bundle not in selected_bundles:
                selected_bundles.append(bundle)

    for bundle in selected_bundles:
        group = groups[bundle]
        downloaded = [download(row, args.force) for row in group]
        if args.restore:
            assembled = assemble_bundle(group, downloaded, args.force)
            restore(assembled, group[0], args.force)


if __name__ == "__main__":
    main()
