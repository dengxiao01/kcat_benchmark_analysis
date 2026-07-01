#!/usr/bin/env python3
"""Build the reproducibility asset bundles intended for Zenodo."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import tarfile
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
RELEASE_DIR = BASE / "release" / "zenodo"
PARTS_DIR = RELEASE_DIR / "parts"
PUBLIC_MANIFEST = BASE / "zenodo_assets_manifest.csv"
PART_SIZE = 128 * 1024 * 1024

TURNUP_ROOT = (
    BASE
    / "external_methods"
    / "AI_file"
    / "turnup"
    / "kcat_prediction_function-main"
    / "kcat_prediction_function-main"
)

SKIP_NAMES = {
    ".git",
    "__pycache__",
    ".DS_Store",
    ".env",
    ".env.vercel",
    "zenodo.txt",
    "esm1b_t33_650M_UR50S.pt",
}
NOTICE = BASE / "THIRD_PARTY_NOTICES.md"
METHOD_SOURCES = BASE / "external_methods" / "METHOD_SOURCES.md"
PUBLIC_METHOD_DIRS = {
    "catapro",
    "catpred",
    "dekp",
    "dlkcat",
    "go_hkp",
    "kcatnet",
    "kinform",
    "pmak",
    "pretkcat",
    "selfprot",
    "turnup",
    "unikp",
}



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare selected large kcat benchmark assets for Zenodo.")
    parser.add_argument("--force", action="store_true", help="Rebuild archives even when they already exist.")
    return parser.parse_args()


def filter_tar_info(info: tarfile.TarInfo) -> tarfile.TarInfo | None:
    parts = Path(info.name).parts
    if any(part in SKIP_NAMES for part in parts):
        return None
    return info


def add_path(tar: tarfile.TarFile, path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)
    tar.add(path, arcname=path.relative_to(BASE).as_posix(), recursive=True, filter=filter_tar_info)


def build_archive(output: Path, paths: list[Path], mode: str, force: bool) -> None:
    if output.exists() and not force:
        print(f"Using existing archive: {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(output, mode) as tar:
        for path in paths:
            print(f"Adding to {output.name}: {path.relative_to(BASE)}", flush=True)
            add_path(tar, path)


def build_catpred_archive(output: Path, force: bool) -> None:
    if output.exists() and not force:
        print(f"Using existing archive: {output}")
        return
    source = BASE / "external_methods" / "CatPred_capsule" / "capsule_data_update.tar.gz"
    prefix = "data/pretrained/production/kcat"
    output.parent.mkdir(parents=True, exist_ok=True)
    selected = 0
    with tarfile.open(source, "r:gz") as src, tarfile.open(output, "w:gz") as dst:
        for original in src:
            if original.name != prefix and not original.name.startswith(prefix + "/"):
                continue
            member = copy.copy(original)
            member.name = f"external_methods/CatPred_capsule/{original.name}"
            fileobj = src.extractfile(original) if original.isfile() else None
            dst.addfile(member, fileobj)
            selected += 1
        add_path(dst, NOTICE)
        add_path(dst, METHOD_SOURCES)
    if not selected:
        raise RuntimeError(f"No CatPred production kcat assets found in: {source}")
    print(f"Selected {selected} CatPred production kcat entries into {output.name}")


def core_result_paths() -> list[Path]:
    paths = [
        BASE / "data" / "final" / "experimental_kcat_truth.csv",
        BASE / "data" / "final" / "benchmark_ready_truth.csv",
        BASE / "data" / "final" / "benchmark_ready_catpred.csv",
        BASE / "README.md",
        BASE / "LICENSE",
        NOTICE,
        METHOD_SOURCES,
        BASE / "reports" / "tables",
        BASE / "reports" / "report_tables",
        BASE / "reports" / "figures",
    ]
    for method_dir in sorted((BASE / "data" / "final").iterdir()):
        if not method_dir.is_dir() or method_dir.name not in PUBLIC_METHOD_DIRS:
            continue
        paths.extend(sorted(method_dir.rglob("*.csv")))
    return list(dict.fromkeys(paths))



def turnup_paths() -> list[Path]:
    return [TURNUP_ROOT / "code" / "data", NOTICE, METHOD_SOURCES]


def other_model_paths() -> list[Path]:
    return [
        BASE / "external_methods" / "CataPro" / "models" / "kcat_models",
        BASE / "external_methods" / "PMAK" / "code" / "save_model",
        BASE / "external_methods" / "KcatNet" / "RESULT" / "model_KcatNet.pt",
        BASE / "external_methods" / "PreTKcat" / "MolGNet.pt",
        BASE / "data" / "final" / "pretkcat" / "pretkcat_extratrees_model.pkl",
        BASE / "external_methods" / "UniKP_official" / "models" / "UniKP for kcat.pkl",
        BASE / "external_methods" / "SELFprot" / "weights",
        BASE / "data" / "final" / "dekp" / "dekp_public_retrained_model.pt",
        NOTICE,
        METHOD_SOURCES,
    ]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def split_asset(bundle_name: str, source: Path, force: bool) -> list[Path]:
    if source.stat().st_size <= PART_SIZE:
        return [source]

    PARTS_DIR.mkdir(parents=True, exist_ok=True)
    existing = sorted(PARTS_DIR.glob(f"{bundle_name}.part[0-9][0-9][0-9]"))
    if existing and not force and sum(path.stat().st_size for path in existing) == source.stat().st_size:
        print(f"Using {len(existing)} existing parts for: {bundle_name}")
        return existing
    for path in existing:
        path.unlink()

    parts = []
    with source.open("rb") as src:
        index = 1
        while chunk := src.read(PART_SIZE):
            part = PARTS_DIR / f"{bundle_name}.part{index:03d}"
            with part.open("wb") as dst:
                dst.write(chunk)
            parts.append(part)
            print(f"Wrote part {index}: {part.name} ({len(chunk)} bytes)", flush=True)
            index += 1
    return parts


def asset_rows(
    bundle_name: str,
    local_path: Path,
    content: str,
    restore_action: str,
    restore_target: str,
    required_for: str,
    force: bool,
) -> list[dict[str, str | int]]:
    if not local_path.exists():
        raise FileNotFoundError(local_path)
    parts = split_asset(bundle_name, local_path, force)
    bundle_sha256 = sha256(local_path)
    rows = []
    for index, part in enumerate(parts, start=1):
        rows.append(
            {
                "asset_name": bundle_name if len(parts) == 1 else part.name,
                "bundle_name": bundle_name,
                "part_index": index,
                "part_count": len(parts),
                "local_path": part.relative_to(BASE).as_posix(),
                "size_bytes": part.stat().st_size,
                "sha256": sha256(part),
                "bundle_size_bytes": local_path.stat().st_size,
                "bundle_sha256": bundle_sha256,
                "content": content,
                "restore_action": restore_action,
                "restore_target": restore_target,
                "required_for": required_for,
                "zenodo_record_id": "",
                "download_url": "",
            }
        )
    return rows


def write_readme(rows: list[dict[str, str | int]]) -> None:
    lines = [
        "# Zenodo Asset Bundle",
        "",
        "These files supplement the GitHub repository with large benchmark results and model assets.",
        "Archives preserve project-relative paths. Extract archive assets from the repository root.",
        "Large bundles are split into 128 MiB parts to avoid proxy timeouts.",
        "The download helper verifies each part, concatenates the original bundle, verifies it again, and restores it.",
        "Third-party model assets retain their upstream licenses and citation requirements.",
        "",
        "| bundle | bundle_size_bytes | parts | restore_action | restore_target | required_for |",
        "| --- | ---: | ---: | --- | --- | --- |",
    ]
    seen = set()
    for row in rows:
        if row["bundle_name"] in seen:
            continue
        seen.add(row["bundle_name"])
        lines.append(
            f"| {row['bundle_name']} | {row['bundle_size_bytes']} | {row['part_count']} | "
            f"{row['restore_action']} | "
            f"{row['restore_target']} | {row['required_for']} |"
        )
    (RELEASE_DIR / "ASSET_README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    RELEASE_DIR.mkdir(parents=True, exist_ok=True)

    core = RELEASE_DIR / "kcat_benchmark_core_data_and_results.tar.gz"
    catpred = RELEASE_DIR / "kcat_benchmark_catpred_kcat_assets.tar.gz"
    turnup = RELEASE_DIR / "kcat_benchmark_turnup_kcat_assets.tar"
    other = RELEASE_DIR / "kcat_benchmark_other_model_assets.tar"

    build_archive(core, core_result_paths(), "w:gz", args.force)
    build_catpred_archive(catpred, args.force)
    build_archive(turnup, turnup_paths(), "w", args.force)
    build_archive(other, other_model_paths(), "w", args.force)

    bundles = [
        (
            core.name,
            core,
            "Unified benchmark CSVs, per-method evaluated tables, reports, figures, and initial analyses.",
            "extract_to_repo_root",
            ".",
            "benchmark inspection and report regeneration",
        ),
        (
            catpred.name,
            catpred,
            "CatPred production kcat ensemble checkpoints selected from the official capsule.",
            "extract_to_repo_root",
            ".",
            "CatPred",
        ),
        (
            "catpred_db.tar.gz",
            BASE / "external_methods" / "CatPred_datas" / "catpred-db.tar.gz",
            "CatPred-DB dataset bundle used for dataset comparison.",
            "place_file",
            "external_methods/CatPred_datas/catpred-db.tar.gz",
            "CatPred dataset comparison",
        ),
        (
            "kinform_results.tar.gz",
            BASE / "external_methods" / "KinForm" / "results.tar.gz",
            "KinForm trained models and cached result assets.",
            "place_file_then_extract",
            "external_methods/KinForm/results.tar.gz",
            "KinForm",
        ),
        (
            turnup.name,
            turnup,
            "TurNuP task model, XGBoost, molecular files, and reaction assets; download base ESM1b separately.",
            "extract_to_repo_root",
            ".",
            "TurNuP-official",
        ),
        (
            other.name,
            other,
            "CataPro, PMAK, KcatNet, PreTKcat, UniKP, SELFprot, and public-retrained DEKP model assets.",
            "extract_to_repo_root",
            ".",
            "remaining benchmark methods",
        ),
    ]
    rows = []
    for bundle in bundles:
        rows.extend(asset_rows(*bundle, force=args.force))

    fieldnames = list(rows[0])
    for path in [RELEASE_DIR / "assets_manifest.csv", PUBLIC_MANIFEST]:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    write_readme(rows)
    print(f"Prepared {len(bundles)} bundles as {len(rows)} Zenodo files under: {RELEASE_DIR}")


if __name__ == "__main__":
    main()
