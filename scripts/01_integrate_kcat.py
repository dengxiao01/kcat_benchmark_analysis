#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1: Integrate all 5 kcat prediction sources into a single JSON file.

Inputs (in project root):
  - GO_HKP.json                         (dict: reaction_id -> kcat value)
  - reaction_kcat_MW_DLKcat.csv         (cols: reactions, kcat)
  - reaction_kcat_MW_MTLKP.csv          (cols: reactions, kcat)
  - reaction_kcat_MW_TurNup.csv         (cols: reactions, kcat)
  - reaction_kcat_MW_UniKP.csv          (cols: reactions, kcat)

Output (in project root):
  - integrated_kcat_simple.json
        list of {"reaction": <id>, "kcat": {<source>: <value>, ...}}
"""

import json
from collections import defaultdict
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_SOURCES = ['GO_HKP', 'DLKcat', 'MTLKP', 'TurNup', 'UniKP']


def read_go_hkp(file_path: Path) -> dict:
    """Read GO_HKP.json: {reaction_id: kcat_value}."""
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return {rxn: {'GO_HKP': val} for rxn, val in data.items()}


def read_csv_data(file_path: Path, source_name: str) -> dict:
    """Read a kcat CSV (cols: reactions, kcat) into {reaction: {source: kcat}}."""
    df = pd.read_csv(file_path)
    result = {}
    for _, row in df.iterrows():
        rxn = row.get('reactions', None)
        val = row.get('kcat', None)
        if rxn and pd.notna(val):
            result.setdefault(rxn, {})[source_name] = val
    return result


def merge_data(data_sources: dict) -> dict:
    merged = defaultdict(dict)
    for source_data in data_sources.values():
        for rxn, kcat_info in source_data.items():
            merged[rxn].update(kcat_info)
    return dict(merged)


def main():
    print("=" * 70)
    print("Step 1: Integrate kcat prediction sources")
    print("=" * 70)

    data_sources = {}
    for source in DATA_SOURCES:
        if source == 'GO_HKP':
            path = BASE / 'GO_HKP.json'
            print(f"  - Reading {path.name} (GO_HKP)...")
            data_sources[source] = read_go_hkp(path)
        else:
            path = BASE / f'reaction_kcat_MW_{source}.csv'
            print(f"  - Reading {path.name}...")
            data_sources[source] = read_csv_data(path, source)
        print(f"    Found {len(data_sources[source])} reactions")

    print("\nMerging data sources...")
    merged = merge_data(data_sources)

    result = [{'reaction': rxn, 'kcat': info} for rxn, info in sorted(merged.items())]

    out_path = BASE / 'integrated_kcat_simple.json'
    print(f"\nWriting {out_path}...")
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print(f"\n✅ Done! Integrated {len(result)} unique reactions.")

    source_counts = defaultdict(int)
    multi_source = 0
    for entry in result:
        if len(entry['kcat']) > 1:
            multi_source += 1
        for src in entry['kcat']:
            source_counts[src] += 1

    print("\nSource coverage:")
    for src, cnt in sorted(source_counts.items()):
        print(f"  - {src}: {cnt} reactions")
    print(f"  - Multi-source reactions: {multi_source}")


if __name__ == '__main__':
    main()
