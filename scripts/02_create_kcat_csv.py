#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2: Convert integrated_kcat_simple.json to a wide CSV.

Output (in project root):
  - kcat_comparison.csv
        columns: reaction, GO_HKP, DLKcat, MTLKP, TurNup, UniKP
"""

import json
from collections import defaultdict
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_SOURCES = ['GO_HKP', 'DLKcat', 'MTLKP', 'TurNup', 'UniKP']


def main():
    print("=" * 70)
    print("Step 2: Convert integrated JSON to wide CSV")
    print("=" * 70)

    src = BASE / 'integrated_kcat_simple.json'
    print(f"Reading {src.name}...")
    with open(src, 'r', encoding='utf-8') as f:
        kcat_data = json.load(f)
    print(f"  - Loaded {len(kcat_data)} reactions")

    print("\nBuilding rows...")
    rows = []
    for entry in kcat_data:
        row = {'reaction': entry['reaction']}
        for source in DATA_SOURCES:
            row[source] = entry['kcat'].get(source)
        rows.append(row)

    df = pd.DataFrame(rows).sort_values('reaction')

    out = BASE / 'kcat_comparison.csv'
    print(f"\nWriting {out.name}...")
    df.to_csv(out, index=False, encoding='utf-8')

    print(f"\n✅ Done!")
    print(f"  - Total rows (reactions): {len(df)}")
    print(f"  - Columns: {len(df.columns)} (1 id + {len(DATA_SOURCES)} sources)")

    print("\nSource coverage:")
    for source in DATA_SOURCES:
        n = df[source].notna().sum()
        print(f"  - {source}: {n} ({n / len(df) * 100:.1f}%)")

    multi = (df[DATA_SOURCES].notna().sum(axis=1) > 1).sum()
    print(f"\n  - Multi-source reactions: {multi}")

    print("\nSample (first 5 rows):")
    print(df.head().to_string())


if __name__ == '__main__':
    main()
