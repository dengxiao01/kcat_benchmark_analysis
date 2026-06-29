#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3: Add GPR (gene-protein-reaction) column from eciML1515.json
and keep only reactions that have a valid GPR.

Inputs (in project root):
  - kcat_comparison.csv
  - eciML1515.json   (COBRA model with reactions[].id and gene_reaction_rule)

Output (in project root):
  - kcat_comparison_with_gpr.csv
"""

import json
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_SOURCES = ['GO_HKP', 'DLKcat', 'MTLKP', 'TurNup', 'UniKP']


def main():
    print("=" * 70)
    print("Step 3: Add GPR column from model and filter")
    print("=" * 70)

    model_path = BASE / 'eciML1515.json'
    print(f"Reading {model_path.name}...")
    with open(model_path, 'r', encoding='utf-8') as f:
        model = json.load(f)

    reaction_gpr = {}
    for rxn in model['reactions']:
        gpr = rxn.get('gene_reaction_rule', '')
        if gpr and isinstance(gpr, str) and gpr.strip():
            reaction_gpr[rxn['id']] = gpr.strip()
    print(f"  - Found {len(reaction_gpr)} reactions with valid GPR")

    src = BASE / 'kcat_comparison.csv'
    print(f"\nReading {src.name}...")
    df = pd.read_csv(src)
    print(f"  - Original: {len(df)} reactions")

    df_f = df[df['reaction'].isin(reaction_gpr.keys())].copy()
    print(f"  - After filter: {len(df_f)} reactions")

    df_f['gpr'] = df_f['reaction'].map(reaction_gpr)

    cols = ['reaction', 'gpr'] + [c for c in df_f.columns if c not in ('reaction', 'gpr')]
    df_f = df_f[cols]

    out = BASE / 'kcat_comparison_with_gpr.csv'
    print(f"\nWriting {out.name}...")
    df_f.to_csv(out, index=False, encoding='utf-8')

    print(f"\n✅ Done!")
    print(f"  - Original: {len(df)}")
    print(f"  - With valid GPR: {len(df_f)} ({len(df_f) / len(df) * 100:.1f}%)")
    print(f"  - Filtered out: {len(df) - len(df_f)}")

    print("\nSource coverage (after filter):")
    for s in DATA_SOURCES:
        n = df_f[s].notna().sum()
        print(f"  - {s}: {n} ({n / len(df_f) * 100:.1f}%)")

    n_simple = (df_f['gpr'].str.contains('and|or', case=False, na=False) == False).sum()
    n_complex = df_f['gpr'].str.contains('and|or', case=False, na=False).sum()
    print(f"\nGPR complexity:")
    print(f"  - Simple (single gene): {n_simple}")
    print(f"  - Complex (contains and/or): {n_complex}")


if __name__ == '__main__':
    main()
