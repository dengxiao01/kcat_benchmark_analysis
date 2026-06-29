#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 4: Append `database` and `fill_method` columns to the comparison CSV.

Input (in project root):
  - kcat_comparison_with_gpr.csv
  - reaction_kcat_MW_databasefill.csv
        cols: kcat, data_type   (index = reaction id)

Output (in project root):
  - kcat_comparison_with_gpr.csv  (overwritten with two extra columns)

Mapping rules (mirroring update_database_fill_v2.py):
  - data_type in {'fill', 'Other_species'}  -> fill_method
  - data_type in {'Database', 'Brenda_SA', ...}  -> database
  - Some rows may have both; the script preserves both columns.
"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA_SOURCES = ['GO_HKP', 'DLKcat', 'MTLKP', 'TurNup', 'UniKP', 'database', 'fill_method']


def main():
    print("=" * 70)
    print("Step 4: Add database and fill_method columns")
    print("=" * 70)

    db_path = BASE / 'reaction_kcat_MW_databasefill.csv'
    print(f"Reading {db_path.name}...")
    df_db = pd.read_csv(db_path, index_col=0)
    print(f"  - Rows: {len(df_db)}")
    print(f"\ndata_type distribution:")
    print(df_db['data_type'].value_counts().to_string())

    database_map, fill_map = {}, {}
    for rxn_id, row in df_db.iterrows():
        kcat, dtype = row['kcat'], row['data_type']
        if dtype in ('fill', 'Other_species'):
            fill_map[rxn_id] = kcat
        else:  # Database, Brenda_SA, ...
            database_map[rxn_id] = kcat
    print(f"\n  - Database type (Database + Brenda_SA...): {len(database_map)}")
    print(f"  - Fill type (fill + Other_species): {len(fill_map)}")

    main_path = BASE / 'kcat_comparison_with_gpr.csv'
    print(f"\nReading {main_path.name}...")
    df = pd.read_csv(main_path)
    for col in ('database', 'fill_method'):
        if col in df.columns:
            df = df.drop(columns=[col])
    print(f"  - Rows: {len(df)}")

    df['database'] = df['reaction'].map(database_map)
    df['fill_method'] = df['reaction'].map(fill_map)

    n_db = df['database'].notna().sum()
    n_fm = df['fill_method'].notna().sum()
    print(f"\n  - database matches: {n_db}")
    print(f"  - fill_method matches: {n_fm}")

    print(f"\nWriting {main_path.name}...")
    df.to_csv(main_path, index=False, encoding='utf-8')

    print(f"\n✅ Done!")
    print(f"  - Total: {len(df)}")
    print(f"  - Has database: {n_db} ({n_db / len(df) * 100:.1f}%)")
    print(f"  - Has fill_method: {n_fm} ({n_fm / len(df) * 100:.1f}%)")
    print(f"  - Both: {((df['database'].notna()) & (df['fill_method'].notna())).sum()}")
    print(f"  - At least one: {((df['database'].notna()) | (df['fill_method'].notna())).sum()}")

    print("\nSource coverage:")
    for s in DATA_SOURCES:
        n = df[s].notna().sum()
        print(f"  - {s}: {n} ({n / len(df) * 100:.1f}%)")


if __name__ == '__main__':
    main()
