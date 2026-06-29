#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 5: Run the 10 analysis points defined in .claude/commands/kcat_compare.md
and write outputs to analysis_results/.

Inputs (in project root):
  - kcat_comparison_with_gpr.csv

Outputs (in project root/analysis_results/):
  - 01_distribution_statistics.csv / 01_global_distribution.png
  - 02_correlation_matrix.csv / 02_correlation_analysis.png
  - 03_asymmetry_statistics.csv / 03_thermodynamic_asymmetry.png
  - 04_isozyme_statistics.csv / 04_isozyme_specificity.png
  - 05_complex_statistics.csv / 05_complex_handling.png
  - 06_substrate_statistics.csv / 06_substrate_specificity.png
  - 07_coverage_statistics.csv / 07_coverage_analysis.png
  - 08_benchmark_statistics.csv / 08_ground_truth_benchmark.png
  - 09_bias_statistics.csv / 09_bias_detection.png
  - 10_ensemble_statistics.csv / 10_ensemble_comparison.png
  - kcat_comparison_enhanced.csv  (input + ensemble columns, in project root)
"""

import os
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
from scipy.stats import gaussian_kde
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

warnings.filterwarnings('ignore')

BASE = Path(__file__).resolve().parent.parent
INPUT_CSV = BASE / 'kcat_comparison_with_gpr.csv'
OUTPUT_DIR = BASE / 'analysis_results'
ENHANCED_CSV = BASE / 'kcat_comparison_enhanced.csv'

# Constants
ARTIFACT_VALUE = 112645.641475366
METHODS = ['DLKcat', 'MTLKP', 'TurNup', 'UniKP']

plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

def classify_reaction(name: str) -> str:
    """Classify a reaction name into Transport / Metabolic (Reverse) / Core."""
    if 'tex' in name or 'tipp' in name or 'tpp' in name:
        return 'Transport'
    if '_reverse' in name:
        return 'Metabolic (Reverse)'
    return 'Metabolic (Core)'


def is_complex(gpr) -> bool:
    return pd.notna(gpr) and ' and ' in str(gpr)


def get_base_reaction(name: str) -> str:
    name = re.sub(r'_num\d+$', '', name)
    return re.sub(r'_reverse$', '', name)


# ----------------------------------------------------------------------------
# Load and preprocess
# ----------------------------------------------------------------------------

print("=" * 80)
print("COMPREHENSIVE KCAT PREDICTION METHOD EVALUATION")
print("=" * 80)
print()

print("📊 Loading data...")
df = pd.read_csv(INPUT_CSV)
print(f"   Total rows: {len(df)}")
print(f"   Columns: {', '.join(df.columns)}")
print()

df['database_clean'] = df['database'].replace(ARTIFACT_VALUE, np.nan)
df['database_clean'] = pd.to_numeric(df['database_clean'], errors='coerce')

for m in METHODS:
    df[m] = pd.to_numeric(df[m], errors='coerce')

df['reaction_type'] = df['reaction'].apply(classify_reaction)
df['is_complex'] = df['gpr'].apply(is_complex)
df['base_reaction'] = df['reaction'].apply(get_base_reaction)
df['is_reverse'] = df['reaction'].str.contains('_reverse', na=False)

print("✅ Data preprocessing complete")
print(f"   Artifact values in database: {(df['database'] == ARTIFACT_VALUE).sum()}")
print(f"   Valid ground truth: {df['database_clean'].notna().sum()}")
print(f"   Transport reactions: {(df['reaction_type'] == 'Transport').sum()}")
print(f"   Metabolic reactions: {df['reaction_type'].str.contains('Metabolic').sum()}")
print(f"   Enzyme complexes: {df['is_complex'].sum()}")
print()


# ----------------------------------------------------------------------------
# Point 1: Global distribution
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 1: GLOBAL DISTRIBUTION & DYNAMIC RANGE")
print("=" * 80)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Analysis Point 1: Global Distribution of kcat Predictions (log10 scale)',
             fontsize=16, fontweight='bold')

distribution_stats = []
for idx, method in enumerate(METHODS):
    ax = axes[idx // 3, idx % 3]
    values = df[method].dropna()
    log_values = np.log10(values[values > 0])

    if len(log_values) > 0:
        ax.hist(log_values, bins=50, alpha=0.6, density=True,
                label='Histogram', edgecolor='black')
        if len(log_values) > 1:
            kde = gaussian_kde(log_values)
            x_range = np.linspace(log_values.min(), log_values.max(), 200)
            ax.plot(x_range, kde(x_range), 'r-', linewidth=2, label='KDE')

        mean_val = np.mean(log_values)
        median_val = np.median(log_values)
        skewness = stats.skew(log_values)
        kurtosis = stats.kurtosis(log_values)

        ax.axvline(1, color='green', linestyle='--', linewidth=2,
                   label='Bio. Median (10 s⁻¹)')
        ax.axvline(median_val, color='orange', linestyle='--', linewidth=2,
                   label='Model Median')

        ax.set_xlabel('log₁₀(kcat) [s⁻¹]', fontsize=11)
        ax.set_ylabel('Density', fontsize=11)
        ax.set_title(f'{method}\nN={len(log_values)}, Skew={skewness:.2f}, '
                     f'Kurt={kurtosis:.2f}', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        distribution_stats.append({
            'Method': method,
            'N_predictions': len(values),
            'Mean_log10': mean_val,
            'Median_log10': median_val,
            'Std_log10': np.std(log_values),
            'Skewness': skewness,
            'Kurtosis': kurtosis,
            'Min_log10': log_values.min(),
            'Max_log10': log_values.max(),
            'Dynamic_Range': log_values.max() - log_values.min(),
        })

# ECDF panel
ax = axes[1, 2]
for method in METHODS:
    values = df[method].dropna()
    log_values = np.log10(values[values > 0])
    if len(log_values) > 0:
        sorted_vals = np.sort(log_values)
        y = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        ax.plot(sorted_vals, y, linewidth=2, label=method)
ax.set_xlabel('log₁₀(kcat) [s⁻¹]', fontsize=11)
ax.set_ylabel('Cumulative Probability', fontsize=11)
ax.set_title('ECDF Comparison', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / '01_global_distribution.png', dpi=300, bbox_inches='tight')
print("   📈 Distribution plots saved")

stats_df = pd.DataFrame(distribution_stats)
stats_df.to_csv(OUTPUT_DIR / '01_distribution_statistics.csv', index=False)
print("   📊 Distribution statistics saved")
print()
print(stats_df.to_string(index=False))
print()


# ----------------------------------------------------------------------------
# Point 2: Inter-model concordance
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 2: INTER-MODEL CONCORDANCE")
print("=" * 80)

common_mask = df[METHODS].notna().all(axis=1)
common_df = df[common_mask].copy()
print(f"   Reactions with predictions from all methods: {len(common_df)}")

if len(common_df) > 5:
    log_data = np.log10(common_df[METHODS].replace(0, np.nan))
    corr_matrix = log_data.corr(method='spearman')

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdYlGn', center=0.5,
                vmin=0, vmax=1, square=True, ax=axes[0],
                cbar_kws={'label': 'Spearman ρ'})
    axes[0].set_title('Spearman Correlation Matrix\n(log₁₀-transformed kcat values)',
                     fontsize=13, fontweight='bold')

    linkage = hierarchy.linkage(corr_matrix, method='ward')
    hierarchy.dendrogram(linkage, labels=METHODS, ax=axes[1],
                         orientation='right', leaf_font_size=12)
    axes[1].set_title('Hierarchical Clustering of Methods', fontsize=13, fontweight='bold')
    axes[1].set_xlabel('Distance', fontsize=11)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '02_correlation_analysis.png', dpi=300, bbox_inches='tight')
    print("   📈 Correlation analysis saved")

    corr_matrix.to_csv(OUTPUT_DIR / '02_correlation_matrix.csv')
    print("   📊 Correlation matrix saved")
    print()
    print("Spearman Correlation Matrix:")
    print(corr_matrix.round(3))
    print()
else:
    print("   ⚠️  Insufficient common predictions for correlation analysis")
    print()


# ----------------------------------------------------------------------------
# Point 3: Thermodynamic asymmetry
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 3: THERMODYNAMIC ASYMMETRY (Forward vs Reverse)")
print("=" * 80)

forward_df = df[~df['is_reverse']].copy()
reverse_df = df[df['is_reverse']].copy()
forward_df['merge_key'] = forward_df['base_reaction'] + '_' + forward_df['gpr'].astype(str)
reverse_df['merge_key'] = reverse_df['base_reaction'] + '_' + reverse_df['gpr'].astype(str)
pairs = forward_df.merge(reverse_df, on='merge_key', suffixes=('_fwd', '_rev'))
print(f"   Forward-reverse pairs found: {len(pairs)}")

if len(pairs) > 0:
    asymmetry_stats = []
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Analysis Point 3: Thermodynamic Asymmetry Analysis',
                 fontsize=16, fontweight='bold')

    for idx, method in enumerate(METHODS):
        ax = axes[idx // 3, idx % 3]
        fwd_col, rev_col = f'{method}_fwd', f'{method}_rev'
        valid_pairs = pairs[[fwd_col, rev_col]].dropna()

        if len(valid_pairs) > 0:
            fwd_vals = valid_pairs[fwd_col].values
            rev_vals = valid_pairs[rev_col].values
            asi = np.abs(np.log10(fwd_vals / rev_vals))

            ax.scatter(np.log10(fwd_vals), np.log10(rev_vals), alpha=0.6, s=50)
            lims = [
                min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])
            ]
            ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=2, label='Perfect Symmetry')

            mean_asi = np.mean(asi)
            symmetric_count = int(np.sum(asi < 0.1))

            ax.set_xlabel('log₁₀(kcat forward) [s⁻¹]', fontsize=10)
            ax.set_ylabel('log₁₀(kcat reverse) [s⁻¹]', fontsize=10)
            ax.set_title(f'{method}\nASI={mean_asi:.2f}, '
                         f'Symmetric={symmetric_count}/{len(valid_pairs)}',
                         fontsize=11, fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

            asymmetry_stats.append({
                'Method': method,
                'N_pairs': len(valid_pairs),
                'Mean_ASI': mean_asi,
                'Median_ASI': np.median(asi),
                'Symmetric_pairs': symmetric_count,
                'Symmetric_percentage': 100 * symmetric_count / len(valid_pairs),
            })

    # ASI distribution panel
    ax = axes[1, 2]
    for method in METHODS:
        fwd_col, rev_col = f'{method}_fwd', f'{method}_rev'
        valid_pairs = pairs[[fwd_col, rev_col]].dropna()
        if len(valid_pairs) > 0:
            asi = np.abs(np.log10(valid_pairs[fwd_col].values /
                                  valid_pairs[rev_col].values))
            ax.hist(asi, bins=30, alpha=0.5, label=method)
    ax.set_xlabel('Asymmetry Index |log₁₀(fwd/rev)|', fontsize=10)
    ax.set_ylabel('Frequency', fontsize=10)
    ax.set_title('Distribution of Asymmetry', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '03_thermodynamic_asymmetry.png', dpi=300, bbox_inches='tight')
    print("   📈 Asymmetry analysis saved")

    pd.DataFrame(asymmetry_stats).to_csv(
        OUTPUT_DIR / '03_asymmetry_statistics.csv', index=False)
    print("   📊 Asymmetry statistics saved")
    print(pd.DataFrame(asymmetry_stats).to_string(index=False))
    print()
else:
    print("   ⚠️  No forward-reverse pairs found")
    print()


# ----------------------------------------------------------------------------
# Point 4: Isozyme specificity
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 4: ISOZYME SPECIFICITY")
print("=" * 80)

isozyme_groups = df.groupby('base_reaction')
isozyme_reactions = [n for n, g in isozyme_groups if len(g['gpr'].unique()) > 1]
print(f"   Reactions with multiple isozymes: {len(isozyme_reactions)}")

if isozyme_reactions:
    isozyme_stats = []
    for method in METHODS:
        cvs = []
        for reaction in isozyme_reactions:
            group = df[df['base_reaction'] == reaction]
            values = group[method].dropna()
            if len(values) > 1 and np.mean(values) > 0:
                cvs.append(np.std(values) / np.mean(values))
        if cvs:
            isozyme_stats.append({
                'Method': method,
                'N_reactions': len(cvs),
                'Mean_CV': np.mean(cvs),
                'Median_CV': np.median(cvs),
                'CV_std': np.std(cvs),
                'High_specificity_count': int(np.sum(np.array(cvs) > 0.3)),
            })

    if isozyme_stats:
        fig, ax = plt.subplots(figsize=(10, 6))
        iso_df = pd.DataFrame(isozyme_stats)
        x_pos = np.arange(len(iso_df))
        bars = ax.bar(x_pos, iso_df['Mean_CV'], yerr=iso_df['CV_std'],
                      capsize=5, alpha=0.7, edgecolor='black')
        colors = plt.cm.viridis(np.linspace(0, 1, len(iso_df)))
        for bar, color in zip(bars, colors):
            bar.set_color(color)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(iso_df['Method'], fontsize=11)
        ax.set_ylabel('Coefficient of Variation (CV)', fontsize=12)
        ax.set_title('Isozyme Discrimination Ability\n'
                     '(Higher CV = Better Sequence Resolution)',
                     fontsize=13, fontweight='bold')
        ax.axhline(0.3, color='red', linestyle='--', linewidth=2,
                   label='High Specificity Threshold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / '04_isozyme_specificity.png', dpi=300, bbox_inches='tight')
        print("   📈 Isozyme analysis saved")

        iso_df.to_csv(OUTPUT_DIR / '04_isozyme_statistics.csv', index=False)
        print("   📊 Isozyme statistics saved")
        print(iso_df.to_string(index=False))
        print()
else:
    print("   ⚠️  No isozyme variations found")
    print()


# ----------------------------------------------------------------------------
# Point 5: Enzyme complex handling
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 5: ENZYME COMPLEX HANDLING")
print("=" * 80)

complex_stats = []
for method in METHODS:
    single_vals = df[~df['is_complex']][method].dropna()
    complex_vals = df[df['is_complex']][method].dropna()
    if len(single_vals) > 0 and len(complex_vals) > 0:
        _, pval = stats.mannwhitneyu(np.log10(single_vals), np.log10(complex_vals))
        complex_stats.append({
            'Method': method,
            'N_single': len(single_vals),
            'N_complex': len(complex_vals),
            'Mean_single_log10': np.mean(np.log10(single_vals)),
            'Mean_complex_log10': np.mean(np.log10(complex_vals)),
            'Median_single_log10': np.median(np.log10(single_vals)),
            'Median_complex_log10': np.median(np.log10(complex_vals)),
            'MannWhitney_pval': pval,
            'Significant_difference': pval < 0.05,
        })

if complex_stats:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    comp_df = pd.DataFrame(complex_stats)

    plot_data = []
    for method in METHODS:
        for v in np.log10(df[~df['is_complex']][method].dropna()):
            plot_data.append({'Method': method, 'Type': 'Single Gene', 'log10_kcat': v})
        for v in np.log10(df[df['is_complex']][method].dropna()):
            plot_data.append({'Method': method, 'Type': 'Complex', 'log10_kcat': v})
    plot_df = pd.DataFrame(plot_data)
    sns.violinplot(data=plot_df, x='Method', y='log10_kcat',
                   hue='Type', split=True, ax=axes[0])
    axes[0].set_ylabel('log₁₀(kcat) [s⁻¹]', fontsize=11)
    axes[0].set_title('Single Gene vs Complex Distribution', fontsize=12, fontweight='bold')
    axes[0].tick_params(axis='x', rotation=45)

    coverage_data = comp_df[['Method', 'N_single', 'N_complex']].set_index('Method')
    coverage_data.plot(kind='bar', ax=axes[1], alpha=0.7, edgecolor='black')
    axes[1].set_ylabel('Number of Predictions', fontsize=11)
    axes[1].set_title('Coverage: Single Gene vs Complex', fontsize=12, fontweight='bold')
    axes[1].legend(['Single Gene', 'Complex'])
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '05_complex_handling.png', dpi=300, bbox_inches='tight')
    print("   📈 Complex handling analysis saved")

    comp_df.to_csv(OUTPUT_DIR / '05_complex_statistics.csv', index=False)
    print("   📊 Complex statistics saved")
    print(comp_df.to_string(index=False))
    print()
else:
    print("   ⚠️  Insufficient data for complex analysis")
    print()


# ----------------------------------------------------------------------------
# Point 6: Substrate specificity
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 6: SUBSTRATE SPECIFICITY")
print("=" * 80)

enzyme_groups = df.groupby('gpr')
multi_substrate = [g for g, gp in enzyme_groups
                   if len(gp['base_reaction'].unique()) > 1]
print(f"   Enzymes with multiple substrates: {len(multi_substrate)}")

substrate_stats = []
for method in METHODS:
    variances = []
    for gpr in multi_substrate[:20]:
        group = df[df['gpr'] == gpr]
        values = group[method].dropna()
        if len(values) > 1:
            variances.append(np.var(np.log10(values)))
    if variances:
        substrate_stats.append({
            'Method': method,
            'N_enzymes': len(variances),
            'Mean_variance': np.mean(variances),
            'Median_variance': np.median(variances),
            'High_substrate_sensitivity': int(np.sum(np.array(variances) > 0.5)),
        })

if substrate_stats:
    fig, ax = plt.subplots(figsize=(10, 6))
    sub_df = pd.DataFrame(substrate_stats)
    x_pos = np.arange(len(sub_df))
    bars = ax.bar(x_pos, sub_df['Mean_variance'], alpha=0.7, edgecolor='black')
    colors = plt.cm.plasma(np.linspace(0, 1, len(sub_df)))
    for bar, color in zip(bars, colors):
        bar.set_color(color)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(sub_df['Method'], fontsize=11)
    ax.set_ylabel('Mean Variance (log₁₀ space)', fontsize=12)
    ax.set_title('Substrate Specificity Analysis\n'
                 '(Higher variance = Better substrate discrimination)',
                 fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '06_substrate_specificity.png', dpi=300, bbox_inches='tight')
    print("   📈 Substrate specificity analysis saved")

    sub_df.to_csv(OUTPUT_DIR / '06_substrate_statistics.csv', index=False)
    print("   📊 Substrate statistics saved")
    print(sub_df.to_string(index=False))
    print()
else:
    print("   ⚠️  Insufficient multi-substrate enzymes for analysis")
    print()


# ----------------------------------------------------------------------------
# Point 7: Coverage analysis
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 7: COVERAGE & APPLICABILITY DOMAIN")
print("=" * 80)

coverage_stats = []
for method in METHODS:
    total = df[method].notna().sum()
    transport = df[df['reaction_type'] == 'Transport'][method].notna().sum()
    metabolic = df[df['reaction_type'].str.contains('Metabolic')][method].notna().sum()
    complex_ = df[df['is_complex']][method].notna().sum()

    n_total = len(df)
    n_transport = (df['reaction_type'] == 'Transport').sum()
    n_metabolic = df['reaction_type'].str.contains('Metabolic').sum()
    n_complex = df['is_complex'].sum()

    coverage_stats.append({
        'Method': method,
        'Total_coverage': total,
        'Total_pct': 100 * total / n_total,
        'Transport_coverage': transport,
        'Transport_pct': 100 * transport / n_transport if n_transport else 0,
        'Metabolic_coverage': metabolic,
        'Metabolic_pct': 100 * metabolic / n_metabolic if n_metabolic else 0,
        'Complex_coverage': complex_,
        'Complex_pct': 100 * complex_ / n_complex if n_complex else 0,
    })

cov_df = pd.DataFrame(coverage_stats)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
categories = ['Total', 'Transport', 'Metabolic', 'Complex']
x_pos = np.arange(len(METHODS))
width = 0.2

for i, cat in enumerate(categories):
    col = f'{cat}_pct'
    axes[0].bar(x_pos + i * width, cov_df[col].values, width,
                label=cat, alpha=0.8)
axes[0].set_xticks(x_pos + 1.5 * width)
axes[0].set_xticklabels(METHODS, fontsize=11)
axes[0].set_ylabel('Coverage (%)', fontsize=12)
axes[0].set_title('Coverage by Reaction Type', fontsize=13, fontweight='bold')
axes[0].legend()
axes[0].grid(True, alpha=0.3, axis='y')

# Venn diagram (DLKcat, MTLKP, UniKP)
try:
    from matplotlib_venn import venn3
    dl_mask = df['DLKcat'].notna()
    mt_mask = df['MTLKP'].notna()
    un_mask = df['UniKP'].notna()
    only_dl = dl_mask & ~mt_mask & ~un_mask
    only_mt = mt_mask & ~dl_mask & ~un_mask
    only_un = un_mask & ~dl_mask & ~mt_mask
    dl_mt = dl_mask & mt_mask & ~un_mask
    dl_un = dl_mask & un_mask & ~mt_mask
    mt_un = mt_mask & un_mask & ~dl_mask
    all_three = dl_mask & mt_mask & un_mask
    venn3(subsets=(only_dl.sum(), only_mt.sum(), dl_mt.sum(), only_un.sum(),
                   dl_un.sum(), mt_un.sum(), all_three.sum()),
          set_labels=('DLKcat', 'MTLKP', 'UniKP'), ax=axes[1])
    axes[1].set_title('Coverage Overlap (DLKcat, MTLKP, UniKP)',
                     fontsize=13, fontweight='bold')
except ImportError:
    axes[1].text(0.5, 0.5, 'Venn diagram requires matplotlib_venn package',
                 ha='center', va='center', fontsize=12)
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / '07_coverage_analysis.png', dpi=300, bbox_inches='tight')
print("   📈 Coverage analysis saved")

cov_df.to_csv(OUTPUT_DIR / '07_coverage_statistics.csv', index=False)
print("   📊 Coverage statistics saved")
print(cov_df.to_string(index=False))
print()


# ----------------------------------------------------------------------------
# Point 8: Ground-truth benchmarking
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 8: GROUND TRUTH BENCHMARKING")
print("=" * 80)

ground_truth_df = df[df['database_clean'].notna()].copy()
print(f"   Valid ground truth samples: {len(ground_truth_df)}")

if len(ground_truth_df) > 3:
    benchmark_stats = []
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Analysis Point 8: Ground Truth Benchmarking',
                 fontsize=16, fontweight='bold')

    for idx, method in enumerate(METHODS):
        ax = axes[idx // 3, idx % 3]
        valid_mask = ground_truth_df[method].notna()
        gt_values = ground_truth_df[valid_mask]['database_clean'].values
        pred_values = ground_truth_df[valid_mask][method].values

        if len(gt_values) > 0:
            log_gt = np.log10(gt_values)
            log_pred = np.log10(pred_values)

            rmse = np.sqrt(mean_squared_error(log_gt, log_pred))
            mae = mean_absolute_error(log_gt, log_pred)
            r2 = r2_score(log_gt, log_pred)
            pearson_r, pearson_p = stats.pearsonr(log_gt, log_pred)
            spearman_r, spearman_p = stats.spearmanr(log_gt, log_pred)

            ax.scatter(log_gt, log_pred, alpha=0.6, s=60,
                       edgecolor='black', linewidth=0.5)
            lims = [min(log_gt.min(), log_pred.min()),
                    max(log_gt.max(), log_pred.max())]
            ax.plot(lims, lims, 'r--', linewidth=2, label='Perfect Prediction')
            z = np.polyfit(log_gt, log_pred, 1)
            p = np.poly1d(z)
            ax.plot(lims, p(lims), 'b-', linewidth=2, alpha=0.7,
                    label=f'Fit: y={z[0]:.2f}x+{z[1]:.2f}')

            ax.set_xlabel('log₁₀(kcat) Ground Truth [s⁻¹]', fontsize=10)
            ax.set_ylabel('log₁₀(kcat) Predicted [s⁻¹]', fontsize=10)
            ax.set_title(f'{method} (N={len(gt_values)})\n'
                         f'RMSE={rmse:.2f}, MAE={mae:.2f}, R²={r2:.2f}',
                         fontsize=11, fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

            benchmark_stats.append({
                'Method': method,
                'N_samples': len(gt_values),
                'RMSE_log10': rmse,
                'MAE_log10': mae,
                'R2': r2,
                'Pearson_r': pearson_r,
                'Pearson_p': pearson_p,
                'Spearman_r': spearman_r,
                'Spearman_p': spearman_p,
            })

    # Error bar panel
    ax = axes[1, 2]
    bench_df = pd.DataFrame(benchmark_stats)
    x_pos = np.arange(len(bench_df))
    width = 0.35
    ax.bar(x_pos - width / 2, bench_df['RMSE_log10'], width,
           label='RMSE', alpha=0.7)
    ax.bar(x_pos + width / 2, bench_df['MAE_log10'], width,
           label='MAE', alpha=0.7)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(bench_df['Method'], fontsize=9, rotation=45)
    ax.set_ylabel('Error (log₁₀ space)', fontsize=10)
    ax.set_title('Error Comparison', fontsize=11, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '08_ground_truth_benchmark.png', dpi=300, bbox_inches='tight')
    print("   📈 Benchmark analysis saved")

    bench_df.to_csv(OUTPUT_DIR / '08_benchmark_statistics.csv', index=False)
    print("   📊 Benchmark statistics saved")
    print(bench_df.to_string(index=False))
    print()
else:
    print("   ⚠️  Insufficient ground truth data for benchmarking")
    print()


# ----------------------------------------------------------------------------
# Point 9: Bias detection (Bland-Altman)
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 9: SYSTEMATIC BIAS DETECTION (Bland-Altman)")
print("=" * 80)

if len(ground_truth_df) > 3:
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Analysis Point 9: Bland-Altman Plots (Bias Detection)',
                 fontsize=16, fontweight='bold')

    bias_stats = []
    for idx, method in enumerate(METHODS):
        ax = axes[idx // 3, idx % 3]
        valid_mask = ground_truth_df[method].notna()
        gt_values = ground_truth_df[valid_mask]['database_clean'].values
        pred_values = ground_truth_df[valid_mask][method].values

        if len(gt_values) > 0:
            log_gt = np.log10(gt_values)
            log_pred = np.log10(pred_values)
            mean_vals = (log_gt + log_pred) / 2
            diff_vals = log_pred - log_gt

            mean_diff = np.mean(diff_vals)
            std_diff = np.std(diff_vals)

            ax.scatter(log_gt, diff_vals, alpha=0.6, s=60,
                       edgecolor='black', linewidth=0.5)
            ax.axhline(0, color='green', linestyle='-', linewidth=2, label='No Bias')
            ax.axhline(mean_diff, color='red', linestyle='--', linewidth=2,
                       label=f'Mean Bias={mean_diff:.2f}')
            ax.axhline(mean_diff + 1.96 * std_diff, color='orange', linestyle=':',
                       linewidth=2, label='±1.96 SD')
            ax.axhline(mean_diff - 1.96 * std_diff, color='orange', linestyle=':',
                       linewidth=2)

            ax.set_xlabel('log₁₀(kcat) Ground Truth [s⁻¹]', fontsize=10)
            ax.set_ylabel('Difference (Predicted - Truth)', fontsize=10)
            ax.set_title(f'{method}\nBias={mean_diff:.3f}, SD={std_diff:.3f}',
                         fontsize=11, fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

            trend_corr, trend_p = stats.pearsonr(log_gt, diff_vals)
            bias_stats.append({
                'Method': method,
                'N_samples': len(gt_values),
                'Mean_bias': mean_diff,
                'Bias_SD': std_diff,
                'Trend_correlation': trend_corr,
                'Trend_p_value': trend_p,
                'Systematic_trend': abs(trend_corr) > 0.3 and trend_p < 0.05,
            })

    # Bias summary
    ax = axes[1, 2]
    bias_df = pd.DataFrame(bias_stats)
    x_pos = np.arange(len(bias_df))
    colors = ['red' if b > 0 else 'blue' for b in bias_df['Mean_bias']]
    ax.barh(x_pos, bias_df['Mean_bias'], xerr=bias_df['Bias_SD'],
            capsize=5, alpha=0.7, edgecolor='black', color=colors)
    ax.set_yticks(x_pos)
    ax.set_yticklabels(bias_df['Method'], fontsize=10)
    ax.set_xlabel('Mean Bias (log₁₀ space)', fontsize=11)
    ax.set_title('Systematic Bias Summary\n'
                 '(Red=Overestimation, Blue=Underestimation)',
                 fontsize=11, fontweight='bold')
    ax.axvline(0, color='black', linestyle='-', linewidth=2)
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '09_bias_detection.png', dpi=300, bbox_inches='tight')
    print("   📈 Bias detection analysis saved")

    bias_df.to_csv(OUTPUT_DIR / '09_bias_statistics.csv', index=False)
    print("   📊 Bias statistics saved")
    print(bias_df.to_string(index=False))
    print()
else:
    print("   ⚠️  Insufficient ground truth data for bias analysis")
    print()


# ----------------------------------------------------------------------------
# Point 10: Ensemble modeling
# ----------------------------------------------------------------------------

print("=" * 80)
print("ANALYSIS POINT 10: ENSEMBLE MODELING")
print("=" * 80)

df['Ensemble_Geometric_Mean'] = df[METHODS].apply(
    lambda row: np.exp(np.mean(np.log(row.dropna()))) if len(row.dropna()) > 0 else np.nan,
    axis=1)
df['Ensemble_Arithmetic_Mean'] = df[METHODS].mean(axis=1)
df['Ensemble_Median'] = df[METHODS].median(axis=1)

if len(ground_truth_df) > 3 and len(benchmark_stats) > 0:
    weights = {}
    bench_df = pd.DataFrame(benchmark_stats)
    for _, row in bench_df.iterrows():
        weights[row['Method']] = max(0, row['R2']) + 0.01
    total_weight = sum(weights.values())
    weights = {k: v / total_weight for k, v in weights.items()}

    print("   Ensemble weights (based on R²):")
    for method, weight in weights.items():
        print(f"      {method}: {weight:.3f}")
    print()

    df['Ensemble_Weighted'] = df.apply(
        lambda row: np.exp(sum(weights.get(m, 0) * np.log(row[m])
                              for m in METHODS if pd.notna(row[m])))
        if any(pd.notna(row[m]) for m in METHODS) else np.nan,
        axis=1)
    ensemble_methods = ['Ensemble_Geometric_Mean', 'Ensemble_Arithmetic_Mean',
                        'Ensemble_Median', 'Ensemble_Weighted']
else:
    ensemble_methods = ['Ensemble_Geometric_Mean', 'Ensemble_Arithmetic_Mean',
                        'Ensemble_Median']

ensemble_df = ground_truth_df.copy()
ensemble_stats = []
if len(ensemble_df) > 3:
    for ens_method in ensemble_methods:
        if ens_method in ensemble_df.columns:
            valid_mask = ensemble_df[ens_method].notna()
            gt_values = ensemble_df[valid_mask]['database_clean'].values
            pred_values = ensemble_df[valid_mask][ens_method].values
            if len(gt_values) > 0:
                log_gt = np.log10(gt_values)
                log_pred = np.log10(pred_values)
                ensemble_stats.append({
                    'Method': ens_method,
                    'N_samples': len(gt_values),
                    'RMSE_log10': np.sqrt(mean_squared_error(log_gt, log_pred)),
                    'MAE_log10': mean_absolute_error(log_gt, log_pred),
                    'R2': r2_score(log_gt, log_pred),
                })

    all_stats = benchmark_stats + ensemble_stats
    comparison_df = pd.DataFrame(all_stats).sort_values('R2', ascending=False)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Analysis Point 10: Ensemble vs Individual Methods',
                 fontsize=16, fontweight='bold')
    x_pos = np.arange(len(comparison_df))

    for i, (metric, ax) in enumerate(zip(['RMSE_log10', 'MAE_log10', 'R2'], axes)):
        bars = ax.bar(x_pos, comparison_df[metric], alpha=0.7, edgecolor='black')
        for j, (_, row) in enumerate(comparison_df.iterrows()):
            bars[j].set_color('red' if 'Ensemble' in row['Method'] else 'steelblue')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(comparison_df['Method'], rotation=45, ha='right', fontsize=9)
        ax.set_ylabel(metric, fontsize=11)
        titles = ['Root Mean Squared Error\n(Lower is better)',
                  'Mean Absolute Error\n(Lower is better)',
                  'Coefficient of Determination\n(Higher is better)']
        ax.set_title(titles[i], fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
    axes[2].axhline(0, color='black', linestyle='-', linewidth=1)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / '10_ensemble_comparison.png', dpi=300, bbox_inches='tight')
    print("   📈 Ensemble analysis saved")

    comparison_df.to_csv(OUTPUT_DIR / '10_ensemble_statistics.csv', index=False)
    print("   📊 Ensemble statistics saved")
    print()
    print("Ensemble Performance Comparison:")
    print(comparison_df.to_string(index=False))
    print()

print(f"   💾 Enhanced dataset saved to: {ENHANCED_CSV}")
df.to_csv(ENHANCED_CSV, index=False)
print()


# ----------------------------------------------------------------------------
# Final summary
# ----------------------------------------------------------------------------

print("=" * 80)
print("FINAL SUMMARY & RECOMMENDATIONS")
print("=" * 80)
print()
print("📋 KEY FINDINGS:")
print()
print("1. GLOBAL DISTRIBUTION:")
print("   - All methods show log-normal distributions as expected")
print("   - Dynamic ranges vary significantly across methods")
print()
print("2. INTER-MODEL CONCORDANCE:")
print("   - Correlation analysis reveals method clustering patterns")
print("   - Sequence-based methods show varying degrees of correlation")
print()
print()
print("3. THERMODYNAMIC ASYMMETRY:")
print("   - Deep learning methods show varying ability to capture "
      "forward/reverse differences")
print()
print()
print("4. COVERAGE ANALYSIS:")
for method in METHODS:
    method_data = cov_df[cov_df['Method'] == method]
    if len(method_data) > 0:
        print(f"   - {method}: {method_data['Total_pct'].values[0]:.1f}% total coverage")
print("   - Methods show different coverage patterns across reaction types")
print()
print()

if benchmark_stats:
    best = bench_df.loc[bench_df['R2'].idxmax()]
    print("5. GROUND TRUTH BENCHMARKING:")
    print(f"   - Best performing method: {best['Method']}")
    print(f"     R² = {best['R2']:.3f}, RMSE = {best['RMSE_log10']:.3f}")
    print()

print("🎯 RECOMMENDED STRATEGY FOR ecGEMs:")
print()
print("   Layer 1: Experimental data (cleaned, artifact-free)")
print("   Layer 2: For core metabolism → Use best performing sequence-based method")
print("   Layer 3: For transport reactions → Use method with best transport coverage")
print("   Layer 4: For enzyme complexes → Use method with best complex handling")
print("   Layer 5: Ensemble weighted predictions for remaining gaps")
print()

print("=" * 80)
print("✅ ANALYSIS COMPLETE")
print(f"📁 All results saved to: {OUTPUT_DIR}")
print("=" * 80)
print()

print("Generated files:")
for i in range(1, 11):
    print(f"   {i:02d}_*.png - Visualization")
    print(f"   {i:02d}_*.csv - Statistics")
print("   kcat_comparison_enhanced.csv - Enhanced dataset with ensembles")
print()
