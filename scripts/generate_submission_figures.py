#!/usr/bin/env python3
"""
generate_submission_figures.py
==============================
Generates all submission-ready figures from master_results.parquet.
All bars include n values and 95% Wilson CI error bars.

Figures:
- Fig 1A: Method comparison bar plot (Top-1 accuracy by method)
- Fig 1B: Mechanism stratification (CODING vs REGULATORY)
- Fig 2: Coverage vs Accuracy scatter
- Fig 3: Score distributions by mechanism
- Extended Data figures as needed

Author: Mechanism-GWAS-Causal-Graphs
Date: 2025-01
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats as scipy_stats
import json
from datetime import datetime

PROJECT_ROOT = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures" / "submission"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Nature Genetics style
plt.style.use('seaborn-v0_8-whitegrid')
NATURE_COLORS = {
    'L2G': '#2171B5',       # Blue
    'cS2G': '#6A51A3',      # Purple
    'NearestGene': '#969696', # Gray
    'FLAMES': '#E6550D',    # Orange (for future)
    'Calibrated': '#31A354', # Green (for future)
    'CODING': '#1B7837',
    'REGULATORY': '#762A83'
}

def wilson_ci(successes: int, n: int, alpha: float = 0.05):
    """Compute Wilson score interval for binomial proportion."""
    if n == 0:
        return 0.0, 0.0, 0.0
    p = successes / n
    z = scipy_stats.norm.ppf(1 - alpha / 2)
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    spread = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denominator
    return p, max(0, center - spread), min(1, center + spread)

def load_master_results():
    """Load the master results table."""
    master_path = RESULTS_DIR / "master_results.parquet"
    df = pd.read_parquet(master_path)
    print(f"Loaded master results: {len(df)} rows")
    return df

def load_statistics():
    """Load pre-computed statistics."""
    stats_path = RESULTS_DIR / "method_statistics.json"
    with open(stats_path) as f:
        return json.load(f)

def fig1a_method_comparison(df, stats):
    """
    Figure 1A: Method comparison bar plot.
    Shows Top-1 accuracy for each method with n values and 95% CI.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    methods = ['L2G', 'cS2G', 'NearestGene']
    available_methods = [m for m in methods if m in stats]
    
    x_pos = range(len(available_methods))
    bars = []
    
    for i, method in enumerate(available_methods):
        method_stats = stats[method]
        accuracy = method_stats['top1_accuracy']
        ci_lo = method_stats['top1_accuracy_ci_lo']
        ci_hi = method_stats['top1_accuracy_ci_hi']
        n = method_stats['n_with_prediction']
        n_correct = method_stats['top1_n']
        
        # Bar with error bar
        bar = ax.bar(i, accuracy, 
                     color=NATURE_COLORS.get(method, '#666666'),
                     edgecolor='black', linewidth=1.5,
                     yerr=[[accuracy - ci_lo], [ci_hi - accuracy]],
                     capsize=5, error_kw={'linewidth': 2})
        bars.append(bar)
        
        # Add n annotation above bar
        ax.annotate(f'n={n}', (i, accuracy + (ci_hi - accuracy) + 3),
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylabel('Top-1 Accuracy (%)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Method', fontsize=12, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(available_methods, fontsize=11)
    ax.set_ylim(0, 110)
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5, label='Random baseline')
    
    ax.set_title('Gene Prioritization Accuracy by Method', fontsize=14, fontweight='bold')
    
    # Add legend for CI
    ax.text(0.02, 0.98, 'Error bars: 95% Wilson CI', transform=ax.transAxes,
            fontsize=9, verticalalignment='top', style='italic')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'Fig1A_method_comparison.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig1A_method_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved Fig1A_method_comparison.pdf/png")

def fig1b_mechanism_stratification(df, stats):
    """
    Figure 1B: CODING vs REGULATORY stratification for L2G.
    Shows performance gap between mechanism classes.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Focus on L2G
    if 'L2G' not in stats:
        print("L2G not in stats - skipping Fig1B")
        return
    
    l2g_stats = stats['L2G']
    
    # Prepare data for grouped bar
    mechanisms = ['CODING', 'REGULATORY']
    x_pos = range(len(mechanisms))
    
    for i, mech in enumerate(mechanisms):
        key = mech.lower()
        acc = l2g_stats[f'{key}_top1']
        ci_lo = l2g_stats[f'{key}_top1_ci_lo']
        ci_hi = l2g_stats[f'{key}_top1_ci_hi']
        n = l2g_stats[f'{key}_n_pred']
        
        bar = ax.bar(i, acc,
                     color=NATURE_COLORS[mech],
                     edgecolor='black', linewidth=1.5,
                     yerr=[[acc - ci_lo], [ci_hi - acc]],
                     capsize=5, error_kw={'linewidth': 2})
        
        ax.annotate(f'n={n}', (i, acc + (ci_hi - acc) + 3),
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylabel('Top-1 Accuracy (%)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Mechanism Class', fontsize=12, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(mechanisms, fontsize=11)
    ax.set_ylim(0, 110)
    
    ax.set_title('L2G Performance by Mechanism Class', fontsize=14, fontweight='bold')
    
    # Add gap annotation
    coding_acc = l2g_stats['coding_top1']
    reg_acc = l2g_stats['regulatory_top1']
    gap = coding_acc - reg_acc
    
    ax.annotate(f'Gap: {gap:.1f} pp', xy=(0.5, min(coding_acc, reg_acc) - 5),
               fontsize=10, ha='center', fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.text(0.02, 0.98, 'Error bars: 95% Wilson CI', transform=ax.transAxes,
            fontsize=9, verticalalignment='top', style='italic')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'Fig1B_mechanism_stratification.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig1B_mechanism_stratification.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved Fig1B_mechanism_stratification.pdf/png")

def fig2_coverage_vs_accuracy(df, stats):
    """
    Figure 2: Coverage vs Accuracy scatter plot.
    Each point is a method, showing trade-off.
    """
    fig, ax = plt.subplots(figsize=(8, 7))
    
    for method, method_stats in stats.items():
        coverage = method_stats['coverage']
        accuracy = method_stats['top1_accuracy']
        n = method_stats['n_with_prediction']
        
        # Point size proportional to n
        size = max(100, n / 2)
        
        ax.scatter(coverage, accuracy,
                   s=size, c=NATURE_COLORS.get(method, '#666666'),
                   edgecolor='black', linewidth=2, alpha=0.8,
                   label=f"{method} (n={n})")
        
        # Add method label
        ax.annotate(method, (coverage + 1.5, accuracy + 1.5),
                   fontsize=10, fontweight='bold')
    
    ax.set_xlabel('Coverage (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Top-1 Accuracy (%)', fontsize=12, fontweight='bold')
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    
    # Add quadrant lines
    ax.axvline(x=80, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=80, color='gray', linestyle='--', alpha=0.5)
    
    ax.set_title('Coverage-Accuracy Trade-off', fontsize=14, fontweight='bold')
    ax.legend(loc='lower left', fontsize=9)
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'Fig2_coverage_vs_accuracy.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig2_coverage_vs_accuracy.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved Fig2_coverage_vs_accuracy.pdf/png")

def fig3_score_distributions(df):
    """
    Figure 3: Score distributions by mechanism class for L2G.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    l2g_df = df[(df['method'] == 'L2G') & (df['has_prediction'] == True)]
    
    # Left: CODING vs REGULATORY distributions
    ax = axes[0]
    for mech in ['CODING', 'REGULATORY']:
        mech_data = l2g_df[l2g_df['mechanism_class'] == mech]['score'].dropna()
        if len(mech_data) > 0:
            ax.hist(mech_data, bins=20, alpha=0.6, 
                   label=f'{mech} (n={len(mech_data)})',
                   color=NATURE_COLORS[mech], edgecolor='black')
    
    ax.set_xlabel('L2G Score', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('L2G Score Distributions', fontsize=12, fontweight='bold')
    ax.legend()
    
    # Right: Correct vs Incorrect predictions
    ax = axes[1]
    correct = l2g_df[l2g_df['top1_correct'] == True]['score'].dropna()
    incorrect = l2g_df[l2g_df['top1_correct'] == False]['score'].dropna()
    
    ax.hist(correct, bins=20, alpha=0.6, label=f'Correct (n={len(correct)})',
           color='#31A354', edgecolor='black')
    if len(incorrect) > 0:
        ax.hist(incorrect, bins=20, alpha=0.6, label=f'Incorrect (n={len(incorrect)})',
               color='#E6550D', edgecolor='black')
    
    ax.set_xlabel('L2G Score', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('L2G Scores: Correct vs Incorrect', fontsize=12, fontweight='bold')
    ax.legend()
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'Fig3_score_distributions.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig3_score_distributions.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved Fig3_score_distributions.pdf/png")

def fig4_nearest_gene_comparison(df, stats):
    """
    Figure 4: NearestGene baseline vs L2G by mechanism.
    Demonstrates value of machine learning over distance baseline.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    methods = ['L2G', 'NearestGene']
    mechanisms = ['CODING', 'REGULATORY']
    
    x = np.arange(len(mechanisms))
    width = 0.35
    
    for i, method in enumerate(methods):
        if method not in stats:
            continue
        method_stats = stats[method]
        
        accs = []
        ci_los = []
        ci_his = []
        ns = []
        
        for mech in mechanisms:
            key = mech.lower()
            acc = method_stats[f'{key}_top1']
            ci_lo = method_stats[f'{key}_top1_ci_lo']
            ci_hi = method_stats[f'{key}_top1_ci_hi']
            n = method_stats[f'{key}_n_pred']
            
            accs.append(acc)
            ci_los.append(acc - ci_lo)
            ci_his.append(ci_hi - acc)
            ns.append(n)
        
        offset = (i - 0.5) * width
        bars = ax.bar(x + offset, accs,
                      width=width - 0.02,
                      label=method,
                      color=NATURE_COLORS.get(method, '#666666'),
                      edgecolor='black', linewidth=1.5,
                      yerr=[ci_los, ci_his],
                      capsize=4, error_kw={'linewidth': 1.5})
        
        # Add n values
        for j, (bar, n) in enumerate(zip(bars, ns)):
            height = bar.get_height()
            ax.annotate(f'n={n}',
                       xy=(bar.get_x() + bar.get_width() / 2, height + ci_his[j] + 2),
                       ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    ax.set_ylabel('Top-1 Accuracy (%)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Mechanism Class', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(mechanisms, fontsize=11)
    ax.set_ylim(0, 115)
    ax.legend(loc='upper right')
    
    ax.set_title('L2G vs Nearest Gene Baseline', fontsize=14, fontweight='bold')
    
    ax.text(0.02, 0.98, 'Error bars: 95% Wilson CI', transform=ax.transAxes,
            fontsize=9, verticalalignment='top', style='italic')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'Fig4_baseline_comparison.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig4_baseline_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Saved Fig4_baseline_comparison.pdf/png")

def generate_summary_table(stats):
    """Generate summary table for manuscript."""
    rows = []
    for method, method_stats in stats.items():
        row = {
            'Method': method,
            'Coverage': f"{method_stats['coverage']:.1f}%",
            'n': method_stats['n_with_prediction'],
            'Top-1 Accuracy': f"{method_stats['top1_accuracy']:.1f}%",
            '95% CI': f"[{method_stats['top1_accuracy_ci_lo']:.1f}%-{method_stats['top1_accuracy_ci_hi']:.1f}%]",
            'CODING': f"{method_stats['coding_top1']:.1f}% (n={method_stats['coding_n_pred']})",
            'REGULATORY': f"{method_stats['regulatory_top1']:.1f}% (n={method_stats['regulatory_n_pred']})"
        }
        rows.append(row)
    
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(FIGURES_DIR / 'Table1_method_summary.tsv', sep='\t', index=False)
    print("Saved Table1_method_summary.tsv")
    return summary_df

def main():
    """Generate all figures."""
    print("=" * 70)
    print("GENERATING SUBMISSION FIGURES")
    print("=" * 70)
    print(f"Output directory: {FIGURES_DIR}")
    print(f"Started at: {datetime.now().isoformat()}")
    
    # Load data
    df = load_master_results()
    stats = load_statistics()
    
    print("\nGenerating figures...")
    print("-" * 50)
    
    # Generate all figures
    fig1a_method_comparison(df, stats)
    fig1b_mechanism_stratification(df, stats)
    fig2_coverage_vs_accuracy(df, stats)
    fig3_score_distributions(df)
    fig4_nearest_gene_comparison(df, stats)
    
    # Generate summary table
    print("\nGenerating summary table...")
    summary = generate_summary_table(stats)
    print(summary.to_string(index=False))
    
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"All figures saved to: {FIGURES_DIR}")

if __name__ == "__main__":
    main()
