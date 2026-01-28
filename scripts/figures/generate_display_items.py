#!/usr/bin/env python3
"""
Generate display items (figures and tables) for Nature Genetics manuscript.

Target: ≤8 display items total
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import json

# Set publication-quality defaults
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

BASE_DIR = Path(__file__).parent.parent
OUTPUT_DIR = BASE_DIR / "figures"
DATA_DIR = BASE_DIR / "data" / "processed" / "baselines"


def load_benchmark_data():
    """Load benchmark and evaluation data."""
    benchmark = pd.read_csv(DATA_DIR / "regulatorybench_v3.tsv", sep="\t")
    candidates = pd.read_csv(DATA_DIR / "evaluation_candidates.tsv", sep="\t")
    loci = pd.read_csv(DATA_DIR / "evaluation_loci.tsv", sep="\t")
    
    with open(DATA_DIR / "regulatorybench_v3_summary.json") as f:
        summary = json.load(f)
    
    return benchmark, candidates, loci, summary


def figure_1_benchmark_overview(benchmark, summary):
    """
    Figure 1: RegulatoryBench v3 Overview
    Panel A: Data sources and composition
    Panel B: Locus-level independence from L2G training
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    # Panel A: Benchmark composition
    ax = axes[0]
    sources = ['CRISPRi\n(ENCODE)', 'MPRA\n(Abell 2022)']
    counts = [summary['evidence_types']['CRISPRi'], summary['evidence_types']['MPRA']]
    colors = ['#2E86AB', '#A23B72']
    
    bars = ax.bar(sources, counts, color=colors, edgecolor='black', linewidth=1)
    ax.set_ylabel('Number of Validated Pairs', fontsize=11)
    ax.set_title('A) Benchmark Composition', fontsize=12, fontweight='bold', loc='left')
    
    # Add value labels
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 100,
                f'{count:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylim(0, max(counts) * 1.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel B: Independence from L2G
    ax = axes[1]
    
    # Create stacked bar showing overlap
    categories = ['Locus\nOverlap', 'Gene\nOverlap']
    overlap = [0.0, 3.2]  # From validation
    independent = [100.0, 96.8]
    
    x = np.arange(len(categories))
    width = 0.6
    
    bars1 = ax.bar(x, independent, width, label='Independent', color='#28A745')
    bars2 = ax.bar(x, overlap, width, bottom=independent, label='Overlap with L2G Training', color='#DC3545')
    
    ax.set_ylabel('Percentage (%)', fontsize=11)
    ax.set_title('B) Independence from L2G Training', fontsize=12, fontweight='bold', loc='left')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylim(0, 105)
    ax.legend(loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add percentage labels
    for i, (ind, ovl) in enumerate(zip(independent, overlap)):
        ax.text(i, ind/2, f'{ind:.1f}%', ha='center', va='center', fontsize=10, fontweight='bold', color='white')
        if ovl > 0:
            ax.text(i, ind + ovl/2, f'{ovl:.1f}%', ha='center', va='center', fontsize=9, color='white')
    
    plt.tight_layout()
    return fig


def figure_2_baseline_comparison(candidates, loci):
    """
    Figure 2: Baseline Method Comparison
    Panel A: ROC curves
    Panel B: Stratified by evidence type
    """
    from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    
    # Panel A: Overall ROC curves
    ax = axes[0]
    
    y_true = candidates['label'].values
    
    methods = {
        'NearestGene': candidates['score_nearest'].values,
        'Within100kb': candidates['score_100kb'].values
    }
    
    colors = {'NearestGene': '#2E86AB', 'Within100kb': '#F18F01'}
    
    for name, scores in methods.items():
        fpr, tpr, _ = roc_curve(y_true, scores)
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, label=f'{name} (AUC={roc_auc:.3f})', color=colors[name], linewidth=2)
    
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('A) Overall ROC Curves', fontsize=12, fontweight='bold', loc='left')
    ax.legend(loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel B: Stratified bar chart
    ax = axes[1]
    
    # Calculate stratified AUCs
    crispri_mask = candidates['evidence_type'] == 'CRISPRi'
    mpra_mask = candidates['evidence_type'] == 'MPRA'
    
    results = []
    for name, scores in methods.items():
        # CRISPRi
        fpr, tpr, _ = roc_curve(y_true[crispri_mask], scores[crispri_mask])
        crispri_auc = auc(fpr, tpr)
        
        # MPRA
        fpr, tpr, _ = roc_curve(y_true[mpra_mask], scores[mpra_mask])
        mpra_auc = auc(fpr, tpr)
        
        results.append({
            'method': name,
            'CRISPRi': crispri_auc,
            'MPRA': mpra_auc
        })
    
    x = np.arange(len(results))
    width = 0.35
    
    crispri_aucs = [r['CRISPRi'] for r in results]
    mpra_aucs = [r['MPRA'] for r in results]
    
    bars1 = ax.bar(x - width/2, crispri_aucs, width, label='CRISPRi', color='#2E86AB')
    bars2 = ax.bar(x + width/2, mpra_aucs, width, label='MPRA', color='#A23B72')
    
    ax.set_ylabel('AUC-ROC')
    ax.set_title('B) Stratified by Evidence Type', fontsize=12, fontweight='bold', loc='left')
    ax.set_xticks(x)
    ax.set_xticklabels([r['method'] for r in results])
    ax.set_ylim(0, 1.1)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.legend(loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.02,
                    f'{height:.2f}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig


def figure_3_distance_distribution(candidates):
    """
    Figure 3: Distance to TSS analysis
    Shows why nearest gene works for CRISPRi but not MPRA
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    # Panel A: Distance distributions by evidence type
    ax = axes[0]
    
    crispri = candidates[candidates['evidence_type'] == 'CRISPRi']
    mpra = candidates[candidates['evidence_type'] == 'MPRA']
    
    crispri_pos = crispri[crispri['label'] == 1]['distance'] / 1000
    mpra_pos = mpra[mpra['label'] == 1]['distance'] / 1000
    
    bins = np.linspace(0, 500, 50)
    
    ax.hist(crispri_pos, bins=bins, alpha=0.7, label=f'CRISPRi (n={len(crispri_pos)})', 
            color='#2E86AB', density=True)
    ax.hist(mpra_pos, bins=bins, alpha=0.7, label=f'MPRA (n={len(mpra_pos)})', 
            color='#A23B72', density=True)
    
    ax.set_xlabel('Distance to TSS (kb)')
    ax.set_ylabel('Density')
    ax.set_title('A) Distance Distribution of True Positives', fontsize=12, fontweight='bold', loc='left')
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel B: Fraction nearest gene by evidence type
    ax = axes[1]
    
    # Calculate fraction that are nearest gene
    crispri_loci = candidates[candidates['evidence_type'] == 'CRISPRi'].groupby('locus_id')
    mpra_loci = candidates[candidates['evidence_type'] == 'MPRA'].groupby('locus_id')
    
    def fraction_nearest_is_positive(group):
        pos_genes = set(group[group['label'] == 1]['gene_symbol'])
        nearest = group.loc[group['distance'].idxmin(), 'gene_symbol']
        return nearest in pos_genes
    
    crispri_nearest = crispri_loci.apply(fraction_nearest_is_positive).mean()
    mpra_nearest = mpra_loci.apply(fraction_nearest_is_positive).mean()
    
    categories = ['CRISPRi', 'MPRA']
    fractions = [crispri_nearest * 100, mpra_nearest * 100]
    colors = ['#2E86AB', '#A23B72']
    
    bars = ax.bar(categories, fractions, color=colors, edgecolor='black', linewidth=1)
    ax.set_ylabel('% Where Nearest Gene is Correct')
    ax.set_title('B) Nearest Gene Accuracy', fontsize=12, fontweight='bold', loc='left')
    ax.set_ylim(0, 100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for bar, frac in zip(bars, fractions):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{frac:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    return fig


def table_1_benchmark_statistics(summary, candidates, loci):
    """Generate summary statistics table."""
    
    stats = {
        'Metric': [
            'Total loci',
            'CRISPRi loci',
            'MPRA loci',
            'Total candidate genes',
            'Positive labels',
            'Negative labels',
            'Positive rate',
            'Mean candidates per locus',
            'Unique genes',
            'L2G training locus overlap',
            'L2G training gene overlap'
        ],
        'Value': [
            f"{len(loci):,}",
            f"{len(loci[loci['evidence_type'] == 'CRISPRi']):,}",
            f"{len(loci[loci['evidence_type'] == 'MPRA']):,}",
            f"{len(candidates):,}",
            f"{candidates['label'].sum():,}",
            f"{(~candidates['label'].astype(bool)).sum():,}",
            f"{candidates['label'].mean()*100:.1f}%",
            f"{len(candidates) / len(loci):.1f}",
            f"{summary['unique_genes']:,}",
            "0.0%",
            "3.2%"
        ]
    }
    
    return pd.DataFrame(stats)


def table_2_baseline_results(candidates):
    """Generate baseline results table."""
    from sklearn.metrics import roc_auc_score, average_precision_score
    
    y_true = candidates['label'].values
    crispri_mask = candidates['evidence_type'] == 'CRISPRi'
    mpra_mask = candidates['evidence_type'] == 'MPRA'
    
    methods = {
        'NearestGene': candidates['score_nearest'].values,
        'Within100kb': candidates['score_100kb'].values
    }
    
    results = []
    for name, scores in methods.items():
        results.append({
            'Method': name,
            'Overall AUC-ROC': f"{roc_auc_score(y_true, scores):.3f}",
            'Overall AUC-PR': f"{average_precision_score(y_true, scores):.3f}",
            'CRISPRi AUC-ROC': f"{roc_auc_score(y_true[crispri_mask], scores[crispri_mask]):.3f}",
            'MPRA AUC-ROC': f"{roc_auc_score(y_true[mpra_mask], scores[mpra_mask]):.3f}"
        })
    
    return pd.DataFrame(results)


def main():
    """Generate all display items."""
    print("Loading data...")
    benchmark, candidates, loci, summary = load_benchmark_data()
    
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Generate figures
    print("\nGenerating Figure 1: Benchmark Overview...")
    fig1 = figure_1_benchmark_overview(benchmark, summary)
    fig1.savefig(OUTPUT_DIR / "figure_1_benchmark_overview.png")
    fig1.savefig(OUTPUT_DIR / "figure_1_benchmark_overview.pdf")
    print(f"  Saved to {OUTPUT_DIR / 'figure_1_benchmark_overview.png'}")
    
    print("\nGenerating Figure 2: Baseline Comparison...")
    fig2 = figure_2_baseline_comparison(candidates, loci)
    fig2.savefig(OUTPUT_DIR / "figure_2_baseline_comparison.png")
    fig2.savefig(OUTPUT_DIR / "figure_2_baseline_comparison.pdf")
    print(f"  Saved to {OUTPUT_DIR / 'figure_2_baseline_comparison.png'}")
    
    print("\nGenerating Figure 3: Distance Analysis...")
    fig3 = figure_3_distance_distribution(candidates)
    fig3.savefig(OUTPUT_DIR / "figure_3_distance_analysis.png")
    fig3.savefig(OUTPUT_DIR / "figure_3_distance_analysis.pdf")
    print(f"  Saved to {OUTPUT_DIR / 'figure_3_distance_analysis.png'}")
    
    # Generate tables
    print("\nGenerating Table 1: Benchmark Statistics...")
    table1 = table_1_benchmark_statistics(summary, candidates, loci)
    table1.to_csv(OUTPUT_DIR / "table_1_benchmark_statistics.csv", index=False)
    print(table1.to_string(index=False))
    
    print("\nGenerating Table 2: Baseline Results...")
    table2 = table_2_baseline_results(candidates)
    table2.to_csv(OUTPUT_DIR / "table_2_baseline_results.csv", index=False)
    print(table2.to_string(index=False))
    
    print("\n" + "="*60)
    print("DISPLAY ITEMS GENERATED")
    print("="*60)
    print(f"\nFigures saved to: {OUTPUT_DIR}")
    print("\nFor Nature Genetics (≤8 display items):")
    print("  Figure 1: Benchmark overview + independence")
    print("  Figure 2: Baseline ROC + stratified performance")
    print("  Figure 3: Distance distribution analysis")
    print("  Table 1: Benchmark statistics")
    print("  Table 2: Baseline results")
    print("\nRemaining slots: 3 (for mechanistic model results)")
    
    plt.close('all')


if __name__ == "__main__":
    main()
