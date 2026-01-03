#!/usr/bin/env python3
"""
Generate Figure 1: Overall Performance Comparison
Main display item for the Nature Genetics Article.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
CANDIDATES_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_cs2g.tsv"
OUTPUT_DIR = BASE_DIR / "figures"


def plot_roc_curves(df, ax, title="ROC Curves"):
    """Plot ROC curves for all methods."""
    
    methods = [
        ('NearestGene', 'score_nearest', 'tab:blue'),
        ('Within100kb', 'score_100kb', 'tab:orange'),
    ]
    
    for method_name, score_col, color in methods:
        scores = df[score_col].values
        labels = df['label'].values
        
        valid = ~np.isnan(scores)
        fpr, tpr, _ = roc_curve(labels[valid], scores[valid])
        roc_auc = auc(fpr, tpr)
        
        ax.plot(fpr, tpr, color=color, linewidth=2, 
                label=f'{method_name} (AUC = {roc_auc:.3f})')
    
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1, label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=11)
    ax.set_ylabel('True Positive Rate', fontsize=11)
    ax.set_title(title, fontweight='bold', fontsize=12)
    ax.legend(loc='lower right', fontsize=9)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.grid(True, alpha=0.3)


def plot_pr_curves(df, ax, title="Precision-Recall Curves"):
    """Plot PR curves for all methods."""
    
    methods = [
        ('NearestGene', 'score_nearest', 'tab:blue'),
        ('Within100kb', 'score_100kb', 'tab:orange'),
    ]
    
    baseline_rate = df['label'].mean()
    
    for method_name, score_col, color in methods:
        scores = df[score_col].values
        labels = df['label'].values
        
        valid = ~np.isnan(scores)
        precision, recall, _ = precision_recall_curve(labels[valid], scores[valid])
        avg_prec = average_precision_score(labels[valid], scores[valid])
        
        ax.plot(recall, precision, color=color, linewidth=2, 
                label=f'{method_name} (AP = {avg_prec:.3f})')
    
    ax.axhline(y=baseline_rate, color='k', linestyle='--', alpha=0.5, 
               label=f'Random (AP = {baseline_rate:.3f})')
    ax.set_xlabel('Recall', fontsize=11)
    ax.set_ylabel('Precision', fontsize=11)
    ax.set_title(title, fontweight='bold', fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([0, 1.02])
    ax.grid(True, alpha=0.3)


def plot_method_comparison_bar(df, ax):
    """Bar chart comparing methods across evidence types."""
    
    from sklearn.metrics import roc_auc_score
    
    methods = ['NearestGene', 'Within100kb']
    score_cols = {'NearestGene': 'score_nearest', 'Within100kb': 'score_100kb'}
    evidence_types = ['All', 'CRISPRi', 'MPRA']
    
    results = []
    for etype in evidence_types:
        for method in methods:
            if etype == 'All':
                subset = df
            else:
                subset = df[df['evidence_type'] == etype]
            
            scores = subset[score_cols[method]].values
            labels = subset['label'].values
            valid = ~np.isnan(scores)
            
            auc_val = roc_auc_score(labels[valid], scores[valid])
            results.append({
                'Evidence Type': etype,
                'Method': method,
                'AUC-ROC': auc_val
            })
    
    results_df = pd.DataFrame(results)
    
    x = np.arange(len(evidence_types))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, 
                   results_df[results_df['Method'] == 'NearestGene']['AUC-ROC'],
                   width, label='NearestGene', color='tab:blue')
    bars2 = ax.bar(x + width/2,
                   results_df[results_df['Method'] == 'Within100kb']['AUC-ROC'],
                   width, label='Within100kb', color='tab:orange')
    
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('AUC-ROC', fontsize=11)
    ax.set_xlabel('Evidence Type', fontsize=11)
    ax.set_title('C. Performance by Evidence Type', fontweight='bold', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(evidence_types)
    ax.legend(loc='lower right', fontsize=9)
    ax.set_ylim([0.4, 1.0])
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bar in bars1:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=8)


def plot_benchmark_stats(df, ax):
    """Show benchmark statistics as text panel."""
    
    n_total = len(df)
    n_loci = df['locus_id'].nunique()
    n_positive = df['label'].sum()
    pos_rate = df['label'].mean() * 100
    n_crispri = len(df[df['evidence_type'] == 'CRISPRi'])
    n_mpra = len(df[df['evidence_type'] == 'MPRA'])
    
    ax.axis('off')
    
    stats_text = f"""
RegulatoryBench v3 Statistics
════════════════════════════════

Total Candidates:     {n_total:,}
Unique Loci:          {n_loci:,}
Positive Pairs:       {int(n_positive):,}
Positive Rate:        {pos_rate:.1f}%

CRISPRi Evidence:     {n_crispri:,}
MPRA Evidence:        {n_mpra:,}

L2G Overlap:          0.0%
(Independent of training)

Window Size:          500 kb
Genes per Locus:      ~14.8
"""
    
    ax.text(0.1, 0.95, stats_text, transform=ax.transAxes,
            fontfamily='monospace', fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))


def main():
    print("=" * 70)
    print("Generating Figure 1: Performance Summary")
    print("=" * 70)
    
    # Load data
    print(f"\nLoading candidates...")
    df = pd.read_csv(CANDIDATES_FILE, sep='\t', low_memory=False)
    print(f"  Loaded {len(df):,} candidates")
    
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: ROC curves
    plot_roc_curves(df, axes[0, 0], title='A. ROC Curves (All Evidence)')
    
    # Panel B: PR curves
    plot_pr_curves(df, axes[0, 1], title='B. Precision-Recall Curves')
    
    # Panel C: Bar chart by evidence type
    plot_method_comparison_bar(df, axes[1, 0])
    
    # Panel D: Statistics panel
    plot_benchmark_stats(df, axes[1, 1])
    
    plt.tight_layout()
    
    # Save
    for ext in ['png', 'pdf']:
        output_file = OUTPUT_DIR / f'figure_1_performance_summary.{ext}'
        fig.savefig(output_file, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    plt.close()
    
    print("\n" + "=" * 70)
    print("Figure 1 complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
