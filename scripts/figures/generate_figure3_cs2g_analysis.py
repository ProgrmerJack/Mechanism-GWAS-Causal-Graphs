#!/usr/bin/env python3
"""
Generate Figure 3: cS2G Analysis
Explains why cS2G shows AUC ≈ 0.50 on RegulatoryBench.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
CANDIDATES_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_cs2g.tsv"
OUTPUT_DIR = BASE_DIR / "figures"


def main():
    print("=" * 70)
    print("Generating Figure 3: cS2G Analysis")
    print("=" * 70)
    
    # Load data
    print(f"\nLoading candidates...")
    df = pd.read_csv(CANDIDATES_FILE, sep='\t', low_memory=False)
    print(f"  Loaded {len(df):,} candidates")
    
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Filter to candidates with cS2G scores
    df_cs2g = df[~df['score_cs2g'].isna()].copy()
    print(f"  With cS2G scores: {len(df_cs2g):,}")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: cS2G score distribution by label
    ax1 = axes[0, 0]
    
    # Separate positive and negative
    pos_scores = df_cs2g[df_cs2g['label'] == 1]['score_cs2g']
    neg_scores = df_cs2g[df_cs2g['label'] == 0]['score_cs2g']
    
    bins = np.linspace(0, 1, 50)
    ax1.hist(neg_scores, bins=bins, alpha=0.6, label=f'Negative (n={len(neg_scores):,})', 
             color='tab:blue', density=True)
    ax1.hist(pos_scores, bins=bins, alpha=0.6, label=f'Positive (n={len(pos_scores):,})', 
             color='tab:red', density=True)
    ax1.set_xlabel('cS2G Score')
    ax1.set_ylabel('Density')
    ax1.set_title('A. cS2G Score Distribution', fontweight='bold', fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Per-locus score analysis
    ax2 = axes[0, 1]
    
    # For each locus, count unique cS2G scores among candidates
    locus_stats = df_cs2g.groupby('locus_id').agg({
        'score_cs2g': lambda x: x.nunique(),
        'gene_symbol': 'count',
        'label': 'sum'
    }).rename(columns={
        'score_cs2g': 'unique_scores',
        'gene_symbol': 'n_candidates',
        'label': 'n_positive'
    })
    
    # Histogram of unique scores per locus
    ax2.hist(locus_stats['unique_scores'], bins=range(1, 20), 
             edgecolor='black', alpha=0.7, color='tab:green')
    ax2.axvline(x=locus_stats['unique_scores'].median(), color='red', 
                linestyle='--', label=f'Median = {locus_stats["unique_scores"].median():.0f}')
    ax2.set_xlabel('Unique cS2G Scores per Locus')
    ax2.set_ylabel('Number of Loci')
    ax2.set_title('B. Score Diversity Within Loci', fontweight='bold', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Panel C: Example locus showing identical scores
    ax3 = axes[1, 0]
    
    # Find a locus with multiple candidates that all have the same score
    single_score_loci = locus_stats[locus_stats['unique_scores'] == 1]
    if len(single_score_loci) > 0:
        example_locus = single_score_loci.index[0]
        example_df = df_cs2g[df_cs2g['locus_id'] == example_locus].copy()
        example_df = example_df.sort_values('distance').head(10)
        
        colors = ['red' if l == 1 else 'blue' for l in example_df['label']]
        bars = ax3.barh(range(len(example_df)), example_df['score_cs2g'], 
                       color=colors, alpha=0.7)
        ax3.set_yticks(range(len(example_df)))
        ax3.set_yticklabels([f"{g[:15]}..." if len(g) > 15 else g 
                            for g in example_df['gene_symbol']], fontsize=8)
        ax3.set_xlabel('cS2G Score')
        ax3.set_title(f'C. Example Locus: All Genes Get Same Score\n(Locus {example_locus})', 
                     fontweight='bold', fontsize=11)
        ax3.axvline(x=example_df['score_cs2g'].iloc[0], color='gray', 
                   linestyle=':', alpha=0.5)
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='red', alpha=0.7, label='True target'),
                          Patch(facecolor='blue', alpha=0.7, label='Non-target')]
        ax3.legend(handles=legend_elements, loc='lower right')
    
    ax3.grid(True, alpha=0.3, axis='x')
    
    # Panel D: Explanation text
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    explanation = """
Why cS2G Shows AUC ≈ 0.50 on RegulatoryBench
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cS2G was designed for GWAS variant-to-gene
assignment, providing a maximum score per gene
aggregated across all linked SNPs.

The Problem:
• cS2G assigns scores to GENES, not POSITIONS
• All candidate genes at the same locus receive
  scores based on their own SNP linkages
• The true target gene is not scored relative
  to the regulatory element position

Result:
• {:.1f}% of loci have candidates with
  identical cS2G scores
• Within-locus discrimination is impossible
• AUC ≈ 0.50 (random performance)

Key Insight:
Methods designed for GWAS-specific positions
cannot be directly applied to arbitrary
regulatory element-to-gene prediction.

This highlights the need for position-aware
methods that score the regulatory→gene link,
not just the gene itself.
""".format(
        100 * (locus_stats['unique_scores'] == 1).mean()
    )
    
    ax4.text(0.05, 0.95, explanation, transform=ax4.transAxes,
            fontfamily='monospace', fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save
    for ext in ['png', 'pdf']:
        output_file = OUTPUT_DIR / f'figure_3_cs2g_analysis.{ext}'
        fig.savefig(output_file, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    plt.close()
    
    # Print statistics
    print("\n" + "=" * 70)
    print("cS2G Statistics")
    print("=" * 70)
    print(f"Loci with only 1 unique score: {(locus_stats['unique_scores'] == 1).sum():,} "
          f"({100 * (locus_stats['unique_scores'] == 1).mean():.1f}%)")
    print(f"Median unique scores per locus: {locus_stats['unique_scores'].median():.0f}")
    print(f"Mean positive score: {pos_scores.mean():.4f}")
    print(f"Mean negative score: {neg_scores.mean():.4f}")


if __name__ == "__main__":
    main()
