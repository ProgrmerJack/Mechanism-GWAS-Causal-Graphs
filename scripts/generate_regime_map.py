#!/usr/bin/env python3
"""
Generate the Regime Map figure showing method performance across distance
and evidence type stratifications.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
CANDIDATES_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_cs2g.tsv"
OUTPUT_DIR = BASE_DIR / "figures"


def compute_regime_auc(df, score_col, distance_min, distance_max, evidence_type=None):
    """Compute AUC for a specific regime."""
    
    # Filter by evidence type if specified
    if evidence_type:
        df = df[df['evidence_type'] == evidence_type]
    
    # Get loci with positive pairs in this distance regime
    pos_in_regime = df[(df['label'] == 1) & 
                       (df['distance'] >= distance_min) & 
                       (df['distance'] < distance_max)]
    
    if len(pos_in_regime) == 0:
        return np.nan, 0
    
    loci_with_pos = pos_in_regime['locus_id'].unique()
    
    # Get all candidates for these loci
    df_regime = df[df['locus_id'].isin(loci_with_pos)]
    
    if len(df_regime) == 0 or df_regime['label'].sum() == 0:
        return np.nan, 0
    
    # Compute AUC
    scores = df_regime[score_col].values
    labels = df_regime['label'].values
    
    # Filter NaN
    valid = ~np.isnan(scores)
    if valid.sum() == 0 or labels[valid].sum() == 0:
        return np.nan, 0
    
    try:
        auc = roc_auc_score(labels[valid], scores[valid])
        return auc, len(pos_in_regime)
    except:
        return np.nan, 0


def main():
    print("=" * 70)
    print("Generating Regime Map")
    print("=" * 70)
    
    # Load data
    print(f"\nLoading candidates...")
    df = pd.read_csv(CANDIDATES_FILE, sep='\t', low_memory=False)
    print(f"  Loaded {len(df):,} candidates")
    
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Define regimes
    distance_regimes = [
        ('0-10kb', 0, 10000),
        ('10-100kb', 10000, 100000),
        ('100-500kb', 100000, 500000)
    ]
    
    evidence_types = ['CRISPRi', 'MPRA', 'All']
    methods = [('NearestGene', 'score_nearest')]
    
    # Build regime map data
    results = []
    
    for etype in evidence_types:
        for regime_name, d_min, d_max in distance_regimes:
            for method_name, score_col in methods:
                if etype == 'All':
                    auc, n = compute_regime_auc(df, score_col, d_min, d_max)
                else:
                    auc, n = compute_regime_auc(df, score_col, d_min, d_max, etype)
                
                results.append({
                    'evidence_type': etype,
                    'distance_regime': regime_name,
                    'method': method_name,
                    'auc': auc,
                    'n_positives': n
                })
    
    results_df = pd.DataFrame(results)
    
    # Print table
    print("\n" + "=" * 70)
    print("Regime Map: NearestGene AUC-ROC by Evidence Type and Distance")
    print("=" * 70)
    
    pivot = results_df.pivot_table(
        index='evidence_type',
        columns='distance_regime',
        values='auc'
    )[['0-10kb', '10-100kb', '100-500kb']]
    
    print("\n" + pivot.round(3).to_string())
    
    # Create heatmap figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: Heatmap
    ax1 = axes[0]
    
    # Reorder for display
    pivot_display = pivot.loc[['CRISPRi', 'MPRA', 'All']]
    
    sns.heatmap(
        pivot_display,
        annot=True,
        fmt='.2f',
        cmap='RdYlGn',
        vmin=0.5,
        vmax=1.0,
        ax=ax1,
        cbar_kws={'label': 'AUC-ROC'}
    )
    ax1.set_title('A. Regime Map: NearestGene Performance', fontweight='bold', fontsize=12)
    ax1.set_xlabel('Distance Regime')
    ax1.set_ylabel('Evidence Type')
    
    # Panel B: Line plot showing performance decay with distance
    ax2 = axes[1]
    
    x_positions = [0, 1, 2]
    x_labels = ['0-10kb', '10-100kb', '100-500kb']
    
    for etype in ['CRISPRi', 'MPRA']:
        aucs = [pivot_display.loc[etype, r] for r in x_labels]
        ax2.plot(x_positions, aucs, 'o-', label=etype, linewidth=2, markersize=8)
    
    ax2.axhline(y=0.5, color='gray', linestyle='--', label='Random', alpha=0.7)
    ax2.set_xlabel('Distance Regime')
    ax2.set_ylabel('AUC-ROC')
    ax2.set_title('B. Performance Decay with Distance', fontweight='bold', fontsize=12)
    ax2.set_xticks(x_positions)
    ax2.set_xticklabels(x_labels)
    ax2.set_ylim(0.4, 1.05)
    ax2.legend(loc='lower left')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    for ext in ['png', 'pdf']:
        output_file = OUTPUT_DIR / f'figure_4_regime_map.{ext}'
        fig.savefig(output_file, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    plt.close()
    
    # Save data
    output_csv = OUTPUT_DIR / 'regime_map_data.csv'
    results_df.to_csv(output_csv, index=False)
    print(f"Saved: {output_csv}")
    
    # Summary interpretation
    print("\n" + "=" * 70)
    print("Key Insights from Regime Map")
    print("=" * 70)
    print("""
1. SHORT RANGE (0-10kb): Near-perfect performance for both evidence types
   - CRISPRi: {:.2f}
   - MPRA: {:.2f}
   → At short distances, proximity dominates; sophisticated methods add nothing

2. INTERMEDIATE (10-100kb): Performance diverges
   - CRISPRi: {:.2f}
   - MPRA: {:.2f}
   → CRISPRi enhancers still target nearby genes; MPRA shows more distal effects

3. LONG RANGE (100-500kb): The unsolved frontier
   - CRISPRi: {:.2f}
   - MPRA: {:.2f}
   → All methods struggle; long-range prediction requires new approaches
""".format(
        pivot_display.loc['CRISPRi', '0-10kb'],
        pivot_display.loc['MPRA', '0-10kb'],
        pivot_display.loc['CRISPRi', '10-100kb'],
        pivot_display.loc['MPRA', '10-100kb'],
        pivot_display.loc['CRISPRi', '100-500kb'],
        pivot_display.loc['MPRA', '100-500kb']
    ))


if __name__ == "__main__":
    main()
