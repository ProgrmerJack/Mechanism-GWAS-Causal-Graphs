"""
Generate stratified regime map showing method performance by distance bins.

This addresses External Critique Attack A: "Gold standards biased toward proximity"
by showing that no single method dominates across all distance regimes.

Distance Bins:
- 0-10kb: Coding variants (proximal)
- 10-100kb: Proximal regulatory
- 100kb-1Mb: Distal regulatory  
- >1Mb: Ultra-distal/complex

Expected Pattern:
- Distance method should dominate at <10kb
- Functional methods (ABC, PoPS) should improve at >50kb
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts for publication
matplotlib.rcParams['font.family'] = 'Arial'
import json
from pathlib import Path

# Load data
print("Loading data...")
# Load locus-gene pairs with distances
pairs_df = pd.read_csv('data/processed/baselines/post2021_locus_gene_pairs_annotated.tsv', sep='\t')

# Load predictions with performance metrics
predictions_df = pd.read_csv('results/baselines/post2021_predictions_all_methods.tsv', sep='\t')

# Merge to get distance information for true genes
merged_df = predictions_df.merge(
    pairs_df[['locus_id', 'gene', 'distance_to_lead']],
    left_on=['locus_id', 'true_gene'],
    right_on=['locus_id', 'gene'],
    how='left'
)

print(f"Merged {len(merged_df)} predictions with distance information")
print(f"Missing distances: {merged_df['distance_to_lead'].isna().sum()}")

# Define distance bins
def classify_distance_regime(distance_bp):
    """Classify distance into regime bins."""
    if pd.isna(distance_bp):
        return 'Unknown'
    elif distance_bp < 10_000:
        return '0-10kb\n(Coding)'
    elif distance_bp < 100_000:
        return '10-100kb\n(Proximal)'
    elif distance_bp < 1_000_000:
        return '100kb-1Mb\n(Distal)'
    else:
        return '>1Mb\n(Ultra-distal)'

merged_df['distance_regime'] = merged_df['distance_to_lead'].apply(classify_distance_regime)

# Remove unknown distances
merged_df = merged_df[merged_df['distance_regime'] != 'Unknown']

print("\nDistance regime distribution:")
print(merged_df.groupby(['distance_regime', 'method']).size().unstack(fill_value=0))

# Calculate performance metrics by regime and method
regime_order = ['0-10kb\n(Coding)', '10-100kb\n(Proximal)', '100kb-1Mb\n(Distal)', '>1Mb\n(Ultra-distal)']
methods = ['Distance', 'ABC_Only', 'PoPS', 'cS2G_LocusAware_max', 'FLAMES', 'eQTL_Only']

# Calculate Top-1 accuracy for each regime × method
results = []
for regime in regime_order:
    for method in methods:
        subset = merged_df[(merged_df['distance_regime'] == regime) & (merged_df['method'] == method)]
        
        if len(subset) == 0:
            continue
        
        top1_accuracy = subset['top1_correct'].mean() * 100
        n = len(subset)
        
        # 95% Confidence interval (Wilson score)
        if n > 0:
            p = top1_accuracy / 100
            z = 1.96
            denominator = 1 + z**2 / n
            center = (p + z**2 / (2*n)) / denominator
            margin = z * np.sqrt((p * (1 - p) / n + z**2 / (4*n**2)) / denominator)
            ci_lower = max(0, (center - margin) * 100)
            ci_upper = min(100, (center + margin) * 100)
        else:
            ci_lower = ci_upper = top1_accuracy
        
        results.append({
            'regime': regime,
            'method': method,
            'top1_accuracy': top1_accuracy,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'n_loci': n
        })

results_df = pd.DataFrame(results)

print("\nPerformance by regime:")
print(results_df.pivot(index='method', columns='regime', values='top1_accuracy').round(1))

# Create publication-quality figure
fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True)

# Color scheme for methods
method_colors = {
    'Distance': '#1f77b4',
    'ABC_Only': '#ff7f0e',
    'PoPS': '#2ca02c',
    'cS2G_LocusAware_max': '#d62728',
    'FLAMES': '#9467bd',
    'eQTL_Only': '#8c564b'
}

method_labels = {
    'Distance': 'Distance',
    'ABC_Only': 'ABC',
    'PoPS': 'PoPS',
    'cS2G_LocusAware_max': 'cS2G',
    'FLAMES': 'FLAMES',
    'eQTL_Only': 'eQTL'
}

for i, regime in enumerate(regime_order):
    ax = axes[i]
    
    regime_data = results_df[results_df['regime'] == regime].sort_values('top1_accuracy', ascending=False)
    
    # Plot bars with error bars
    x_positions = np.arange(len(regime_data))
    bars = ax.bar(
        x_positions,
        regime_data['top1_accuracy'],
        color=[method_colors[m] for m in regime_data['method']],
        alpha=0.8,
        edgecolor='black',
        linewidth=0.5
    )
    
    # Add 95% CI error bars
    errors_lower = regime_data['top1_accuracy'] - regime_data['ci_lower']
    errors_upper = regime_data['ci_upper'] - regime_data['top1_accuracy']
    ax.errorbar(
        x_positions,
        regime_data['top1_accuracy'],
        yerr=[errors_lower, errors_upper],
        fmt='none',
        ecolor='black',
        capsize=3,
        capthick=1,
        linewidth=1
    )
    
    # Add sample size annotations
    for j, (idx, row) in enumerate(regime_data.iterrows()):
        ax.text(
            j, row['top1_accuracy'] + 5,
            f"n={row['n_loci']}",
            ha='center',
            va='bottom',
            fontsize=7
        )
    
    # Labels
    ax.set_title(regime, fontsize=11, fontweight='bold')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([method_labels[m] for m in regime_data['method']], rotation=45, ha='right', fontsize=9)
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    if i == 0:
        ax.set_ylabel('Top-1 Accuracy (%)', fontsize=11, fontweight='bold')
    
    # Add reference line at 50%
    ax.axhline(50, color='gray', linestyle=':', linewidth=1, alpha=0.5)

plt.tight_layout()

# Save figure
output_path = Path('figures/baselines/Extended_Data_Figure_4_Stratified_Regime_Map.pdf')
output_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nStratified regime map saved to: {output_path}")

# Save metrics to JSON
metrics = {
    'regime_performance': results_df.to_dict('records'),
    'key_findings': {
        'coding_regime_winner': results_df[results_df['regime'] == '0-10kb\n(Coding)'].nlargest(1, 'top1_accuracy')['method'].values[0] if len(results_df[results_df['regime'] == '0-10kb\n(Coding)']) > 0 else 'N/A',
        'distal_regime_winner': results_df[results_df['regime'] == '100kb-1Mb\n(Distal)'].nlargest(1, 'top1_accuracy')['method'].values[0] if len(results_df[results_df['regime'] == '100kb-1Mb\n(Distal)']) > 0 else 'N/A',
        'total_loci_analyzed': len(merged_df['locus_id'].unique())
    }
}

metrics_path = Path('data/processed/baselines/stratified_regime_metrics.json')
with open(metrics_path, 'w') as f:
    json.dump(metrics, f, indent=2)

print(f"Metrics saved to: {metrics_path}")

# Print summary
print("\n=== SUMMARY ===")
print(f"Total loci analyzed: {len(merged_df['locus_id'].unique())}")
print("\nDistance regime distribution:")
regime_counts = merged_df.groupby('distance_regime')['locus_id'].nunique()
for regime in regime_order:
    count = regime_counts.get(regime, 0)
    print(f"  {regime}: {count} loci")

print("\nKey findings:")
for regime in regime_order:
    regime_results = results_df[results_df['regime'] == regime].nlargest(1, 'top1_accuracy')
    if len(regime_results) > 0:
        winner = regime_results.iloc[0]
        print(f"  {regime}: {method_labels[winner['method']]} ({winner['top1_accuracy']:.1f}%, n={winner['n_loci']})")
