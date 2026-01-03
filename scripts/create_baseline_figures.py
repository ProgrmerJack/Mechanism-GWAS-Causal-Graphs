"""
Create publication-ready figures for baseline comparison results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set publication style
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9

# Paths
BASE_DIR = Path(__file__).parent.parent
RESULTS_DIR = BASE_DIR / "results" / "baselines"
FIGURES_DIR = BASE_DIR / "figures" / "baselines"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load all result files"""
    metrics = pd.read_csv(RESULTS_DIR / "post2021_comparison_metrics.tsv", sep='\t')
    by_tier = pd.read_csv(RESULTS_DIR / "post2021_performance_by_tier.tsv", sep='\t')
    by_trait = pd.read_csv(RESULTS_DIR / "post2021_performance_by_trait.tsv", sep='\t')
    predictions = pd.read_csv(RESULTS_DIR / "post2021_predictions_all_methods.tsv", sep='\t')
    
    return metrics, by_tier, by_trait, predictions

def create_overall_performance_barplot(metrics):
    """Figure 1: Bar plot of Top-1/3/5/10 accuracy"""
    
    # Prepare data
    methods = metrics['method'].values
    top1 = metrics['top1_accuracy'].values * 100
    top3 = metrics['top3_accuracy'].values * 100
    top5 = metrics['top5_accuracy'].values * 100
    top10 = metrics['top10_accuracy'].values * 100
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(methods))
    width = 0.2
    
    bars1 = ax.bar(x - 1.5*width, top1, width, label='Top-1', color='#d62728')
    bars2 = ax.bar(x - 0.5*width, top3, width, label='Top-3', color='#ff7f0e')
    bars3 = ax.bar(x + 0.5*width, top5, width, label='Top-5', color='#2ca02c')
    bars4 = ax.bar(x + 1.5*width, top10, width, label='Top-10', color='#1f77b4')
    
    # Formatting
    ax.set_ylabel('Accuracy (%)', fontweight='bold')
    ax.set_xlabel('Method', fontweight='bold')
    ax.set_title('Baseline Method Performance on Post-2021 Benchmark (n=63 loci)', 
                 fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha='right')
    ax.legend(loc='upper right', framealpha=0.9)
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels on bars
    for bars in [bars1, bars2, bars3, bars4]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.1f}',
                       ha='center', va='bottom', fontsize=7)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Figure1_Overall_Performance.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "Figure1_Overall_Performance.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 1 saved: Overall performance bar plot")

def create_mrr_comparison(metrics):
    """Figure 2: Mean Reciprocal Rank comparison"""
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    methods = metrics['method'].values
    mrr = metrics['mrr'].values
    
    # Color code by performance
    colors = ['#d62728' if m == 'Distance' else '#1f77b4' for m in methods]
    
    bars = ax.barh(methods, mrr, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Formatting
    ax.set_xlabel('Mean Reciprocal Rank (MRR)', fontweight='bold')
    ax.set_title('Mean Reciprocal Rank Comparison\nHigher = Better Gene Prioritization', 
                 fontweight='bold', pad=15)
    ax.set_xlim(0, 1.0)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, mrr)):
        ax.text(value + 0.02, bar.get_y() + bar.get_height()/2., 
               f'{value:.3f}',
               ha='left', va='center', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Figure2_MRR_Comparison.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "Figure2_MRR_Comparison.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 2 saved: MRR comparison")

def create_performance_by_tier_heatmap(by_tier):
    """Figure 3: Performance by evidence tier heatmap"""
    
    # Pivot for heatmap (methods × tiers, value = Top-1 accuracy)
    pivot = by_tier.pivot(index='method', columns='evidence_tier', values='top1_accuracy')
    
    # Reorder tiers
    tier_order = ['Tier1_Mendelian', 'Tier1_Coding', 'Tier1_CRISPR', 'Tier1_Drug', 'Tier2_MultiEvidence']
    pivot = pivot[[t for t in tier_order if t in pivot.columns]]
    
    # Reorder methods (Distance first)
    method_order = ['Distance', 'eQTL_Only', 'PoPS', 'CS2G_Proxy', 'FLAMES', 'ABC_Only']
    pivot = pivot.loc[[m for m in method_order if m in pivot.index]]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    sns.heatmap(pivot * 100, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': 'Top-1 Accuracy (%)'},
                linewidths=0.5, linecolor='gray',
                vmin=0, vmax=100, ax=ax)
    
    ax.set_title('Baseline Performance by Evidence Tier\n(Top-1 Accuracy %)', 
                 fontweight='bold', pad=15)
    ax.set_xlabel('Evidence Tier', fontweight='bold')
    ax.set_ylabel('Method', fontweight='bold')
    
    # Rotate labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.setp(ax.get_yticklabels(), rotation=0)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Figure3_Performance_By_Tier.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "Figure3_Performance_By_Tier.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 3 saved: Performance by evidence tier heatmap")

def create_topk_accuracy_curves(metrics, predictions):
    """Figure 4: Top-k accuracy curves"""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    methods = metrics['method'].unique()
    
    # Line styles and colors
    colors = {'Distance': '#d62728', 'ABC_Only': '#ff7f0e', 'eQTL_Only': '#2ca02c',
              'PoPS': '#1f77b4', 'CS2G_Proxy': '#9467bd', 'FLAMES': '#8c564b'}
    
    markers = {'Distance': 'o', 'ABC_Only': 's', 'eQTL_Only': '^',
               'PoPS': 'D', 'CS2G_Proxy': 'v', 'FLAMES': 'p'}
    
    k_values = [1, 3, 5, 10]
    
    for method in methods:
        row = metrics[metrics['method'] == method].iloc[0]
        accuracies = [
            row['top1_accuracy'] * 100,
            row['top3_accuracy'] * 100,
            row['top5_accuracy'] * 100,
            row['top10_accuracy'] * 100
        ]
        
        ax.plot(k_values, accuracies, 
               marker=markers.get(method, 'o'),
               color=colors.get(method, 'gray'),
               label=method, linewidth=2, markersize=8)
    
    # Formatting
    ax.set_xlabel('Top-k', fontweight='bold')
    ax.set_ylabel('Accuracy (%)', fontweight='bold')
    ax.set_title('Top-k Accuracy Curves for Baseline Methods', 
                 fontweight='bold', pad=15)
    ax.set_xticks(k_values)
    ax.set_xticklabels([f'Top-{k}' for k in k_values])
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower right', framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Figure4_TopK_Curves.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "Figure4_TopK_Curves.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 4 saved: Top-k accuracy curves")

def create_rank_distribution_boxplot(predictions):
    """Figure 5: True gene rank distribution boxplot"""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data
    methods = predictions['method'].unique()
    
    data = [predictions[predictions['method'] == m]['true_gene_rank'].values 
            for m in methods]
    
    # Create boxplot
    bp = ax.boxplot(data, labels=methods, patch_artist=True,
                    showfliers=True, notch=True)
    
    # Color boxes
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd', '#8c564b']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Formatting
    ax.set_ylabel('True Gene Rank', fontweight='bold')
    ax.set_xlabel('Method', fontweight='bold')
    ax.set_title('Distribution of True Gene Rankings\n(Lower = Better)', 
                 fontweight='bold', pad=15)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_yscale('log')
    
    # Add median lines
    medians = [np.median(d) for d in data]
    for i, median in enumerate(medians):
        ax.text(i+1.1, median, f'{median:.0f}', 
               ha='left', va='center', fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Figure5_Rank_Distribution.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "Figure5_Rank_Distribution.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 5 saved: Rank distribution boxplot")

def create_trait_category_heatmap(by_trait):
    """Supplementary Figure: Performance by trait category"""
    
    # Filter to top trait categories (≥3 loci)
    trait_counts = by_trait.groupby('trait_category')['n_loci'].first()
    top_traits = trait_counts[trait_counts >= 3].sort_values(ascending=False).index
    
    filtered = by_trait[by_trait['trait_category'].isin(top_traits)]
    
    # Pivot
    pivot = filtered.pivot(index='method', columns='trait_category', values='top1_accuracy')
    
    # Reorder methods
    method_order = ['Distance', 'eQTL_Only', 'PoPS', 'CS2G_Proxy', 'FLAMES', 'ABC_Only']
    pivot = pivot.loc[[m for m in method_order if m in pivot.index]]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    sns.heatmap(pivot * 100, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Top-1 Accuracy (%)'},
                linewidths=0.5, linecolor='gray',
                vmin=0, vmax=100, ax=ax)
    
    ax.set_title('Baseline Performance by Trait Category\n(Top-1 Accuracy %, n≥3 loci)', 
                 fontweight='bold', pad=15)
    ax.set_xlabel('Trait Category', fontweight='bold')
    ax.set_ylabel('Method', fontweight='bold')
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.setp(ax.get_yticklabels(), rotation=0)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "FigureS1_Performance_By_Trait.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "FigureS1_Performance_By_Trait.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Supplementary Figure 1 saved: Performance by trait category")

def create_method_comparison_scatter(predictions):
    """Supplementary Figure: Distance vs ABC/eQTL scatter"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Distance vs ABC
    distance_ranks = predictions[predictions['method'] == 'Distance']['true_gene_rank'].values
    abc_ranks = predictions[predictions['method'] == 'ABC_Only']['true_gene_rank'].values
    
    ax1.scatter(distance_ranks, abc_ranks, alpha=0.6, s=50, color='#1f77b4', edgecolor='black', linewidth=0.5)
    ax1.plot([1, max(distance_ranks.max(), abc_ranks.max())], 
            [1, max(distance_ranks.max(), abc_ranks.max())], 
            'k--', alpha=0.5, linewidth=1)
    ax1.set_xlabel('Distance Rank', fontweight='bold')
    ax1.set_ylabel('ABC Rank', fontweight='bold')
    ax1.set_title('Distance vs ABC Prioritization', fontweight='bold')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    
    # Distance vs eQTL
    eqtl_ranks = predictions[predictions['method'] == 'eQTL_Only']['true_gene_rank'].values
    
    ax2.scatter(distance_ranks, eqtl_ranks, alpha=0.6, s=50, color='#2ca02c', edgecolor='black', linewidth=0.5)
    ax2.plot([1, max(distance_ranks.max(), eqtl_ranks.max())], 
            [1, max(distance_ranks.max(), eqtl_ranks.max())], 
            'k--', alpha=0.5, linewidth=1)
    ax2.set_xlabel('Distance Rank', fontweight='bold')
    ax2.set_ylabel('eQTL Rank', fontweight='bold')
    ax2.set_title('Distance vs eQTL Prioritization', fontweight='bold')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "FigureS2_Method_Comparison_Scatter.pdf", bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "FigureS2_Method_Comparison_Scatter.png", bbox_inches='tight')
    plt.close()
    
    print("✓ Supplementary Figure 2 saved: Method comparison scatter plots")

def main():
    """Create all figures"""
    
    print("Loading data...")
    metrics, by_tier, by_trait, predictions = load_data()
    
    print(f"\nCreating figures...")
    print(f"Output directory: {FIGURES_DIR}")
    print()
    
    # Main figures
    create_overall_performance_barplot(metrics)
    create_mrr_comparison(metrics)
    create_performance_by_tier_heatmap(by_tier)
    create_topk_accuracy_curves(metrics, predictions)
    create_rank_distribution_boxplot(predictions)
    
    # Supplementary figures
    create_trait_category_heatmap(by_trait)
    create_method_comparison_scatter(predictions)
    
    print("\n" + "="*70)
    print("✓ All figures created successfully!")
    print(f"✓ Figures saved to: {FIGURES_DIR}")
    print("="*70)
    
    print("\nFigure Summary:")
    print("  Main Figures:")
    print("    Figure 1: Overall performance bar plot (Top-1/3/5/10)")
    print("    Figure 2: MRR comparison")
    print("    Figure 3: Performance by evidence tier heatmap")
    print("    Figure 4: Top-k accuracy curves")
    print("    Figure 5: Rank distribution boxplot")
    print("  Supplementary Figures:")
    print("    Figure S1: Performance by trait category heatmap")
    print("    Figure S2: Method comparison scatter plots")

if __name__ == "__main__":
    main()
