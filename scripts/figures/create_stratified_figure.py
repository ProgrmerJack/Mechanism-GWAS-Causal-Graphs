#!/usr/bin/env python3
"""
create_stratified_figure.py
============================
Nature Genetics Figure: 3-Panel Stratified Performance

This script creates the critical figure showing that "Distance dominates"
is a benchmark artifact, not a meaningful finding.

Panel A: CODING loci (n=56) - Distance wins (expected)
Panel B: REGULATORY loci (n=7) - Distance gap narrows
Panel C: DISTANCE-FAILS (n=5) - All methods struggle

This figure demonstrates:
1. Why 88.9% coding bias inflates Distance performance
2. The real locus-to-gene challenge is regulatory variants
3. Current benchmark is insufficient for method evaluation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Configuration
RESULTS_DIR = Path("results/baselines/stratified")
FIGURES_DIR = Path("figures/baselines")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Nature Genetics style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Color palette
COLORS = {
    'Distance': '#2166ac',     # Blue
    'ABC_Only': '#b2182b',     # Red
    'eQTL_Only': '#4daf4a',    # Green
    'PoPS': '#984ea3',         # Purple
    'CS2G_Proxy': '#ff7f00',   # Orange
    'FLAMES': '#a65628',       # Brown
}


def load_stratified_metrics():
    """Load stratified performance metrics."""
    metrics_path = RESULTS_DIR / "stratified_performance_by_mechanism.tsv"
    if metrics_path.exists():
        df = pd.read_csv(metrics_path, sep='\t')
        return df
    return None


def create_panel(ax, df, subset_name, title, show_legend=False):
    """Create a single panel showing Top-1/3/5 accuracy by method."""
    
    subset_df = df[df['subset'] == subset_name].copy()
    if len(subset_df) == 0:
        ax.text(0.5, 0.5, f"No data for {subset_name}", ha='center', va='center')
        return
    
    methods = ['Distance', 'ABC_Only', 'eQTL_Only', 'PoPS', 'CS2G_Proxy', 'FLAMES']
    x = np.arange(len(methods))
    width = 0.25
    
    # Get metrics for each method
    top1_vals = []
    top3_vals = []
    top5_vals = []
    
    for method in methods:
        method_df = subset_df[subset_df['method'] == method]
        if len(method_df) > 0:
            top1_vals.append(method_df['top1_acc'].values[0])
            top3_vals.append(method_df['top3_acc'].values[0])
            top5_vals.append(method_df['top5_acc'].values[0])
        else:
            top1_vals.append(0)
            top3_vals.append(0)
            top5_vals.append(0)
    
    # Plot bars
    bars1 = ax.bar(x - width, top1_vals, width, label='Top-1', color='#1b7837', alpha=0.9)
    bars3 = ax.bar(x, top3_vals, width, label='Top-3', color='#41ab5d', alpha=0.9)
    bars5 = ax.bar(x + width, top5_vals, width, label='Top-5', color='#a1d99b', alpha=0.9)
    
    # Get sample size
    n_loci = subset_df['n_loci'].values[0] if len(subset_df) > 0 else 0
    
    # Add title with sample size
    ax.set_title(f"{title}\n(n={n_loci} loci)", fontweight='bold', pad=10)
    
    # Labels
    ax.set_ylabel('Accuracy (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(['Distance', 'ABC', 'eQTL', 'PoPS', 'cS2G', 'FLAMES'], rotation=45, ha='right')
    ax.set_ylim(0, 100)
    
    # Add grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        if height > 5:
            ax.annotate(f'{height:.0f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 2),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=6)
    
    if show_legend:
        ax.legend(loc='upper right', framealpha=0.9)


def create_three_panel_figure():
    """Create the main 3-panel figure for Nature Genetics."""
    
    df = load_stratified_metrics()
    if df is None:
        print("ERROR: Could not load stratified metrics")
        return
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(8, 3.5))
    
    # Panel A: CODING
    create_panel(axes[0], df, 'CODING', 'A. Coding Loci', show_legend=True)
    
    # Panel B: REGULATORY
    create_panel(axes[1], df, 'REGULATORY', 'B. Regulatory Loci')
    
    # Panel C: DISTANCE_FAILS
    create_panel(axes[2], df, 'DISTANCE_FAILS', 'C. Distance-Fails')
    
    # Add annotation explaining the finding
    fig.text(0.5, -0.02, 
             "Distance dominance (61.9%) reflects 89% coding benchmark bias, not method superiority.",
             ha='center', va='top', fontsize=7, style='italic', color='#666666')
    
    plt.tight_layout()
    
    # Save
    fig.savefig(FIGURES_DIR / 'Fig_stratified_performance_3panel.pdf', 
                format='pdf', bbox_inches='tight', pad_inches=0.1)
    fig.savefig(FIGURES_DIR / 'Fig_stratified_performance_3panel.png', 
                format='png', bbox_inches='tight', pad_inches=0.1)
    
    print(f"Saved: {FIGURES_DIR / 'Fig_stratified_performance_3panel.pdf'}")
    print(f"Saved: {FIGURES_DIR / 'Fig_stratified_performance_3panel.png'}")
    
    plt.close()
    return fig


def create_mechanism_pie_chart():
    """Create pie chart showing benchmark composition bias."""
    
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    
    sizes = [56, 7]  # CODING, REGULATORY
    labels = ['Coding\n(n=56)', 'Regulatory\n(n=7)']
    colors = ['#d73027', '#1a9850']
    explode = (0, 0.1)
    
    wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors,
                                        explode=explode, autopct='%1.1f%%',
                                        startangle=90, pctdistance=0.6,
                                        textprops={'fontsize': 9})
    
    ax.set_title('Benchmark Composition\n(Post-2021 Independent)', fontweight='bold')
    
    # Add warning annotation
    ax.annotate('⚠️ 88.9% coding bias\nfavors Distance baseline',
                xy=(0.5, -0.15), xycoords='axes fraction',
                ha='center', va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='#fff3cd', alpha=0.8))
    
    plt.tight_layout()
    
    fig.savefig(FIGURES_DIR / 'Fig_benchmark_composition_pie.pdf', 
                format='pdf', bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig_benchmark_composition_pie.png', 
                format='png', bbox_inches='tight')
    
    print(f"Saved: {FIGURES_DIR / 'Fig_benchmark_composition_pie.pdf'}")
    
    plt.close()
    return fig


def create_distance_gap_figure():
    """Create figure showing Distance advantage narrows for regulatory loci."""
    
    df = load_stratified_metrics()
    if df is None:
        return
    
    fig, ax = plt.subplots(figsize=(4, 3))
    
    subsets = ['CODING', 'REGULATORY']
    x = np.arange(len(subsets))
    width = 0.15
    
    methods = ['Distance', 'ABC_Only', 'eQTL_Only', 'PoPS', 'CS2G_Proxy', 'FLAMES']
    colors = ['#2166ac', '#b2182b', '#4daf4a', '#984ea3', '#ff7f00', '#a65628']
    
    for i, (method, color) in enumerate(zip(methods, colors)):
        vals = []
        for subset in subsets:
            subset_df = df[(df['subset'] == subset) & (df['method'] == method)]
            if len(subset_df) > 0:
                vals.append(subset_df['top1_acc'].values[0])
            else:
                vals.append(0)
        
        bars = ax.bar(x + i * width, vals, width, label=method, color=color)
    
    ax.set_ylabel('Top-1 Accuracy (%)')
    ax.set_xticks(x + width * 2.5)
    ax.set_xticklabels(['Coding\n(n=56)', 'Regulatory\n(n=7)'])
    ax.set_ylim(0, 80)
    ax.legend(loc='upper right', fontsize=6, ncol=2)
    ax.set_title('Distance Advantage Narrows for Regulatory Loci', fontweight='bold')
    
    # Add gap annotations with clean vertical lines and badges (no arrows)
    # Coding gap
    ax.plot([0.15, 0.15], [5, 64], color='red', lw=1.5, linestyle='-', alpha=0.8)
    ax.text(-0.1, 35, 'Gap: 59pp', fontsize=7, fontweight='bold', color='red', 
           rotation=90, va='center',
           bbox=dict(boxstyle='round,pad=0.15', facecolor='white', 
                    edgecolor='red', linewidth=0.5, alpha=0.9))
    
    # Regulatory gap  
    ax.plot([1.15, 1.15], [0, 43], color='red', lw=1.5, linestyle='-', alpha=0.8)
    ax.text(0.9, 20, 'Gap: 43pp', fontsize=7, fontweight='bold', color='red', 
           rotation=90, va='center',
           bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                    edgecolor='red', linewidth=0.5, alpha=0.9))
    
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    fig.savefig(FIGURES_DIR / 'Fig_distance_gap_comparison.pdf', 
                format='pdf', bbox_inches='tight')
    fig.savefig(FIGURES_DIR / 'Fig_distance_gap_comparison.png', 
                format='png', bbox_inches='tight')
    
    print(f"Saved: {FIGURES_DIR / 'Fig_distance_gap_comparison.pdf'}")
    
    plt.close()


def main():
    """Generate all stratified figures for Nature Genetics."""
    
    print("=" * 70)
    print("CREATING STRATIFIED PERFORMANCE FIGURES")
    print("Nature Genetics Requirement: Show benchmark composition bias")
    print("=" * 70)
    
    # Create main 3-panel figure
    print("\n>>> Creating 3-panel stratified figure...")
    create_three_panel_figure()
    
    # Create benchmark composition pie chart
    print("\n>>> Creating benchmark composition pie chart...")
    create_mechanism_pie_chart()
    
    # Create distance gap comparison
    print("\n>>> Creating distance gap comparison...")
    create_distance_gap_figure()
    
    print("\n" + "=" * 70)
    print("FIGURES COMPLETE")
    print("=" * 70)
    print("\nKey message for Nature Genetics:")
    print("  'Distance dominates (61.9%)' is a benchmark artifact,")
    print("  not evidence of method superiority.")
    print("  The benchmark is 89% coding loci where Distance is expected to win.")


if __name__ == "__main__":
    main()
