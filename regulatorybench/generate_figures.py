#!/usr/bin/env python3
"""
Generate Publication Figures for Nature Genetics Analysis

Creates 6 display items as required by NG Analysis format:
1. Figure 1: Task taxonomy schematic (conceptual - just creates data)
2. Figure 2: Task A results - Distance dominance
3. Figure 3: Task B results - ABC vs baselines
4. Figure 4: Leakage audit summary
5. Table 1: Method comparison by task
6. Table 2: Controlled evaluation results
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

# Use script directory for relative paths
SCRIPT_DIR = Path(__file__).parent.resolve()

# Set publication style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

OUTPUT_DIR = SCRIPT_DIR / "figures"
OUTPUT_DIR.mkdir(exist_ok=True)


def load_results():
    """Load all evaluation results."""
    results = {}
    
    task_a_path = SCRIPT_DIR / "benchmarks/task_a_evaluation_results.json"
    if task_a_path.exists():
        with open(task_a_path) as f:
            results['task_a'] = json.load(f)
    
    task_b_path = SCRIPT_DIR / "benchmarks/task_b_baseline_results.json"
    if task_b_path.exists():
        with open(task_b_path) as f:
            results['task_b'] = json.load(f)
    
    leakage_b = SCRIPT_DIR / "benchmarks/task_b_leakage_audit.json"
    if leakage_b.exists():
        with open(leakage_b) as f:
            results['leakage_b'] = json.load(f)
    
    return results


def figure2_task_a_results(results):
    """
    Figure 2: Task A (GWAS→Gene) Results
    
    Shows distance dominance in GWAS locus-to-gene mapping.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    
    methods = results['task_a']['methods']
    
    # Sort by AUROC and take top 10
    methods_sorted = sorted(methods, key=lambda x: x['auroc'], reverse=True)[:12]
    
    names = [m['name'] for m in methods_sorted]
    aurocs = [m['auroc'] for m in methods_sorted]
    categories = [m['category'] for m in methods_sorted]
    
    # Color by category
    category_colors = {
        'Distance': '#E74C3C',
        'GWAS-derived': '#3498DB',
        'ABC': '#2ECC71',
        'Enhancer': '#9B59B6',
        'Regulatory': '#F39C12',
        'eQTL': '#1ABC9C',
        'scATAC': '#E91E63',
        'Expression': '#00BCD4',
        'Functional': '#607D8B',
        'Hi-C': '#FF5722'
    }
    
    colors = [category_colors.get(c, '#95A5A6') for c in categories]
    
    bars = ax.barh(range(len(names)), aurocs, color=colors, edgecolor='black', linewidth=0.5)
    
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.set_xlabel('AUROC', fontweight='bold')
    ax.set_title('Task A: GWAS Credible Set → Causal Gene\n(14,016 pairs, 560 loci)', fontweight='bold')
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5, label='Random')
    ax.set_xlim(0.4, 1.0)
    ax.invert_yaxis()
    
    # Add value labels
    for bar, auroc in zip(bars, aurocs):
        ax.text(auroc + 0.01, bar.get_y() + bar.get_height()/2, 
                f'{auroc:.3f}', va='center', fontsize=8)
    
    # Legend for categories
    unique_cats = list(dict.fromkeys(categories))  # Preserve order
    legend_patches = [mpatches.Patch(color=category_colors.get(c, '#95A5A6'), label=c) 
                      for c in unique_cats]
    ax.legend(handles=legend_patches, loc='lower right', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'figure2_task_a_results.png')
    fig.savefig(OUTPUT_DIR / 'figure2_task_a_results.pdf')
    plt.close()
    print(f"Saved Figure 2: Task A results")


def figure3_task_b_results(results):
    """
    Figure 3: Task B (Enhancer→Gene) Results
    
    Shows ABC vs distance comparison.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    # Panel A: Main comparison
    methods = results['task_b']['methods']
    names = [m['name'].replace('_', ' ') for m in methods]
    aurocs = [m['auroc'] for m in methods if m['auroc'] is not None]
    auprcs = [m['auprc'] for m in methods if m['auprc'] is not None]
    
    x = np.arange(len(names))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, aurocs, width, label='AUROC', color='#3498DB', edgecolor='black')
    bars2 = ax1.bar(x + width/2, auprcs, width, label='AUPRC', color='#E74C3C', edgecolor='black')
    
    ax1.set_ylabel('Score', fontweight='bold')
    ax1.set_title('Task B: Enhancer → Gene\n(19,825 pairs)', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, rotation=15, ha='right')
    ax1.legend()
    ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax1.set_ylim(0, 1)
    
    # Add value labels
    for bar, val in zip(bars1, aurocs):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                f'{val:.3f}', ha='center', va='bottom', fontsize=8)
    
    # Panel B: Distance stratification
    strata = results['task_b']['distance_strata']
    strata_names = list(strata.keys())
    strata_aurocs = [strata[s]['auroc'] if strata[s]['auroc'] else 0.5 for s in strata_names]
    strata_n = [strata[s]['n_pairs'] for s in strata_names]
    
    colors = ['#E74C3C', '#F39C12', '#3498DB']
    bars = ax2.bar(strata_names, strata_aurocs, color=colors, edgecolor='black')
    
    ax2.set_ylabel('Distance Baseline AUROC', fontweight='bold')
    ax2.set_xlabel('Distance to TSS', fontweight='bold')
    ax2.set_title('Distance Baseline Performance\nby Distance Stratum', fontweight='bold')
    ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax2.set_ylim(0.5, 0.8)
    
    # Add sample sizes
    for bar, n, auroc in zip(bars, strata_n, strata_aurocs):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'n={n:,}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'figure3_task_b_results.png')
    fig.savefig(OUTPUT_DIR / 'figure3_task_b_results.pdf')
    plt.close()
    print(f"Saved Figure 3: Task B results")


def figure4_leakage_audit(results):
    """
    Figure 4: Leakage Audit Summary
    
    Shows data contamination issues in Task B benchmark.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Panel A: Leakage types (only show first 3 audits for clarity)
    audits = results['leakage_b']['audits']
    # Filter to only the 3 main leakage types for visualization
    main_audits = [a for a in audits if 'Gene ID' not in a.get('name', '')][:3]
    audit_names = ['ABC Training\nOverlap', 'Gene Overlap\nAcross Sources', 'Genomic\nProximity']
    affected_pct = [a['affected'] / a['total'] * 100 for a in main_audits]
    
    colors = ['#E74C3C' if pct > 50 else '#F39C12' if pct > 20 else '#2ECC71' for pct in affected_pct]
    bars = ax1.bar(audit_names, affected_pct, color=colors, edgecolor='black')
    
    ax1.set_ylabel('Affected Pairs (%)', fontweight='bold')
    ax1.set_title('Task B Leakage Audit Results\n(N=19,825 pairs)', fontweight='bold')
    ax1.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% threshold')
    ax1.set_ylim(0, 100)
    
    for bar, pct in zip(bars, affected_pct):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, 
                f'{pct:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Panel B: Leave-one-study-out splits
    loso = results['leakage_b']['loso_splits']
    study_names = ['K562\n(held out)', 'Other ENCODE\n(held out)', 'Fulco 2019\n(held out)']
    test_sizes = [loso['holdout_ENCODE_CRISPRi_K562']['test_size'],
                  loso['holdout_ENCODE_CRISPRi_heldout']['test_size'],
                  loso['holdout_Fulco_2019']['test_size']]
    test_pos = [loso['holdout_ENCODE_CRISPRi_K562']['test_pos'],
                loso['holdout_ENCODE_CRISPRi_heldout']['test_pos'],
                loso['holdout_Fulco_2019']['test_pos']]
    
    x = np.arange(len(study_names))
    width = 0.35
    
    bars1 = ax2.bar(x - width/2, test_sizes, width, label='Test pairs', color='#3498DB', edgecolor='black')
    bars2 = ax2.bar(x + width/2, test_pos, width, label='Positives', color='#E74C3C', edgecolor='black')
    
    ax2.set_ylabel('Count', fontweight='bold')
    ax2.set_title('Leave-One-Study-Out Splits\nfor Unbiased Evaluation', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(study_names)
    ax2.legend()
    
    # Add positive rates
    for i, (total, pos) in enumerate(zip(test_sizes, test_pos)):
        rate = pos / total * 100
        ax2.text(i, max(total, pos) + 300, f'{rate:.1f}%', ha='center', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'figure4_leakage_audit.png')
    fig.savefig(OUTPUT_DIR / 'figure4_leakage_audit.pdf')
    plt.close()
    print(f"Saved Figure 4: Leakage audit")


def table1_method_comparison(results):
    """
    Table 1: Method Comparison by Task Type
    """
    rows = []
    
    # Task A methods
    for m in results['task_a']['methods'][:10]:
        rows.append({
            'Task': 'A (GWAS→Gene)',
            'Method': m['name'],
            'Category': m['category'],
            'AUROC': f"{m['auroc']:.3f}",
            'Coverage': f"{m['coverage']*100:.0f}%",
            'Applicable': 'Yes'
        })
    
    # Task B methods
    for m in results['task_b']['methods']:
        if m['auroc'] is not None:
            rows.append({
                'Task': 'B (Enhancer→Gene)',
                'Method': m['name'].replace('_', ' '),
                'Category': 'Baseline',
                'AUROC': f"{m['auroc']:.3f}",
                'Coverage': f"{m['coverage']*100:.0f}%",
                'Applicable': 'Yes'
            })
    
    df = pd.DataFrame(rows)
    
    # Save as CSV
    df.to_csv(OUTPUT_DIR / 'table1_method_comparison.csv', index=False)
    
    print(f"Saved Table 1: Method comparison ({len(df)} rows)")
    return df


def table2_controlled_evaluation(results):
    """
    Table 2: Leakage-Controlled Evaluation Recommendations
    """
    rows = [
        {
            'Issue': 'ABC Training Overlap',
            'Affected': '78%',
            'Risk Level': 'CRITICAL',
            'Mitigation': 'Exclude K562/Fulco when evaluating ABC',
            'Impact': 'ABC performance likely inflated'
        },
        {
            'Issue': 'Gene Overlap',
            'Affected': '64%',
            'Risk Level': 'HIGH',
            'Mitigation': 'Leave-one-study-out CV',
            'Impact': 'Gene-specific features may transfer'
        },
        {
            'Issue': 'Genomic Proximity',
            'Affected': '92%',
            'Risk Level': 'MODERATE',
            'Mitigation': 'Chromosome holdout (chr8,9)',
            'Impact': 'Nearby loci may share causal genes'
        }
    ]
    
    df = pd.DataFrame(rows)
    
    df.to_csv(OUTPUT_DIR / 'table2_leakage_controls.csv', index=False)
    
    print(f"Saved Table 2: Leakage controls")
    return df


def figure1_task_taxonomy():
    """
    Figure 1: Task Taxonomy Schematic
    
    Conceptual diagram showing Task A/B/C definitions.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Draw boxes for each task type
    box_props = dict(boxstyle='round,pad=0.5', facecolor='lightblue', edgecolor='black', linewidth=2)
    
    # Task A
    ax.add_patch(plt.Rectangle((0.05, 0.65), 0.28, 0.3, facecolor='#E8F4FD', edgecolor='#2C3E50', linewidth=2, zorder=1))
    ax.text(0.19, 0.90, 'TASK A', fontsize=14, fontweight='bold', ha='center', va='center')
    ax.text(0.19, 0.82, 'GWAS → Gene', fontsize=11, ha='center', va='center')
    ax.text(0.19, 0.73, 'Input: Credible Set\nOutput: Causal Gene\nMethods: L2G, PoPS, FLAMES', 
            fontsize=8, ha='center', va='center', linespacing=1.5)
    
    # Task B
    ax.add_patch(plt.Rectangle((0.36, 0.65), 0.28, 0.3, facecolor='#E8FDF4', edgecolor='#27AE60', linewidth=2, zorder=1))
    ax.text(0.50, 0.90, 'TASK B', fontsize=14, fontweight='bold', ha='center', va='center')
    ax.text(0.50, 0.82, 'Enhancer → Gene', fontsize=11, ha='center', va='center')
    ax.text(0.50, 0.73, 'Input: Enhancer Region\nOutput: Target Gene(s)\nMethods: ABC, rE2G, Hi-C', 
            fontsize=8, ha='center', va='center', linespacing=1.5)
    
    # Task C
    ax.add_patch(plt.Rectangle((0.67, 0.65), 0.28, 0.3, facecolor='#FDE8F4', edgecolor='#8E44AD', linewidth=2, zorder=1))
    ax.text(0.81, 0.90, 'TASK C', fontsize=14, fontweight='bold', ha='center', va='center')
    ax.text(0.81, 0.82, 'Variant → Gene', fontsize=11, ha='center', va='center')
    ax.text(0.81, 0.73, 'Input: SNP Position\nOutput: Affected Gene(s)\nMethods: eQTL, VEP', 
            fontsize=8, ha='center', va='center', linespacing=1.5)
    
    # Arrows showing relationships
    ax.annotate('', xy=(0.36, 0.80), xytext=(0.33, 0.80),
                arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    ax.annotate('', xy=(0.67, 0.80), xytext=(0.64, 0.80),
                arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    
    # Key insight box
    ax.add_patch(plt.Rectangle((0.15, 0.15), 0.70, 0.35, facecolor='#FFF3E0', edgecolor='#E65100', linewidth=2, zorder=1))
    ax.text(0.50, 0.42, 'KEY INSIGHT', fontsize=12, fontweight='bold', ha='center', va='center', color='#E65100')
    ax.text(0.50, 0.30, 
            'Published benchmarks often conflate these distinct tasks.\n'
            'Task A (GWAS→Gene): Distance dominates (AUROC=0.93)\n'
            'Task B (E2G): ABC outperforms distance (AUROC=0.88 vs 0.87)\n'
            'Methods should only be evaluated on their applicable task.',
            fontsize=9, ha='center', va='center', linespacing=1.8)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('RegulatoryBench v4: Task-Stratified Evaluation Framework', fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'figure1_task_taxonomy.png')
    fig.savefig(OUTPUT_DIR / 'figure1_task_taxonomy.pdf')
    plt.close()
    print(f"Saved Figure 1: Task taxonomy")


def main():
    """Generate all publication figures."""
    print("\nGenerating Publication Figures for Nature Genetics Analysis")
    print("=" * 60)
    
    results = load_results()
    
    if not results:
        print("ERROR: No results files found. Run evaluation scripts first.")
        return
    
    # Generate figures
    figure1_task_taxonomy()
    
    if 'task_a' in results:
        figure2_task_a_results(results)
    
    if 'task_b' in results:
        figure3_task_b_results(results)
    
    if 'leakage_b' in results:
        figure4_leakage_audit(results)
    
    # Generate tables
    if 'task_a' in results and 'task_b' in results:
        table1_method_comparison(results)
    
    if 'leakage_b' in results:
        table2_controlled_evaluation(results)
    
    print("\n" + "=" * 60)
    print(f"All figures saved to: {OUTPUT_DIR}")
    print("\nDisplay items for NG Analysis:")
    print("  1. Figure 1: Task taxonomy (conceptual)")
    print("  2. Figure 2: Task A results")
    print("  3. Figure 3: Task B results")
    print("  4. Figure 4: Leakage audit")
    print("  5. Table 1: Method comparison")
    print("  6. Table 2: Leakage controls")


if __name__ == "__main__":
    main()
