#!/usr/bin/env python3
"""
Generate Calibration Curve Figure (Extended Data Fig. 3)
Shows whether predicted scores correspond to true probabilities for L2G and PoPS
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json

# Set publication-quality style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for Nature

def calculate_calibration(y_true, y_pred, n_bins=10):
    """
    Calculate calibration curve using equal-frequency binning
    Returns: bin_means (predicted prob), bin_frequencies (observed prob), bin_counts
    """
    # Remove any NaN values
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    
    if len(y_pred) == 0:
        return np.array([]), np.array([]), np.array([])
    
    # Sort by predicted probability
    sorted_indices = np.argsort(y_pred)
    y_true_sorted = y_true[sorted_indices]
    y_pred_sorted = y_pred[sorted_indices]
    
    # Create equal-frequency bins
    bin_edges = np.percentile(y_pred_sorted, np.linspace(0, 100, n_bins + 1))
    bin_edges[-1] = max(y_pred_sorted.max(), bin_edges[-1])  # Ensure last edge covers max
    
    bin_means = []
    bin_frequencies = []
    bin_counts = []
    
    for i in range(n_bins):
        mask = (y_pred_sorted >= bin_edges[i]) & (y_pred_sorted < bin_edges[i + 1])
        if i == n_bins - 1:  # Include upper edge for last bin
            mask = (y_pred_sorted >= bin_edges[i]) & (y_pred_sorted <= bin_edges[i + 1])
        
        if np.sum(mask) > 0:
            bin_means.append(np.mean(y_pred_sorted[mask]))
            bin_frequencies.append(np.mean(y_true_sorted[mask]))
            bin_counts.append(np.sum(mask))
    
    return np.array(bin_means), np.array(bin_frequencies), np.array(bin_counts)

def calculate_ece(y_true, y_pred, n_bins=10):
    """Calculate Expected Calibration Error"""
    bin_means, bin_frequencies, bin_counts = calculate_calibration(y_true, y_pred, n_bins)
    if len(bin_counts) == 0:
        return np.nan
    total = np.sum(bin_counts)
    ece = np.sum(bin_counts * np.abs(bin_means - bin_frequencies)) / total
    return ece

# Load predictions file that has all methods and scores
print("Loading predictions data...")
predictions = pd.read_csv('results/baselines/post2021_predictions_all_methods.tsv', sep='\t')

print(f"Available methods: {predictions['method'].unique()}")

# For calibration curves, we'll use ABC_Only (supervised method) and PoPS (polygenic method)
# These represent different approaches: ABC uses regulatory networks, PoPS uses gene-level features

# Calculate calibration for ABC (supervised functional method)
abc_data = predictions[predictions['method'] == 'ABC_Only'].copy()
abc_data_clean = abc_data[abc_data['true_gene_score'].notna()].copy()

y_true_abc = abc_data_clean['top1_correct'].values
y_pred_abc = abc_data_clean['true_gene_score'].values

print(f"ABC: {len(abc_data_clean)} loci with scores")

# Calculate calibration for ABC
bin_means_abc, bin_freq_abc, bin_counts_abc = calculate_calibration(y_true_abc, y_pred_abc, n_bins=10)
ece_abc = calculate_ece(y_true_abc, y_pred_abc, n_bins=10)

print(f"ABC ECE: {ece_abc:.4f}")

# Calculate calibration for PoPS
pops_data = predictions[predictions['method'] == 'PoPS'].copy()
pops_data_clean = pops_data[pops_data['true_gene_score'].notna()].copy()

y_true_pops = pops_data_clean['top1_correct'].values
y_pred_pops = pops_data_clean['true_gene_score'].values

print(f"PoPS: {len(pops_data_clean)} loci with scores")

# Calculate calibration for PoPS
bin_means_pops, bin_freq_pops, bin_counts_pops = calculate_calibration(y_true_pops, y_pred_pops, n_bins=10)
ece_pops = calculate_ece(y_true_pops, y_pred_pops, n_bins=10)

print(f"PoPS ECE: {ece_pops:.4f}")

has_pops = len(pops_data_clean) > 0

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))

# Plot ABC calibration
ax = axes[0]
ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Perfect calibration')
ax.plot(bin_means_abc, bin_freq_abc, 'o-', color='#2E86AB', linewidth=2, 
        markersize=6, label=f'ABC (ECE={ece_abc:.3f})')

# Add confidence bars based on bin size
for i in range(len(bin_means_abc)):
    n = bin_counts_abc[i]
    p = bin_freq_abc[i]
    if p > 0 and p < 1:
        se = np.sqrt(p * (1 - p) / n)
        ax.errorbar(bin_means_abc[i], bin_freq_abc[i], yerr=1.96*se, 
                    color='#2E86AB', alpha=0.3, capsize=3)

ax.set_xlabel('Predicted probability (ABC score)', fontsize=9)
ax.set_ylabel('Observed frequency (fraction correct)', fontsize=9)
ax.set_title('ABC Calibration', fontsize=10, fontweight='bold')
ax.legend(loc='upper left', frameon=False, fontsize=8)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(alpha=0.2, linewidth=0.5)

# Plot PoPS calibration
ax = axes[1]
ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Perfect calibration')
ax.plot(bin_means_pops, bin_freq_pops, 'o-', color='#A23B72', linewidth=2,
        markersize=6, label=f'PoPS (ECE={ece_pops:.3f})')

for i in range(len(bin_means_pops)):
    n = bin_counts_pops[i]
    p = bin_freq_pops[i]
    if p > 0 and p < 1:
        se = np.sqrt(p * (1 - p) / n)
        ax.errorbar(bin_means_pops[i], bin_freq_pops[i], yerr=1.96*se,
                    color='#A23B72', alpha=0.3, capsize=3)

ax.set_xlabel('Predicted probability (PoPS score)', fontsize=9)
ax.set_ylabel('Observed frequency (fraction correct)', fontsize=9)
ax.set_title('PoPS Calibration', fontsize=10, fontweight='bold')
ax.legend(loc='upper left', frameon=False, fontsize=8)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(alpha=0.2, linewidth=0.5)

plt.tight_layout()

# Save figure
output_path = Path('figures/baselines/Extended_Data_Figure_3.pdf')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nCalibration curve saved to: {output_path}")

# Also save as PNG for quick viewing
plt.savefig(output_path.with_suffix('.png'), dpi=150, bbox_inches='tight')

plt.show()

# Save calibration metrics to JSON
calibration_metrics = {
    'abc': {
        'ece': float(ece_abc),
        'n_bins': 10,
        'n_pairs': int(len(abc_data_clean)),
        'bin_means': bin_means_abc.tolist(),
        'bin_frequencies': bin_freq_abc.tolist(),
        'bin_counts': bin_counts_abc.tolist()
    },
    'pops': {
        'ece': float(ece_pops),
        'n_bins': 10,
        'n_pairs': int(len(pops_data_clean)),
        'bin_means': bin_means_pops.tolist(),
        'bin_frequencies': bin_freq_pops.tolist(),
        'bin_counts': bin_counts_pops.tolist()
    }
}

with open('data/processed/baselines/calibration_metrics.json', 'w') as f:
    json.dump(calibration_metrics, f, indent=2)

print("Calibration metrics saved to: data/processed/baselines/calibration_metrics.json")
print("\nSummary:")
print(f"ABC: ECE = {ece_abc:.4f} (well-calibrated if ECE < 0.05)")
print(f"PoPS: ECE = {ece_pops:.4f}")
