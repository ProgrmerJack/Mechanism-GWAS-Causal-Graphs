#!/usr/bin/env python3
"""
PROPER APPLES-TO-APPLES CALIBRATION COMPARISON

This script addresses the reviewer's core criticism:
- ECE = 0.683 (L2G v25.12 gene-level collapse) is INVALID
- L2G v25.12 is defined per GWAS credible set, not per gene-disease
- We must use study-locus → disease aligned evaluation

We compute:
1. Full benchmark (n=14,016): Our method vs L2G v22.09
2. Matched subset (n=6,344): Our method vs L2G v25.12 (semantic-aligned)

This produces REVIEWER-PROOF comparisons.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json

def compute_ece(y_true, y_prob, n_bins=10):
    """Expected Calibration Error with quantile binning"""
    y_true = np.array(y_true)
    y_prob = np.array(y_prob)
    
    # Handle edge cases
    if len(y_true) == 0 or len(y_prob) == 0:
        return np.nan
    
    # Use quantile binning for robustness
    try:
        bin_edges = np.percentile(y_prob, np.linspace(0, 100, n_bins + 1))
        bin_edges = np.unique(bin_edges)
    except:
        bin_edges = np.linspace(0, 1, n_bins + 1)
    
    ece = 0.0
    total_samples = len(y_true)
    
    for i in range(len(bin_edges) - 1):
        if i == len(bin_edges) - 2:
            mask = (y_prob >= bin_edges[i]) & (y_prob <= bin_edges[i+1])
        else:
            mask = (y_prob >= bin_edges[i]) & (y_prob < bin_edges[i+1])
        
        if mask.sum() > 0:
            bin_accuracy = y_true[mask].mean()
            bin_confidence = y_prob[mask].mean()
            bin_size = mask.sum()
            ece += (bin_size / total_samples) * abs(bin_accuracy - bin_confidence)
    
    return ece

def compute_brier(y_true, y_prob):
    """Brier score"""
    return np.mean((np.array(y_prob) - np.array(y_true)) ** 2)

def main():
    print("=" * 80)
    print("PROPER APPLES-TO-APPLES CALIBRATION COMPARISON")
    print("Addressing reviewer criticism: L2G v25.12 gene-level collapse is INVALID")
    print("=" * 80)
    
    # Load the semantic-matched predictions file
    matched_file = Path("l2g_semantic_matched_predictions.tsv")
    if not matched_file.exists():
        print("ERROR: l2g_semantic_matched_predictions.tsv not found")
        return
    
    print("\n1. Loading semantic-matched predictions...")
    df = pd.read_csv(matched_file, sep='\t')
    print(f"   Total rows: {len(df)}")
    
    # Identify which rows have L2G v25.12 predictions (semantic-aligned)
    # These have non-empty l2g_max values
    df['has_l2g_v25'] = df['l2g_max'].notna() & (df['l2g_max'] != '')
    
    # Convert columns to proper types
    df['truth'] = df['truth'].astype(bool)
    df['calibrated_prob'] = pd.to_numeric(df['calibrated_prob'], errors='coerce')
    df['l2g_max'] = pd.to_numeric(df['l2g_max'], errors='coerce')
    df['l2g_mean'] = pd.to_numeric(df['l2g_mean'], errors='coerce')
    
    # Get the matched subset (rows where L2G v25.12 has predictions)
    matched_subset = df[df['has_l2g_v25']].copy()
    print(f"   Matched subset (with L2G v25.12): {len(matched_subset)} rows")
    print(f"   Match rate: {100 * len(matched_subset) / len(df):.1f}%")
    
    # Compute metrics on full benchmark
    print("\n" + "=" * 80)
    print("PANEL A: FULL BENCHMARK (n={})".format(len(df)))
    print("Methods with gene-disease native outputs")
    print("=" * 80)
    
    full_valid = df['calibrated_prob'].notna()
    our_full_ece = compute_ece(df.loc[full_valid, 'truth'], df.loc[full_valid, 'calibrated_prob'])
    our_full_brier = compute_brier(df.loc[full_valid, 'truth'], df.loc[full_valid, 'calibrated_prob'])
    our_full_n = full_valid.sum()
    
    print(f"\n   Our Method (Mechanism Graphs):")
    print(f"     ECE = {our_full_ece:.4f}")
    print(f"     Brier = {our_full_brier:.4f}")
    print(f"     n = {our_full_n}")
    
    # Load L2G v22.09 from the original benchmark file
    benchmark_file = Path("validation_bundle/calibration/gene_predictions_with_l2g_versions.tsv")
    if benchmark_file.exists():
        bench_df = pd.read_csv(benchmark_file, sep='\t')
        v22_valid = bench_df['l2g_v22_max'].notna()
        if v22_valid.any():
            v22_ece = compute_ece(bench_df.loc[v22_valid, 'truth'], bench_df.loc[v22_valid, 'l2g_v22_max'])
            v22_brier = compute_brier(bench_df.loc[v22_valid, 'truth'], bench_df.loc[v22_valid, 'l2g_v22_max'])
            v22_n = v22_valid.sum()
            
            print(f"\n   L2G v22.09 (gene-level max):")
            print(f"     ECE = {v22_ece:.4f}")
            print(f"     Brier = {v22_brier:.4f}")
            print(f"     n = {v22_n}")
            print(f"\n   RATIO: Our method is {v22_ece/our_full_ece:.1f}x better calibrated than L2G v22.09")
    
    # Compute metrics on matched subset (APPLES-TO-APPLES for v25.12)
    print("\n" + "=" * 80)
    print("PANEL B: MATCHED SUBSET (n={})".format(len(matched_subset)))
    print("Semantic-aligned: study-locus → study → disease mapping")
    print("This is the ONLY valid comparison for L2G v25.12")
    print("=" * 80)
    
    # Our method on matched subset
    matched_valid = matched_subset['calibrated_prob'].notna()
    our_matched_ece = compute_ece(
        matched_subset.loc[matched_valid, 'truth'], 
        matched_subset.loc[matched_valid, 'calibrated_prob']
    )
    our_matched_brier = compute_brier(
        matched_subset.loc[matched_valid, 'truth'],
        matched_subset.loc[matched_valid, 'calibrated_prob']
    )
    our_matched_n = matched_valid.sum()
    
    print(f"\n   Our Method (on matched subset):")
    print(f"     ECE = {our_matched_ece:.4f}")
    print(f"     Brier = {our_matched_brier:.4f}")
    print(f"     n = {our_matched_n}")
    
    # L2G v25.12 semantic-aligned (max aggregation)
    l2g_max_valid = matched_subset['l2g_max'].notna()
    l2g_max_ece = compute_ece(
        matched_subset.loc[l2g_max_valid, 'truth'],
        matched_subset.loc[l2g_max_valid, 'l2g_max']
    )
    l2g_max_brier = compute_brier(
        matched_subset.loc[l2g_max_valid, 'truth'],
        matched_subset.loc[l2g_max_valid, 'l2g_max']
    )
    l2g_max_n = l2g_max_valid.sum()
    
    print(f"\n   L2G v25.12 (max, semantic-aligned):")
    print(f"     ECE = {l2g_max_ece:.4f}")
    print(f"     Brier = {l2g_max_brier:.4f}")
    print(f"     n = {l2g_max_n}")
    
    # L2G v25.12 semantic-aligned (mean aggregation)
    l2g_mean_valid = matched_subset['l2g_mean'].notna()
    l2g_mean_ece = compute_ece(
        matched_subset.loc[l2g_mean_valid, 'truth'],
        matched_subset.loc[l2g_mean_valid, 'l2g_mean']
    )
    l2g_mean_brier = compute_brier(
        matched_subset.loc[l2g_mean_valid, 'truth'],
        matched_subset.loc[l2g_mean_valid, 'l2g_mean']
    )
    l2g_mean_n = l2g_mean_valid.sum()
    
    print(f"\n   L2G v25.12 (mean, semantic-aligned):")
    print(f"     ECE = {l2g_mean_ece:.4f}")
    print(f"     Brier = {l2g_mean_brier:.4f}")
    print(f"     n = {l2g_mean_n}")
    
    print("\n" + "=" * 80)
    print("COMPARISON RATIOS (on matched subset)")
    print("=" * 80)
    print(f"\n   Our method vs L2G v25.12 (max):  {l2g_max_ece/our_matched_ece:.1f}x better ECE")
    print(f"   Our method vs L2G v25.12 (mean): {l2g_mean_ece/our_matched_ece:.1f}x better ECE")
    
    print("\n" + "=" * 80)
    print("PITFALL DEMONSTRATION: Why gene-level collapse is INVALID")
    print("=" * 80)
    
    # Load the gene-level collapsed results
    if benchmark_file.exists():
        bench_df = pd.read_csv(benchmark_file, sep='\t')
        v25_max_valid = bench_df['l2g_v25_max'].notna()
        if v25_max_valid.any():
            v25_collapsed_ece = compute_ece(
                bench_df.loc[v25_max_valid, 'truth'],
                bench_df.loc[v25_max_valid, 'l2g_v25_max']
            )
            print(f"\n   L2G v25.12 (gene-level collapse, INVALID):")
            print(f"     ECE = {v25_collapsed_ece:.4f}")
            print(f"     n = {v25_max_valid.sum()}")
            print(f"\n   This is INVALID because L2G v25.12 is defined per GWAS credible set.")
            print(f"   Collapsing across study-loci conflates disease-specific signals.")
            print(f"   The correct ECE (semantic-aligned) is {l2g_max_ece:.4f} (max) or {l2g_mean_ece:.4f} (mean).")
            print(f"\n   INFLATION FACTOR: {v25_collapsed_ece/l2g_mean_ece:.1f}x worse due to semantic mismatch")
    
    # Summary table for manuscript
    print("\n" + "=" * 80)
    print("MANUSCRIPT TABLE FORMAT")
    print("=" * 80)
    
    print("\n   PANEL A: Full Benchmark (n=14,016) - Gene-Disease Native Methods")
    print("   " + "-" * 60)
    print(f"   {'Method':<35} {'ECE':>8} {'Brier':>8} {'n':>8}")
    print("   " + "-" * 60)
    print(f"   {'Mechanism Graphs (ours)':<35} {our_full_ece:>8.4f} {our_full_brier:>8.4f} {our_full_n:>8}")
    if benchmark_file.exists():
        print(f"   {'L2G v22.09':<35} {v22_ece:>8.4f} {v22_brier:>8.4f} {v22_n:>8}")
    print("   " + "-" * 60)
    
    print(f"\n   PANEL B: Matched Subset (n={len(matched_subset)}) - Semantic-Aligned L2G v25.12")
    print("   " + "-" * 60)
    print(f"   {'Method':<35} {'ECE':>8} {'Brier':>8} {'n':>8}")
    print("   " + "-" * 60)
    print(f"   {'Mechanism Graphs (ours)':<35} {our_matched_ece:>8.4f} {our_matched_brier:>8.4f} {our_matched_n:>8}")
    print(f"   {'L2G v25.12 (mean aggregation)':<35} {l2g_mean_ece:>8.4f} {l2g_mean_brier:>8.4f} {l2g_mean_n:>8}")
    print(f"   {'L2G v25.12 (max aggregation)':<35} {l2g_max_ece:>8.4f} {l2g_max_brier:>8.4f} {l2g_max_n:>8}")
    print("   " + "-" * 60)
    
    print("\n   PITFALL (show separately, clearly labeled INVALID):")
    print("   " + "-" * 60)
    if benchmark_file.exists():
        print(f"   {'L2G v25.12 (gene-level collapse)':<35} {v25_collapsed_ece:>8.4f} {'N/A':>8} {v25_max_valid.sum():>8}")
        print(f"   {'↳ INVALID: conflates disease contexts':<35}")
    
    # Save results to JSON for manuscript updates
    results = {
        "panel_a_full_benchmark": {
            "n": int(our_full_n),
            "our_method": {"ECE": float(our_full_ece), "Brier": float(our_full_brier)},
            "l2g_v22_09": {"ECE": float(v22_ece), "Brier": float(v22_brier), "n": int(v22_n)} if benchmark_file.exists() else None,
        },
        "panel_b_matched_subset": {
            "n": int(our_matched_n),
            "our_method": {"ECE": float(our_matched_ece), "Brier": float(our_matched_brier)},
            "l2g_v25_12_mean": {"ECE": float(l2g_mean_ece), "Brier": float(l2g_mean_brier), "n": int(l2g_mean_n)},
            "l2g_v25_12_max": {"ECE": float(l2g_max_ece), "Brier": float(l2g_max_brier), "n": int(l2g_max_n)},
        },
        "pitfall_invalid": {
            "l2g_v25_12_gene_collapse": {"ECE": float(v25_collapsed_ece)} if benchmark_file.exists() else None,
        },
        "ratios": {
            "our_vs_v22_09": float(v22_ece/our_full_ece) if benchmark_file.exists() else None,
            "our_vs_v25_12_mean_matched": float(l2g_mean_ece/our_matched_ece),
            "our_vs_v25_12_max_matched": float(l2g_max_ece/our_matched_ece),
        }
    }
    
    with open("proper_calibration_comparison.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 80)
    print("RESULTS SAVED TO: proper_calibration_comparison.json")
    print("=" * 80)

if __name__ == "__main__":
    main()
