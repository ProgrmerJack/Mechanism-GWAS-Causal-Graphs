#!/usr/bin/env python3
"""
Evaluate all baselines including cS2G on RegulatoryBench.
Computes stratified AUC-ROC and AUC-PR by evidence type.
"""

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
CANDIDATES_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_cs2g.tsv"


def compute_metrics(labels, scores, score_name=""):
    """Compute AUC-ROC and AUC-PR for a set of predictions."""
    # Filter out NaN scores
    valid = ~np.isnan(scores)
    labels_valid = labels[valid]
    scores_valid = scores[valid]
    
    if len(labels_valid) == 0 or labels_valid.sum() == 0:
        return {'auc_roc': np.nan, 'auc_pr': np.nan, 'n': 0, 'n_pos': 0}
    
    try:
        auc_roc = roc_auc_score(labels_valid, scores_valid)
        auc_pr = average_precision_score(labels_valid, scores_valid)
    except ValueError:
        auc_roc = np.nan
        auc_pr = np.nan
    
    return {
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
        'n': len(labels_valid),
        'n_pos': labels_valid.sum(),
        'coverage': valid.sum() / len(valid)
    }


def main():
    print("=" * 70)
    print("RegulatoryBench Baseline Evaluation (with cS2G)")
    print("=" * 70)
    
    # Load candidates
    print(f"\nLoading candidates from {CANDIDATES_FILE}...")
    df = pd.read_csv(CANDIDATES_FILE, sep='\t', low_memory=False)
    print(f"  Loaded {len(df):,} candidates")
    print(f"  Positive rate: {df['label'].mean()*100:.2f}%")
    
    # Define methods to evaluate
    methods = {
        'NearestGene': 'score_nearest',
        'Within100kb': 'score_100kb',
        'cS2G': 'score_cs2g'
    }
    
    results = []
    
    # Overall evaluation
    print("\n" + "=" * 70)
    print("OVERALL RESULTS")
    print("=" * 70)
    print(f"\n{'Method':<15} {'AUC-ROC':>10} {'AUC-PR':>10} {'Coverage':>10} {'N':>12} {'N+':>8}")
    print("-" * 70)
    
    for method_name, col_name in methods.items():
        if col_name not in df.columns:
            print(f"{method_name:<15} {'N/A':>10}")
            continue
            
        scores = df[col_name].values
        labels = df['label'].values
        
        metrics = compute_metrics(labels, scores, method_name)
        
        print(f"{method_name:<15} {metrics['auc_roc']:>10.4f} {metrics['auc_pr']:>10.4f} "
              f"{metrics['coverage']*100:>9.1f}% {metrics['n']:>12,} {metrics['n_pos']:>8,}")
        
        results.append({
            'method': method_name,
            'evidence_type': 'All',
            **metrics
        })
    
    # Stratified by evidence type
    print("\n" + "=" * 70)
    print("STRATIFIED BY EVIDENCE TYPE")
    print("=" * 70)
    
    for etype in df['evidence_type'].unique():
        mask = df['evidence_type'] == etype
        df_sub = df[mask]
        
        print(f"\n{etype} (n={len(df_sub):,}, {df_sub['label'].sum():,} positives)")
        print(f"{'Method':<15} {'AUC-ROC':>10} {'AUC-PR':>10} {'Coverage':>10}")
        print("-" * 50)
        
        for method_name, col_name in methods.items():
            if col_name not in df.columns:
                continue
                
            scores = df_sub[col_name].values
            labels = df_sub['label'].values
            
            metrics = compute_metrics(labels, scores, method_name)
            
            print(f"{method_name:<15} {metrics['auc_roc']:>10.4f} {metrics['auc_pr']:>10.4f} "
                  f"{metrics['coverage']*100:>9.1f}%")
            
            results.append({
                'method': method_name,
                'evidence_type': etype,
                **metrics
            })
    
    # Create regime map (by distance)
    print("\n" + "=" * 70)
    print("REGIME MAP: Performance by Distance")
    print("=" * 70)
    
    distance_regimes = [
        ('0-10kb', 0, 10000),
        ('10-100kb', 10000, 100000),
        ('100-500kb', 100000, 500000)
    ]
    
    print(f"\n{'Regime':<12} {'N pairs':>10} {'NearestGene':>12} {'cS2G':>12}")
    print("-" * 50)
    
    for regime_name, d_min, d_max in distance_regimes:
        mask = (df['distance'] >= d_min) & (df['distance'] < d_max) & (df['label'] == 1)
        n_positives = mask.sum()
        
        # Get all candidates in loci that have positives in this regime
        loci_with_pos = df.loc[mask, 'locus_id'].unique()
        mask_regime = df['locus_id'].isin(loci_with_pos)
        df_regime = df[mask_regime]
        
        if len(df_regime) == 0:
            print(f"{regime_name:<12} {n_positives:>10,} {'N/A':>12} {'N/A':>12}")
            continue
        
        nearest_auc = compute_metrics(df_regime['label'].values, df_regime['score_nearest'].values)['auc_roc']
        cs2g_auc = compute_metrics(df_regime['label'].values, df_regime['score_cs2g'].values)['auc_roc'] if 'score_cs2g' in df.columns else np.nan
        
        print(f"{regime_name:<12} {n_positives:>10,} {nearest_auc:>12.3f} {cs2g_auc:>12.3f}")
    
    # Save results
    results_df = pd.DataFrame(results)
    output_file = BASE_DIR / "data" / "processed" / "baselines" / "baseline_evaluation_results.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
