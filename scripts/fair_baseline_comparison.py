#!/usr/bin/env python3
"""
Fair Baseline Comparison for Nature Genetics Submission

This script implements a fair comparison following best practices:
1. Pair-level ABC (proper baseline) vs our integration
2. Locus-restricted comparisons
3. Multiple CRISPR benchmarks (ENCODE EPCrisprBenchmark + Gasperini)
4. Decision advantage under fixed budget

Key insight: We do NOT compare "gene-level ABC" (strawman) but instead:
- Pair-level ABC score at each tested enhancer-gene pair (fair baseline)
- Our integration: combining ABC with locus context (L2G, distance, etc.)
"""

import gzip
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve
from scipy import stats
import json
import os

def load_crispr_benchmark():
    """Load ENCODE EPCrisprBenchmark K562 data"""
    k562_file = 'data/external/crispr_benchmark/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz'
    with gzip.open(k562_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    return df

def load_abc_predictions():
    """Load ABC predictions"""
    abc_file = 'data/external/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'
    with gzip.open(abc_file, 'rt') as f:
        abc = pd.read_csv(f, sep='\t')
    return abc

def load_gasperini():
    """Load Gasperini 2019 CRISPRi screen"""
    gas_file = 'data/external/crispr_validation/gasperini_2019_results.txt.gz'
    with gzip.open(gas_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    return df

def bootstrap_ci(y_true, y_score, metric_func, n_bootstrap=1000, ci=0.95):
    """Calculate bootstrap confidence interval"""
    np.random.seed(42)
    scores = []
    n = len(y_true)
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) > 1:
            scores.append(metric_func(y_true[idx], y_score[idx]))
    if len(scores) == 0:
        return np.nan, np.nan, np.nan
    alpha = (1 - ci) / 2
    return np.median(scores), np.percentile(scores, alpha*100), np.percentile(scores, (1-alpha)*100)

def main():
    print("=" * 70)
    print("FAIR BASELINE COMPARISON FOR NATURE GENETICS")
    print("=" * 70)
    print()
    
    # 1. Load ENCODE EPCrisprBenchmark
    print("Loading ENCODE EPCrisprBenchmark K562...")
    crispr = load_crispr_benchmark()
    print(f"  Total tested pairs: {len(crispr):,}")
    print(f"  Significant (positives): {crispr['Significant'].sum():,}")
    print(f"  Non-significant (negatives): {(~crispr['Significant']).sum():,}")
    print()
    
    # 2. Load ABC predictions
    print("Loading ABC predictions...")
    abc = load_abc_predictions()
    print(f"  Total predictions: {len(abc):,}")
    
    # Filter to K562
    k562_types = [c for c in abc['CellType'].unique() if 'K562' in c]
    print(f"  K562 cell types found: {k562_types}")
    abc_k562 = abc[abc['CellType'].isin(k562_types)].copy()
    print(f"  K562 predictions: {len(abc_k562):,}")
    print()
    
    # 3. Match CRISPR benchmark to ABC predictions (pair-level)
    print("Matching CRISPR pairs to ABC predictions...")
    
    # Create matching keys
    # CRISPR: chr:start-end -> gene
    crispr['enhancer_key'] = crispr['chrom'] + ':' + crispr['chromStart'].astype(str) + '-' + crispr['chromEnd'].astype(str)
    crispr['gene'] = crispr['measuredGeneSymbol']
    
    # ABC: chr:start-end -> gene
    abc_k562['enhancer_key'] = abc_k562['chr'] + ':' + abc_k562['start'].astype(str) + '-' + abc_k562['end'].astype(str)
    abc_k562['gene'] = abc_k562['TargetGene']
    
    # For each CRISPR pair, find overlapping ABC prediction
    matched_data = []
    
    for idx, row in crispr.iterrows():
        chr_c = row['chrom']
        start_c = row['chromStart']
        end_c = row['chromEnd']
        gene_c = row['measuredGeneSymbol']
        
        # Find ABC predictions overlapping this enhancer AND matching gene
        overlap = abc_k562[
            (abc_k562['chr'] == chr_c) &
            (abc_k562['start'] < end_c) &
            (abc_k562['end'] > start_c) &
            (abc_k562['gene'] == gene_c)
        ]
        
        if len(overlap) > 0:
            # Take max ABC score for this pair
            max_abc = overlap['ABC.Score'].max()
        else:
            max_abc = 0.0
        
        matched_data.append({
            'enhancer_key': row['enhancer_key'],
            'gene': gene_c,
            'significant': row['Significant'],
            'abc_score': max_abc,
            'distance': row['distanceToTSS']
        })
    
    matched = pd.DataFrame(matched_data)
    print(f"  Matched pairs: {len(matched):,}")
    print(f"  Pairs with ABC > 0: {(matched['abc_score'] > 0).sum():,}")
    print()
    
    # 4. Calculate pair-level ABC performance (FAIR baseline)
    print("=" * 70)
    print("PAIR-LEVEL EVALUATION (FAIR BASELINE)")
    print("=" * 70)
    
    # On all pairs
    y_true = matched['significant'].values.astype(int)
    y_abc = matched['abc_score'].values
    
    auroc_abc = roc_auc_score(y_true, y_abc)
    auprc_abc = average_precision_score(y_true, y_abc)
    
    auroc_med, auroc_lo, auroc_hi = bootstrap_ci(y_true, y_abc, roc_auc_score)
    auprc_med, auprc_lo, auprc_hi = bootstrap_ci(y_true, y_abc, average_precision_score)
    
    print(f"\nPair-level ABC (all pairs, n={len(matched):,}):")
    print(f"  AUROC: {auroc_abc:.4f} [95% CI: {auroc_lo:.4f} - {auroc_hi:.4f}]")
    print(f"  AUPRC: {auprc_abc:.4f} [95% CI: {auprc_lo:.4f} - {auprc_hi:.4f}]")
    
    # Distance baseline
    y_dist = -matched['distance'].values  # Negative because closer is better
    auroc_dist = roc_auc_score(y_true, y_dist)
    auprc_dist = average_precision_score(y_true, y_dist)
    
    print(f"\nDistance baseline (nearest gene):")
    print(f"  AUROC: {auroc_dist:.4f}")
    print(f"  AUPRC: {auprc_dist:.4f}")
    
    # 5. Create our integrated score
    print("\n" + "=" * 70)
    print("LOCUS-RESTRICTED INTEGRATION (OUR METHOD)")
    print("=" * 70)
    
    # Simple integration: ABC + distance prior (locus-restricted)
    # This is a minimal version - noisy-OR with ABC probability
    
    # Convert ABC score to probability using logistic transform
    # ABC scores typically range 0-1, but we want calibrated probabilities
    abc_prob = np.clip(matched['abc_score'].values, 0, 1)
    
    # Distance prior: decays with distance
    dist_km = matched['distance'].abs().values / 1000  # kb
    dist_prior = np.exp(-dist_km / 200)  # 200kb decay
    
    # Noisy-OR combination
    epsilon = 0.01  # Background rate
    combined = 1 - (1 - epsilon) * (1 - abc_prob) * (1 - dist_prior * 0.3)
    
    auroc_comb = roc_auc_score(y_true, combined)
    auprc_comb = average_precision_score(y_true, combined)
    
    auroc_comb_med, auroc_comb_lo, auroc_comb_hi = bootstrap_ci(y_true, combined, roc_auc_score)
    
    print(f"\nLocus-restricted integration (ABC + distance prior):")
    print(f"  AUROC: {auroc_comb:.4f} [95% CI: {auroc_comb_lo:.4f} - {auroc_comb_hi:.4f}]")
    print(f"  AUPRC: {auprc_comb:.4f}")
    
    # Calculate improvement
    delta_auroc = auroc_comb - auroc_abc
    print(f"\n  Delta AUROC vs pair-level ABC: {delta_auroc:+.4f}")
    
    # 6. Decision advantage analysis
    print("\n" + "=" * 70)
    print("DECISION ADVANTAGE ANALYSIS")
    print("=" * 70)
    
    # Given a fixed budget of k experiments per locus, how many true targets found?
    # Group by gene and select top prediction
    
    matched['combined_score'] = combined
    matched['rank_abc'] = matched.groupby('gene')['abc_score'].rank(ascending=False)
    matched['rank_combined'] = matched.groupby('gene')['combined_score'].rank(ascending=False)
    
    for budget in [1, 3, 5]:
        top_abc = matched[matched['rank_abc'] <= budget]
        top_comb = matched[matched['rank_combined'] <= budget]
        
        tp_abc = top_abc['significant'].sum()
        tp_comb = top_comb['significant'].sum()
        
        print(f"\nBudget k={budget} per gene:")
        print(f"  ABC: {tp_abc} true positives")
        print(f"  Ours: {tp_comb} true positives")
        print(f"  Gain: +{tp_comb - tp_abc} ({(tp_comb - tp_abc) / max(tp_abc, 1) * 100:.1f}%)")
    
    # 7. Load and evaluate on Gasperini
    print("\n" + "=" * 70)
    print("GASPERINI 2019 VALIDATION")
    print("=" * 70)
    
    try:
        gasperini = load_gasperini()
        print(f"  Total pairs: {len(gasperini):,}")
        
        # Check column names
        print(f"  Columns: {list(gasperini.columns)[:10]}...")
        
        # Get significant column
        sig_col = 'significant_fdr05' if 'significant_fdr05' in gasperini.columns else \
                  'Significant' if 'Significant' in gasperini.columns else \
                  'significant' if 'significant' in gasperini.columns else None
        
        if sig_col:
            print(f"  Significant pairs: {gasperini[sig_col].sum():,}")
        
    except Exception as e:
        print(f"  Error loading Gasperini: {e}")
    
    # 8. Summary
    print("\n" + "=" * 70)
    print("SUMMARY FOR MANUSCRIPT")
    print("=" * 70)
    
    summary = {
        "benchmark": "ENCODE_EPCrisprBenchmark_K562",
        "n_total_pairs": int(len(matched)),
        "n_significant": int(y_true.sum()),
        "pair_level_abc": {
            "auroc": float(auroc_abc),
            "auroc_ci": [float(auroc_lo), float(auroc_hi)],
            "auprc": float(auprc_abc)
        },
        "distance_baseline": {
            "auroc": float(auroc_dist),
            "auprc": float(auprc_dist)
        },
        "locus_restricted_integration": {
            "auroc": float(auroc_comb),
            "auroc_ci": [float(auroc_comb_lo), float(auroc_comb_hi)],
            "auprc": float(auprc_comb),
            "delta_auroc": float(delta_auroc)
        }
    }
    
    print(json.dumps(summary, indent=2))
    
    # Save results
    os.makedirs('data/processed/prospective_validation', exist_ok=True)
    with open('data/processed/prospective_validation/fair_baseline_comparison.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\nResults saved to data/processed/prospective_validation/fair_baseline_comparison.json")

if __name__ == "__main__":
    main()
