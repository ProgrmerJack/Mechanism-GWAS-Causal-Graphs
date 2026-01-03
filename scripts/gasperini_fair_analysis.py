#!/usr/bin/env python3
"""
Comprehensive Gasperini 2019 Analysis with Fair Baselines
"""

import gzip
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats
import json
import os

def bootstrap_auroc(y_true, y_score, n_boot=1000):
    """Bootstrap AUROC with CI"""
    np.random.seed(42)
    scores = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) > 1:
            scores.append(roc_auc_score(y_true[idx], y_score[idx]))
    return np.median(scores), np.percentile(scores, 2.5), np.percentile(scores, 97.5)

def main():
    print("=" * 70)
    print("GASPERINI 2019 - COMPREHENSIVE FAIR BASELINE ANALYSIS")
    print("=" * 70)
    
    # Load Gasperini
    gas_file = 'data/external/crispr_validation/gasperini_2019_results.txt.gz'
    with gzip.open(gas_file, 'rt', errors='replace') as f:
        gas = pd.read_csv(f, sep='\t', low_memory=False)
    
    print(f"\nTotal tested pairs: {len(gas):,}")
    print(f"Columns: {list(gas.columns)}")
    
    # Convert p-value to numeric
    gas['pvalue.empirical.adjusted'] = pd.to_numeric(gas['pvalue.empirical.adjusted'], errors='coerce')
    
    # Define significance
    gas['significant'] = gas['pvalue.empirical.adjusted'] < 0.05
    n_sig = gas['significant'].sum()
    print(f"\nSignificant (padj < 0.05): {n_sig:,}")
    print(f"Non-significant: {(~gas['significant']).sum():,}")
    
    # Parse enhancer coordinates from actual columns
    print("\nParsing enhancer coordinates...")
    
    # Use target_site columns which contain the enhancer coordinates
    if 'target_site.chr' in gas.columns:
        gas['chrom'] = gas['target_site.chr']
        gas['start'] = pd.to_numeric(gas['target_site.start'], errors='coerce')
        gas['end'] = pd.to_numeric(gas['target_site.stop'], errors='coerce')
        print("Using target_site columns for coordinates")
    else:
        # Fallback to parsing pairs4merge
        print(f"Sample pairs4merge: {gas['pairs4merge'].iloc[:3].tolist()}")
        
        def parse_coords(s):
            try:
                parts = str(s).split('|')
                if len(parts) >= 1:
                    coord = parts[0]
                    if ':' in coord and '-' in coord:
                        chrom = coord.split(':')[0]
                        start_end = coord.split(':')[1]
                        start = int(start_end.split('-')[0])
                        end = int(start_end.split('-')[1])
                        return chrom, start, end
            except:
                pass
            return None, None, None
        
        coords = gas['pairs4merge'].apply(parse_coords)
        gas['chrom'] = [c[0] for c in coords]
        gas['start'] = [c[1] for c in coords]
        gas['end'] = [c[2] for c in coords]
    
    # Drop rows without coordinates
    gas = gas.dropna(subset=['chrom', 'start', 'end'])
    print(f"Pairs with valid coordinates: {len(gas):,}")
    
    # Load ABC predictions
    print("\nLoading ABC predictions...")
    abc_file = 'data/external/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'
    with gzip.open(abc_file, 'rt') as f:
        abc = pd.read_csv(f, sep='\t')
    
    print(f"Total ABC predictions: {len(abc):,}")
    
    # Filter to K562 (matches Gasperini cell type)
    k562_types = [c for c in abc['CellType'].unique() if 'K562' in c]
    abc_k562 = abc[abc['CellType'].isin(k562_types)].copy()
    print(f"K562 ABC predictions: {len(abc_k562):,}")
    
    # Also get blood-related cell types for broader coverage
    blood_types = ['K562-Roadmap', 'erythroblast-Corces2016', 'GM12878-Roadmap']
    abc_blood = abc[abc['CellType'].isin(blood_types)].copy()
    print(f"Blood cell types ABC predictions: {len(abc_blood):,}")
    
    # Match Gasperini pairs to ABC predictions
    print("\nMatching Gasperini pairs to ABC predictions...")
    
    matched_data = []
    for idx, row in gas.iterrows():
        if idx % 50000 == 0:
            print(f"  Processing {idx:,}/{len(gas):,}...")
        
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        gene = row['gene_short_name']
        
        # Find overlapping ABC predictions for this gene
        overlap = abc_blood[
            (abc_blood['chr'] == chrom) &
            (abc_blood['start'] < end) &
            (abc_blood['end'] > start) &
            (abc_blood['TargetGene'] == gene)
        ]
        
        abc_score = overlap['ABC.Score'].max() if len(overlap) > 0 else 0.0
        
        matched_data.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'gene': gene,
            'significant': row['significant'],
            'abc_score': abc_score if pd.notna(abc_score) else 0.0
        })
    
    matched = pd.DataFrame(matched_data)
    print(f"\nMatched pairs: {len(matched):,}")
    print(f"Pairs with ABC > 0: {(matched['abc_score'] > 0).sum():,}")
    print(f"Significant with ABC > 0: {matched[matched['abc_score'] > 0]['significant'].sum():,}")
    
    # Subset to pairs with ABC coverage for fair comparison
    covered = matched[matched['abc_score'] > 0].copy()
    print(f"\n=== PAIRS WITH ABC COVERAGE (n={len(covered):,}) ===")
    print(f"Significant: {covered['significant'].sum():,}")
    print(f"Non-significant: {(~covered['significant']).sum():,}")
    
    if len(covered) > 1 and covered['significant'].sum() > 0 and (~covered['significant']).sum() > 0:
        y_true = covered['significant'].values.astype(int)
        y_abc = covered['abc_score'].values
        
        auroc = roc_auc_score(y_true, y_abc)
        auroc_med, auroc_lo, auroc_hi = bootstrap_auroc(y_true, y_abc)
        
        print(f"\nPair-level ABC AUROC: {auroc:.4f} [95% CI: {auroc_lo:.4f} - {auroc_hi:.4f}]")
    
    # Distance analysis
    print("\n" + "=" * 70)
    print("DISTANCE BASELINE ANALYSIS")
    print("=" * 70)
    
    # Get gene TSS information
    # Calculate approximate distance to gene
    # For now, use enhancer midpoint as proxy
    
    print("\nSummary saved.")
    
    # Save results
    results = {
        "dataset": "Gasperini_2019",
        "n_total": len(gas),
        "n_significant": int(gas['significant'].sum()),
        "n_with_abc": int((matched['abc_score'] > 0).sum()),
        "n_sig_with_abc": int(matched[matched['abc_score'] > 0]['significant'].sum())
    }
    
    if len(covered) > 1:
        results["auroc_abc_pair_level"] = float(auroc)
        results["auroc_ci"] = [float(auroc_lo), float(auroc_hi)]
    
    os.makedirs('data/processed/prospective_validation', exist_ok=True)
    with open('data/processed/prospective_validation/gasperini_fair_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
