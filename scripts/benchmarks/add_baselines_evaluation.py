#!/usr/bin/env python3
"""
Add L2G and Additional Baselines to Evaluation
===============================================

This script adds Open Targets L2G scores and other baseline methods
to the evaluation on RegulatoryBench v3.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
EXTERNAL_DIR = DATA_DIR / "external"
OUTPUT_DIR = DATA_DIR / "processed" / "baselines"


def load_candidates() -> pd.DataFrame:
    """Load evaluation candidates."""
    candidates_file = OUTPUT_DIR / "evaluation_candidates.tsv"
    df = pd.read_csv(candidates_file, sep='\t')
    print(f"Loaded {len(df)} candidates from {len(df['locus_id'].unique())} loci")
    return df


def load_l2g_scores() -> pd.DataFrame:
    """Load L2G scores from Open Targets."""
    l2g_file = EXTERNAL_DIR / "opentargets_l2g" / "processed" / "l2g_processed.parquet"
    
    if not l2g_file.exists():
        print(f"L2G file not found: {l2g_file}")
        return pd.DataFrame()
    
    df = pd.read_parquet(l2g_file)
    print(f"Loaded {len(df)} L2G predictions")
    
    # Create lookup key: chr_pos_gene
    df['chr'] = df['chrom'].astype(str)
    df['pos'] = df['pos'].astype(int)
    df['gene'] = df['ensembl_gene_id']
    
    return df


def load_cs2g_scores() -> pd.DataFrame:
    """Load cS2G scores if available."""
    cs2g_files = [
        EXTERNAL_DIR / "cs2g" / "cs2g_scores.tsv",
        EXTERNAL_DIR / "cs2g" / "cs2g_predictions.tsv",
        DATA_DIR / "raw" / "cs2g_scores.tsv",
    ]
    
    for f in cs2g_files:
        if f.exists():
            df = pd.read_csv(f, sep='\t')
            print(f"Loaded cS2G scores from {f}")
            return df
    
    print("cS2G scores not found")
    return pd.DataFrame()


def add_l2g_scores(candidates_df: pd.DataFrame, l2g_df: pd.DataFrame) -> pd.DataFrame:
    """Add L2G scores to candidates."""
    print("\n=== Adding L2G Scores ===")
    
    if l2g_df.empty:
        candidates_df['score_l2g'] = np.nan
        return candidates_df
    
    # Create lookup dict: (chr, pos) -> {gene: score}
    l2g_lookup = {}
    for _, row in l2g_df.iterrows():
        key = (str(row['chr']).replace('chr', ''), int(row['pos']))
        if key not in l2g_lookup:
            l2g_lookup[key] = {}
        l2g_lookup[key][row['gene']] = row['l2g_score']
    
    print(f"  Created lookup with {len(l2g_lookup)} positions")
    
    # Match candidates to L2G scores
    scores = []
    matched = 0
    
    for _, cand in candidates_df.iterrows():
        chr_norm = str(cand['chr']).replace('chr', '')
        pos = int(cand['locus_pos'])
        gene = cand['gene_ensembl']
        
        key = (chr_norm, pos)
        if key in l2g_lookup and gene in l2g_lookup[key]:
            scores.append(l2g_lookup[key][gene])
            matched += 1
        else:
            # Try nearby positions (±1bp for rounding)
            found = False
            for delta in [-1, 1, -2, 2]:
                alt_key = (chr_norm, pos + delta)
                if alt_key in l2g_lookup and gene in l2g_lookup[alt_key]:
                    scores.append(l2g_lookup[alt_key][gene])
                    matched += 1
                    found = True
                    break
            if not found:
                scores.append(np.nan)
    
    candidates_df['score_l2g'] = scores
    print(f"  Matched {matched} candidates ({matched/len(candidates_df)*100:.1f}%)")
    print(f"  Missing: {candidates_df['score_l2g'].isna().sum()}")
    
    return candidates_df


def evaluate_method(candidates_df: pd.DataFrame, score_col: str, 
                   method_name: str) -> dict:
    """Evaluate a single method."""
    # Filter to candidates with scores
    valid = candidates_df[candidates_df[score_col].notna()].copy()
    
    if len(valid) == 0:
        return {'method': method_name, 'auc_roc': np.nan, 'auc_pr': np.nan,
                'precision_at_1': np.nan, 'n_candidates': 0, 'n_positives': 0}
    
    y_true = valid['label']
    y_score = valid[score_col]
    
    # Handle edge cases
    if len(y_true.unique()) < 2:
        auc_roc = 0.5
        auc_pr = y_true.mean()
    else:
        auc_roc = roc_auc_score(y_true, y_score)
        auc_pr = average_precision_score(y_true, y_score)
    
    # Precision@1 by locus
    correct = 0
    total = 0
    for locus_id, group in valid.groupby('locus_id'):
        if group['label'].sum() == 0:
            continue
        total += 1
        top_gene = group.loc[group[score_col].idxmax()]
        if top_gene['label'] == 1:
            correct += 1
    
    precision_at_1 = correct / total if total > 0 else np.nan
    
    return {
        'method': method_name,
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
        'precision_at_1': precision_at_1,
        'n_candidates': len(valid),
        'n_positives': int(y_true.sum()),
        'n_loci': total
    }


def evaluate_by_evidence_type(candidates_df: pd.DataFrame, score_col: str,
                              method_name: str) -> dict:
    """Evaluate stratified by evidence type."""
    results = {}
    
    for evidence_type in ['CRISPRi', 'MPRA']:
        subset = candidates_df[candidates_df['evidence_type'] == evidence_type]
        if len(subset) == 0:
            continue
        
        result = evaluate_method(subset, score_col, f"{method_name}_{evidence_type}")
        results[evidence_type] = result
    
    return results


def main():
    print("="*60)
    print("ADDING ADDITIONAL BASELINES TO EVALUATION")
    print("="*60)
    
    # Load candidates
    candidates_df = load_candidates()
    
    # Load baseline scores
    l2g_df = load_l2g_scores()
    cs2g_df = load_cs2g_scores()
    
    # Add L2G scores
    candidates_df = add_l2g_scores(candidates_df, l2g_df)
    
    # Evaluate all methods
    print("\n=== Evaluating Methods ===")
    
    methods = [
        ('score_nearest', 'NearestGene'),
        ('score_l2g', 'L2G'),
    ]
    
    results = []
    
    for score_col, method_name in methods:
        if score_col not in candidates_df.columns:
            continue
        
        # Add score column if missing
        if score_col == 'score_nearest' and score_col not in candidates_df.columns:
            candidates_df['score_nearest'] = 1.0 / (candidates_df['distance'].abs() + 1)
        
        print(f"\n  {method_name}:")
        result = evaluate_method(candidates_df, score_col, method_name)
        results.append(result)
        
        print(f"    AUC-ROC: {result['auc_roc']:.3f}")
        print(f"    AUC-PR: {result['auc_pr']:.3f}")
        print(f"    Precision@1: {result['precision_at_1']:.3f}")
        print(f"    Coverage: {result['n_candidates']}/{len(candidates_df)} candidates")
        
        # Stratified results
        strat_results = evaluate_by_evidence_type(candidates_df, score_col, method_name)
        for evidence_type, strat_result in strat_results.items():
            print(f"    {evidence_type}: AUC-ROC={strat_result['auc_roc']:.3f}, "
                  f"P@1={strat_result['precision_at_1']:.3f}")
    
    # Save results
    results_df = pd.DataFrame(results)
    results_file = OUTPUT_DIR / "baseline_results_extended.tsv"
    results_df.to_csv(results_file, sep='\t', index=False)
    print(f"\nSaved: {results_file}")
    
    # Save updated candidates
    candidates_file = OUTPUT_DIR / "evaluation_candidates_with_scores.tsv"
    candidates_df.to_csv(candidates_file, sep='\t', index=False)
    print(f"Saved: {candidates_file}")
    
    # Print comparison table
    print("\n" + "="*60)
    print("BASELINE COMPARISON")
    print("="*60)
    print(results_df.to_string(index=False))
    
    return results_df


if __name__ == "__main__":
    results = main()
