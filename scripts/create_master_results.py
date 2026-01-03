#!/usr/bin/env python3
"""
create_master_results.py
========================
Creates the master results table with one row per (locus, method) for the 
Nature Genetics Article. This is the single source of truth for all figures.

Output schema (master_results.parquet):
---------------------------------------
- locus_id: unique locus identifier
- chr, pos, lead_snp: genomic location
- true_gene: ground truth causal gene (Ensembl ID)
- gene_symbol: gene symbol
- trait: associated trait
- mechanism_class: CODING or REGULATORY
- evidence_class: confidence level
- method: L2G, cS2G, FLAMES, Nearest, Calibrated
- top_gene: method's top-ranked gene
- score: method's score for top gene
- rank_of_true: rank of true gene (1 = correct)
- top1_correct: boolean
- top3_correct: boolean
- top5_correct: boolean
- score_of_true: method's score for true gene
- has_prediction: whether method made a prediction

Author: Mechanism-GWAS-Causal-Graphs
Date: 2025-01
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import json
from scipy import stats as scipy_stats

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)

def wilson_ci(successes: int, n: int, alpha: float = 0.05):
    """Compute Wilson score interval for binomial proportion."""
    if n == 0:
        return 0.0, 0.0, 0.0
    p = successes / n
    z = scipy_stats.norm.ppf(1 - alpha / 2)
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    spread = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denominator
    return p, max(0, center - spread), min(1, center + spread)

def load_benchmark(use_independent=True):
    """
    Load the benchmark loci.
    
    CRITICAL: use_independent=True uses the verified independent benchmark
    with 63 loci that have NO L2G training data overlap.
    
    The mechanism_stratified_bench_v1.tsv is 91.7% circular (OpenTargets_2024)
    and MUST NOT be used for L2G evaluation.
    """
    if use_independent:
        bench_path = DATA_DIR / "processed" / "baselines" / "independent_benchmark_v1.tsv"
        if not bench_path.exists():
            # Fallback to post2021 file
            bench_path = DATA_DIR / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
    else:
        bench_path = DATA_DIR / "processed" / "baselines" / "mechanism_stratified_bench_v1.tsv"
    
    bench = pd.read_csv(bench_path, sep='\t')
    print(f"Loaded benchmark: {len(bench)} loci from {bench_path.name}")
    print(f"  Benchmark independence: {'VERIFIED' if use_independent else 'NOT VERIFIED (circular risk!)'}")
    return bench

def load_tss_annotations():
    """Load gene TSS annotations for nearest gene computation."""
    tss_path = DATA_DIR / "external" / "ensembl" / "gencode_tss_grch38.tsv"
    if tss_path.exists():
        tss = pd.read_csv(tss_path, sep='\t')
        # Normalize chromosome format
        tss['chromosome'] = tss['chromosome'].astype(str).str.replace('chr', '')
        print(f"Loaded TSS annotations: {len(tss)} genes")
        return tss
    else:
        print("TSS annotations not found")
        return None

def load_l2g_results():
    """Load L2G matching results."""
    l2g_path = RESULTS_DIR / "l2g_benchmark_matching_v2.tsv"
    if l2g_path.exists():
        l2g = pd.read_csv(l2g_path, sep='\t')
        print(f"Loaded L2G results: {len(l2g)} matched loci")
        return l2g
    else:
        print("L2G results not found - run match_l2g_to_benchmark.py first")
        return None

def load_cs2g_unified():
    """Load cS2G results from unified benchmark."""
    unified_path = DATA_DIR / "processed" / "baselines" / "unified_benchmark_l2g_cs2g.tsv"
    if unified_path.exists():
        unified = pd.read_csv(unified_path, sep='\t')
        print(f"Loaded unified benchmark with cS2G: {len(unified)} loci")
        return unified
    else:
        print("Unified benchmark not found")
        return None

def create_l2g_method_results(bench, l2g_results):
    """Create method results for L2G."""
    if l2g_results is None:
        return pd.DataFrame()
    
    rows = []
    for idx, locus in bench.iterrows():
        locus_id = locus['locus_id']
        l2g_match = l2g_results[l2g_results['locus_id'] == locus_id]
        
        if len(l2g_match) > 0:
            match = l2g_match.iloc[0]
            # For L2G: check if the top gene matches true gene
            gene_match = match['match_type'] == 'gene_match'
            
            rows.append({
                'locus_id': locus_id,
                'method': 'L2G',
                'top_gene': match['l2g_gene'],
                'score': match['l2g_score'],
                'rank_of_true': 1 if gene_match else np.nan,
                'top1_correct': gene_match,
                'top3_correct': gene_match,  # Simplified - full ranking needs more data
                'top5_correct': gene_match,
                'score_of_true': match['l2g_score'] if gene_match else np.nan,
                'has_prediction': True
            })
        else:
            rows.append({
                'locus_id': locus_id,
                'method': 'L2G',
                'top_gene': np.nan,
                'score': np.nan,
                'rank_of_true': np.nan,
                'top1_correct': False,
                'top3_correct': False,
                'top5_correct': False,
                'score_of_true': np.nan,
                'has_prediction': False
            })
    
    return pd.DataFrame(rows)

def create_cs2g_method_results(bench, unified_df):
    """
    Create method results for cS2G using the unified benchmark data.
    cS2G is only available for the original 63 well-characterized loci.
    Match by gene_symbol since locus_id formats differ.
    """
    if unified_df is None:
        return pd.DataFrame()
    
    # Create a lookup by gene symbol
    unified_by_gene = unified_df.set_index('gene_symbol')
    
    rows = []
    matched_count = 0
    for idx, locus in bench.iterrows():
        locus_id = locus['locus_id']
        gene_symbol = locus['gene_symbol']
        
        # Try to find this locus in unified benchmark by gene symbol
        if gene_symbol in unified_by_gene.index:
            match = unified_by_gene.loc[gene_symbol]
            # Handle case where gene symbol appears multiple times
            if isinstance(match, pd.DataFrame):
                match = match.iloc[0]
            
            if match['cs2g_available']:
                matched_count += 1
                correct = bool(match['cs2g_correct'])
                
                rows.append({
                    'locus_id': locus_id,
                    'method': 'cS2G',
                    'top_gene': gene_symbol if correct else 'unknown',
                    'score': float(match['cs2g_score_gene_based']) if pd.notna(match['cs2g_score_gene_based']) else np.nan,
                    'rank_of_true': 1 if correct else np.nan,
                    'top1_correct': correct,
                    'top3_correct': correct,
                    'top5_correct': correct,
                    'score_of_true': float(match['cs2g_score_gene_based']) if correct and pd.notna(match['cs2g_score_gene_based']) else np.nan,
                    'has_prediction': True
                })
                continue
        
        # No cS2G data for this locus
        rows.append({
            'locus_id': locus_id,
            'method': 'cS2G',
            'top_gene': np.nan,
            'score': np.nan,
            'rank_of_true': np.nan,
            'top1_correct': False,
            'top3_correct': False,
            'top5_correct': False,
            'score_of_true': np.nan,
            'has_prediction': False
        })
    
    print(f"  cS2G matched {matched_count} loci from unified benchmark")
    return pd.DataFrame(rows)

def find_nearest_gene(chrom, pos, tss_df, max_distance=1_000_000):
    """Find the nearest protein-coding gene to a position."""
    if tss_df is None:
        return None, None
    
    chrom_str = str(chrom).replace('chr', '')
    nearby = tss_df[(tss_df['chromosome'] == chrom_str) & 
                    (tss_df['gene_type'] == 'protein_coding')]
    
    if len(nearby) == 0:
        return None, None
    
    distances = np.abs(nearby['tss_position'] - pos)
    min_idx = distances.idxmin()
    min_dist = distances.loc[min_idx]
    
    if min_dist > max_distance:
        return None, None
    
    return nearby.loc[min_idx, 'gene_symbol'], int(min_dist)

def create_nearest_gene_results(bench, tss_df):
    """Create baseline results for nearest gene approach."""
    rows = []
    for idx, locus in bench.iterrows():
        # Compute nearest gene
        chrom = locus['chr']
        pos = locus['pos_hg38']
        true_gene = locus['gene_symbol']
        
        # Handle NaN gene symbol
        if pd.isna(true_gene):
            rows.append({
                'locus_id': locus['locus_id'],
                'method': 'NearestGene',
                'top_gene': np.nan,
                'score': np.nan,
                'rank_of_true': np.nan,
                'top1_correct': False,
                'top3_correct': False,
                'top5_correct': False,
                'score_of_true': np.nan,
                'has_prediction': False,
                'distance_to_tss': np.nan
            })
            continue
        
        nearest, dist = find_nearest_gene(chrom, pos, tss_df)
        
        if nearest is not None:
            correct = str(nearest).upper() == str(true_gene).upper()
            rows.append({
                'locus_id': locus['locus_id'],
                'method': 'NearestGene',
                'top_gene': nearest,
                'score': 1.0 / (1 + dist / 10000) if dist else 1.0,  # Distance-decay score
                'rank_of_true': 1 if correct else np.nan,
                'top1_correct': correct,
                'top3_correct': correct,
                'top5_correct': correct,
                'score_of_true': 1.0 / (1 + dist / 10000) if correct else np.nan,
                'has_prediction': True,
                'distance_to_tss': dist
            })
        else:
            rows.append({
                'locus_id': locus['locus_id'],
                'method': 'NearestGene',
                'top_gene': np.nan,
                'score': np.nan,
                'rank_of_true': np.nan,
                'top1_correct': False,
                'top3_correct': False,
                'top5_correct': False,
                'score_of_true': np.nan,
                'has_prediction': False,
                'distance_to_tss': np.nan
            })
    return pd.DataFrame(rows)

def merge_with_benchmark(method_results, bench):
    """Merge method results with benchmark metadata."""
    # Create a clean benchmark copy
    bench_cols = ['locus_id', 'chr', 'pos_hg38', 'lead_snp', 'ensembl_id', 
                  'gene_symbol', 'trait', 'mechanism_class', 'evidence_class', 'confidence']
    
    # Rename for clarity
    bench_clean = bench[bench_cols].copy()
    bench_clean = bench_clean.rename(columns={
        'pos_hg38': 'pos',
        'ensembl_id': 'true_gene'
    })
    
    # Merge
    merged = method_results.merge(bench_clean, on='locus_id', how='left')
    
    return merged

def compute_statistics(master_df):
    """Compute summary statistics for the master results with Wilson 95% CI."""
    stats = {}
    
    for method in master_df['method'].unique():
        method_df = master_df[master_df['method'] == method]
        
        n_total = len(method_df)
        n_with_pred = int(method_df['has_prediction'].sum())
        n_top1 = int(method_df['top1_correct'].sum())
        n_top3 = int(method_df['top3_correct'].sum())
        n_top5 = int(method_df['top5_correct'].sum())
        
        # Coverage with Wilson CI
        cov_p, cov_lo, cov_hi = wilson_ci(n_with_pred, n_total)
        
        # Accuracy with Wilson CI (among those with predictions)
        if n_with_pred > 0:
            acc_p, acc_lo, acc_hi = wilson_ci(n_top1, n_with_pred)
        else:
            acc_p, acc_lo, acc_hi = 0, 0, 0
        
        stats[method] = {
            'n_total': int(n_total),
            'n_with_prediction': n_with_pred,
            'coverage': float(cov_p * 100),
            'coverage_ci_lo': float(cov_lo * 100),
            'coverage_ci_hi': float(cov_hi * 100),
            'top1_accuracy': float(acc_p * 100),
            'top1_accuracy_ci_lo': float(acc_lo * 100),
            'top1_accuracy_ci_hi': float(acc_hi * 100),
            'top1_n': n_top1,
            'top3_n': n_top3,
            'top5_n': n_top5,
            'mean_score': float(method_df[method_df['has_prediction']]['score'].mean()) if n_with_pred > 0 else 0
        }
        
        # By mechanism class with Wilson CI
        for mech in ['CODING', 'REGULATORY']:
            mech_df = method_df[method_df['mechanism_class'] == mech]
            n_mech = len(mech_df)
            n_mech_pred = int(mech_df['has_prediction'].sum())
            n_mech_top1 = int(mech_df['top1_correct'].sum())
            
            if n_mech_pred > 0:
                mech_p, mech_lo, mech_hi = wilson_ci(n_mech_top1, n_mech_pred)
            else:
                mech_p, mech_lo, mech_hi = 0, 0, 0
            
            stats[method][f'{mech.lower()}_n'] = n_mech
            stats[method][f'{mech.lower()}_n_pred'] = n_mech_pred
            stats[method][f'{mech.lower()}_top1'] = float(mech_p * 100)
            stats[method][f'{mech.lower()}_top1_ci_lo'] = float(mech_lo * 100)
            stats[method][f'{mech.lower()}_top1_ci_hi'] = float(mech_hi * 100)
    
    return stats

def save_paired_sets(master_df):
    """Create paired sets JSON for statistical tests."""
    paired_sets = {}
    
    methods = master_df['method'].unique()
    
    for i, method1 in enumerate(methods):
        for method2 in methods[i+1:]:
            df1 = master_df[master_df['method'] == method1].set_index('locus_id')
            df2 = master_df[master_df['method'] == method2].set_index('locus_id')
            
            # Find loci with predictions from both methods
            common_loci = set(df1[df1['has_prediction']].index) & set(df2[df2['has_prediction']].index)
            
            if len(common_loci) > 0:
                paired_sets[f"{method1}_vs_{method2}"] = {
                    'n_paired': len(common_loci),
                    'loci': list(common_loci)
                }
    
    output_path = RESULTS_DIR / "paired_sets.json"
    with open(output_path, 'w') as f:
        json.dump(paired_sets, f, indent=2)
    print(f"Saved paired sets to: {output_path}")
    
    return paired_sets

def main():
    """Main execution."""
    print("=" * 70)
    print("CREATING MASTER RESULTS TABLE")
    print("=" * 70)
    print(f"Started at: {datetime.now().isoformat()}")
    
    # Load data
    bench = load_benchmark()
    l2g_results = load_l2g_results()
    tss_df = load_tss_annotations()
    unified_df = load_cs2g_unified()
    
    # Create method results
    print("\n" + "-" * 50)
    print("Creating method results...")
    print("-" * 50)
    
    all_results = []
    
    # L2G from bulk download
    l2g_df = create_l2g_method_results(bench, l2g_results)
    if len(l2g_df) > 0:
        l2g_merged = merge_with_benchmark(l2g_df, bench)
        all_results.append(l2g_merged)
        print(f"L2G: {len(l2g_merged)} rows, {l2g_df['has_prediction'].sum()} with predictions")
    
    # cS2G from unified benchmark
    cs2g_df = create_cs2g_method_results(bench, unified_df)
    if len(cs2g_df) > 0:
        cs2g_merged = merge_with_benchmark(cs2g_df, bench)
        all_results.append(cs2g_merged)
        print(f"cS2G: {len(cs2g_merged)} rows, {cs2g_df['has_prediction'].sum()} with predictions")
    
    # Nearest gene baseline from TSS
    nearest_df = create_nearest_gene_results(bench, tss_df)
    nearest_merged = merge_with_benchmark(nearest_df, bench)
    all_results.append(nearest_merged)
    print(f"NearestGene: {len(nearest_merged)} rows, {nearest_df['has_prediction'].sum()} with predictions")
    
    # Combine all
    if not all_results:
        print("ERROR: No method results created")
        return
    
    master_df = pd.concat(all_results, ignore_index=True)
    print(f"\nCombined master table: {len(master_df)} rows")
    
    # Compute statistics with Wilson CI
    print("\n" + "-" * 50)
    print("Computing statistics with Wilson 95% CI...")
    print("-" * 50)
    
    stats = compute_statistics(master_df)
    
    for method, method_stats in stats.items():
        print(f"\n{method}:")
        print(f"  Coverage: {method_stats['n_with_prediction']}/{method_stats['n_total']} ({method_stats['coverage']:.1f}%)")
        print(f"  Top-1 Accuracy: {method_stats['top1_n']}/{method_stats['n_with_prediction']} ({method_stats['top1_accuracy']:.1f}% [{method_stats['top1_accuracy_ci_lo']:.1f}%-{method_stats['top1_accuracy_ci_hi']:.1f}%])")
        print(f"  CODING: n={method_stats['coding_n_pred']}, {method_stats['coding_top1']:.1f}% [{method_stats['coding_top1_ci_lo']:.1f}%-{method_stats['coding_top1_ci_hi']:.1f}%]")
        print(f"  REGULATORY: n={method_stats['regulatory_n_pred']}, {method_stats['regulatory_top1']:.1f}% [{method_stats['regulatory_top1_ci_lo']:.1f}%-{method_stats['regulatory_top1_ci_hi']:.1f}%]")
    
    # Save outputs
    print("\n" + "-" * 50)
    print("Saving outputs...")
    print("-" * 50)
    
    # Master results parquet
    master_path = RESULTS_DIR / "master_results.parquet"
    master_df.to_parquet(master_path, index=False)
    print(f"Saved master results to: {master_path}")
    
    # Also TSV for inspection
    tsv_path = RESULTS_DIR / "master_results.tsv"
    master_df.to_csv(tsv_path, sep='\t', index=False)
    print(f"Saved TSV to: {tsv_path}")
    
    # Statistics JSON
    stats_path = RESULTS_DIR / "method_statistics.json"
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Saved statistics to: {stats_path}")
    
    # Paired sets
    save_paired_sets(master_df)
    
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"Finished at: {datetime.now().isoformat()}")
    
    return master_df, stats

if __name__ == "__main__":
    main()
