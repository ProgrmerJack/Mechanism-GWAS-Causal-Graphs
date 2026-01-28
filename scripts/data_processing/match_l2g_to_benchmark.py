#!/usr/bin/env python3
"""
match_l2g_to_benchmark.py
=========================
Improved matching of L2G data to benchmark loci using multiple strategies:
1. Exact position matching (chr:pos)
2. Range matching (within ±500kb of lead variant)
3. Gene-based matching (Ensembl ID)

The key insight: L2G provides scores for study_id + variant + gene combinations.
We need to match our benchmark loci to relevant L2G studies/variants.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

def load_data():
    """Load benchmark and L2G data."""
    # Load benchmark
    bench_path = DATA_DIR / "processed" / "baselines" / "mechanism_stratified_bench_v1.tsv"
    bench = pd.read_csv(bench_path, sep='\t')
    print(f"Loaded benchmark: {len(bench)} loci")
    
    # Load L2G data
    l2g_path = DATA_DIR / "external" / "opentargets_l2g" / "processed" / "l2g_processed.parquet"
    l2g = pd.read_parquet(l2g_path)
    print(f"Loaded L2G: {len(l2g)} entries")
    
    return bench, l2g

def normalize_data(bench, l2g):
    """Normalize chromosome and position formats."""
    # Normalize chromosome format (remove 'chr' prefix if present)
    bench['chr_norm'] = bench['chr'].astype(str).str.replace('chr', '', regex=False)
    l2g['chr_norm'] = l2g['chrom'].astype(str).str.replace('chr', '', regex=False)
    
    # Ensure positions are integers
    bench['pos'] = bench['pos_hg38'].astype(int)
    l2g['pos'] = l2g['pos'].astype(int)
    
    # Extract Ensembl gene ID without version
    bench['ensembl_clean'] = bench['ensembl_id'].str.split('.').str[0]
    l2g['ensembl_clean'] = l2g['ensembl_gene_id'].str.split('.').str[0]
    
    return bench, l2g

def match_by_position(bench, l2g, window_kb=500):
    """Match benchmark loci to L2G by position (within window)."""
    window = window_kb * 1000  # Convert to bp
    
    matches = []
    matched_loci = set()
    
    for idx, locus in bench.iterrows():
        chrom = locus['chr_norm']
        pos = locus['pos']
        gene = locus['ensembl_clean']
        
        # Filter L2G to same chromosome and within window
        l2g_chrom = l2g[l2g['chr_norm'] == chrom]
        l2g_nearby = l2g_chrom[
            (l2g_chrom['pos'] >= pos - window) & 
            (l2g_chrom['pos'] <= pos + window)
        ]
        
        if len(l2g_nearby) > 0:
            matched_loci.add(locus['locus_id'])
            
            # Find L2G entries for the same gene
            l2g_gene = l2g_nearby[l2g_nearby['ensembl_clean'] == gene]
            
            if len(l2g_gene) > 0:
                # Get highest L2G score for this gene
                best_match = l2g_gene.loc[l2g_gene['l2g_score'].idxmax()]
                matches.append({
                    'locus_id': locus['locus_id'],
                    'bench_chr': chrom,
                    'bench_pos': pos,
                    'bench_gene': gene,
                    'bench_mechanism': locus['mechanism_class'],
                    'l2g_variant_id': best_match['variant_id'],
                    'l2g_gene': best_match['ensembl_clean'],
                    'l2g_score': best_match['l2g_score'],
                    'l2g_study': best_match['study_id'],
                    'match_type': 'gene_match',
                    'distance': abs(best_match['pos'] - pos)
                })
            else:
                # No gene match, get best overall score
                best_match = l2g_nearby.loc[l2g_nearby['l2g_score'].idxmax()]
                matches.append({
                    'locus_id': locus['locus_id'],
                    'bench_chr': chrom,
                    'bench_pos': pos,
                    'bench_gene': gene,
                    'bench_mechanism': locus['mechanism_class'],
                    'l2g_variant_id': best_match['variant_id'],
                    'l2g_gene': best_match['ensembl_clean'],
                    'l2g_score': best_match['l2g_score'],
                    'l2g_study': best_match['study_id'],
                    'match_type': 'position_only',
                    'distance': abs(best_match['pos'] - pos)
                })
    
    return pd.DataFrame(matches), matched_loci

def evaluate_l2g_performance(matches_df, bench):
    """Evaluate L2G's ability to identify correct causal gene."""
    results = {
        'total_loci': len(bench),
        'matched_loci': len(matches_df),
        'coverage': len(matches_df) / len(bench) * 100 if len(bench) > 0 else 0
    }
    
    # For loci with gene match, check if L2G correctly prioritized the causal gene
    gene_matches = matches_df[matches_df['match_type'] == 'gene_match']
    results['gene_matches'] = len(gene_matches)
    
    # Check if L2G score > 0.5 (good prediction)
    high_score = gene_matches[gene_matches['l2g_score'] > 0.5]
    results['high_score_gene_matches'] = len(high_score)
    
    # By mechanism class
    if len(matches_df) > 0:
        by_mechanism = matches_df.groupby('bench_mechanism').agg({
            'locus_id': 'count',
            'l2g_score': 'mean'
        }).rename(columns={'locus_id': 'n', 'l2g_score': 'mean_l2g_score'})
        results['by_mechanism'] = by_mechanism.to_dict()
    
    return results

def main():
    print("=" * 70)
    print("L2G TO BENCHMARK MATCHING")
    print("=" * 70)
    
    # Load data
    bench, l2g = load_data()
    bench, l2g = normalize_data(bench, l2g)
    
    print(f"\nBenchmark chromosomes: {sorted(bench['chr_norm'].unique())}")
    print(f"L2G chromosomes (sample): {sorted(l2g['chr_norm'].unique())[:10]}...")
    print(f"L2G unique variants: {l2g['variant_id'].nunique()}")
    print(f"L2G unique genes: {l2g['ensembl_clean'].nunique()}")
    
    # Strategy 1: Position matching with 500kb window
    print("\n" + "-" * 50)
    print("STRATEGY 1: Position matching (±500kb)")
    print("-" * 50)
    
    matches, matched_loci = match_by_position(bench, l2g, window_kb=500)
    
    coverage = len(matched_loci) / len(bench) * 100
    print(f"Matched loci: {len(matched_loci)}/{len(bench)} ({coverage:.1f}%)")
    
    if len(matches) > 0:
        print(f"\nMatches by type:")
        print(matches['match_type'].value_counts())
        
        print(f"\nL2G score distribution for matches:")
        print(matches['l2g_score'].describe())
        
        print(f"\nMatches by mechanism:")
        print(matches.groupby('bench_mechanism')['l2g_score'].agg(['count', 'mean', 'std']))
        
        # Save results
        output_path = RESULTS_DIR / "l2g_benchmark_matching_v2.tsv"
        matches.to_csv(output_path, sep='\t', index=False)
        print(f"\nSaved matches to: {output_path}")
    
    # Evaluate performance
    print("\n" + "-" * 50)
    print("PERFORMANCE EVALUATION")
    print("-" * 50)
    
    results = evaluate_l2g_performance(matches, bench)
    for key, value in results.items():
        if key != 'by_mechanism':
            print(f"{key}: {value}")
    
    # Report unmatched loci
    unmatched = bench[~bench['locus_id'].isin(matched_loci)]
    if len(unmatched) > 0:
        print(f"\nUnmatched loci: {len(unmatched)}")
        print("Sample unmatched:")
        print(unmatched[['locus_id', 'chr_norm', 'pos', 'ensembl_clean', 'mechanism_class']].head(10))
    
    return matches, results

if __name__ == "__main__":
    main()
