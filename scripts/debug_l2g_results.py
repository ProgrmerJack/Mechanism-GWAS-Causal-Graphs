#!/usr/bin/env python3
"""Debug L2G results - understand why 38.7% accuracy."""

import pandas as pd
import json

def main():
    # Load master results
    master = pd.read_csv('results/master_results_independent.tsv', sep='\t')
    l2g = master[master['method'] == 'L2G'].copy()
    nearest = master[master['method'] == 'NearestGene'].copy()
    
    print("L2G WRONG CASES:")
    print("="*70)
    wrong = l2g[~l2g['top1_correct']]
    for idx, row in wrong.iterrows():
        print(f"  {row['gene_symbol']}: L2G predicted '{row['top_gene']}' (score={row['score']:.3f})")
    
    print(f"\nTotal wrong: {len(wrong)}/62")
    
    print("\n\nL2G CORRECT CASES:")
    print("="*70)
    correct = l2g[l2g['top1_correct']]
    for idx, row in correct.iterrows():
        print(f"  {row['gene_symbol']}: correct (score={row['score']:.3f})")
    
    print(f"\nTotal correct: {len(correct)}/62")
    
    # Check overlap - where does L2G disagree with NearestGene?
    print("\n\nCOMPARISON: L2G vs NearestGene")
    print("="*70)
    
    l2g_set = set(l2g[l2g['top1_correct']]['gene_symbol'])
    nearest_set = set(nearest[nearest['top1_correct']]['gene_symbol'])
    
    both_correct = l2g_set & nearest_set
    l2g_only = l2g_set - nearest_set
    nearest_only = nearest_set - l2g_set
    both_wrong = set(l2g['gene_symbol']) - l2g_set - nearest_set
    
    print(f"Both correct: {len(both_correct)}")
    print(f"L2G only correct: {len(l2g_only)} -> {l2g_only}")
    print(f"NearestGene only correct: {len(nearest_only)} -> {nearest_only}")
    print(f"Both wrong: {len(both_wrong)} -> {both_wrong}")
    
    # Load benchmark to check mechanism types
    bench = pd.read_csv('data/processed/baselines/independent_benchmark_v1.tsv', sep='\t')
    
    print("\n\nACCURACY BY MECHANISM TYPE:")
    print("="*70)
    
    for mech in ['REGULATORY', 'CODING', 'AMBIGUOUS']:
        mech_genes = set(bench[bench['mechanism'] == mech]['gene_symbol'])
        l2g_correct_mech = mech_genes & l2g_set
        nearest_correct_mech = mech_genes & nearest_set
        
        total = len(mech_genes)
        if total > 0:
            l2g_pct = 100 * len(l2g_correct_mech) / total
            near_pct = 100 * len(nearest_correct_mech) / total
            print(f"{mech}: L2G {len(l2g_correct_mech)}/{total} ({l2g_pct:.1f}%), NearestGene {len(nearest_correct_mech)}/{total} ({near_pct:.1f}%)")

if __name__ == "__main__":
    main()
