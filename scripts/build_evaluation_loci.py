#!/usr/bin/env python3
"""
Build Evaluation Loci with Positive and Negative Candidates
============================================================

This script creates proper evaluation loci by:
1. Taking each experimental positive (CRISPRi/MPRA)
2. Finding all genes within a window (e.g., 500kb)
3. Labeling positives vs negatives
4. Creating locus-level evaluation sets
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
EXTERNAL_DIR = DATA_DIR / "external"
OUTPUT_DIR = DATA_DIR / "processed" / "baselines"

# Parameters
WINDOW_SIZE = 500_000  # 500kb window around each locus


def load_gene_coordinates() -> pd.DataFrame:
    """Load gene coordinates from genome annotation."""
    gene_file = EXTERNAL_DIR / "crispr_benchmark" / "resources" / "genome_annotations" / "CollapsedGeneBounds.hg38.bed"
    
    if not gene_file.exists():
        raise FileNotFoundError(f"Gene coordinates file not found: {gene_file}")
    
    # Read BED file (first row is header)
    df = pd.read_csv(gene_file, sep='\t', skiprows=1, header=None,
                     names=['chr', 'start', 'end', 'symbol', 'score', 'strand', 'ensembl_id', 'gene_type'])
    
    # Filter to protein-coding genes only
    df = df[df['gene_type'] == 'protein_coding'].copy()
    
    # Calculate TSS
    df['tss'] = df.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1)
    
    print(f"Loaded {len(df)} protein-coding genes")
    return df


def load_benchmark() -> pd.DataFrame:
    """Load the RegulatoryBench v3 benchmark."""
    benchmark_file = OUTPUT_DIR / "regulatorybench_v3.tsv"
    
    if not benchmark_file.exists():
        raise FileNotFoundError(f"Benchmark not found: {benchmark_file}")
    
    df = pd.read_csv(benchmark_file, sep='\t')
    print(f"Loaded benchmark: {len(df)} entries")
    return df


def find_genes_in_window(chrom: str, pos: int, gene_df: pd.DataFrame, 
                         window: int = WINDOW_SIZE) -> pd.DataFrame:
    """Find all genes within a window of a position."""
    # Normalize chromosome format
    chrom_norm = f"chr{chrom}" if not str(chrom).startswith('chr') else str(chrom)
    
    # Filter to chromosome
    chr_genes = gene_df[gene_df['chr'] == chrom_norm].copy()
    
    # Find genes with TSS within window
    chr_genes['distance'] = np.abs(chr_genes['tss'] - pos)
    nearby = chr_genes[chr_genes['distance'] <= window].copy()
    
    return nearby.sort_values('distance')


def create_evaluation_loci(benchmark_df: pd.DataFrame, 
                          gene_df: pd.DataFrame,
                          window: int = WINDOW_SIZE) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create evaluation loci with positive and negative candidates.
    
    Returns:
        loci_df: DataFrame of loci with metadata
        candidates_df: DataFrame of all candidates with labels
    """
    print(f"\n=== Creating Evaluation Loci (window={window/1000:.0f}kb) ===")
    
    loci = []
    all_candidates = []
    
    # Group benchmark by locus
    # For CRISPRi: use locus_id (enhancer element)
    # For MPRA: use (chr, pos)
    
    # Process CRISPRi loci
    crispr_data = benchmark_df[benchmark_df['evidence_type'] == 'CRISPRi']
    crispr_loci = crispr_data.groupby('locus_id')
    
    for locus_id, locus_entries in crispr_loci:
        chrom = locus_entries['chr'].iloc[0]
        pos = int(locus_entries['pos'].iloc[0])
        positive_genes = set(locus_entries['gene_ensembl'].dropna())
        positive_symbols = set(locus_entries['gene_symbol'].dropna())
        
        # Find all candidate genes in window
        nearby_genes = find_genes_in_window(chrom, pos, gene_df, window)
        
        if len(nearby_genes) == 0:
            continue
            
        # Create candidate entries
        for _, gene in nearby_genes.iterrows():
            is_positive = gene['ensembl_id'] in positive_genes or gene['symbol'] in positive_symbols
            
            all_candidates.append({
                'locus_id': locus_id,
                'chr': chrom,
                'locus_pos': pos,
                'gene_ensembl': gene['ensembl_id'],
                'gene_symbol': gene['symbol'],
                'gene_tss': gene['tss'],
                'distance': gene['distance'],
                'label': 1 if is_positive else 0,
                'evidence_type': 'CRISPRi'
            })
        
        loci.append({
            'locus_id': locus_id,
            'chr': chrom,
            'pos': pos,
            'evidence_type': 'CRISPRi',
            'n_positives': len(positive_genes),
            'n_candidates': len(nearby_genes),
            'positive_rate': len(positive_genes) / len(nearby_genes) if len(nearby_genes) > 0 else 0
        })
    
    print(f"  CRISPRi: {len([l for l in loci if l['evidence_type'] == 'CRISPRi'])} loci")
    
    # Process MPRA loci
    mpra_data = benchmark_df[benchmark_df['evidence_type'] == 'MPRA']
    mpra_loci = mpra_data.groupby(['chr', 'pos'])
    
    for (chrom, pos), locus_entries in mpra_loci:
        locus_id = f"mpra_{chrom}_{int(pos)}"
        positive_genes = set(locus_entries['gene_ensembl'].dropna())
        
        # Find all candidate genes in window
        nearby_genes = find_genes_in_window(chrom, int(pos), gene_df, window)
        
        if len(nearby_genes) == 0:
            continue
        
        # Create candidate entries
        for _, gene in nearby_genes.iterrows():
            is_positive = gene['ensembl_id'] in positive_genes
            
            all_candidates.append({
                'locus_id': locus_id,
                'chr': chrom,
                'locus_pos': int(pos),
                'gene_ensembl': gene['ensembl_id'],
                'gene_symbol': gene['symbol'],
                'gene_tss': gene['tss'],
                'distance': gene['distance'],
                'label': 1 if is_positive else 0,
                'evidence_type': 'MPRA'
            })
        
        loci.append({
            'locus_id': locus_id,
            'chr': chrom,
            'pos': int(pos),
            'evidence_type': 'MPRA',
            'n_positives': len(positive_genes),
            'n_candidates': len(nearby_genes),
            'positive_rate': len(positive_genes) / len(nearby_genes) if len(nearby_genes) > 0 else 0
        })
    
    print(f"  MPRA: {len([l for l in loci if l['evidence_type'] == 'MPRA'])} loci")
    
    loci_df = pd.DataFrame(loci)
    candidates_df = pd.DataFrame(all_candidates)
    
    # Summary statistics
    print(f"\n  Total loci: {len(loci_df)}")
    print(f"  Total candidates: {len(candidates_df)}")
    print(f"  Total positives: {candidates_df['label'].sum()}")
    print(f"  Total negatives: {len(candidates_df) - candidates_df['label'].sum()}")
    print(f"  Positive rate: {candidates_df['label'].mean():.3f}")
    
    return loci_df, candidates_df


def run_baseline_evaluation(candidates_df: pd.DataFrame) -> pd.DataFrame:
    """
    Run baseline evaluation methods on the candidates.
    """
    print("\n=== Running Baseline Evaluations ===")
    
    results = []
    
    # Method 1: Nearest Gene (inverse distance)
    print("\n  Method: NearestGene (inverse distance)")
    candidates_df['score_nearest'] = 1.0 / (candidates_df['distance'].abs() + 1)
    
    # Normalize per locus
    for locus_id, group in candidates_df.groupby('locus_id'):
        total = group['score_nearest'].sum()
        if total > 0:
            candidates_df.loc[group.index, 'score_nearest'] = group['score_nearest'] / total
    
    # Calculate metrics
    from sklearn.metrics import roc_auc_score, average_precision_score
    
    y_true = candidates_df['label']
    y_score = candidates_df['score_nearest']
    
    if len(y_true.unique()) > 1:
        auc_roc = roc_auc_score(y_true, y_score)
        auc_pr = average_precision_score(y_true, y_score)
    else:
        auc_roc = 0.5
        auc_pr = y_true.mean()
    
    print(f"    AUC-ROC: {auc_roc:.3f}")
    print(f"    AUC-PR: {auc_pr:.3f}")
    
    results.append({
        'method': 'NearestGene',
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
        'n_candidates': len(candidates_df),
        'n_positives': int(y_true.sum())
    })
    
    # Method 2: Within 100kb
    print("\n  Method: Within100kb (binary)")
    candidates_df['score_100kb'] = (candidates_df['distance'] <= 100_000).astype(float)
    
    y_score = candidates_df['score_100kb']
    if len(y_true.unique()) > 1:
        auc_roc = roc_auc_score(y_true, y_score)
        auc_pr = average_precision_score(y_true, y_score)
    else:
        auc_roc = 0.5
        auc_pr = y_true.mean()
    
    print(f"    AUC-ROC: {auc_roc:.3f}")
    print(f"    AUC-PR: {auc_pr:.3f}")
    
    results.append({
        'method': 'Within100kb',
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
        'n_candidates': len(candidates_df),
        'n_positives': int(y_true.sum())
    })
    
    # Method 3: Precision@1 by locus
    print("\n  Precision@1 by locus:")
    
    for method in ['NearestGene', 'Within100kb']:
        score_col = 'score_nearest' if method == 'NearestGene' else 'score_100kb'
        
        correct = 0
        total = 0
        for locus_id, group in candidates_df.groupby('locus_id'):
            if group['label'].sum() == 0:
                continue  # Skip loci without positives
            total += 1
            # Get top-scored gene
            top_gene = group.loc[group[score_col].idxmax()]
            if top_gene['label'] == 1:
                correct += 1
        
        precision_at_1 = correct / total if total > 0 else 0
        print(f"    {method}: {precision_at_1:.3f} ({correct}/{total} loci)")
        
        # Update results
        for r in results:
            if r['method'] == method:
                r['precision_at_1'] = precision_at_1
    
    return pd.DataFrame(results)


def main():
    print("="*60)
    print("BUILDING EVALUATION LOCI WITH CANDIDATES")
    print("="*60)
    
    # Load data
    gene_df = load_gene_coordinates()
    benchmark_df = load_benchmark()
    
    # Create evaluation loci
    loci_df, candidates_df = create_evaluation_loci(benchmark_df, gene_df)
    
    # Run baseline evaluation
    results_df = run_baseline_evaluation(candidates_df)
    
    # Save outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    loci_file = OUTPUT_DIR / "evaluation_loci.tsv"
    loci_df.to_csv(loci_file, sep='\t', index=False)
    print(f"\nSaved loci: {loci_file}")
    
    candidates_file = OUTPUT_DIR / "evaluation_candidates.tsv"
    candidates_df.to_csv(candidates_file, sep='\t', index=False)
    print(f"Saved candidates: {candidates_file}")
    
    results_file = OUTPUT_DIR / "baseline_results.tsv"
    results_df.to_csv(results_file, sep='\t', index=False)
    print(f"Saved results: {results_file}")
    
    # Print final summary
    print("\n" + "="*60)
    print("EVALUATION LOCI SUMMARY")
    print("="*60)
    print(f"Total loci: {len(loci_df)}")
    print(f"  CRISPRi: {len(loci_df[loci_df['evidence_type'] == 'CRISPRi'])}")
    print(f"  MPRA: {len(loci_df[loci_df['evidence_type'] == 'MPRA'])}")
    print(f"\nTotal candidates: {len(candidates_df)}")
    print(f"  Positives: {candidates_df['label'].sum()}")
    print(f"  Negatives: {len(candidates_df) - candidates_df['label'].sum()}")
    print(f"  Positive rate: {candidates_df['label'].mean():.3%}")
    print("\nBaseline Results:")
    print(results_df.to_string(index=False))
    
    return loci_df, candidates_df, results_df


if __name__ == "__main__":
    loci, candidates, results = main()
