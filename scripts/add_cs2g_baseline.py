#!/usr/bin/env python3
"""
Add cS2G baseline scores to the evaluation candidates.

This script loads precomputed cS2G scores (Gazal et al. 2022, Nat Genet)
and matches them to our RegulatoryBench evaluation candidates.

cS2G scores are SNP-gene pairs with scores from 0-1 representing the
combined probability that a SNP regulates a gene.
"""

import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent
CS2G_DIR = BASE_DIR / "data" / "external" / "cs2g" / "cS2G_1000GEUR"
CANDIDATES_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates.tsv"
OUTPUT_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_cs2g.tsv"


def load_cs2g_scores():
    """Load all cS2G scores from chromosome files."""
    print("Loading cS2G scores...")
    
    all_scores = []
    
    for chrom in range(1, 23):
        filepath = CS2G_DIR / f"cS2G.{chrom}.SGscore.gz"
        
        if not filepath.exists():
            print(f"  Warning: {filepath} not found, skipping")
            continue
            
        print(f"  Loading chromosome {chrom}...")
        
        # Read the gzipped file
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
        df['chrom'] = chrom
        all_scores.append(df)
    
    if not all_scores:
        raise ValueError("No cS2G score files found!")
    
    scores = pd.concat(all_scores, ignore_index=True)
    print(f"  Loaded {len(scores):,} SNP-gene pairs")
    
    return scores


def parse_snp_position(snp_id):
    """
    Parse SNP ID to extract chromosome and position.
    cS2G uses format like 'rs12345' or 'chr1:12345:A:G'
    """
    if snp_id.startswith('rs'):
        # Will need to use a SNP lookup table - skip for now
        return None, None
    elif ':' in snp_id:
        parts = snp_id.split(':')
        chrom = parts[0].replace('chr', '')
        try:
            pos = int(parts[1])
            return chrom, pos
        except ValueError:
            return None, None
    return None, None


def match_scores_to_candidates(candidates, cs2g_scores):
    """
    Match cS2G scores to evaluation candidates.
    
    cS2G scores are keyed by (SNP, GENE), so we need to:
    1. Find SNPs near each locus position
    2. Match genes in our candidate set
    """
    print("\nMatching cS2G scores to candidates...")
    
    # Create a lookup dictionary: gene -> list of (snp, score) tuples
    gene_scores = {}
    for _, row in cs2g_scores.iterrows():
        gene = row.get('GENE', row.get('gene', ''))
        snp = row.get('SNP', row.get('snp', ''))
        score = row.get('cS2G', row.get('score', 0))
        
        if gene not in gene_scores:
            gene_scores[gene] = []
        gene_scores[gene].append((snp, score))
    
    print(f"  Built index for {len(gene_scores):,} genes")
    
    # Match to candidates
    # For simplicity, we use the max cS2G score for each gene
    candidates['score_cs2g'] = candidates['gene_symbol'].map(
        lambda g: max([s for _, s in gene_scores.get(g, [(None, np.nan)])], default=np.nan)
    )
    
    matched = candidates['score_cs2g'].notna().sum()
    print(f"  Matched {matched:,} / {len(candidates):,} candidates ({100*matched/len(candidates):.1f}%)")
    
    return candidates


def main():
    print("=" * 60)
    print("Adding cS2G Baseline to Evaluation Candidates")
    print("=" * 60)
    
    # Check if cS2G data exists
    if not CS2G_DIR.exists():
        print(f"\nError: cS2G data not found at {CS2G_DIR}")
        print("Please download from: https://zenodo.org/records/7754032")
        return
    
    # Load cS2G scores
    cs2g_scores = load_cs2g_scores()
    
    # Examine the data structure
    print(f"\ncS2G columns: {list(cs2g_scores.columns)}")
    print(f"Sample rows:")
    print(cs2g_scores.head())
    
    # Load evaluation candidates
    print(f"\nLoading candidates from {CANDIDATES_FILE}...")
    candidates = pd.read_csv(CANDIDATES_FILE, sep='\t')
    print(f"  Loaded {len(candidates):,} candidates")
    
    # Match scores
    candidates = match_scores_to_candidates(candidates, cs2g_scores)
    
    # Save output
    print(f"\nSaving to {OUTPUT_FILE}...")
    candidates.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("Summary Statistics")
    print("=" * 60)
    
    has_cs2g = candidates['score_cs2g'].notna()
    print(f"\nOverall coverage: {has_cs2g.sum():,} / {len(candidates):,} ({100*has_cs2g.mean():.1f}%)")
    
    if 'evidence_type' in candidates.columns:
        for etype in candidates['evidence_type'].unique():
            mask = candidates['evidence_type'] == etype
            coverage = has_cs2g[mask].mean()
            print(f"  {etype}: {100*coverage:.1f}% coverage")
    
    if has_cs2g.sum() > 0:
        print(f"\ncS2G score distribution (where available):")
        print(f"  Mean: {candidates.loc[has_cs2g, 'score_cs2g'].mean():.4f}")
        print(f"  Median: {candidates.loc[has_cs2g, 'score_cs2g'].median():.4f}")
        print(f"  Min: {candidates.loc[has_cs2g, 'score_cs2g'].min():.4f}")
        print(f"  Max: {candidates.loc[has_cs2g, 'score_cs2g'].max():.4f}")


if __name__ == "__main__":
    main()
