"""
Generate Locus-Gene Pairs for Post-2021 Benchmark
==================================================

Creates all candidate locus-gene pairs for the 63-locus post-2021 benchmark.

For each locus:
- Finds all genes within ±1Mb
- Annotates with ABC scores
- Annotates with eQTL scores  
- Computes distance scores

Output: post2021_locus_gene_pairs_annotated.tsv

Author: Mechanism-GWAS-Causal-Graphs team
Date: December 2025
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import gzip
from typing import Dict
import logging

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_benchmark() -> pd.DataFrame:
    """Load post-2021 independent benchmark."""
    benchmark_file = project_root / "data" / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
    df = pd.read_csv(benchmark_file, sep='\t')
    logger.info(f"Loaded benchmark: {len(df)} loci")
    return df


def load_gene_tss() -> pd.DataFrame:
    """Load gene TSS positions from GENCODE."""
    tss_file = project_root / "data" / "external" / "ensembl" / "gencode_tss_grch38.tsv"
    
    if not tss_file.exists():
        logger.error(f"TSS file not found: {tss_file}")
        return pd.DataFrame()
    
    df = pd.read_csv(tss_file, sep='\t')
    df = df.rename(columns={
        'gene_symbol': 'gene',
        'chromosome': 'chr',
        'tss_position': 'tss'
    })
    
    logger.info(f"Loaded TSS for {len(df):,} genes")
    return df


def create_locus_gene_pairs(benchmark_df: pd.DataFrame, gene_tss_df: pd.DataFrame, window_kb: int = 1000) -> pd.DataFrame:
    """
    Create all locus-gene pairs within window.
    
    Args:
        benchmark_df: Post-2021 benchmark
        gene_tss_df: Gene TSS annotations
        window_kb: Window size around lead SNP (default 1 Mb)
    
    Returns:
        DataFrame with locus-gene pairs
    """
    pairs = []
    window_bp = window_kb * 1000
    
    for _, locus in benchmark_df.iterrows():
        locus_chr = str(locus['chr'])
        locus_pos = locus['pos_hg38']
        locus_id = locus['locus_id']
        
        # Filter genes on same chromosome
        chr_genes = gene_tss_df[
            (gene_tss_df['chr'] == locus_chr) | 
            (gene_tss_df['chr'] == f"chr{locus_chr}") |
            (gene_tss_df['chr'].astype(str).str.replace('chr', '') == locus_chr)
        ]
        
        # Calculate distances and filter by window
        for _, gene in chr_genes.iterrows():
            distance = abs(int(gene['tss']) - locus_pos)
            
            if distance <= window_bp:
                pairs.append({
                    'locus_id': locus_id,
                    'gene': gene['gene'],
                    'chr': locus_chr,
                    'tss': gene['tss'],
                    'locus_pos': locus_pos,
                    'distance_to_lead': distance,
                    'within_1mb': True
                })
    
    pairs_df = pd.DataFrame(pairs)
    logger.info(f"Created {len(pairs_df):,} locus-gene pairs for {pairs_df['locus_id'].nunique()} loci")
    logger.info(f"Median genes per locus: {pairs_df.groupby('locus_id').size().median():.0f}")
    
    return pairs_df


def load_abc_predictions() -> pd.DataFrame:
    """Load ABC enhancer-gene predictions."""
    abc_file = project_root / "data" / "external" / "abc" / "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
    
    if not abc_file.exists():
        logger.warning(f"ABC file not found: {abc_file}")
        return pd.DataFrame()
    
    logger.info("Loading ABC predictions...")
    with gzip.open(abc_file, 'rt') as f:
        abc_df = pd.read_csv(f, sep='\t')
    
    logger.info(f"Loaded {len(abc_df):,} ABC predictions")
    return abc_df


def annotate_with_abc(pairs_df: pd.DataFrame, abc_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate pairs with ABC scores."""
    if abc_df.empty:
        pairs_df['abc_score'] = 0.0
        return pairs_df
    
    # Create gene → max ABC score mapping
    gene_abc = {}
    for _, row in abc_df.iterrows():
        gene = row['TargetGene']
        score = row['ABC.Score.Numerator']  # Using raw numerator as in original script
        
        if gene not in gene_abc or score > gene_abc[gene]:
            gene_abc[gene] = score
    
    # Annotate pairs
    pairs_df['abc_score'] = pairs_df['gene'].map(gene_abc).fillna(0.0)
    
    n_with_abc = (pairs_df['abc_score'] > 0).sum()
    logger.info(f"ABC annotation: {n_with_abc}/{len(pairs_df)} pairs ({100*n_with_abc/len(pairs_df):.1f}%) with ABC > 0")
    logger.info(f"ABC score range: {pairs_df['abc_score'].min():.4f} - {pairs_df['abc_score'].max():.4f}")
    
    return pairs_df


def load_eqtl_scores() -> pd.DataFrame:
    """Load GTEx eQTL scores."""
    eqtl_file = project_root / "data" / "external" / "gtex" / "gtex_gene_eqtl_scores.tsv"
    
    if not eqtl_file.exists():
        logger.warning(f"eQTL file not found: {eqtl_file}")
        return pd.DataFrame()
    
    df = pd.read_csv(eqtl_file, sep='\t')
    logger.info(f"Loaded eQTL scores for {len(df):,} genes")
    
    return df


def annotate_with_eqtl(pairs_df: pd.DataFrame, eqtl_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate pairs with eQTL scores."""
    if eqtl_df.empty:
        pairs_df['eqtl_score'] = 0.0
        return pairs_df
    
    # Create gene → eQTL score mapping
    gene_eqtl = dict(zip(eqtl_df['gene_symbol'], eqtl_df['eqtl_score']))
    
    # Annotate pairs
    pairs_df['eqtl_score'] = pairs_df['gene'].map(gene_eqtl).fillna(0.0)
    
    n_with_eqtl = (pairs_df['eqtl_score'] > 0).sum()
    logger.info(f"eQTL annotation: {n_with_eqtl}/{len(pairs_df)} pairs ({100*n_with_eqtl/len(pairs_df):.1f}%) with eQTL > 0")
    logger.info(f"eQTL score range: {pairs_df['eqtl_score'].min():.4f} - {pairs_df['eqtl_score'].max():.4f}")
    
    return pairs_df


def compute_distance_scores(pairs_df: pd.DataFrame) -> pd.DataFrame:
    """Compute distance-based scores (higher = closer)."""
    max_dist = pairs_df['distance_to_lead'].max()
    pairs_df['distance_score'] = 1.0 - (pairs_df['distance_to_lead'] / max_dist)
    
    logger.info(f"Distance score range: {pairs_df['distance_score'].min():.4f} - {pairs_df['distance_score'].max():.4f}")
    
    return pairs_df


def main():
    """Main execution."""
    logger.info("="*80)
    logger.info("GENERATING LOCUS-GENE PAIRS FOR POST-2021 BENCHMARK")
    logger.info("="*80)
    
    # Load data
    logger.info("\n1. Loading benchmark and gene annotations...")
    benchmark_df = load_benchmark()
    gene_tss_df = load_gene_tss()
    
    if gene_tss_df.empty:
        logger.error("Cannot proceed without gene TSS annotations")
        return
    
    # Create pairs
    logger.info("\n2. Creating locus-gene pairs...")
    pairs_df = create_locus_gene_pairs(benchmark_df, gene_tss_df)
    
    # Load functional annotations
    logger.info("\n3. Loading ABC predictions...")
    abc_df = load_abc_predictions()
    
    logger.info("\n4. Loading eQTL scores...")
    eqtl_df = load_eqtl_scores()
    
    # Annotate pairs
    logger.info("\n5. Annotating with ABC scores...")
    pairs_df = annotate_with_abc(pairs_df, abc_df)
    
    logger.info("\n6. Annotating with eQTL scores...")
    pairs_df = annotate_with_eqtl(pairs_df, eqtl_df)
    
    logger.info("\n7. Computing distance scores...")
    pairs_df = compute_distance_scores(pairs_df)
    
    # Save results
    output_file = project_root / "data" / "processed" / "baselines" / "post2021_locus_gene_pairs_annotated.tsv"
    pairs_df.to_csv(output_file, sep='\t', index=False)
    
    logger.info("\n" + "="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    logger.info(f"Total pairs: {len(pairs_df):,}")
    logger.info(f"Unique loci: {pairs_df['locus_id'].nunique()}")
    logger.info(f"Unique genes: {pairs_df['gene'].nunique()}")
    logger.info(f"ABC coverage: {(pairs_df['abc_score'] > 0).sum():,} pairs ({100*(pairs_df['abc_score'] > 0).sum()/len(pairs_df):.1f}%)")
    logger.info(f"eQTL coverage: {(pairs_df['eqtl_score'] > 0).sum():,} pairs ({100*(pairs_df['eqtl_score'] > 0).sum()/len(pairs_df):.1f}%)")
    logger.info(f"Median genes/locus: {pairs_df.groupby('locus_id').size().median():.0f}")
    logger.info(f"\nOutput: {output_file}")
    logger.info("="*80)


if __name__ == "__main__":
    main()
