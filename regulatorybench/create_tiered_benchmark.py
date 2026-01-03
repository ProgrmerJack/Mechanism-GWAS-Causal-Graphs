#!/usr/bin/env python3
"""
Create Tiered Gold Standard Benchmark for Nature Genetics Article

Combines:
- Tier-0 (Platinum): G2P-only Mendelian genes (123 genes, zero contamination)
- Tier-1 (Clinical): ClinVar pathogenic variants (2,040 genes)
- Tier-2 (Experimental): STING-seq/CRISPRi validation (when available)

Author: Generated for Nature Genetics v3 Platinum manuscript
Date: 2025-12-21
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Dict, Set

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def create_tiered_benchmark(
    g2p_file: str = "benchmarks/evidence_curation/mendelian_genes_g2p_only.tsv",
    clinvar_file: str = "benchmarks/mendelian_genes_curated.tsv",
    output_dir: str = "benchmarks"
) -> None:
    """
    Create tiered benchmark structure.
    
    Parameters
    ----------
    g2p_file : str
        Path to G2P-only genes (Tier-0 Platinum)
    clinvar_file : str
        Path to ClinVar curated genes (Tier-1)
    output_dir : str
        Output directory for tiered benchmarks
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("="*60)
    logger.info("Creating Tiered Gold Standard")
    logger.info("="*60)
    
    # Load Tier-0 (Platinum G2P)
    logger.info("\nTier-0 (Platinum Mendelian): Loading G2P genes...")
    tier0_df = pd.read_csv(g2p_file, sep='\t')
    tier0_df['tier'] = 'tier0_platinum'
    tier0_df['tier_label'] = 'Platinum Mendelian (G2P)'
    logger.info(f"  Loaded {len(tier0_df)} Tier-0 genes")
    
    # Load Tier-1 (ClinVar Clinical)
    logger.info("\nTier-1 (Clinical Genetics): Loading ClinVar genes...")
    tier1_df = pd.read_csv(clinvar_file, sep='\t')
    tier1_df['tier'] = 'tier1_clinical'
    tier1_df['tier_label'] = 'Clinical Genetics (ClinVar)'
    logger.info(f"  Loaded {len(tier1_df)} Tier-1 genes")
    
    # Identify overlap
    tier0_genes = set(tier0_df['gene_symbol'])
    tier1_genes = set(tier1_df['gene_symbol'])
    overlap = tier0_genes & tier1_genes
    logger.info(f"\n  Overlap between Tier-0 and Tier-1: {len(overlap)} genes")
    
    # Create combined Tier 0+1 (keeping tier labels)
    combined_df = pd.concat([tier0_df, tier1_df], ignore_index=True)
    logger.info(f"\nCombined Tier 0+1: {len(combined_df)} total evidence records")
    
    # Save individual tiers
    tier0_output = output_dir / "task_a_evidence_tier0_platinum.tsv"
    tier0_df.to_csv(tier0_output, sep='\t', index=False)
    logger.info(f"\nSaved Tier-0 to: {tier0_output}")
    
    tier1_output = output_dir / "task_a_evidence_tier1_clinical.tsv"
    tier1_df.to_csv(tier1_output, sep='\t', index=False)
    logger.info(f"Saved Tier-1 to: {tier1_output}")
    
    # Save combined tiered benchmark
    combined_output = output_dir / "task_a_evidence_tiered.tsv"
    combined_df.to_csv(combined_output, sep='\t', index=False)
    logger.info(f"Saved combined tiered evidence to: {combined_output}")
    
    # Generate summary statistics
    logger.info("\n" + "="*60)
    logger.info("TIERED BENCHMARK SUMMARY")
    logger.info("="*60)
    logger.info(f"\nTier-0 (Platinum): {len(tier0_genes)} unique genes")
    logger.info(f"Tier-1 (Clinical): {len(tier1_genes)} unique genes")
    logger.info(f"Overlap: {len(overlap)} genes ({100*len(overlap)/len(tier0_genes):.1f}% of Tier-0)")
    logger.info(f"Total unique genes: {len(tier0_genes | tier1_genes)}")
    
    logger.info("\nEvidence source breakdown:")
    for source in combined_df['evidence_source'].unique():
        count = (combined_df['evidence_source'] == source).sum()
        logger.info(f"  {source}: {count}")
    
    logger.info("\nConfidence breakdown:")
    for conf in combined_df['confidence'].value_counts().items():
        logger.info(f"  {conf[0]}: {conf[1]}")
    
    logger.info("\n" + "="*60)
    logger.info("Tiered benchmark creation complete!")
    logger.info("="*60)
    
    # Return summary for verification
    summary = {
        'tier0_genes': len(tier0_genes),
        'tier1_genes': len(tier1_genes),
        'overlap': len(overlap),
        'total_unique': len(tier0_genes | tier1_genes),
        'tier0_file': str(tier0_output),
        'tier1_file': str(tier1_output),
        'combined_file': str(combined_output)
    }
    
    return summary


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Create tiered gold standard benchmark")
    parser.add_argument("--g2p-file", default="benchmarks/evidence_curation/mendelian_genes_g2p_only.tsv",
                       help="Path to G2P-only genes (Tier-0)")
    parser.add_argument("--clinvar-file", default="benchmarks/mendelian_genes_curated.tsv",
                       help="Path to ClinVar curated genes (Tier-1)")
    parser.add_argument("--output-dir", default="benchmarks",
                       help="Output directory")
    
    args = parser.parse_args()
    
    summary = create_tiered_benchmark(
        g2p_file=args.g2p_file,
        clinvar_file=args.clinvar_file,
        output_dir=args.output_dir
    )
    
    print("\nSummary:", summary)
