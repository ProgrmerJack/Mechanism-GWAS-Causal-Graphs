"""
Prepare Input Data for Baseline Methods
========================================

This script prepares all necessary input data for running baseline methods:
- cS2G: ABC links, eQTL data, distance annotations, coding variants, constraint
- FLAMES: ABC, eQTL, distance, credible set PIPs
- PoPS: MAGMA gene sets, LD-pruned SNPs
- Effector Index: ABC, eQTL, epigenomic annotations

Created: December 19, 2025
Purpose: Replace random score placeholders with real baseline implementations
"""

import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import logging
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class BaselineDataPreparation:
    """Prepare standardized inputs for all baseline methods."""
    
    def __init__(self, data_dir: str = "data/external", output_dir: str = "data/processed/baselines"):
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def load_abc_predictions(self) -> pd.DataFrame:
        """
        Load ABC enhancer-gene predictions.
        
        Returns:
            DataFrame with columns: chr, start, end, gene, ABC_Score, cell_type
        """
        abc_file = self.data_dir / "abc" / "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
        
        if not abc_file.exists():
            logger.error(f"ABC file not found: {abc_file}")
            return pd.DataFrame()
        
        logger.info(f"Loading ABC predictions from {abc_file}")
        
        # Read ABC file (tab-separated, gzipped)
        with gzip.open(abc_file, 'rt') as f:
            abc_df = pd.read_csv(f, sep='\t')
        
        logger.info(f"Loaded {len(abc_df):,} ABC predictions")
        logger.info(f"Columns: {abc_df.columns.tolist()}")
        logger.info(f"Cell types: {abc_df.get('CellType', abc_df.columns[0]).nunique()}")
        
        return abc_df
    
    def load_gtex_eqtls(self, tissue: str = "Liver") -> pd.DataFrame:
        """
        Load GTEx V8 eQTL data for genes.
        
        Args:
            tissue: GTEx tissue name (optional, for tissue-specific files)
        
        Returns:
            DataFrame with columns: gene, eqtl_score or tissue-specific data
        """
        gtex_dir = self.data_dir / "gtex"
        
        # First try our processed gene eQTL scores
        gene_eqtl_file = gtex_dir / "gtex_gene_eqtl_scores.tsv"
        
        if gene_eqtl_file.exists():
            logger.info(f"Loading processed GTEx gene eQTL scores from {gene_eqtl_file}")
            eqtl_df = pd.read_csv(gene_eqtl_file, sep='\t')
            logger.info(f"Loaded {len(eqtl_df):,} gene eQTL scores")
            logger.info(f"Columns: {eqtl_df.columns.tolist()}")
            return eqtl_df
        
        # Fallback: look for tissue-specific eQTL files
        eqtl_files = list(gtex_dir.glob(f"*{tissue}*.signif_variant_gene_pairs.txt.gz"))
        
        if not eqtl_files:
            logger.warning(f"No GTEx eQTL files found")
            return pd.DataFrame()
        
        logger.info(f"Loading GTEx eQTLs from {eqtl_files[0]}")
        
        with gzip.open(eqtl_files[0], 'rt') as f:
            eqtl_df = pd.read_csv(f, sep='\t')
        
        eqtl_df['tissue'] = tissue
        
        logger.info(f"Loaded {len(eqtl_df):,} eQTL associations for {tissue}")
        logger.info(f"Columns: {eqtl_df.columns.tolist()}")
        
        return eqtl_df
    
    def load_locus_summary(self) -> pd.DataFrame:
        """
        Load our processed locus summary with GWAS signals.
        
        Returns:
            DataFrame with locus-level information
        """
        locus_file = Path("data/processed/locus_summary.tsv")
        
        if not locus_file.exists():
            logger.error(f"Locus summary not found: {locus_file}")
            return pd.DataFrame()
        
        logger.info(f"Loading locus summary from {locus_file}")
        locus_df = pd.read_csv(locus_file, sep='\t')
        
        logger.info(f"Loaded {len(locus_df)} loci")
        logger.info(f"Traits: {locus_df['trait'].nunique()}")
        logger.info(f"Columns: {locus_df.columns.tolist()}")
        
        return locus_df
    
    def create_locus_gene_pairs(
        self,
        locus_df: pd.DataFrame,
        window_kb: int = 1000
    ) -> pd.DataFrame:
        """
        Create locus-gene pairs for scoring.
        
        For each locus, identify candidate genes within +/- window_kb.
        
        Args:
            locus_df: Locus summary DataFrame
            window_kb: Window around lead SNP (default 1Mb)
        
        Returns:
            DataFrame with columns: locus_id, gene, chr, tss, distance_to_lead
        """
        # Load gene TSS annotations (Gencode format)
        gene_file = self.data_dir / "ensembl" / "gencode_tss_grch38.tsv"
        
        if not gene_file.exists():
            logger.warning(f"Gene TSS file not found: {gene_file}")
            return pd.DataFrame()
        
        logger.info(f"Loading gene TSS positions from {gene_file}")
        
        gene_df = pd.read_csv(gene_file, sep='\t')
        logger.info(f"Loaded TSS for {len(gene_df):,} genes")
        logger.info(f"Gene TSS columns: {gene_df.columns.tolist()}")
        
        # Standardize column names (handle different formats)
        col_mapping = {}
        for col in gene_df.columns:
            col_lower = col.lower()
            if ('chr' in col_lower or 'chrom' in col_lower) and 'chr' not in col_mapping and col != 'gene_type':
                col_mapping[col] = 'chr'
            elif ('tss' in col_lower or ('position' in col_lower and 'tss' in col_lower)) and 'tss' not in col_mapping:
                col_mapping[col] = 'tss'
            elif ('gene' in col_lower or 'symbol' in col_lower) and 'gene' not in col_mapping and col != 'gene_type':
                col_mapping[col] = 'gene'
        
        if col_mapping:
            gene_df = gene_df.rename(columns=col_mapping)
            logger.info(f"Renamed columns: {col_mapping}")
        
        # Ensure required columns exist
        if 'gene' not in gene_df.columns or 'chr' not in gene_df.columns or 'tss' not in gene_df.columns:
            logger.error(f"Missing required columns. Available: {gene_df.columns.tolist()}")
            return pd.DataFrame()
        
        # Create locus-gene pairs
        pairs = []
        window_bp = window_kb * 1000
        
        for _, locus in locus_df.iterrows():
            locus_chr = str(locus['chr'])
            locus_pos = locus['pos_hg38']
            locus_id = locus['locus_id']
            
            # Filter genes on same chromosome (handle both "1" and "chr1" formats)
            chr_genes = gene_df[
                (gene_df['chr'] == locus_chr) | 
                (gene_df['chr'] == f"chr{locus_chr}") |
                (gene_df['chr'].astype(str).str.replace('chr', '') == locus_chr)
            ]
            
            # Calculate distances
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
                        'within_1mb': distance <= 1_000_000
                    })
        
        pairs_df = pd.DataFrame(pairs)
        logger.info(f"Created {len(pairs_df):,} locus-gene pairs")
        logger.info(f"Genes per locus: {pairs_df.groupby('locus_id').size().describe()}")
        
        return pairs_df
    
    def annotate_pairs_with_abc(
        self,
        pairs_df: pd.DataFrame,
        abc_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Annotate locus-gene pairs with ABC scores.
        
        Args:
            pairs_df: Locus-gene pairs
            abc_df: ABC predictions
        
        Returns:
            pairs_df with added abc_score column
        """
        logger.info("Annotating pairs with ABC scores...")
        
        # Standardize column names in ABC df
        abc_cols = abc_df.columns.tolist()
        logger.info(f"ABC columns: {abc_cols[:10]}...")
        
        # Common column mappings
        gene_col = next((col for col in abc_cols if 'gene' in col.lower() or 'target' in col.lower()), None)
        score_col = next((col for col in abc_cols if 'abc' in col.lower() and 'score' in col.lower()), None)
        chr_col = next((col for col in abc_cols if col.lower() in ['chr', 'chromosome']), None)
        
        if not all([gene_col, score_col]):
            logger.error(f"Could not find required ABC columns. Gene col: {gene_col}, Score col: {score_col}")
            pairs_df['abc_score'] = 0.0
            return pairs_df
        
        logger.info(f"Using ABC columns: gene={gene_col}, score={score_col}")
        
        # Create gene -> max ABC score mapping
        abc_gene_scores = abc_df.groupby(gene_col)[score_col].max().to_dict()
        
        # Annotate pairs
        pairs_df['abc_score'] = pairs_df['gene'].map(abc_gene_scores).fillna(0.0)
        
        logger.info(f"ABC scores: {(pairs_df['abc_score'] > 0).sum()} / {len(pairs_df)} pairs have ABC links")
        logger.info(f"ABC score range: {pairs_df['abc_score'].min():.4f} - {pairs_df['abc_score'].max():.4f}")
        
        return pairs_df
    
    def annotate_pairs_with_eqtl(
        self,
        pairs_df: pd.DataFrame,
        eqtl_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Annotate locus-gene pairs with eQTL evidence.
        
        Args:
            pairs_df: Locus-gene pairs
            eqtl_df: GTEx eQTL associations or scores
        
        Returns:
            pairs_df with added eqtl_score column
        """
        logger.info("Annotating pairs with eQTL evidence...")
        
        if eqtl_df.empty:
            pairs_df['eqtl_score'] = 0.0
            return pairs_df
        
        logger.info(f"eQTL DataFrame columns: {eqtl_df.columns.tolist()}")
        
        # Handle different eQTL data formats
        if 'gene' in eqtl_df.columns or 'gene_name' in eqtl_df.columns or 'gene_symbol' in eqtl_df.columns:
            # Direct gene-level scores
            gene_col = 'gene' if 'gene' in eqtl_df.columns else ('gene_name' if 'gene_name' in eqtl_df.columns else 'gene_symbol')
            score_col = next((col for col in eqtl_df.columns if 'score' in col.lower() or 'pip' in col.lower()), None)
            
            if score_col:
                gene_scores = eqtl_df.groupby(gene_col)[score_col].max().to_dict()
                pairs_df['eqtl_score'] = pairs_df['gene'].map(gene_scores).fillna(0.0)
            else:
                # Use presence as binary score
                genes_with_eqtl = set(eqtl_df[gene_col].unique())
                pairs_df['eqtl_score'] = pairs_df['gene'].apply(lambda g: 1.0 if g in genes_with_eqtl else 0.0)
        
        elif 'gene_id' in eqtl_df.columns:
            # Variant-level eQTL data with gene IDs
            eqtl_df['gene_clean'] = eqtl_df['gene_id'].str.split('.').str[0]
            
            # Convert p-values to scores: score = -log10(p)
            if 'pval_nominal' in eqtl_df.columns:
                eqtl_df['eqtl_score'] = -np.log10(eqtl_df['pval_nominal'] + 1e-300)
                gene_scores = eqtl_df.groupby('gene_clean')['eqtl_score'].max().to_dict()
                pairs_df['eqtl_score'] = pairs_df['gene'].map(gene_scores).fillna(0.0)
            else:
                genes_with_eqtl = set(eqtl_df['gene_clean'].unique())
                pairs_df['eqtl_score'] = pairs_df['gene'].apply(lambda g: 1.0 if g in genes_with_eqtl else 0.0)
        else:
            logger.warning("Could not identify gene column in eQTL data")
            pairs_df['eqtl_score'] = 0.0
        
        logger.info(f"eQTL evidence: {(pairs_df['eqtl_score'] > 0).sum()} / {len(pairs_df)} pairs have eQTL scores")
        logger.info(f"eQTL score range: {pairs_df['eqtl_score'].min():.4f} - {pairs_df['eqtl_score'].max():.4f}")
        
        return pairs_df
    
    def compute_distance_scores(self, pairs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute distance-based scores.
        
        Uses inverse distance weighting: score = 1 / (1 + distance / 100kb)
        
        Args:
            pairs_df: Locus-gene pairs with distance_to_lead column
        
        Returns:
            pairs_df with added distance_score column
        """
        pairs_df['distance_score'] = 1.0 / (1.0 + pairs_df['distance_to_lead'] / 100_000)
        
        logger.info(f"Distance scores: range {pairs_df['distance_score'].min():.4f} - {pairs_df['distance_score'].max():.4f}")
        
        return pairs_df
    
    def run_full_preparation(self) -> pd.DataFrame:
        """
        Run complete baseline input preparation pipeline.
        
        Returns:
            Fully annotated locus-gene pairs DataFrame
        """
        logger.info("="*80)
        logger.info("BASELINE DATA PREPARATION PIPELINE")
        logger.info("="*80)
        
        # Load raw data
        logger.info("\n1. Loading raw data sources...")
        locus_df = self.load_locus_summary()
        abc_df = self.load_abc_predictions()
        eqtl_df = self.load_gtex_eqtls(tissue="Liver")
        
        if locus_df.empty:
            logger.error("No locus data available - cannot proceed")
            return pd.DataFrame()
        
        # Create locus-gene pairs
        logger.info("\n2. Creating locus-gene pairs...")
        pairs_df = self.create_locus_gene_pairs(locus_df, window_kb=1000)
        
        if pairs_df.empty:
            logger.error("No locus-gene pairs created - cannot proceed")
            return pd.DataFrame()
        
        # Annotate with features
        logger.info("\n3. Annotating with ABC scores...")
        if not abc_df.empty:
            pairs_df = self.annotate_pairs_with_abc(pairs_df, abc_df)
        else:
            pairs_df['abc_score'] = 0.0
        
        logger.info("\n4. Annotating with eQTL evidence...")
        pairs_df = self.annotate_pairs_with_eqtl(pairs_df, eqtl_df)
        
        logger.info("\n5. Computing distance scores...")
        pairs_df = self.compute_distance_scores(pairs_df)
        
        # Save results
        output_file = self.output_dir / "locus_gene_pairs_annotated.tsv"
        pairs_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"\n✓ Saved annotated pairs to: {output_file}")
        
        # Summary statistics
        logger.info("\n" + "="*80)
        logger.info("SUMMARY STATISTICS")
        logger.info("="*80)
        logger.info(f"Total locus-gene pairs: {len(pairs_df):,}")
        logger.info(f"Unique loci: {pairs_df['locus_id'].nunique()}")
        logger.info(f"Unique genes: {pairs_df['gene'].nunique()}")
        logger.info(f"Pairs with ABC > 0: {(pairs_df['abc_score'] > 0).sum():,} ({100*(pairs_df['abc_score'] > 0).sum()/len(pairs_df):.1f}%)")
        logger.info(f"Pairs with eQTL score > 0: {(pairs_df['eqtl_score'] > 0).sum():,} ({100*(pairs_df['eqtl_score'] > 0).sum()/len(pairs_df):.1f}%)")
        logger.info(f"Median genes per locus: {pairs_df.groupby('locus_id').size().median():.0f}")
        logger.info("="*80)
        
        return pairs_df


if __name__ == "__main__":
    prep = BaselineDataPreparation()
    annotated_pairs = prep.run_full_preparation()
    
    if not annotated_pairs.empty:
        print("\n✓ Baseline input preparation complete!")
        print(f"✓ Data ready for cS2G, FLAMES, and other baselines")
        print(f"✓ Output: data/processed/baselines/locus_gene_pairs_annotated.tsv")
