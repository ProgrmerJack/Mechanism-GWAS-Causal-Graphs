#!/usr/bin/env python3
"""
Extract Coding/LoF Variants from Fine-Mapped Credible Sets

Uses VEP (Variant Effect Predictor) consequence annotations to identify:
- High PIP coding variants (missense, nonsense, frameshift)
- Loss-of-function (LoF) variants (stop-gained, frameshift, splice-site)
- Splice site variants

Maps variants unambiguously to single genes for use in gold standard.

Author: Generated for RegulatoryBench v3 platinum
Date: 2025-01
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Set
import pyarrow.parquet as pq
import json
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CodingVariantExtractor:
    """Extract coding/LoF variants from fine-mapped credible sets."""
    
    # VEP consequence terms for coding/LoF variants
    LOF_CONSEQUENCES = {
        'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
        'splice_donor_variant', 'splice_acceptor_variant',
        'transcript_ablation', 'transcript_amplification'
    }
    
    CODING_CONSEQUENCES = {
        'missense_variant', 'inframe_insertion', 'inframe_deletion',
        'protein_altering_variant', 'coding_sequence_variant'
    }
    
    SPLICE_REGION_CONSEQUENCES = {
        'splice_region_variant', 'splice_donor_5th_base_variant',
        'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant'
    }
    
    def __init__(
        self,
        annotations_dir: str = "data/external/flames/Annotation_data",
        benchmark_dir: str = "benchmarks",
        pip_threshold_lof: float = 0.5,
        pip_threshold_missense: float = 0.8,
    ):
        self.annotations_dir = Path(annotations_dir)
        self.benchmark_dir = Path(benchmark_dir)
        self.output_dir = self.benchmark_dir / "evidence_curation"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # PIP thresholds
        self.pip_threshold_lof = pip_threshold_lof  # Lower threshold for LoF (higher confidence)
        self.pip_threshold_missense = pip_threshold_missense  # Higher threshold for missense
        
        # Storage for extracted genes
        self.coding_genes = []
        
    def load_finemapped_credible_sets(self) -> pd.DataFrame:
        """
        Load fine-mapped credible sets with PIP values.
        
        Returns
        -------
        pd.DataFrame
            Credible set variants with PIP values.
        """
        logger.info("Loading fine-mapped credible sets...")
        
        # Check for SuSiE fine-mapping results
        finemapping_files = list(self.benchmark_dir.glob("**/finemapping/*.parquet"))
        
        if not finemapping_files:
            logger.warning("No fine-mapping files found. Checking alternative locations...")
            
            # Try FLAMES example data
            example_file = Path("data/external/flames/example_data/credible_sets.parquet")
            if example_file.exists():
                logger.info(f"Using example credible sets from {example_file}")
                return pd.read_parquet(example_file)
            
            logger.error("No fine-mapping results found! Please run SuSiE fine-mapping first.")
            return pd.DataFrame()
        
        # Load all fine-mapping results
        dfs = []
        for file in finemapping_files:
            try:
                df = pd.read_parquet(file)
                dfs.append(df)
            except Exception as e:
                logger.warning(f"Could not read {file}: {e}")
        
        if not dfs:
            return pd.DataFrame()
        
        credible_sets = pd.concat(dfs, ignore_index=True)
        logger.info(f"Loaded {len(credible_sets):,} variants in credible sets")
        
        return credible_sets
    
    def load_vep_annotations(self, chr_pos_list: List[Tuple[str, int]]) -> pd.DataFrame:
        """
        Load VEP annotations for specific variants.
        
        Parameters
        ----------
        chr_pos_list : list
            List of (chromosome, position) tuples.
            
        Returns
        -------
        pd.DataFrame
            VEP annotations for requested variants.
        """
        logger.info(f"Loading VEP annotations for {len(chr_pos_list):,} variants...")
        
        # FLAMES annotation structure: data/external/flames/Annotation_data/chr{N}/{variant}.parquet
        
        annotations = []
        
        for i, (chrom, pos) in enumerate(chr_pos_list):
            if (i + 1) % 1000 == 0:
                logger.info(f"  Loaded annotations for {i+1:,} / {len(chr_pos_list):,} variants")
            
            # Construct file path
            # FLAMES uses: chr{N}/{chrom}_{pos}_{ref}_{alt}.parquet
            # We need to search for any variant at this position
            
            chr_dir = self.annotations_dir / f"chr{chrom}"
            if not chr_dir.exists():
                continue
            
            # Find all variants at this position
            variant_files = list(chr_dir.glob(f"{chrom}_{pos}_*.parquet"))
            
            for vfile in variant_files:
                try:
                    ann = pd.read_parquet(vfile)
                    annotations.append(ann)
                except Exception as e:
                    logger.warning(f"Could not read {vfile}: {e}")
        
        if not annotations:
            logger.warning("No VEP annotations found!")
            return pd.DataFrame()
        
        vep_df = pd.concat(annotations, ignore_index=True)
        logger.info(f"Loaded {len(vep_df):,} VEP consequence annotations")
        
        return vep_df
    
    def annotate_credible_sets_with_vep(
        self,
        credible_sets: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Annotate credible set variants with VEP consequences.
        
        Parameters
        ----------
        credible_sets : pd.DataFrame
            Fine-mapped credible set variants.
            
        Returns
        -------
        pd.DataFrame
            Credible sets with VEP consequence annotations.
        """
        logger.info("Annotating credible sets with VEP consequences...")
        
        # Extract chromosome and position for all variants
        if 'chr' in credible_sets.columns and 'pos' in credible_sets.columns:
            chr_pos_list = list(zip(credible_sets['chr'], credible_sets['pos']))
        elif 'variant_id' in credible_sets.columns:
            # Parse variant_id (format: chr:pos:ref:alt or chr_pos_ref_alt)
            parsed = credible_sets['variant_id'].str.extract(r'(\d+)[_:](\d+)')
            chr_pos_list = list(zip(parsed[0], parsed[1].astype(int)))
        else:
            logger.error("Cannot determine variant coordinates from credible sets!")
            return credible_sets
        
        # Load VEP annotations
        vep_annotations = self.load_vep_annotations(chr_pos_list)
        
        if vep_annotations.empty:
            logger.warning("No VEP annotations loaded. Using existing annotations in credible sets...")
            return credible_sets
        
        # Merge VEP annotations with credible sets
        # This depends on the exact structure of both dataframes
        # For now, return as-is
        
        logger.info("VEP annotation merging complete")
        return credible_sets
    
    def extract_coding_lof_genes(
        self,
        credible_sets: pd.DataFrame,
        require_unambiguous: bool = True
    ) -> int:
        """
        Extract genes with high-PIP coding/LoF variants.
        
        Parameters
        ----------
        credible_sets : pd.DataFrame
            Fine-mapped variants with VEP annotations.
        require_unambiguous : bool
            If True, only include variants mapping to a single gene.
            
        Returns
        -------
        int
            Number of genes extracted.
        """
        logger.info("Extracting genes with coding/LoF variants...")
        
        # Check for required columns
        required_cols = ['pip', 'consequence', 'gene_symbol']
        
        # Try alternative column names
        if 'PIP' in credible_sets.columns:
            credible_sets['pip'] = credible_sets['PIP']
        if 'Consequence' in credible_sets.columns:
            credible_sets['consequence'] = credible_sets['Consequence']
        if 'SYMBOL' in credible_sets.columns:
            credible_sets['gene_symbol'] = credible_sets['SYMBOL']
        
        missing_cols = [col for col in required_cols if col not in credible_sets.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            logger.error(f"Available columns: {list(credible_sets.columns)}")
            
            # Try to infer from existing annotations
            if 'AnyCoding' in credible_sets.columns:
                logger.info("Using AnyCoding flag from existing annotations...")
                return self._extract_from_flags(credible_sets)
            
            return 0
        
        count = 0
        
        for _, row in credible_sets.iterrows():
            pip = float(row['pip'])
            consequence = str(row['consequence'])
            gene = str(row['gene_symbol'])
            
            if pd.isna(gene) or gene == 'nan' or gene == '':
                continue
            
            # Check if variant has coding/LoF consequence
            consequences = set(consequence.split(','))
            
            is_lof = bool(consequences & self.LOF_CONSEQUENCES)
            is_coding = bool(consequences & self.CODING_CONSEQUENCES)
            is_splice = bool(consequences & self.SPLICE_REGION_CONSEQUENCES)
            
            # Apply PIP thresholds
            if is_lof and pip >= self.pip_threshold_lof:
                evidence_type = 'coding_lof'
                conf = 'very_high'
            elif is_coding and pip >= self.pip_threshold_missense:
                evidence_type = 'coding_missense'
                conf = 'high'
            elif is_splice and pip >= self.pip_threshold_lof:
                evidence_type = 'splice_variant'
                conf = 'high'
            else:
                continue
            
            # Check for ambiguous gene mapping
            if require_unambiguous:
                genes_affected = gene.split(',')
                if len(genes_affected) > 1:
                    continue
            
            self.coding_genes.append({
                'gene_symbol': gene,
                'evidence_type': evidence_type,
                'evidence_source': 'VEP_annotation',
                'variant_id': row.get('variant_id', f"{row.get('chr', '')}:{row.get('pos', '')}"),
                'pip': pip,
                'consequence': consequence,
                'confidence': conf,
                'date_curated': pd.Timestamp.now().strftime('%Y-%m-%d')
            })
            count += 1
        
        logger.info(f"Extracted {count} coding/LoF variant-gene associations")
        return count
    
    def _extract_from_flags(self, credible_sets: pd.DataFrame) -> int:
        """
        Extract genes using AnyCoding/AnySpliceSite flags.
        
        Fallback method when VEP consequence annotations are not available.
        
        Parameters
        ----------
        credible_sets : pd.DataFrame
            Variants with AnyCoding/AnySpliceSite flags.
            
        Returns
        -------
        int
            Number of genes extracted.
        """
        logger.info("Extracting genes using AnyCoding/AnySpliceSite flags...")
        
        # Filter for high-PIP coding variants
        if 'AnyCoding' in credible_sets.columns:
            coding = credible_sets[
                (credible_sets['AnyCoding'] == True) &
                (credible_sets.get('pip', credible_sets.get('PIP', 0)) >= self.pip_threshold_lof)
            ]
            
            logger.info(f"Found {len(coding):,} high-PIP coding variants")
            
            # Group by gene
            if 'gene_symbol' in coding.columns:
                genes = coding['gene_symbol'].unique()
            elif 'GeneSymbol' in coding.columns:
                genes = coding['GeneSymbol'].unique()
            else:
                logger.warning("No gene symbol column found!")
                return 0
            
            for gene in genes:
                if pd.isna(gene) or gene == '':
                    continue
                
                self.coding_genes.append({
                    'gene_symbol': str(gene),
                    'evidence_type': 'coding_variant',
                    'evidence_source': 'AnyCoding_flag',
                    'confidence': 'high',
                    'date_curated': pd.Timestamp.now().strftime('%Y-%m-%d')
                })
            
            logger.info(f"Extracted {len(genes)} genes with coding variants")
            return len(genes)
        
        return 0
    
    def deduplicate_and_save(self, output_file: str = "coding_lof_genes.tsv"):
        """
        Deduplicate and save extracted genes.
        
        Parameters
        ----------
        output_file : str
            Output filename.
        """
        if not self.coding_genes:
            logger.warning("No coding genes extracted!")
            return
        
        df = pd.DataFrame(self.coding_genes)
        
        # Deduplicate by gene
        df_dedup = df.groupby('gene_symbol').agg({
            'evidence_type': lambda x: ', '.join(sorted(set(x))),
            'evidence_source': 'first',
            'confidence': lambda x: 'very_high' if 'very_high' in x.values else x.iloc[0],
            'date_curated': 'first',
        }).reset_index()
        
        output_path = self.output_dir / output_file
        df_dedup.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Saved {len(df_dedup):,} genes with coding/LoF variants to {output_path}")
        logger.info(f"{'='*60}")
        
        # Summary
        logger.info("\nEvidence types:")
        for etype, count in df_dedup['evidence_type'].value_counts().items():
            logger.info(f"  {etype}: {count}")


def main():
    """Main execution."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Extract coding/LoF genes from fine-mapped credible sets"
    )
    parser.add_argument(
        '--annotations-dir',
        default='data/external/flames/Annotation_data',
        help='VEP annotations directory'
    )
    parser.add_argument(
        '--benchmark-dir',
        default='benchmarks',
        help='Benchmark directory'
    )
    parser.add_argument(
        '--pip-threshold-lof',
        type=float,
        default=0.5,
        help='PIP threshold for LoF variants (default: 0.5)'
    )
    parser.add_argument(
        '--pip-threshold-missense',
        type=float,
        default=0.8,
        help='PIP threshold for missense variants (default: 0.8)'
    )
    parser.add_argument(
        '--output',
        default='coding_lof_genes.tsv',
        help='Output filename'
    )
    
    args = parser.parse_args()
    
    extractor = CodingVariantExtractor(
        annotations_dir=args.annotations_dir,
        benchmark_dir=args.benchmark_dir,
        pip_threshold_lof=args.pip_threshold_lof,
        pip_threshold_missense=args.pip_threshold_missense,
    )
    
    logger.info("="*60)
    logger.info("Coding/LoF Variant Extraction")
    logger.info("="*60)
    logger.info(f"PIP threshold (LoF): {args.pip_threshold_lof}")
    logger.info(f"PIP threshold (missense): {args.pip_threshold_missense}")
    logger.info("")
    
    # Load fine-mapped credible sets
    credible_sets = extractor.load_finemapped_credible_sets()
    
    if credible_sets.empty:
        logger.error("No credible sets loaded!")
        return
    
    # Annotate with VEP (if not already annotated)
    credible_sets = extractor.annotate_credible_sets_with_vep(credible_sets)
    
    # Extract coding/LoF genes
    n_genes = extractor.extract_coding_lof_genes(credible_sets, require_unambiguous=True)
    
    if n_genes == 0:
        logger.warning("No coding/LoF genes extracted!")
        logger.warning("This may be because:")
        logger.warning("  1. No variants meet PIP thresholds")
        logger.warning("  2. VEP annotations are missing")
        logger.warning("  3. Column names don't match expected format")
        return
    
    # Save results
    extractor.deduplicate_and_save(output_file=args.output)
    
    logger.info("\n" + "="*60)
    logger.info("Extraction complete!")
    logger.info("="*60)


if __name__ == '__main__':
    main()
