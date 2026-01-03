#!/usr/bin/env python3
"""
Parse GTEx V8 eQTL data to extract gene-level eQTL scores.

From GTEx_Analysis_v8_eQTL.tar (1.5 GB), extracts:
- Significant variant-gene pairs across 49 tissues
- Maximum eQTL significance per gene (min p-value across tissues)
- Converts to eQTL scores: -log10(p-value)

Output: data/external/gtex/gtex_gene_eqtl_scores.tsv
Columns: gene_id, gene_symbol, best_tissue, min_pvalue, eqtl_score, n_tissues
"""

import sys
import tarfile
import gzip
from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Ensembl ID to gene symbol mapping (for tier1 genes)
# In production, would use full Ensembl GTF annotation
ENSEMBL_TO_SYMBOL = {
    'ENSG00000130203': 'APOE',
    'ENSG00000134243': 'SORT1',
    'ENSG00000148737': 'TCF7L2',
    'ENSG00000130164': 'LDLR',
    'ENSG00000169174': 'PCSK9',
    'ENSG00000084674': 'APOB',
    'ENSG00000087237': 'CETP',
    'ENSG00000113161': 'HMGCR',
    'ENSG00000175906': 'LPL',
    'ENSG00000132170': 'PPARG',
    'ENSG00000254647': 'INS',
    # Add more as needed
}


def load_ensembl_gene_mapping():
    """
    Load Ensembl ID to gene symbol mapping from GENCODE GTF.
    
    Returns dict: {ensembl_id -> gene_symbol}
    """
    logger.info("Loading Ensembl ID to gene symbol mapping from GENCODE...")
    
    base_dir = Path(__file__).parent.parent
    gtf_path = base_dir / 'data' / 'external' / 'ensembl' / 'gencode.v38.annotation.gtf.gz'
    
    if not gtf_path.exists():
        logger.warning(f"GENCODE GTF not found: {gtf_path}")
        logger.warning("Using minimal hardcoded mapping for tier1 genes")
        return ENSEMBL_TO_SYMBOL.copy()
    
    id_to_symbol = {}
    
    with gzip.open(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            # Only process gene features
            if fields[2] != 'gene':
                continue
            
            # Parse attributes
            attrs = {}
            for item in fields[8].strip(';').split(';'):
                item = item.strip()
                if not item:
                    continue
                try:
                    key, value = item.split(' ', 1)
                    attrs[key] = value.strip('"')
                except ValueError:
                    continue
            
            gene_id = attrs.get('gene_id', '').split('.')[0]  # Remove version
            gene_name = attrs.get('gene_name', '')
            gene_type = attrs.get('gene_type', '')
            
            # Only keep protein-coding genes
            if gene_type == 'protein_coding' and gene_id and gene_name:
                id_to_symbol[gene_id] = gene_name
    
    logger.info(f"Loaded {len(id_to_symbol):,} Ensembl ID mappings")
    return id_to_symbol


def parse_gtex_tar(tar_path, ensembl_to_symbol):
    """
    Parse GTEx V8 eQTL tar file and extract gene-level scores.
    
    For each tissue, reads *.signif_variant_gene_pairs.txt.gz
    Aggregates across tissues to get max eQTL significance per gene.
    """
    logger.info(f"Parsing GTEx tar: {tar_path}")
    
    # Store best (minimum) p-value per gene across all tissues
    gene_eqtl = defaultdict(lambda: {
        'min_pvalue': 1.0,
        'best_tissue': None,
        'n_tissues': 0,
        'gene_symbol': None
    })
    
    tissues_processed = 0
    total_pairs = 0
    
    with tarfile.open(tar_path, 'r') as tar:
        for member in tar.getmembers():
            # Only process significant variant-gene pairs
            if not member.name.endswith('.signif_variant_gene_pairs.txt.gz'):
                continue
            
            # Extract tissue name from filename
            # Format: GTEx_Analysis_v8_eQTL/{Tissue}.v8.signif_variant_gene_pairs.txt.gz
            parts = member.name.split('/')
            if len(parts) < 2:
                continue
            
            tissue = parts[-1].replace('.v8.signif_variant_gene_pairs.txt.gz', '')
            tissues_processed += 1
            
            logger.info(f"  Processing tissue {tissues_processed}: {tissue}...")
            
            # Extract and read the file
            f = tar.extractfile(member)
            df = pd.read_csv(f, sep='\t', compression='gzip')
            
            logger.info(f"    Significant pairs: {len(df):,}")
            total_pairs += len(df)
            
            # Process each variant-gene pair
            for _, row in df.iterrows():
                gene_id = row['gene_id'].split('.')[0]  # Remove version
                pval = row['pval_nominal']
                
                # Get gene symbol
                if gene_id not in ensembl_to_symbol:
                    continue
                
                gene_symbol = ensembl_to_symbol[gene_id]
                
                # Update if this is the best p-value for this gene
                if pval < gene_eqtl[gene_id]['min_pvalue']:
                    gene_eqtl[gene_id]['min_pvalue'] = pval
                    gene_eqtl[gene_id]['best_tissue'] = tissue
                    gene_eqtl[gene_id]['gene_symbol'] = gene_symbol
                
                if gene_eqtl[gene_id]['gene_symbol'] is None:
                    gene_eqtl[gene_id]['gene_symbol'] = gene_symbol
                    gene_eqtl[gene_id]['n_tissues'] += 1
    
    logger.info(f"Processed {tissues_processed} tissues")
    logger.info(f"Total significant pairs: {total_pairs:,}")
    logger.info(f"Unique genes with eQTLs: {len(gene_eqtl):,}")
    
    return dict(gene_eqtl)


def convert_to_dataframe(gene_eqtl):
    """
    Convert gene eQTL dict to pandas DataFrame.
    
    Adds eQTL score = -log10(min_pvalue)
    """
    logger.info("Converting to DataFrame...")
    
    rows = []
    for gene_id, info in gene_eqtl.items():
        # Calculate eQTL score
        min_pval = info['min_pvalue']
        eqtl_score = -np.log10(min_pval + 1e-300)  # Avoid log(0)
        
        rows.append({
            'gene_id': gene_id,
            'gene_symbol': info['gene_symbol'],
            'best_tissue': info['best_tissue'],
            'min_pvalue': min_pval,
            'eqtl_score': eqtl_score,
            'n_tissues': info['n_tissues']
        })
    
    df = pd.DataFrame(rows)
    
    # Sort by eQTL score (descending)
    df = df.sort_values('eqtl_score', ascending=False)
    
    logger.info(f"Created DataFrame with {len(df):,} genes")
    logger.info(f"  eQTL score range: {df['eqtl_score'].min():.2f} - {df['eqtl_score'].max():.2f}")
    logger.info(f"  Median eQTL score: {df['eqtl_score'].median():.2f}")
    
    return df


def main():
    """Main entry point."""
    print("=" * 80)
    print("GTEx V8 eQTL Gene Score Extraction")
    print("=" * 80)
    print()
    
    # Define paths
    base_dir = Path(__file__).parent.parent
    tar_path = base_dir / 'data' / 'external' / 'gtex' / 'GTEx_Analysis_v8_eQTL.tar'
    output_path = base_dir / 'data' / 'external' / 'gtex' / 'gtex_gene_eqtl_scores.tsv'
    
    # Check if tar exists
    if not tar_path.exists():
        print(f"ERROR: GTEx tar not found: {tar_path}")
        print("Please download first:")
        print("  curl -L -o GTEx_Analysis_v8_eQTL.tar https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar")
        sys.exit(1)
    
    # Check tar file size
    tar_size_gb = tar_path.stat().st_size / 1e9
    print(f"GTEx tar file: {tar_path}")
    print(f"Size: {tar_size_gb:.2f} GB")
    
    if tar_size_gb < 1.0:
        print(f"ERROR: Tar file too small ({tar_size_gb:.2f} GB)")
        print("Expected ~1.5 GB, likely incomplete download")
        print("Please re-download GTEx data")
        sys.exit(1)
    
    print()
    
    # Load Ensembl ID mapping
    ensembl_to_symbol = load_ensembl_gene_mapping()
    print()
    
    # Parse GTEx tar
    gene_eqtl = parse_gtex_tar(tar_path, ensembl_to_symbol)
    print()
    
    # Convert to DataFrame
    df = convert_to_dataframe(gene_eqtl)
    print()
    
    # Write output
    print(f"Writing results to {output_path}...")
    df.to_csv(output_path, sep='\t', index=False)
    
    print(f"Wrote {len(df):,} genes")
    print()
    print("SUCCESS: GTEx eQTL parsing complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
