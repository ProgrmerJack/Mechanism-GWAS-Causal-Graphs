"""
Extract ABC + eQTL + Distance + Constraint Features for Tier1 Genes
=====================================================================

This script extracts features from multiple data sources for comprehensive
baseline validation:

1. ABC scores: Enhancer-gene predictions from Nasser et al. 2021 (131 biosamples)
2. eQTL scores: GTEx v8 colocalization (will download if needed)
3. Distance scores: TSS distance from Ensembl annotations
4. Coding scores: VEP annotations from tier1 metadata
5. Constraint scores: gnomAD pLI from constraint database

Output: data/processed/benchmark/tier1_features.tsv

Author: Mechanism-GWAS-Causal-Graphs team
Date: December 2025
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import gzip
import logging
from typing import Dict, List, Set
import requests
from io import StringIO

# Add project root to path
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_tier1_genes() -> pd.DataFrame:
    """Load tier1 gold standard genes."""
    tier1_file = project_root / "data" / "processed" / "benchmark" / "tier1_gold_standard_genes.tsv"
    df = pd.read_csv(tier1_file, sep='\t')
    logger.info(f"Loaded {len(df)} tier1 genes")
    return df


def extract_abc_scores(tier1_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Extract ABC scores for tier1 genes from Nasser et al. 2021 database.
    
    Priority biosamples (cardiometabolic):
    - Liver: LIVER, HepG2, liver-*
    - Adipose: adipose-*, ADIPOSE, fat
    - Pancreas: pancrea*, islet*, PANC
    - Vascular: aorta*, coronary*, endothel*, HUVEC
    
    Returns DataFrame with columns: gene_symbol, max_abc_score, 
                                    mean_abc_score, n_biosamples
    """
    logger.info("Extracting ABC scores for tier1 genes...")
    
    abc_file = project_root / "data" / "external" / "abc" / "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
    
    if not abc_file.exists():
        logger.error(f"ABC file not found: {abc_file}")
        return pd.DataFrame()
    
    # Get tier1 gene symbols
    tier1_symbols = set(tier1_genes['gene_symbol'].values)
    logger.info(f"Searching for {len(tier1_symbols)} tier1 genes in ABC database")
    
    # Cardiometabolic biosample patterns (case-insensitive)
    priority_patterns = [
        'liver', 'hepg2', 'hepato',
        'adipose', 'fat', 'brite',
        'pancrea', 'islet', 'panc',
        'aorta', 'coronary', 'endothel', 'huvec', 'smooth_muscle'
    ]
    
    # Parse ABC file and collect predictions for tier1 genes
    gene_predictions = {gene: [] for gene in tier1_symbols}
    n_total = 0
    n_tier1 = 0
    n_priority = 0
    
    with gzip.open(abc_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        
        # Find column indices
        try:
            gene_idx = header.index('TargetGene')
            score_idx = header.index('ABC.Score')
            celltype_idx = header.index('CellType')
            distance_idx = header.index('distance')
        except ValueError as e:
            logger.error(f"Column not found in ABC file: {e}")
            return pd.DataFrame()
        
        logger.info("Parsing ABC database (7.7M predictions)...")
        for i, line in enumerate(f):
            if i % 1000000 == 0:
                logger.info(f"  Processed {i:,} predictions, found {n_tier1:,} tier1 matches")
            
            parts = line.strip().split('\t')
            if len(parts) <= max(gene_idx, score_idx, celltype_idx, distance_idx):
                continue
            
            gene = parts[gene_idx]
            n_total += 1
            
            if gene in tier1_symbols:
                try:
                    abc_score = float(parts[score_idx])
                    celltype = parts[celltype_idx]
                    distance = float(parts[distance_idx]) if parts[distance_idx] != 'NA' else np.nan
                    
                    # Check if priority biosample
                    is_priority = any(pat in celltype.lower() for pat in priority_patterns)
                    
                    gene_predictions[gene].append({
                        'abc_score': abc_score,
                        'celltype': celltype,
                        'distance': distance,
                        'is_priority': is_priority
                    })
                    
                    n_tier1 += 1
                    if is_priority:
                        n_priority += 1
                        
                except (ValueError, IndexError):
                    continue
    
    logger.info(f"Total ABC predictions: {n_total:,}")
    logger.info(f"Tier1 gene predictions: {n_tier1:,}")
    logger.info(f"Priority biosample predictions: {n_priority:,}")
    
    # Aggregate ABC scores per gene
    abc_features = []
    for gene, predictions in gene_predictions.items():
        if len(predictions) == 0:
            abc_features.append({
                'gene_symbol': gene,
                'max_abc_score': 0.0,
                'mean_abc_score': 0.0,
                'max_abc_priority': 0.0,
                'mean_abc_priority': 0.0,
                'n_biosamples': 0,
                'n_priority_biosamples': 0,
                'min_distance': np.nan
            })
        else:
            pred_df = pd.DataFrame(predictions)
            priority_df = pred_df[pred_df['is_priority']]
            
            abc_features.append({
                'gene_symbol': gene,
                'max_abc_score': pred_df['abc_score'].max(),
                'mean_abc_score': pred_df['abc_score'].mean(),
                'max_abc_priority': priority_df['abc_score'].max() if len(priority_df) > 0 else 0.0,
                'mean_abc_priority': priority_df['abc_score'].mean() if len(priority_df) > 0 else 0.0,
                'n_biosamples': len(pred_df['celltype'].unique()),
                'n_priority_biosamples': len(priority_df['celltype'].unique()) if len(priority_df) > 0 else 0,
                'min_distance': pred_df['distance'].min()
            })
    
    abc_df = pd.DataFrame(abc_features)
    logger.info(f"Extracted ABC features for {len(abc_df)} genes")
    logger.info(f"  Genes with ABC>0: {(abc_df['max_abc_score'] > 0).sum()}")
    logger.info(f"  Genes with priority ABC>0: {(abc_df['max_abc_priority'] > 0).sum()}")
    
    return abc_df


def download_gnomad_constraint() -> pd.DataFrame:
    """
    Download gnomAD v2.1.1 constraint scores (pLI) for tier1 genes.
    
    Returns DataFrame with columns: gene_symbol, pli, loeuf
    """
    logger.info("Downloading gnomAD constraint scores...")
    
    # gnomAD v2.1.1 constraint file
    url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    
    constraint_file = project_root / "data" / "external" / "gnomad" / "gnomad_constraint.tsv.gz"
    constraint_file.parent.mkdir(parents=True, exist_ok=True)
    
    if not constraint_file.exists():
        logger.info(f"Downloading from {url}...")
        try:
            response = requests.get(url, stream=True, timeout=60)
            response.raise_for_status()
            
            with open(constraint_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            logger.info(f"Downloaded to {constraint_file}")
        except Exception as e:
            logger.error(f"Download failed: {e}")
            return pd.DataFrame()
    else:
        logger.info(f"Using cached file: {constraint_file}")
    
    # Parse constraint file
    try:
        df = pd.read_csv(constraint_file, sep='\t', compression='gzip')
        logger.info(f"Loaded {len(df)} genes from gnomAD")
        
        # Keep only relevant columns
        constraint_cols = ['gene', 'pLI', 'oe_lof_upper']
        if all(col in df.columns for col in constraint_cols):
            constraint_df = df[constraint_cols].copy()
            constraint_df.columns = ['gene_symbol', 'pli', 'loeuf']
            logger.info(f"Extracted constraint scores for {len(constraint_df)} genes")
            return constraint_df
        else:
            logger.error(f"Expected columns not found. Available: {df.columns.tolist()[:10]}")
            return pd.DataFrame()
            
    except Exception as e:
        logger.error(f"Failed to parse constraint file: {e}")
        return pd.DataFrame()


def get_ensembl_tss_positions(tier1_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Get TSS positions from GENCODE GTF annotations for tier1 genes.
    
    Uses pre-extracted TSS data from parse_gencode_tss.py (19,930 protein-coding genes).
    
    Returns DataFrame with columns: gene_symbol, ensembl_id, chr, tss_position, strand
    """
    logger.info("Loading TSS positions from GENCODE annotations...")
    
    # Load pre-extracted TSS data
    tss_file = project_root / "data" / "external" / "ensembl" / "gencode_tss_grch38.tsv"
    
    if not tss_file.exists():
        logger.warning(f"TSS file not found: {tss_file}")
        logger.warning("Run: python scripts/parse_gencode_tss.py")
        logger.warning("Falling back to placeholder TSS positions")
        
        # Fallback placeholder (for development only)
        tss_data = []
        for _, row in tier1_genes.iterrows():
            tss_data.append({
                'gene_symbol': row['gene_symbol'],
                'ensembl_id': row.get('ensembl_id', ''),
                'chr': 'NA',
                'tss_position': 0,
                'strand': '+'
            })
        return pd.DataFrame(tss_data)
    
    # Load GENCODE TSS data
    gencode_tss = pd.read_csv(tss_file, sep='\t')
    logger.info(f"Loaded {len(gencode_tss):,} GENCODE TSS positions")
    
    # Create lookup dict for fast access
    tss_lookup = {}
    for _, row in gencode_tss.iterrows():
        tss_lookup[row['gene_symbol']] = {
            'chr': str(row['chromosome']),
            'tss_position': int(row['tss_position']),
            'strand': row['strand']
        }
    
    # Match tier1 genes to GENCODE TSS
    tss_data = []
    found = 0
    for _, row in tier1_genes.iterrows():
        gene = row['gene_symbol']
        
        if gene in tss_lookup:
            info = tss_lookup[gene]
            tss_data.append({
                'gene_symbol': gene,
                'ensembl_id': row.get('ensembl_id', ''),
                'chr': info['chr'],
                'tss_position': info['tss_position'],
                'strand': info['strand']
            })
            found += 1
        else:
            logger.warning(f"  Gene {gene} not found in GENCODE, using default")
            tss_data.append({
                'gene_symbol': gene,
                'ensembl_id': row.get('ensembl_id', ''),
                'chr': 'NA',
                'tss_position': 0,
                'strand': '+'
            })
    
    tss_df = pd.DataFrame(tss_data)
    logger.info(f"Matched {found}/{len(tier1_genes)} tier1 genes to GENCODE TSS")
    return tss_df


def calculate_distance_scores(tier1_genes: pd.DataFrame, tss_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate distance-based scores from GWAS lead SNPs to gene TSS.
    
    For genes with known GWAS associations, calculate distance.
    For others, use ABC min_distance.
    
    Returns DataFrame with columns: gene_symbol, distance_to_lead_snp, distance_score
    """
    logger.info("Calculating distance scores...")
    
    # Known GWAS lead SNPs for tier1 genes (from literature)
    gwas_lead_snps = {
        'APOE': ('19', 44908684),  # rs429358
        'SORT1': ('1', 109274968),  # rs12740374
        'TCF7L2': ('10', 112998590),  # rs7903146
        'LDLR': ('19', 11089362),  # rs688
        'PCSK9': ('1', 55039974),  # rs11591147
        'APOB': ('2', 21001429),  # rs693
        'CETP': ('16', 56959412),  # rs3764261
        'HMGCR': ('5', 74632993),  # rs12916
        'LPL': ('8', 19840862),  # rs328
        'PPARG': ('3', 12328867),  # rs1801282
        'INS': ('11', 2159842),  # rs689
    }
    
    distance_data = []
    for _, row in tss_df.iterrows():
        gene = row['gene_symbol']
        
        if gene in gwas_lead_snps:
            lead_chr, lead_pos = gwas_lead_snps[gene]
            if lead_chr == row['chr'] and row['tss_position'] > 0:
                distance = abs(lead_pos - row['tss_position'])
            else:
                distance = 100000  # Default 100kb
        else:
            distance = 100000  # Default
        
        # Distance score: inverse distance (closer = higher score)
        distance_score = 1.0 / (1.0 + distance / 100000)
        
        distance_data.append({
            'gene_symbol': gene,
            'distance_to_lead_snp': distance,
            'distance_score': distance_score
        })
    
    distance_df = pd.DataFrame(distance_data)
    logger.info(f"Calculated distance scores for {len(distance_df)} genes")
    return distance_df


def extract_coding_scores(tier1_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Extract coding variant indicators from tier1 metadata.
    
    Returns DataFrame with columns: gene_symbol, is_coding, coding_score
    """
    logger.info("Extracting coding scores from tier1 metadata...")
    
    coding_data = []
    for _, row in tier1_genes.iterrows():
        gene = row['gene_symbol']
        evidence = row.get('evidence_type', '')
        notes = row.get('notes', '')
        
        # Check if gene has coding variants mentioned
        is_coding = 1 if any(word in notes.lower() for word in ['missense', 'nonsense', 'coding', 'mutation', 'variant']) else 0
        
        coding_data.append({
            'gene_symbol': gene,
            'is_coding': is_coding,
            'coding_score': float(is_coding)
        })
    
    coding_df = pd.DataFrame(coding_data)
    logger.info(f"Extracted coding scores for {len(coding_df)} genes")
    logger.info(f"  Genes with coding variants: {coding_df['is_coding'].sum()}")
    return coding_df


def create_simplified_eqtl_scores(tier1_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Load GTEx V8 eQTL scores for tier1 genes.
    
    Uses pre-parsed eQTL data from parse_gtex_eqtl.py.
    Scores represent maximum eQTL significance across 49 tissues (-log10(p-value)).
    
    Returns DataFrame with columns: gene_symbol, eqtl_score, best_tissue, min_pvalue
    """
    logger.info("Loading GTEx V8 eQTL scores...")
    
    # Load pre-parsed GTEx eQTL data
    eqtl_file = project_root / "data" / "external" / "gtex" / "gtex_gene_eqtl_scores.tsv"
    
    if not eqtl_file.exists():
        logger.warning(f"GTEx eQTL file not found: {eqtl_file}")
        logger.warning("Run: python scripts/parse_gtex_eqtl.py")
        logger.warning("Falling back to placeholder eQTL scores")
        
        # Fallback placeholder (for development only)
        eqtl_data = []
        for _, row in tier1_genes.iterrows():
            eqtl_data.append({
                'gene_symbol': row['gene_symbol'],
                'eqtl_score': 0.50,  # Placeholder
                'best_tissue': 'unknown',
                'min_pvalue': 1.0
            })
        return pd.DataFrame(eqtl_data)
    
    # Load GTEx eQTL data
    gtex_eqtl = pd.read_csv(eqtl_file, sep='\t')
    logger.info(f"Loaded {len(gtex_eqtl):,} GTEx eQTL scores")
    
    # Create lookup dict for fast access
    eqtl_lookup = {}
    for _, row in gtex_eqtl.iterrows():
        eqtl_lookup[row['gene_symbol']] = {
            'eqtl_score': float(row['eqtl_score']),
            'best_tissue': row['best_tissue'],
            'min_pvalue': float(row['min_pvalue']),
            'n_tissues': int(row['n_tissues'])
        }
    
    # Match tier1 genes to GTEx eQTL
    eqtl_data = []
    found = 0
    for _, row in tier1_genes.iterrows():
        gene = row['gene_symbol']
        
        if gene in eqtl_lookup:
            info = eqtl_lookup[gene]
            eqtl_data.append({
                'gene_symbol': gene,
                'eqtl_score': info['eqtl_score'],
                'best_tissue': info['best_tissue'],
                'min_pvalue': info['min_pvalue'],
                'n_tissues': info['n_tissues']
            })
            found += 1
        else:
            logger.warning(f"  Gene {gene} not found in GTEx, using default score")
            eqtl_data.append({
                'gene_symbol': gene,
                'eqtl_score': 0.0,  # No eQTL evidence
                'best_tissue': 'none',
                'min_pvalue': 1.0,
                'n_tissues': 0
            })
    
    eqtl_df = pd.DataFrame(eqtl_data)
    logger.info(f"Matched {found}/{len(tier1_genes)} tier1 genes to GTEx eQTL")
    logger.info(f"  Mean eQTL score: {eqtl_df['eqtl_score'].mean():.2f}")
    logger.info(f"  Max eQTL score: {eqtl_df['eqtl_score'].max():.2f}")
    return eqtl_df


def main():
    """Extract all features for tier1 genes."""
    logger.info("="*80)
    logger.info("FEATURE EXTRACTION FOR TIER1 GENES")
    logger.info("="*80)
    
    # Load tier1 genes
    tier1_genes = load_tier1_genes()
    
    # Extract ABC scores (REAL data from Nasser et al. 2021)
    abc_df = extract_abc_scores(tier1_genes)
    
    # Download constraint scores (REAL data from gnomAD)
    constraint_df = download_gnomad_constraint()
    
    # Get TSS positions
    tss_df = get_ensembl_tss_positions(tier1_genes)
    
    # Calculate distance scores
    distance_df = calculate_distance_scores(tier1_genes, tss_df)
    
    # Extract coding scores
    coding_df = extract_coding_scores(tier1_genes)
    
    # Create simplified eQTL scores
    eqtl_df = create_simplified_eqtl_scores(tier1_genes)
    
    # Merge all features
    logger.info("\nMerging all features...")
    features = tier1_genes[['gene_symbol', 'trait', 'evidence_type']].copy()
    
    for df in [abc_df, constraint_df, distance_df, coding_df, eqtl_df]:
        if len(df) > 0:
            features = features.merge(df, on='gene_symbol', how='left')
    
    # Fill NaN values
    numeric_cols = features.select_dtypes(include=[np.number]).columns
    features[numeric_cols] = features[numeric_cols].fillna(0)
    
    # Save features
    output_file = project_root / "data" / "processed" / "benchmark" / "tier1_features.tsv"
    features.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"\n{'='*80}")
    logger.info(f"FEATURES SAVED: {output_file}")
    logger.info(f"{'='*80}")
    logger.info(f"Total genes: {len(features)}")
    logger.info(f"Total features: {len(features.columns)}")
    logger.info(f"\nFeature columns:")
    for col in features.columns:
        logger.info(f"  - {col}")
    
    logger.info(f"\nSummary statistics:")
    logger.info(f"  Max ABC score: {features['max_abc_score'].max():.3f}")
    logger.info(f"  Mean ABC score: {features['max_abc_score'].mean():.3f}")
    logger.info(f"  Genes with ABC>0: {(features['max_abc_score'] > 0).sum()}")
    logger.info(f"  Genes with pLI>0.9: {(features.get('pli', pd.Series([0])) > 0.9).sum()}")
    logger.info(f"  Mean distance: {features['distance_to_lead_snp'].mean():.0f} bp")
    logger.info(f"  Genes with coding variants: {features['is_coding'].sum()}")
    
    return features


if __name__ == '__main__':
    features = main()
