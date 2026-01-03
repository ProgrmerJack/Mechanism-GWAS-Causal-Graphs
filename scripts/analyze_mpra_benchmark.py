#!/usr/bin/env python3
"""
Analyze Abell et al. 2022 MPRA data for RegulatoryBench v3 construction.

This script processes the MPRA results from:
Abell et al. "Multiple Causal Variants Underlie Genetic Associations in Humans"
Science 2022 - 49,256 variants tested in lymphoblastoid cells (LCLs)

Key output: Variants with significant allelic effects that can serve as 
experimental ground truth for variant-to-gene benchmark.
"""

import pandas as pd
import os
from pathlib import Path

# Define paths
BASE_DIR = Path("c:/Users/Jack0/GitHub/Mechanism-GWAS-Causal-Graphs")
MPRA_DIR = BASE_DIR / "data/external/mpra_abell_2022/nsabell-mpra-v2-d62fca8"
OUTPUT_DIR = BASE_DIR / "data/processed"

def load_mpra_results():
    """Load the main MPRA summary statistics."""
    sumstats_file = MPRA_DIR / "sumstats/1KG_novaSeq_DESeq2_Love_Base.txt"
    
    print(f"Loading MPRA results from {sumstats_file}")
    df = pd.read_csv(sumstats_file, sep='\t')
    print(f"  Total variants tested: {len(df):,}")
    print(f"  Columns: {list(df.columns)}")
    
    return df

def load_eqtl_mapping():
    """Load the eQTL variant-gene mapping (Supplementary Table S1)."""
    eqtl_file = MPRA_DIR / "supptables/SupplementaryTableS1.txt"
    
    print(f"\nLoading eQTL mapping from {eqtl_file}")
    df = pd.read_csv(eqtl_file, sep='\t')
    print(f"  Total eQTL-gene pairs: {len(df):,}")
    print(f"  Unique variants: {df['SNP'].nunique():,}")
    print(f"  Unique genes: {df['gene'].str.split('_').str[0].nunique():,}")
    
    return df

def analyze_mpra_hits(df, allele_padj_threshold=0.05, expr_padj_threshold=0.05):
    """Identify significant MPRA hits."""
    print(f"\n=== MPRA Hit Analysis ===")
    print(f"Thresholds: allele padj < {allele_padj_threshold}, expression padj < {expr_padj_threshold}")
    
    # Filter for variants with eQTL gene annotations (not controls)
    eqtl_variants = df[df['geneID'].notna() & (df['geneID'] != 'NA')].copy()
    print(f"\nVariants with eQTL gene annotation: {len(eqtl_variants):,}")
    
    # Expression effects (allele-independent)
    expr_hits = eqtl_variants[eqtl_variants['padj_expr'] < expr_padj_threshold]
    print(f"Expression hits (padj < {expr_padj_threshold}): {len(expr_hits):,}")
    
    # Allelic effects (variant-specific regulatory effects)
    allele_hits = eqtl_variants[eqtl_variants['padj_allele'] < allele_padj_threshold]
    print(f"Allelic hits (padj < {allele_padj_threshold}): {len(allele_hits):,}")
    
    # Combined: both expression AND allelic effects
    both_hits = eqtl_variants[
        (eqtl_variants['padj_expr'] < expr_padj_threshold) & 
        (eqtl_variants['padj_allele'] < allele_padj_threshold)
    ]
    print(f"Both expression AND allelic hits: {len(both_hits):,}")
    
    # Strong allelic hits (more stringent threshold)
    strong_allele = eqtl_variants[eqtl_variants['padj_allele'] < 0.01]
    print(f"Strong allelic hits (padj < 0.01): {len(strong_allele):,}")
    
    very_strong = eqtl_variants[eqtl_variants['padj_allele'] < 0.001]
    print(f"Very strong allelic hits (padj < 0.001): {len(very_strong):,}")
    
    return {
        'all_eqtl': eqtl_variants,
        'expr_hits': expr_hits,
        'allele_hits': allele_hits,
        'both_hits': both_hits,
        'strong_allele': strong_allele
    }

def extract_benchmark_variants(hits_dict, eqtl_df):
    """Extract variants suitable for variant-to-gene benchmark."""
    print("\n=== Extracting Benchmark Variants ===")
    
    allele_hits = hits_dict['allele_hits'].copy()
    
    # Parse chromosome and position
    allele_hits['chr_pos'] = allele_hits['chrom'].str.replace('chr', '') + '_' + allele_hits['pos'].astype(str)
    
    # Extract gene information
    allele_hits['ensembl_id'] = allele_hits['geneID'].str.split('_').str[0]
    allele_hits['gene_symbol'] = allele_hits['geneID'].str.split('_').str[1]
    
    # Merge with eQTL data to get effect direction
    # Create matching key for eQTL data
    eqtl_df['chr_pos'] = eqtl_df['SNP']
    eqtl_df['ensembl_id'] = eqtl_df['gene'].str.split('_').str[0]
    
    # Key statistics
    print(f"\nAllelic hits with gene annotations: {len(allele_hits):,}")
    print(f"  Unique chromosomal positions: {allele_hits['chr_pos'].nunique():,}")
    print(f"  Unique genes: {allele_hits['ensembl_id'].nunique():,}")
    
    # Effect size distribution
    print(f"\nAllelic effect size distribution:")
    print(f"  Min: {allele_hits['log2FoldChange_allele'].min():.3f}")
    print(f"  Max: {allele_hits['log2FoldChange_allele'].max():.3f}")
    print(f"  Mean: {allele_hits['log2FoldChange_allele'].mean():.3f}")
    print(f"  Median: {allele_hits['log2FoldChange_allele'].median():.3f}")
    
    # Direction of effects
    positive = (allele_hits['log2FoldChange_allele'] > 0).sum()
    negative = (allele_hits['log2FoldChange_allele'] < 0).sum()
    print(f"\nEffect direction:")
    print(f"  Positive (alt > ref): {positive:,}")
    print(f"  Negative (alt < ref): {negative:,}")
    
    return allele_hits

def create_variant_gene_benchmark(allele_hits):
    """Create the final variant-to-gene benchmark file."""
    print("\n=== Creating Variant-Gene Benchmark ===")
    
    # Select key columns for benchmark
    benchmark_cols = [
        'VarID', 'chrom', 'pos', 'ref', 'alt',
        'geneID', 'ensembl_id', 'gene_symbol',
        'log2FoldChange_allele', 'padj_allele',
        'log2FoldChange_expr', 'padj_expr'
    ]
    
    # Handle missing columns
    available_cols = [c for c in benchmark_cols if c in allele_hits.columns]
    benchmark = allele_hits[available_cols].copy()
    
    # Add metadata
    benchmark['source'] = 'Abell_Science_2022'
    benchmark['cell_type'] = 'LCL'
    benchmark['assay'] = 'MPRA'
    
    # Save to file
    output_file = OUTPUT_DIR / "mpra_abell_benchmark.tsv"
    benchmark.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved benchmark file: {output_file}")
    print(f"  Total variant-gene pairs: {len(benchmark):,}")
    
    return benchmark

def analyze_gene_coverage(allele_hits, eqtl_df):
    """Analyze which genes have MPRA-validated variants."""
    print("\n=== Gene Coverage Analysis ===")
    
    # Get unique genes from MPRA hits
    mpra_genes = allele_hits['ensembl_id'].unique()
    print(f"Genes with MPRA-validated variants: {len(mpra_genes):,}")
    
    # Number of validated variants per gene
    variants_per_gene = allele_hits.groupby('ensembl_id').size()
    print(f"\nVariants per gene:")
    print(f"  Mean: {variants_per_gene.mean():.1f}")
    print(f"  Median: {variants_per_gene.median():.0f}")
    print(f"  Max: {variants_per_gene.max()}")
    print(f"  Genes with 1 variant: {(variants_per_gene == 1).sum()}")
    print(f"  Genes with 2+ variants: {(variants_per_gene >= 2).sum()}")
    print(f"  Genes with 5+ variants: {(variants_per_gene >= 5).sum()}")
    
    # Top genes by number of variants
    print(f"\nTop 10 genes by number of validated variants:")
    top_genes = variants_per_gene.nlargest(10)
    for gene, count in top_genes.items():
        gene_symbol = allele_hits[allele_hits['ensembl_id'] == gene]['gene_symbol'].iloc[0] if 'gene_symbol' in allele_hits.columns else 'NA'
        print(f"  {gene} ({gene_symbol}): {count} variants")
    
    return mpra_genes

def main():
    print("=" * 60)
    print("MPRA Benchmark Analysis - Abell et al. 2022")
    print("=" * 60)
    
    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load data
    mpra_df = load_mpra_results()
    eqtl_df = load_eqtl_mapping()
    
    # Analyze hits at different thresholds
    hits_dict = analyze_mpra_hits(mpra_df)
    
    # Extract benchmark variants
    allele_hits = extract_benchmark_variants(hits_dict, eqtl_df)
    
    # Create variant-gene benchmark
    benchmark = create_variant_gene_benchmark(allele_hits)
    
    # Analyze gene coverage
    mpra_genes = analyze_gene_coverage(allele_hits, eqtl_df)
    
    # Summary statistics for paper
    print("\n" + "=" * 60)
    print("SUMMARY FOR REGULATORYBENCH V3")
    print("=" * 60)
    print(f"""
MPRA Dataset: Abell et al. Science 2022
- Total variants tested: {len(mpra_df):,}
- Variants with significant allelic effects (FDR < 0.05): {len(hits_dict['allele_hits']):,}
- Unique genes with MPRA-validated variants: {len(mpra_genes):,}
- Cell type: Lymphoblastoid cells (LCLs)
- Assay: Massively Parallel Reporter Assay (MPRA)
- Ground truth type: Allelic regulatory effects (ref vs alt)

This provides variant→gene pairs where:
1. The variant was experimentally tested in an MPRA
2. The variant shows significant allele-specific regulatory activity
3. The target gene is the eQTL gene (expression QTL)

Key advantage: These are DIRECT experimental measurements of variant
regulatory effects, not inferred from population genetics.
""")
    
    return benchmark, mpra_genes

if __name__ == "__main__":
    benchmark, genes = main()
