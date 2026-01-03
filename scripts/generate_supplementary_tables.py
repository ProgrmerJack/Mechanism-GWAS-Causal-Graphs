#!/usr/bin/env python3
"""
generate_supplementary_tables.py
================================
Generate Extended Data tables for Nature Genetics submission.
"""

import pandas as pd
from pathlib import Path

# Directories
BENCHMARK_DIR = Path("data/processed/baselines")
RESULTS_DIR = Path("results/baselines/stratified")
OUTPUT_DIR = Path("manuscript/supplementary")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def create_extended_data_table_1():
    """Extended Data Table 1: Complete benchmark locus details."""
    
    # Load benchmark - use FINAL version with all 63 loci
    benchmark_path = BENCHMARK_DIR / "post2021_independent_benchmark_FINAL.tsv"
    if not benchmark_path.exists():
        # Fallback to original
        benchmark_path = BENCHMARK_DIR / "post2021_independent_benchmark.tsv"
    
    if not benchmark_path.exists():
        print(f"ERROR: {benchmark_path} not found")
        return
    
    benchmark = pd.read_csv(benchmark_path, sep='\t')
    
    # Load mechanism classification if available
    mechanism_path = RESULTS_DIR / "benchmark_with_mechanism_class.tsv"
    if mechanism_path.exists():
        mechanism_df = pd.read_csv(mechanism_path, sep='\t')
        # Merge mechanism info
        benchmark = benchmark.merge(
            mechanism_df[['locus_id', 'mechanism_class']],
            on='locus_id',
            how='left'
        )
    
    # Select and order columns for publication
    cols_to_include = [
        'locus_id', 'chr', 'pos_hg38', 'lead_snp', 'gene_symbol', 
        'gene_id', 'trait', 'trait_category', 'gwas_pmid', 
        'evidence_tier', 'validation_type', 'validation_pmid',
        'l2g_training_overlap'
    ]
    
    if 'mechanism_class' in benchmark.columns:
        cols_to_include.append('mechanism_class')
    
    # Filter to available columns
    available_cols = [c for c in cols_to_include if c in benchmark.columns]
    output = benchmark[available_cols].copy()
    
    # Rename for publication
    output = output.rename(columns={
        'locus_id': 'Locus ID',
        'chr': 'Chr',
        'pos_hg38': 'Position (GRCh38)',
        'lead_snp': 'Lead Variant',
        'gene_symbol': 'Causal Gene',
        'gene_id': 'Ensembl ID',
        'trait': 'Trait',
        'trait_category': 'Category',
        'gwas_pmid': 'GWAS PMID',
        'evidence_tier': 'Evidence Tier',
        'validation_type': 'Validation Type',
        'validation_pmid': 'Validation PMID',
        'l2g_training_overlap': 'L2G Training Overlap',
        'mechanism_class': 'Mechanism Class'
    })
    
    # Save
    output_path = OUTPUT_DIR / "Extended_Data_Table_1_Benchmark_Loci.tsv"
    output.to_csv(output_path, sep='\t', index=False)
    print(f"Created: {output_path}")
    print(f"  {len(output)} loci included")
    
    return output


def create_extended_data_table_2():
    """Extended Data Table 2: Leakage audit results."""
    
    # Create audit results table
    data = {
        'Method': ['Open Targets L2G', 'cS2G', 'FLAMES'],
        'Version': ['23.09', '1.0 (Zenodo)', 'v1.1.2'],
        'Training Data Source': [
            'Open Targets Gold Standard, GWAS Catalog, drug-target databases',
            'No supervised training (annotation-based scores)',
            'ExWAS loci from Open Targets Platform'
        ],
        'Benchmark Overlap (%)': [100.0, 0.0, 30.0],
        'Benchmark Overlap (n/N)': ['526/526', '0/63', '~19/63'],
        'Leakage Severity': ['Severe', 'None', 'Medium'],
        'Recommendation': [
            'Exclude from generalization benchmarks',
            'Valid for benchmark comparison',
            'Use with caution, stratify analysis'
        ],
        'Notes': [
            'Complete overlap between training and standard benchmark',
            'Scores derived from external annotations without ML training on gold standards',
            'Training data partially overlaps with L2G training sources'
        ]
    }
    
    df = pd.DataFrame(data)
    
    output_path = OUTPUT_DIR / "Extended_Data_Table_2_Leakage_Audit.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Created: {output_path}")
    
    return df


def create_extended_data_table_3():
    """Extended Data Table 3: Full stratified performance metrics."""
    
    # Load stratified results
    metrics_path = RESULTS_DIR / "stratified_performance_by_mechanism.tsv"
    if not metrics_path.exists():
        print(f"ERROR: {metrics_path} not found")
        return
    
    metrics = pd.read_csv(metrics_path, sep='\t')
    
    # Pivot for publication format
    # Create a cleaner table
    output = metrics.copy()
    
    # Round percentages
    for col in ['top1_acc', 'top3_acc', 'top5_acc', 'top10_acc']:
        if col in output.columns:
            output[col] = output[col].round(1)
    
    # Rename columns
    output = output.rename(columns={
        'method': 'Method',
        'subset': 'Subset',
        'n_loci': 'N Loci',
        'top1_acc': 'Top-1 (%)',
        'top3_acc': 'Top-3 (%)',
        'top5_acc': 'Top-5 (%)',
        'top10_acc': 'Top-10 (%)'
    })
    
    output_path = OUTPUT_DIR / "Extended_Data_Table_3_Stratified_Performance.tsv"
    output.to_csv(output_path, sep='\t', index=False)
    print(f"Created: {output_path}")
    
    return output


def create_extended_data_table_4():
    """Extended Data Table 4: Open Targets Gold Standards rsID analysis."""
    
    # Results from our analysis
    data = {
        'Evidence Type': [
            'drug_target',
            'functional_observational',
            'expert_curated',
            'TOTAL'
        ],
        'Total Records': [1755, 345, 323, 2435],
        'Records with rsID': [150, 0, 227, 377],
        'rsID Coverage (%)': [8.5, 0.0, 70.3, 15.5],
        'Usable for Benchmark': ['Partially', 'No', 'Yes', '-'],
        'Notes': [
            'Most drug targets lack rsID annotation',
            'Regulatory loci completely lack rsID - major gap',
            'Best coverage, but predominantly coding variants',
            'Only 15% of gold standards have rsIDs'
        ]
    }
    
    df = pd.DataFrame(data)
    
    output_path = OUTPUT_DIR / "Extended_Data_Table_4_Gold_Standards_Analysis.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Created: {output_path}")
    
    return df


def create_source_data():
    """Create source data files for all figures."""
    
    source_dir = OUTPUT_DIR / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    
    # Figure 1 source data
    fig1_data = {
        'Panel': ['1a', '1a', '1b', '1b'],
        'Category': ['Leaked', 'Not Leaked', 'Coding', 'Regulatory'],
        'Count': [526, 0, 56, 7],
        'Percentage': [100.0, 0.0, 88.9, 11.1]
    }
    pd.DataFrame(fig1_data).to_csv(source_dir / "Fig1_source_data.csv", index=False)
    
    # Figure 2 source data (load from stratified results)
    metrics_path = RESULTS_DIR / "stratified_performance_by_mechanism.tsv"
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path, sep='\t')
        metrics.to_csv(source_dir / "Fig2_source_data.csv", index=False)
    
    print(f"Created source data in: {source_dir}")


def main():
    """Generate all supplementary tables."""
    
    print("=" * 70)
    print("GENERATING SUPPLEMENTARY TABLES FOR NATURE GENETICS")
    print("=" * 70)
    
    print("\n>>> Extended Data Table 1: Benchmark loci...")
    create_extended_data_table_1()
    
    print("\n>>> Extended Data Table 2: Leakage audit...")
    create_extended_data_table_2()
    
    print("\n>>> Extended Data Table 3: Stratified performance...")
    create_extended_data_table_3()
    
    print("\n>>> Extended Data Table 4: Gold standards analysis...")
    create_extended_data_table_4()
    
    print("\n>>> Source data for figures...")
    create_source_data()
    
    print("\n" + "=" * 70)
    print("SUPPLEMENTARY TABLES COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
