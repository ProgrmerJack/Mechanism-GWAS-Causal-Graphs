#!/usr/bin/env python3
"""
STING-seq Full Framework Validation with REAL ABC Integration

This script computes path-probability using the COMPLETE noisy-OR framework:
    P = 1 - (1-ε)(1-P_L2G)(1-P_ABC)(1-P_eQTL)(1-P_coloc)

For STING-seq, we can integrate:
- L2G: Locus-to-Gene scores (machine learning based)
- ABC: Activity-by-Contact enhancer-gene predictions (experimental + Hi-C based)

The key innovation: ABC provides INDEPENDENT evidence that should BOOST performance.
"""

import pandas as pd
import numpy as np
import gzip
import json
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score
from typing import Dict, Tuple, Optional

# Constants
LEAK_PROB = 0.01  # Background probability ε

BASE_DIR = Path(__file__).parent if '__file__' in dir() else Path('.')
if not (BASE_DIR / 'data').exists():
    BASE_DIR = Path('c:/Users/Jack0/GitHub/Mechanism-GWAS-Causal-Graphs')

DATA_DIR = BASE_DIR / 'data'
OUTPUT_DIR = DATA_DIR / 'processed' / 'prospective_validation'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_gene_symbol_to_ensembl_mapping() -> Dict[str, str]:
    """Load gene symbol to Ensembl ID mapping from ENSG annotation file."""
    ensg_path = DATA_DIR / 'external' / 'flames' / 'Annotation_data' / 'ENSG' / 'ENSG.v102.genes.parquet'
    
    print(f"[Mapping] Loading gene mapping from ENSG annotations...")
    df = pd.read_parquet(ensg_path, columns=['ensembl_gene_id', 'external_gene_name', 'hgnc_symbol'])
    
    # Create multiple mappings (external_gene_name and hgnc_symbol can differ)
    mapping = {}
    for _, row in df.iterrows():
        ens_id = row['ensembl_gene_id']
        if pd.notna(row['external_gene_name']):
            mapping[row['external_gene_name']] = ens_id
        if pd.notna(row['hgnc_symbol']):
            mapping[row['hgnc_symbol']] = ens_id
    
    print(f"[Mapping] Created {len(mapping):,} gene symbol -> Ensembl ID mappings")
    return mapping


def load_abc_predictions() -> pd.DataFrame:
    """Load ABC enhancer-gene predictions."""
    abc_path = DATA_DIR / 'external' / 'abc' / 'AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'
    
    print(f"[ABC] Loading ABC predictions from {abc_path}...")
    
    # Load only needed columns
    df = pd.read_csv(abc_path, sep='\t', 
                     usecols=['TargetGene', 'ABC.Score', 'CellType', 'chr', 'start', 'end'])
    
    print(f"[ABC] Loaded {len(df):,} predictions across {df['CellType'].nunique()} cell types")
    print(f"[ABC] {df['TargetGene'].nunique():,} unique genes")
    
    return df


def compute_max_abc_per_gene(abc_df: pd.DataFrame, gene_mapping: Dict[str, str]) -> Dict[str, float]:
    """Compute maximum ABC score per gene (Ensembl ID) across all cell types."""
    print("[ABC] Computing max ABC score per gene...")
    
    # Map gene symbols to Ensembl IDs
    abc_df = abc_df.copy()
    abc_df['gene_id'] = abc_df['TargetGene'].map(gene_mapping)
    
    # Filter to genes with valid Ensembl IDs
    abc_mapped = abc_df.dropna(subset=['gene_id'])
    print(f"[ABC] {len(abc_mapped):,} predictions with valid Ensembl IDs ({abc_df['gene_id'].notna().mean()*100:.1f}%)")
    
    # Compute max ABC per gene (taking maximum across all cell types and enhancers)
    max_abc = abc_mapped.groupby('gene_id')['ABC.Score'].max().to_dict()
    print(f"[ABC] Computed max ABC for {len(max_abc):,} genes")
    print(f"[ABC] ABC score range: {min(max_abc.values()):.4f} - {max(max_abc.values()):.4f}")
    
    return max_abc


def load_sting_seq_benchmark() -> pd.DataFrame:
    """Load STING-seq benchmark data."""
    sting_path = DATA_DIR / 'external' / 'sting_seq' / 'sting_seq_cre_gene_pairs.tsv'
    
    # Read TSV, skipping comment lines
    df = pd.read_csv(sting_path, sep='\t', comment='#')
    
    # Rename columns for consistency
    df = df.rename(columns={
        'target_gene_ensembl': 'gene_id',
        'target_gene': 'gene_symbol'
    })
    
    # Filter to valid gene pairs (remove NO_TARGET rows)
    df = df[df['gene_id'].notna() & (df['gene_id'] != 'NA')]
    
    print(f"[STING-seq] Loaded {len(df)} validated CRE-gene pairs")
    print(f"  - {df['gene_id'].nunique()} unique genes")
    print(f"  - {df['rsid'].nunique()} unique variants")
    
    return df


def load_l2g_predictions() -> pd.DataFrame:
    """Load L2G predictions."""
    l2g_path = DATA_DIR / 'external' / 'opentargets_l2g' / 'processed' / 'l2g_processed.parquet'
    
    df = pd.read_parquet(l2g_path)
    
    # Rename columns for consistency
    df = df.rename(columns={'ensembl_gene_id': 'gene_id'})
    
    print(f"[L2G] Loaded {len(df):,} predictions")
    print(f"[L2G] Unique genes: {df['gene_id'].nunique():,}")
    return df


def compute_path_probability_with_abc(l2g_score: float, abc_score: float) -> float:
    """
    Compute path-probability using noisy-OR with L2G and ABC.
    
    P = 1 - (1-ε)(1-P_L2G)(1-P_ABC)
    
    This assumes L2G and ABC provide INDEPENDENT evidence.
    """
    # Handle missing values
    l2g_score = float(l2g_score) if pd.notna(l2g_score) else 0.0
    abc_score = float(abc_score) if pd.notna(abc_score) else 0.0
    
    # Noisy-OR aggregation
    prob_no_path = (1 - LEAK_PROB) * (1 - l2g_score) * (1 - abc_score)
    path_prob = 1 - prob_no_path
    
    return path_prob


def prepare_validation_dataset(l2g_df: pd.DataFrame, 
                                sting_seq_df: pd.DataFrame,
                                max_abc_scores: Dict[str, float]) -> pd.DataFrame:
    """Prepare unified validation dataset with path-probability."""
    
    # Get validated genes
    if 'gene_id' in sting_seq_df.columns:
        validated_genes = set(sting_seq_df['gene_id'].unique())
    elif 'gene_ensembl' in sting_seq_df.columns:
        validated_genes = set(sting_seq_df['gene_ensembl'].unique())
    else:
        # Map gene symbols
        gene_mapping = load_gene_symbol_to_ensembl_mapping()
        validated_genes = set(sting_seq_df['target_gene'].map(gene_mapping).dropna().unique())
    
    print(f"[Validation] {len(validated_genes)} validated genes from STING-seq")
    
    # Prepare L2G dataframe
    df = l2g_df.copy()
    df = df.rename(columns={'y_proba': 'l2g_score'}) if 'y_proba' in df.columns else df
    
    # Ensure gene_id column
    if 'gene_id' not in df.columns:
        if 'gene_ensembl' in df.columns:
            df['gene_id'] = df['gene_ensembl']
    
    # Add ABC scores
    df['abc_score'] = df['gene_id'].map(max_abc_scores).fillna(0.0)
    
    # Report ABC coverage
    has_abc = (df['abc_score'] > 0).sum()
    print(f"[ABC Coverage] {has_abc:,} of {len(df):,} predictions have ABC scores ({has_abc/len(df)*100:.1f}%)")
    
    # Check ABC coverage on validated genes
    validated_with_abc = df[df['gene_id'].isin(validated_genes) & (df['abc_score'] > 0)]
    print(f"[ABC on Validated] {validated_with_abc['gene_id'].nunique()} of {len(validated_genes)} validated genes have ABC")
    
    # Compute path-probability with ABC integration
    print("[PathProb] Computing path-probability with ABC integration...")
    df['path_probability_l2g_only'] = df['l2g_score'].apply(
        lambda x: 1 - (1 - LEAK_PROB) * (1 - (x if pd.notna(x) else 0))
    )
    df['path_probability_with_abc'] = df.apply(
        lambda row: compute_path_probability_with_abc(row['l2g_score'], row['abc_score']),
        axis=1
    )
    
    # Create validation labels
    df['is_validated'] = df['gene_id'].isin(validated_genes).astype(int)
    
    print(f"[Dataset] Final dataset: {len(df):,} predictions, {df['is_validated'].sum():,} validated")
    
    return df


def compute_auroc_with_ci(y_true: np.ndarray, y_score: np.ndarray, 
                           n_bootstrap: int = 1000, seed: int = 42) -> Tuple[float, float, float]:
    """Compute AUROC with bootstrap confidence interval."""
    np.random.seed(seed)
    
    auroc = roc_auc_score(y_true, y_score)
    
    # Bootstrap CI
    n = len(y_true)
    aurocs = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, n, replace=True)
        try:
            aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
        except ValueError:
            continue
    
    ci_low = np.percentile(aurocs, 2.5)
    ci_high = np.percentile(aurocs, 97.5)
    
    return auroc, ci_low, ci_high


def compute_delta_auroc_with_ci(y_true: np.ndarray, y_score1: np.ndarray, y_score2: np.ndarray,
                                 n_bootstrap: int = 1000, seed: int = 42) -> Tuple[float, float, float, float]:
    """Compute ΔAUROC with bootstrap CI and p-value."""
    np.random.seed(seed)
    
    auroc1 = roc_auc_score(y_true, y_score1)
    auroc2 = roc_auc_score(y_true, y_score2)
    delta = auroc1 - auroc2
    
    # Bootstrap
    n = len(y_true)
    deltas = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, n, replace=True)
        try:
            a1 = roc_auc_score(y_true[idx], y_score1[idx])
            a2 = roc_auc_score(y_true[idx], y_score2[idx])
            deltas.append(a1 - a2)
        except ValueError:
            continue
    
    ci_low = np.percentile(deltas, 2.5)
    ci_high = np.percentile(deltas, 97.5)
    
    # P-value: proportion of bootstrap samples where delta <= 0
    p_value = np.mean(np.array(deltas) <= 0)
    
    return delta, ci_low, ci_high, p_value


def run_full_framework_validation():
    """Run complete validation comparing Path-prob with ABC vs L2G alone."""
    
    print("=" * 70)
    print("STING-seq FULL FRAMEWORK VALIDATION WITH ABC INTEGRATION")
    print("True Multi-Source Evidence Integration")
    print("=" * 70)
    print()
    
    # Step 1: Load gene mapping
    print("[1/5] Loading gene symbol -> Ensembl ID mapping...")
    gene_mapping = load_gene_symbol_to_ensembl_mapping()
    
    # Step 2: Load ABC predictions
    print("\n[2/5] Loading ABC enhancer-gene predictions...")
    abc_df = load_abc_predictions()
    max_abc_scores = compute_max_abc_per_gene(abc_df, gene_mapping)
    
    # Step 3: Load STING-seq benchmark
    print("\n[3/5] Loading STING-seq benchmark...")
    sting_seq_df = load_sting_seq_benchmark()
    
    # Step 4: Load L2G predictions
    print("\n[4/5] Loading L2G predictions...")
    l2g_df = load_l2g_predictions()
    
    # Step 5: Prepare validation dataset
    print("\n[5/5] Preparing unified validation dataset...")
    df = prepare_validation_dataset(l2g_df, sting_seq_df, max_abc_scores)
    
    # Run validation
    print("\n" + "=" * 70)
    print("VALIDATION RESULTS")
    print("=" * 70)
    
    y_true = df['is_validated'].values
    y_l2g = df['l2g_score'].fillna(0).values
    y_path_l2g_only = df['path_probability_l2g_only'].values
    y_path_with_abc = df['path_probability_with_abc'].values
    
    # Compute AUROCs
    print("\n--- METHOD PERFORMANCE ---\n")
    
    auroc_l2g, ci_l2g_low, ci_l2g_high = compute_auroc_with_ci(y_true, y_l2g)
    auroc_path_l2g, ci_path_l2g_low, ci_path_l2g_high = compute_auroc_with_ci(y_true, y_path_l2g_only)
    auroc_path_abc, ci_path_abc_low, ci_path_abc_high = compute_auroc_with_ci(y_true, y_path_with_abc)
    
    print(f"  L2G alone:              AUROC {auroc_l2g:.4f} [{ci_l2g_low:.3f}, {ci_l2g_high:.3f}]")
    print(f"  Path-prob (L2G only):   AUROC {auroc_path_l2g:.4f} [{ci_path_l2g_low:.3f}, {ci_path_l2g_high:.3f}]")
    print(f"  Path-prob (L2G + ABC):  AUROC {auroc_path_abc:.4f} [{ci_path_abc_low:.3f}, {ci_path_abc_high:.3f}]")
    
    # Statistical comparison
    print("\n--- PAIRWISE COMPARISONS ---\n")
    
    delta1, ci1_low, ci1_high, p1 = compute_delta_auroc_with_ci(y_true, y_path_with_abc, y_l2g)
    delta2, ci2_low, ci2_high, p2 = compute_delta_auroc_with_ci(y_true, y_path_with_abc, y_path_l2g_only)
    
    sig1 = "***" if p1 < 0.001 else "**" if p1 < 0.01 else "*" if p1 < 0.05 else ""
    sig2 = "***" if p2 < 0.001 else "**" if p2 < 0.01 else "*" if p2 < 0.05 else ""
    
    print(f"  Path-prob (L2G+ABC) vs L2G:           ΔAUROC {delta1:+.4f} [{ci1_low:+.4f}, {ci1_high:+.4f}]  p={p1:.4f} {sig1}")
    print(f"  Path-prob (L2G+ABC) vs Path-prob(L2G): ΔAUROC {delta2:+.4f} [{ci2_low:+.4f}, {ci2_high:+.4f}]  p={p2:.4f} {sig2}")
    
    # ABC contribution analysis
    print("\n--- ABC CONTRIBUTION ANALYSIS ---\n")
    
    df_with_abc = df[df['abc_score'] > 0]
    df_without_abc = df[df['abc_score'] == 0]
    
    print(f"  Genes WITH ABC scores: {len(df_with_abc):,} ({len(df_with_abc)/len(df)*100:.1f}%)")
    print(f"  Genes WITHOUT ABC scores: {len(df_without_abc):,} ({len(df_without_abc)/len(df)*100:.1f}%)")
    
    # Performance on genes with ABC
    if df_with_abc['is_validated'].sum() > 10:
        auroc_abc_subset, _, _ = compute_auroc_with_ci(
            df_with_abc['is_validated'].values,
            df_with_abc['path_probability_with_abc'].values
        )
        auroc_l2g_subset, _, _ = compute_auroc_with_ci(
            df_with_abc['is_validated'].values,
            df_with_abc['l2g_score'].fillna(0).values
        )
        print(f"\n  On genes WITH ABC:")
        print(f"    L2G AUROC: {auroc_l2g_subset:.4f}")
        print(f"    Path-prob (L2G+ABC) AUROC: {auroc_abc_subset:.4f}")
        print(f"    Improvement: {auroc_abc_subset - auroc_l2g_subset:+.4f}")
    
    # Correlation analysis
    corr_l2g_abc = df[df['abc_score'] > 0][['l2g_score', 'abc_score']].corr().iloc[0, 1]
    print(f"\n  L2G-ABC correlation (on shared genes): {corr_l2g_abc:.4f}")
    
    # Save results
    results = {
        'auroc_l2g': auroc_l2g,
        'auroc_path_l2g_only': auroc_path_l2g,
        'auroc_path_with_abc': auroc_path_abc,
        'delta_auroc_vs_l2g': delta1,
        'delta_auroc_vs_l2g_ci': [ci1_low, ci1_high],
        'delta_p_value': p1,
        'abc_coverage_pct': len(df_with_abc) / len(df) * 100,
        'n_validated': int(df['is_validated'].sum()),
        'n_total': len(df)
    }
    
    results_path = OUTPUT_DIR / 'full_framework_with_abc_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {results_path}")
    
    print("\n" + "=" * 70)
    print("FULL FRAMEWORK VALIDATION WITH ABC COMPLETE")
    print("=" * 70)
    
    return results


if __name__ == '__main__':
    run_full_framework_validation()
