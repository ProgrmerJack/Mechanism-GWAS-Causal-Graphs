#!/usr/bin/env python3
"""
STING-seq Full Framework Validation with VARIANT-SPECIFIC ABC Integration

CRITICAL METHODOLOGICAL FIX:
The previous approaches used GENE-LEVEL ABC (max ABC for each gene globally).
This is WRONG because it promotes genes with any enhancer contact anywhere.

The CORRECT approach: For each STING-seq GWAS variant, find ABC enhancers that
OVERLAP the variant position, then use those specific enhancer-gene scores.

This is TRUE locus-specific multi-source evidence integration.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score
from typing import Dict, Tuple, List, Set

# Constants
LEAK_PROB = 0.01  # Background probability ε
WINDOW_SIZE = 10000  # bp window around variant (10kb for maximum coverage)

# Blood-relevant cell types (K562 is exact match for STING-seq)
BLOOD_CELL_TYPES = [
    'K562-Roadmap',
    'erythroblast-Corces2016',
    'CD14-positive_monocyte-Novakovic2016',
    'CD14-positive_monocyte-ENCODE',
    'CD4-positive_helper_T_cell-Corces2016',
    'CD8-positive_alpha-beta_T_cell-Corces2016',
    'B_cell-ENCODE',
    'GM12878-Roadmap',
    'natural_killer_cell-Corces2016',
]

BASE_DIR = Path(__file__).parent.parent if '__file__' in dir() else Path('.')
if not (BASE_DIR / 'data').exists():
    BASE_DIR = Path('c:/Users/Jack0/GitHub/Mechanism-GWAS-Causal-Graphs')

DATA_DIR = BASE_DIR / 'data'
OUTPUT_DIR = DATA_DIR / 'processed' / 'prospective_validation'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_gene_symbol_to_ensembl_mapping() -> Dict[str, str]:
    """Load gene symbol to Ensembl ID mapping."""
    ensg_path = DATA_DIR / 'external' / 'flames' / 'Annotation_data' / 'ENSG' / 'ENSG.v102.genes.parquet'
    
    print(f"[Mapping] Loading gene mapping...")
    df = pd.read_parquet(ensg_path, columns=['ensembl_gene_id', 'external_gene_name', 'hgnc_symbol'])
    
    mapping = {}
    for _, row in df.iterrows():
        ens_id = row['ensembl_gene_id']
        if pd.notna(row['external_gene_name']):
            mapping[row['external_gene_name']] = ens_id
        if pd.notna(row['hgnc_symbol']):
            mapping[row['hgnc_symbol']] = ens_id
    
    print(f"[Mapping] {len(mapping):,} symbol -> Ensembl mappings")
    return mapping


def load_sting_seq_variants() -> pd.DataFrame:
    """Load STING-seq variants with positions."""
    path = DATA_DIR / 'external' / 'sting_seq' / 'sting_seq_cre_gene_pairs.tsv'
    
    df = pd.read_csv(path, sep='\t', comment='#')
    df = df.rename(columns={
        'target_gene_ensembl': 'gene_id',
        'target_gene': 'gene_symbol'
    })
    df = df[df['gene_id'].notna() & (df['gene_id'] != 'NA')]
    
    # Parse chromosome and position
    df['chr'] = 'chr' + df['chr'].astype(str)
    df['position'] = df['position'].astype(int)
    
    print(f"[STING-seq] {len(df)} CRE-gene pairs at {df['rsid'].nunique()} unique variants")
    print(f"[STING-seq] Validated genes: {df['gene_id'].nunique()}")
    
    return df


def load_abc_predictions_blood() -> pd.DataFrame:
    """Load ABC predictions for blood cell types."""
    path = DATA_DIR / 'external' / 'abc' / 'AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'
    
    print(f"[ABC] Loading ABC predictions (blood cells)...")
    
    df = pd.read_csv(path, sep='\t', 
                     usecols=['chr', 'start', 'end', 'TargetGene', 'ABC.Score', 'CellType'])
    
    # Filter to blood cell types
    df = df[df['CellType'].isin(BLOOD_CELL_TYPES)]
    
    print(f"[ABC] {len(df):,} blood cell predictions")
    print(f"[ABC] Cell types: {df['CellType'].nunique()}")
    
    return df


def find_variant_overlapping_abc(variant_chr: str, variant_pos: int, 
                                  abc_df: pd.DataFrame, 
                                  gene_mapping: Dict[str, str],
                                  window: int = 500) -> Dict[str, float]:
    """
    Find ABC enhancer-gene predictions that overlap a specific variant.
    
    Returns dict of gene_id -> max ABC score for enhancers overlapping variant.
    """
    # Filter to chromosome
    chr_abc = abc_df[abc_df['chr'] == variant_chr]
    
    # Find enhancers overlapping variant ± window
    overlapping = chr_abc[
        (chr_abc['start'] <= variant_pos + window) & 
        (chr_abc['end'] >= variant_pos - window)
    ]
    
    if len(overlapping) == 0:
        return {}
    
    # Map gene symbols to Ensembl IDs
    overlapping = overlapping.copy()
    overlapping['gene_id'] = overlapping['TargetGene'].map(gene_mapping)
    overlapping = overlapping.dropna(subset=['gene_id'])
    
    if len(overlapping) == 0:
        return {}
    
    # Take max ABC per gene
    result = overlapping.groupby('gene_id')['ABC.Score'].max().to_dict()
    
    return result


def load_l2g_predictions() -> pd.DataFrame:
    """Load L2G predictions."""
    path = DATA_DIR / 'external' / 'opentargets_l2g' / 'processed' / 'l2g_processed.parquet'
    
    df = pd.read_parquet(path)
    df = df.rename(columns={'ensembl_gene_id': 'gene_id'})
    
    print(f"[L2G] {len(df):,} predictions")
    return df


def compute_path_probability(l2g_score: float, abc_score: float) -> float:
    """Noisy-OR: P = 1 - (1-ε)(1-L2G)(1-ABC)"""
    l2g = float(l2g_score) if pd.notna(l2g_score) else 0.0
    abc = float(abc_score) if pd.notna(abc_score) else 0.0
    
    return 1 - (1 - LEAK_PROB) * (1 - l2g) * (1 - abc)


def compute_auroc_ci(y_true: np.ndarray, y_score: np.ndarray, 
                     n_boot: int = 2000, seed: int = 42) -> Tuple[float, float, float]:
    """AUROC with bootstrap CI."""
    np.random.seed(seed)
    auroc = roc_auc_score(y_true, y_score)
    
    n = len(y_true)
    aurocs = []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        try:
            aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
        except:
            continue
    
    return auroc, np.percentile(aurocs, 2.5), np.percentile(aurocs, 97.5)


def compute_delta_auroc(y_true: np.ndarray, y1: np.ndarray, y2: np.ndarray,
                        n_boot: int = 2000, seed: int = 42) -> Tuple[float, float, float, float]:
    """ΔAUROC with CI and p-value."""
    np.random.seed(seed)
    delta = roc_auc_score(y_true, y1) - roc_auc_score(y_true, y2)
    
    n = len(y_true)
    deltas = []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        try:
            d = roc_auc_score(y_true[idx], y1[idx]) - roc_auc_score(y_true[idx], y2[idx])
            deltas.append(d)
        except:
            continue
    
    return delta, np.percentile(deltas, 2.5), np.percentile(deltas, 97.5), np.mean(np.array(deltas) <= 0)


def run_variant_specific_validation():
    """Run validation with VARIANT-SPECIFIC ABC integration."""
    
    print("=" * 80)
    print("STING-seq FULL FRAMEWORK VALIDATION")
    print("VARIANT-SPECIFIC ABC INTEGRATION (Enhancers overlapping GWAS variants)")
    print("=" * 80)
    print()
    
    # Load data
    gene_mapping = load_gene_symbol_to_ensembl_mapping()
    sting_seq_df = load_sting_seq_variants()
    abc_df = load_abc_predictions_blood()
    l2g_df = load_l2g_predictions()
    
    # Get unique STING-seq variants
    variant_df = sting_seq_df[['rsid', 'chr', 'position']].drop_duplicates()
    print(f"\n[Processing] {len(variant_df)} unique STING-seq variants")
    
    # For each variant, find overlapping ABC predictions
    print("[ABC] Computing variant-specific ABC scores...")
    variant_abc_scores = {}  # variant -> {gene_id: abc_score}
    
    for _, row in variant_df.iterrows():
        rsid = row['rsid']
        scores = find_variant_overlapping_abc(
            row['chr'], row['position'], abc_df, gene_mapping, WINDOW_SIZE
        )
        if scores:
            variant_abc_scores[rsid] = scores
    
    variants_with_abc = len(variant_abc_scores)
    total_genes_with_abc = sum(len(v) for v in variant_abc_scores.values())
    print(f"[ABC] {variants_with_abc}/{len(variant_df)} variants have overlapping ABC enhancers")
    print(f"[ABC] {total_genes_with_abc} variant-gene ABC scores computed")
    
    # Get validated gene set per variant
    validated = sting_seq_df.groupby('rsid')['gene_id'].apply(set).to_dict()
    all_validated_genes = set(sting_seq_df['gene_id'].unique())
    
    # Build prediction dataset per variant
    # For fair comparison, we evaluate at the GENE level across all L2G predictions
    # But add variant-specific ABC only where we have overlapping enhancers
    
    print("\n[Building] Creating locus-aware prediction dataset...")
    
    # Approach: For each variant, get L2G predictions for nearby genes
    # Add variant-specific ABC where available
    records = []
    
    for rsid, var_info in variant_df.groupby('rsid').first().iterrows():
        chr_num = var_info['chr']
        pos = var_info['position']
        
        # Get validated genes for this variant
        valid_genes = validated.get(rsid, set())
        
        # Get L2G predictions for this locus (genes within 1Mb)
        # Note: L2G doesn't have locus info, so we use all predictions
        # In practice, L2G is computed per-locus by Open Targets
        
        # Get variant-specific ABC scores
        abc_scores_here = variant_abc_scores.get(rsid, {})
        
        # For genes that have ABC scores at this variant, create records
        for gene_id, abc_score in abc_scores_here.items():
            # Look up L2G score for this gene
            l2g_row = l2g_df[l2g_df['gene_id'] == gene_id]
            if len(l2g_row) > 0:
                l2g_score = l2g_row['y_proba'].iloc[0] if 'y_proba' in l2g_row.columns else l2g_row['l2g_score'].iloc[0]
            else:
                l2g_score = 0.0
            
            is_validated = gene_id in valid_genes
            
            path_prob = compute_path_probability(l2g_score, abc_score)
            
            records.append({
                'rsid': rsid,
                'gene_id': gene_id,
                'l2g_score': l2g_score,
                'abc_score': abc_score,
                'path_probability': path_prob,
                'is_validated': int(is_validated)
            })
    
    if not records:
        print("ERROR: No variant-gene pairs with both L2G and ABC scores!")
        return None
    
    df = pd.DataFrame(records)
    print(f"\n[Dataset] {len(df)} variant-gene pairs with ABC overlap")
    print(f"[Dataset] {df['is_validated'].sum()} validated pairs")
    print(f"[Dataset] {(df['l2g_score'] > 0).sum()} pairs with L2G > 0")
    
    # Compute metrics
    print("\n" + "=" * 80)
    print("RESULTS: VARIANT-SPECIFIC NOISY-OR AGGREGATION")
    print("=" * 80)
    
    y_true = df['is_validated'].values
    y_l2g = df['l2g_score'].values
    y_path = df['path_probability'].values
    
    if y_true.sum() < 5:
        print("WARNING: Too few validated pairs for robust AUROC")
        return None
    
    print("\n--- AUROC PERFORMANCE ---\n")
    
    auroc_l2g, ci_l2g_lo, ci_l2g_hi = compute_auroc_ci(y_true, y_l2g)
    auroc_path, ci_p_lo, ci_p_hi = compute_auroc_ci(y_true, y_path)
    
    print(f"  L2G alone:                    AUROC {auroc_l2g:.4f} [{ci_l2g_lo:.3f}, {ci_l2g_hi:.3f}]")
    print(f"  Path-prob (L2G + var-ABC):    AUROC {auroc_path:.4f} [{ci_p_lo:.3f}, {ci_p_hi:.3f}]")
    
    delta, ci_lo, ci_hi, pval = compute_delta_auroc(y_true, y_path, y_l2g)
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
    print(f"\n  ΔAUROC (Path-prob vs L2G): {delta:+.4f} [{ci_lo:+.4f}, {ci_hi:+.4f}]  p={pval:.4f} {sig}")
    
    # Correlation
    corr = df[['l2g_score', 'abc_score']].corr().iloc[0, 1]
    print(f"\n  L2G-ABC correlation: r = {corr:.4f}")
    
    # Additional diagnostics
    print("\n--- DIAGNOSTICS ---\n")
    
    print(f"  Mean L2G score (validated): {df[df['is_validated']==1]['l2g_score'].mean():.4f}")
    print(f"  Mean L2G score (non-valid): {df[df['is_validated']==0]['l2g_score'].mean():.4f}")
    print(f"  Mean ABC score (validated): {df[df['is_validated']==1]['abc_score'].mean():.4f}")
    print(f"  Mean ABC score (non-valid): {df[df['is_validated']==0]['abc_score'].mean():.4f}")
    
    # Save results
    results = {
        'method': 'variant_specific_abc',
        'window_size_bp': WINDOW_SIZE,
        'n_variants_with_abc': variants_with_abc,
        'n_variant_gene_pairs': len(df),
        'n_validated_pairs': int(df['is_validated'].sum()),
        'auroc_l2g': auroc_l2g,
        'auroc_l2g_ci': [ci_l2g_lo, ci_l2g_hi],
        'auroc_path_prob': auroc_path,
        'auroc_path_prob_ci': [ci_p_lo, ci_p_hi],
        'delta_auroc': delta,
        'delta_auroc_ci': [ci_lo, ci_hi],
        'p_value': pval,
        'l2g_abc_correlation': corr,
        'blood_cell_types_used': BLOOD_CELL_TYPES
    }
    
    outpath = OUTPUT_DIR / 'variant_specific_abc_validation_results.json'
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n\nResults saved to: {outpath}")
    
    # Save detailed data for inspection
    df.to_csv(OUTPUT_DIR / 'variant_specific_abc_validation_data.csv', index=False)
    
    print("\n" + "=" * 80)
    print("VARIANT-SPECIFIC ABC VALIDATION COMPLETE")
    print("=" * 80)
    
    return results


if __name__ == '__main__':
    run_variant_specific_validation()
