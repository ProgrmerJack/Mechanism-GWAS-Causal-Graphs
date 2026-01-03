#!/usr/bin/env python3
"""
STING-seq Full Framework Validation with BLOOD-MATCHED ABC Integration

The critical fix: Filter ABC predictions to BLOOD-RELEVANT cell types only.
STING-seq uses K562 cells, so we must match the tissue context.

This is the methodologically correct implementation of multi-source evidence integration.
"""

import pandas as pd
import numpy as np
import gzip
import json
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score
from typing import Dict, Tuple, List, Optional

# Constants
LEAK_PROB = 0.01  # Background probability ε

# Blood-relevant cell types from ABC dataset
# These are cell types relevant to blood traits and STING-seq validation
BLOOD_CELL_TYPES = [
    # Exact STING-seq cell line
    'K562-Roadmap',
    
    # Erythroid lineage (most relevant for blood traits)
    'erythroblast-Corces2016',
    
    # Monocytes/Macrophages (blood cell types)
    'CD14-positive_monocyte-Novakovic2016',
    'CD14-positive_monocyte-ENCODE',
    'CD14-positive_monocytes-Roadmap',
    'CD14-positive_monocyte_treated_with_RPMI_4h-Novakovic2016',
    'CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016',
    
    # T cells
    'CD4-positive_helper_T_cell-Corces2016',
    'CD4-positive_helper_T_cell-ENCODE',
    'CD8-positive_alpha-beta_T_cell-Corces2016',
    'CD8-positive_alpha-beta_T_cell-ENCODE',
    'CD3-positive_T_cell-Roadmap',
    'T-cell-ENCODE',
    
    # B cells
    'B_cell-ENCODE',
    'CD19-positive_B_cell-Roadmap',
    
    # NK cells
    'natural_killer_cell-Corces2016',
    'CD56-positive_natural_killer_cells-Roadmap',
    
    # Lymphoblastoid (blood-derived)
    'GM12878-Roadmap',
    
    # Dendritic cells
    'dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_6_hour-Garber2017',
    'dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_4_hour-Garber2017',
    
    # Other relevant myeloid/lymphoid
    'MM.1S-ENCODE',
    'OCI-LY7-ENCODE',
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
    
    print(f"[Mapping] Loading gene mapping from ENSG annotations...")
    df = pd.read_parquet(ensg_path, columns=['ensembl_gene_id', 'external_gene_name', 'hgnc_symbol'])
    
    mapping = {}
    for _, row in df.iterrows():
        ens_id = row['ensembl_gene_id']
        if pd.notna(row['external_gene_name']):
            mapping[row['external_gene_name']] = ens_id
        if pd.notna(row['hgnc_symbol']):
            mapping[row['hgnc_symbol']] = ens_id
    
    print(f"[Mapping] Created {len(mapping):,} gene symbol -> Ensembl ID mappings")
    return mapping


def load_abc_predictions_blood_only() -> pd.DataFrame:
    """Load ABC predictions filtered to blood-relevant cell types ONLY."""
    abc_path = DATA_DIR / 'external' / 'abc' / 'AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'
    
    print(f"[ABC] Loading ABC predictions (BLOOD CELL TYPES ONLY)...")
    print(f"[ABC] Filtering to {len(BLOOD_CELL_TYPES)} blood-relevant cell types")
    
    # Load full file
    df = pd.read_csv(abc_path, sep='\t', 
                     usecols=['TargetGene', 'ABC.Score', 'CellType', 'chr', 'start', 'end'])
    
    total_before = len(df)
    n_celltypes_before = df['CellType'].nunique()
    
    # Filter to blood cell types
    df = df[df['CellType'].isin(BLOOD_CELL_TYPES)]
    
    total_after = len(df)
    n_celltypes_after = df['CellType'].nunique()
    
    print(f"[ABC] Filtered: {total_before:,} -> {total_after:,} predictions")
    print(f"[ABC] Cell types: {n_celltypes_before} -> {n_celltypes_after}")
    print(f"[ABC] Blood cell types found: {sorted(df['CellType'].unique())}")
    print(f"[ABC] {df['TargetGene'].nunique():,} unique genes in blood cells")
    
    return df


def compute_max_abc_per_gene(abc_df: pd.DataFrame, gene_mapping: Dict[str, str]) -> Dict[str, float]:
    """Compute maximum ABC score per gene across BLOOD cell types only."""
    print("[ABC] Computing max ABC score per gene (blood cells only)...")
    
    abc_df = abc_df.copy()
    abc_df['gene_id'] = abc_df['TargetGene'].map(gene_mapping)
    
    abc_mapped = abc_df.dropna(subset=['gene_id'])
    print(f"[ABC] {len(abc_mapped):,} predictions with valid Ensembl IDs")
    
    max_abc = abc_mapped.groupby('gene_id')['ABC.Score'].max().to_dict()
    print(f"[ABC] Computed max blood-ABC for {len(max_abc):,} genes")
    
    if max_abc:
        print(f"[ABC] Score range: {min(max_abc.values()):.4f} - {max(max_abc.values()):.4f}")
    
    return max_abc


def load_sting_seq_benchmark() -> pd.DataFrame:
    """Load STING-seq benchmark data."""
    sting_path = DATA_DIR / 'external' / 'sting_seq' / 'sting_seq_cre_gene_pairs.tsv'
    
    df = pd.read_csv(sting_path, sep='\t', comment='#')
    df = df.rename(columns={
        'target_gene_ensembl': 'gene_id',
        'target_gene': 'gene_symbol'
    })
    df = df[df['gene_id'].notna() & (df['gene_id'] != 'NA')]
    
    print(f"[STING-seq] Loaded {len(df)} validated CRE-gene pairs")
    print(f"  - {df['gene_id'].nunique()} unique genes")
    
    return df


def load_l2g_predictions() -> pd.DataFrame:
    """Load L2G predictions."""
    l2g_path = DATA_DIR / 'external' / 'opentargets_l2g' / 'processed' / 'l2g_processed.parquet'
    
    df = pd.read_parquet(l2g_path)
    df = df.rename(columns={'ensembl_gene_id': 'gene_id'})
    
    print(f"[L2G] Loaded {len(df):,} predictions")
    return df


def compute_path_probability(l2g_score: float, abc_score: float) -> float:
    """Noisy-OR path probability: P = 1 - (1-ε)(1-P_L2G)(1-P_ABC)"""
    l2g = float(l2g_score) if pd.notna(l2g_score) else 0.0
    abc = float(abc_score) if pd.notna(abc_score) else 0.0
    
    prob_no_path = (1 - LEAK_PROB) * (1 - l2g) * (1 - abc)
    return 1 - prob_no_path


def compute_auroc_with_ci(y_true: np.ndarray, y_score: np.ndarray, 
                           n_bootstrap: int = 2000, seed: int = 42) -> Tuple[float, float, float]:
    """AUROC with bootstrap 95% CI."""
    np.random.seed(seed)
    
    auroc = roc_auc_score(y_true, y_score)
    
    n = len(y_true)
    aurocs = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, n, replace=True)
        try:
            aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
        except ValueError:
            continue
    
    return auroc, np.percentile(aurocs, 2.5), np.percentile(aurocs, 97.5)


def compute_delta_auroc(y_true: np.ndarray, y_score1: np.ndarray, y_score2: np.ndarray,
                        n_bootstrap: int = 2000, seed: int = 42) -> Tuple[float, float, float, float]:
    """ΔAUROC with bootstrap CI and p-value."""
    np.random.seed(seed)
    
    delta = roc_auc_score(y_true, y_score1) - roc_auc_score(y_true, y_score2)
    
    n = len(y_true)
    deltas = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, n, replace=True)
        try:
            d = roc_auc_score(y_true[idx], y_score1[idx]) - roc_auc_score(y_true[idx], y_score2[idx])
            deltas.append(d)
        except ValueError:
            continue
    
    ci_low = np.percentile(deltas, 2.5)
    ci_high = np.percentile(deltas, 97.5)
    p_value = np.mean(np.array(deltas) <= 0)
    
    return delta, ci_low, ci_high, p_value


def run_blood_matched_validation():
    """Run validation with blood-matched ABC integration."""
    
    print("=" * 80)
    print("STING-seq FULL FRAMEWORK VALIDATION")
    print("TISSUE-MATCHED ABC INTEGRATION (Blood Cell Types Only)")
    print("=" * 80)
    print()
    
    # Load data
    print("[1/5] Loading gene mapping...")
    gene_mapping = load_gene_symbol_to_ensembl_mapping()
    
    print("\n[2/5] Loading ABC predictions (BLOOD CELLS ONLY)...")
    abc_df = load_abc_predictions_blood_only()
    max_abc = compute_max_abc_per_gene(abc_df, gene_mapping)
    
    print("\n[3/5] Loading STING-seq benchmark...")
    sting_seq_df = load_sting_seq_benchmark()
    validated_genes = set(sting_seq_df['gene_id'].unique())
    
    print("\n[4/5] Loading L2G predictions...")
    l2g_df = load_l2g_predictions()
    
    print("\n[5/5] Preparing validation dataset...")
    df = l2g_df.copy()
    if 'y_proba' in df.columns:
        df = df.rename(columns={'y_proba': 'l2g_score'})
    
    # Add blood-matched ABC
    df['abc_score'] = df['gene_id'].map(max_abc).fillna(0.0)
    
    # Compute path probabilities
    df['path_prob_l2g_only'] = df['l2g_score'].apply(
        lambda x: 1 - (1 - LEAK_PROB) * (1 - (x if pd.notna(x) else 0))
    )
    df['path_prob_with_abc'] = df.apply(
        lambda r: compute_path_probability(r['l2g_score'], r['abc_score']), axis=1
    )
    
    # Labels
    df['is_validated'] = df['gene_id'].isin(validated_genes).astype(int)
    
    # Report coverage
    has_abc = (df['abc_score'] > 0).sum()
    validated_with_abc = df[df['is_validated'] == 1]['abc_score'] > 0
    print(f"\n[Coverage] {has_abc:,} of {len(df):,} predictions have blood-ABC ({has_abc/len(df)*100:.1f}%)")
    print(f"[Coverage] {validated_with_abc.sum()} of {df['is_validated'].sum()} validated genes have blood-ABC")
    
    # Validation
    print("\n" + "=" * 80)
    print("RESULTS: TISSUE-MATCHED NOISY-OR AGGREGATION")
    print("=" * 80)
    
    y_true = df['is_validated'].values
    y_l2g = df['l2g_score'].fillna(0).values
    y_path_l2g = df['path_prob_l2g_only'].values
    y_path_abc = df['path_prob_with_abc'].values
    
    # Nearest gene baseline
    if 'distance_to_tss' in df.columns:
        y_nearest = 1 / (1 + df['distance_to_tss'].fillna(1e6).values / 1e6)
    else:
        y_nearest = None
    
    print("\n--- AUROC PERFORMANCE ---\n")
    
    auroc_l2g, ci_l2g_lo, ci_l2g_hi = compute_auroc_with_ci(y_true, y_l2g)
    auroc_path_l2g, ci_pl_lo, ci_pl_hi = compute_auroc_with_ci(y_true, y_path_l2g)
    auroc_path_abc, ci_pa_lo, ci_pa_hi = compute_auroc_with_ci(y_true, y_path_abc)
    
    print(f"  L2G alone:                    AUROC {auroc_l2g:.4f} [{ci_l2g_lo:.3f}, {ci_l2g_hi:.3f}]")
    print(f"  Path-prob (L2G only):         AUROC {auroc_path_l2g:.4f} [{ci_pl_lo:.3f}, {ci_pl_hi:.3f}]")
    print(f"  Path-prob (L2G + Blood-ABC):  AUROC {auroc_path_abc:.4f} [{ci_pa_lo:.3f}, {ci_pa_hi:.3f}]")
    
    if y_nearest is not None:
        auroc_near, ci_n_lo, ci_n_hi = compute_auroc_with_ci(y_true, y_nearest)
        print(f"  Nearest gene baseline:        AUROC {auroc_near:.4f} [{ci_n_lo:.3f}, {ci_n_hi:.3f}]")
    
    # Statistical comparisons
    print("\n--- PAIRWISE COMPARISONS ---\n")
    
    delta, ci_lo, ci_hi, pval = compute_delta_auroc(y_true, y_path_abc, y_l2g)
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
    print(f"  Path-prob(L2G+ABC) vs L2G:  ΔAUROC = {delta:+.4f} [{ci_lo:+.4f}, {ci_hi:+.4f}]  p={pval:.4f} {sig}")
    
    delta2, ci2_lo, ci2_hi, pval2 = compute_delta_auroc(y_true, y_path_abc, y_path_l2g)
    sig2 = "***" if pval2 < 0.001 else "**" if pval2 < 0.01 else "*" if pval2 < 0.05 else ""
    print(f"  Path-prob(L2G+ABC) vs L2G-only:  ΔAUROC = {delta2:+.4f} [{ci2_lo:+.4f}, {ci2_hi:+.4f}]  p={pval2:.4f} {sig2}")
    
    # Subset analysis
    print("\n--- ABC CONTRIBUTION ANALYSIS ---\n")
    
    df_with_abc = df[df['abc_score'] > 0]
    print(f"  Predictions with blood-ABC: {len(df_with_abc):,} ({len(df_with_abc)/len(df)*100:.1f}%)")
    
    if df_with_abc['is_validated'].sum() >= 10:
        auroc_abc_sub, _, _ = compute_auroc_with_ci(
            df_with_abc['is_validated'].values,
            df_with_abc['path_prob_with_abc'].values
        )
        auroc_l2g_sub, _, _ = compute_auroc_with_ci(
            df_with_abc['is_validated'].values,
            df_with_abc['l2g_score'].fillna(0).values
        )
        print(f"\n  On genes WITH blood-ABC:")
        print(f"    L2G AUROC: {auroc_l2g_sub:.4f}")
        print(f"    Path-prob (L2G+ABC): {auroc_abc_sub:.4f}")
        print(f"    Δ = {auroc_abc_sub - auroc_l2g_sub:+.4f}")
    
    # Correlation analysis
    corr_df = df[df['abc_score'] > 0][['l2g_score', 'abc_score']].dropna()
    if len(corr_df) > 10:
        corr = corr_df.corr().iloc[0, 1]
        print(f"\n  L2G-ABC correlation: r = {corr:.4f}")
    
    # Save results
    results = {
        'method': 'blood_tissue_matched_abc',
        'n_blood_cell_types': len(BLOOD_CELL_TYPES),
        'blood_cell_types_found': sorted(abc_df['CellType'].unique().tolist()),
        'auroc_l2g': auroc_l2g,
        'auroc_l2g_ci': [ci_l2g_lo, ci_l2g_hi],
        'auroc_path_l2g_only': auroc_path_l2g,
        'auroc_path_l2g_only_ci': [ci_pl_lo, ci_pl_hi],
        'auroc_path_with_blood_abc': auroc_path_abc,
        'auroc_path_with_blood_abc_ci': [ci_pa_lo, ci_pa_hi],
        'delta_auroc_vs_l2g': delta,
        'delta_auroc_vs_l2g_ci': [ci_lo, ci_hi],
        'delta_p_value': pval,
        'abc_coverage_pct': len(df_with_abc) / len(df) * 100,
        'n_validated': int(df['is_validated'].sum()),
        'n_total': len(df),
        'improved_over_all_celltypes': auroc_path_abc > 0.6232  # Previous result
    }
    
    outpath = OUTPUT_DIR / 'blood_matched_abc_validation_results.json'
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n\nResults saved to: {outpath}")
    
    print("\n" + "=" * 80)
    print("BLOOD-MATCHED ABC VALIDATION COMPLETE")
    print("=" * 80)
    
    return results


if __name__ == '__main__':
    run_blood_matched_validation()
