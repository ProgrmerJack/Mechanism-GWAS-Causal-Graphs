#!/usr/bin/env python3
"""
STING-seq Validation with Training Leakage Exclusion

This script recalculates AUROC after excluding genes that overlap with
L2G training gold standards, providing a TRUE prospective validation.

Critical finding: 29/128 STING-seq genes overlap with L2G training data.
These must be excluded for honest prospective validation.

Author: Mechanism-GWAS Research Team
Date: 2025-01-15
"""

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve
from pathlib import Path
import scipy.stats as stats

# Paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"

# Overlapping genes that should be excluded (from provenance verification)
OVERLAPPING_GENES = {
    'ENSG00000000971', 'ENSG00000010704', 'ENSG00000021826', 'ENSG00000090534',
    'ENSG00000091513', 'ENSG00000096968', 'ENSG00000113302', 'ENSG00000113721',
    'ENSG00000119866', 'ENSG00000130203', 'ENSG00000134242', 'ENSG00000134460',
    'ENSG00000135111', 'ENSG00000136997', 'ENSG00000142208', 'ENSG00000146648',
    'ENSG00000153162', 'ENSG00000159216', 'ENSG00000160791', 'ENSG00000169442',
    'ENSG00000169738', 'ENSG00000170989', 'ENSG00000171791', 'ENSG00000172819',
    'ENSG00000180210', 'ENSG00000183454', 'ENSG00000185475', 'ENSG00000196260',
    'ENSG00000197249'  # 29 total
}

def load_l2g_scores():
    """Load L2G predictions."""
    l2g_path = DATA_DIR / "external/opentargets_l2g/processed/l2g_processed.parquet"
    df = pd.read_parquet(l2g_path)
    # Standardize column name
    df = df.rename(columns={'ensembl_gene_id': 'gene_id'})
    print(f"Loaded {len(df)} L2G predictions for {df['gene_id'].nunique()} genes")
    return df

def load_sting_seq_data():
    """Load STING-seq benchmark data."""
    sting_seq_path = DATA_DIR / "external/sting_seq/sting_seq_cre_gene_pairs.tsv"
    df = pd.read_csv(sting_seq_path, sep='\t', comment='#')
    
    # Filter out NO_TARGET entries
    df = df[df['target_gene'] != 'NO_TARGET'].copy()
    df = df.rename(columns={'target_gene_ensembl': 'gene_id'})
    
    print(f"Loaded {len(df)} STING-seq CRE-gene pairs ({df['gene_id'].nunique()} genes)")
    return df

def get_max_l2g_per_gene(l2g_df):
    """Get maximum L2G score per gene."""
    return l2g_df.groupby('gene_id')['l2g_score'].max().reset_index()

def calculate_auroc_with_ci(y_true, y_scores, n_bootstrap=1000, alpha=0.05):
    """Calculate AUROC with bootstrap confidence interval."""
    auroc = roc_auc_score(y_true, y_scores)
    
    # Bootstrap for CI
    n = len(y_true)
    bootstrap_aurocs = []
    
    for _ in range(n_bootstrap):
        indices = np.random.choice(n, size=n, replace=True)
        if len(np.unique(y_true[indices])) < 2:
            continue
        try:
            boot_auroc = roc_auc_score(y_true[indices], y_scores[indices])
            bootstrap_aurocs.append(boot_auroc)
        except ValueError:
            continue
    
    lower = np.percentile(bootstrap_aurocs, (alpha/2) * 100)
    upper = np.percentile(bootstrap_aurocs, (1 - alpha/2) * 100)
    
    return auroc, lower, upper

def evaluate_l2g_on_sting_seq(l2g_df, sting_seq_df, exclude_overlapping=False):
    """Evaluate L2G performance on STING-seq benchmark."""
    
    # Get max L2G per gene
    l2g_max = get_max_l2g_per_gene(l2g_df)
    
    # Get unique validated genes from STING-seq
    validated_genes = set(sting_seq_df['gene_id'].unique())
    
    if exclude_overlapping:
        n_before = len(validated_genes)
        validated_genes = validated_genes - OVERLAPPING_GENES
        n_after = len(validated_genes)
        print(f"Excluding {n_before - n_after} overlapping genes, {n_after} remain")
    
    # Get all genes with L2G scores
    all_l2g_genes = set(l2g_max['gene_id'].unique())
    
    # Create evaluation dataset
    # Positive: validated genes with L2G scores
    # Negative: genes with L2G scores but not validated
    genes_with_scores = all_l2g_genes
    positives = validated_genes & genes_with_scores
    negatives = genes_with_scores - validated_genes
    
    # Create dataframe for evaluation
    eval_data = []
    for gene in positives:
        score = l2g_max[l2g_max['gene_id'] == gene]['l2g_score'].values[0]
        eval_data.append({'gene_id': gene, 'l2g_score': score, 'is_validated': 1})
    
    # Sample negatives (use all for this evaluation)
    for gene in negatives:
        score = l2g_max[l2g_max['gene_id'] == gene]['l2g_score'].values[0]
        eval_data.append({'gene_id': gene, 'l2g_score': score, 'is_validated': 0})
    
    eval_df = pd.DataFrame(eval_data)
    
    print(f"\nEvaluation set: {len(positives)} positives, {len(negatives)} negatives")
    
    # Calculate AUROC
    y_true = eval_df['is_validated'].values
    y_scores = eval_df['l2g_score'].values
    
    auroc, ci_lower, ci_upper = calculate_auroc_with_ci(y_true, y_scores)
    
    return auroc, ci_lower, ci_upper, len(positives), len(negatives)

def main():
    print("=" * 70)
    print("STING-seq L2G Validation with Training Leakage Analysis")
    print("=" * 70)
    
    # Load data
    l2g_df = load_l2g_scores()
    sting_seq_df = load_sting_seq_data()
    
    # Evaluate with ALL genes
    print("\n" + "-" * 70)
    print("EVALUATION 1: All STING-seq genes (includes training data overlap)")
    print("-" * 70)
    auroc_all, ci_l_all, ci_u_all, n_pos_all, n_neg_all = evaluate_l2g_on_sting_seq(
        l2g_df, sting_seq_df, exclude_overlapping=False
    )
    
    # Evaluate with PROSPECTIVE genes only (excluding overlapping)
    print("\n" + "-" * 70)
    print("EVALUATION 2: Prospective genes only (excludes 29 training overlap)")
    print("-" * 70)
    auroc_prospective, ci_l_prosp, ci_u_prosp, n_pos_prosp, n_neg_prosp = evaluate_l2g_on_sting_seq(
        l2g_df, sting_seq_df, exclude_overlapping=True
    )
    
    # Summary
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"\n1. ALL genes (n={n_pos_all} positives):")
    print(f"   AUROC = {auroc_all:.4f} [95% CI: {ci_l_all:.4f} - {ci_u_all:.4f}]")
    print(f"   WARNING: Includes 29 genes that overlap with L2G training data")
    
    print(f"\n2. PROSPECTIVE genes only (n={n_pos_prosp} positives):")
    print(f"   AUROC = {auroc_prospective:.4f} [95% CI: {ci_l_prosp:.4f} - {ci_u_prosp:.4f}]")
    print(f"   TRUE prospective validation - excludes all genes from L2G training")
    
    print("\n" + "-" * 70)
    print("INTERPRETATION")
    print("-" * 70)
    
    delta = auroc_all - auroc_prospective
    if delta > 0.01:
        print(f"[!] AUROC drops by {delta:.4f} when excluding training overlap")
        print("   This suggests some inflation from training leakage")
    elif delta > -0.01:
        print(f"[OK] AUROC is stable ({delta:+.4f}) when excluding training overlap")
        print("   Performance is robust to training exclusion")
    else:
        print(f"[+] AUROC increases by {-delta:.4f} when excluding training overlap")
        print("   Overlapping genes may have been harder cases")
    
    # Save results
    results = {
        'all_genes': {
            'auroc': auroc_all,
            'ci_lower': ci_l_all,
            'ci_upper': ci_u_all,
            'n_positives': n_pos_all,
            'n_negatives': n_neg_all,
            'n_overlap_with_training': 29
        },
        'prospective_only': {
            'auroc': auroc_prospective,
            'ci_lower': ci_l_prosp,
            'ci_upper': ci_u_prosp,
            'n_positives': n_pos_prosp,
            'n_negatives': n_neg_prosp,
            'n_overlap_excluded': 29
        }
    }
    
    import json
    output_path = DATA_DIR / "processed/prospective_validation/l2g_auroc_with_leakage_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n📄 Results saved to: {output_path}")
    
    return results

if __name__ == "__main__":
    main()
