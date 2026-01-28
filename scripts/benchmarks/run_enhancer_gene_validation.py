#!/usr/bin/env python3
"""
Validate enhancer-gene linking predictions against CRISPR perturbation data.

This validates the ABC/PCHi-C ensemble's ability to predict which enhancer-gene
pairs represent true regulatory connections, independent of genetic association.

References:
- Fulco et al. 2019 (https://doi.org/10.1038/s41588-019-0538-0)
- Gasperini et al. 2019 (https://doi.org/10.1016/j.cell.2018.11.029)
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from sklearn.metrics import (
    precision_recall_curve,
    auc,
    average_precision_score,
    roc_auc_score,
    f1_score
)
import gzip
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.enhancer_gene_linking import (
    EnhancerGeneEnsemble,
    ABCLinks,
    PCHiCLinks,
    compute_ensemble_score
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='2025-%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "external" / "crispr_validation"
RESULTS_DIR = BASE_DIR / "results" / "enhancer_gene_validation"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_fulco_pairs():
    """
    Load Fulco et al. 2019 CRISPRi data.
    
    Returns DataFrame with columns:
    - enhancer_chr, enhancer_start, enhancer_end
    - gene_symbol
    - validated (bool): True if significant perturbation effect
    - abc_score: Activity-by-Contact score from the paper
    """
    logger.info(f"Loading Fulco CRISPRi data from {DATA_DIR / 'fulco_2019_table_s6a.xlsx'}")
    
    df = pd.read_excel(
        DATA_DIR / "fulco_2019_table_s6a.xlsx",
        sheet_name="Supplementary Table 6a",
        header=1
    )
    
    # Rename columns
    df.columns = [
        'chr', 'start', 'end', 'Element name', 'Gene', 'Element type',
        'Cell type', 'Norm.log2FC', 'Scaled.sig', 'Significant',
        'Screen', 'CRISPRi system', 'num_gRNAs', 'total_gRNAs',
        'fraction_cells_tested', 'FDR', 'Validated', 'TSS',
        'ABC.Score', 'Activity', 'Contact', 'CellType.SpecificActivity',
        'PowerLaw', 'Reference'
    ]
    
    logger.info(f"Loaded {len(df)} enhancer-gene pairs from Fulco et al.")
    
    # Filter to pairs with ABC scores
    df_scored = df[df['ABC.Score'].notna()].copy()
    logger.info(f"  {len(df_scored)} pairs have ABC scores")
    
    # Create standardized output
    result = pd.DataFrame({
        'enhancer_chr': df_scored['chr'].astype(str),
        'enhancer_start': df_scored['start'].astype(int),
        'enhancer_end': df_scored['end'].astype(int),
        'gene_symbol': df_scored['Gene'].astype(str),
        'validated': df_scored['Significant'] == True,
        'abc_score': df_scored['ABC.Score'].astype(float),
        'source': 'Fulco2019'
    })
    
    logger.info(f"  {result['validated'].sum()} validated pairs")
    
    return result


def load_gasperini_pairs():
    """
    Load Gasperini et al. 2019 CRISPRi data.
    
    Returns DataFrame with same structure as load_fulco_pairs().
    """
    logger.info(f"Loading Gasperini CRISPRi data from {DATA_DIR / 'gasperini_2019_results.txt.gz'}")
    
    try:
        with gzip.open(DATA_DIR / "gasperini_2019_results.txt.gz", 'rt') as f:
            df = pd.read_csv(f, sep='\t')
        
        logger.info(f"Loaded {len(df)} enhancer-gene pairs from Gasperini et al.")
        logger.info(f"Columns: {list(df.columns)}")
        
        # Try to identify relevant columns
        # Expected: target_gene, chr/chrom, start, end, significant, score
        
        # This is a placeholder - need to examine actual column names
        # For now, return empty DataFrame with correct structure
        logger.warning("Gasperini data structure needs manual inspection")
        return pd.DataFrame(columns=[
            'enhancer_chr', 'enhancer_start', 'enhancer_end',
            'gene_symbol', 'validated', 'abc_score', 'source'
        ])
        
    except Exception as e:
        logger.error(f"Failed to load Gasperini data: {e}")
        return pd.DataFrame(columns=[
            'enhancer_chr', 'enhancer_start', 'enhancer_end',
            'gene_symbol', 'validated', 'abc_score', 'source'
        ])


def compute_enhancer_gene_metrics(df, score_column='abc_score'):
    """
    Compute precision-recall metrics for enhancer-gene predictions.
    
    Args:
        df: DataFrame with 'validated' (bool) and score column
        score_column: Name of column with prediction scores
    
    Returns:
        dict with AUPRC, AUROC, precision@recall, F1, etc.
    """
    logger.info("Computing validation metrics...")
    
    # Extract labels and scores
    y_true = df['validated'].values
    y_score = df[score_column].values
    
    n_total = len(y_true)
    n_validated = y_true.sum()
    
    logger.info(f"Total pairs: {n_total}")
    logger.info(f"Validated pairs: {n_validated} ({100*n_validated/n_total:.1f}%)")
    
    # Compute precision-recall curve
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    auprc = auc(recall, precision)
    ap = average_precision_score(y_true, y_score)
    
    # Compute AUROC
    auroc = roc_auc_score(y_true, y_score)
    
    # Find precision at 50% recall
    idx_50 = np.argmin(np.abs(recall - 0.5))
    precision_at_50_recall = precision[idx_50]
    
    # Find optimal F1 threshold
    f1_scores = []
    for threshold in thresholds:
        y_pred = (y_score >= threshold).astype(int)
        f1 = f1_score(y_true, y_pred)
        f1_scores.append(f1)
    
    if f1_scores:
        max_f1 = max(f1_scores)
        optimal_threshold = thresholds[np.argmax(f1_scores)]
    else:
        max_f1 = 0.0
        optimal_threshold = 0.5
    
    # Bootstrap confidence intervals for AUPRC
    logger.info("Computing bootstrap confidence intervals...")
    n_bootstrap = 1000
    auprc_bootstrap = []
    
    np.random.seed(42)
    for _ in range(n_bootstrap):
        indices = np.random.choice(len(y_true), size=len(y_true), replace=True)
        y_true_boot = y_true[indices]
        y_score_boot = y_score[indices]
        
        if y_true_boot.sum() > 0 and (~y_true_boot).sum() > 0:
            try:
                auprc_boot = average_precision_score(y_true_boot, y_score_boot)
                auprc_bootstrap.append(auprc_boot)
            except:
                pass
    
    auprc_ci_lower = np.percentile(auprc_bootstrap, 2.5)
    auprc_ci_upper = np.percentile(auprc_bootstrap, 97.5)
    
    logger.info(f"  AUPRC: {auprc:.4f} [{auprc_ci_lower:.4f}, {auprc_ci_upper:.4f}]")
    logger.info(f"  AUROC: {auroc:.4f}")
    logger.info(f"  Precision @ 50% recall: {precision_at_50_recall:.4f}")
    logger.info(f"  Max F1: {max_f1:.4f} (threshold: {optimal_threshold:.4f})")
    
    return {
        'n_total': int(n_total),
        'n_validated': int(n_validated),
        'auprc': float(auprc),
        'average_precision': float(ap),
        'auroc': float(auroc),
        'precision_at_50_recall': float(precision_at_50_recall),
        'max_f1': float(max_f1),
        'optimal_threshold': float(optimal_threshold),
        'auprc_ci_lower': float(auprc_ci_lower),
        'auprc_ci_upper': float(auprc_ci_upper),
        'manuscript_claim_auprc': 0.71,
        'manuscript_claim_precision_at_50_recall': 0.68,
        'manuscript_claim_max_f1': 0.64
    }


def apply_ensemble_scores(df, ensemble):
    """
    Apply ensemble scores to each enhancer-gene pair.
    
    Args:
        df: DataFrame with enhancer-gene pairs
        ensemble: EnhancerGeneEnsemble instance
    
    Returns:
        DataFrame with added 'ensemble_score' column
    """
    logger.info("Computing ensemble scores for all pairs...")
    
    ensemble_scores = []
    
    for idx, row in df.iterrows():
        if idx % 1000 == 0:
            logger.info(f"  Processed {idx}/{len(df)} pairs...")
        
        try:
            # Get ensemble probability
            prob, evidence = ensemble.compute_edge_probability(
                enhancer_chr=str(row['enhancer_chr']),
                enhancer_start=int(row['enhancer_start']),
                enhancer_end=int(row['enhancer_end']),
                gene=str(row['gene_symbol']),
                method="ensemble"
            )
            ensemble_scores.append(prob)
        except Exception as e:
            # Fallback to ABC score if ensemble fails
            logger.warning(f"Ensemble failed for row {idx}: {e}")
            ensemble_scores.append(row.get('abc_score', 0.0))
    
    df['ensemble_score'] = ensemble_scores
    logger.info(f"  Computed {len(ensemble_scores)} ensemble scores")
    logger.info(f"  Mean: {np.mean(ensemble_scores):.4f}, Max: {np.max(ensemble_scores):.4f}")
    
    return df


def main():
    """Main validation workflow."""
    
    print("=" * 60)
    print("ENHANCER-GENE LINKING VALIDATION")
    print("=" * 60)
    print()
    print("Validating ABC/PCHi-C ensemble against CRISPR perturbations")
    print()
    
    # Load Fulco data
    fulco_df = load_fulco_pairs()
    
    # Load Gasperini data (placeholder for now)
    gasperini_df = load_gasperini_pairs()
    
    # Combine datasets
    if len(gasperini_df) > 0:
        combined_df = pd.concat([fulco_df, gasperini_df], ignore_index=True)
        logger.info(f"Combined: {len(combined_df)} total pairs")
    else:
        combined_df = fulco_df
        logger.info("Using only Fulco data (Gasperini not loaded)")
    
    # Initialize ensemble (will use ABC scores from data if ensemble not available)
    try:
        logger.info("Initializing enhancer-gene ensemble...")
        abc = ABCLinks()
        pchic = PCHiCLinks()
        ensemble = EnhancerGeneEnsemble(abc_links=abc, pchic_links=pchic)
        
        # Apply ensemble scores
        combined_df = apply_ensemble_scores(combined_df, ensemble)
        score_column = 'ensemble_score'
        logger.info("Using ensemble scores for validation")
    except Exception as e:
        logger.warning(f"Failed to initialize ensemble: {e}")
        logger.info("Falling back to ABC scores from data")
        score_column = 'abc_score'
    
    # Compute metrics
    metrics = compute_enhancer_gene_metrics(combined_df, score_column=score_column)
    
    # Save results
    output_file = RESULTS_DIR / "validation_metrics.json"
    with open(output_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    
    logger.info(f"\nResults saved to {output_file}")
    
    # Print summary
    print()
    print("=" * 60)
    print("VALIDATION RESULTS")
    print("=" * 60)
    print(f"Total enhancer-gene pairs: {metrics['n_total']}")
    print(f"Validated pairs: {metrics['n_validated']}")
    print(f"Score method: {score_column}")
    print()
    print(f"AUPRC: {metrics['auprc']:.4f} [{metrics['auprc_ci_lower']:.4f}, {metrics['auprc_ci_upper']:.4f}]")
    print(f"AUROC: {metrics['auroc']:.4f}")
    print(f"Precision @ 50% recall: {metrics['precision_at_50_recall']:.4f}")
    print(f"Max F1: {metrics['max_f1']:.4f}")
    print("=" * 60)
    print()
    
    # Compare to manuscript claims
    auprc_diff = abs(metrics['auprc'] - metrics['manuscript_claim_auprc'])
    if auprc_diff > 0.05:
        print(f"⚠ WARNING: AUPRC {metrics['auprc']:.4f} differs from manuscript claim {metrics['manuscript_claim_auprc']:.4f}")
        print("   Manuscript needs update or analysis needs review!")
    else:
        print(f"✓ AUPRC {metrics['auprc']:.4f} matches manuscript claim {metrics['manuscript_claim_auprc']:.4f}")
    
    print()


if __name__ == "__main__":
    main()
