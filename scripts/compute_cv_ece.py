#!/usr/bin/env python3
"""
Compute Proper 5-Fold Cross-Validated ECE

This script implements the 5-fold cross-validation with locus-level stratification
as described in the manuscript methodology.
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.isotonic import IsotonicRegression
from pathlib import Path


def compute_ece(y_true, y_prob, n_bins=10):
    """Compute Expected Calibration Error."""
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    
    bins = np.linspace(0, 1, n_bins + 1)
    bin_indices = np.digitize(y_prob, bins[:-1]) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    ece = 0.0
    n = len(y_true)
    
    for i in range(n_bins):
        mask = bin_indices == i
        if mask.sum() > 0:
            bin_acc = y_true[mask].mean()
            bin_conf = y_prob[mask].mean()
            bin_weight = mask.sum() / n
            ece += bin_weight * abs(bin_acc - bin_conf)
    
    return ece


def load_predictions():
    """Load gene predictions with calibration."""
    # Try validation_bundle first
    vb_path = Path('validation_bundle/calibration/gene_predictions_with_calibration.tsv')
    if vb_path.exists():
        return pd.read_csv(vb_path, sep='\t')
    
    # Try results directory
    res_path = Path('results/calibration_validation/gene_predictions_with_calibration.tsv')
    if res_path.exists():
        return pd.read_csv(res_path, sep='\t')
    
    raise FileNotFoundError("Could not find gene predictions file")


def main():
    print("=" * 70)
    print("5-FOLD CROSS-VALIDATED ECE COMPUTATION")
    print("=" * 70)
    
    # Load predictions
    print("\n1. Loading predictions...")
    df = load_predictions()
    print(f"   Loaded {len(df)} predictions")
    print(f"   Columns: {list(df.columns)}")
    
    # Check for required columns
    required_cols = ['probability', 'is_true_positive']
    if 'disease' in df.columns:
        required_cols.append('disease')
    
    for col in required_cols[:2]:
        if col not in df.columns:
            # Try alternative column names
            if 'predicted_probability' in df.columns:
                df['probability'] = df['predicted_probability']
            if 'true_positive' in df.columns:
                df['is_true_positive'] = df['true_positive']
            if 'label' in df.columns:
                df['is_true_positive'] = df['label']
    
    print(f"\n   Available columns: {list(df.columns)}")
    
    # Map actual column names
    # truth column is boolean, calibrated_prob is the probability
    if 'calibrated_prob' in df.columns:
        y_prob = df['calibrated_prob'].values.astype(float)
    elif 'probability' in df.columns:
        y_prob = df['probability'].values.astype(float)
    else:
        raise ValueError(f"Could not find probability column. Available: {list(df.columns)}")
    
    if 'truth' in df.columns:
        y_true = df['truth'].values.astype(int)
    elif 'is_true_positive' in df.columns:
        y_true = df['is_true_positive'].values.astype(int)
    else:
        raise ValueError(f"Could not find truth column. Available: {list(df.columns)}")
    
    # Locus stratification - use Disease column
    if 'Disease' in df.columns:
        locus_ids = df['Disease'].values
    elif 'disease' in df.columns:
        locus_ids = df['disease'].values
    elif 'locus' in df.columns:
        locus_ids = df['locus'].values
    else:
        locus_ids = np.arange(len(df))
    
    print(f"\n2. Data summary:")
    print(f"   Total predictions: {len(y_prob)}")
    print(f"   Positive rate: {y_true.mean():.3%}")
    print(f"   Unique loci/diseases: {len(np.unique(locus_ids))}")
    
    # Pre-isotonic ECE
    pre_ece = compute_ece(y_true, y_prob)
    print(f"\n3. Pre-isotonic ECE: {pre_ece:.4f}")
    
    # 5-fold cross-validated ECE with isotonic calibration
    print("\n4. 5-Fold CV ECE with Isotonic Calibration...")
    
    # Create unique locus mapping for stratification
    unique_loci = np.unique(locus_ids)
    locus_to_idx = {l: i for i, l in enumerate(unique_loci)}
    locus_indices = np.array([locus_to_idx[l] for l in locus_ids])
    
    # Create locus-level labels for stratification
    locus_labels = np.array([y_true[locus_ids == l].max() for l in unique_loci])
    
    # 5-fold CV
    n_folds = 5
    fold_eces = []
    all_calibrated_probs = np.zeros_like(y_prob)
    all_fold_assignments = np.zeros(len(y_prob), dtype=int)
    
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    
    for fold, (train_locus_idx, test_locus_idx) in enumerate(skf.split(unique_loci, locus_labels)):
        # Get train/test loci
        train_loci = set(unique_loci[train_locus_idx])
        test_loci = set(unique_loci[test_locus_idx])
        
        # Get train/test masks at prediction level
        train_mask = np.array([l in train_loci for l in locus_ids])
        test_mask = np.array([l in test_loci for l in locus_ids])
        
        # Fit isotonic on training fold
        iso = IsotonicRegression(out_of_bounds='clip')
        iso.fit(y_prob[train_mask], y_true[train_mask])
        
        # Transform test fold
        calibrated_probs = iso.transform(y_prob[test_mask])
        all_calibrated_probs[test_mask] = calibrated_probs
        all_fold_assignments[test_mask] = fold
        
        # Compute ECE on test fold
        fold_ece = compute_ece(y_true[test_mask], calibrated_probs)
        fold_eces.append(fold_ece)
        
        print(f"   Fold {fold+1}: ECE = {fold_ece:.4f} (n_test={test_mask.sum()})")
    
    # Average ECE across folds
    mean_ece = np.mean(fold_eces)
    std_ece = np.std(fold_eces)
    
    # Bootstrap CI
    np.random.seed(42)
    n_bootstrap = 1000
    boot_eces = []
    
    for _ in range(n_bootstrap):
        # Sample with replacement at locus level
        boot_loci = np.random.choice(unique_loci, size=len(unique_loci), replace=True)
        boot_mask = np.isin(locus_ids, boot_loci)
        boot_ece = compute_ece(y_true[boot_mask], all_calibrated_probs[boot_mask])
        boot_eces.append(boot_ece)
    
    boot_eces = np.array(boot_eces)
    ci_lower = np.percentile(boot_eces, 2.5)
    ci_upper = np.percentile(boot_eces, 97.5)
    
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"\nPre-isotonic ECE: {pre_ece:.4f}")
    print(f"\nPost-isotonic 5-Fold CV ECE:")
    print(f"  Mean:   {mean_ece:.4f}")
    print(f"  Std:    {std_ece:.4f}")
    print(f"  95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
    print(f"\nManuscript format: ECE = {mean_ece:.3f} [{ci_lower:.3f}--{ci_upper:.3f}]")
    
    # Also compute weighted average across diseases for comparison
    if 'Disease' in df.columns or 'disease' in df.columns:
        disease_col = 'Disease' if 'Disease' in df.columns else 'disease'
        print("\n" + "-" * 70)
        print("Comparison with disease-level weighted average:")
        disease_eces = []
        disease_weights = []
        for disease in df[disease_col].unique():
            dmask = df[disease_col] == disease
            d_ece = compute_ece(y_true[dmask], all_calibrated_probs[dmask])
            disease_eces.append(d_ece)
            disease_weights.append(dmask.sum())
        
        weighted_ece = np.average(disease_eces, weights=disease_weights)
        print(f"  Disease-weighted ECE: {weighted_ece:.4f}")
    
    # Save results
    results = {
        'pre_isotonic_ece': float(pre_ece),
        'cv_mean_ece': float(mean_ece),
        'cv_std_ece': float(std_ece),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
        'fold_eces': [float(e) for e in fold_eces],
        'n_predictions': int(len(y_prob)),
        'n_true_positives': int(y_true.sum()),
        'n_loci': int(len(unique_loci)),
        'n_folds': n_folds
    }
    
    import json
    output_path = Path('results/calibration_validation/cv_ece_results.json')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_path}")


if __name__ == '__main__':
    main()
