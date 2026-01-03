#!/usr/bin/env python3
"""
Bootstrap Confidence Intervals for All Benchmarks
==================================================

Calculates 95% confidence intervals using bootstrap resampling for:
1. Post-2021 Independent Benchmark (n=63)
2. Expanded Regulatory Benchmark (n=440)
3. cS2G comparison (n=21)

This addresses Nature Genetics reviewer concerns about statistical rigor
with formal uncertainty quantification.

Author: Mechanism-GWAS-Causal-Graphs
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def bootstrap_accuracy(
    predictions: np.ndarray, 
    truths: np.ndarray, 
    n_bootstrap: int = 10000,
    confidence_level: float = 0.95
) -> Dict:
    """
    Calculate bootstrap confidence interval for accuracy.
    
    For binary predictions against binary truth labels:
    - predictions: 1 = predicted positive, 0 = predicted negative
    - truths: 1 = actually positive, 0 = actually negative
    
    Returns accuracy CI.
    """
    n = len(predictions)
    if n == 0:
        return {'accuracy': np.nan, 'ci_lower': np.nan, 'ci_upper': np.nan, 'n': 0}
    
    # Calculate observed accuracy
    observed_accuracy = np.mean(predictions == truths)
    
    # Bootstrap resampling
    rng = np.random.default_rng(42)
    bootstrap_accuracies = []
    
    for _ in range(n_bootstrap):
        indices = rng.choice(n, size=n, replace=True)
        boot_pred = predictions[indices]
        boot_truth = truths[indices]
        boot_acc = np.mean(boot_pred == boot_truth)
        bootstrap_accuracies.append(boot_acc)
    
    bootstrap_accuracies = np.array(bootstrap_accuracies)
    
    # Calculate percentile CI
    alpha = 1 - confidence_level
    ci_lower = np.percentile(bootstrap_accuracies, 100 * alpha / 2)
    ci_upper = np.percentile(bootstrap_accuracies, 100 * (1 - alpha / 2))
    
    return {
        'accuracy': float(observed_accuracy),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
        'n': n,
        'se': float(np.std(bootstrap_accuracies))
    }


def bootstrap_top1_accuracy(
    scores: pd.DataFrame,
    score_column: str,
    truth_column: str = 'is_truth_gene',
    group_column: str = 'locus_id',
    n_bootstrap: int = 10000,
    confidence_level: float = 0.95
) -> Dict:
    """
    Calculate bootstrap CI for top-1 accuracy (does highest-scored gene match truth?).
    
    This is how L2G/cS2G accuracy is typically measured in gene prioritization.
    """
    # For each locus, check if top-scored gene is truth gene
    locus_correct = []
    
    for locus_id, group in scores.groupby(group_column):
        if group[score_column].isna().all():
            continue
        
        # Get top-scored gene
        top_gene_idx = group[score_column].idxmax()
        top_is_truth = group.loc[top_gene_idx, truth_column]
        locus_correct.append(int(top_is_truth))
    
    if len(locus_correct) == 0:
        return {'accuracy': np.nan, 'ci_lower': np.nan, 'ci_upper': np.nan, 'n': 0}
    
    correct = np.array(locus_correct)
    truths = np.ones(len(correct))  # All should be 1 for correct top-1
    
    # Bootstrap
    return bootstrap_accuracy(correct, truths, n_bootstrap, confidence_level)


def calculate_regulatory_benchmark_ci(data_dir: Path, n_bootstrap: int = 10000) -> Dict:
    """
    Calculate CI for expanded regulatory benchmark.
    
    This benchmark has 440 credible sets from ABC model silver-standard.
    """
    logger.info("Calculating CI for Expanded Regulatory Benchmark (n=440)...")
    
    accuracy_file = data_dir / "regulatory_benchmark_accuracy.tsv"
    if not accuracy_file.exists():
        logger.warning(f"File not found: {accuracy_file}")
        return {}
    
    df = pd.read_csv(accuracy_file, sep='\t')
    
    # The 'correct' column indicates if L2G top gene matched truth
    # Assuming 'correct' is 1 if L2G correctly predicted truth gene
    predictions = df['correct'].values.astype(int)
    truths = np.ones(len(predictions))  # All benchmark genes are true positives
    
    ci = bootstrap_accuracy(predictions, truths, n_bootstrap)
    
    logger.info(f"  Accuracy: {ci['accuracy']:.1%} [{ci['ci_lower']:.1%}, {ci['ci_upper']:.1%}]")
    
    return {
        'benchmark': 'Expanded Regulatory (ABC Silver Standard)',
        'n_loci': ci['n'],
        **ci
    }


def calculate_post2021_benchmark_ci(data_dir: Path, n_bootstrap: int = 10000) -> Dict:
    """
    Calculate CI for post-2021 independent benchmark (n=63 loci).
    """
    logger.info("Calculating CI for Post-2021 Independent Benchmark (n=63)...")
    
    # Load L2G scores
    l2g_file = data_dir / "l2g_benchmark_scores.tsv"
    if not l2g_file.exists():
        logger.warning(f"File not found: {l2g_file}")
        return {}
    
    l2g_df = pd.read_csv(l2g_file, sep='\t')
    
    # Filter to matched loci
    matched = l2g_df[l2g_df['l2g_score'].notna()].copy()
    
    # L2G prediction correct if score > 0.5
    predictions = (matched['l2g_score'] > 0.5).astype(int).values
    truths = np.ones(len(predictions))  # All are true positives
    
    ci = bootstrap_accuracy(predictions, truths, n_bootstrap)
    
    logger.info(f"  Coverage: {len(matched)}/{len(l2g_df)} ({len(matched)/len(l2g_df):.1%})")
    logger.info(f"  Accuracy: {ci['accuracy']:.1%} [{ci['ci_lower']:.1%}, {ci['ci_upper']:.1%}]")
    
    return {
        'benchmark': 'Post-2021 Independent',
        'n_total': len(l2g_df),
        'n_matched': len(matched),
        'coverage': len(matched) / len(l2g_df),
        **ci
    }


def calculate_cs2g_benchmark_ci(data_dir: Path, n_bootstrap: int = 10000) -> Dict:
    """
    Calculate CI for cS2G benchmark results.
    """
    logger.info("Calculating CI for cS2G Benchmark...")
    
    cs2g_file = data_dir / "cs2g_benchmark_liftover.tsv"
    if not cs2g_file.exists():
        logger.warning(f"File not found: {cs2g_file}")
        return {}
    
    cs2g_df = pd.read_csv(cs2g_file, sep='\t')
    
    # Filter to matched loci
    matched = cs2g_df[cs2g_df['cs2g_score'].notna()].copy()
    
    # cS2G prediction correct if score > 0.5
    predictions = (matched['cs2g_score'] > 0.5).astype(int).values
    truths = np.ones(len(predictions))  # All are true positives
    
    ci = bootstrap_accuracy(predictions, truths, n_bootstrap)
    
    logger.info(f"  Coverage: {len(matched)}/{len(cs2g_df)} ({len(matched)/len(cs2g_df):.1%})")
    logger.info(f"  Accuracy: {ci['accuracy']:.1%} [{ci['ci_lower']:.1%}, {ci['ci_upper']:.1%}]")
    
    return {
        'benchmark': 'cS2G (Gazal et al. 2022)',
        'n_total': len(cs2g_df),
        'n_matched': len(matched),
        'coverage': len(matched) / len(cs2g_df),
        **ci
    }


def compare_methods_bootstrap(data_dir: Path, n_bootstrap: int = 10000) -> Dict:
    """
    Statistical comparison between L2G and cS2G using paired bootstrap.
    """
    logger.info("Comparing L2G vs cS2G with paired bootstrap...")
    
    l2g_file = data_dir / "l2g_benchmark_scores.tsv"
    cs2g_file = data_dir / "cs2g_benchmark_liftover.tsv"
    
    if not l2g_file.exists() or not cs2g_file.exists():
        return {}
    
    l2g_df = pd.read_csv(l2g_file, sep='\t')
    cs2g_df = pd.read_csv(cs2g_file, sep='\t')
    
    # Merge on locus
    merged = l2g_df.merge(cs2g_df, on='locus_id', suffixes=('_l2g', '_cs2g'))
    
    # Filter to loci with both scores
    both_matched = merged[
        merged['l2g_score'].notna() & 
        merged['cs2g_score'].notna()
    ].copy()
    
    if len(both_matched) < 5:
        logger.warning(f"Too few paired observations: {len(both_matched)}")
        return {}
    
    # Calculate predictions
    l2g_correct = (both_matched['l2g_score'] > 0.5).astype(int).values
    cs2g_correct = (both_matched['cs2g_score'] > 0.5).astype(int).values
    
    # Bootstrap difference
    rng = np.random.default_rng(42)
    n = len(l2g_correct)
    diffs = []
    
    for _ in range(n_bootstrap):
        indices = rng.choice(n, size=n, replace=True)
        l2g_acc = np.mean(l2g_correct[indices])
        cs2g_acc = np.mean(cs2g_correct[indices])
        diffs.append(l2g_acc - cs2g_acc)
    
    diffs = np.array(diffs)
    
    observed_diff = np.mean(l2g_correct) - np.mean(cs2g_correct)
    ci_lower = np.percentile(diffs, 2.5)
    ci_upper = np.percentile(diffs, 97.5)
    
    # P-value approximation (proportion of bootstrap samples where diff <= 0)
    p_value = np.mean(diffs <= 0)
    
    logger.info(f"  Paired loci: {n}")
    logger.info(f"  L2G accuracy: {np.mean(l2g_correct):.1%}")
    logger.info(f"  cS2G accuracy: {np.mean(cs2g_correct):.1%}")
    logger.info(f"  Difference: {observed_diff:.1%} [{ci_lower:.1%}, {ci_upper:.1%}]")
    logger.info(f"  P-value (L2G > cS2G): {1-p_value:.3f}")
    
    return {
        'n_paired': n,
        'l2g_accuracy': float(np.mean(l2g_correct)),
        'cs2g_accuracy': float(np.mean(cs2g_correct)),
        'difference': float(observed_diff),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
        'p_value': float(p_value)
    }


def main():
    """Run bootstrap CI analysis for all benchmarks."""
    # Find project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    data_dir = project_root / "data" / "processed" / "baselines"
    
    n_bootstrap = 10000
    
    logger.info("=" * 70)
    logger.info("BOOTSTRAP CONFIDENCE INTERVALS FOR ALL BENCHMARKS")
    logger.info(f"Number of bootstrap samples: {n_bootstrap}")
    logger.info("=" * 70)
    
    results = {}
    
    # 1. Expanded Regulatory Benchmark
    results['regulatory'] = calculate_regulatory_benchmark_ci(data_dir, n_bootstrap)
    
    # 2. Post-2021 Independent Benchmark
    results['post2021'] = calculate_post2021_benchmark_ci(data_dir, n_bootstrap)
    
    # 3. cS2G Benchmark
    results['cs2g'] = calculate_cs2g_benchmark_ci(data_dir, n_bootstrap)
    
    # 4. Method comparison
    results['comparison'] = compare_methods_bootstrap(data_dir, n_bootstrap)
    
    # Save results
    output_file = data_dir / "bootstrap_confidence_intervals.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    logger.info(f"Saved results to {output_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("MANUSCRIPT SUMMARY: BOOTSTRAP 95% CONFIDENCE INTERVALS")
    print("=" * 80)
    
    print("\n## Main Results")
    print("-" * 60)
    
    if results.get('post2021'):
        r = results['post2021']
        print(f"\n**Post-2021 Independent Benchmark (L2G)**")
        print(f"  Coverage: {r.get('n_matched', 'N/A')}/{r.get('n_total', 'N/A')} ({r.get('coverage', 0):.1%})")
        print(f"  Accuracy: {r['accuracy']:.1%} (95% CI: {r['ci_lower']:.1%}–{r['ci_upper']:.1%})")
    
    if results.get('regulatory'):
        r = results['regulatory']
        print(f"\n**Expanded Regulatory Benchmark (L2G)**")
        print(f"  N: {r['n']}")
        print(f"  Accuracy: {r['accuracy']:.1%} (95% CI: {r['ci_lower']:.1%}–{r['ci_upper']:.1%})")
    
    if results.get('cs2g'):
        r = results['cs2g']
        print(f"\n**cS2G Baseline (Gazal et al. 2022)**")
        print(f"  Coverage: {r.get('n_matched', 'N/A')}/{r.get('n_total', 'N/A')} ({r.get('coverage', 0):.1%})")
        print(f"  Accuracy: {r['accuracy']:.1%} (95% CI: {r['ci_lower']:.1%}–{r['ci_upper']:.1%})")
    
    if results.get('comparison') and results['comparison'].get('n_paired', 0) > 0:
        r = results['comparison']
        print(f"\n**Method Comparison (paired loci)**")
        print(f"  N paired: {r['n_paired']}")
        print(f"  L2G: {r['l2g_accuracy']:.1%}")
        print(f"  cS2G: {r['cs2g_accuracy']:.1%}")
        print(f"  Difference: {r['difference']:.1%} (95% CI: {r['ci_lower']:.1%}–{r['ci_upper']:.1%})")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
