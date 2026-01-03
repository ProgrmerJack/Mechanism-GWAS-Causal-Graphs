#!/usr/bin/env python3
"""
Task B Baseline Evaluation

Evaluates multiple baselines on the Task B (Enhancer→Gene) benchmark:
1. Distance to TSS (baseline)
2. ABC model (direct from benchmark data)
3. FLAMES ABC features (where available)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
from dataclasses import dataclass
from typing import Dict, List, Optional
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Use script directory for relative paths
SCRIPT_DIR = Path(__file__).parent.resolve()


@dataclass
class MethodResult:
    """Results for a single method evaluation."""
    method_name: str
    task_type: str
    auroc: float
    auprc: float
    coverage: float
    n_pairs_scored: int
    n_pairs_total: int
    source_breakdown: Optional[Dict[str, float]] = None


def evaluate_distance_baseline(bench: pd.DataFrame) -> MethodResult:
    """
    Evaluate distance-to-TSS baseline.
    
    Uses negative distance (closer = higher score) as prediction.
    """
    logger.info("Evaluating Distance-to-TSS baseline...")
    
    if 'distanceToTSS' not in bench.columns:
        # Calculate distance from enhancer to TSS
        bench['distanceToTSS'] = np.abs(
            ((bench['chromStart'] + bench['chromEnd']) / 2) - bench['startTSS']
        )
    
    # Filter to pairs with valid distances and labels
    valid = bench[bench['distanceToTSS'].notna() & bench['is_positive'].notna()].copy()
    
    if len(valid) < 10:
        return MethodResult(
            method_name='Distance_to_TSS',
            task_type='B',
            auroc=np.nan,
            auprc=np.nan,
            coverage=0.0,
            n_pairs_scored=0,
            n_pairs_total=len(bench)
        )
    
    # Score: negative distance (closer = better)
    scores = -valid['distanceToTSS'].values
    labels = valid['is_positive'].astype(int).values
    
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    coverage = len(valid) / len(bench)
    
    # Breakdown by source
    source_breakdown = {}
    for source in valid['source'].unique():
        subset = valid[valid['source'] == source]
        if len(subset) > 10 and subset['is_positive'].sum() > 0:
            s = -subset['distanceToTSS'].values
            l = subset['is_positive'].astype(int).values
            try:
                source_breakdown[source] = roc_auc_score(l, s)
            except:
                pass
    
    logger.info(f"Distance baseline: AUROC={auroc:.3f}, AUPRC={auprc:.3f}, coverage={coverage:.1%}")
    
    return MethodResult(
        method_name='Distance_to_TSS',
        task_type='B',
        auroc=auroc,
        auprc=auprc,
        coverage=coverage,
        n_pairs_scored=len(valid),
        n_pairs_total=len(bench),
        source_breakdown=source_breakdown
    )


def evaluate_abc_from_fulco(bench: pd.DataFrame) -> MethodResult:
    """
    Evaluate ABC scores directly from Fulco 2019 data.
    
    Fulco data includes pre-computed ABC scores.
    """
    logger.info("Evaluating ABC scores (Fulco 2019 subset)...")
    
    # Filter to Fulco data which has ABC scores
    fulco = bench[bench['source'] == 'Fulco_2019'].copy()
    
    if 'ABC_score' not in fulco.columns or len(fulco) == 0:
        logger.warning("No ABC scores available in Fulco data")
        return MethodResult(
            method_name='ABC_Fulco2019',
            task_type='B',
            auroc=np.nan,
            auprc=np.nan,
            coverage=0.0,
            n_pairs_scored=0,
            n_pairs_total=len(bench)
        )
    
    # Filter to valid pairs
    valid = fulco[fulco['ABC_score'].notna() & fulco['is_positive'].notna()].copy()
    
    if len(valid) < 10 or valid['is_positive'].sum() == 0:
        logger.warning(f"Insufficient valid pairs: {len(valid)}, positives: {valid['is_positive'].sum() if len(valid) > 0 else 0}")
        return MethodResult(
            method_name='ABC_Fulco2019',
            task_type='B',
            auroc=np.nan,
            auprc=np.nan,
            coverage=0.0,
            n_pairs_scored=len(valid),
            n_pairs_total=len(bench)
        )
    
    scores = valid['ABC_score'].values
    labels = valid['is_positive'].astype(int).values
    
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    coverage = len(valid) / len(bench)
    
    logger.info(f"ABC (Fulco 2019): AUROC={auroc:.3f}, AUPRC={auprc:.3f}, coverage={coverage:.1%}")
    
    return MethodResult(
        method_name='ABC_Fulco2019',
        task_type='B',
        auroc=auroc,
        auprc=auprc,
        coverage=coverage,
        n_pairs_scored=len(valid),
        n_pairs_total=len(bench)
    )


def evaluate_by_distance_strata(bench: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """
    Evaluate distance baseline stratified by distance ranges.
    
    Ranges: 0-10kb, 10-100kb, 100-500kb
    """
    logger.info("Computing distance-stratified metrics...")
    
    if 'distanceToTSS' not in bench.columns:
        bench['distanceToTSS'] = np.abs(
            ((bench['chromStart'] + bench['chromEnd']) / 2) - bench['startTSS']
        )
    
    strata = {
        '0-10kb': (0, 10000),
        '10-100kb': (10000, 100000),
        '100-500kb': (100000, 500000)
    }
    
    results = {}
    for stratum_name, (min_dist, max_dist) in strata.items():
        subset = bench[
            (bench['distanceToTSS'] >= min_dist) & 
            (bench['distanceToTSS'] < max_dist) &
            bench['is_positive'].notna()
        ]
        
        if len(subset) < 20 or subset['is_positive'].sum() == 0:
            results[stratum_name] = {
                'n_pairs': len(subset),
                'n_positive': int(subset['is_positive'].sum()) if len(subset) > 0 else 0,
                'auroc': np.nan,
                'auprc': np.nan
            }
            continue
            
        scores = -subset['distanceToTSS'].values
        labels = subset['is_positive'].astype(int).values
        
        try:
            auroc = roc_auc_score(labels, scores)
            auprc = average_precision_score(labels, scores)
        except:
            auroc, auprc = np.nan, np.nan
        
        results[stratum_name] = {
            'n_pairs': len(subset),
            'n_positive': int(subset['is_positive'].sum()),
            'auroc': auroc,
            'auprc': auprc
        }
        
        logger.info(f"  {stratum_name}: n={len(subset)}, pos={int(subset['is_positive'].sum())}, AUC={auroc:.3f}")
    
    return results


def main():
    """Run Task B baseline evaluation."""
    
    benchmark_path = SCRIPT_DIR / "benchmarks/task_b_enhancer_to_gene.parquet"
    if not benchmark_path.exists():
        logger.error(f"Benchmark not found: {benchmark_path}")
        return
    
    logger.info(f"Loading benchmark: {benchmark_path}")
    bench = pd.read_parquet(benchmark_path)
    
    # Ensure label column
    if 'is_positive' not in bench.columns and 'Regulated' in bench.columns:
        bench['is_positive'] = bench['Regulated']
    
    logger.info(f"Benchmark size: {len(bench)} pairs")
    logger.info(f"Positives: {bench['is_positive'].sum()} ({bench['is_positive'].mean():.1%})")
    logger.info(f"Sources: {bench['source'].unique().tolist()}")
    
    results = []
    
    # 1. Distance baseline
    dist_result = evaluate_distance_baseline(bench)
    results.append(dist_result)
    
    # 2. ABC from Fulco
    abc_result = evaluate_abc_from_fulco(bench)
    results.append(abc_result)
    
    # 3. Distance-stratified analysis
    strata_results = evaluate_by_distance_strata(bench)
    
    # Print summary
    print("\n" + "=" * 70)
    print("TASK B (ENHANCER → GENE) BASELINE EVALUATION")
    print("=" * 70)
    
    print("\n### Method Comparison ###\n")
    print(f"{'Method':<25} {'AUROC':>8} {'AUPRC':>8} {'Coverage':>10} {'N Scored':>12}")
    print("-" * 70)
    for r in results:
        auroc_str = f"{r.auroc:.3f}" if not np.isnan(r.auroc) else "N/A"
        auprc_str = f"{r.auprc:.3f}" if not np.isnan(r.auprc) else "N/A"
        print(f"{r.method_name:<25} {auroc_str:>8} {auprc_str:>8} {r.coverage:>9.1%} {r.n_pairs_scored:>12,}")
    
    if dist_result.source_breakdown:
        print("\n### Distance Baseline by Source ###\n")
        for source, auroc in dist_result.source_breakdown.items():
            print(f"  {source}: AUROC={auroc:.3f}")
    
    print("\n### Distance-Stratified Analysis ###\n")
    for stratum, metrics in strata_results.items():
        auroc_str = f"{metrics['auroc']:.3f}" if not np.isnan(metrics['auroc']) else "N/A"
        print(f"  {stratum}: n={metrics['n_pairs']}, pos={metrics['n_positive']}, AUROC={auroc_str}")
    
    # Save results
    output = {
        'methods': [
            {
                'name': r.method_name,
                'task': r.task_type,
                'auroc': float(r.auroc) if not np.isnan(r.auroc) else None,
                'auprc': float(r.auprc) if not np.isnan(r.auprc) else None,
                'coverage': float(r.coverage),
                'n_scored': r.n_pairs_scored,
                'n_total': r.n_pairs_total
            }
            for r in results
        ],
        'distance_strata': {k: {kk: float(vv) if isinstance(vv, (float, np.floating)) and not np.isnan(vv) else vv 
                                for kk, vv in v.items()} 
                           for k, v in strata_results.items()},
        'benchmark_info': {
            'n_pairs': len(bench),
            'n_positive': int(bench['is_positive'].sum()),
            'positive_rate': float(bench['is_positive'].mean()),
            'sources': bench['source'].unique().tolist()
        }
    }
    
    output_path = SCRIPT_DIR / "benchmarks/task_b_baseline_results.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    logger.info(f"Saved results to {output_path}")


if __name__ == "__main__":
    main()
