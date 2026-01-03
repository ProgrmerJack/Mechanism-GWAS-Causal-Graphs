#!/usr/bin/env python3
"""
Distance-Stratified Evaluation for Regime Map Analysis

Evaluates GWAS-to-gene methods across distance bins to show where
proximity baseline is sufficient vs where integration adds value.

Key insight for Nature Genetics: Performance varies by TSS-to-SNP distance

Author: Generated for Nature Genetics v3 Platinum manuscript  
Date: 2025-12-21
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class DistanceStratifiedEvaluator:
    """Evaluate methods across distance regimes."""
    
    # Distance bins (TSS to lead SNP)
    DISTANCE_BINS = [
        (0, 10_000, "0-10kb (proximal)"),
        (10_000, 100_000, "10-100kb (near)"),
        (100_000, 500_000, "100-500kb (distal)"),
        (500_000, float('inf'), ">500kb (very distal)")
    ]
    
    def __init__(self, benchmark_path: str, results_dir: str = "benchmarks/results"):
        """
        Initialize evaluator.
        
        Parameters
        ----------
        benchmark_path : str
            Path to benchmark parquet file with predictions
        results_dir : str
            Output directory for results
        """
        self.benchmark_path = Path(benchmark_path)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        self.df = None
        self.label_col = None
        self.method_cols = []
        
    def load_data(self) -> None:
        """Load benchmark with predictions."""
        logger.info(f"Loading benchmark from {self.benchmark_path}")
        self.df = pd.read_parquet(self.benchmark_path)
        logger.info(f"Loaded {len(self.df)} gene-locus pairs")
        
        # Identify label column
        label_candidates = [c for c in self.df.columns if 'label' in c.lower() or 'truth' in c.lower()]
        if label_candidates:
            self.label_col = label_candidates[0]
        else:
            # Use TruthDistanceRank as binary label
            self.label_col = 'is_causal'
            self.df['is_causal'] = self.df['TruthDistanceRank'] == 1
            
        logger.info(f"Using label column: {self.label_col}")
        logger.info(f"Positive labels: {self.df[self.label_col].sum()}")
        
        # Identify method prediction columns (rank-based)
        rank_cols = [c for c in self.df.columns if 'Rank' in c and c != 'TruthDistanceRank' and c != 'DistanceRank']
        self.method_cols = rank_cols[:10]  # Top 10 methods
        logger.info(f"Found {len(self.method_cols)} prediction methods")
        
        # Calculate TSS-to-SNP distances if not present
        if 'tss_distance' not in self.df.columns:
            # Use GeneBodyDistanceToBestSNP or PromoterDistanceToBestSNP
            if 'GeneBodyDistanceToBestSNP' in self.df.columns:
                self.df['tss_distance'] = self.df['GeneBodyDistanceToBestSNP'].abs()
            elif 'PromoterDistanceToBestSNP' in self.df.columns:
                self.df['tss_distance'] = self.df['PromoterDistanceToBestSNP'].abs()
            else:
                logger.error("No distance column found")
                raise ValueError("Need TSS-to-SNP distance information")
                
        logger.info(f"Distance range: {self.df['tss_distance'].min():.0f} to {self.df['tss_distance'].max():.0f} bp")
        
    def assign_distance_bins(self) -> None:
        """Assign each gene-locus pair to distance bin."""
        self.df['distance_bin'] = None
        self.df['distance_label'] = None
        
        for min_dist, max_dist, label in self.DISTANCE_BINS:
            mask = (self.df['tss_distance'] >= min_dist) & (self.df['tss_distance'] < max_dist)
            self.df.loc[mask, 'distance_bin'] = f"{min_dist}-{max_dist}"
            self.df.loc[mask, 'distance_label'] = label
            
        logger.info("\nDistance bin distribution:")
        for label in self.df['distance_label'].value_counts().items():
            logger.info(f"  {label[0]}: {label[1]} pairs")
            
    def calculate_topk_precision(self, method_col: str, k: int = 1) -> Dict[str, float]:
        """
        Calculate Top-K precision per distance bin.
        
        Parameters
        ----------
        method_col : str
            Method prediction column (rank-based)
        k : int
            Top-K threshold
            
        Returns
        -------
        dict
            {distance_label: top-k precision}
        """
        results = {}
        
        for _, _, label in self.DISTANCE_BINS:
            bin_mask = self.df['distance_label'] == label
            bin_df = self.df[bin_mask]
            
            if len(bin_df) == 0:
                results[label] = np.nan
                continue
                
            # For each locus, check if top-k prediction is correct
            topk_correct = []
            for locus_id, locus_df in bin_df.groupby('Disease'):
                # Get true causal genes
                true_genes = set(locus_df[locus_df[self.label_col]]['TargetGene'])
                
                if len(true_genes) == 0:
                    continue
                    
                # Get top-k predicted genes
                topk_pred = locus_df.nsmallest(k, method_col)['TargetGene']
                
                # Check if any true gene in top-k
                correct = any(gene in true_genes for gene in topk_pred)
                topk_correct.append(correct)
                
            results[label] = np.mean(topk_correct) if topk_correct else np.nan
            
        return results
        
    def evaluate_all_methods(self) -> pd.DataFrame:
        """
        Evaluate all methods across distance bins.
        
        Returns
        -------
        DataFrame
            Results with columns: method, distance_bin, top1_precision, top3_precision, top5_precision
        """
        logger.info("\nEvaluating all methods across distance bins...")
        
        results = []
        
        for method_col in self.method_cols:
            method_name = method_col.replace('Rank', '').replace('ConnectionStrengthRank', 'ABC')
            
            # Calculate Top-1, Top-3, Top-5
            top1 = self.calculate_topk_precision(method_col, k=1)
            top3 = self.calculate_topk_precision(method_col, k=3)
            top5 = self.calculate_topk_precision(method_col, k=5)
            
            for dist_label in [label for _, _, label in self.DISTANCE_BINS]:
                results.append({
                    'method': method_name,
                    'distance_bin': dist_label,
                    'top1_precision': top1.get(dist_label, np.nan),
                    'top3_precision': top3.get(dist_label, np.nan),
                    'top5_precision': top5.get(dist_label, np.nan)
                })
                
        results_df = pd.DataFrame(results)
        
        # Add baseline (proximity = DistanceRank)
        if 'DistanceRank' in self.df.columns:
            top1_prox = self.calculate_topk_precision('DistanceRank', k=1)
            top3_prox = self.calculate_topk_precision('DistanceRank', k=3)
            top5_prox = self.calculate_topk_precision('DistanceRank', k=5)
            
            for dist_label in [label for _, _, label in self.DISTANCE_BINS]:
                results_df = pd.concat([results_df, pd.DataFrame([{
                    'method': 'Proximity (baseline)',
                    'distance_bin': dist_label,
                    'top1_precision': top1_prox.get(dist_label, np.nan),
                    'top3_precision': top3_prox.get(dist_label, np.nan),
                    'top5_precision': top5_prox.get(dist_label, np.nan)
                }])], ignore_index=True)
                
        return results_df
        
    def compute_delta_from_baseline(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """Add delta from proximity baseline."""
        baseline_df = results_df[results_df['method'] == 'Proximity (baseline)']
        
        delta_rows = []
        for _, row in results_df.iterrows():
            if row['method'] == 'Proximity (baseline)':
                continue
                
            baseline_row = baseline_df[baseline_df['distance_bin'] == row['distance_bin']]
            if len(baseline_row) == 0:
                continue
                
            delta_rows.append({
                'method': row['method'],
                'distance_bin': row['distance_bin'],
                'delta_top1': row['top1_precision'] - baseline_row['top1_precision'].values[0],
                'delta_top3': row['top3_precision'] - baseline_row['top3_precision'].values[0],
                'delta_top5': row['top5_precision'] - baseline_row['top5_precision'].values[0]
            })
            
        delta_df = pd.DataFrame(delta_rows)
        
        # Merge with original results
        results_with_delta = results_df.merge(delta_df, on=['method', 'distance_bin'], how='left')
        
        return results_with_delta
        
    def run(self) -> pd.DataFrame:
        """Execute full distance-stratified evaluation."""
        logger.info("="*60)
        logger.info("Distance-Stratified Evaluation (Regime Map)")
        logger.info("="*60)
        
        self.load_data()
        self.assign_distance_bins()
        results_df = self.evaluate_all_methods()
        results_with_delta = self.compute_delta_from_baseline(results_df)
        
        # Save results
        output_file = self.results_dir / "distance_stratified_performance.csv"
        results_with_delta.to_csv(output_file, index=False)
        logger.info(f"\nSaved results to: {output_file}")
        
        # Print summary
        logger.info("\n" + "="*60)
        logger.info("DISTANCE REGIME SUMMARY")
        logger.info("="*60)
        
        for dist_label in [label for _, _, label in self.DISTANCE_BINS]:
            logger.info(f"\n{dist_label}:")
            bin_results = results_with_delta[results_with_delta['distance_bin'] == dist_label]
            
            logger.info("  Top-1 Precision:")
            for _, row in bin_results.nsmallest(3, 'top1_precision', keep='all').iterrows():
                logger.info(f"    {row['method']}: {row['top1_precision']:.3f}")
                
        return results_with_delta


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Distance-stratified evaluation")
    parser.add_argument("--benchmark", default="benchmarks/task_a_gwas_to_gene_v3_platinum.parquet",
                       help="Path to benchmark file")
    parser.add_argument("--results-dir", default="benchmarks/results",
                       help="Output directory")
    
    args = parser.parse_args()
    
    evaluator = DistanceStratifiedEvaluator(
        benchmark_path=args.benchmark,
        results_dir=args.results_dir
    )
    
    results = evaluator.run()
    print("\nDistance stratification complete!")
