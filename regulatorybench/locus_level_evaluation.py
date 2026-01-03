#!/usr/bin/env python3
"""
Locus-Level Gene Prioritization Metrics

Peer reviewers correctly noted that pair-level AUROC is not the right metric
for gene prioritization. What matters is: at each GWAS locus, does the method
rank the true causal gene(s) above false candidate genes?

LOCUS-LEVEL METRICS:
- Top-1 Accuracy: Is the #1 predicted gene the causal gene?
- Top-3 Accuracy: Is the causal gene in the top 3 predictions?
- Mean Reciprocal Rank (MRR): Average of 1/rank of causal gene
- Mean Average Precision (MAP): mAP across loci

BOOTSTRAP CONFIDENCE INTERVALS:
- Locus-level bootstrap (not pair-level)
- 1000 resamples
- 95% percentile CIs

Author: RegulatoryBench v2.0
Date: December 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from typing import Dict, List, Tuple
from collections import defaultdict
import warnings
from scipy import stats

warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).parent.resolve()


class LocusLevelEvaluator:
    """Evaluate gene prioritization at locus level, not pair level."""
    
    def __init__(self, benchmark_v2: bool = True):
        self.benchmarks_dir = SCRIPT_DIR / "benchmarks"
        self.results_dir = SCRIPT_DIR / "results"
        self.use_v2 = benchmark_v2
        
    def load_benchmark(self) -> pd.DataFrame:
        """Load Task A benchmark (v2 with independent labels)."""
        if self.use_v2:
            path = self.benchmarks_dir / "task_a_gwas_to_gene_v2.parquet"
        else:
            path = self.benchmarks_dir / "task_a_gwas_to_gene.parquet"
        
        df = pd.read_parquet(path)
        print(f"Loaded benchmark: {len(df)} pairs, {df.locus_id.nunique()} loci, "
              f"{df.is_positive.sum()} positives ({df.is_positive.mean():.2%})")
        return df
    
    def compute_locus_ranks(self, df: pd.DataFrame, score_col: str) -> pd.DataFrame:
        """
        For each locus, rank genes by prediction score.
        
        Returns DataFrame with locus_id, gene_symbol, rank, is_positive.
        """
        # Rank within each locus (higher score = better rank)
        df = df.copy()
        
        # Special handling for distance - lower is better
        if 'distance' in score_col.lower():
            df['rank'] = df.groupby('locus_id')[score_col].rank(ascending=True, method='min')
        else:
            df['rank'] = df.groupby('locus_id')[score_col].rank(ascending=False, method='min')
        
        return df[['locus_id', 'gene_symbol', score_col, 'rank', 'is_positive']]
    
    def top_k_accuracy(self, ranks_df: pd.DataFrame, k: int) -> float:
        """
        Top-k accuracy: fraction of loci where ≥1 causal gene is in top k.
        
        For loci with no positives, ignore (test set loci only).
        """
        locus_results = []
        
        for locus_id, locus_df in ranks_df.groupby('locus_id'):
            positives = locus_df[locus_df['is_positive'] == True]
            
            if len(positives) == 0:
                continue  # Locus has no causal genes (negative locus)
            
            # Is any causal gene in top k?
            top_k_genes = locus_df.nsmallest(k, 'rank')
            has_causal_in_top_k = (top_k_genes['is_positive'] == True).any()
            locus_results.append(has_causal_in_top_k)
        
        return np.mean(locus_results) if locus_results else 0.0
    
    def mean_reciprocal_rank(self, ranks_df: pd.DataFrame) -> float:
        """
        MRR: Average of 1/rank for causal genes.
        
        If a locus has multiple causal genes, take best (min) rank.
        """
        reciprocal_ranks = []
        
        for locus_id, locus_df in ranks_df.groupby('locus_id'):
            positives = locus_df[locus_df['is_positive'] == True]
            
            if len(positives) == 0:
                continue
            
            # Best rank among causal genes
            best_rank = positives['rank'].min()
            reciprocal_ranks.append(1.0 / best_rank)
        
        return np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
    
    def mean_average_precision(self, ranks_df: pd.DataFrame) -> float:
        """
        MAP: Mean average precision across loci.
        
        AP = mean precision at each recall point (each true positive).
        """
        aps = []
        
        for locus_id, locus_df in ranks_df.groupby('locus_id'):
            # Sort by rank
            locus_df = locus_df.sort_values('rank')
            
            positives = locus_df['is_positive'].values
            n_positives = positives.sum()
            
            if n_positives == 0:
                continue
            
            # Precision at each position where we find a positive
            precisions = []
            n_correct = 0
            for i, is_pos in enumerate(positives, start=1):
                if is_pos:
                    n_correct += 1
                    precisions.append(n_correct / i)
            
            # Average precision for this locus
            ap = np.mean(precisions) if precisions else 0.0
            aps.append(ap)
        
        return np.mean(aps) if aps else 0.0
    
    def evaluate_method(self, df: pd.DataFrame, score_col: str, 
                       method_name: str) -> Dict:
        """Compute all locus-level metrics for a method."""
        # Compute ranks
        ranks_df = self.compute_locus_ranks(df, score_col)
        
        # Compute metrics
        metrics = {
            'method': method_name,
            'score_column': score_col,
            'top1_accuracy': self.top_k_accuracy(ranks_df, k=1),
            'top3_accuracy': self.top_k_accuracy(ranks_df, k=3),
            'top5_accuracy': self.top_k_accuracy(ranks_df, k=5),
            'mrr': self.mean_reciprocal_rank(ranks_df),
            'map': self.mean_average_precision(ranks_df),
            'n_loci_evaluated': ranks_df[ranks_df['is_positive'] == True]['locus_id'].nunique()
        }
        
        return metrics
    
    def bootstrap_locus_level(self, df: pd.DataFrame, score_col: str,
                              n_bootstrap: int = 1000) -> Dict:
        """
        Bootstrap locus-level metrics by resampling loci (not pairs).
        
        This gives proper confidence intervals for locus-level performance.
        """
        # Get loci with positives
        positive_loci = df[df['is_positive'] == True]['locus_id'].unique()
        n_loci = len(positive_loci)
        
        print(f"Bootstrapping {n_bootstrap} samples from {n_loci} loci...")
        
        bootstrap_results = {
            'top1': [],
            'top3': [],
            'top5': [],
            'mrr': [],
            'map': []
        }
        
        for i in range(n_bootstrap):
            # Resample loci with replacement
            resampled_loci = np.random.choice(positive_loci, size=n_loci, replace=True)
            resampled_df = df[df['locus_id'].isin(resampled_loci)]
            
            # Compute metrics
            ranks_df = self.compute_locus_ranks(resampled_df, score_col)
            
            bootstrap_results['top1'].append(self.top_k_accuracy(ranks_df, k=1))
            bootstrap_results['top3'].append(self.top_k_accuracy(ranks_df, k=3))
            bootstrap_results['top5'].append(self.top_k_accuracy(ranks_df, k=5))
            bootstrap_results['mrr'].append(self.mean_reciprocal_rank(ranks_df))
            bootstrap_results['map'].append(self.mean_average_precision(ranks_df))
        
        # Compute 95% CIs
        cis = {}
        for metric, values in bootstrap_results.items():
            cis[f'{metric}_ci_lower'] = np.percentile(values, 2.5)
            cis[f'{metric}_ci_upper'] = np.percentile(values, 97.5)
            cis[f'{metric}_std'] = np.std(values)
        
        return cis
    
    def compare_methods(self, df: pd.DataFrame, 
                       methods: List[Tuple[str, str]]) -> pd.DataFrame:
        """
        Compare multiple methods with locus-level metrics.
        
        Args:
            methods: List of (score_column, method_name) tuples
        """
        results = []
        
        for score_col, method_name in methods:
            print(f"\nEvaluating {method_name}...")
            
            # Main metrics
            metrics = self.evaluate_method(df, score_col, method_name)
            
            # Bootstrap CIs
            print(f"  Computing bootstrap CIs...")
            cis = self.bootstrap_locus_level(df, score_col, n_bootstrap=1000)
            
            # Combine
            metrics.update(cis)
            results.append(metrics)
        
        results_df = pd.DataFrame(results)
        
        # Sort by MRR (best overall metric)
        results_df = results_df.sort_values('mrr', ascending=False)
        
        return results_df
    
    def save_results(self, results_df: pd.DataFrame, output_name: str = "locus_level_results"):
        """Save locus-level evaluation results."""
        output_path = self.results_dir / f"{output_name}.parquet"
        results_df.to_parquet(output_path, index=False)
        print(f"\nSaved results: {output_path}")
        
        # Also save human-readable version
        csv_path = self.results_dir / f"{output_name}.csv"
        results_df.to_csv(csv_path, index=False)
        print(f"Saved CSV: {csv_path}")
    
    def print_comparison_table(self, results_df: pd.DataFrame):
        """Print formatted comparison table."""
        print("\n" + "="*100)
        print("LOCUS-LEVEL GENE PRIORITIZATION PERFORMANCE")
        print("="*100)
        print(f"{'Method':<20} {'Top-1':>8} {'Top-3':>8} {'Top-5':>8} {'MRR':>8} {'MAP':>8} {'N Loci':>8}")
        print("-"*100)
        
        for _, row in results_df.iterrows():
            print(f"{row['method']:<20} "
                  f"{row['top1_accuracy']:>8.3f} "
                  f"{row['top3_accuracy']:>8.3f} "
                  f"{row['top5_accuracy']:>8.3f} "
                  f"{row['mrr']:>8.3f} "
                  f"{row['map']:>8.3f} "
                  f"{row['n_loci_evaluated']:>8.0f}")
            
            # Print CIs
            print(f"{'  95% CI':<20} "
                  f"±{row['top1_std']:>7.3f} "
                  f"±{row['top3_std']:>7.3f} "
                  f"±{row['top5_std']:>7.3f} "
                  f"±{row['mrr_std']:>7.3f} "
                  f"±{row['map_std']:>7.3f}")
        
        print("="*100)


def main():
    """Run locus-level evaluation on Task A v2 benchmark."""
    evaluator = LocusLevelEvaluator(benchmark_v2=True)
    df = evaluator.load_benchmark()
    
    # Define methods to compare
    methods = [
        ('MaxABC', 'ABC'),
        ('POPS.Score', 'PoPS'),
        ('GeneBodyDistanceToBestSNP', 'Distance Baseline'),
    ]
    
    # Remove methods without scores
    available_methods = [(col, name) for col, name in methods if col in df.columns]
    print(f"\nEvaluating {len(available_methods)} methods: {[m[1] for m in available_methods]}")
    
    # Run comparison
    results_df = evaluator.compare_methods(df, available_methods)
    
    # Print results
    evaluator.print_comparison_table(results_df)
    
    # Save results
    evaluator.save_results(results_df, output_name="locus_level_results_v2")
    
    print("\n" + "="*100)
    print("LOCUS-LEVEL EVALUATION COMPLETE!")
    print("="*100)


if __name__ == '__main__':
    main()
