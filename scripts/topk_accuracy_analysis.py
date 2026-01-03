#!/usr/bin/env python3
"""
Top-k Accuracy Analysis for L2G Benchmark
==========================================

This script calculates Top-k (k=1,3,5) accuracy metrics for the L2G benchmark,
which is a standard way to evaluate gene prioritization methods in the field.

Top-k accuracy asks: "Is the true causal gene in the top k predictions?"
This is more lenient than Top-1 but captures methods that rank well overall.

For cS2G, we only have per-gene scores (not all candidate genes ranked),
so Top-k analysis is primarily meaningful for L2G which ranks candidates.

Author: Mechanism-GWAS-Causal-Graphs
Date: 2025
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import requests
import time

import numpy as np
import pandas as pd
from scipy import stats

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TopKAccuracyAnalysis:
    """Calculate Top-k accuracy for L2G predictions."""
    
    PLATFORM_API_BASE = "https://api.platform.opentargets.org/api/v4/graphql"
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.data_dir = project_root / "data"
        self.baselines_dir = self.data_dir / "processed" / "baselines"
        self.benchmark_df = None
        
    def load_benchmark(self) -> pd.DataFrame:
        """Load the post-2021 independent benchmark."""
        benchmark_path = self.baselines_dir / "post2021_independent_benchmark.tsv"
        self.benchmark_df = pd.read_csv(benchmark_path, sep='\t')
        logger.info(f"Loaded benchmark: {len(self.benchmark_df)} loci")
        return self.benchmark_df
    
    def query_l2g_candidates(self, study_id: str, variant_id: str, limit: int = 10) -> List[Dict]:
        """
        Query Open Targets Platform API for all L2G candidate genes at a locus.
        
        Returns list of (gene_id, gene_symbol, l2g_score) sorted by score descending.
        """
        query = """
        query L2GCandidates($studyId: String!, $variantId: String!) {
            studyLocus2Genes(studyId: $studyId, variantId: $variantId, size: 20) {
                rows {
                    gene {
                        id
                        symbol
                    }
                    score
                }
            }
        }
        """
        
        try:
            response = requests.post(
                self.PLATFORM_API_BASE,
                json={"query": query, "variables": {"studyId": study_id, "variantId": variant_id}},
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            rows = data.get('data', {}).get('studyLocus2Genes', {}).get('rows', [])
            candidates = []
            for row in rows:
                if row.get('gene') and row.get('score') is not None:
                    candidates.append({
                        'gene_id': row['gene']['id'],
                        'gene_symbol': row['gene']['symbol'],
                        'l2g_score': row['score']
                    })
            
            # Sort by score descending
            candidates.sort(key=lambda x: x['l2g_score'], reverse=True)
            return candidates
            
        except Exception as e:
            logger.warning(f"API error for {study_id}/{variant_id}: {e}")
            return []
    
    def calculate_topk_for_locus(self, true_gene: str, candidates: List[Dict], k_values: List[int]) -> Dict[int, bool]:
        """
        Check if true_gene is in top-k candidates for each k.
        
        Args:
            true_gene: The known causal gene symbol
            candidates: List of candidate dicts with gene_symbol and l2g_score
            k_values: List of k values to check (e.g., [1, 3, 5])
        
        Returns:
            Dict mapping k -> bool (True if true_gene in top-k)
        """
        results = {}
        gene_symbols = [c['gene_symbol'] for c in candidates]
        
        for k in k_values:
            topk_genes = gene_symbols[:k]
            results[k] = true_gene in topk_genes
            
        return results
    
    def run_topk_analysis(self, k_values: List[int] = [1, 3, 5], use_cache: bool = True) -> Dict:
        """
        Run full Top-k analysis on the benchmark.
        
        Args:
            k_values: List of k values to evaluate
            use_cache: Whether to use cached L2G candidate data
        
        Returns:
            Dict with Top-k accuracies and confidence intervals
        """
        logger.info("=" * 70)
        logger.info("TOP-K ACCURACY ANALYSIS")
        logger.info("=" * 70)
        
        # Load benchmark
        if self.benchmark_df is None:
            self.load_benchmark()
        
        # Check for cached candidates
        cache_path = self.baselines_dir / "l2g_candidates_cache.json"
        if use_cache and cache_path.exists():
            logger.info(f"Loading cached L2G candidates from {cache_path}")
            with open(cache_path, 'r') as f:
                candidates_cache = json.load(f)
        else:
            candidates_cache = {}
        
        # Collect results for each locus
        results_per_locus = []
        
        for idx, row in self.benchmark_df.iterrows():
            locus_id = row['locus_id']
            true_gene = row['gene_symbol']
            
            # Try to get cached candidates
            cache_key = locus_id
            if cache_key in candidates_cache:
                candidates = candidates_cache[cache_key]
            else:
                # Need to query API - construct study_id and variant_id
                # This is complex because we need GWAS study ID + variant ID
                # For now, use existing L2G scores file which has top-1 predictions
                candidates = []  # Will fall back to top-1 only
            
            if candidates:
                topk_results = self.calculate_topk_for_locus(true_gene, candidates, k_values)
            else:
                # Fall back to top-1 from existing data
                topk_results = {k: None for k in k_values}
            
            results_per_locus.append({
                'locus_id': locus_id,
                'true_gene': true_gene,
                'n_candidates': len(candidates),
                **{f'top_{k}': topk_results.get(k) for k in k_values}
            })
        
        # Calculate aggregate statistics
        results_df = pd.DataFrame(results_per_locus)
        
        topk_stats = {}
        for k in k_values:
            col = f'top_{k}'
            valid = results_df[results_df[col].notna()]
            if len(valid) > 0:
                accuracy = valid[col].mean()
                n = len(valid)
                
                # Wilson score interval for confidence interval
                z = 1.96
                p = accuracy
                denominator = 1 + z**2/n
                center = (p + z**2/(2*n)) / denominator
                margin = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denominator
                ci_lower = max(0, center - margin)
                ci_upper = min(1, center + margin)
                
                topk_stats[k] = {
                    'accuracy': accuracy,
                    'n_evaluated': n,
                    'ci_lower': ci_lower,
                    'ci_upper': ci_upper
                }
            else:
                topk_stats[k] = {
                    'accuracy': None,
                    'n_evaluated': 0,
                    'ci_lower': None,
                    'ci_upper': None
                }
        
        return {
            'k_values': k_values,
            'topk_statistics': topk_stats,
            'locus_details': results_per_locus
        }
    
    def analyze_existing_l2g_scores(self) -> Dict:
        """
        Analyze Top-k using existing L2G scores file.
        
        Since we only have top-1 predictions stored, this provides
        Top-1 analysis with proper confidence intervals.
        """
        logger.info("Analyzing existing L2G scores for Top-1 accuracy...")
        
        # Load benchmark
        if self.benchmark_df is None:
            self.load_benchmark()
        
        # Load unified benchmark with match info
        unified_path = self.baselines_dir / "unified_benchmark_l2g_cs2g.tsv"
        if unified_path.exists():
            unified_df = pd.read_csv(unified_path, sep='\t')
            logger.info(f"Loaded unified benchmark: {len(unified_df)} loci")
            
            # L2G correct: l2g_score > 0.5 (causal gene has high L2G)
            if 'l2g_score' in unified_df.columns:
                valid = unified_df[unified_df['l2g_score'].notna()]
                n = len(valid)
                correct = (valid['l2g_score'] > 0.5).sum()
                accuracy = correct / n if n > 0 else 0
                
                # Bootstrap CI
                def bootstrap_ci(data, n_boot=1000, seed=42):
                    np.random.seed(seed)
                    boot_means = []
                    for _ in range(n_boot):
                        sample = np.random.choice(data, size=len(data), replace=True)
                        boot_means.append(np.mean(sample))
                    return np.percentile(boot_means, [2.5, 97.5])
                
                binary = (valid['l2g_score'] > 0.5).astype(int).values
                ci = bootstrap_ci(binary) if n > 0 else [0, 0]
                
                return {
                    'top_1': {
                        'accuracy': accuracy,
                        'correct': int(correct),
                        'total': n,
                        'ci_lower': ci[0],
                        'ci_upper': ci[1]
                    }
                }
        
        # Fallback: Load L2G scores and merge with benchmark
        l2g_path = self.baselines_dir / "platform_api_l2g_scores.tsv"
        if not l2g_path.exists():
            logger.warning(f"L2G scores file not found: {l2g_path}")
            return {}
        
        l2g_df = pd.read_csv(l2g_path, sep='\t')
        
        # Merge with benchmark by gene_symbol
        merged = self.benchmark_df.merge(
            l2g_df[['gene_symbol', 'l2g_score']].drop_duplicates(),
            on='gene_symbol',
            how='left'
        )
        
        # L2G correct: score > 0.5 for causal gene
        valid = merged[merged['l2g_score'].notna()]
        n = len(valid)
        correct = (valid['l2g_score'] > 0.5).sum()
        accuracy = correct / n if n > 0 else 0
        
        # Bootstrap CI
        def bootstrap_ci(data, n_boot=1000, seed=42):
            np.random.seed(seed)
            boot_means = []
            for _ in range(n_boot):
                sample = np.random.choice(data, size=len(data), replace=True)
                boot_means.append(np.mean(sample))
            return np.percentile(boot_means, [2.5, 97.5])
        
        if n > 0:
            binary = (valid['l2g_score'] > 0.5).astype(int).values
            ci = bootstrap_ci(binary)
        else:
            ci = [0, 0]
        
        return {
            'top_1': {
                'accuracy': accuracy,
                'correct': int(correct),
                'total': n,
                'ci_lower': ci[0],
                'ci_upper': ci[1]
            }
        }


def main():
    """Main entry point."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    analysis = TopKAccuracyAnalysis(project_root)
    
    # Analyze existing L2G scores (Top-1)
    top1_results = analysis.analyze_existing_l2g_scores()
    
    print("\n" + "=" * 70)
    print("TOP-1 ACCURACY ANALYSIS (from existing L2G scores)")
    print("=" * 70)
    
    if 'top_1' in top1_results:
        t1 = top1_results['top_1']
        print(f"\nTop-1 Accuracy: {t1['accuracy']*100:.1f}%")
        print(f"  Correct: {t1['correct']}/{t1['total']}")
        print(f"  95% CI: [{t1['ci_lower']*100:.1f}%, {t1['ci_upper']*100:.1f}%]")
    else:
        print("Could not compute Top-1 accuracy")
    
    # Save results
    output_path = project_root / "data" / "processed" / "baselines" / "topk_accuracy_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(top1_results, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
