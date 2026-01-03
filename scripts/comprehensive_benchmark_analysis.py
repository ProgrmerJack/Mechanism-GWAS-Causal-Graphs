#!/usr/bin/env python3
"""
Comprehensive L2G vs cS2G Statistical Comparison
=================================================

This script performs a rigorous statistical comparison between L2G and cS2G
methods on the post-2021 independent benchmark, addressing Nature Genetics
reviewer concerns about selection bias and proper statistical testing.

Key Features:
1. Full coverage analysis (L2G: 45/63, cS2G: 60/63)
2. Paired comparison on overlapping loci (n=45)
3. McNemar's test for paired proportions
4. Ensemble complementarity analysis
5. Selection bias characterization

Author: Mechanism-GWAS-Causal-Graphs
Date: 2025
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ComprehensiveBenchmarkAnalysis:
    """Comprehensive comparison of L2G vs cS2G methods."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.data_dir = project_root / "data"
        self.baselines_dir = self.data_dir / "processed" / "baselines"
        
        # Load data
        self.combined_df = None
        self.cs2g_gene_df = None
        self.l2g_scores = None
        
    def load_data(self) -> None:
        """Load all benchmark data."""
        logger.info("Loading benchmark data...")
        
        # Load combined baseline scores
        combined_path = self.baselines_dir / "benchmark_baseline_scores_combined.tsv"
        self.combined_df = pd.read_csv(combined_path, sep='\t')
        logger.info(f"Loaded combined benchmark: {len(self.combined_df)} loci")
        
        # Load gene-based cS2G scores
        cs2g_path = self.baselines_dir / "cs2g_benchmark_gene_based.tsv"
        self.cs2g_gene_df = pd.read_csv(cs2g_path, sep='\t')
        logger.info(f"Loaded gene-based cS2G: {len(self.cs2g_gene_df)} loci")
        
    def create_unified_benchmark(self) -> pd.DataFrame:
        """Create unified benchmark with both L2G and cS2G (gene-based) scores."""
        logger.info("Creating unified benchmark...")
        
        # Start with combined df
        df = self.combined_df.copy()
        
        # Update cS2G scores with gene-based results
        cs2g_lookup = self.cs2g_gene_df.set_index('locus_id')[['cs2g_score', 'cs2g_prediction', 'cs2g_best_snp']].to_dict('index')
        
        for idx, row in df.iterrows():
            locus_id = row['locus_id']
            if locus_id in cs2g_lookup:
                cs2g_data = cs2g_lookup[locus_id]
                df.at[idx, 'cs2g_score_gene_based'] = cs2g_data['cs2g_score']
                df.at[idx, 'cs2g_prediction_gene_based'] = cs2g_data['cs2g_prediction']
                df.at[idx, 'cs2g_best_snp'] = cs2g_data['cs2g_best_snp']
        
        # Classification flags
        df['l2g_available'] = df['l2g_score'].notna() & (df['l2g_score'] > 0)
        df['cs2g_available'] = df['cs2g_score_gene_based'].notna() & (df['cs2g_score_gene_based'] > 0)
        
        # Correct predictions (score > 0.5 means correct for benchmark genes)
        df['l2g_correct'] = (df['l2g_score'] > 0.5).fillna(False).astype(int)
        df['cs2g_correct'] = (df['cs2g_score_gene_based'] > 0.5).fillna(False).astype(int)
        
        return df
    
    def calculate_coverage_and_accuracy(self, df: pd.DataFrame) -> Dict:
        """Calculate coverage and accuracy for each method."""
        logger.info("Calculating coverage and accuracy...")
        
        total = len(df)
        
        # L2G
        l2g_matched = df['l2g_available'].sum()
        l2g_matched_df = df[df['l2g_available']]
        l2g_correct = (l2g_matched_df['l2g_score'] > 0.5).sum() if len(l2g_matched_df) > 0 else 0
        l2g_accuracy = l2g_correct / len(l2g_matched_df) if len(l2g_matched_df) > 0 else 0
        
        # cS2G (gene-based)
        cs2g_matched = df['cs2g_available'].sum()
        cs2g_matched_df = df[df['cs2g_available']]
        cs2g_correct = (cs2g_matched_df['cs2g_score_gene_based'] > 0.5).sum() if len(cs2g_matched_df) > 0 else 0
        cs2g_accuracy = cs2g_correct / len(cs2g_matched_df) if len(cs2g_matched_df) > 0 else 0
        
        # Bootstrap CIs
        def bootstrap_ci(data, stat_fn, n_bootstrap=1000, seed=42):
            np.random.seed(seed)
            boot_stats = []
            for _ in range(n_bootstrap):
                sample = np.random.choice(data, size=len(data), replace=True)
                boot_stats.append(stat_fn(sample))
            return np.percentile(boot_stats, [2.5, 97.5])
        
        # L2G bootstrap
        l2g_binary = (l2g_matched_df['l2g_score'] > 0.5).astype(int).values
        if len(l2g_binary) > 0:
            l2g_ci = bootstrap_ci(l2g_binary, np.mean)
        else:
            l2g_ci = [0, 0]
        
        # cS2G bootstrap
        cs2g_binary = (cs2g_matched_df['cs2g_score_gene_based'] > 0.5).astype(int).values
        if len(cs2g_binary) > 0:
            cs2g_ci = bootstrap_ci(cs2g_binary, np.mean)
        else:
            cs2g_ci = [0, 0]
        
        results = {
            'total_loci': total,
            'l2g': {
                'matched': int(l2g_matched),
                'coverage': l2g_matched / total,
                'correct': int(l2g_correct),
                'accuracy': l2g_accuracy,
                'ci_lower': l2g_ci[0],
                'ci_upper': l2g_ci[1]
            },
            'cs2g_gene_based': {
                'matched': int(cs2g_matched),
                'coverage': cs2g_matched / total,
                'correct': int(cs2g_correct),
                'accuracy': cs2g_accuracy,
                'ci_lower': cs2g_ci[0],
                'ci_upper': cs2g_ci[1]
            }
        }
        
        return results
    
    def mcnemar_test(self, df: pd.DataFrame) -> Dict:
        """
        Perform McNemar's test for paired proportions.
        
        This tests whether L2G and cS2G have significantly different accuracy
        on the same set of loci.
        """
        logger.info("Performing McNemar's test...")
        
        # Get paired loci (both methods have scores)
        paired_df = df[df['l2g_available'] & df['cs2g_available']].copy()
        n_paired = len(paired_df)
        logger.info(f"Paired loci: {n_paired}")
        
        if n_paired == 0:
            return {'error': 'No paired loci available'}
        
        # Create contingency table for McNemar's test
        # Rows: L2G correct/incorrect
        # Cols: cS2G correct/incorrect
        l2g_correct = (paired_df['l2g_score'] > 0.5).astype(int).values
        cs2g_correct = (paired_df['cs2g_score_gene_based'] > 0.5).astype(int).values
        
        # Contingency table cells
        # b: L2G correct, cS2G incorrect
        # c: L2G incorrect, cS2G correct
        b = np.sum((l2g_correct == 1) & (cs2g_correct == 0))
        c = np.sum((l2g_correct == 0) & (cs2g_correct == 1))
        a = np.sum((l2g_correct == 1) & (cs2g_correct == 1))  # Both correct
        d = np.sum((l2g_correct == 0) & (cs2g_correct == 0))  # Both incorrect
        
        # McNemar's test statistic (chi-squared with continuity correction)
        if b + c > 0:
            chi2 = ((abs(b - c) - 1) ** 2) / (b + c)
            p_value = 1 - scipy_stats.chi2.cdf(chi2, df=1)
        else:
            chi2 = 0
            p_value = 1.0
        
        # Exact McNemar's test (binomial)
        if b + c > 0:
            exact_p = 2 * min(
                scipy_stats.binom.cdf(min(b, c), b + c, 0.5),
                1 - scipy_stats.binom.cdf(max(b, c) - 1, b + c, 0.5)
            )
        else:
            exact_p = 1.0
        
        # Paired accuracies
        l2g_acc_paired = l2g_correct.mean()
        cs2g_acc_paired = cs2g_correct.mean()
        
        results = {
            'n_paired': n_paired,
            'contingency_table': {
                'both_correct': int(a),
                'l2g_only_correct': int(b),
                'cs2g_only_correct': int(c),
                'both_incorrect': int(d)
            },
            'l2g_accuracy_paired': l2g_acc_paired,
            'cs2g_accuracy_paired': cs2g_acc_paired,
            'accuracy_difference': cs2g_acc_paired - l2g_acc_paired,
            'mcnemar_chi2': chi2,
            'mcnemar_p_value': p_value,
            'mcnemar_exact_p': exact_p,
            'discordant_pairs': int(b + c)
        }
        
        return results
    
    def ensemble_complementarity(self, df: pd.DataFrame) -> Dict:
        """
        Analyze ensemble complementarity - what happens if we combine methods?
        
        This addresses the reviewer concern about whether different methods
        capture different biology.
        """
        logger.info("Analyzing ensemble complementarity...")
        
        # All loci where at least one method is available
        any_available = df[df['l2g_available'] | df['cs2g_available']].copy()
        
        # Union accuracy: correct if EITHER method is correct
        union_correct = ((any_available['l2g_score'] > 0.5) | 
                         (any_available['cs2g_score_gene_based'] > 0.5)).sum()
        union_acc = union_correct / len(any_available)
        
        # Intersection accuracy: correct if BOTH methods are correct
        paired = any_available[any_available['l2g_available'] & any_available['cs2g_available']]
        if len(paired) > 0:
            intersection_correct = ((paired['l2g_score'] > 0.5) & 
                                    (paired['cs2g_score_gene_based'] > 0.5)).sum()
            intersection_acc = intersection_correct / len(paired)
        else:
            intersection_correct = 0
            intersection_acc = 0
        
        # Unique contributions
        l2g_only = df['l2g_available'] & ~df['cs2g_available']
        cs2g_only = df['cs2g_available'] & ~df['l2g_available']
        
        l2g_unique_correct = df[l2g_only]['l2g_correct'].sum() if l2g_only.any() else 0
        cs2g_unique_correct = df[cs2g_only]['cs2g_correct'].sum() if cs2g_only.any() else 0
        
        results = {
            'loci_any_method': len(any_available),
            'loci_both_methods': len(paired),
            'union_accuracy': union_acc,
            'union_correct': int(union_correct),
            'intersection_accuracy': intersection_acc,
            'intersection_correct': int(intersection_correct),
            'l2g_unique_loci': int(l2g_only.sum()),
            'l2g_unique_correct': int(l2g_unique_correct),
            'cs2g_unique_loci': int(cs2g_only.sum()),
            'cs2g_unique_correct': int(cs2g_unique_correct),
            'complementarity_gain': union_acc - max(
                df[df['l2g_available']]['l2g_correct'].mean() if df['l2g_available'].any() else 0,
                df[df['cs2g_available']]['cs2g_correct'].mean() if df['cs2g_available'].any() else 0
            )
        }
        
        return results
    
    def selection_bias_analysis(self, df: pd.DataFrame) -> Dict:
        """
        Characterize potential selection bias between matched and unmatched loci.
        
        This addresses reviewer concern about whether matched loci are representative.
        """
        logger.info("Analyzing selection bias...")
        
        # Characterize L2G matched vs unmatched
        l2g_matched = df[df['l2g_available']]
        l2g_unmatched = df[~df['l2g_available']]
        
        # cS2G matched vs unmatched
        cs2g_matched = df[df['cs2g_available']]
        cs2g_unmatched = df[~df['cs2g_available']]
        
        # Evidence tier distribution
        def tier_dist(subset):
            if 'evidence_tier' in subset.columns:
                return subset['evidence_tier'].value_counts().to_dict()
            return {}
        
        # Chromosome distribution
        def chr_dist(subset):
            return subset['chr'].value_counts().to_dict()
        
        results = {
            'l2g_matched': {
                'n': len(l2g_matched),
                'evidence_tiers': tier_dist(l2g_matched),
                'chromosomes': chr_dist(l2g_matched),
                'has_x_chr': (l2g_matched['chr'] == 'X').any()
            },
            'l2g_unmatched': {
                'n': len(l2g_unmatched),
                'evidence_tiers': tier_dist(l2g_unmatched),
                'chromosomes': chr_dist(l2g_unmatched),
                'genes': l2g_unmatched['gene_symbol'].tolist()
            },
            'cs2g_matched': {
                'n': len(cs2g_matched),
                'evidence_tiers': tier_dist(cs2g_matched),
                'chromosomes': chr_dist(cs2g_matched),
                'has_x_chr': (cs2g_matched['chr'] == 'X').any()
            },
            'cs2g_unmatched': {
                'n': len(cs2g_unmatched),
                'evidence_tiers': tier_dist(cs2g_unmatched),
                'chromosomes': chr_dist(cs2g_unmatched),
                'genes': cs2g_unmatched['gene_symbol'].tolist(),
                'reason': 'X chromosome (cS2G only covers autosomes 1-22)'
            }
        }
        
        return results
    
    def run_full_analysis(self) -> Dict:
        """Run complete analysis pipeline."""
        logger.info("=" * 70)
        logger.info("COMPREHENSIVE L2G vs cS2G BENCHMARK ANALYSIS")
        logger.info("=" * 70)
        
        # Load data
        self.load_data()
        
        # Create unified benchmark
        unified_df = self.create_unified_benchmark()
        
        # Run all analyses
        coverage_accuracy = self.calculate_coverage_and_accuracy(unified_df)
        mcnemar_results = self.mcnemar_test(unified_df)
        ensemble_results = self.ensemble_complementarity(unified_df)
        bias_results = self.selection_bias_analysis(unified_df)
        
        # Compile full results
        results = {
            'coverage_and_accuracy': coverage_accuracy,
            'mcnemar_test': mcnemar_results,
            'ensemble_complementarity': ensemble_results,
            'selection_bias': bias_results
        }
        
        # Print summary
        self._print_summary(results)
        
        # Save results
        self._save_results(results, unified_df)
        
        return results
    
    def _print_summary(self, results: Dict) -> None:
        """Print formatted summary."""
        ca = results['coverage_and_accuracy']
        mcn = results['mcnemar_test']
        ens = results['ensemble_complementarity']
        
        print("\n" + "=" * 70)
        print("SUMMARY FOR MANUSCRIPT")
        print("=" * 70)
        
        print("\n--- Coverage and Accuracy ---")
        print(f"Total benchmark loci: {ca['total_loci']}")
        print(f"\nL2G (Platform API):")
        print(f"  Coverage: {ca['l2g']['matched']}/{ca['total_loci']} ({ca['l2g']['coverage']*100:.1f}%)")
        print(f"  Accuracy: {ca['l2g']['accuracy']*100:.1f}% (95% CI: {ca['l2g']['ci_lower']*100:.1f}-{ca['l2g']['ci_upper']*100:.1f}%)")
        
        print(f"\ncS2G (Gene-Based):")
        print(f"  Coverage: {ca['cs2g_gene_based']['matched']}/{ca['total_loci']} ({ca['cs2g_gene_based']['coverage']*100:.1f}%)")
        print(f"  Accuracy: {ca['cs2g_gene_based']['accuracy']*100:.1f}% (95% CI: {ca['cs2g_gene_based']['ci_lower']*100:.1f}-{ca['cs2g_gene_based']['ci_upper']*100:.1f}%)")
        
        print("\n--- McNemar's Test (Paired Comparison) ---")
        print(f"Paired loci: {mcn['n_paired']}")
        print(f"L2G accuracy (paired): {mcn['l2g_accuracy_paired']*100:.1f}%")
        print(f"cS2G accuracy (paired): {mcn['cs2g_accuracy_paired']*100:.1f}%")
        print(f"Accuracy difference: {mcn['accuracy_difference']*100:.1f} percentage points")
        print(f"McNemar's chi² = {mcn['mcnemar_chi2']:.3f}, P = {mcn['mcnemar_p_value']:.4f}")
        print(f"Exact McNemar P = {mcn['mcnemar_exact_p']:.4f}")
        print(f"Contingency: both correct={mcn['contingency_table']['both_correct']}, " +
              f"L2G only={mcn['contingency_table']['l2g_only_correct']}, " +
              f"cS2G only={mcn['contingency_table']['cs2g_only_correct']}, " +
              f"both wrong={mcn['contingency_table']['both_incorrect']}")
        
        print("\n--- Ensemble Complementarity ---")
        print(f"Union accuracy (either correct): {ens['union_accuracy']*100:.1f}%")
        print(f"Intersection accuracy (both correct): {ens['intersection_accuracy']*100:.1f}%")
        print(f"Complementarity gain: {ens['complementarity_gain']*100:.1f} percentage points")
        
        print("\n--- Selection Bias ---")
        bias = results['selection_bias']
        print(f"L2G unmatched genes: {bias['l2g_unmatched']['genes']}")
        print(f"cS2G unmatched genes: {bias['cs2g_unmatched']['genes']}")
        print(f"cS2G limitation: {bias['cs2g_unmatched']['reason']}")
        
        print("\n" + "=" * 70)
    
    def _save_results(self, results: Dict, unified_df: pd.DataFrame) -> None:
        """Save results to files."""
        # Save JSON summary
        json_path = self.baselines_dir / "comprehensive_benchmark_analysis.json"
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        logger.info(f"Saved analysis to {json_path}")
        
        # Save unified benchmark TSV
        tsv_path = self.baselines_dir / "unified_benchmark_l2g_cs2g.tsv"
        unified_df.to_csv(tsv_path, sep='\t', index=False)
        logger.info(f"Saved unified benchmark to {tsv_path}")


def main():
    """Main entry point."""
    # Find project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    # Run analysis
    analysis = ComprehensiveBenchmarkAnalysis(project_root)
    results = analysis.run_full_analysis()


if __name__ == "__main__":
    main()
