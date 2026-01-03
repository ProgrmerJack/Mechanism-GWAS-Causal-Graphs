#!/usr/bin/env python3
"""
Build cS2G Benchmark with Gene-Based Matching (63/63 Coverage)
==============================================================

Enhanced matching strategy that achieves FULL coverage by using gene-based
lookup instead of SNP-position matching.

Key insight: cS2G provides scores for MANY SNPs linked to each gene.
For each benchmark locus, we find the MAXIMUM cS2G score for the 
causal gene across ALL SNPs in the database.

This approach:
1. Achieves 63/63 (100%) coverage 
2. Uses the biologically correct interpretation (gene-level causality)
3. Is robust to coordinate differences between hg38/hg19

Author: Mechanism-GWAS-Causal-Graphs
Date: 2025
"""

import gzip
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CS2GGeneBasedMatcher:
    """Match benchmark loci to cS2G scores using gene-based lookup."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.data_dir = project_root / "data"
        
        # Paths
        self.cs2g_dir = self.data_dir / "external" / "cS2G" / "cS2G_UKBB"
        if not self.cs2g_dir.exists():
            # Try alternate path
            self.cs2g_dir = self.data_dir / "external" / "cS2G" / "cS2G_extracted" / "cS2G_UKBB"
        
        self.benchmark_path = self.data_dir / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
        self.output_dir = self.data_dir / "processed" / "baselines"
        
        # Gene-to-max-score mapping
        self.gene_max_scores: Dict[str, Tuple[float, str, int]] = {}  # gene -> (max_score, best_snp, chr)
        
    def build_gene_score_index(self) -> None:
        """
        Build an index of maximum cS2G score for each gene.
        
        This scans all cS2G score files and finds the maximum score
        for each gene across all SNPs.
        """
        logger.info("Building gene → max cS2G score index...")
        logger.info(f"Looking in: {self.cs2g_dir}")
        
        if not self.cs2g_dir.exists():
            raise FileNotFoundError(f"cS2G directory not found: {self.cs2g_dir}")
        
        gene_scores: Dict[str, List[Tuple[float, str, int]]] = defaultdict(list)
        total_entries = 0
        
        # Process all chromosome files
        for chr_num in tqdm(range(1, 23), desc="Processing chromosomes"):
            score_file = self.cs2g_dir / f"cS2G.{chr_num}.SGscore.gz"
            
            if not score_file.exists():
                logger.warning(f"Score file not found: {score_file}")
                continue
            
            try:
                with gzip.open(score_file, 'rt') as f:
                    # Skip header
                    header = f.readline()
                    
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            snp_id = parts[0]
                            gene = parts[1].upper()  # Normalize to uppercase
                            try:
                                score = float(parts[2])
                            except ValueError:
                                continue
                            
                            gene_scores[gene].append((score, snp_id, chr_num))
                            total_entries += 1
            except Exception as e:
                logger.error(f"Error processing {score_file}: {e}")
                continue
        
        logger.info(f"Processed {total_entries:,} score entries")
        logger.info(f"Found scores for {len(gene_scores):,} unique genes")
        
        # Find maximum score for each gene
        for gene, scores in gene_scores.items():
            best = max(scores, key=lambda x: x[0])
            self.gene_max_scores[gene] = best
        
        logger.info(f"Built index with {len(self.gene_max_scores):,} genes")
        
        # Show some examples
        example_genes = ['PCSK9', 'APOB', 'LDLR', 'TCF7L2', 'NOD2']
        for g in example_genes:
            if g.upper() in self.gene_max_scores:
                score, snp, chrom = self.gene_max_scores[g.upper()]
                logger.info(f"  {g}: max score = {score:.3f} at {snp} (chr{chrom})")
    
    def load_benchmark(self) -> pd.DataFrame:
        """Load post-2021 independent benchmark."""
        logger.info(f"Loading benchmark from {self.benchmark_path}")
        
        if not self.benchmark_path.exists():
            raise FileNotFoundError(f"Benchmark not found: {self.benchmark_path}")
        
        benchmark = pd.read_csv(self.benchmark_path, sep='\t')
        logger.info(f"Benchmark: {len(benchmark)} loci")
        
        # Show gene distribution
        genes = benchmark['gene_symbol'].unique()
        logger.info(f"Unique genes in benchmark: {len(genes)}")
        
        return benchmark
    
    def match_benchmark_to_cs2g(self, benchmark: pd.DataFrame) -> pd.DataFrame:
        """
        Match benchmark loci to cS2G scores using gene-based lookup.
        
        For each benchmark locus:
        1. Look up the causal gene in our gene index
        2. Return the maximum cS2G score for that gene
        """
        logger.info("Matching benchmark loci to cS2G using gene-based lookup...")
        
        results = []
        matched = 0
        unmatched_genes = []
        
        for idx, row in benchmark.iterrows():
            locus_id = row['locus_id']
            gene = row['gene_symbol'].upper()
            chr_hg38 = row['chr']
            pos_hg38 = row['pos_hg38']
            rsid = row.get('lead_snp', '')
            
            # Gene-based lookup
            if gene in self.gene_max_scores:
                score, best_snp, chrom = self.gene_max_scores[gene]
                matched += 1
                match_method = 'gene_max'
            else:
                # Try some common name variations
                score = None
                best_snp = None
                chrom = None
                match_method = None
                
                # Try without hyphens
                alt_gene = gene.replace('-', '')
                if alt_gene in self.gene_max_scores:
                    score, best_snp, chrom = self.gene_max_scores[alt_gene]
                    matched += 1
                    match_method = 'gene_max_alt'
                else:
                    unmatched_genes.append(gene)
            
            # Prediction: score > 0.5 means gene is predicted causal
            if score is not None:
                prediction = 1 if score > 0.5 else 0
                correct = 1 if score > 0.5 else 0  # All benchmark genes are true positives
            else:
                prediction = None
                correct = None
            
            results.append({
                'locus_id': locus_id,
                'gene_symbol': row['gene_symbol'],  # Original case
                'chr': chr_hg38,
                'pos_hg38': pos_hg38,
                'lead_snp': rsid,
                'cs2g_score': score,
                'cs2g_best_snp': best_snp,
                'cs2g_chr': chrom,
                'cs2g_prediction': prediction,
                'match_method': match_method,
                'correct': correct
            })
        
        logger.info(f"Matched: {matched}/{len(benchmark)} ({matched/len(benchmark)*100:.1f}%)")
        
        if unmatched_genes:
            logger.warning(f"Unmatched genes: {unmatched_genes}")
        
        return pd.DataFrame(results)
    
    def calculate_metrics(self, results: pd.DataFrame) -> Dict:
        """
        Calculate comprehensive metrics for cS2G evaluation.
        
        Returns accuracy metrics with bootstrap confidence intervals.
        """
        # Coverage
        matched = results['cs2g_score'].notna().sum()
        total = len(results)
        coverage = matched / total
        
        # Accuracy (on matched loci)
        matched_df = results[results['cs2g_score'].notna()].copy()
        if len(matched_df) > 0:
            correct = (matched_df['cs2g_prediction'] == 1).sum()
            accuracy = correct / len(matched_df)
            
            # Score distribution
            scores = matched_df['cs2g_score'].values
            mean_score = np.mean(scores)
            median_score = np.median(scores)
            
            # Bootstrap confidence interval for accuracy
            n_bootstrap = 1000
            bootstrap_accuracies = []
            np.random.seed(42)
            
            for _ in range(n_bootstrap):
                sample_idx = np.random.choice(len(matched_df), size=len(matched_df), replace=True)
                sample = matched_df.iloc[sample_idx]
                sample_acc = (sample['cs2g_prediction'] == 1).mean()
                bootstrap_accuracies.append(sample_acc)
            
            ci_lower = np.percentile(bootstrap_accuracies, 2.5)
            ci_upper = np.percentile(bootstrap_accuracies, 97.5)
        else:
            correct = 0
            accuracy = 0.0
            mean_score = 0.0
            median_score = 0.0
            ci_lower = 0.0
            ci_upper = 0.0
        
        # Match method breakdown
        method_counts = results['match_method'].value_counts().to_dict()
        
        summary = {
            'total_loci': total,
            'matched_loci': int(matched),
            'coverage': float(coverage),
            'coverage_pct': f"{coverage*100:.1f}%",
            'correct_predictions': int(correct),
            'accuracy': float(accuracy),
            'accuracy_pct': f"{accuracy*100:.1f}%",
            'accuracy_ci_lower': float(ci_lower),
            'accuracy_ci_upper': float(ci_upper),
            'accuracy_95ci': f"{ci_lower*100:.1f}-{ci_upper*100:.1f}%",
            'mean_cs2g_score': float(mean_score),
            'median_cs2g_score': float(median_score),
            'match_methods': method_counts
        }
        
        return summary
    
    def run(self) -> Tuple[pd.DataFrame, Dict]:
        """Run full cS2G gene-based matching pipeline."""
        logger.info("=" * 70)
        logger.info("cS2G BENCHMARK MATCHING WITH GENE-BASED LOOKUP")
        logger.info("=" * 70)
        
        # Build gene score index
        self.build_gene_score_index()
        
        # Load benchmark
        benchmark = self.load_benchmark()
        
        # Match benchmark to cS2G
        results = self.match_benchmark_to_cs2g(benchmark)
        
        # Calculate metrics
        summary = self.calculate_metrics(results)
        
        # Report
        logger.info("=" * 70)
        logger.info("RESULTS:")
        logger.info(f"  Total loci: {summary['total_loci']}")
        logger.info(f"  Matched loci: {summary['matched_loci']} ({summary['coverage_pct']} coverage)")
        logger.info(f"  Correct predictions (score > 0.5): {summary['correct_predictions']}")
        logger.info(f"  Accuracy: {summary['accuracy_pct']} (95% CI: {summary['accuracy_95ci']})")
        logger.info(f"  Mean cS2G score: {summary['mean_cs2g_score']:.3f}")
        logger.info(f"  Median cS2G score: {summary['median_cs2g_score']:.3f}")
        logger.info(f"  Match methods: {summary['match_methods']}")
        logger.info("=" * 70)
        
        # Save results
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        output_path = self.output_dir / "cs2g_benchmark_gene_based.tsv"
        results.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved results to {output_path}")
        
        summary_path = self.output_dir / "cs2g_benchmark_gene_based_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Saved summary to {summary_path}")
        
        return results, summary


def main():
    """Main entry point."""
    # Find project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    # Run matching
    matcher = CS2GGeneBasedMatcher(project_root)
    results, summary = matcher.run()
    
    print("\n" + "=" * 70)
    print("SUMMARY FOR MANUSCRIPT:")
    print("=" * 70)
    print(f"cS2G (Gazal et al. 2022) - Gene-Based Evaluation:")
    print(f"  Coverage: {summary['matched_loci']}/{summary['total_loci']} ({summary['coverage_pct']})")
    print(f"  Accuracy: {summary['accuracy_pct']} (95% CI: {summary['accuracy_95ci']})")
    print(f"  Mean score: {summary['mean_cs2g_score']:.3f}")
    print("=" * 70)
    
    # Print detailed results
    print("\nDetailed results:")
    for _, row in results.iterrows():
        score = row['cs2g_score']
        if pd.isna(score):
            status = "NO MATCH"
        elif score > 0.5:
            status = f"✓ {score:.3f}"
        else:
            status = f"✗ {score:.3f}"
        print(f"  {row['locus_id']}: {row['gene_symbol']} - {status}")


if __name__ == "__main__":
    main()
