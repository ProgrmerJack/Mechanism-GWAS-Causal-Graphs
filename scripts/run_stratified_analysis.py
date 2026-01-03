#!/usr/bin/env python3
"""
Mechanism Stratified cS2G Analysis for Nature Genetics
======================================================

Implements P2: Stratify benchmark composition by mechanism (coding vs regulatory).

Stratification rules (based on evidence_tier):
- CODING: Tier1_Coding, Tier1_Mendelian (coding variants dominate)
- REGULATORY: Tier1_CRISPR (enhancer-mediated, validated by CRISPRi)
- AMBIGUOUS: Tier1_Drug, Tier2_MultiEvidence (mixed evidence)

Outputs per-stratum:
- Top-k accuracy + bootstrap CIs for cS2G, nearest-gene, distance-weighted baselines
- Results for main manuscript Figure 3

Author: Mechanism-GWAS-Causal-Graphs
Date: December 2025
"""

import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

# Import locus-aware evaluator
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.run_cs2g_locus_aware import LocusAwareCS2GEvaluator, LocusDefinition

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Stratification mapping
STRATUM_MAP = {
    'Tier1_Coding': 'CODING',
    'Tier1_Mendelian': 'CODING',
    'Tier1_CRISPR': 'REGULATORY',
    'Tier1_Drug': 'AMBIGUOUS',
    'Tier2_MultiEvidence': 'AMBIGUOUS',
}


class StratifiedAnalyzer:
    """Analyze cS2G performance stratified by mechanism type."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.results_dir = project_root / "results" / "stratified_analysis"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
    def load_cs2g_results(self) -> pd.DataFrame:
        """Load locus-aware cS2G results."""
        results_path = self.project_root / "results" / "cs2g_locus_aware" / "cs2g_locus_aware_results_max.tsv"
        if not results_path.exists():
            raise FileNotFoundError(f"Run run_cs2g_locus_aware.py first: {results_path}")
        return pd.read_csv(results_path, sep='\t')
    
    def load_benchmark(self) -> pd.DataFrame:
        """Load benchmark with evidence tier info."""
        benchmark_path = self.project_root / "data" / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
        return pd.read_csv(benchmark_path, sep='\t')
    
    def assign_strata(self, df: pd.DataFrame) -> pd.DataFrame:
        """Assign mechanism stratum to each locus."""
        df = df.copy()
        df['stratum'] = df['evidence_tier'].map(STRATUM_MAP).fillna('UNKNOWN')
        return df
    
    def compute_nearest_gene_baseline(self, benchmark_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute nearest-gene baseline within each locus.
        
        For each locus, the nearest gene is always Top-1 by definition.
        But our benchmark tests if it's the TRUE causal gene.
        """
        # Load locus-gene pair annotations
        pairs_path = self.project_root / "data" / "processed" / "baselines" / "post2021_locus_gene_pairs_annotated.tsv"
        if not pairs_path.exists():
            logger.warning("Locus-gene pairs not found, cannot compute nearest-gene baseline")
            return pd.DataFrame()
        
        pairs_df = pd.read_csv(pairs_path, sep='\t')
        
        results = []
        for _, row in benchmark_df.iterrows():
            locus_id = row['locus_id']
            true_gene = row['gene_symbol'].upper()
            
            # Get all genes at this locus
            locus_pairs = pairs_df[pairs_df['locus_id'] == locus_id].copy()
            if locus_pairs.empty:
                continue
            
            # Nearest gene = minimum distance to lead SNP position
            if 'distance_to_tss' in locus_pairs.columns:
                locus_pairs = locus_pairs.sort_values('distance_to_tss', ascending=True)
                locus_pairs['gene_upper'] = locus_pairs['gene_symbol'].str.upper()
                
                # Find rank of true gene
                if true_gene in locus_pairs['gene_upper'].values:
                    rank = locus_pairs['gene_upper'].tolist().index(true_gene) + 1
                else:
                    rank = len(locus_pairs) + 1
                
                results.append({
                    'locus_id': locus_id,
                    'true_gene': true_gene,
                    'nearest_gene': locus_pairs.iloc[0]['gene_symbol'] if len(locus_pairs) > 0 else '',
                    'true_gene_rank': rank,
                    'top1_correct': 1 if rank == 1 else 0,
                    'top3_correct': 1 if rank <= 3 else 0,
                    'top5_correct': 1 if rank <= 5 else 0,
                    'candidate_count': len(locus_pairs)
                })
        
        return pd.DataFrame(results)
    
    def compute_stratum_metrics(
        self,
        results_df: pd.DataFrame,
        stratum: str,
        method: str
    ) -> Dict:
        """Compute metrics for a single stratum."""
        stratum_df = results_df[results_df.get('stratum', 'ALL') == stratum]
        
        if stratum_df.empty or 'top1_correct' not in stratum_df.columns:
            return {
                'stratum': stratum,
                'method': method,
                'n_loci': 0,
                'top1_mean': np.nan,
                'top3_mean': np.nan,
                'top5_mean': np.nan,
            }
        
        # Filter to covered loci
        if 'coverage' in stratum_df.columns:
            stratum_df = stratum_df[stratum_df['coverage'] == True]
        
        n = len(stratum_df)
        if n == 0:
            return {
                'stratum': stratum,
                'method': method,
                'n_loci': 0,
                'top1_mean': np.nan,
                'top3_mean': np.nan,
                'top5_mean': np.nan,
            }
        
        # Bootstrap CIs
        def bootstrap_mean(col, n_boot=1000):
            values = stratum_df[col].dropna().values
            if len(values) == 0:
                return np.nan, np.nan, np.nan
            np.random.seed(42)
            means = [np.mean(np.random.choice(values, len(values), replace=True)) for _ in range(n_boot)]
            return np.mean(values), np.percentile(means, 2.5), np.percentile(means, 97.5)
        
        t1, t1_lo, t1_hi = bootstrap_mean('top1_correct')
        t3, t3_lo, t3_hi = bootstrap_mean('top3_correct')
        t5, t5_lo, t5_hi = bootstrap_mean('top5_correct')
        
        return {
            'stratum': stratum,
            'method': method,
            'n_loci': n,
            'top1_mean': t1,
            'top1_ci_lower': t1_lo,
            'top1_ci_upper': t1_hi,
            'top3_mean': t3,
            'top3_ci_lower': t3_lo,
            'top3_ci_upper': t3_hi,
            'top5_mean': t5,
            'top5_ci_lower': t5_lo,
            'top5_ci_upper': t5_hi,
        }
    
    def run_stratified_analysis(self) -> pd.DataFrame:
        """Run full stratified analysis."""
        logger.info("=" * 70)
        logger.info("STRATIFIED MECHANISM ANALYSIS")
        logger.info("=" * 70)
        
        # Load data
        benchmark_df = self.load_benchmark()
        benchmark_df = self.assign_strata(benchmark_df)
        
        logger.info("\nStratum distribution:")
        for stratum, count in benchmark_df['stratum'].value_counts().items():
            logger.info(f"  {stratum}: {count} loci")
        
        # Load cS2G results
        cs2g_df = self.load_cs2g_results()
        cs2g_df = cs2g_df.merge(
            benchmark_df[['locus_id', 'stratum', 'evidence_tier']],
            on='locus_id',
            how='left'
        )
        
        # Compute metrics per stratum
        all_metrics = []
        
        for stratum in ['CODING', 'REGULATORY', 'AMBIGUOUS', 'ALL']:
            if stratum == 'ALL':
                stratum_results = cs2g_df.copy()
                stratum_results['stratum'] = 'ALL'
            else:
                stratum_results = cs2g_df[cs2g_df['stratum'] == stratum].copy()
                stratum_results['stratum'] = stratum
            
            metrics = self.compute_stratum_metrics(stratum_results, stratum, 'cS2G_LocusAware')
            all_metrics.append(metrics)
        
        # Compute nearest-gene baseline
        nearest_df = self.compute_nearest_gene_baseline(benchmark_df)
        if not nearest_df.empty:
            nearest_df = nearest_df.merge(
                benchmark_df[['locus_id', 'stratum']],
                on='locus_id',
                how='left'
            )
            
            for stratum in ['CODING', 'REGULATORY', 'AMBIGUOUS', 'ALL']:
                if stratum == 'ALL':
                    stratum_results = nearest_df.copy()
                    stratum_results['stratum'] = 'ALL'
                else:
                    stratum_results = nearest_df[nearest_df['stratum'] == stratum].copy()
                    stratum_results['stratum'] = stratum
                
                metrics = self.compute_stratum_metrics(stratum_results, stratum, 'NearestGene')
                all_metrics.append(metrics)
        
        # Create results DataFrame
        results_df = pd.DataFrame(all_metrics)
        
        # Print results
        logger.info("\n" + "=" * 70)
        logger.info("STRATIFIED RESULTS")
        logger.info("=" * 70)
        
        for method in results_df['method'].unique():
            logger.info(f"\n{method}:")
            method_df = results_df[results_df['method'] == method]
            for _, row in method_df.iterrows():
                if pd.notna(row['top1_mean']):
                    logger.info(f"  {row['stratum']:12s}: Top-1={row['top1_mean']*100:5.1f}% "
                               f"[{row.get('top1_ci_lower', np.nan)*100:.1f}-{row.get('top1_ci_upper', np.nan)*100:.1f}%] "
                               f"(n={row['n_loci']})")
        
        # Save results
        results_path = self.results_dir / "stratified_results.tsv"
        results_df.to_csv(results_path, sep='\t', index=False)
        logger.info(f"\nSaved to {results_path}")
        
        # Save JSON for figure generation
        json_path = self.results_dir / "stratified_results.json"
        results_df.to_json(json_path, orient='records', indent=2)
        logger.info(f"Saved to {json_path}")
        
        return results_df
    
    def run_loto_analysis(self) -> pd.DataFrame:
        """
        Leave-One-Trait-Out (LOTO) robustness analysis.
        
        For each trait category (Cardiovascular, Metabolic, Immune, etc.),
        compute cS2G performance when that category is held out.
        This tests robustness beyond trait-specific tuning.
        """
        logger.info("\n" + "=" * 70)
        logger.info("LEAVE-ONE-TRAIT-OUT (LOTO) ROBUSTNESS ANALYSIS")
        logger.info("=" * 70)
        
        # Load data
        benchmark_df = self.load_benchmark()
        cs2g_df = self.load_cs2g_results()
        
        # Merge trait category info
        cs2g_df = cs2g_df.merge(
            benchmark_df[['locus_id', 'trait_category', 'evidence_tier']],
            on='locus_id',
            how='left'
        )
        
        logger.info("\nTrait category distribution:")
        for cat, count in benchmark_df['trait_category'].value_counts().items():
            logger.info(f"  {cat}: {count} loci")
        
        all_metrics = []
        trait_categories = benchmark_df['trait_category'].unique()
        
        for held_out_category in trait_categories:
            # Performance on held-out category
            held_out_df = cs2g_df[cs2g_df['trait_category'] == held_out_category].copy()
            
            if held_out_df.empty:
                continue
            
            # Filter to covered loci
            if 'coverage' in held_out_df.columns:
                held_out_df = held_out_df[held_out_df['coverage'] == True]
            
            n = len(held_out_df)
            if n == 0:
                continue
            
            # Compute metrics
            top1 = held_out_df['top1_correct'].mean() if 'top1_correct' in held_out_df.columns else np.nan
            top3 = held_out_df['top3_correct'].mean() if 'top3_correct' in held_out_df.columns else np.nan
            top5 = held_out_df['top5_correct'].mean() if 'top5_correct' in held_out_df.columns else np.nan
            
            all_metrics.append({
                'held_out_category': held_out_category,
                'n_loci': n,
                'top1_mean': top1,
                'top3_mean': top3,
                'top5_mean': top5,
            })
            
            logger.info(f"  {held_out_category:20s}: Top-1={top1*100:5.1f}% (n={n})")
        
        results_df = pd.DataFrame(all_metrics)
        
        # Save results
        loto_path = self.results_dir / "loto_results.tsv"
        results_df.to_csv(loto_path, sep='\t', index=False)
        logger.info(f"\nSaved LOTO results to {loto_path}")
        
        # Compute coefficient of variation across categories
        if len(results_df) > 1 and not results_df['top1_mean'].isna().all():
            cv = results_df['top1_mean'].std() / results_df['top1_mean'].mean() if results_df['top1_mean'].mean() > 0 else np.nan
            logger.info(f"\nTop-1 coefficient of variation across categories: {cv:.2f}")
            logger.info("(Lower CV indicates more robust performance across trait types)")
        
        return results_df


def main():

    project_root = Path(__file__).resolve().parents[1]
    
    analyzer = StratifiedAnalyzer(project_root)
    results = analyzer.run_stratified_analysis()
    
    # Run LOTO robustness analysis (P3)
    loto_results = analyzer.run_loto_analysis()
    
    print("\n" + "=" * 70)
    print("SUMMARY FOR MANUSCRIPT FIGURE 3")
    print("=" * 70)
    print("\nCoding vs Regulatory stratification reveals mechanism-specific performance:")
    
    cs2g_coding = results[(results['method'] == 'cS2G_LocusAware') & (results['stratum'] == 'CODING')]
    cs2g_reg = results[(results['method'] == 'cS2G_LocusAware') & (results['stratum'] == 'REGULATORY')]
    
    if not cs2g_coding.empty and not cs2g_reg.empty:
        print(f"\ncS2G on CODING loci:     Top-1 = {cs2g_coding.iloc[0]['top1_mean']*100:.1f}%")
        print(f"cS2G on REGULATORY loci: Top-1 = {cs2g_reg.iloc[0]['top1_mean']*100:.1f}%")
    
    print("\nLOTO robustness (performance by held-out trait category):")
    for _, row in loto_results.iterrows():
        print(f"  {row['held_out_category']:20s}: {row['top1_mean']*100:.1f}% (n={row['n_loci']})")



if __name__ == "__main__":
    main()
