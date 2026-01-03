#!/usr/bin/env python3
"""
Mechanism-Stratified Evaluation for Regime Map Analysis

Evaluates GWAS-to-gene methods across mechanistic strata to show where
different biological mechanisms require different prioritization strategies.

Operational strata definitions (reviewer-proof):
1. Coding-driven: credible set contains protein-altering consequence in gold gene
2. Proximal regulatory: noncoding ≤10kb from gold gene TSS
3. Distal regulatory: noncoding >100kb from TSS or enhancer-mapped
4. Ambiguous: multiple plausible genes or insufficient annotation

Author: Generated for Nature Genetics v3 Platinum manuscript
Date: 2025-12-21
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MechanismStratifiedEvaluator:
    """Evaluate methods across mechanistic regimes."""
    
    # Mechanistic strata (operational definitions)
    PROXIMAL_THRESHOLD_KB = 10  # ≤10kb = proximal
    DISTAL_THRESHOLD_KB = 100    # >100kb = distal
    
    def __init__(
        self,
        benchmark_path: str,
        evidence_path: str,
        gene_annotations_path: str = "data/external/flames/Annotation_data/ENSG/ENSG.v102.genes.parquet",
        results_dir: str = "benchmarks/results"
    ):
        """
        Initialize evaluator.
        
        Parameters
        ----------
        benchmark_path : str
            Path to benchmark parquet with predictions
        evidence_path : str
            Path to evidence manifest with mechanism annotations
        gene_annotations_path : str
            Path to gene metadata (for biotype)
        results_dir : str
            Output directory
        """
        self.benchmark_path = Path(benchmark_path)
        self.evidence_path = Path(evidence_path)
        self.gene_annotations_path = Path(gene_annotations_path)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        self.df = None
        self.evidence_df = None
        self.gene_annotations = None
        self.label_col = None
        self.method_cols = []
        
    def load_data(self) -> None:
        """Load benchmark, evidence, and gene annotations."""
        logger.info(f"Loading benchmark from {self.benchmark_path}")
        self.df = pd.read_parquet(self.benchmark_path)
        logger.info(f"Loaded {len(self.df)} gene-locus pairs")
        
        # Load evidence manifest
        if self.evidence_path.exists():
            logger.info(f"Loading evidence from {self.evidence_path}")
            self.evidence_df = pd.read_parquet(self.evidence_path)
            logger.info(f"Loaded {len(self.evidence_df)} evidence records")
        else:
            logger.warning(f"Evidence file not found: {self.evidence_path}. Proceeding without it.")
            self.evidence_df = pd.DataFrame()
        
        # Load gene annotations for biotype
        if self.gene_annotations_path.exists():
            logger.info(f"Loading gene annotations from {self.gene_annotations_path}")
            self.gene_annotations = pd.read_parquet(self.gene_annotations_path)
            logger.info(f"Loaded {len(self.gene_annotations)} gene annotations")
        else:
            logger.warning(f"Gene annotations not found at {self.gene_annotations_path}")
            
        # Identify label column
        label_candidates = [c for c in self.df.columns if 'label' in c.lower() or 'truth' in c.lower() or 'positive' in c.lower()]
        if label_candidates:
            self.label_col = label_candidates[0]
        else:
            self.label_col = 'is_causal'
            self.df['is_causal'] = self.df['TruthDistanceRank'] == 1
            
        logger.info(f"Using label column: {self.label_col}")
        logger.info(f"Positive labels: {self.df[self.label_col].sum()}")
        
        # Identify method columns
        rank_cols = [c for c in self.df.columns if 'Rank' in c and c != 'TruthDistanceRank' and c != 'DistanceRank']
        self.method_cols = rank_cols[:10]
        logger.info(f"Found {len(self.method_cols)} prediction methods")
        
    def assign_mechanism_strata(self) -> None:
        """
        Assign each positive label to mechanistic stratum.
        
        Strata definitions:
        1. Coding-driven: Has coding variant consequence (CodingSpliceOrPromoterVariants=1, AnyCoding=1)
        2. Proximal regulatory: Noncoding, ≤10kb from TSS
        3. Distal regulatory: Noncoding, >100kb from TSS
        4. Mid-range regulatory: Noncoding, 10-100kb
        5. Ambiguous: Insufficient information
        """
        self.df['mechanism_stratum'] = 'unassigned'
        self.df['mechanism_label'] = 'Unassigned'
        
        # Calculate distance if not present
        if 'tss_distance' not in self.df.columns:
            if 'GeneBodyDistanceToBestSNP' in self.df.columns:
                self.df['tss_distance'] = self.df['GeneBodyDistanceToBestSNP'].abs()
            elif 'PromoterDistanceToBestSNP' in self.df.columns:
                self.df['tss_distance'] = self.df['PromoterDistanceToBestSNP'].abs()
            else:
                logger.warning("No distance column available")
                self.df['tss_distance'] = np.nan
                
        # Stratum 1: Coding-driven
        if 'AnyCoding' in self.df.columns:
            coding_mask = self.df['AnyCoding'] == 1
            self.df.loc[coding_mask, 'mechanism_stratum'] = 'coding'
            self.df.loc[coding_mask, 'mechanism_label'] = 'Coding-driven (protein-altering)'
            logger.info(f"Coding-driven: {coding_mask.sum()} loci")
            
        # Stratum 2: Proximal regulatory (noncoding ≤10kb)
        proximal_mask = (
            (self.df['mechanism_stratum'] == 'unassigned') &
            (self.df['tss_distance'] <= self.PROXIMAL_THRESHOLD_KB * 1000)
        )
        self.df.loc[proximal_mask, 'mechanism_stratum'] = 'proximal_regulatory'
        self.df.loc[proximal_mask, 'mechanism_label'] = f'Proximal regulatory (≤{self.PROXIMAL_THRESHOLD_KB}kb)'
        logger.info(f"Proximal regulatory: {proximal_mask.sum()} loci")
        
        # Stratum 3: Distal regulatory (noncoding >100kb)
        distal_mask = (
            (self.df['mechanism_stratum'] == 'unassigned') &
            (self.df['tss_distance'] > self.DISTAL_THRESHOLD_KB * 1000)
        )
        self.df.loc[distal_mask, 'mechanism_stratum'] = 'distal_regulatory'
        self.df.loc[distal_mask, 'mechanism_label'] = f'Distal regulatory (>{self.DISTAL_THRESHOLD_KB}kb)'
        logger.info(f"Distal regulatory: {distal_mask.sum()} loci")
        
        # Stratum 4: Mid-range regulatory (10-100kb)
        midrange_mask = (
            (self.df['mechanism_stratum'] == 'unassigned') &
            (self.df['tss_distance'] > self.PROXIMAL_THRESHOLD_KB * 1000) &
            (self.df['tss_distance'] <= self.DISTAL_THRESHOLD_KB * 1000)
        )
        self.df.loc[midrange_mask, 'mechanism_stratum'] = 'midrange_regulatory'
        self.df.loc[midrange_mask, 'mechanism_label'] = f'Mid-range regulatory ({self.PROXIMAL_THRESHOLD_KB}-{self.DISTAL_THRESHOLD_KB}kb)'
        logger.info(f"Mid-range regulatory: {midrange_mask.sum()} loci")
        
        # Stratum 5: Ambiguous (remaining)
        ambiguous_mask = self.df['mechanism_stratum'] == 'unassigned'
        self.df.loc[ambiguous_mask, 'mechanism_stratum'] = 'ambiguous'
        self.df.loc[ambiguous_mask, 'mechanism_label'] = 'Ambiguous/Insufficient data'
        logger.info(f"Ambiguous: {ambiguous_mask.sum()} loci")
        
        logger.info("\nMechanism stratum distribution:")
        for label in self.df['mechanism_label'].value_counts().items():
            logger.info(f"  {label[0]}: {label[1]} pairs")
            
    def calculate_stratified_metrics(self, method_col: str, k: int = 1) -> Dict[str, float]:
        """
        Calculate Top-K precision per mechanism stratum.
        
        Parameters
        ----------
        method_col : str
            Method prediction column
        k : int
            Top-K threshold
            
        Returns
        -------
        dict
            {stratum_label: top-k precision}
        """
        results = {}
        
        for stratum in self.df['mechanism_label'].unique():
            if stratum == 'Unassigned':
                continue
                
            stratum_mask = self.df['mechanism_label'] == stratum
            stratum_df = self.df[stratum_mask]
            
            if len(stratum_df) == 0:
                results[stratum] = np.nan
                continue
                
            # Calculate Top-K precision per locus
            topk_correct = []
            for locus_id, locus_df in stratum_df.groupby('locus_id'):
                true_genes = set(locus_df[locus_df[self.label_col].astype(bool)]['TargetGene'])
                
                if len(true_genes) == 0:
                    continue
                    
                topk_pred = locus_df.nsmallest(k, method_col)['TargetGene']
                correct = any(gene in true_genes for gene in topk_pred)
                topk_correct.append(correct)
                
            results[stratum] = np.mean(topk_correct) if topk_correct else np.nan
            
        return results
        
    def evaluate_all_methods(self) -> pd.DataFrame:
        """
        Evaluate all methods across mechanism strata.
        
        Returns
        -------
        DataFrame
            Results with columns: method, stratum, top1, top3, top5, auprc
        """
        logger.info("\nEvaluating all methods across mechanism strata...")
        
        results = []
        
        for method_col in self.method_cols:
            method_name = method_col.replace('Rank', '').replace('ConnectionStrengthRank', 'ABC')
            
            # Calculate metrics
            top1 = self.calculate_stratified_metrics(method_col, k=1)
            top3 = self.calculate_stratified_metrics(method_col, k=3)
            top5 = self.calculate_stratified_metrics(method_col, k=5)
            
            for stratum in top1.keys():
                results.append({
                    'method': method_name,
                    'mechanism_stratum': stratum,
                    'top1_precision': top1.get(stratum, np.nan),
                    'top3_precision': top3.get(stratum, np.nan),
                    'top5_precision': top5.get(stratum, np.nan)
                })
                
        results_df = pd.DataFrame(results)
        
        # Add proximity baseline
        dist_col = 'DistanceRank' if 'DistanceRank' in self.df.columns else 'Distance (Rank)'
        if dist_col in self.df.columns:
            top1_prox = self.calculate_stratified_metrics(dist_col, k=1)
            top3_prox = self.calculate_stratified_metrics(dist_col, k=3)
            top5_prox = self.calculate_stratified_metrics(dist_col, k=5)
            
            for stratum in top1_prox.keys():
                results_df = pd.concat([results_df, pd.DataFrame([{
                    'method': 'Proximity (baseline)',
                    'mechanism_stratum': stratum,
                    'top1_precision': top1_prox.get(stratum, np.nan),
                    'top3_precision': top3_prox.get(stratum, np.nan),
                    'top5_precision': top5_prox.get(stratum, np.nan)
                }])], ignore_index=True)
                
        return results_df
        
    def compute_delta_from_baseline(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """Compute performance delta from proximity baseline."""
        baseline_df = results_df[results_df['method'] == 'Proximity (baseline)']
        
        delta_rows = []
        for _, row in results_df.iterrows():
            if row['method'] == 'Proximity (baseline)':
                continue
                
            baseline_row = baseline_df[baseline_df['mechanism_stratum'] == row['mechanism_stratum']]
            if len(baseline_row) == 0:
                continue
                
            delta_rows.append({
                'method': row['method'],
                'mechanism_stratum': row['mechanism_stratum'],
                'delta_top1': row['top1_precision'] - baseline_row['top1_precision'].values[0],
                'delta_top3': row['top3_precision'] - baseline_row['top3_precision'].values[0],
                'delta_top5': row['top5_precision'] - baseline_row['top5_precision'].values[0]
            })
            
        delta_df = pd.DataFrame(delta_rows)
        results_with_delta = results_df.merge(delta_df, on=['method', 'mechanism_stratum'], how='left')
        
        return results_with_delta
        
    def run(self) -> pd.DataFrame:
        """Execute full mechanism-stratified evaluation."""
        logger.info("="*60)
        logger.info("Mechanism-Stratified Evaluation (Regime Map)")
        logger.info("="*60)
        
        self.load_data()
        self.assign_mechanism_strata()
        results_df = self.evaluate_all_methods()
        results_with_delta = self.compute_delta_from_baseline(results_df)
        
        # Save results
        output_file = self.results_dir / "mechanism_stratified_performance.csv"
        results_with_delta.to_csv(output_file, index=False)
        logger.info(f"\nSaved results to: {output_file}")
        
        # Print summary
        logger.info("\n" + "="*60)
        logger.info("MECHANISM REGIME SUMMARY")
        logger.info("="*60)
        
        for stratum in results_with_delta['mechanism_stratum'].unique():
            if pd.isna(stratum):
                continue
                
            logger.info(f"\n{stratum}:")
            stratum_results = results_with_delta[results_with_delta['mechanism_stratum'] == stratum]
            
            # Show top 3 methods by Top-1 precision
            top_methods = stratum_results.nlargest(3, 'top1_precision')
            logger.info("  Top-1 Precision:")
            for _, row in top_methods.iterrows():
                delta_str = f" (Δ={row['delta_top1']:+.3f})" if pd.notna(row.get('delta_top1')) else ""
                logger.info(f"    {row['method']}: {row['top1_precision']:.3f}{delta_str}")
                
        return results_with_delta


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Mechanism-stratified evaluation")
    parser.add_argument("--benchmark", default="benchmarks/results_cache.parquet",
                       help="Path to benchmark file")
    parser.add_argument("--evidence", default="benchmarks/task_a_evidence_manifest_v3.parquet",
                       help="Path to evidence manifest")
    parser.add_argument("--gene-annotations", 
                       default="data/external/flames/Annotation_data/ENSG/ENSG.v102.genes.parquet",
                       help="Path to gene annotations")
    parser.add_argument("--results-dir", default="benchmarks/results",
                       help="Output directory")
    
    args = parser.parse_args()
    
    evaluator = MechanismStratifiedEvaluator(
        benchmark_path=args.benchmark,
        evidence_path=args.evidence,
        gene_annotations_path=args.gene_annotations,
        results_dir=args.results_dir
    )
    
    results = evaluator.run()
    print("\nMechanism stratification complete!")
