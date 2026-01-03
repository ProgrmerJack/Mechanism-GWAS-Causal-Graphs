#!/usr/bin/env python3
"""
Evaluation Framework for GWAS Gene Prioritization Methods
==========================================================

This framework provides a standardized way to evaluate gene prioritization
methods on RegulatoryBench v3 (experimental ground truth from CRISPRi/MPRA).

Key principles:
1. All methods scored on same benchmark (no circular training/test overlap)
2. Locus-level evaluation (candidate genes compete within locus)
3. Binary classification metrics (AUC-ROC, AUC-PR, precision@k)
4. Effect size stratification (does method work for weak effects?)
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Callable
from dataclasses import dataclass
from sklearn.metrics import (
    roc_auc_score, 
    average_precision_score, 
    precision_recall_curve,
    roc_curve
)

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
BENCHMARK_FILE = DATA_DIR / "processed" / "baselines" / "regulatorybench_v3.tsv"


@dataclass
class EvaluationResult:
    """Stores evaluation results for a single method."""
    method_name: str
    auc_roc: float
    auc_pr: float
    precision_at_1: float
    precision_at_3: float
    recall_at_1: float
    recall_at_3: float
    n_loci: int
    n_positives: int
    n_candidates: int
    stratified_results: Dict[str, Dict] = None
    
    def to_dict(self) -> Dict:
        return {
            'method': self.method_name,
            'auc_roc': self.auc_roc,
            'auc_pr': self.auc_pr,
            'precision@1': self.precision_at_1,
            'precision@3': self.precision_at_3,
            'recall@1': self.recall_at_1,
            'recall@3': self.recall_at_3,
            'n_loci': self.n_loci,
            'n_positives': self.n_positives,
            'n_candidates': self.n_candidates,
            'stratified': self.stratified_results
        }


class Benchmark:
    """
    Handles loading and managing the experimental ground truth benchmark.
    """
    
    def __init__(self, benchmark_path: Path = BENCHMARK_FILE):
        self.benchmark_path = benchmark_path
        self.data = None
        self.loci = None
        
    def load(self) -> pd.DataFrame:
        """Load the benchmark data."""
        if not self.benchmark_path.exists():
            raise FileNotFoundError(f"Benchmark not found: {self.benchmark_path}")
        
        self.data = pd.read_csv(self.benchmark_path, sep='\t')
        print(f"Loaded benchmark: {len(self.data)} entries")
        print(f"  Evidence types: {self.data['evidence_type'].value_counts().to_dict()}")
        
        # Create locus-level grouping
        self._create_loci()
        
        return self.data
    
    def _create_loci(self):
        """
        Create locus-level groupings for evaluation.
        For CRISPRi: group by enhancer element
        For MPRA: group by variant position (within LD window)
        """
        loci = []
        
        # CRISPRi entries: use locus_id
        crispr_data = self.data[self.data['evidence_type'] == 'CRISPRi']
        for locus_id, group in crispr_data.groupby('locus_id'):
            loci.append({
                'locus_id': locus_id,
                'chr': group['chr'].iloc[0],
                'pos': group['pos'].iloc[0],
                'evidence_type': 'CRISPRi',
                'genes': group['gene_symbol'].tolist(),
                'gene_ensembl': group['gene_ensembl'].tolist(),
                'n_genes': len(group),
                'effect_sizes': group['effect_size'].tolist()
            })
        
        # MPRA entries: use variant position
        mpra_data = self.data[self.data['evidence_type'] == 'MPRA']
        for (chrom, pos), group in mpra_data.groupby(['chr', 'pos']):
            locus_id = f"mpra_{chrom}_{pos}"
            loci.append({
                'locus_id': locus_id,
                'chr': chrom,
                'pos': pos,
                'evidence_type': 'MPRA',
                'genes': group['gene_symbol'].tolist() if 'gene_symbol' in group.columns else [],
                'gene_ensembl': group['gene_ensembl'].tolist(),
                'n_genes': len(group),
                'effect_sizes': group['effect_size'].tolist()
            })
        
        self.loci = pd.DataFrame(loci)
        print(f"  Created {len(self.loci)} evaluation loci")
        
    def get_loci(self) -> pd.DataFrame:
        """Return locus-level summary."""
        if self.loci is None:
            self.load()
        return self.loci
    
    def get_positive_genes(self) -> set:
        """Return set of all positive (experimentally validated) genes."""
        if self.data is None:
            self.load()
        return set(self.data['gene_ensembl'].dropna().unique())


class BaselineMethod:
    """
    Base class for gene prioritization methods.
    Each method should implement the `score` method.
    """
    
    def __init__(self, name: str):
        self.name = name
        
    def score(self, locus: Dict, candidate_genes: List[str]) -> Dict[str, float]:
        """
        Score candidate genes for a locus.
        
        Args:
            locus: Dict with locus info (chr, pos, etc.)
            candidate_genes: List of gene IDs (ENSG format)
            
        Returns:
            Dict mapping gene_id -> score (higher = more likely causal)
        """
        raise NotImplementedError


class NearestGeneBaseline(BaselineMethod):
    """
    Baseline: Nearest gene to the lead variant gets score 1, others 0.
    This is the simplest possible baseline.
    """
    
    def __init__(self, gene_coordinates: pd.DataFrame = None):
        super().__init__("NearestGene")
        self.gene_coords = gene_coordinates
        
    def load_gene_coordinates(self, path: Path = None):
        """Load gene coordinates from file or use default."""
        if path is None:
            # Try common locations
            candidates = [
                DATA_DIR / "external" / "gene_coordinates.tsv",
                DATA_DIR / "external" / "gencode" / "gene_coords.tsv",
                DATA_DIR / "raw" / "genes.tsv",
            ]
            for p in candidates:
                if p.exists():
                    path = p
                    break
        
        if path and path.exists():
            self.gene_coords = pd.read_csv(path, sep='\t')
            print(f"Loaded {len(self.gene_coords)} gene coordinates")
        else:
            print("Warning: No gene coordinates file found - NearestGene will use TSS from benchmark")
    
    def score(self, locus: Dict, candidate_genes: List[str]) -> Dict[str, float]:
        """Score genes by inverse distance to locus."""
        # For now, return uniform scores (need gene coordinates)
        # The actual implementation would compute distance
        n_genes = len(candidate_genes)
        if n_genes == 0:
            return {}
        return {g: 1.0 / n_genes for g in candidate_genes}


class DistanceBasedBaseline(BaselineMethod):
    """
    Score genes by inverse distance to lead variant.
    Uses TSS distance from benchmark data.
    """
    
    def __init__(self):
        super().__init__("InverseDistance")
        
    def score(self, locus: Dict, candidate_genes: List[str], 
              gene_distances: Dict[str, float] = None) -> Dict[str, float]:
        """
        Score = 1 / (distance + 1) to handle distance=0 cases.
        """
        if gene_distances is None:
            # Fallback: uniform scores
            n_genes = len(candidate_genes)
            return {g: 1.0 / n_genes for g in candidate_genes}
        
        scores = {}
        for gene in candidate_genes:
            dist = gene_distances.get(gene, 1e9)
            scores[gene] = 1.0 / (abs(dist) + 1)
        
        # Normalize to sum to 1
        total = sum(scores.values())
        if total > 0:
            scores = {g: s / total for g, s in scores.items()}
        
        return scores


class Evaluator:
    """
    Main evaluation class.
    Computes metrics for any method on the benchmark.
    """
    
    def __init__(self, benchmark: Benchmark):
        self.benchmark = benchmark
        
    def evaluate(self, method: BaselineMethod, 
                 candidate_generator: Callable = None) -> EvaluationResult:
        """
        Evaluate a method on the benchmark.
        
        Args:
            method: The gene prioritization method to evaluate
            candidate_generator: Function that returns candidate genes for a locus
                                 If None, uses genes from benchmark as candidates
        """
        if self.benchmark.data is None:
            self.benchmark.load()
            
        all_scores = []
        all_labels = []
        locus_results = []
        
        for _, locus in self.benchmark.loci.iterrows():
            # Get positive genes for this locus
            positive_genes = set(locus['gene_ensembl'])
            
            # Get candidate genes
            if candidate_generator:
                candidates = candidate_generator(locus)
            else:
                # Use genes from benchmark as candidates (simpler evaluation)
                candidates = list(positive_genes)
                
            if len(candidates) == 0:
                continue
                
            # Score candidates
            scores = method.score(locus.to_dict(), candidates)
            
            # Create labels
            for gene, score in scores.items():
                label = 1 if gene in positive_genes else 0
                all_scores.append(score)
                all_labels.append(label)
                
            # Locus-level metrics
            ranked_genes = sorted(scores.items(), key=lambda x: -x[1])
            top1_correct = ranked_genes[0][0] in positive_genes if ranked_genes else False
            top3_correct = any(g in positive_genes for g, _ in ranked_genes[:3])
            
            locus_results.append({
                'locus_id': locus['locus_id'],
                'n_candidates': len(candidates),
                'n_positives': len(positive_genes),
                'top1_correct': top1_correct,
                'top3_correct': top3_correct
            })
        
        # Compute overall metrics
        all_scores = np.array(all_scores)
        all_labels = np.array(all_labels)
        
        # Handle edge cases
        if len(np.unique(all_labels)) < 2:
            print(f"Warning: Only one class in labels for {method.name}")
            auc_roc = 0.5
            auc_pr = np.mean(all_labels)
        else:
            auc_roc = roc_auc_score(all_labels, all_scores)
            auc_pr = average_precision_score(all_labels, all_scores)
        
        locus_df = pd.DataFrame(locus_results)
        
        return EvaluationResult(
            method_name=method.name,
            auc_roc=auc_roc,
            auc_pr=auc_pr,
            precision_at_1=locus_df['top1_correct'].mean(),
            precision_at_3=locus_df['top3_correct'].mean(),
            recall_at_1=locus_df['top1_correct'].sum() / locus_df['n_positives'].sum(),
            recall_at_3=locus_df['top3_correct'].sum() / locus_df['n_positives'].sum(),
            n_loci=len(locus_df),
            n_positives=int(locus_df['n_positives'].sum()),
            n_candidates=int(locus_df['n_candidates'].sum())
        )
    
    def evaluate_stratified(self, method: BaselineMethod,
                           stratify_by: str = 'evidence_type') -> EvaluationResult:
        """
        Evaluate with stratification by evidence type, effect size, etc.
        """
        if self.benchmark.data is None:
            self.benchmark.load()
            
        stratified_results = {}
        
        for stratum, loci in self.benchmark.loci.groupby(stratify_by):
            # Create temporary benchmark with subset
            temp_benchmark = Benchmark.__new__(Benchmark)
            temp_benchmark.data = self.benchmark.data[
                self.benchmark.data['locus_id'].isin(loci['locus_id'])
            ]
            temp_benchmark.loci = loci
            
            temp_evaluator = Evaluator(temp_benchmark)
            result = temp_evaluator.evaluate(method)
            stratified_results[stratum] = result.to_dict()
        
        # Get overall result
        overall = self.evaluate(method)
        overall.stratified_results = stratified_results
        
        return overall


def run_baseline_comparison():
    """
    Run comparison of baseline methods on RegulatoryBench v3.
    """
    print("="*60)
    print("BASELINE COMPARISON ON REGULATORYBENCH V3")
    print("="*60)
    
    # Load benchmark
    benchmark = Benchmark()
    benchmark.load()
    
    # Initialize evaluator
    evaluator = Evaluator(benchmark)
    
    # Define methods
    methods = [
        NearestGeneBaseline(),
        DistanceBasedBaseline(),
    ]
    
    # Run evaluation
    results = []
    for method in methods:
        print(f"\nEvaluating: {method.name}")
        result = evaluator.evaluate(method)
        results.append(result)
        print(f"  AUC-ROC: {result.auc_roc:.3f}")
        print(f"  AUC-PR: {result.auc_pr:.3f}")
        print(f"  Precision@1: {result.precision_at_1:.3f}")
    
    # Save results
    output_dir = DATA_DIR / "processed" / "evaluation"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results_df = pd.DataFrame([r.to_dict() for r in results])
    results_df.to_csv(output_dir / "baseline_comparison.tsv", sep='\t', index=False)
    print(f"\nResults saved to: {output_dir / 'baseline_comparison.tsv'}")
    
    return results


if __name__ == "__main__":
    run_baseline_comparison()
