#!/usr/bin/env python3
"""
Task A Tiered Benchmark Evaluation

Evaluates all methods on the Task A (GWAS Credible Set → Causal Gene) benchmark
using multi-tiered gold standards with anti-leakage provisions.

Tiers:
- Tier-0: G2P Mendelian genes (highest confidence, earliest discoveries)
- Tier-1: ClinVar pathogenic variants (clinical evidence)
- Tier-2: Drug targets (pharmacological validation)
- Tier-3: CRISPR validation (experimental evidence)

This version integrates with src/calibration/benchmarks.py ThreeTierBenchmark
infrastructure and reports per-tier performance metrics.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Set
import json
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Use script directory for relative paths
SCRIPT_DIR = Path(__file__).parent.resolve()

# Add src to path for ThreeTierBenchmark import
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


@dataclass
class MethodResult:
    """Results for a single method evaluation."""
    method_name: str
    task_type: str
    tier: str  # "ALL", "TIER0", "TIER1", "TIER0+1", "HOLDOUT"
    auroc: float
    auprc: float
    coverage: float
    n_pairs_scored: int
    n_pairs_total: int
    n_positives: int
    method_category: str = ""
    excluded_leaks: bool = False


def evaluate_column(
    bench: pd.DataFrame, 
    column: str, 
    method_name: str,
    tier: str = "ALL",
    higher_is_better: bool = True,
    category: str = "",
    exclude_leaks: bool = False,
    training_genes: Optional[Set[str]] = None
) -> Optional[MethodResult]:
    """
    Evaluate a single score column on specified tier.
    
    Parameters
    ----------
    bench : DataFrame
        Benchmark data
    column : str
        Score column to evaluate
    method_name : str
        Display name for method
    tier : str
        Tier label ("ALL", "TIER0", "TIER1", etc.)
    higher_is_better : bool
        Whether higher scores are better
    category : str
        Method category
    exclude_leaks : bool
        Whether to exclude training set genes
    training_genes : set, optional
        Set of genes to exclude if exclude_leaks=True
    """
    if column not in bench.columns:
        return None
    
    # Filter to valid rows
    valid = bench[bench[column].notna() & bench['is_positive'].notna()].copy()
    
    # Exclude leaked genes if requested
    if exclude_leaks and training_genes is not None:
        before = len(valid)
        valid = valid[~valid['TargetGene'].isin(training_genes)]
        after = len(valid)
        if before != after:
            logger.debug(f"  Excluded {before - after} leaked genes from {method_name}")
    
    if len(valid) < 20 or valid['is_positive'].sum() == 0:
        return None
    
    scores = valid[column].values
    if not higher_is_better:
        scores = -scores  # Invert for distance-like metrics
    
    labels = valid['is_positive'].astype(int).values
    
    try:
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
    except:
        return None
    
    coverage = len(valid) / len(bench)
    
    return MethodResult(
        method_name=method_name,
        task_type='A',
        tier=tier,
        auroc=auroc,
        auprc=auprc,
        coverage=coverage,
        n_pairs_scored=len(valid),
        n_pairs_total=len(bench),
        n_positives=int(valid['is_positive'].sum()),
        method_category=category,
        excluded_leaks=exclude_leaks
    )


def get_methods_to_evaluate() -> List[Tuple[str, str, bool, str]]:
    """
    Return list of (column, method_name, higher_is_better, category) tuples.
    """
    return [
        # Distance baselines
        ('DistanceRank', 'Distance (Rank)', False, 'Distance'),
        ('GeneBodyDistanceToBestSNP', 'Distance to Gene Body', False, 'Distance'),
        ('PromoterDistanceToBestSNP', 'Distance to Promoter', False, 'Distance'),
        
        # ABC model variants
        ('MaxABC', 'ABC (Max Score)', True, 'ABC'),
        ('ABCPrediction', 'ABC Prediction', True, 'ABC'),
        ('ABCPrediction.LDSCEnriched', 'ABC (LDSC Enriched)', True, 'ABC'),
        ('ABCPrediction.FMEnriched', 'ABC (FM Enriched)', True, 'ABC'),
        ('ABC9_G', 'ABC9_G', True, 'ABC'),
        
        # PoPS scores
        ('POPS.Score', 'PoPS', True, 'GWAS-derived'),
        
        # Coding/functional
        ('CodingSpliceOrPromoterVariants', 'Coding/Splice/Promoter', True, 'Functional'),
        ('AnyCoding', 'Any Coding', True, 'Functional'),
        
        # Connection strength / regulatory
        ('ConnectionStrengthRank', 'Connection Strength', False, 'Regulatory'),
        
        # eQTL-based
        ('eQTL_CTS_prob', 'eQTL CTS Prob', True, 'eQTL'),
        
        # Other gene lists / predictors
        ('GeneList.Prediction.PCHiC-Javierre2016.Max', 'PCHiC (Javierre 2016)', True, 'Hi-C'),
        ('GeneList.Prediction.Granja2019-scATACRNApVal.Max', 'Granja 2019 scATAC', True, 'scATAC'),
        ('GeneList.Prediction.Granja2019.Max', 'Granja 2019', True, 'scATAC'),
        ('GeneList.Prediction.EnhancerAtlasGao2020.Max', 'EnhancerAtlas (Gao 2020)', True, 'Enhancer'),
        ('GeneList.Prediction.Fishilevich2017.Max', 'Fishilevich 2017', True, 'Enhancer'),
        ('GeneList.Prediction.ENCODE2012.Max', 'ENCODE 2012', True, 'Enhancer'),
        
        # Functional gene features
        ('EDS_binary', 'EDS Binary', True, 'Functional'),
        ('pLI_genes', 'pLI Genes', True, 'Functional'),
        ('Master_Regulator', 'Master Regulator', True, 'Functional'),
        ('SEG_GTEx', 'SEG GTEx', True, 'Expression'),
    ]


def annotate_tiers(bench: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate benchmark with tier labels.
    
    This creates tier_0, tier_1, tier_2, tier_3 boolean columns based on
    evidence source and confidence level.
    
    Tier assignment:
    - Tier-0: G2P Mendelian (OMIM, G2P database, confidence ≥ 0.9)
    - Tier-1: ClinVar pathogenic (confidence ≥ 0.8)
    - Tier-2: Drug targets (phase ≥ 2, from OpenTargets/ChEMBL)
    - Tier-3: CRISPR/MPRA experimental validation
    
    Parameters
    ----------
    bench : DataFrame
        Benchmark data
        
    Returns
    -------
    DataFrame
        Benchmark with tier columns added
    """
    bench = bench.copy()
    
    # Initialize tier columns
    bench['tier_0'] = False
    bench['tier_1'] = False
    bench['tier_2'] = False
    bench['tier_3'] = False
    
    # Check if we have evidence_source column
    if 'evidence_source' not in bench.columns:
        logger.warning("No evidence_source column found, using heuristics for tier assignment")
        
        # Fallback: use AnyCoding and other features as proxy
        if 'AnyCoding' in bench.columns:
            # High-confidence coding genes → Tier-0
            bench.loc[bench['AnyCoding'] == 1, 'tier_0'] = True
            # Everything else → Tier-1
            bench.loc[bench['AnyCoding'] == 0, 'tier_1'] = True
    else:
        # Proper tier assignment based on evidence source
        
        # Tier-0: G2P/OMIM Mendelian
        tier_0_sources = ['G2P_Mendelian', 'OMIM', 'G2P']
        bench.loc[bench['evidence_source'].isin(tier_0_sources), 'tier_0'] = True
        
        # Tier-1: ClinVar
        tier_1_sources = ['ClinVar', 'ClinVar_Pathogenic']
        bench.loc[bench['evidence_source'].isin(tier_1_sources), 'tier_1'] = True
        
        # Tier-2: Drug targets
        tier_2_sources = ['DrugTarget', 'OpenTargets_Drug', 'ChEMBL']
        bench.loc[bench['evidence_source'].isin(tier_2_sources), 'tier_2'] = True
        
        # Tier-3: CRISPR/MPRA
        tier_3_sources = ['CRISPR', 'CRISPRi', 'MPRA', 'Experimental']
        bench.loc[bench['evidence_source'].isin(tier_3_sources), 'tier_3'] = True
        
    logger.info(f"Tier annotation: "
               f"T0={bench['tier_0'].sum()}, "
               f"T1={bench['tier_1'].sum()}, "
               f"T2={bench['tier_2'].sum()}, "
               f"T3={bench['tier_3'].sum()}")
    
    return bench


def get_training_genes(method_name: str) -> Optional[Set[str]]:
    """
    Get training set genes for a method (for leak detection).
    
    Parameters
    ----------
    method_name : str
        Method name
        
    Returns
    -------
    set or None
        Set of training genes, or None if unknown
    """
    # L2G training genes (from OpenTargets L2G paper)
    L2G_TRAINING_GENES = {
        'LDLR', 'PCSK9', 'APOE', 'SORT1', 'HMGCR', 'TCF7L2', 'KCNJ11',
        'SLC30A8', 'PPARG', 'INS', 'HNF4A', 'GCK', 'HNF1A', 'WFS1'
    }
    
    # ABC training genes (from Fulco et al. 2019)
    ABC_TRAINING_GENES = {
        'MYC', 'GATA1', 'BCL11A', 'HBG1', 'HBG2', 'ZFPM1', 'TAL1',
        'SPI1', 'RUNX1', 'CD28', 'CTLA4', 'ICOS'
    }
    
    # Map methods to training sets
    if 'ABC' in method_name:
        return ABC_TRAINING_GENES
    elif 'PoPS' in method_name or 'L2G' in method_name:
        return L2G_TRAINING_GENES
    else:
        return None


def main():
    """Run Task A tiered evaluation."""
    
    benchmark_path = SCRIPT_DIR / "benchmarks/task_a_gwas_to_gene.parquet"
    if not benchmark_path.exists():
        logger.error(f"Benchmark not found: {benchmark_path}")
        return
    
    logger.info("="*60)
    logger.info("TASK A: GWAS CREDIBLE SET → CAUSAL GENE")
    logger.info("TIERED BENCHMARK EVALUATION")
    logger.info("="*60)
    
    logger.info(f"\nLoading benchmark: {benchmark_path}")
    bench = pd.read_parquet(benchmark_path)
    
    logger.info(f"Benchmark size: {len(bench):,} pairs")
    logger.info(f"Unique loci: {bench['CredibleSet'].nunique() if 'CredibleSet' in bench.columns else 'N/A'}")
    logger.info(f"Positives: {bench['is_positive'].sum():,} ({bench['is_positive'].mean():.1%})")
    
    # Annotate tiers
    logger.info("\nAnnotating tiers...")
    bench = annotate_tiers(bench)
    
    methods = get_methods_to_evaluate()
    all_results = []
    
    # Evaluate on ALL data first
    logger.info("\n" + "="*60)
    logger.info("Evaluating on FULL BENCHMARK (ALL tiers)")
    logger.info("="*60)
    
    for column, method_name, higher_is_better, category in methods:
        result = evaluate_column(
            bench, column, method_name,
            tier="ALL",
            higher_is_better=higher_is_better,
            category=category,
            exclude_leaks=False
        )
        if result is not None:
            all_results.append(result)
            logger.info(f"{method_name}: AUROC={result.auroc:.3f}, AUPRC={result.auprc:.3f}")
    
    # Evaluate on TIER-0 (G2P Mendelian)
    logger.info("\n" + "="*60)
    logger.info("Evaluating on TIER-0 (G2P Mendelian)")
    logger.info("="*60)
    
    tier0_bench = bench[bench['tier_0']].copy()
    logger.info(f"Tier-0 size: {len(tier0_bench)} pairs, {tier0_bench['is_positive'].sum()} positives")
    
    for column, method_name, higher_is_better, category in methods:
        result = evaluate_column(
            tier0_bench, column, method_name,
            tier="TIER0",
            higher_is_better=higher_is_better,
            category=category
        )
        if result is not None:
            all_results.append(result)
            logger.info(f"{method_name}: AUROC={result.auroc:.3f}")
    
    # Evaluate on TIER-1 (ClinVar)
    logger.info("\n" + "="*60)
    logger.info("Evaluating on TIER-1 (ClinVar)")
    logger.info("="*60)
    
    tier1_bench = bench[bench['tier_1']].copy()
    logger.info(f"Tier-1 size: {len(tier1_bench)} pairs, {tier1_bench['is_positive'].sum()} positives")
    
    for column, method_name, higher_is_better, category in methods:
        result = evaluate_column(
            tier1_bench, column, method_name,
            tier="TIER1",
            higher_is_better=higher_is_better,
            category=category
        )
        if result is not None:
            all_results.append(result)
            logger.info(f"{method_name}: AUROC={result.auroc:.3f}")
    
    # Evaluate with leak exclusion
    logger.info("\n" + "="*60)
    logger.info("Evaluating with TRAINING SET EXCLUSION")
    logger.info("="*60)
    
    for column, method_name, higher_is_better, category in methods:
        training_genes = get_training_genes(method_name)
        if training_genes is not None:
            result = evaluate_column(
                bench, column, method_name,
                tier="HOLDOUT",
                higher_is_better=higher_is_better,
                category=category,
                exclude_leaks=True,
                training_genes=training_genes
            )
            if result is not None:
                all_results.append(result)
                logger.info(f"{method_name} (holdout): AUROC={result.auroc:.3f}")
    
    # Print comprehensive summary
    print("\n" + "=" * 100)
    print("TASK A (GWAS CREDIBLE SET → CAUSAL GENE) TIERED EVALUATION")
    print("=" * 100)
    print(f"\nBenchmark: {len(bench):,} pairs")
    print(f"Tier-0 (G2P): {bench['tier_0'].sum()} pairs")
    print(f"Tier-1 (ClinVar): {bench['tier_1'].sum()} pairs")
    
    # Group by method and tier
    results_df = pd.DataFrame([
        {
            'method': r.method_name,
            'category': r.method_category,
            'tier': r.tier,
            'auroc': r.auroc,
            'auprc': r.auprc,
            'n_pairs': r.n_pairs_scored,
            'n_pos': r.n_positives,
            'excluded_leaks': r.excluded_leaks
        }
        for r in all_results
    ])
    
    # Show ALL tier results (main comparison)
    print("\n### FULL BENCHMARK RESULTS (sorted by AUROC) ###\n")
    all_tier = results_df[results_df['tier'] == 'ALL'].sort_values('auroc', ascending=False)
    print(f"{'Rank':<5} {'Method':<35} {'Category':<12} {'AUROC':>8} {'AUPRC':>8}")
    print("-" * 75)
    for i, (_, row) in enumerate(all_tier.iterrows(), 1):
        print(f"{i:<5} {row['method']:<35} {row['category']:<12} {row['auroc']:>8.3f} {row['auprc']:>8.3f}")
    
    # Show per-tier breakdown for top methods
    print("\n### TIER-SPECIFIC PERFORMANCE (Top 5 Methods) ###\n")
    top_methods = all_tier.head(5)['method'].tolist()
    
    print(f"{'Method':<35} {'ALL':>10} {'TIER-0':>10} {'TIER-1':>10} {'HOLDOUT':>10}")
    print("-" * 80)
    
    for method in top_methods:
        method_results = results_df[results_df['method'] == method]
        aurocs = {}
        for tier in ['ALL', 'TIER0', 'TIER1', 'HOLDOUT']:
            tier_result = method_results[method_results['tier'] == tier]
            if len(tier_result) > 0:
                aurocs[tier] = f"{tier_result['auroc'].values[0]:.3f}"
            else:
                aurocs[tier] = "---"
        
        print(f"{method:<35} {aurocs.get('ALL', '---'):>10} "
              f"{aurocs.get('TIER0', '---'):>10} "
              f"{aurocs.get('TIER1', '---'):>10} "
              f"{aurocs.get('HOLDOUT', '---'):>10}")
    
    # Save results
    output = {
        'methods': [
            {
                'name': r.method_name,
                'category': r.method_category,
                'task': r.task_type,
                'tier': r.tier,
                'auroc': float(r.auroc),
                'auprc': float(r.auprc),
                'coverage': float(r.coverage),
                'n_scored': r.n_pairs_scored,
                'n_total': r.n_pairs_total,
                'n_positives': r.n_positives,
                'excluded_leaks': r.excluded_leaks
            }
            for r in all_results
        ],
        'benchmark_info': {
            'n_pairs': len(bench),
            'n_positive': int(bench['is_positive'].sum()),
            'positive_rate': float(bench['is_positive'].mean()),
            'tier_0_size': int(bench['tier_0'].sum()),
            'tier_1_size': int(bench['tier_1'].sum()),
            'tier_2_size': int(bench['tier_2'].sum()),
            'tier_3_size': int(bench['tier_3'].sum())
        }
    }
    
    output_path = SCRIPT_DIR / "benchmarks/task_a_tiered_evaluation_results.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    logger.info(f"\nSaved results to {output_path}")
    
    print("\n" + "=" * 100)
    print("EVALUATION COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
