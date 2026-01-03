#!/usr/bin/env python3
"""
Benchmark gene prioritization methods against Open Targets gold standards.

Implements proper metrics:
- Top-1 accuracy (was #1 ranked gene correct?)
- Mean Reciprocal Rank (MRR)
- Recall@k (Recall@5, Recall@10, Recall@20)
- Expected Calibration Error (ECE) - probability calibration
- Brier score - calibration metric
- Per-trait stratification
- Per-confidence stratification

This is the CORRECT way to benchmark, not just "did we find it in top 20".
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import json
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Paths
GOLD_STANDARDS_FILE = Path("data/external/open_targets/curated_gold_standards.tsv")
LIPID_GOLD_STANDARDS = Path("data/external/open_targets/lipid_gold_standards.tsv")
TIER1_RESULTS = Path("outputs/tier1/tier1_scores.tsv")  # From previous work
OUTPUT_DIR = Path("outputs/benchmarking")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

def load_gold_standards(confidence_levels=None):
    """Load gold standard gene-locus pairs."""
    print(f"Loading gold standards from {GOLD_STANDARDS_FILE}")
    gs = pd.read_csv(GOLD_STANDARDS_FILE, sep='\t')
    
    if confidence_levels:
        gs = gs[gs['confidence'].isin(confidence_levels)]
    
    print(f"Loaded {len(gs)} gold standards")
    print(f"  High confidence: {len(gs[gs['confidence']=='High'])}")
    print(f"  Medium confidence: {len(gs[gs['confidence']=='Medium'])}")
    
    return gs

def calculate_top_k_accuracy(ranks, k):
    """Calculate what fraction of genes were found in top-k."""
    return np.mean(ranks <= k)

def calculate_mrr(ranks):
    """
    Calculate Mean Reciprocal Rank.
    MRR = average of (1 / rank of correct gene)
    
    If correct gene ranked #1: contributes 1.0
    If correct gene ranked #2: contributes 0.5
    If correct gene ranked #10: contributes 0.1
    """
    reciprocal_ranks = 1.0 / ranks
    return np.mean(reciprocal_ranks)

def calculate_ece(predicted_probs, binary_outcomes, n_bins=10):
    """
    Calculate Expected Calibration Error.
    
    For predictions grouped into bins by probability:
    ECE = sum over bins of: |bin accuracy - bin avg confidence| * bin_size
    
    A well-calibrated model should have: if it predicts 70% probability, 
    70% of those predictions should be correct.
    """
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    ece = 0.0
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        in_bin = (predicted_probs > bin_lower) & (predicted_probs <= bin_upper)
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            accuracy_in_bin = binary_outcomes[in_bin].mean()
            avg_confidence_in_bin = predicted_probs[in_bin].mean()
            ece += np.abs(avg_confidence_in_bin - accuracy_in_bin) * prop_in_bin
    
    return ece

def calculate_brier_score(predicted_probs, binary_outcomes):
    """
    Calculate Brier score: mean squared error of probability predictions.
    Lower is better (0 = perfect calibration).
    """
    return np.mean((predicted_probs - binary_outcomes) ** 2)

def rank_genes_for_locus(method_scores: pd.DataFrame, gene_id: str) -> Tuple[int, float]:
    """
    Given scores for all genes in a locus, find rank of target gene.
    
    Returns:
        (rank, score) where rank=1 means top-ranked gene
    """
    # Sort by score descending
    sorted_genes = method_scores.sort_values('score', ascending=False)
    
    # Find rank of target gene (1-indexed)
    try:
        rank = sorted_genes.index.tolist().index(gene_id) + 1
        score = sorted_genes.loc[gene_id, 'score']
        return rank, score
    except (ValueError, KeyError):
        # Gene not scored by this method
        return np.inf, 0.0

def benchmark_method(method_name: str, 
                     gold_standards: pd.DataFrame,
                     method_scores_dir: Path) -> Dict:
    """
    Benchmark a single method against gold standards.
    
    For each gold standard locus:
    1. Load all gene scores for that locus
    2. Find rank of gold standard gene
    3. Record score and rank
    
    Then calculate aggregate metrics.
    """
    print(f"\n{'='*80}")
    print(f"Benchmarking: {method_name}")
    print(f"{'='*80}")
    
    ranks = []
    scores = []
    binary_outcomes = []  # 1 if top-ranked, 0 otherwise
    
    # Track per-trait performance
    trait_ranks = defaultdict(list)
    
    for idx, row in gold_standards.iterrows():
        gene_id = row['gene_id']
        trait = row['trait_name']
        locus_key = f"{row['chromosome_GRCh38']}_{row['position_GRCh38']}"
        
        # Load scores for this locus
        # (This assumes each locus has a file with scores for all genes)
        # You'll need to adapt this to your actual file structure
        locus_file = method_scores_dir / f"{locus_key}_scores.tsv"
        
        if not locus_file.exists():
            # Method didn't score this locus
            continue
        
        locus_scores = pd.read_csv(locus_file, sep='\t', index_col=0)
        rank, score = rank_genes_for_locus(locus_scores, gene_id)
        
        if rank != np.inf:
            ranks.append(rank)
            scores.append(score)
            binary_outcomes.append(1 if rank == 1 else 0)
            trait_ranks[trait].append(rank)
    
    if len(ranks) == 0:
        print(f"⚠ No loci scored by {method_name}")
        return {}
    
    ranks = np.array(ranks)
    scores = np.array(scores)
    binary_outcomes = np.array(binary_outcomes)
    
    # Calculate metrics
    metrics = {
        'method': method_name,
        'n_loci': len(ranks),
        'top1_accuracy': calculate_top_k_accuracy(ranks, 1),
        'top5_accuracy': calculate_top_k_accuracy(ranks, 5),
        'top10_accuracy': calculate_top_k_accuracy(ranks, 10),
        'top20_accuracy': calculate_top_k_accuracy(ranks, 20),
        'mrr': calculate_mrr(ranks),
        'ece': calculate_ece(scores, binary_outcomes),
        'brier_score': calculate_brier_score(scores, binary_outcomes),
        'median_rank': np.median(ranks),
        'mean_rank': np.mean(ranks)
    }
    
    # Print results
    print(f"\nLoci evaluated: {metrics['n_loci']}")
    print(f"\nAccuracy Metrics:")
    print(f"  Top-1 Accuracy:  {metrics['top1_accuracy']:.3f}  (correct gene ranked #1)")
    print(f"  Top-5 Accuracy:  {metrics['top5_accuracy']:.3f}  (correct gene in top 5)")
    print(f"  Top-10 Accuracy: {metrics['top10_accuracy']:.3f}  (correct gene in top 10)")
    print(f"  Top-20 Accuracy: {metrics['top20_accuracy']:.3f}  (correct gene in top 20)")
    print(f"\nRanking Metrics:")
    print(f"  Mean Reciprocal Rank (MRR): {metrics['mrr']:.3f}")
    print(f"  Median rank: {metrics['median_rank']:.1f}")
    print(f"  Mean rank: {metrics['mean_rank']:.1f}")
    print(f"\nCalibration Metrics:")
    print(f"  Expected Calibration Error (ECE): {metrics['ece']:.3f}")
    print(f"  Brier Score: {metrics['brier_score']:.3f}")
    
    # Per-trait breakdown for top traits
    print(f"\nPer-Trait Performance (Top 10 traits by # loci):")
    trait_summary = []
    for trait, trait_rank_list in sorted(trait_ranks.items(), 
                                         key=lambda x: len(x[1]), 
                                         reverse=True)[:10]:
        trait_arr = np.array(trait_rank_list)
        trait_summary.append({
            'trait': trait,
            'n_loci': len(trait_arr),
            'top1_acc': calculate_top_k_accuracy(trait_arr, 1),
            'mrr': calculate_mrr(trait_arr)
        })
        print(f"  {trait:40s}: n={len(trait_arr):3d}, Top-1={trait_summary[-1]['top1_acc']:.3f}, MRR={trait_summary[-1]['mrr']:.3f}")
    
    metrics['per_trait'] = trait_summary
    
    return metrics

def create_comparison_table(all_metrics: List[Dict]) -> pd.DataFrame:
    """Create comparison table across methods."""
    comparison = []
    
    for m in all_metrics:
        comparison.append({
            'Method': m['method'],
            'N Loci': m['n_loci'],
            'Top-1': f"{m['top1_accuracy']:.3f}",
            'Top-5': f"{m['top5_accuracy']:.3f}",
            'Top-10': f"{m['top10_accuracy']:.3f}",
            'Top-20': f"{m['top20_accuracy']:.3f}",
            'MRR': f"{m['mrr']:.3f}",
            'ECE': f"{m['ece']:.3f}",
            'Brier': f"{m['brier_score']:.3f}",
            'Median Rank': f"{m['median_rank']:.1f}"
        })
    
    return pd.DataFrame(comparison)

def main():
    print("="*80)
    print("GENE PRIORITIZATION BENCHMARKING")
    print("Open Targets Gold Standards (Mountjoy et al. 2021)")
    print("="*80)
    
    # Load gold standards
    gs_all = load_gold_standards()
    gs_high = load_gold_standards(confidence_levels=['High'])
    gs_lipid = pd.read_csv(LIPID_GOLD_STANDARDS, sep='\t')
    
    print(f"\n{len(gs_lipid)} gold standards for lipid traits")
    
    # TODO: Benchmark each method
    # For now, this is a template showing what needs to be done
    
    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("""
To complete this benchmarking, we need to:

1. Format our method's scores to match expected structure:
   - One file per locus: {chr}_{pos}_scores.tsv
   - Columns: gene_id, score
   - All candidate genes in locus ranked by score

2. Load and format baseline scores similarly:
   - Open Targets L2G (if we can get official scores)
   - Our cS2G-inspired proxy
   - Our FLAMES approximation
   - PoPS
   - Effector Index

3. Run benchmark_method() for each:
   - Our probabilistic causal graph method
   - L2G (official if available, else skip)
   - cS2G-inspired proxy
   - FLAMES approximation
   - PoPS
   - Effector Index
   - Simple nearest-gene baseline

4. Create comparison table and visualizations

5. Stratify by:
   - Confidence level (High vs Medium)
   - Trait type (lipid vs disease vs other)
   - Evidence class (drug target vs functional vs expert curated)
   - Locus complexity (# candidate genes)

This will give us Nature Biotechnology-grade benchmarking with:
✓ Proper metrics (not just Recall@20)
✓ Calibration assessment (ECE, Brier)
✓ Stratified analysis
✓ Fair comparison to official published methods
    """)
    
    # Save template for method comparison
    template = {
        'benchmark_structure': {
            'gold_standards': str(GOLD_STANDARDS_FILE),
            'n_total_loci': len(gs_all),
            'n_high_confidence': len(gs_high),
            'n_lipid_loci': len(gs_lipid),
            'required_metrics': [
                'top1_accuracy', 'top5_accuracy', 'top10_accuracy', 'top20_accuracy',
                'mrr', 'ece', 'brier_score', 'median_rank'
            ],
            'stratifications': [
                'by_confidence', 'by_trait', 'by_evidence_class'
            ]
        }
    }
    
    template_file = OUTPUT_DIR / "benchmark_template.json"
    with open(template_file, 'w') as f:
        json.dump(template, f, indent=2)
    
    print(f"\n✓ Saved benchmark template to {template_file}")

if __name__ == '__main__':
    main()
