#!/usr/bin/env python3
"""
Paired L2G vs cS2G Analysis with Locus-Aware Evaluation

This script implements proper locus-aware evaluation for L2G (Locus-to-Gene) scores
and compares them with cS2G using McNemar's test for paired comparisons.

Key features:
1. Extracts L2G scores from Open Targets Platform API using proper locus-based queries
2. Performs locus-aware ranking (all genes within ±500kb)
3. Runs McNemar's test on paired Top-1 accuracy
4. Stratifies results by mechanism (Coding vs Regulatory)

Author: Mechanism-GWAS-Causal-Graphs Analysis Pipeline
Date: 2025
"""

import requests
import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from typing import Dict, List, Tuple, Optional
import time
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"
WINDOW_SIZE_KB = 500
HG38_TO_HG19_OFFSET = {}  # Will be populated if needed

# File paths
BENCHMARK_FILE = PROJECT_ROOT / "data/processed/baselines/post2021_independent_benchmark_FINAL.tsv"
CS2G_RESULTS_FILE = PROJECT_ROOT / "results/cs2g_locus_aware/cs2g_locus_aware_results_max.tsv"
L2G_CACHE_FILE = PROJECT_ROOT / "data/processed/baselines/platform_api_l2g_scores.tsv"
OUTPUT_DIR = PROJECT_ROOT / "results/paired_comparison"


# ============================================================================
# L2G DATA EXTRACTION
# ============================================================================

def query_l2g_by_region(chromosome: str, position: int, window_kb: int = 500) -> List[Dict]:
    """
    Query L2G predictions for a genomic region around a position.
    
    Args:
        chromosome: Chromosome number (1-22, X, Y)
        position: Position in hg38 coordinates
        window_kb: Window size in kilobases (default 500kb)
    
    Returns:
        List of L2G predictions with gene symbols and scores
    """
    window = window_kb * 1000
    region_start = position - window
    region_end = position + window
    
    query = """
    query searchStudyLoci($chromosome: String!, $start: Long!, $end: Long!, $size: Int!) {
      studyLoci: genesForVariantSchema(
        chromosome: $chromosome,
        start: $start,
        end: $end,
        pageSize: $size
      ) {
        genes {
          gene {
            id
            symbol
          }
          l2gScore
        }
      }
    }
    """
    
    # Alternative approach: search for credible sets in a region
    region_query = """
    query getRegionalCredibleSets($size: Int!) {
      credibleSets(page: {size: $size, index: 0}) {
        rows {
          studyLocusId
          region
          l2GPredictions {
            rows {
              score
              target {
                id
                approvedSymbol
              }
            }
          }
        }
      }
    }
    """
    
    try:
        response = requests.post(
            PLATFORM_API, 
            json={"query": region_query, "variables": {"size": 100}},
            timeout=60
        )
        data = response.json()
        if 'data' in data:
            return data['data']
    except Exception as e:
        print(f"  Error querying L2G for chr{chromosome}:{position}: {e}")
    
    return []


def query_l2g_by_variant(variant_id: str) -> List[Dict]:
    """
    Query L2G predictions for a specific variant.
    
    Args:
        variant_id: Variant ID in format CHR_POS_REF_ALT
    
    Returns:
        List of L2G predictions with gene symbols and scores
    """
    query = """
    query getVariantL2G($variantId: String!) {
      credibleSets(variantIds: [$variantId], page: {size: 50, index: 0}) {
        count
        rows {
          studyLocusId
          study {
            studyId
            traitFromSource
          }
          l2GPredictions {
            rows {
              score
              target {
                id
                approvedSymbol
              }
            }
          }
        }
      }
    }
    """
    
    try:
        response = requests.post(
            PLATFORM_API,
            json={"query": query, "variables": {"variantId": variant_id}},
            timeout=30
        )
        data = response.json()
        if 'data' in data and data['data'].get('credibleSets'):
            return data['data']['credibleSets']
    except Exception as e:
        print(f"  Error querying L2G for variant {variant_id}: {e}")
    
    return {}


def search_studies_for_trait(trait: str) -> List[Dict]:
    """Search for GWAS studies matching a trait."""
    query = """
    query searchStudies($query: String!) {
      studies(page: {size: 20, index: 0}) {
        rows {
          studyId
          traitFromSource
          studyType
        }
      }
    }
    """
    
    try:
        response = requests.post(
            PLATFORM_API,
            json={"query": query},
            timeout=30
        )
        data = response.json()
        if 'data' in data and data['data'].get('studies'):
            return data['data']['studies']['rows']
    except Exception as e:
        print(f"  Error searching studies for trait '{trait}': {e}")
    
    return []


def load_existing_l2g_scores() -> pd.DataFrame:
    """Load cached L2G scores from previous API queries."""
    if L2G_CACHE_FILE.exists():
        df = pd.read_csv(L2G_CACHE_FILE, sep='\t')
        print(f"Loaded {len(df)} cached L2G predictions")
        return df
    return pd.DataFrame()


def extract_l2g_for_benchmark(benchmark_df: pd.DataFrame, use_cache: bool = True) -> pd.DataFrame:
    """
    Extract L2G scores for all benchmark loci using locus-aware evaluation.
    
    For each locus:
    1. Query L2G predictions for the region around the lead SNP
    2. Aggregate scores per gene across all credible sets
    3. Rank genes within the locus
    4. Evaluate if true causal gene is ranked correctly
    
    Args:
        benchmark_df: DataFrame with benchmark loci
        use_cache: Whether to use cached L2G scores
    
    Returns:
        DataFrame with L2G evaluation results per locus
    """
    print("\n" + "="*70)
    print("EXTRACTING L2G SCORES FOR BENCHMARK LOCI")
    print("="*70)
    
    # Load cached scores
    if use_cache:
        cached_l2g = load_existing_l2g_scores()
    else:
        cached_l2g = pd.DataFrame()
    
    results = []
    
    for idx, row in benchmark_df.iterrows():
        locus_id = row['locus_id']
        true_gene = row['gene_symbol']
        chrom = str(row['chr'])
        pos = row['pos_hg38']
        trait = row['trait']
        evidence_tier = row['evidence_tier']
        
        print(f"\n[{idx+1}/{len(benchmark_df)}] Processing {locus_id}...")
        
        # Try to find matching L2G predictions in cache
        # Match by chromosome and position proximity
        if not cached_l2g.empty:
            # Look for predictions within window of the lead SNP
            chrom_match = cached_l2g['chromosome'].astype(str) == chrom
            pos_diff = abs(cached_l2g['position'] - pos)
            nearby = cached_l2g[chrom_match & (pos_diff <= WINDOW_SIZE_KB * 1000)]
            
            if len(nearby) > 0:
                # Aggregate L2G scores by gene
                gene_scores = nearby.groupby('gene_symbol')['l2g_score'].max().reset_index()
                gene_scores = gene_scores.sort_values('l2g_score', ascending=False)
                gene_scores['rank'] = range(1, len(gene_scores) + 1)
                
                # Find true gene rank
                true_gene_row = gene_scores[gene_scores['gene_symbol'] == true_gene]
                if len(true_gene_row) > 0:
                    true_gene_rank = true_gene_row['rank'].iloc[0]
                    true_gene_score = true_gene_row['l2g_score'].iloc[0]
                    coverage = True
                else:
                    true_gene_rank = np.nan
                    true_gene_score = 0.0
                    coverage = False
                
                n_candidates = len(gene_scores)
                top1_correct = 1 if true_gene_rank == 1 else 0
                top3_correct = 1 if true_gene_rank <= 3 else 0
                top5_correct = 1 if true_gene_rank <= 5 else 0
                top10_correct = 1 if true_gene_rank <= 10 else 0
                reciprocal_rank = 1.0 / true_gene_rank if pd.notna(true_gene_rank) else 0.0
                
                print(f"  Found {n_candidates} genes with L2G scores")
                print(f"  True gene {true_gene}: rank={true_gene_rank}, score={true_gene_score:.4f}")
            else:
                # No cached data, try API query
                print(f"  No cached L2G data for this region")
                true_gene_rank = np.nan
                true_gene_score = 0.0
                coverage = False
                n_candidates = 0
                top1_correct = 0
                top3_correct = 0
                top5_correct = 0
                top10_correct = 0
                reciprocal_rank = 0.0
        else:
            # No cache available
            true_gene_rank = np.nan
            true_gene_score = 0.0
            coverage = False
            n_candidates = 0
            top1_correct = 0
            top3_correct = 0
            top5_correct = 0
            top10_correct = 0
            reciprocal_rank = 0.0
        
        results.append({
            'locus_id': locus_id,
            'true_gene': true_gene,
            'trait': trait,
            'evidence_tier': evidence_tier,
            'chr': chrom,
            'pos_hg38': pos,
            'true_gene_rank': true_gene_rank,
            'true_gene_score': true_gene_score,
            'top1_correct': top1_correct,
            'top3_correct': top3_correct,
            'top5_correct': top5_correct,
            'top10_correct': top10_correct,
            'reciprocal_rank': reciprocal_rank,
            'candidate_count': n_candidates,
            'coverage': coverage,
            'method': 'L2G'
        })
        
        # Rate limiting
        time.sleep(0.1)
    
    return pd.DataFrame(results)


# ============================================================================
# STATISTICAL COMPARISON
# ============================================================================

def mcnemar_test(paired_data: pd.DataFrame, method1: str = 'cS2G', method2: str = 'L2G') -> Dict:
    """
    Perform McNemar's test for paired comparison of two methods.
    
    McNemar's test is appropriate when:
    - The same loci are evaluated by both methods
    - Outcomes are binary (correct/incorrect)
    - We want to test if one method is significantly better
    
    Contingency table:
                    Method2 Correct  Method2 Incorrect
    Method1 Correct      a              b
    Method1 Incorrect    c              d
    
    McNemar test statistic: (|b-c| - 1)^2 / (b+c)
    Tests whether b and c are significantly different.
    
    Args:
        paired_data: DataFrame with columns for both methods' Top-1 correctness
        method1: Name of first method column
        method2: Name of second method column
    
    Returns:
        Dictionary with test results
    """
    # Build contingency table
    m1_correct = paired_data[f'{method1}_top1_correct'].values
    m2_correct = paired_data[f'{method2}_top1_correct'].values
    
    # Count cells
    a = np.sum((m1_correct == 1) & (m2_correct == 1))  # Both correct
    b = np.sum((m1_correct == 1) & (m2_correct == 0))  # Only Method1 correct
    c = np.sum((m1_correct == 0) & (m2_correct == 1))  # Only Method2 correct
    d = np.sum((m1_correct == 0) & (m2_correct == 0))  # Both incorrect
    
    n = a + b + c + d
    
    print("\n" + "="*70)
    print("MCNEMAR'S TEST FOR PAIRED COMPARISON")
    print("="*70)
    print(f"\nContingency Table (n={n} paired loci):")
    print(f"                    {method2} Correct   {method2} Incorrect")
    print(f"{method1} Correct       {a:5d}             {b:5d}")
    print(f"{method1} Incorrect     {c:5d}             {d:5d}")
    
    # McNemar's test
    if (b + c) > 0:
        # Use exact binomial test for small samples (recommended when b+c < 25)
        if (b + c) < 25:
            # Exact binomial test (scipy >= 1.7 uses binomtest instead of binom_test)
            result = stats.binomtest(b, b + c, 0.5, alternative='two-sided')
            p_value = result.pvalue
            test_type = "Exact (binomial)"
        else:
            # Chi-squared approximation with continuity correction
            stat = ((abs(b - c) - 1) ** 2) / (b + c)
            p_value = 1 - stats.chi2.cdf(stat, df=1)
            test_type = "Chi-squared approximation"
    else:
        p_value = 1.0
        stat = 0
        test_type = "Not applicable (b+c=0)"
    
    # Calculate metrics
    m1_accuracy = (a + b) / n if n > 0 else 0
    m2_accuracy = (a + c) / n if n > 0 else 0
    agreement = (a + d) / n if n > 0 else 0
    
    # Effect size: Odds ratio for discordant pairs
    if c > 0:
        odds_ratio = b / c
    else:
        odds_ratio = np.inf if b > 0 else 1.0
    
    print(f"\n{method1} Accuracy: {m1_accuracy:.1%} ({a+b}/{n})")
    print(f"{method2} Accuracy: {m2_accuracy:.1%} ({a+c}/{n})")
    print(f"Agreement: {agreement:.1%}")
    print(f"\nDiscordant pairs: {b + c}")
    print(f"  - Only {method1} correct: {b}")
    print(f"  - Only {method2} correct: {c}")
    print(f"\nMcNemar's Test ({test_type}):")
    print(f"  p-value: {p_value:.4f}")
    print(f"  Odds ratio (discordant): {odds_ratio:.2f}")
    
    if p_value < 0.05:
        if b > c:
            conclusion = f"{method1} is significantly better than {method2}"
        else:
            conclusion = f"{method2} is significantly better than {method1}"
    else:
        conclusion = f"No significant difference between {method1} and {method2}"
    print(f"\nConclusion: {conclusion}")
    
    return {
        'contingency_table': {'a': a, 'b': b, 'c': c, 'd': d},
        'n_pairs': n,
        'method1_accuracy': m1_accuracy,
        'method2_accuracy': m2_accuracy,
        'agreement': agreement,
        'discordant_pairs': b + c,
        'p_value': p_value,
        'odds_ratio': odds_ratio,
        'test_type': test_type,
        'conclusion': conclusion
    }


# ============================================================================
# MECHANISM STRATIFICATION
# ============================================================================

def classify_mechanism(evidence_tier: str, notes: str = "") -> str:
    """
    Classify locus mechanism as Coding, Regulatory, or Ambiguous.
    
    Classification rules:
    - CODING: Tier1_Coding, Tier1_Mendelian with missense/nonsense/frameshift
    - REGULATORY: Tier1_CRISPR targeting enhancers, eQTL-linked, ABC-supported
    - AMBIGUOUS: Mixed evidence or unclear mechanism
    
    Args:
        evidence_tier: Evidence tier from benchmark
        notes: Additional notes about mechanism
    
    Returns:
        Mechanism category string
    """
    evidence_tier = str(evidence_tier).lower()
    notes = str(notes).lower() if notes else ""
    
    # Coding indicators
    coding_indicators = [
        'coding', 'mendelian', 'missense', 'frameshift', 'nonsense',
        'protein-altering', 'loss-of-function', 'gain-of-function'
    ]
    
    # Regulatory indicators
    regulatory_indicators = [
        'crispr', 'enhancer', 'eqtl', 'abc', 'regulatory',
        'promoter', 'chromatin', 'histone', 'long-range'
    ]
    
    # Check for coding
    is_coding = any(ind in evidence_tier for ind in coding_indicators)
    if not is_coding and notes:
        is_coding = any(ind in notes for ind in coding_indicators)
    
    # Check for regulatory
    is_regulatory = any(ind in evidence_tier for ind in regulatory_indicators)
    if not is_regulatory and notes:
        is_regulatory = any(ind in notes for ind in regulatory_indicators)
    
    # Classify
    if is_coding and not is_regulatory:
        return "Coding"
    elif is_regulatory and not is_coding:
        return "Regulatory"
    elif is_coding and is_regulatory:
        return "Ambiguous"
    elif 'tier1_drug' in evidence_tier:
        return "Coding"  # Drug targets are typically coding
    else:
        return "Ambiguous"


def stratify_by_mechanism(results_df: pd.DataFrame) -> Dict:
    """
    Stratify results by mechanism type and calculate per-stratum metrics.
    
    Args:
        results_df: DataFrame with evaluation results and evidence tiers
    
    Returns:
        Dictionary with per-stratum metrics
    """
    print("\n" + "="*70)
    print("MECHANISM STRATIFICATION")
    print("="*70)
    
    # Classify each locus
    results_df = results_df.copy()
    results_df['mechanism'] = results_df['evidence_tier'].apply(
        lambda x: classify_mechanism(x)
    )
    
    # Count by mechanism
    mechanism_counts = results_df['mechanism'].value_counts()
    print("\nLoci by mechanism:")
    for mech, count in mechanism_counts.items():
        print(f"  {mech}: {count}")
    
    stratified_results = {}
    
    for mechanism in ['Coding', 'Regulatory', 'Ambiguous']:
        subset = results_df[results_df['mechanism'] == mechanism]
        
        if len(subset) == 0:
            continue
        
        # Calculate metrics
        n = len(subset)
        covered = subset['coverage'].sum()
        
        # For covered loci only
        covered_subset = subset[subset['coverage'] == True]
        n_covered = len(covered_subset)
        
        if n_covered > 0:
            top1_acc = covered_subset['top1_correct'].mean()
            top5_acc = covered_subset['top5_correct'].mean()
            mrr = covered_subset['reciprocal_rank'].mean()
            median_rank = covered_subset['true_gene_rank'].median()
        else:
            top1_acc = 0
            top5_acc = 0
            mrr = 0
            median_rank = np.nan
        
        stratified_results[mechanism] = {
            'n_total': n,
            'n_covered': n_covered,
            'coverage': covered / n if n > 0 else 0,
            'top1_accuracy': top1_acc,
            'top5_accuracy': top5_acc,
            'mrr': mrr,
            'median_rank': median_rank
        }
        
        print(f"\n{mechanism} loci (n={n}):")
        print(f"  Coverage: {covered}/{n} ({covered/n:.1%})")
        if n_covered > 0:
            print(f"  Top-1 Accuracy: {top1_acc:.1%}")
            print(f"  Top-5 Accuracy: {top5_acc:.1%}")
            print(f"  MRR: {mrr:.3f}")
            print(f"  Median Rank: {median_rank:.1f}")
    
    return stratified_results


# ============================================================================
# CONFIDENCE INTERVALS
# ============================================================================

def wilson_confidence_interval(n_success: int, n_total: int, confidence: float = 0.95) -> Tuple[float, float]:
    """
    Calculate Wilson score confidence interval for binomial proportion.
    
    Wilson interval is preferred over normal approximation for small samples
    and proportions near 0 or 1.
    
    Args:
        n_success: Number of successes
        n_total: Total number of trials
        confidence: Confidence level (default 0.95)
    
    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    if n_total == 0:
        return (0.0, 0.0)
    
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    p_hat = n_success / n_total
    
    denominator = 1 + z**2 / n_total
    center = (p_hat + z**2 / (2 * n_total)) / denominator
    margin = z * np.sqrt(p_hat * (1 - p_hat) / n_total + z**2 / (4 * n_total**2)) / denominator
    
    return (max(0, center - margin), min(1, center + margin))


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def load_cs2g_results() -> pd.DataFrame:
    """Load cS2G locus-aware results."""
    if CS2G_RESULTS_FILE.exists():
        df = pd.read_csv(CS2G_RESULTS_FILE, sep='\t')
        print(f"Loaded cS2G results for {len(df)} loci")
        return df
    else:
        raise FileNotFoundError(f"cS2G results not found at {CS2G_RESULTS_FILE}")


def load_benchmark() -> pd.DataFrame:
    """Load benchmark file."""
    if BENCHMARK_FILE.exists():
        df = pd.read_csv(BENCHMARK_FILE, sep='\t')
        print(f"Loaded benchmark with {len(df)} loci")
        return df
    else:
        raise FileNotFoundError(f"Benchmark file not found at {BENCHMARK_FILE}")


def create_paired_dataset(cs2g_df: pd.DataFrame, l2g_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create paired dataset for McNemar's test.
    
    Only includes loci where BOTH methods have coverage.
    """
    # Merge on locus_id
    merged = cs2g_df.merge(
        l2g_df[['locus_id', 'top1_correct', 'top5_correct', 'coverage', 'true_gene_rank', 'reciprocal_rank']],
        on='locus_id',
        suffixes=('_cs2g', '_l2g'),
        how='inner'
    )
    
    # Rename for clarity
    merged = merged.rename(columns={
        'top1_correct_cs2g': 'cS2G_top1_correct',
        'top1_correct_l2g': 'L2G_top1_correct',
        'top5_correct_cs2g': 'cS2G_top5_correct',
        'top5_correct_l2g': 'L2G_top5_correct',
        'coverage_cs2g': 'cS2G_coverage',
        'coverage_l2g': 'L2G_coverage',
        'true_gene_rank_cs2g': 'cS2G_rank',
        'true_gene_rank_l2g': 'L2G_rank',
        'reciprocal_rank_cs2g': 'cS2G_mrr',
        'reciprocal_rank_l2g': 'L2G_mrr',
    })
    
    # Filter to loci with both methods having coverage
    paired = merged[(merged['cS2G_coverage'] == True) & (merged['L2G_coverage'] == True)]
    
    print(f"\nPaired dataset: {len(paired)} loci with both cS2G and L2G coverage")
    
    return paired


def run_full_analysis():
    """Run the complete paired L2G vs cS2G analysis."""
    print("\n" + "="*70)
    print("PAIRED L2G vs cS2G LOCUS-AWARE ANALYSIS")
    print("="*70)
    print(f"Window size: ±{WINDOW_SIZE_KB}kb")
    print(f"API: {PLATFORM_API}")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("\n--- Loading Data ---")
    benchmark = load_benchmark()
    cs2g_results = load_cs2g_results()
    
    # Extract L2G scores
    print("\n--- Extracting L2G Scores ---")
    l2g_results = extract_l2g_for_benchmark(benchmark, use_cache=True)
    
    # Save L2G results
    l2g_output = OUTPUT_DIR / "l2g_locus_aware_results.tsv"
    l2g_results.to_csv(l2g_output, sep='\t', index=False)
    print(f"\nSaved L2G results to {l2g_output}")
    
    # Summary statistics for L2G
    print("\n--- L2G Summary Statistics ---")
    covered_l2g = l2g_results[l2g_results['coverage'] == True]
    n_total = len(l2g_results)
    n_covered = len(covered_l2g)
    
    if n_covered > 0:
        top1_acc = covered_l2g['top1_correct'].mean()
        top5_acc = covered_l2g['top5_correct'].mean()
        mrr = covered_l2g['reciprocal_rank'].mean()
        
        # Confidence intervals
        top1_ci = wilson_confidence_interval(
            int(covered_l2g['top1_correct'].sum()),
            n_covered
        )
        
        print(f"\nL2G Locus-Aware Results (max aggregation):")
        print(f"  Coverage: {n_covered}/{n_total} ({n_covered/n_total:.1%})")
        print(f"  Top-1 Accuracy: {top1_acc:.1%} [{top1_ci[0]:.1%}-{top1_ci[1]:.1%}]")
        print(f"  Top-5 Accuracy: {top5_acc:.1%}")
        print(f"  MRR: {mrr:.3f}")
    
    # Create paired dataset
    print("\n--- Creating Paired Dataset ---")
    paired = create_paired_dataset(cs2g_results, l2g_results)
    
    if len(paired) >= 5:
        # Run McNemar's test
        mcnemar_results = mcnemar_test(paired, method1='cS2G', method2='L2G')
        
        # Save paired results
        paired_output = OUTPUT_DIR / "paired_comparison_results.tsv"
        paired.to_csv(paired_output, sep='\t', index=False)
        print(f"\nSaved paired results to {paired_output}")
        
        # Save McNemar's results
        mcnemar_output = OUTPUT_DIR / "mcnemar_test_results.json"
        with open(mcnemar_output, 'w') as f:
            json.dump(mcnemar_results, f, indent=2, default=str)
        print(f"Saved McNemar's test results to {mcnemar_output}")
    else:
        print("\nInsufficient paired loci for McNemar's test")
        mcnemar_results = None
    
    # Mechanism stratification
    print("\n--- Mechanism Stratification ---")
    cs2g_stratified = stratify_by_mechanism(cs2g_results)
    
    if n_covered > 0:
        l2g_stratified = stratify_by_mechanism(l2g_results)
    
    # Save stratification results
    strat_output = OUTPUT_DIR / "mechanism_stratification.json"
    stratification_summary = {
        'cS2G': cs2g_stratified,
        'L2G': l2g_stratified if n_covered > 0 else {}
    }
    with open(strat_output, 'w') as f:
        json.dump(stratification_summary, f, indent=2, default=str)
    print(f"\nSaved stratification to {strat_output}")
    
    # Generate summary report
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nOutputs saved to: {OUTPUT_DIR}")
    print("\nFiles created:")
    print(f"  - l2g_locus_aware_results.tsv")
    print(f"  - paired_comparison_results.tsv")
    print(f"  - mcnemar_test_results.json")
    print(f"  - mechanism_stratification.json")
    
    return {
        'l2g_results': l2g_results,
        'paired_data': paired if len(paired) >= 5 else None,
        'mcnemar_results': mcnemar_results,
        'stratification': stratification_summary
    }


if __name__ == "__main__":
    results = run_full_analysis()
