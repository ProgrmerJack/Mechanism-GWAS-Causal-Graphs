#!/usr/bin/env python3
"""
Build Expanded Regulatory Benchmark from E2G Benchmarking Data

This script processes the EngreitzLab E2G benchmarking dataset to create
a robust regulatory evaluation set for L2G and cS2G baselines.

The E2G benchmark contains silver-standard variant-gene links based on:
1. ABC (Activity-by-Contact) enhancer-gene predictions
2. POPS (Polygenic Priority Score) rankings  
3. Various regulatory evidence types (eQTL, CRISPR, PCHiC, etc.)

Key columns in UKBiobank.ABCGene.anyabc.tsv:
- CredibleSet: chr:start-end-id format
- Disease: Trait abbreviation
- TargetGene: Gene symbol
- truth: Whether this is a validated link (TRUE/FALSE)
- MaxABC: Maximum ABC score across cell types
- ABCPrediction: Binary ABC prediction
- Various method predictions

We filter for high-confidence regulatory loci where:
1. Truth label is available
2. Non-coding mechanism (not coding/splice variants)
3. Has ABC or other regulatory predictions

Author: Analysis Pipeline
Date: 2025-01-19
"""

import pandas as pd
import numpy as np
import requests
import json
import time
from pathlib import Path
from collections import defaultdict

# Paths
DATA_DIR = Path(__file__).parent.parent / "data"
EXTERNAL_DIR = DATA_DIR / "external" / "E2G_benchmarking" / "resources"
OUTPUT_DIR = DATA_DIR / "processed" / "baselines"

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"


def parse_credible_set_position(cs_str: str) -> dict:
    """
    Parse credible set string to extract chromosome and position.
    Format: chr11:117072267-120072267-2
    Returns: {'chrom': '11', 'start': 117072267, 'end': 120072267, 'id': '2'}
    """
    try:
        chrom_part, rest = cs_str.split(':')
        chrom = chrom_part.replace('chr', '')
        parts = rest.split('-')
        start = int(parts[0])
        end = int(parts[1])
        cs_id = parts[2] if len(parts) > 2 else '1'
        return {
            'chrom': chrom,
            'start': start,
            'end': end,
            'center': (start + end) // 2,
            'id': cs_id
        }
    except Exception as e:
        print(f"Error parsing credible set: {cs_str}: {e}")
        return None


def load_e2g_benchmark() -> pd.DataFrame:
    """Load and parse the E2G benchmarking data."""
    filepath = EXTERNAL_DIR / "UKBiobank.ABCGene.anyabc.tsv"
    print(f"Loading E2G benchmark from {filepath}")
    
    df = pd.read_csv(filepath, sep='\t')
    print(f"Loaded {len(df)} rows with {len(df.columns)} columns")
    
    # Parse credible set positions
    positions = []
    for cs in df['CredibleSet']:
        pos = parse_credible_set_position(cs)
        positions.append(pos)
    
    df['parsed_pos'] = positions
    df['chrom'] = df['parsed_pos'].apply(lambda x: x['chrom'] if x else None)
    df['center_pos'] = df['parsed_pos'].apply(lambda x: x['center'] if x else None)
    
    return df


def filter_regulatory_loci(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for high-confidence regulatory (non-coding) loci.
    
    Criteria:
    1. Has truth label
    2. Not coding/splice variants (CodingSpliceOrPromoterVariants == FALSE)
    3. Has some regulatory evidence (ABC or other)
    """
    # Filter for non-missing truth labels
    df_truth = df[df['truth'].notna()].copy()
    print(f"  Records with truth labels: {len(df_truth)}")
    
    # Filter for non-coding (excluding coding variants)
    df_noncoding = df_truth[df_truth['CodingSpliceOrPromoterVariants'] == False].copy()
    print(f"  Non-coding variants: {len(df_noncoding)}")
    
    # Get true positives (validated regulatory links)
    df_true = df_noncoding[df_noncoding['truth'] == True].copy()
    print(f"  True regulatory links: {len(df_true)}")
    
    return df_noncoding, df_true


def get_unique_loci(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get unique locus-gene pairs with their truth labels.
    Each credible set can have multiple genes tested.
    """
    # Create unique locus ID
    df['locus_id'] = df['CredibleSet'] + '_' + df['TargetGene']
    
    # Group by credible set to understand locus structure
    loci = df.groupby('CredibleSet').agg({
        'Disease': 'first',
        'TargetGene': lambda x: list(x.unique()),
        'truth': lambda x: list(x),
        'chrom': 'first',
        'center_pos': 'first',
        'MaxABC': 'max',
        'ABCPrediction': 'any'
    }).reset_index()
    
    print(f"  Unique credible sets (loci): {len(loci)}")
    
    return loci


def query_l2g_for_locus(chrom: str, pos: int, window: int = 500000) -> list:
    """
    Query L2G scores for a genomic region using Platform API.
    Uses credibleSets endpoint to find overlapping fine-mapping data.
    """
    query = """
    query credibleSetsByRegion($chrom: String!, $start: Long!, $end: Long!) {
      credibleSets(
        chromosome: $chrom
        start: $start
        end: $end
      ) {
        rows {
          studyLocusId
          variant {
            id
            chromosome
            position
          }
          l2GPredictions {
            rows {
              target {
                id
                approvedSymbol
              }
              score
            }
          }
        }
      }
    }
    """
    
    variables = {
        'chrom': str(chrom),
        'start': max(0, pos - window),
        'end': pos + window
    }
    
    try:
        response = requests.post(
            PLATFORM_API,
            json={'query': query, 'variables': variables},
            headers={'Content-Type': 'application/json'},
            timeout=30
        )
        
        if response.status_code == 200:
            data = response.json()
            if 'data' in data and data['data'] and 'credibleSets' in data['data']:
                return data['data']['credibleSets']['rows']
        
        return []
    except Exception as e:
        print(f"  API error for {chrom}:{pos}: {e}")
        return []


def extract_l2g_for_gene(credible_sets: list, target_gene: str) -> dict:
    """
    Extract L2G score for a specific target gene from credible set data.
    Returns the maximum L2G score for the gene across all overlapping credible sets.
    """
    max_score = None
    best_variant = None
    
    for cs in credible_sets:
        if not cs.get('l2GPredictions') or not cs['l2GPredictions'].get('rows'):
            continue
            
        for pred in cs['l2GPredictions']['rows']:
            gene_symbol = pred.get('target', {}).get('approvedSymbol', '')
            if gene_symbol.upper() == target_gene.upper():
                score = pred.get('score', 0)
                if max_score is None or score > max_score:
                    max_score = score
                    if cs.get('variant'):
                        best_variant = cs['variant'].get('id')
    
    return {
        'l2g_score': max_score,
        'best_variant': best_variant
    }


def build_regulatory_benchmark(df: pd.DataFrame, sample_size: int = None) -> pd.DataFrame:
    """
    Build benchmark by querying L2G for each regulatory locus.
    
    Args:
        df: DataFrame with regulatory loci
        sample_size: If specified, sample this many loci for testing
    """
    # Get true positives only for validation
    df_true = df[df['truth'] == True].copy()
    
    if sample_size:
        df_true = df_true.sample(n=min(sample_size, len(df_true)), random_state=42)
    
    print(f"\nBuilding benchmark for {len(df_true)} true regulatory links...")
    
    results = []
    cache = {}  # Cache API calls by region
    
    for idx, row in df_true.iterrows():
        chrom = row['chrom']
        pos = row['center_pos']
        gene = row['TargetGene']
        disease = row['Disease']
        cs = row['CredibleSet']
        abc_score = row.get('MaxABC', 0)
        
        # Create cache key
        cache_key = f"{chrom}:{pos}"
        
        # Query L2G (with caching)
        if cache_key not in cache:
            credible_sets = query_l2g_for_locus(chrom, pos)
            cache[cache_key] = credible_sets
            time.sleep(0.2)  # Rate limiting
        else:
            credible_sets = cache[cache_key]
        
        # Extract L2G for target gene
        l2g_result = extract_l2g_for_gene(credible_sets, gene)
        
        results.append({
            'credible_set': cs,
            'chrom': chrom,
            'position': pos,
            'disease': disease,
            'target_gene': gene,
            'truth': True,
            'abc_score': abc_score,
            'l2g_score': l2g_result['l2g_score'],
            'l2g_variant': l2g_result['best_variant'],
            'n_credible_sets': len(credible_sets)
        })
        
        if len(results) % 50 == 0:
            print(f"  Processed {len(results)} loci...")
    
    return pd.DataFrame(results)


def evaluate_l2g_performance(benchmark_df: pd.DataFrame) -> dict:
    """
    Evaluate L2G performance on the regulatory benchmark.
    """
    # Filter to those with L2G scores
    df_scored = benchmark_df[benchmark_df['l2g_score'].notna()].copy()
    
    total = len(benchmark_df)
    covered = len(df_scored)
    
    # At various thresholds
    thresholds = [0.1, 0.3, 0.5, 0.7, 0.9]
    performance = {}
    
    for thresh in thresholds:
        correct = (df_scored['l2g_score'] >= thresh).sum()
        performance[f'correct_at_{thresh}'] = correct
        performance[f'accuracy_at_{thresh}'] = correct / covered if covered > 0 else 0
    
    performance['total_loci'] = total
    performance['l2g_coverage'] = covered
    performance['coverage_rate'] = covered / total if total > 0 else 0
    performance['mean_l2g'] = df_scored['l2g_score'].mean() if covered > 0 else None
    performance['median_l2g'] = df_scored['l2g_score'].median() if covered > 0 else None
    
    return performance


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("BUILDING EXPANDED REGULATORY BENCHMARK")
    print("=" * 70)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load E2G benchmark
    print("\n1. Loading E2G benchmarking data...")
    df = load_e2g_benchmark()
    
    # Filter for regulatory loci
    print("\n2. Filtering for regulatory (non-coding) loci...")
    df_noncoding, df_true = filter_regulatory_loci(df)
    
    # Get unique loci statistics
    print("\n3. Analyzing locus structure...")
    loci = get_unique_loci(df_noncoding)
    
    # Build benchmark - start with sample for testing
    print("\n4. Building regulatory benchmark (querying L2G)...")
    
    # First do a quick test with 20 loci
    print("\n   Testing with 20 loci first...")
    test_benchmark = build_regulatory_benchmark(df_noncoding, sample_size=20)
    
    # Save test results
    test_output = OUTPUT_DIR / "regulatory_benchmark_test.tsv"
    test_benchmark.to_csv(test_output, sep='\t', index=False)
    print(f"   Saved test benchmark to {test_output}")
    
    # Evaluate test performance
    print("\n5. Evaluating L2G performance on test set...")
    test_perf = evaluate_l2g_performance(test_benchmark)
    
    print(f"\n   TEST RESULTS:")
    print(f"   - Total loci: {test_perf['total_loci']}")
    print(f"   - L2G coverage: {test_perf['l2g_coverage']} ({test_perf['coverage_rate']:.1%})")
    if test_perf['mean_l2g']:
        print(f"   - Mean L2G score: {test_perf['mean_l2g']:.3f}")
        print(f"   - Accuracy at 0.5: {test_perf['accuracy_at_0.5']:.1%}")
    
    # If test looks good, build full benchmark
    # Note: Full run would take ~30 minutes for 200+ loci
    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    
    # Summary statistics
    print(f"\nE2G Dataset Statistics:")
    print(f"  - Total variant-gene pairs: {len(df)}")
    print(f"  - Non-coding pairs: {len(df_noncoding)}")
    print(f"  - True regulatory links: {len(df_true)}")
    print(f"  - Unique diseases: {df['Disease'].nunique()}")
    print(f"  - Unique credible sets: {df['CredibleSet'].nunique()}")
    
    # List diseases
    print(f"\nDisease categories in benchmark:")
    disease_counts = df_true.groupby('Disease').size().sort_values(ascending=False)
    for disease, count in disease_counts.head(20).items():
        print(f"    {disease}: {count} true links")
    
    return test_benchmark, test_perf


if __name__ == "__main__":
    benchmark, performance = main()
