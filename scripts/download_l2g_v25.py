#!/usr/bin/env python3
"""
Download Open Targets Platform L2G v25.12 predictions via GraphQL API.

This script queries the current Open Targets Platform (v25.12) for L2G predictions
corresponding to our benchmark variants, allowing us to compare calibration between
L2G v22.09 (bulk download) and L2G v25.12 (Platform API).
"""

import requests
import pandas as pd
import json
from pathlib import Path
from tqdm import tqdm
import time

# Open Targets Platform GraphQL API
PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

def query_credible_set_l2g(study_locus_id: str) -> dict:
    """
    Query L2G predictions for a specific credible set (study-locus) from Platform API.
    
    Args:
        study_locus_id: StudyLocus ID in format like "GCST004131_1_205723522"
        
    Returns:
        Dict with L2G predictions including scores and features
    """
    query = """
    query getL2GPredictions($studyLocusId: String!) {
      credibleSet(studyLocusId: $studyLocusId) {
        studyLocusId
        studyId
        chromosome
        position
        variant {
          id
          rsIds
        }
        l2GPredictions {
          count
          rows {
            score
            target {
              id
              approvedSymbol
            }
            features {
              name
              value
              shapValue
            }
            shapBaseValue
          }
        }
      }
    }
    """
    
    variables = {"studyLocusId": study_locus_id}
    
    try:
        response = requests.post(
            PLATFORM_API,
            json={"query": query, "variables": variables},
            timeout=30
        )
        response.raise_for_status()
        data = response.json()
        
        if "errors" in data:
            print(f"GraphQL errors for {study_locus_id}: {data['errors']}")
            return None
            
        return data.get("data", {}).get("credibleSet")
        
    except Exception as e:
        print(f"Error querying {study_locus_id}: {e}")
        return None


def load_benchmark_predictions():
    """Load our v22.09 benchmark L2G predictions with ground truth labels."""
    
    # Path to our processed L2G data
    l2g_path = Path("data/external/opentargets_l2g/processed/l2g_processed.parquet")
    
    if not l2g_path.exists():
        raise FileNotFoundError(f"L2G data not found at {l2g_path}")
    
    print(f"Loading v22.09 benchmark from {l2g_path}")
    df = pd.read_parquet(l2g_path)
    
    print(f"Loaded {len(df):,} L2G predictions")
    print(f"Columns: {df.columns.tolist()}")
    print(f"\nFirst few rows:")
    print(df.head())
    
    return df


def construct_study_locus_id(row) -> str:
    """
    Construct study-locus ID from benchmark data.
    
    Platform StudyLocus IDs follow format: STUDYID_CHR_POS
    Example: GCST004131_1_205723522
    
    Args:
        row: DataFrame row with study_id, chr, pos
        
    Returns:
        Study-locus ID string
    """
    # Need to determine how to construct this from our data
    # Our v22.09 data has: study_id, variant_id, gene_id, y2g_score
    # Platform expects: studyLocusId in format STUDYID_CHR_POS
    
    # For now, extract chr and pos from variant_id (format: CHR_POS_REF_ALT)
    variant_parts = row['variant_id'].split('_')
    if len(variant_parts) >= 2:
        chr_num = variant_parts[0]
        pos = variant_parts[1]
        return f"{row['study_id']}_{chr_num}_{pos}"
    else:
        return None


def download_platform_l2g():
    """
    Download L2G v25.12 predictions from Platform API for our benchmark variants.
    
    Returns:
        DataFrame with Platform L2G predictions matched to our benchmark
    """
    
    # Load our v22.09 benchmark
    benchmark_df = load_benchmark_predictions()
    
    # Get unique study-locus combinations to query
    print("\nConstructing study-locus IDs...")
    benchmark_df['study_locus_id'] = benchmark_df.apply(construct_study_locus_id, axis=1)
    
    # Remove rows where study_locus_id construction failed
    valid_df = benchmark_df[benchmark_df['study_locus_id'].notna()].copy()
    print(f"Valid study-locus IDs: {len(valid_df):,}")
    
    # Get unique study-locus IDs (one per credible set)
    unique_loci = valid_df['study_locus_id'].unique()
    print(f"Unique credible sets to query: {len(unique_loci):,}")
    
    # Query Platform API for each credible set
    print(f"\nQuerying Platform API for L2G v25.12 predictions...")
    platform_predictions = []
    
    for study_locus_id in tqdm(unique_loci[:100]):  # Start with first 100 for testing
        result = query_credible_set_l2g(study_locus_id)
        
        if result and result.get('l2GPredictions'):
            l2g_preds = result['l2GPredictions']
            
            for pred in l2g_preds.get('rows', []):
                platform_predictions.append({
                    'study_locus_id': study_locus_id,
                    'study_id': result['studyId'],
                    'chromosome': result['chromosome'],
                    'position': result['position'],
                    'variant_id': result['variant']['id'] if result.get('variant') else None,
                    'gene_id': pred['target']['id'],
                    'gene_symbol': pred['target']['approvedSymbol'],
                    'l2g_score_v25': pred['score'],
                    'shap_base_value': pred.get('shapBaseValue'),
                    'n_features': len(pred.get('features', []))
                })
        
        # Rate limiting - be nice to the API
        time.sleep(0.1)
    
    # Convert to DataFrame
    platform_df = pd.DataFrame(platform_predictions)
    
    print(f"\nDownloaded {len(platform_df):,} L2G v25.12 predictions")
    
    # Save results
    output_path = Path("data/external/opentargets_l2g/processed/l2g_v25_platform.parquet")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    platform_df.to_parquet(output_path)
    print(f"Saved to {output_path}")
    
    return platform_df


def merge_with_benchmark():
    """
    Merge Platform v25.12 predictions with v22.09 benchmark for comparison.
    
    Returns:
        DataFrame with both v22.09 and v25.12 scores for the same gene-locus pairs
    """
    
    # Load v22.09 benchmark
    benchmark_df = load_benchmark_predictions()
    
    # Load v25.12 Platform predictions
    platform_path = Path("data/external/opentargets_l2g/processed/l2g_v25_platform.parquet")
    
    if not platform_path.exists():
        print("Platform v25.12 data not found. Run download_platform_l2g() first.")
        return None
    
    platform_df = pd.read_parquet(platform_path)
    
    # Merge on study_locus_id + gene_id
    print(f"\nMerging v22.09 and v25.12 predictions...")
    
    # Add study_locus_id to benchmark if not present
    if 'study_locus_id' not in benchmark_df.columns:
        benchmark_df['study_locus_id'] = benchmark_df.apply(construct_study_locus_id, axis=1)
    
    # Merge
    merged_df = benchmark_df.merge(
        platform_df[['study_locus_id', 'gene_id', 'l2g_score_v25', 'shap_base_value']],
        on=['study_locus_id', 'gene_id'],
        how='inner',
        suffixes=('_v22', '_v25')
    )
    
    print(f"Matched {len(merged_df):,} gene-locus pairs between v22.09 and v25.12")
    
    # Calculate score differences
    merged_df['score_diff'] = merged_df['l2g_score_v25'] - merged_df['y2g_score']
    merged_df['score_ratio'] = merged_df['l2g_score_v25'] / merged_df['y2g_score'].replace(0, 1e-10)
    
    # Save merged comparison
    output_path = Path("data/external/opentargets_l2g/processed/l2g_v22_v25_comparison.parquet")
    merged_df.to_parquet(output_path)
    print(f"Saved comparison to {output_path}")
    
    return merged_df


if __name__ == "__main__":
    print("=" * 80)
    print("Downloading Open Targets Platform L2G v25.12 predictions")
    print("=" * 80)
    
    # Step 1: Download Platform v25.12 predictions
    platform_df = download_platform_l2g()
    
    # Step 2: Merge with v22.09 benchmark
    comparison_df = merge_with_benchmark()
    
    if comparison_df is not None:
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        print(f"Total matched predictions: {len(comparison_df):,}")
        print(f"v22.09 score range: {comparison_df['y2g_score'].min():.4f} - {comparison_df['y2g_score'].max():.4f}")
        print(f"v25.12 score range: {comparison_df['l2g_score_v25'].min():.4f} - {comparison_df['l2g_score_v25'].max():.4f}")
        print(f"Mean score difference: {comparison_df['score_diff'].mean():.4f}")
        print(f"Median score difference: {comparison_df['score_diff'].median():.4f}")
        print(f"Correlation: {comparison_df[['y2g_score', 'l2g_score_v25']].corr().iloc[0, 1]:.4f}")
