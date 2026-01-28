#!/usr/bin/env python3
"""
Build Expanded Regulatory Benchmark with Variant-Level L2G Queries

This script:
1. Loads the E2G benchmarking data with true regulatory links
2. Maps credible sets to their lead variants from the SuSiE files
3. Queries Platform API for L2G scores using variant IDs
4. Creates a robust benchmark for regulatory mechanism evaluation

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
import gzip

# Paths
DATA_DIR = Path(__file__).parent.parent / "data"
EXTERNAL_DIR = DATA_DIR / "external" / "E2G_benchmarking" / "resources"
SUSIE_DIR = EXTERNAL_DIR / "191010_UKBB_SuSiE_hg38_liftover"
OUTPUT_DIR = DATA_DIR / "processed" / "baselines"

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"


def load_e2g_benchmark() -> pd.DataFrame:
    """Load the E2G benchmarking data."""
    filepath = EXTERNAL_DIR / "UKBiobank.ABCGene.anyabc.tsv"
    print(f"Loading E2G benchmark from {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    print(f"  Loaded {len(df)} rows")
    return df


def get_true_regulatory_links(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for true regulatory (non-coding) links.
    These are silver-standard validated enhancer-gene connections.
    """
    # Filter for non-coding true positives
    df_true = df[
        (df['truth'] == True) & 
        (df['CodingSpliceOrPromoterVariants'] == False)
    ].copy()
    
    print(f"  True non-coding regulatory links: {len(df_true)}")
    return df_true


def load_lead_variants_for_trait(trait: str) -> pd.DataFrame:
    """Load the SuSiE variant data for a given trait."""
    trait_dir = SUSIE_DIR / trait
    variant_file = trait_dir / "variant.list.txt"
    
    if not variant_file.exists():
        return None
    
    try:
        df = pd.read_csv(variant_file, sep='\t')
        return df
    except Exception as e:
        print(f"  Warning: Could not load {variant_file}: {e}")
        return None


def get_lead_variant_for_credible_set(variants_df: pd.DataFrame, cs_id: str) -> dict:
    """
    Get the lead (highest PIP) variant for a credible set.
    Returns variant info including chr_pos_ref_alt format for API.
    """
    if variants_df is None:
        return None
    
    # Filter to this credible set
    cs_variants = variants_df[variants_df['CredibleSet'] == cs_id].copy()
    
    if len(cs_variants) == 0:
        return None
    
    # Get variant with highest posterior probability
    lead_idx = cs_variants['PosteriorProb'].idxmax()
    lead = cs_variants.loc[lead_idx]
    
    # Build variant ID in Platform format: CHR_POS_REF_ALT
    chrom = str(lead['chromosome']).replace('chr', '')
    pos = int(lead['position'])
    ref = str(lead['allele1'])
    alt = str(lead['allele2'])
    rsid = str(lead['rsid'])
    pip = float(lead['PosteriorProb'])
    
    variant_id = f"{chrom}_{pos}_{ref}_{alt}"
    
    return {
        'variant_id': variant_id,
        'rsid': rsid,
        'chrom': chrom,
        'position': pos,
        'ref': ref,
        'alt': alt,
        'pip': pip
    }


def query_l2g_for_variant(variant_id: str) -> dict:
    """
    Query L2G scores for a variant using the Platform API.
    Returns dict with gene -> score mapping.
    """
    query = """
    query getVariantL2G($variantId: String!) {
      variant(variantId: $variantId) {
        id
        rsIds
        credibleSets {
          count
          rows {
            studyLocusId
            study {
              id
              traitFromSourceMappedIds
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
    }
    """
    
    try:
        response = requests.post(
            PLATFORM_API,
            json={'query': query, 'variables': {'variantId': variant_id}},
            headers={'Content-Type': 'application/json'},
            timeout=30
        )
        
        if response.status_code != 200:
            return {'error': f'HTTP {response.status_code}'}
        
        data = response.json()
        
        if 'errors' in data:
            return {'error': data['errors'][0]['message'][:100]}
        
        if not data.get('data') or not data['data'].get('variant'):
            return {'error': 'Variant not found'}
        
        variant = data['data']['variant']
        credible_sets = variant.get('credibleSets', {}).get('rows', [])
        
        # Aggregate L2G scores across all credible sets
        gene_scores = defaultdict(float)
        n_credible_sets = len(credible_sets)
        
        for cs in credible_sets:
            l2g_preds = cs.get('l2GPredictions', {}).get('rows', [])
            for pred in l2g_preds:
                gene = pred.get('target', {}).get('approvedSymbol', '')
                score = pred.get('score', 0)
                if gene and score > gene_scores[gene]:
                    gene_scores[gene] = score
        
        return {
            'gene_scores': dict(gene_scores),
            'n_credible_sets': n_credible_sets,
            'rsids': variant.get('rsIds', [])
        }
        
    except Exception as e:
        return {'error': str(e)[:100]}


def build_regulatory_benchmark(benchmark_df: pd.DataFrame, max_loci: int = None) -> pd.DataFrame:
    """
    Build the regulatory benchmark by querying L2G for each true link.
    """
    results = []
    
    # Cache variant data per trait
    trait_variants = {}
    
    # Track processing
    processed = 0
    found = 0
    matched = 0
    
    # Group by credible set to avoid duplicate queries
    grouped = benchmark_df.groupby('CredibleSet').first().reset_index()
    
    if max_loci:
        grouped = grouped.head(max_loci)
    
    print(f"\nProcessing {len(grouped)} unique credible sets...")
    
    for idx, row in grouped.iterrows():
        cs_id = row['CredibleSet']
        disease = row['Disease']
        
        # Load variant data for this trait if not cached
        if disease not in trait_variants:
            trait_variants[disease] = load_lead_variants_for_trait(disease)
        
        # Get lead variant
        lead_var = get_lead_variant_for_credible_set(trait_variants[disease], cs_id)
        
        if lead_var is None:
            continue
        
        found += 1
        
        # Query L2G
        l2g_result = query_l2g_for_variant(lead_var['variant_id'])
        time.sleep(0.15)  # Rate limiting
        
        # Get all true genes for this credible set
        cs_links = benchmark_df[benchmark_df['CredibleSet'] == cs_id]
        
        for _, link in cs_links.iterrows():
            target_gene = link['TargetGene']
            abc_score = link.get('MaxABC', 0)
            
            # Get L2G score for target gene
            l2g_score = None
            if 'gene_scores' in l2g_result:
                l2g_score = l2g_result['gene_scores'].get(target_gene)
                
                # Try case-insensitive match
                if l2g_score is None:
                    for gene, score in l2g_result['gene_scores'].items():
                        if gene.upper() == target_gene.upper():
                            l2g_score = score
                            break
            
            results.append({
                'credible_set': cs_id,
                'disease': disease,
                'target_gene': target_gene,
                'truth': True,
                'abc_score': abc_score,
                'lead_variant_id': lead_var['variant_id'],
                'lead_rsid': lead_var['rsid'],
                'lead_pip': lead_var['pip'],
                'l2g_score': l2g_score,
                'l2g_n_credible_sets': l2g_result.get('n_credible_sets', 0),
                'l2g_error': l2g_result.get('error')
            })
            
            if l2g_score is not None:
                matched += 1
        
        processed += 1
        if processed % 20 == 0:
            print(f"  Processed {processed}/{len(grouped)} credible sets...")
    
    print(f"\n  Total credible sets: {len(grouped)}")
    print(f"  Found lead variants: {found}")
    print(f"  Gene-L2G matches: {matched}")
    
    return pd.DataFrame(results)


def evaluate_benchmark(df: pd.DataFrame) -> dict:
    """Evaluate L2G performance on the benchmark."""
    # Filter to those with L2G scores
    df_scored = df[df['l2g_score'].notna()].copy()
    
    total = len(df)
    covered = len(df_scored)
    
    results = {
        'total_links': total,
        'l2g_coverage': covered,
        'coverage_rate': covered / total if total > 0 else 0,
    }
    
    if covered > 0:
        results['mean_l2g'] = df_scored['l2g_score'].mean()
        results['median_l2g'] = df_scored['l2g_score'].median()
        
        # Accuracy at thresholds
        for thresh in [0.1, 0.3, 0.5, 0.7]:
            correct = (df_scored['l2g_score'] >= thresh).sum()
            results[f'accuracy_at_{thresh}'] = correct / covered
    
    return results


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("BUILDING EXPANDED REGULATORY BENCHMARK V2")
    print("(Using variant-level L2G queries)")
    print("=" * 70)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load benchmark
    print("\n1. Loading E2G benchmark data...")
    df = load_e2g_benchmark()
    
    # Get true regulatory links
    print("\n2. Filtering for true regulatory links...")
    df_true = get_true_regulatory_links(df)
    
    # Show disease distribution
    print("\n3. Disease distribution:")
    diseases = df_true.groupby('Disease').size().sort_values(ascending=False)
    for d, n in diseases.head(10).items():
        print(f"    {d}: {n} links")
    
    # Build benchmark - test first with 100 loci
    print("\n4. Building benchmark (100 loci test)...")
    benchmark = build_regulatory_benchmark(df_true, max_loci=100)
    
    # Save results
    output_file = OUTPUT_DIR / "regulatory_benchmark_expanded.tsv"
    benchmark.to_csv(output_file, sep='\t', index=False)
    print(f"\n   Saved to {output_file}")
    
    # Evaluate
    print("\n5. Evaluating L2G performance:")
    metrics = evaluate_benchmark(benchmark)
    
    print(f"\n   RESULTS:")
    print(f"   - Total gene-locus links: {metrics['total_links']}")
    print(f"   - L2G coverage: {metrics['l2g_coverage']} ({metrics['coverage_rate']:.1%})")
    
    if 'mean_l2g' in metrics:
        print(f"   - Mean L2G (true positives): {metrics['mean_l2g']:.3f}")
        print(f"   - Median L2G (true positives): {metrics['median_l2g']:.3f}")
        print(f"   - Recall at L2G > 0.5: {metrics.get('accuracy_at_0.5', 0):.1%}")
        print(f"   - Recall at L2G > 0.3: {metrics.get('accuracy_at_0.3', 0):.1%}")
    
    # Summary
    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    
    print(f"""
This benchmark evaluates L2G's ability to identify true enhancer-gene
regulatory connections identified through the ABC (Activity-by-Contact)
model. These are non-coding GWAS loci where functional genomics evidence
supports specific gene targets.

Key insights:
- Coverage indicates what fraction of regulatory loci L2G can score
- Recall at threshold X shows what fraction of true genes have L2G > X
- This is a more challenging test than standard GWAS loci because
  regulatory mechanisms often act at distance without obvious features
""")
    
    return benchmark, metrics


if __name__ == "__main__":
    benchmark, metrics = main()
