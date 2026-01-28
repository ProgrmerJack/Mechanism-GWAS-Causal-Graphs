#!/usr/bin/env python3
"""
Query Open Targets Platform API for L2G scores and variant data.
Uses correct GraphQL schema based on API introspection (v25.4.4).

Schema notes:
- Study type has 'id' not 'studyId'
- CredibleSet has l2GPredictions field
- Variant search uses variantId format: CHROM_POS_REF_ALT
"""

import requests
import json
import pandas as pd
import time
from pathlib import Path

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

def get_variant_info(variant_id):
    """Query variant info and associated credible sets with L2G data."""
    query = """
    query getVariant($variantId: String!) {
      variant(variantId: $variantId) {
        id
        chromosome
        position
        refAllele
        altAllele
        gnomadVariantId
        mostSevereConsequence
        inSilicoPredictors {
          method
          assessment
          score
        }
      }
      credibleSets(variantIds: [$variantId], page: {size: 50, index: 0}) {
        count
        rows {
          studyLocusId
          region
          chromosome
          position
          finemappingMethod
          pValueMantissa
          pValueExponent
          beta
          study {
            id
            traitFromSource
            studyType
            nSamples
          }
          l2GPredictions {
            count
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
    variables = {"variantId": variant_id}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=60)
    return response.json()

def search_for_study(search_term, limit=20):
    """Search for GWAS studies by trait name."""
    query = """
    query searchStudies($term: String!, $size: Int!) {
      search(queryString: $term, entityNames: ["study"], page: {size: $size, index: 0}) {
        total
        hits {
          ... on Study {
            id
            traitFromSource
            studyType
            nSamples
            nCases
          }
        }
      }
    }
    """
    variables = {"term": search_term, "size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    return response.json()

def get_credible_sets_for_study(study_id, limit=100):
    """Get all credible sets and L2G predictions for a study."""
    query = """
    query getStudyCredibleSets($studyId: String!, $size: Int!) {
      credibleSets(studyIds: [$studyId], page: {size: $size, index: 0}) {
        count
        rows {
          studyLocusId
          region
          chromosome
          position
          finemappingMethod
          pValueMantissa
          pValueExponent
          study {
            id
            traitFromSource
          }
          locus {
            variantId
            posteriorProbability
          }
          l2GPredictions {
            count
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
    variables = {"studyId": study_id, "size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=120)
    return response.json()

def rsid_to_variant_id(rsid):
    """Convert rsID to Open Targets variant ID format (CHROM_POS_REF_ALT).
    
    Note: This uses a simple lookup API. For production use, consider dbSNP queries.
    """
    # Use Open Targets search API for rsID lookup
    query = """
    query searchVariant($term: String!) {
      search(queryString: $term, entityNames: ["variant"], page: {size: 1, index: 0}) {
        total
        hits {
          ... on Variant {
            id
            chromosome
            position
            refAllele
            altAllele
          }
        }
      }
    }
    """
    variables = {"term": rsid}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    result = response.json()
    if result.get('data', {}).get('search', {}).get('hits'):
        return result['data']['search']['hits'][0].get('id')
    return None

def get_l2g_for_benchmark_variants(benchmark_df):
    """Query L2G scores for benchmark variants.
    
    Args:
        benchmark_df: DataFrame with columns including 'lead_snp' or 'chr', 'pos_hg38'
    
    Returns:
        DataFrame with L2G scores for each variant
    """
    results = []
    
    for idx, row in benchmark_df.iterrows():
        # Try rsID first
        variant_id = None
        if 'lead_snp' in row and pd.notna(row['lead_snp']):
            rsid = row['lead_snp']
            print(f"Looking up {rsid}...", end=" ")
            variant_id = rsid_to_variant_id(rsid)
            if variant_id:
                print(f"Found: {variant_id}")
            else:
                print("Not found")
        
        # If we have the variant ID, query L2G
        if variant_id:
            var_data = get_variant_info(variant_id)
            if 'data' in var_data:
                cs_data = var_data['data'].get('credibleSets', {})
                for cs in cs_data.get('rows', []):
                    l2g_preds = cs.get('l2GPredictions', {}).get('rows', [])
                    study = cs.get('study', {})
                    for pred in l2g_preds:
                        results.append({
                            'benchmark_locus': row.get('locus_id', idx),
                            'rsid': row.get('lead_snp', ''),
                            'variant_id': variant_id,
                            'study_id': study.get('id', ''),
                            'trait': study.get('traitFromSource', ''),
                            'gene_symbol': pred['target']['approvedSymbol'],
                            'gene_id': pred['target']['id'],
                            'l2g_score': pred['score'],
                            'region': cs.get('region', ''),
                        })
            time.sleep(0.5)  # Rate limiting
    
    return pd.DataFrame(results)

def main():
    print("=" * 70)
    print("OPEN TARGETS PLATFORM API - L2G Data Extraction")
    print("=" * 70)
    
    # Test 1: Look up a known variant (SORT1 locus rs12740374)
    print("\n1. Testing variant lookup for rs12740374 (SORT1 locus):")
    variant_id = rsid_to_variant_id("rs12740374")
    print(f"   Variant ID: {variant_id}")
    
    if variant_id:
        var_data = get_variant_info(variant_id)
        if 'data' in var_data and var_data['data'].get('variant'):
            v = var_data['data']['variant']
            print(f"   Position: chr{v['chromosome']}:{v['position']}")
            print(f"   Alleles: {v['refAllele']}/{v['altAllele']}")
            print(f"   Consequence: {v.get('mostSevereConsequence', 'N/A')}")
        
        cs_data = var_data['data'].get('credibleSets', {})
        print(f"\n   Credible sets containing variant: {cs_data.get('count', 0)}")
        for cs in cs_data.get('rows', [])[:3]:
            study = cs.get('study', {})
            print(f"\n   Study: {study.get('id', 'N/A')}")
            print(f"   Trait: {study.get('traitFromSource', 'N/A')}")
            l2g_preds = cs.get('l2GPredictions', {}).get('rows', [])
            print(f"   L2G predictions ({len(l2g_preds)}):")
            for pred in l2g_preds[:5]:
                print(f"     {pred['target']['approvedSymbol']}: {pred['score']:.4f}")
    
    # Test 2: Search for LDL cholesterol GWAS
    print("\n" + "=" * 70)
    print("2. Searching for LDL cholesterol GWAS studies:")
    search_results = search_for_study("LDL cholesterol", limit=10)
    if 'data' in search_results and search_results['data'].get('search'):
        hits = search_results['data']['search'].get('hits', [])
        print(f"   Found {search_results['data']['search']['total']} studies")
        for hit in hits[:5]:
            if hit:
                print(f"   - {hit.get('id', 'N/A')}: {hit.get('traitFromSource', 'N/A')} (n={hit.get('nSamples', 'N/A')})")
    
    # Test 3: Get credible sets for a specific GWAS study
    print("\n" + "=" * 70)
    print("3. Getting credible sets for an LDL GWAS study:")
    if 'data' in search_results and search_results['data'].get('search', {}).get('hits'):
        first_study = None
        for hit in search_results['data']['search']['hits']:
            if hit and hit.get('studyType') == 'gwas':
                first_study = hit
                break
        
        if first_study:
            study_id = first_study['id']
            print(f"   Study: {study_id}")
            cs_results = get_credible_sets_for_study(study_id, limit=10)
            if 'data' in cs_results and cs_results['data'].get('credibleSets'):
                cs = cs_results['data']['credibleSets']
                print(f"   Total credible sets: {cs.get('count', 0)}")
                for row in cs.get('rows', [])[:3]:
                    print(f"\n   Region: {row.get('region', 'N/A')}")
                    print(f"   P-value: {row.get('pValueMantissa', 'N/A')}e{row.get('pValueExponent', 'N/A')}")
                    l2g_preds = row.get('l2GPredictions', {}).get('rows', [])
                    if l2g_preds:
                        print(f"   Top L2G genes:")
                        for pred in l2g_preds[:3]:
                            print(f"     {pred['target']['approvedSymbol']}: {pred['score']:.4f}")
    
    # Test 4: Load benchmark and query L2G
    print("\n" + "=" * 70)
    print("4. Testing benchmark variant L2G lookup:")
    benchmark_path = Path(r"C:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs\data\processed\baselines\post2021_independent_benchmark_FINAL.tsv")
    if benchmark_path.exists():
        benchmark_df = pd.read_csv(benchmark_path, sep='\t')
        print(f"   Loaded {len(benchmark_df)} benchmark loci")
        
        # Query first 3 for testing
        sample_df = benchmark_df.head(3)
        l2g_results = get_l2g_for_benchmark_variants(sample_df)
        
        if not l2g_results.empty:
            print(f"\n   Retrieved {len(l2g_results)} L2G predictions:")
            print(l2g_results[['benchmark_locus', 'gene_symbol', 'l2g_score']].head(10).to_string())
    
    print("\n" + "=" * 70)
    print("Platform API L2G extraction complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
