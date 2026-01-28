#!/usr/bin/env python3
"""
Query Open Targets Platform API v25 for L2G scores.
Updated to use correct GraphQL schema discovered via introspection.
"""

import requests
import json
import pandas as pd
import time
from pathlib import Path

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

def get_studies(limit=20):
    """Query available studies."""
    query = """
    query getStudies($size: Int!) {
      studies(page: {size: $size, index: 0}) {
        count
        rows {
          studyId
          traitFromSource
          nCases
          nSamples
          studyType
        }
      }
    }
    """
    variables = {"size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    return response.json()

def get_credible_sets(study_ids=None, variant_ids=None, limit=50):
    """Get credible sets with L2G predictions."""
    query = """
    query getCredibleSets($studyIds: [String!], $variantIds: [String!], $size: Int!) {
      credibleSets(studyIds: $studyIds, variantIds: $variantIds, page: {size: $size, index: 0}) {
        count
        rows {
          studyLocusId
          region
          pValueMantissa
          pValueExponent
          ldSet {
            tagVariantId
            r2Overall
          }
          locus {
            variantId
            posteriorProbability
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
    variables = {"studyIds": study_ids, "variantIds": variant_ids, "size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=60)
    return response.json()

def get_credible_set_by_id(study_locus_id):
    """Get single credible set by ID with full L2G data."""
    query = """
    query getCredibleSet($studyLocusId: String!) {
      credibleSet(studyLocusId: $studyLocusId) {
        studyLocusId
        region
        study {
          studyId
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
    """
    variables = {"studyLocusId": study_locus_id}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    return response.json()

def search_gwas_studies_by_trait(trait_query, limit=10):
    """Search for GWAS studies by trait name."""
    query = """
    query searchStudies($query: String!, $size: Int!) {
      search(queryString: $query, entityNames: ["study"], page: {size: $size, index: 0}) {
        total
        hits {
          ... on Study {
            studyId
            traitFromSource
            studyType
            nSamples
            nCases
          }
        }
      }
    }
    """
    variables = {"query": trait_query, "size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    return response.json()

def query_variant_l2g(variant_id):
    """Query L2G data for a specific variant."""
    query = """
    query getVariantL2G($variantId: String!) {
      variant(variantId: $variantId) {
        id
        chromosome
        position
        refAllele
        altAllele
        gnomadVariantId
      }
      credibleSets(variantIds: [$variantId], page: {size: 20, index: 0}) {
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
    variables = {"variantId": variant_id}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    return response.json()

def main():
    """Main execution - explore API and get sample L2G data."""
    print("=" * 70)
    print("OPEN TARGETS PLATFORM API - L2G Data Query v2")
    print("=" * 70)
    
    # Step 1: Get sample studies
    print("\n1. Sample GWAS studies:")
    studies = get_studies(limit=10)
    if 'data' in studies and studies['data'].get('studies'):
        study_data = studies['data']['studies']
        print(f"   Total studies available: {study_data['count']}")
        gwas_studies = [s for s in study_data['rows'] if s.get('studyType') == 'gwas']
        print(f"   GWAS studies in sample: {len(gwas_studies)}")
        for s in gwas_studies[:5]:
            print(f"   - {s['studyId']}: {s['traitFromSource']} (n={s.get('nSamples', 'N/A')})")
    else:
        print(f"   Error: {json.dumps(studies, indent=2)[:500]}")
        return
    
    # Step 2: Search for lipid-related studies
    print("\n2. Searching for lipid GWAS studies:")
    lipid_search = search_gwas_studies_by_trait("LDL cholesterol", limit=5)
    if 'data' in lipid_search and lipid_search['data'].get('search'):
        hits = lipid_search['data']['search']['hits']
        print(f"   Found {lipid_search['data']['search']['total']} studies")
        for hit in hits[:5]:
            if hit:  # Some hits may be None
                print(f"   - {hit.get('studyId', 'N/A')}: {hit.get('traitFromSource', 'N/A')}")
    
    # Step 3: Get credible sets for a GWAS study
    if gwas_studies:
        study_id = gwas_studies[0]['studyId']
        print(f"\n3. Getting credible sets for study: {study_id}")
        cs_data = get_credible_sets(study_ids=[study_id], limit=10)
        if 'data' in cs_data and cs_data['data'].get('credibleSets'):
            cs = cs_data['data']['credibleSets']
            print(f"   Total credible sets: {cs['count']}")
            for row in cs.get('rows', [])[:5]:
                locus_id = row.get('studyLocusId', 'N/A')
                region = row.get('region', 'N/A')
                print(f"\n   Locus: {locus_id}")
                print(f"   Region: {region}")
                
                # L2G predictions
                l2g_preds = row.get('l2GPredictions', {}).get('rows', [])
                if l2g_preds:
                    print("   L2G Predictions:")
                    for pred in l2g_preds[:5]:
                        gene = pred['target']['approvedSymbol']
                        score = pred['score']
                        print(f"     {gene}: {score:.4f}")
                else:
                    print("   No L2G predictions available")
        else:
            print(f"   Error: {json.dumps(cs_data, indent=2)[:800]}")
    
    # Step 4: Query L2G for a specific well-known variant (SORT1 locus)
    print("\n4. Querying L2G for SORT1 locus variant (rs12740374):")
    # Format: CHROM_POS_REF_ALT
    # rs12740374 is at chr1:109274968 (GRCh38), G>T
    variant_id = "1_109274968_G_T"
    variant_data = query_variant_l2g(variant_id)
    if 'data' in variant_data:
        if variant_data['data'].get('variant'):
            v = variant_data['data']['variant']
            print(f"   Variant: {v['id']}")
            print(f"   Position: chr{v['chromosome']}:{v['position']}")
        
        cs = variant_data['data'].get('credibleSets', {})
        if cs.get('rows'):
            print(f"   Credible sets containing this variant: {cs['count']}")
            for row in cs['rows'][:3]:
                study = row.get('study', {})
                print(f"\n   Study: {study.get('studyId')} - {study.get('traitFromSource')}")
                l2g_preds = row.get('l2GPredictions', {}).get('rows', [])
                for pred in l2g_preds[:5]:
                    gene = pred['target']['approvedSymbol']
                    score = pred['score']
                    print(f"     L2G: {gene} = {score:.4f}")
    else:
        print(f"   Error or no data: {json.dumps(variant_data, indent=2)[:500]}")
    
    print("\n" + "=" * 70)
    print("Platform API L2G query complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
