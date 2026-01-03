#!/usr/bin/env python3
"""
Query Open Targets Platform API for L2G scores.
Uses the new Platform API endpoint (api.platform.opentargets.org)
since the Genetics API has been merged into Platform.

This script:
1. Queries the Platform GraphQL API for L2G data
2. Retrieves credibleSets with L2G scores
3. Maps to our benchmark loci for official baseline comparison
"""

import requests
import json
import pandas as pd
import time
from pathlib import Path

# API endpoint
PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

def get_api_info():
    """Get API version and available types."""
    query = """
    {
      meta {
        name
        apiVersion {
          x
          y
          z
        }
        dataVersion {
          year
          month
          iteration
        }
      }
    }
    """
    response = requests.post(PLATFORM_API, json={"query": query}, timeout=30)
    return response.json()

def get_available_studies(limit=10):
    """Query available GWAS studies."""
    query = """
    query getStudies($size: Int!) {
      gwasStudies(page: {size: $size, index: 0}) {
        count
        rows {
          studyId
          traitFromSource
          nCases
          nSamples
        }
      }
    }
    """
    variables = {"size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=30)
    return response.json()

def get_credible_sets_for_study(study_id, limit=50):
    """Get credible sets and L2G data for a specific study."""
    query = """
    query getCredibleSets($studyId: String!, $size: Int!) {
      gwasCredibleSets(studyId: $studyId, page: {size: $size, index: 0}) {
        count
        rows {
          studyLocusId
          variant {
            id
            chromosome
            position
            refAllele
            altAllele
          }
          locus {
            variantId
            posteriorProbability
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
    variables = {"studyId": study_id, "size": limit}
    response = requests.post(PLATFORM_API, json={"query": query, "variables": variables}, timeout=60)
    return response.json()

def introspect_schema():
    """Introspect the GraphQL schema to find L2G-related types."""
    query = """
    {
      __schema {
        types {
          name
          description
          fields {
            name
            description
          }
        }
      }
    }
    """
    response = requests.post(PLATFORM_API, json={"query": query}, timeout=30)
    data = response.json()
    
    # Filter for types related to L2G or credible sets
    l2g_related = []
    for t in data.get('data', {}).get('__schema', {}).get('types', []):
        name = t.get('name', '')
        desc = t.get('description', '') or ''
        if any(kw in name.lower() or kw in desc.lower() for kw in ['l2g', 'locus', 'gene', 'credible', 'gwas', 'variant']):
            if not name.startswith('__'):
                l2g_related.append({
                    'name': name,
                    'description': desc[:200] if desc else '',
                    'fields': [f['name'] for f in (t.get('fields') or [])][:10]
                })
    
    return l2g_related

def explore_gwas_queries():
    """Explore available GWAS-related queries."""
    query = """
    {
      __type(name: "Query") {
        fields {
          name
          description
          args {
            name
            type {
              name
              kind
            }
          }
        }
      }
    }
    """
    response = requests.post(PLATFORM_API, json={"query": query}, timeout=30)
    data = response.json()
    
    # Filter for GWAS-related queries
    gwas_queries = []
    for f in data.get('data', {}).get('__type', {}).get('fields', []):
        name = f.get('name', '')
        desc = f.get('description', '') or ''
        if any(kw in name.lower() or kw in desc.lower() for kw in ['gwas', 'locus', 'variant', 'credible', 'l2g', 'study']):
            gwas_queries.append({
                'name': name,
                'description': desc[:150] if desc else '',
                'args': [a['name'] for a in f.get('args', [])]
            })
    
    return gwas_queries

def main():
    """Main execution."""
    print("=" * 70)
    print("OPEN TARGETS PLATFORM API - L2G Data Query")
    print("=" * 70)
    
    # Step 1: Get API info
    print("\n1. API Information:")
    api_info = get_api_info()
    if 'data' in api_info:
        meta = api_info['data']['meta']
        print(f"   Name: {meta['name']}")
        v = meta['apiVersion']
        print(f"   API Version: {v['x']}.{v['y']}.{v['z']}")
        dv = meta['dataVersion']
        print(f"   Data Version: {dv['year']}.{dv['month']}.{dv['iteration']}")
    else:
        print(f"   Error: {api_info}")
        return
    
    # Step 2: Explore available queries
    print("\n2. Available GWAS-related queries:")
    gwas_queries = explore_gwas_queries()
    for q in gwas_queries[:15]:
        print(f"   - {q['name']}: {q['description'][:60]}...")
        if q['args']:
            print(f"     Args: {', '.join(q['args'][:5])}")
    
    # Step 3: Get sample studies
    print("\n3. Sample GWAS studies:")
    studies = get_available_studies(limit=5)
    if 'data' in studies and 'gwasStudies' in studies['data']:
        study_data = studies['data']['gwasStudies']
        print(f"   Total studies: {study_data['count']}")
        for row in study_data['rows']:
            print(f"   - {row['studyId']}: {row['traitFromSource']} (n={row['nSamples']})")
    else:
        print(f"   Response: {json.dumps(studies, indent=2)[:500]}")
    
    # Step 4: Try to get credible sets for a cardiometabolic study
    print("\n4. Querying credible sets for sample study...")
    # Try to find a lipid-related study
    if 'data' in studies and studies['data'].get('gwasStudies', {}).get('rows'):
        study_id = studies['data']['gwasStudies']['rows'][0]['studyId']
        print(f"   Using study: {study_id}")
        
        cs_data = get_credible_sets_for_study(study_id, limit=5)
        if 'data' in cs_data:
            print(f"   Response keys: {list(cs_data['data'].keys())}")
            if cs_data['data'].get('gwasCredibleSets'):
                cs = cs_data['data']['gwasCredibleSets']
                print(f"   Total credible sets: {cs.get('count', 'N/A')}")
                for row in cs.get('rows', [])[:3]:
                    print(f"   - Locus: {row.get('studyLocusId', 'N/A')}")
                    if row.get('l2GPredictions'):
                        l2g = row['l2GPredictions']['rows']
                        for pred in l2g[:3]:
                            gene = pred['target']['approvedSymbol']
                            score = pred['score']
                            print(f"     L2G: {gene} = {score:.3f}")
        else:
            print(f"   Error: {json.dumps(cs_data, indent=2)[:500]}")
    
    # Step 5: Introspect schema for L2G types
    print("\n5. L2G-related schema types:")
    l2g_types = introspect_schema()
    for t in l2g_types[:10]:
        print(f"   - {t['name']}: {t['description'][:50]}...")
        if t['fields']:
            print(f"     Fields: {', '.join(t['fields'][:5])}")
    
    print("\n" + "=" * 70)
    print("API exploration complete. See above for available data.")
    print("=" * 70)

if __name__ == "__main__":
    main()
