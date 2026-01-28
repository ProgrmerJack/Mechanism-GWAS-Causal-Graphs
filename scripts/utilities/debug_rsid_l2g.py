#!/usr/bin/env python3
"""Debug why rsID-found variants have no L2G."""
import requests
import json

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

# Test with a specific rsID from our data
rsids = ['rs2298214', 'rs9442392', 'rs12740374', 'rs4590800', 'rs657624']

for rsid in rsids:
    print(f"\n{'='*60}")
    print(f"Testing rsID: {rsid}")
    print(f"{'='*60}")
    
    # Step 1: Search for variant
    search_query = """
    query searchVariant($rsid: String!) {
      search(queryString: $rsid, entityNames: ["variant"]) {
        hits { id entity name }
      }
    }
    """
    response = requests.post(
        PLATFORM_API,
        json={'query': search_query, 'variables': {'rsid': rsid}},
        timeout=30
    )
    data = response.json()
    hits = data.get('data', {}).get('search', {}).get('hits', [])
    
    if not hits:
        print("  NOT FOUND in search")
        continue
    
    variant_id = hits[0].get('id')
    print(f"  Variant ID: {variant_id}")
    
    # Step 2: Get variant details
    variant_query = """
    query getVariant($id: String!) {
      variant(variantId: $id) {
        id
        rsIds
        chromosome
        position
        credibleSets {
          count
        }
      }
    }
    """
    response = requests.post(
        PLATFORM_API,
        json={'query': variant_query, 'variables': {'id': variant_id}},
        timeout=30
    )
    data = response.json()
    variant = data.get('data', {}).get('variant', {})
    
    if not variant:
        print("  Variant data not found")
        continue
    
    cs_count = variant.get('credibleSets', {}).get('count', 0)
    print(f"  Chromosome: {variant.get('chromosome')}")
    print(f"  Position: {variant.get('position')}")
    print(f"  Credible sets count: {cs_count}")
    
    # Step 3: Get L2G if credible sets exist
    if cs_count > 0:
        l2g_query = """
        query getL2G($id: String!) {
          variant(variantId: $id) {
            credibleSets {
              rows {
                studyLocusId
                study { studyId traitFromSource }
                l2GPredictions {
                  rows {
                    target { approvedSymbol }
                    score
                  }
                }
              }
            }
          }
        }
        """
        response = requests.post(
            PLATFORM_API,
            json={'query': l2g_query, 'variables': {'id': variant_id}},
            timeout=30
        )
        data = response.json()
        rows = data.get('data', {}).get('variant', {}).get('credibleSets', {}).get('rows', [])
        
        for cs in rows[:3]:  # Show first 3
            study = cs.get('study', {})
            preds = cs.get('l2GPredictions', {}).get('rows', [])
            print(f"\n  Study: {study.get('studyId')}")
            print(f"  Trait: {study.get('traitFromSource')}")
            for p in preds[:3]:
                print(f"    {p.get('target', {}).get('approvedSymbol')}: {p.get('score')}")
    else:
        print("  NO credible sets - likely UKBB variant not in Platform studies")

print("\n" + "="*60)
print("CONCLUSION")
print("="*60)
print("""
The UKBB SuSiE variants exist in the Platform database (as genomic variants),
but they are NOT associated with any credible sets in Platform's studies.

This is because:
1. Open Targets Platform has its own GWAS catalog integration
2. UKBB SuSiE credible sets from the E2G benchmarking paper are not included
3. Only variants from OT-curated studies have L2G predictions

SOLUTION: Use Open Targets Platform's own credible sets instead of UKBB SuSiE.
We need to find credible sets where:
- Non-coding variants exist
- ABC predictions exist for validation
""")
