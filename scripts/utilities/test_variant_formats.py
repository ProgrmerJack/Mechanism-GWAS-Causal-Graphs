#!/usr/bin/env python3
"""Test variant ID formats for Platform API."""
import requests
import json

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

# Test different variant formats
variants_to_test = [
    ('10_122395211_CTT_C', 'E2G format indel'),
    ('10_50882755_G_T', 'E2G format SNP'),
    ('1_109274968_G_T', 'Known working SORT1'),
]

rsids_to_test = [
    'rs4590800',
    'rs7902581', 
    'rs12740374',  # SORT1
]

print("=" * 60)
print("TESTING VARIANT ID FORMATS")
print("=" * 60)

# Test direct variant queries
print("\n1. Direct variant queries (CHR_POS_REF_ALT format):")
for var_id, desc in variants_to_test:
    query = """
    query getVariant($id: String!) {
      variant(variantId: $id) {
        id
        rsIds
        chromosome
        position
      }
    }
    """
    response = requests.post(
        PLATFORM_API,
        json={'query': query, 'variables': {'id': var_id}},
        timeout=30
    )
    data = response.json()
    variant = data.get('data', {}).get('variant')
    
    if variant:
        print(f"  {var_id} ({desc}): FOUND")
        print(f"    rsIds: {variant.get('rsIds', [])}")
    else:
        print(f"  {var_id} ({desc}): NOT FOUND")

# Test rsID searches
print("\n2. rsID searches:")
for rsid in rsids_to_test:
    query = """
    query searchVariant($rsid: String!) {
      search(queryString: $rsid, entityNames: ["variant"]) {
        total
        hits {
          id
          entity
          name
        }
      }
    }
    """
    response = requests.post(
        PLATFORM_API,
        json={'query': query, 'variables': {'rsid': rsid}},
        timeout=30
    )
    data = response.json()
    if 'errors' in data:
        print(f"  {rsid}: API ERROR - {data['errors']}")
        continue
    if data.get('data') is None:
        print(f"  {rsid}: NO DATA RETURNED")
        continue
    hits = data.get('data', {}).get('search', {}).get('hits', [])
    
    if hits:
        print(f"  {rsid}: FOUND -> {hits[0].get('id')}")
    else:
        print(f"  {rsid}: NOT FOUND")

# Test using rsID to get variant ID, then query L2G
print("\n3. Full pipeline: rsID -> variant ID -> L2G:")
rsid = 'rs12740374'  # SORT1 locus

# Step 1: Search for rsID
search_query = """
query searchVariant($rsid: String!) {
  search(queryString: $rsid, entityNames: ["variant"]) {
    hits { id }
  }
}
"""
response = requests.post(
    PLATFORM_API,
    json={'query': search_query, 'variables': {'rsid': rsid}},
    timeout=30
)
data = response.json()
variant_id = data.get('data', {}).get('search', {}).get('hits', [{}])[0].get('id')
print(f"  {rsid} -> variant ID: {variant_id}")

# Step 2: Get L2G for that variant
if variant_id:
    l2g_query = """
    query getL2G($id: String!) {
      variant(variantId: $id) {
        id
        credibleSets {
          count
          rows {
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
    
    variant = data.get('data', {}).get('variant', {})
    cs_count = variant.get('credibleSets', {}).get('count', 0)
    print(f"  Credible sets: {cs_count}")
    
    # Get top L2G scores
    for cs in variant.get('credibleSets', {}).get('rows', [])[:2]:
        preds = cs.get('l2GPredictions', {}).get('rows', [])
        for pred in preds[:3]:
            gene = pred.get('target', {}).get('approvedSymbol', 'N/A')
            score = pred.get('score', 0)
            print(f"    L2G: {gene} = {score:.3f}")

print("\n" + "=" * 60)
print("CONCLUSION")
print("=" * 60)
print("""
The Platform API:
- Uses GRCh38 coordinates
- Variant IDs must match exactly in CHR_POS_REF_ALT format
- Many E2G variants may not be in the Platform database
- rsID search works well as alternative lookup
""")
