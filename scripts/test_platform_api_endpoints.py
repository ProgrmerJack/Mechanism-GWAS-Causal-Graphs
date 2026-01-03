#!/usr/bin/env python3
"""Test Platform API endpoints for credible sets."""
import requests
import json

# Test 1: Get schema to see available queries for credibleSets
introspect_query = """
query {
  __type(name: "Query") {
    fields {
      name
      args {
        name
        type { name kind ofType { name kind } }
      }
    }
  }
}
"""

print("=" * 70)
print("TESTING PLATFORM API FOR CREDIBLE SETS")
print("=" * 70)

response = requests.post(
    'https://api.platform.opentargets.org/api/v4/graphql',
    json={'query': introspect_query},
    headers={'Content-Type': 'application/json'},
    timeout=30
)

print("\n1. Available Query Fields:")
data = response.json()
if 'data' in data and data['data']:
    for field in data['data']['__type']['fields']:
        if 'credible' in field['name'].lower() or 'study' in field['name'].lower():
            args = [f"{a['name']}" for a in field.get('args', [])]
            print(f"  - {field['name']}({', '.join(args)})")

# Test 2: Try credibleSets by region
print("\n2. Testing credibleSets by region query:")

region_query = """
query testRegion {
  credibleSets(chromosome: "1", start: 109000000, end: 110000000) {
    count
  }
}
"""

response = requests.post(
    'https://api.platform.opentargets.org/api/v4/graphql',
    json={'query': region_query},
    headers={'Content-Type': 'application/json'},
    timeout=30
)

print(f"  Status: {response.status_code}")
print(f"  Response: {json.dumps(response.json(), indent=2)[:500]}")

# Test 3: Try studyLocusId based query
print("\n3. Testing alternative: search for variants directly:")

variant_query = """
query testVariant {
  search(queryString: "rs12740374", entityNames: ["variant"]) {
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
    'https://api.platform.opentargets.org/api/v4/graphql',
    json={'query': variant_query},
    headers={'Content-Type': 'application/json'},
    timeout=30
)

print(f"  Status: {response.status_code}")
result = response.json()
print(f"  Search results: {json.dumps(result, indent=2)[:1000]}")

# Test 4: Query variant directly
print("\n4. Testing direct variant query (SORT1 locus):")

direct_query = """
query getVariant {
  variant(variantId: "1_109274968_G_T") {
    id
    rsIds
    chromosome
    position
    credibleSets {
      count
      rows {
        studyLocusId
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
    'https://api.platform.opentargets.org/api/v4/graphql',
    json={'query': direct_query},
    headers={'Content-Type': 'application/json'},
    timeout=30
)

print(f"  Status: {response.status_code}")
result = response.json()
if 'data' in result and result['data'] and 'variant' in result['data']:
    var = result['data']['variant']
    if var:
        print(f"  Found variant: {var.get('id')}")
        cs = var.get('credibleSets', {})
        print(f"  Credible sets count: {cs.get('count', 0)}")
        if cs.get('rows'):
            for row in cs['rows'][:3]:
                print(f"    - {row.get('studyLocusId')}")
                l2g = row.get('l2GPredictions', {}).get('rows', [])
                for pred in l2g[:3]:
                    gene = pred.get('target', {}).get('approvedSymbol', 'N/A')
                    score = pred.get('score', 0)
                    print(f"      L2G: {gene} = {score:.3f}")
else:
    print(f"  Error: {json.dumps(result, indent=2)[:500]}")

# Test 5: Look for study-based credible sets 
print("\n5. Testing gwasCredibleSet query (by study):")

study_query = """
query testStudyCS {
  study(studyId: "GCST006612") {
    id
    traitFromSourceMappedIds
    credibleSets {
      count
    }
  }
}
"""

response = requests.post(
    'https://api.platform.opentargets.org/api/v4/graphql',
    json={'query': study_query},
    headers={'Content-Type': 'application/json'},
    timeout=30
)

print(f"  Status: {response.status_code}")
print(f"  Response: {json.dumps(response.json(), indent=2)[:1000]}")

print("\n" + "=" * 70)
print("CONCLUSION: Use variant-based queries with known lead variants")
print("=" * 70)
