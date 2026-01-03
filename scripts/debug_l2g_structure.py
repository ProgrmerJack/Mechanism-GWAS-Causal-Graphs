#!/usr/bin/env python3
"""Debug L2G query structure."""
import requests
import json

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

variant_id = "1_999842_C_A"  # rs2298214

# Full debug query
query = """
query getL2G($id: String!) {
  variant(variantId: $id) {
    id
    credibleSets(page: {size: 5, index: 0}) {
      count
      rows {
        studyLocusId
        study { studyId traitFromSource }
        l2GPredictions(page: {size: 10, index: 0}) {
          count
          rows {
            target { approvedSymbol id }
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
    json={'query': query, 'variables': {'id': variant_id}},
    timeout=30
)
data = response.json()

print("Response:")
print(json.dumps(data, indent=2))

# Check if we need different structure
variant = data.get('data', {}).get('variant', {})
if variant:
    print(f"\nCredible sets count: {variant.get('credibleSets', {}).get('count', 0)}")
    rows = variant.get('credibleSets', {}).get('rows', [])
    for i, cs in enumerate(rows):
        print(f"\nCredible Set {i+1}:")
        print(f"  StudyLocusId: {cs.get('studyLocusId')}")
        print(f"  Study: {cs.get('study')}")
        l2g = cs.get('l2GPredictions', {})
        print(f"  L2G count: {l2g.get('count', 0)}")
        for p in l2g.get('rows', [])[:3]:
            gene = p.get('target', {}).get('approvedSymbol', 'N/A')
            score = p.get('score', 0)
            print(f"    {gene}: {score}")
