#!/usr/bin/env python3
"""Debug L2G query structure - fix API schema."""
import requests
import json

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

# First introspect the Study type
introspect_query = """
{
  __type(name: "Study") {
    fields {
      name
      type { name kind ofType { name } }
    }
  }
}
"""

response = requests.post(
    PLATFORM_API,
    json={'query': introspect_query},
    timeout=30
)
data = response.json()
print("Study type fields:")
for field in data.get('data', {}).get('__type', {}).get('fields', []):
    print(f"  {field['name']}: {field['type']}")

# Now query with correct fields
variant_id = "1_999842_C_A"  # rs2298214

query = """
query getL2G($id: String!) {
  variant(variantId: $id) {
    id
    credibleSets(page: {size: 5, index: 0}) {
      count
      rows {
        studyLocusId
        study { id studyType traitFromSource }
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

if 'errors' in data:
    print("\nErrors:", data['errors'])
else:
    variant = data.get('data', {}).get('variant', {})
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
