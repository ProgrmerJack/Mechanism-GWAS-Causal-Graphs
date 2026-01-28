#!/usr/bin/env python3
"""Introspect Platform API schema to understand correct field names."""

import requests
import json

API = "https://api.platform.opentargets.org/api/v4/graphql"

def introspect_type(type_name):
    query = '{ __type(name: "%s") { name fields { name type { name kind ofType { name } } } } }' % type_name
    r = requests.post(API, json={"query": query})
    return r.json()

# Introspect key types
print("=" * 60)
print("VARIANT TYPE FIELDS:")
print("=" * 60)
variant_schema = introspect_type("Variant")
if 'data' in variant_schema and variant_schema['data'].get('__type'):
    for field in variant_schema['data']['__type']['fields'][:30]:
        print(f"  {field['name']}: {field['type']}")

print("\n" + "=" * 60)
print("SEARCH RESULT TYPE FIELDS:")
print("=" * 60)
search_schema = introspect_type("SearchResult")
if 'data' in search_schema and search_schema['data'].get('__type'):
    for field in search_schema['data']['__type'].get('fields', [])[:20]:
        print(f"  {field['name']}: {field['type']}")

# Try a simpler variant query
print("\n" + "=" * 60)
print("TESTING VARIANT QUERY:")
print("=" * 60)
variant_query = """
query {
  variant(variantId: "1_109274968_G_T") {
    id
    chromosome
    position
  }
}
"""
r = requests.post(API, json={"query": variant_query})
print(json.dumps(r.json(), indent=2))

# Try search with correct return type
print("\n" + "=" * 60)
print("TESTING SEARCH:")
print("=" * 60)
search_query = """
query {
  search(queryString: "LDL cholesterol", entityNames: ["study"], page: {size: 5, index: 0}) {
    total
    hits {
      id
      name
      entity
      score
    }
  }
}
"""
r2 = requests.post(API, json={"query": search_query})
print(json.dumps(r2.json(), indent=2))

# Query credible sets directly
print("\n" + "=" * 60)
print("TESTING CREDIBLE SETS QUERY:")
print("=" * 60)
cs_query = """
query {
  credibleSets(page: {size: 3, index: 0}) {
    count
    rows {
      studyLocusId
      region
      chromosome
      position
      study {
        id
        traitFromSource
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
r3 = requests.post(API, json={"query": cs_query})
print(json.dumps(r3.json(), indent=2))
