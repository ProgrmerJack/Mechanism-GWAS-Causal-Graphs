#!/usr/bin/env python3
"""Test Platform API variant lookup."""

import requests
import json

API = "https://api.platform.opentargets.org/api/v4/graphql"

# rs12740374 is at chr1:109817590 (hg19) or chr1:109274968 (hg38)
# Platform API uses hg38, format: CHROM_POS_REF_ALT

# Try hg38 position for rs12740374
variant_id_38 = "1_109274968_G_T"

query = """
query getVariant($variantId: String!) {
  variant(variantId: $variantId) {
    id
    chromosome
    position
    refAllele
    altAllele
    gnomadVariantId
  }
}
"""

print(f"Testing variant ID: {variant_id_38}")
r = requests.post(API, json={"query": query, "variables": {"variantId": variant_id_38}})
result = r.json()
print(f"Response: {json.dumps(result, indent=2)}")

# Also try introspecting the variant search endpoint
print("\n" + "=" * 50)
print("Testing search API:")
search_query = """
query {
  search(queryString: "rs12740374", entityNames: ["variant"], page: {size: 5, index: 0}) {
    total
    hits {
      ... on Variant {
        id
        chromosome
        position
      }
    }
  }
}
"""
r2 = requests.post(API, json={"query": search_query})
print(json.dumps(r2.json(), indent=2))

# Test study search
print("\n" + "=" * 50)
print("Testing study search:")
study_query = """
query {
  search(queryString: "LDL cholesterol", entityNames: ["study"], page: {size: 5, index: 0}) {
    total
    hits {
      ... on Study {
        id
        traitFromSource
        studyType
      }
    }
  }
}
"""
r3 = requests.post(API, json={"query": study_query})
print(json.dumps(r3.json(), indent=2))
