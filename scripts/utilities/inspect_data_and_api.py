#!/usr/bin/env python3
"""Inspect cS2G coordinate system and Platform API schema."""

import requests
import gzip
from pathlib import Path

# 1. Get Platform API Study fields
print("=" * 60)
print("PLATFORM API SCHEMA")
print("=" * 60)

query = '{ __type(name: "Study") { fields { name } } }'
r = requests.post('https://api.platform.opentargets.org/api/v4/graphql', 
                  json={'query': query})
fields = [f['name'] for f in r.json()['data']['__type']['fields']]
print("Study fields:", fields[:15], "...")
print()

# 2. Get CredibleSet fields
query2 = '{ __type(name: "CredibleSet") { fields { name } } }'
r2 = requests.post('https://api.platform.opentargets.org/api/v4/graphql', 
                   json={'query': query2})
cs_fields = [f['name'] for f in r2.json()['data']['__type']['fields']]
print("CredibleSet fields:", cs_fields)
print()

# 3. Read cS2G file header
print("=" * 60)
print("cS2G FILE FORMAT (chr1)")
print("=" * 60)
cs2g_path = Path(r'C:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs\data\external\cS2G\cS2G_extracted\cS2G_UKBB\cS2G.1.SGscore.gz')
if cs2g_path.exists():
    print(f"File: {cs2g_path.name}")
    with gzip.open(cs2g_path, 'rt') as f:
        for i in range(6):
            line = f.readline().strip()
            print(f"Line {i}: {line[:100]}{'...' if len(line) > 100 else ''}")
else:
    print("cS2G file not found")

# 4. Read allsnps.txt header  
print("\n" + "=" * 60)
print("cS2G allsnps.txt (coordinate reference)")
print("=" * 60)
allsnps_path = Path(r'C:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs\data\external\cS2G\cS2G_extracted\cS2G_UKBB\allsnps.txt.gz')
if allsnps_path.exists():
    print(f"File: {allsnps_path.name}")
    with gzip.open(allsnps_path, 'rt') as f:
        for i in range(6):
            line = f.readline().strip()
            print(f"Line {i}: {line[:120]}{'...' if len(line) > 120 else ''}")
else:
    print("allsnps.txt not found")

# 5. Check CS2G_INFO.md for coordinate system documentation
print("\n" + "=" * 60)
print("cS2G INFO FILE")
print("=" * 60)
info_path = Path(r'C:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs\data\external\cS2G\CS2G_INFO.md')
if info_path.exists():
    print(info_path.read_text()[:500])
else:
    print("Info file not found")

# 6. Test Platform API with correct field names
print("\n" + "=" * 60)
print("PLATFORM API TEST QUERY")
print("=" * 60)
test_query = """
query {
  studies(page: {size: 3, index: 0}) {
    count
    rows {
      id
      traitFromSource
      studyType
    }
  }
}
"""
r3 = requests.post('https://api.platform.opentargets.org/api/v4/graphql',
                   json={'query': test_query})
print("Response:", r3.json())
