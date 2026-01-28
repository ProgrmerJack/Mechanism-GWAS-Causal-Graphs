#!/usr/bin/env python3
"""
Check Zenodo record and deposit status.
"""

import os
import requests
import json


def main():
    print("=" * 60)
    print("ZENODO STATUS CHECK")
    print("=" * 60)
    
    # Check the old record
    record_id = "17880201"
    print(f"\n1. Checking record {record_id}...")
    r = requests.get(f"https://zenodo.org/api/records/{record_id}")
    print(f"   Status: {r.status_code}")
    if r.status_code == 410:
        print("   RESULT: Record has been DELETED")
    elif r.status_code == 200:
        data = r.json()
        print(f"   Title: {data.get('metadata', {}).get('title', 'N/A')}")
        print(f"   DOI: {data.get('doi', 'N/A')}")
    else:
        print(f"   Response: {r.text[:200]}")
    
    # Search for related records
    print("\n2. Searching for related records...")
    search_terms = ["Mechanism GWAS Causal", "gene-mechanism-causal-link"]
    for term in search_terms:
        r = requests.get("https://zenodo.org/api/records", params={"q": term, "size": 5})
        if r.status_code == 200:
            hits = r.json().get("hits", {}).get("hits", [])
            print(f"   Search '{term}': {len(hits)} results")
            for h in hits:
                rec_id = h.get("id")
                title = h.get("metadata", {}).get("title", "N/A")[:60]
                print(f"      - {rec_id}: {title}")
    
    # Check if token is available
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("\n3. ZENODO_TOKEN not set - cannot check deposits")
        print("   Set it with: $env:ZENODO_TOKEN = 'your_token'")
        return
    
    print("\n3. Checking your deposits...")
    headers = {"Authorization": f"Bearer {token}"}
    r = requests.get("https://zenodo.org/api/deposit/depositions", headers=headers)
    
    if r.status_code == 200:
        deposits = r.json()
        print(f"   Found {len(deposits)} deposit(s)")
        for d in deposits[:20]:
            dep_id = d.get("id")
            title = d.get("title", "Untitled")[:50]
            state = d.get("state", "unknown")
            submitted = d.get("submitted", False)
            doi = d.get("metadata", {}).get("prereserve_doi", {}).get("doi", "N/A")
            print(f"   - ID {dep_id}: '{title}' (state={state}, submitted={submitted})")
            print(f"     DOI: {doi}")
            
            # Check if this is related to the deleted record
            conceptrecid = d.get("conceptrecid")
            if conceptrecid:
                print(f"     Concept ID: {conceptrecid}")
    else:
        print(f"   Error: {r.status_code} - {r.text[:200]}")
    
    print("\n" + "=" * 60)
    print("RECOMMENDATION:")
    print("Since record 17880201 has been deleted, you need to either:")
    print("  1. Create a NEW deposit (fresh DOI)")
    print("  2. Use an EXISTING draft deposit if you have one")
    print("=" * 60)


if __name__ == "__main__":
    main()
