#!/usr/bin/env python3
"""Check files in Zenodo deposit."""
import requests

TOKEN = "v0vwEqX8u9dw6MUFZqAQJSGjwcqA3JImFA5zQbPJx4MIJrhlfQgVp77jJz7p"
DEPOSIT_ID = 18165020

r = requests.get(
    f"https://zenodo.org/api/deposit/depositions/{DEPOSIT_ID}",
    params={"access_token": TOKEN}
)

if r.status_code != 200:
    print(f"Error: {r.status_code}")
    print(r.text)
else:
    d = r.json()
    files = d.get("files", [])
    print(f"Total files: {len(files)}")
    print("\nAll files:")
    for f in files:
        size = f.get("filesize", 0)
        size_mb = size / (1024 * 1024)
        print(f"  - {f['filename']} ({size_mb:.2f} MB)")
