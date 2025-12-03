#!/usr/bin/env python3
"""
Create Zenodo v1.2.0 - Clean version WITHOUT manuscript/submission files
This creates a new version that only contains code, data, and reproducibility materials
"""

import requests
import json
import os
import zipfile
import tempfile
from datetime import datetime

ZENODO_API_TOKEN = "v0vwEqX8u9dw6MUFZqAQJSGjwcqA3JImFA5zQbPJx4MIJrhlfQgVp77jJz7p"
ZENODO_API_URL = "https://zenodo.org/api/deposit/depositions"

# The concept DOI from v1.0.0 (links all versions)
EXISTING_DEPOSIT_ID = "17798899"

headers = {
    "Authorization": f"Bearer {ZENODO_API_TOKEN}",
}

def create_new_version(deposit_id):
    """Create a new version of an existing deposit"""
    response = requests.post(
        f"{ZENODO_API_URL}/{deposit_id}/actions/newversion",
        headers=headers
    )
    if response.status_code == 201:
        return response.json()
    else:
        print(f"Error creating new version: {response.status_code}")
        print(response.text)
        return None

def get_draft_deposit(deposit_id):
    """Get the draft deposit info"""
    response = requests.get(
        f"{ZENODO_API_URL}/{deposit_id}",
        headers=headers
    )
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error getting deposit: {response.status_code}")
        print(response.text)
        return None

def delete_file(deposit_id, file_id):
    """Delete a file from a draft deposit"""
    response = requests.delete(
        f"{ZENODO_API_URL}/{deposit_id}/files/{file_id}",
        headers=headers
    )
    if response.status_code == 204:
        return True
    else:
        print(f"Error deleting file: {response.status_code}")
        return False

def upload_file(bucket_url, filepath, filename=None):
    """Upload a file to a Zenodo deposit"""
    if filename is None:
        filename = os.path.basename(filepath)
    
    with open(filepath, 'rb') as f:
        data = f.read()
    
    response = requests.put(
        f"{bucket_url}/{filename}",
        headers={"Authorization": f"Bearer {ZENODO_API_TOKEN}"},
        data=data
    )
    
    if response.status_code in [200, 201]:
        print(f"  ✓ Uploaded: {filename}")
        return True
    else:
        print(f"  ✗ Failed to upload {filename}: {response.status_code}")
        print(response.text)
        return False

def update_metadata(deposit_id, metadata):
    """Update the metadata of a deposit"""
    response = requests.put(
        f"{ZENODO_API_URL}/{deposit_id}",
        headers={**headers, "Content-Type": "application/json"},
        json={"metadata": metadata}
    )
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error updating metadata: {response.status_code}")
        print(response.text)
        return None

def publish_deposit(deposit_id):
    """Publish the deposit"""
    response = requests.post(
        f"{ZENODO_API_URL}/{deposit_id}/actions/publish",
        headers=headers
    )
    if response.status_code == 202:
        return response.json()
    else:
        print(f"Error publishing: {response.status_code}")
        print(response.text)
        return None

def main():
    print("\n" + "="*70)
    print("Creating Zenodo v1.2.0 - Clean Version (No Manuscript)")
    print("="*70)
    
    base_path = os.path.dirname(os.path.abspath(__file__))
    
    # Step 1: Create new version
    print("\n[1/6] Creating new version from existing deposit...")
    new_version = create_new_version(EXISTING_DEPOSIT_ID)
    
    if not new_version:
        print("Failed to create new version. Exiting.")
        return None
    
    # Get the draft deposit URL
    draft_url = new_version.get('links', {}).get('latest_draft')
    if draft_url:
        draft_id = draft_url.split('/')[-1]
    else:
        print("Could not find draft URL. Exiting.")
        return None
    
    print(f"  ✓ New draft created: {draft_id}")
    
    # Step 2: Get draft details
    print("\n[2/6] Getting draft deposit details...")
    draft = get_draft_deposit(draft_id)
    
    if not draft:
        print("Failed to get draft. Exiting.")
        return None
    
    bucket_url = draft['links']['bucket']
    existing_files = draft.get('files', [])
    
    print(f"  Found {len(existing_files)} existing files")
    
    # Step 3: Delete manuscript-related files
    print("\n[3/6] Removing manuscript and submission files...")
    
    files_to_delete = [
        'manuscript_main.tex',
        'manuscript_main_v2.tex',
        'manuscript_references.bib',
        'submission_cover_letter.tex',
        'submission_reporting_summary.md'
    ]
    
    for f in existing_files:
        filename = f.get('filename', '')
        if any(d in filename for d in files_to_delete) or 'manuscript' in filename.lower() or 'submission' in filename.lower():
            print(f"  Deleting: {filename}")
            delete_file(draft_id, f['id'])
    
    # Step 4: Update metadata
    print("\n[4/6] Updating metadata for v1.2.0...")
    
    metadata = {
        "title": "Mechanism-First Causal Graphs for Noncoding GWAS: Code and Data",
        "description": """This repository contains the **code, configuration, and reproducibility artifacts** for the paper "Path-Probability Models Outperform Point-Estimate Scores for Noncoding GWAS Gene Prioritization".

**Note:** This version (v1.2.0) contains ONLY code and data - no manuscript files. The manuscript is available separately for journal review.

## Key Features
- SuSiE-based fine-mapping with multi-signal support
- coloc.susie colocalization for multi-causal variant handling
- ABC (Activity-by-Contact) and PCHi-C enhancer-gene linking
- Noisy-OR probabilistic inference model
- Per-module probability calibration (ECE < 0.05)
- eQTL Catalogue cross-study replication (78% replication rate)
- Three-tier benchmark system with anti-leakage provisions

## Results
- 76% recall at rank 20 vs 58% for Open Targets L2G
- Effect size correlation r=0.89 with eQTL Catalogue
- Calibrated probabilities at each module

## Reproduction
Full pipeline reproducible via Snakemake workflow. See REPRODUCE.md and START_HERE.md for instructions.""",
        "upload_type": "software",
        "creators": [
            {
                "name": "Ashuraliyev, Abduxoliq",
                "affiliation": "Independent Researcher, Tashkent, Uzbekistan"
            }
        ],
        "keywords": [
            "GWAS", "fine-mapping", "colocalization", "causal graphs",
            "gene prioritization", "eQTL", "regulatory genomics",
            "probabilistic inference", "mechanism graphs", "SuSiE",
            "ABC model", "PCHi-C", "coloc.susie", "path-probability"
        ],
        "license": "MIT",
        "access_right": "open",
        "related_identifiers": [
            {
                "identifier": "https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs",
                "relation": "isSupplementTo",
                "scheme": "url"
            }
        ],
        "version": "1.2.0",
        "language": "eng",
        "notes": "v1.2.0: Clean version containing code and reproducibility materials only. Manuscript removed for journal review process."
    }
    
    update_result = update_metadata(draft_id, metadata)
    if update_result:
        print("  ✓ Metadata updated")
    
    # Step 5: Upload updated files
    print("\n[5/6] Uploading clean files...")
    
    # Files to upload (NO manuscript/submission)
    files_to_upload = [
        'README.md',
        'LICENSE',
        'CITATION.cff',
        '.zenodo.json',
        'environment.yml',
        'REPRODUCE.md',
        'START_HERE.md',
    ]
    
    for rel_path in files_to_upload:
        filepath = os.path.join(base_path, rel_path)
        if os.path.exists(filepath):
            upload_file(bucket_url, filepath)
        else:
            print(f"  ! File not found: {rel_path}")
    
    # Upload source code archive (no manuscript)
    print("\n  Creating clean source archive...")
    with tempfile.TemporaryDirectory() as tmpdir:
        archive_path = os.path.join(tmpdir, "mechanism-gwas-source-v1.2.0.zip")
        
        with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Add src/ directory
            src_path = os.path.join(base_path, 'src')
            if os.path.exists(src_path):
                for root, dirs, files in os.walk(src_path):
                    dirs[:] = [d for d in dirs if d != '__pycache__']
                    for file in files:
                        if not file.endswith('.pyc'):
                            filepath = os.path.join(root, file)
                            arcname = os.path.relpath(filepath, base_path)
                            zipf.write(filepath, arcname)
            
            # Add workflow/ directory
            workflow_path = os.path.join(base_path, 'workflow')
            if os.path.exists(workflow_path):
                for root, dirs, files in os.walk(workflow_path):
                    for file in files:
                        filepath = os.path.join(root, file)
                        arcname = os.path.relpath(filepath, base_path)
                        zipf.write(filepath, arcname)
            
            # Add config/ directory
            config_path = os.path.join(base_path, 'config')
            if os.path.exists(config_path):
                for root, dirs, files in os.walk(config_path):
                    for file in files:
                        filepath = os.path.join(root, file)
                        arcname = os.path.relpath(filepath, base_path)
                        zipf.write(filepath, arcname)
            
            # Add data/manifests/ directory
            manifests_path = os.path.join(base_path, 'data', 'manifests')
            if os.path.exists(manifests_path):
                for file in os.listdir(manifests_path):
                    filepath = os.path.join(manifests_path, file)
                    if os.path.isfile(filepath):
                        arcname = os.path.join('data', 'manifests', file)
                        zipf.write(filepath, arcname)
            
            # Add tests/ directory
            tests_path = os.path.join(base_path, 'tests')
            if os.path.exists(tests_path):
                for root, dirs, files in os.walk(tests_path):
                    dirs[:] = [d for d in dirs if d != '__pycache__']
                    for file in files:
                        if not file.endswith('.pyc'):
                            filepath = os.path.join(root, file)
                            arcname = os.path.relpath(filepath, base_path)
                            zipf.write(filepath, arcname)
            
            # Add scripts/ directory
            scripts_path = os.path.join(base_path, 'scripts')
            if os.path.exists(scripts_path):
                for root, dirs, files in os.walk(scripts_path):
                    dirs[:] = [d for d in dirs if d != '__pycache__']
                    for file in files:
                        if not file.endswith('.pyc'):
                            filepath = os.path.join(root, file)
                            arcname = os.path.relpath(filepath, base_path)
                            zipf.write(filepath, arcname)
        
        size_mb = os.path.getsize(archive_path) / (1024 * 1024)
        print(f"    Archive size: {size_mb:.2f} MB")
        
        upload_file(bucket_url, archive_path, "mechanism-gwas-source-v1.2.0.zip")
    
    # Step 6: Publish
    print("\n[6/6] Publishing v1.2.0...")
    published = publish_deposit(draft_id)
    
    if published:
        doi = published.get('doi', 'N/A')
        doi_url = published.get('doi_url', 'N/A')
        
        print("\n" + "="*70)
        print("SUCCESS! Clean version v1.2.0 published.")
        print("="*70)
        print(f"\n  DOI:     {doi}")
        print(f"  URL:     {doi_url}")
        print(f"  ID:      {draft_id}")
        
        # Save results
        results = {
            'doi': doi,
            'doi_url': doi_url,
            'id': draft_id,
            'version': '1.2.0',
            'upload_date': datetime.now().isoformat(),
            'note': 'Clean version - no manuscript/submission files'
        }
        
        results_path = os.path.join(base_path, 'ZENODO_DOI.json')
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Results saved to: {results_path}")
        
        print("\n" + "-"*70)
        print("NOTE: Previous versions (v1.0.0, v1.1.0) still exist and contain")
        print("manuscript files. To fully remove them, contact Zenodo support.")
        print("-"*70)
        
        return results
    else:
        print("\nFailed to publish. Draft remains at:")
        print(f"  https://zenodo.org/deposit/{draft_id}")
        return {'id': draft_id, 'status': 'draft'}

if __name__ == "__main__":
    main()
