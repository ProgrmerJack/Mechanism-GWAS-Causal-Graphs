#!/usr/bin/env python3
"""
Zenodo Upload Script for Mechanism-GWAS-Causal-Graphs
Uploads the repository to Zenodo and returns a DOI
"""

import requests
import json
import os
import sys
import zipfile
import tempfile
from datetime import datetime

# Zenodo API configuration
ZENODO_API_TOKEN = "v0vwEqX8u9dw6MUFZqAQJSGjwcqA3JImFA5zQbPJx4MIJrhlfQgVp77jJz7p"
ZENODO_API_URL = "https://zenodo.org/api/deposit/depositions"

headers = {
    "Authorization": f"Bearer {ZENODO_API_TOKEN}",
    "Content-Type": "application/json"
}

def create_deposit(metadata):
    """Create a new Zenodo deposit"""
    response = requests.post(
        ZENODO_API_URL,
        headers=headers,
        json={"metadata": metadata}
    )
    if response.status_code == 201:
        return response.json()
    else:
        print(f"Error creating deposit: {response.status_code}")
        print(response.text)
        return None

def upload_file(deposit_id, bucket_url, filepath, filename=None):
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

def publish_deposit(deposit_id):
    """Publish the deposit to get a DOI"""
    response = requests.post(
        f"{ZENODO_API_URL}/{deposit_id}/actions/publish",
        headers={"Authorization": f"Bearer {ZENODO_API_TOKEN}"}
    )
    if response.status_code == 202:
        return response.json()
    else:
        print(f"Error publishing: {response.status_code}")
        print(response.text)
        return None

def create_archive(base_path, output_path):
    """Create a zip archive of the repository"""
    print("Creating archive...")
    
    # Files and directories to include
    include_patterns = [
        'src/',
        'workflow/',
        'manuscript/',
        'data/manifests/',
        'tests/',
        'environment.yml',
        'requirements.txt',
        'setup.py',
        'README.md',
        'LICENSE',
        'CITATION.cff',
        '.zenodo.json',
    ]
    
    # Files to exclude
    exclude_patterns = [
        '__pycache__',
        '.git',
        '.pytest_cache',
        '*.pyc',
        '.DS_Store',
    ]
    
    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(base_path):
            # Skip excluded directories
            dirs[:] = [d for d in dirs if d not in ['__pycache__', '.git', '.pytest_cache', 'node_modules']]
            
            for file in files:
                if file.endswith('.pyc') or file == '.DS_Store':
                    continue
                    
                filepath = os.path.join(root, file)
                arcname = os.path.relpath(filepath, base_path)
                
                # Include if matches any include pattern or is a top-level file
                for pattern in include_patterns:
                    if arcname.startswith(pattern) or arcname == pattern.rstrip('/'):
                        zipf.write(filepath, arcname)
                        break
    
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"  Archive created: {output_path} ({size_mb:.2f} MB)")
    return output_path

def main():
    print("\n" + "="*70)
    print("Mechanism-GWAS-Causal-Graphs: Zenodo Upload")
    print("="*70)
    
    # Metadata from .zenodo.json
    metadata = {
        "title": "Mechanism-First Causal Graphs for Noncoding GWAS: Code and Data",
        "description": """This repository contains the code, configuration, and reproducibility artifacts for the paper "Path-Probability Models Outperform Point-Estimate Scores for Noncoding GWAS Gene Prioritization".

The framework implements probabilistic mechanism graphs with path-probability inference for GWAS gene prioritization, demonstrating that explicit causal paths outperform point-estimate locus-to-gene scores.

Key features:
- SuSiE-based fine-mapping with multi-signal support
- coloc.susie colocalization for multi-causal variant handling
- ABC (Activity-by-Contact) and PCHi-C enhancer-gene linking
- Noisy-OR probabilistic inference model
- Per-module probability calibration (ECE < 0.05)
- eQTL Catalogue cross-study replication (78% replication rate)
- Three-tier benchmark system with anti-leakage provisions

Results:
- 76% recall at rank 20 vs 58% for Open Targets L2G
- Effect size correlation r=0.89 with eQTL Catalogue
- Calibrated probabilities at each module

Full pipeline reproducible via Snakemake workflow.""",
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
        "version": "1.0.0",
        "language": "eng"
    }
    
    # Create deposit
    print("\nCreating Zenodo deposit...")
    deposit = create_deposit(metadata)
    
    if not deposit:
        print("Failed to create deposit. Exiting.")
        return None
    
    deposit_id = deposit['id']
    bucket_url = deposit['links']['bucket']
    print(f"✓ Created deposit ID: {deposit_id}")
    
    # Base path
    base_path = os.path.dirname(os.path.abspath(__file__))
    
    # Upload key files individually
    print("\nUploading files...")
    
    files_to_upload = [
        ('README.md', None),
        ('LICENSE', None),
        ('CITATION.cff', None),
        ('.zenodo.json', None),
        ('environment.yml', None),
        ('requirements.txt', None),
        ('setup.py', None),
        ('manuscript/main.tex', 'manuscript_main.tex'),
        ('manuscript/references.bib', 'manuscript_references.bib'),
        ('workflow/Snakefile', 'workflow_Snakefile'),
    ]
    
    for rel_path, custom_name in files_to_upload:
        filepath = os.path.join(base_path, rel_path)
        if os.path.exists(filepath):
            upload_file(deposit_id, bucket_url, filepath, custom_name)
        else:
            print(f"  ! File not found: {rel_path}")
    
    # Upload source code as zip
    print("\nCreating source code archive...")
    with tempfile.TemporaryDirectory() as tmpdir:
        archive_path = os.path.join(tmpdir, "mechanism-gwas-source-v1.0.0.zip")
        
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
        
        size_mb = os.path.getsize(archive_path) / (1024 * 1024)
        print(f"  Archive size: {size_mb:.2f} MB")
        
        upload_file(deposit_id, bucket_url, archive_path, "mechanism-gwas-source-v1.0.0.zip")
    
    # Publish
    print("\nPublishing deposit...")
    published = publish_deposit(deposit_id)
    
    if published:
        doi = published.get('doi', 'N/A')
        doi_url = published.get('doi_url', 'N/A')
        
        print("\n" + "="*70)
        print("SUCCESS! Deposit published.")
        print("="*70)
        print(f"\n  DOI:     {doi}")
        print(f"  URL:     {doi_url}")
        print(f"  ID:      {deposit_id}")
        
        # Save results
        results = {
            'doi': doi,
            'doi_url': doi_url,
            'id': deposit_id,
            'upload_date': datetime.now().isoformat()
        }
        
        results_path = os.path.join(base_path, 'ZENODO_DOI.json')
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Results saved to: {results_path}")
        
        # Instructions for updating manuscript
        print("\n" + "-"*70)
        print("NEXT STEPS:")
        print("-"*70)
        print(f"1. Update manuscript/main.tex:")
        print(f"   Replace 'zenodo.XXXXXXX' with '{doi}'")
        print(f"\n2. Update CITATION.cff:")
        print(f"   Replace 'zenodo.XXXXXXX' with '{doi}'")
        print(f"\n3. Commit and push to GitHub")
        print("-"*70)
        
        return results
    else:
        print("\nFailed to publish. Deposit remains in draft state.")
        print(f"You can edit it at: https://zenodo.org/deposit/{deposit_id}")
        return {'id': deposit_id, 'status': 'draft'}

if __name__ == "__main__":
    main()
