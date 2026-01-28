#!/usr/bin/env python3
"""
Comprehensive Zenodo Upload Script - Version 2

This script creates a new version (v2) of an existing Zenodo deposit and uploads
all data files and benchmarks. It handles large files, many files, and provides
robust error handling with progress tracking.

Target DOI: 10.5281/zenodo.17880201
New Version: 2

Usage:
    # Set your Zenodo access token (get from https://zenodo.org/account/settings/applications/)
    export ZENODO_TOKEN="your_token_here"
    
    # Run the upload
    python scripts/upload_zenodo_v2.py
    
    # Dry-run mode (just scan files, don't upload)
    python scripts/upload_zenodo_v2.py --dry-run
    
    # Upload specific folders only
    python scripts/upload_zenodo_v2.py --folders data benchmarks

Requirements:
    pip install requests tqdm

Author: Generated for Mechanism-GWAS-Causal-Graphs project
"""

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

try:
    import requests
except ImportError:
    print("ERROR: 'requests' package required. Install with: pip install requests")
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print("WARNING: 'tqdm' not installed. Install for progress bars: pip install tqdm")
    tqdm = None


class ZenodoUploader:
    """Handles uploading files to Zenodo with version management."""
    
    # Zenodo API endpoints
    ZENODO_API = "https://zenodo.org/api"
    SANDBOX_API = "https://sandbox.zenodo.org/api"
    
    def __init__(
        self,
        access_token: str,
        concept_doi: str = "10.5281/zenodo.17880201",
        sandbox: bool = False,
        chunk_size: int = 100 * 1024 * 1024,  # 100 MB chunks
    ):
        """
        Initialize Zenodo uploader.
        
        Parameters
        ----------
        access_token : str
            Zenodo personal access token
        concept_doi : str
            Concept DOI (version-independent) of the deposit
        sandbox : bool
            Use sandbox server for testing
        chunk_size : int
            Size of chunks for large file uploads (bytes)
        """
        self.access_token = access_token
        self.concept_doi = concept_doi
        self.base_url = self.SANDBOX_API if sandbox else self.ZENODO_API
        self.chunk_size = chunk_size
        self.session = requests.Session()
        self.session.headers.update({
            "Authorization": f"Bearer {self.access_token}",
        })
        
        # Extract record ID from DOI
        self.record_id = self.concept_doi.split(".")[-1]
        
        print(f"Initialized ZenodoUploader:")
        print(f"  API: {self.base_url}")
        print(f"  Concept DOI: {self.concept_doi}")
        print(f"  Record ID: {self.record_id}")
        print()
    
    def get_latest_version(self) -> Dict:
        """Get the latest version of the deposit."""
        url = f"{self.base_url}/records/{self.record_id}"
        
        print(f"Fetching latest version from: {url}")
        response = self.session.get(url)
        
        if response.status_code == 200:
            data = response.json()
            print(f"  ✓ Found deposit: {data.get('metadata', {}).get('title', 'N/A')}")
            print(f"  Version: {data.get('metadata', {}).get('version', 'N/A')}")
            return data
        else:
            print(f"  ✗ Error fetching deposit: {response.status_code}")
            print(f"  Response: {response.text}")
            raise Exception(f"Failed to fetch deposit: {response.status_code}")
    
    def create_new_version(self, latest_deposit: Dict) -> Dict:
        """Create a new version of the deposit."""
        # Get the deposit ID from the latest version
        deposit_id = latest_deposit.get("id")
        
        # Check if deposit is published (required to create new version)
        if latest_deposit.get("state") != "done":
            print("  ⚠ Deposit is not published. Publishing first...")
            publish_url = f"{self.base_url}/deposit/depositions/{deposit_id}/actions/publish"
            pub_response = self.session.post(publish_url)
            if pub_response.status_code not in [200, 201, 202]:
                raise Exception(f"Failed to publish deposit: {pub_response.status_code}")
            time.sleep(2)  # Wait for publish
        
        # Create new version
        newversion_url = f"{self.base_url}/deposit/depositions/{deposit_id}/actions/newversion"
        
        print(f"Creating new version...")
        print(f"  URL: {newversion_url}")
        
        response = self.session.post(newversion_url)
        
        if response.status_code in [200, 201]:
            data = response.json()
            # Get the draft URL
            draft_url = data.get("links", {}).get("latest_draft")
            if draft_url:
                # Fetch the draft
                draft_response = self.session.get(draft_url)
                if draft_response.status_code == 200:
                    draft_data = draft_response.json()
                    new_id = draft_data.get("id")
                    print(f"  ✓ Created new version draft with ID: {new_id}")
                    return draft_data
            
            raise Exception("Could not find draft URL in response")
        else:
            print(f"  ✗ Error creating new version: {response.status_code}")
            print(f"  Response: {response.text}")
            raise Exception(f"Failed to create new version: {response.status_code}")
    
    def update_metadata(self, deposit_id: str, version: str = "2") -> None:
        """Update deposit metadata with new version number."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}"
        
        # Get current metadata
        response = self.session.get(url)
        if response.status_code != 200:
            raise Exception(f"Failed to fetch deposit metadata: {response.status_code}")
        
        current_data = response.json()
        metadata = current_data.get("metadata", {})
        
        # Update version
        metadata["version"] = version
        
        # Update metadata
        print(f"Updating metadata with version {version}...")
        update_response = self.session.put(
            url,
            json={"metadata": metadata},
            headers={"Content-Type": "application/json"}
        )
        
        if update_response.status_code == 200:
            print(f"  ✓ Metadata updated")
        else:
            print(f"  ⚠ Warning: Metadata update returned {update_response.status_code}")
            print(f"  Response: {update_response.text}")
    
    def delete_existing_files(self, deposit_id: str) -> None:
        """Delete existing files from the draft."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}"
        response = self.session.get(url)
        
        if response.status_code != 200:
            print(f"  ⚠ Could not fetch existing files")
            return
        
        data = response.json()
        files = data.get("files", [])
        
        if not files:
            print("  No existing files to delete")
            return
        
        print(f"Deleting {len(files)} existing files...")
        for file_entry in files:
            file_id = file_entry.get("id")
            file_name = file_entry.get("filename")
            delete_url = f"{self.base_url}/deposit/depositions/{deposit_id}/files/{file_id}"
            
            del_response = self.session.delete(delete_url)
            if del_response.status_code in [200, 204]:
                print(f"  ✓ Deleted: {file_name}")
            else:
                print(f"  ✗ Failed to delete {file_name}: {del_response.status_code}")
    
    def compute_checksum(self, filepath: Path) -> str:
        """Compute MD5 checksum of a file."""
        md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                md5.update(chunk)
        return md5.hexdigest()
    
    def upload_file(
        self,
        deposit_id: str,
        filepath: Path,
        remote_path: str,
    ) -> bool:
        """
        Upload a single file to the deposit.
        
        Parameters
        ----------
        deposit_id : str
            Deposit ID
        filepath : Path
            Local file path
        remote_path : str
            Path in Zenodo (relative path for organization)
        
        Returns
        -------
        bool
            True if upload successful
        """
        file_size = filepath.stat().st_size
        file_size_mb = file_size / (1024 * 1024)
        
        # Create bucket URL
        bucket_url = f"{self.base_url}/deposit/depositions/{deposit_id}/files"
        
        # Prepare file upload
        filename = remote_path.replace("\\", "/")  # Ensure forward slashes
        
        print(f"  Uploading: {filename} ({file_size_mb:.2f} MB)")
        
        # Upload with progress bar
        if tqdm and file_size > 10 * 1024 * 1024:  # Use progress bar for files > 10MB
            with open(filepath, "rb") as f:
                with tqdm(total=file_size, unit='B', unit_scale=True, desc=filename) as pbar:
                    # Read file in chunks and update progress
                    data = b""
                    for chunk in iter(lambda: f.read(self.chunk_size), b""):
                        data += chunk
                        pbar.update(len(chunk))
                    
                    # Upload
                    files = {"file": (filename, data)}
                    response = self.session.post(bucket_url, files=files)
        else:
            # Small file - upload directly
            with open(filepath, "rb") as f:
                files = {"file": (filename, f)}
                response = self.session.post(bucket_url, files=files)
        
        if response.status_code in [200, 201]:
            print(f"  ✓ Uploaded successfully")
            return True
        else:
            print(f"  ✗ Upload failed: {response.status_code}")
            print(f"  Response: {response.text[:500]}")
            return False
    
    def publish_deposit(self, deposit_id: str) -> bool:
        """Publish the deposit."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}/actions/publish"
        
        print("Publishing deposit...")
        response = self.session.post(url)
        
        if response.status_code in [200, 201, 202]:
            data = response.json()
            doi = data.get("doi", "N/A")
            print(f"  ✓ Published successfully!")
            print(f"  DOI: {doi}")
            return True
        else:
            print(f"  ✗ Publishing failed: {response.status_code}")
            print(f"  Response: {response.text}")
            return False


def scan_files(root_dir: Path, folders: List[str]) -> List[Tuple[Path, str]]:
    """
    Scan folders and collect all files to upload.
    
    Parameters
    ----------
    root_dir : Path
        Project root directory
    folders : List[str]
        Folders to scan
    
    Returns
    -------
    List[Tuple[Path, str]]
        List of (filepath, remote_path) tuples
    """
    files_to_upload = []
    
    for folder in folders:
        folder_path = root_dir / folder
        
        if not folder_path.exists():
            print(f"  ⚠ Folder not found: {folder}")
            continue
        
        print(f"Scanning: {folder}/")
        
        # Recursively find all files
        all_files = list(folder_path.rglob("*"))
        file_count = sum(1 for f in all_files if f.is_file())
        
        print(f"  Found {file_count} files")
        
        for filepath in all_files:
            if filepath.is_file():
                # Create relative path for Zenodo
                relative_path = filepath.relative_to(root_dir)
                remote_path = str(relative_path).replace("\\", "/")
                
                files_to_upload.append((filepath, remote_path))
    
    return files_to_upload


def format_size(size_bytes: int) -> str:
    """Format bytes as human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"


def main():
    parser = argparse.ArgumentParser(
        description="Upload files to Zenodo (create new version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Upload everything (data, benchmarks, regulatorybench/benchmarks)
  python scripts/upload_zenodo_v2.py
  
  # Dry-run (scan files only, no upload)
  python scripts/upload_zenodo_v2.py --dry-run
  
  # Upload specific folders only
  python scripts/upload_zenodo_v2.py --folders data benchmarks
  
  # Use sandbox for testing
  python scripts/upload_zenodo_v2.py --sandbox --dry-run

Environment:
  ZENODO_TOKEN    Zenodo personal access token (required)
                  Get from: https://zenodo.org/account/settings/applications/
        """
    )
    
    parser.add_argument(
        "--token",
        type=str,
        default=os.environ.get("ZENODO_TOKEN"),
        help="Zenodo access token (or set ZENODO_TOKEN env var)"
    )
    
    parser.add_argument(
        "--doi",
        type=str,
        default="10.5281/zenodo.17880201",
        help="Concept DOI of the deposit (default: 10.5281/zenodo.17880201)"
    )
    
    parser.add_argument(
        "--version",
        type=str,
        default="2",
        help="Version number for the new version (default: 2)"
    )
    
    parser.add_argument(
        "--folders",
        nargs="+",
        default=["data", "benchmarks", "regulatorybench/benchmarks"],
        help="Folders to upload (default: data benchmarks regulatorybench/benchmarks)"
    )
    
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).parent.parent,
        help="Project root directory"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Scan files only, don't upload"
    )
    
    parser.add_argument(
        "--sandbox",
        action="store_true",
        help="Use Zenodo sandbox for testing"
    )
    
    parser.add_argument(
        "--skip-delete",
        action="store_true",
        help="Don't delete existing files in the new version draft"
    )
    
    args = parser.parse_args()
    
    # Check token
    if not args.token:
        print("ERROR: Zenodo access token required!")
        print()
        print("Set environment variable:")
        print("  export ZENODO_TOKEN='your_token_here'")
        print()
        print("Or use --token argument:")
        print("  python scripts/upload_zenodo_v2.py --token 'your_token_here'")
        print()
        print("Get your token from:")
        print("  https://zenodo.org/account/settings/applications/")
        sys.exit(1)
    
    root_dir = args.root.resolve()
    print("=" * 70)
    print("Zenodo Upload Script - Version 2")
    print("=" * 70)
    print(f"Project root: {root_dir}")
    print(f"Target DOI: {args.doi}")
    print(f"New version: {args.version}")
    print(f"Folders: {', '.join(args.folders)}")
    if args.dry_run:
        print("Mode: DRY RUN (no upload)")
    if args.sandbox:
        print("Using: SANDBOX (test server)")
    print()
    
    # Step 1: Scan files
    print("=" * 70)
    print("STEP 1: Scanning files")
    print("=" * 70)
    
    files_to_upload = scan_files(root_dir, args.folders)
    
    if not files_to_upload:
        print("No files found to upload!")
        sys.exit(1)
    
    # Calculate total size
    total_size = sum(f[0].stat().st_size for f in files_to_upload)
    
    print()
    print(f"Summary:")
    print(f"  Total files: {len(files_to_upload)}")
    print(f"  Total size: {format_size(total_size)}")
    print()
    
    if args.dry_run:
        print("Dry-run mode: Listing first 20 files...")
        for filepath, remote_path in files_to_upload[:20]:
            size = filepath.stat().st_size
            print(f"  {remote_path} ({format_size(size)})")
        if len(files_to_upload) > 20:
            print(f"  ... and {len(files_to_upload) - 20} more files")
        print()
        print("Dry-run complete. Use without --dry-run to upload.")
        sys.exit(0)
    
    # Step 2: Initialize Zenodo uploader
    print("=" * 70)
    print("STEP 2: Connecting to Zenodo")
    print("=" * 70)
    
    uploader = ZenodoUploader(
        access_token=args.token,
        concept_doi=args.doi,
        sandbox=args.sandbox,
    )
    
    # Step 3: Get latest version
    print()
    print("=" * 70)
    print("STEP 3: Fetching latest version")
    print("=" * 70)
    
    try:
        latest_deposit = uploader.get_latest_version()
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Step 4: Create new version
    print()
    print("=" * 70)
    print("STEP 4: Creating new version")
    print("=" * 70)
    
    try:
        new_draft = uploader.create_new_version(latest_deposit)
        deposit_id = new_draft.get("id")
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Step 5: Update metadata
    print()
    print("=" * 70)
    print("STEP 5: Updating metadata")
    print("=" * 70)
    
    uploader.update_metadata(deposit_id, args.version)
    
    # Step 6: Delete existing files (optional)
    if not args.skip_delete:
        print()
        print("=" * 70)
        print("STEP 6: Deleting existing files")
        print("=" * 70)
        
        uploader.delete_existing_files(deposit_id)
    else:
        print()
        print("=" * 70)
        print("STEP 6: Skipping file deletion (--skip-delete)")
        print("=" * 70)
    
    # Step 7: Upload files
    print()
    print("=" * 70)
    print("STEP 7: Uploading files")
    print("=" * 70)
    print(f"Uploading {len(files_to_upload)} files...")
    print()
    
    upload_success = 0
    upload_failed = 0
    
    for i, (filepath, remote_path) in enumerate(files_to_upload, 1):
        print(f"[{i}/{len(files_to_upload)}]", end=" ")
        
        try:
            success = uploader.upload_file(deposit_id, filepath, remote_path)
            if success:
                upload_success += 1
            else:
                upload_failed += 1
        except Exception as e:
            print(f"  ✗ Error: {e}")
            upload_failed += 1
        
        # Small delay to avoid rate limiting
        time.sleep(0.5)
    
    print()
    print(f"Upload complete:")
    print(f"  Success: {upload_success}")
    print(f"  Failed: {upload_failed}")
    print()
    
    if upload_failed > 0:
        print("⚠ Some uploads failed. Review errors above.")
        print("  You can re-run this script to retry failed uploads.")
        print()
    
    # Step 8: Publish (optional, ask user)
    print("=" * 70)
    print("STEP 8: Publish deposit")
    print("=" * 70)
    print()
    print("The new version draft is ready.")
    print("Review it at: https://zenodo.org/deposit/{0}".format(deposit_id))
    print()
    
    if upload_failed == 0:
        response = input("Publish now? (yes/no): ").strip().lower()
        
        if response in ["yes", "y"]:
            uploader.publish_deposit(deposit_id)
        else:
            print("Not publishing. You can publish manually from the Zenodo website.")
    else:
        print("Not publishing due to upload failures.")
        print("Fix errors and re-run, or publish manually after reviewing.")
    
    print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
