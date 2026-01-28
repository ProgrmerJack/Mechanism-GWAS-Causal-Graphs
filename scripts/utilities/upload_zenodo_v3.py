#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zenodo Upload Script v3 - Production-Ready Large Dataset Upload

This script properly handles Zenodo's hard limits:
  - 50 GB per record (can request up to 200 GB quota increase)
  - 100 files per record maximum

SOLUTION: Creates ZIP archives of large file collections to stay under 100 files.

Key features:
  1. Automatically creates ZIP archives for directories with many files
  2. Uses Zenodo's new bucket API for large files (50GB per file supported)
  3. Proper rate limiting with exponential backoff
  4. Resume capability - tracks uploaded files
  5. Handles connection timeouts gracefully
  6. Streaming ZIP creation (memory efficient)

Target DOI: 10.5281/zenodo.17880201

Usage:
    # Set your Zenodo access token
    export ZENODO_TOKEN="your_token_here"
    
    # Step 1: Create ZIP archives (do this first!)
    python scripts/upload_zenodo_v3.py --create-zips
    
    # Step 2: Upload (after ZIPs are created)
    python scripts/upload_zenodo_v3.py --upload
    
    # Or do both:
    python scripts/upload_zenodo_v3.py --create-zips --upload

Author: Generated for Mechanism-GWAS-Causal-Graphs project
Date: 2025
"""

import argparse
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time
import zipfile
import shutil
import io
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from datetime import datetime
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Fix Windows console encoding for Unicode characters
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except ImportError:
    print("ERROR: 'requests' package required. Install with: pip install requests")
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print("WARNING: 'tqdm' not installed. Progress bars disabled.")
    tqdm = None


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class ZenodoConfig:
    """Zenodo configuration and limits."""
    MAX_FILES_PER_RECORD: int = 100
    MAX_SIZE_PER_RECORD_GB: int = 50  # Default, can request 200GB
    RATE_LIMIT_REQUESTS_PER_MINUTE: int = 100
    RATE_LIMIT_REQUESTS_PER_HOUR: int = 5000
    
    # Safety margins
    MAX_FILES_SAFE: int = 80  # Leave room for additional files
    MAX_SIZE_PER_RECORD_SAFE_GB: int = 45  # Leave margin for safety
    
    # Retry configuration
    MAX_RETRIES: int = 5
    INITIAL_BACKOFF_SECONDS: float = 2.0
    MAX_BACKOFF_SECONDS: float = 120.0
    
    # Connection configuration
    CHUNK_SIZE: int = 10 * 1024 * 1024  # 10 MB chunks
    LARGE_FILE_CHUNK_SIZE: int = 50 * 1024 * 1024  # 50 MB chunks for files > 1GB
    LARGE_FILE_THRESHOLD: int = 1 * 1024 * 1024 * 1024  # 1 GB
    
    # Timeout configuration (connect_timeout, read_timeout)
    CONNECT_TIMEOUT: int = 30  # 30 seconds to establish connection
    READ_TIMEOUT: int = 3600  # 1 hour for reading response (large uploads)
    LARGE_FILE_READ_TIMEOUT: Optional[int] = None  # No timeout for files > 1GB
    
    # Rate limiting
    MIN_DELAY_BETWEEN_REQUESTS: float = 0.7  # ~85 req/min (under 100)


@dataclass
class FileEntry:
    """Represents a file to upload."""
    local_path: Path
    remote_name: str
    size_bytes: int
    is_zip: bool = False
    source_dir: Optional[str] = None


@dataclass
class UploadProgress:
    """Tracks upload progress for resume capability."""
    uploaded_files: Set[str] = field(default_factory=set)
    failed_files: Set[str] = field(default_factory=set)
    total_bytes_uploaded: int = 0
    start_time: Optional[datetime] = None
    deposit_id: Optional[str] = None  # Store deposit ID for proper resume
    bucket_url: Optional[str] = None  # Store bucket URL for direct upload
    
    def save(self, path: Path):
        """Save progress to file."""
        data = {
            "uploaded_files": list(self.uploaded_files),
            "failed_files": list(self.failed_files),
            "total_bytes_uploaded": self.total_bytes_uploaded,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "last_update": datetime.now().isoformat(),
            "deposit_id": self.deposit_id,
            "bucket_url": self.bucket_url
        }
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
    
    @classmethod
    def load(cls, path: Path) -> 'UploadProgress':
        """Load progress from file."""
        if not path.exists():
            return cls()
        
        try:
            with open(path, 'r') as f:
                data = json.load(f)
            
            progress = cls()
            progress.uploaded_files = set(data.get("uploaded_files", []))
            progress.failed_files = set(data.get("failed_files", []))
            progress.total_bytes_uploaded = data.get("total_bytes_uploaded", 0)
            progress.deposit_id = data.get("deposit_id")  # Load saved deposit ID
            progress.bucket_url = data.get("bucket_url")  # Load saved bucket URL
            
            start_str = data.get("start_time")
            if start_str:
                progress.start_time = datetime.fromisoformat(start_str)
            
            return progress
        except Exception as e:
            print(f"Warning: Could not load progress file: {e}")
            return cls()


# =============================================================================
# ZIP ARCHIVE CREATION
# =============================================================================

class ZipArchiveManager:
    """Manages creation of ZIP archives for large directories."""
    
    def __init__(self, root_dir: Path, output_dir: Path, config: ZenodoConfig):
        self.root_dir = root_dir
        self.output_dir = output_dir
        self.config = config
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def should_zip_directory(self, dir_path: Path, file_count: int) -> bool:
        """Determine if a directory should be zipped."""
        # Zip if there are many files
        return file_count > 10
    
    def create_archive_for_directory(
        self,
        source_dir: Path,
        archive_name: str,
        compression: int = zipfile.ZIP_DEFLATED,
        compresslevel: int = 6  # Balanced speed/size
    ) -> Optional[Path]:
        """
        Create a ZIP archive for a directory.
        
        Parameters
        ----------
        source_dir : Path
            Directory to archive
        archive_name : str
            Name of the archive file (without .zip extension)
        compression : int
            ZIP compression method
        compresslevel : int
            Compression level (1-9)
        
        Returns
        -------
        Optional[Path]
            Path to created archive, or None if failed
        """
        archive_path = self.output_dir / f"{archive_name}.zip"
        
        # Skip if already exists and recent
        if archive_path.exists():
            # Check if source is newer
            archive_mtime = archive_path.stat().st_mtime
            source_mtime = source_dir.stat().st_mtime
            
            # Also check for any newer file in source
            any_newer = False
            for f in source_dir.rglob("*"):
                if f.is_file() and f.stat().st_mtime > archive_mtime:
                    any_newer = True
                    break
            
            if not any_newer:
                print(f"  ⏭ Archive exists and is up-to-date: {archive_path.name}")
                return archive_path
            else:
                print(f"  ↻ Source modified, recreating: {archive_path.name}")
        
        print(f"  📦 Creating: {archive_path.name}")
        
        # Get all files to add
        all_files = list(source_dir.rglob("*"))
        files_to_add = [f for f in all_files if f.is_file()]
        
        if not files_to_add:
            print(f"    No files to archive in {source_dir.name}")
            return None
        
        total_size = sum(f.stat().st_size for f in files_to_add)
        
        # Create archive with progress
        temp_path = archive_path.with_suffix('.zip.tmp')
        
        try:
            with zipfile.ZipFile(
                temp_path, 'w',
                compression=compression,
                compresslevel=compresslevel
            ) as zf:
                if tqdm:
                    pbar = tqdm(
                        files_to_add,
                        desc=f"    {archive_name}",
                        unit="file"
                    )
                else:
                    pbar = files_to_add
                    print(f"    Adding {len(files_to_add)} files...")
                
                for file_path in pbar:
                    # Calculate relative path within archive
                    arcname = file_path.relative_to(source_dir)
                    zf.write(file_path, arcname)
            
            # Rename temp to final
            shutil.move(temp_path, archive_path)
            
            final_size = archive_path.stat().st_size
            compression_ratio = (1 - final_size / total_size) * 100 if total_size > 0 else 0
            
            print(f"    ✓ Created: {format_size(final_size)} ({compression_ratio:.1f}% compression)")
            return archive_path
            
        except Exception as e:
            print(f"    ✗ Error creating archive: {e}")
            if temp_path.exists():
                temp_path.unlink()
            return None
    
    def create_all_archives(self, folders_config: Dict[str, dict]) -> List[FileEntry]:
        """
        Create ZIP archives for all configured folders.
        
        Parameters
        ----------
        folders_config : Dict[str, dict]
            Configuration for each folder. Keys are relative paths,
            values contain:
            - 'zip' (bool): Create ZIP of entire directory
            - 'split_zip' (bool): Create separate ZIPs per subdirectory
            - 'flat_only' (bool): Upload only root-level files, skip subdirs
            - 'skip_subdirs' (list): Names of subdirectories to skip
        
        Returns
        -------
        List[FileEntry]
            List of files to upload (ZIPs and individual files)
        """
        files_to_upload = []
        
        for rel_path, config in folders_config.items():
            full_path = self.root_dir / rel_path
            
            if not full_path.exists():
                print(f"⚠ Folder not found: {rel_path}")
                continue
            
            print(f"\n📁 Processing: {rel_path}")
            
            # NEW: Handle flat_only mode - upload only root-level files
            if config.get('flat_only', False):
                skip_subdirs = set(config.get('skip_subdirs', []))
                file_count = 0
                total_size = 0
                
                for item in sorted(full_path.iterdir()):
                    # Skip directories or skip specified subdirectories
                    if item.is_dir():
                        if item.name in skip_subdirs:
                            print(f"    ⏭ Skipping extracted dir: {item.name}")
                        continue
                    
                    # It's a file - add it
                    files_to_upload.append(FileEntry(
                        local_path=item,
                        remote_name=f"{rel_path}/{item.name}".replace('\\', '/'),
                        size_bytes=item.stat().st_size,
                        is_zip=False,
                        source_dir=rel_path
                    ))
                    file_count += 1
                    total_size += item.stat().st_size
                
                print(f"    📄 Added {file_count} root-level files ({format_size(total_size)})")
            
            elif config.get('zip', False):
                # Create a single ZIP for this folder
                archive_name = rel_path.replace('/', '_').replace('\\', '_')
                archive_path = self.create_archive_for_directory(
                    full_path, archive_name
                )
                
                if archive_path:
                    files_to_upload.append(FileEntry(
                        local_path=archive_path,
                        remote_name=archive_path.name,
                        size_bytes=archive_path.stat().st_size,
                        is_zip=True,
                        source_dir=rel_path
                    ))
            
            elif config.get('split_zip', False):
                # Create separate ZIPs for each subdirectory
                for subdir in sorted(full_path.iterdir()):
                    if subdir.is_dir():
                        # Count files in this subdir
                        file_count = sum(1 for _ in subdir.rglob("*") if _.is_file())
                        
                        if file_count > 0:
                            archive_name = f"{rel_path.replace('/', '_')}_{subdir.name}"
                            archive_path = self.create_archive_for_directory(
                                subdir, archive_name
                            )
                            
                            if archive_path:
                                files_to_upload.append(FileEntry(
                                    local_path=archive_path,
                                    remote_name=archive_path.name,
                                    size_bytes=archive_path.stat().st_size,
                                    is_zip=True,
                                    source_dir=f"{rel_path}/{subdir.name}"
                                ))
                    elif subdir.is_file():
                        # Individual files at root of folder
                        files_to_upload.append(FileEntry(
                            local_path=subdir,
                            remote_name=f"{rel_path}/{subdir.name}".replace('\\', '/'),
                            size_bytes=subdir.stat().st_size,
                            is_zip=False,
                            source_dir=rel_path
                        ))
            
            else:
                # Add files individually (for small folders)
                skip_subdirs = set(config.get('skip_subdirs', []))
                
                for file_path in full_path.rglob("*"):
                    if file_path.is_file():
                        # Check if file is in a skipped subdirectory
                        try:
                            rel_to_folder = file_path.relative_to(full_path)
                            parts = rel_to_folder.parts
                            if parts and parts[0] in skip_subdirs:
                                continue  # Skip this file
                        except ValueError:
                            pass
                        
                        rel_file_path = file_path.relative_to(self.root_dir)
                        files_to_upload.append(FileEntry(
                            local_path=file_path,
                            remote_name=str(rel_file_path).replace('\\', '/'),
                            size_bytes=file_path.stat().st_size,
                            is_zip=False,
                            source_dir=rel_path
                        ))
        
        return files_to_upload


# =============================================================================
# ZENODO API CLIENT
# =============================================================================

class ZenodoClient:
    """Zenodo API client with robust error handling."""
    
    ZENODO_API = "https://zenodo.org/api"
    SANDBOX_API = "https://sandbox.zenodo.org/api"
    
    def __init__(
        self,
        access_token: str,
        sandbox: bool = False,
        config: ZenodoConfig = None
    ):
        self.access_token = access_token
        self.base_url = self.SANDBOX_API if sandbox else self.ZENODO_API
        self.config = config or ZenodoConfig()
        
        # Create session with retry logic
        self.session = requests.Session()
        self.session.headers.update({
            "Authorization": f"Bearer {self.access_token}",
        })
        
        # Configure retries with connection pooling optimized for large uploads
        # CRITICAL: Adapter retries are DISABLED for PUT requests to prevent file corruption
        # See urllib3 issue #459: "urllib3 fails to seek file object to beginning when it retries"
        # When adapter retries a file upload, the file generator/object is NOT reset,
        # causing subsequent retries to send wrong bytes or incomplete data.
        # Solution: Disable adapter-level retries for PUT, handle all retries at application level
        retry_strategy = Retry(
            total=0,  # DISABLE total retries (we handle retries manually at application level)
            connect=3,  # Only retry connection establishment (before file is sent)
            read=0,  # Don't retry read errors (would corrupt upload)
            status_forcelist=[502, 503, 504],  # Only for GET/HEAD requests
            allowed_methods=["HEAD", "GET", "OPTIONS"],  # Explicitly EXCLUDE PUT/POST from auto-retry
            backoff_factor=1
        )
        # Optimized connection pool for large file uploads
        # - pool_connections: Keep multiple persistent connections
        # - pool_maxsize: Allow many concurrent connections
        # - max_retries: Apply retry strategy
        # - pool_block: Don't block waiting for connection from pool
        adapter = HTTPAdapter(
            pool_connections=10,
            pool_maxsize=20,
            max_retries=retry_strategy,
            pool_block=False
        )
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        
        # Rate limiting state
        self._last_request_time = 0
        self._request_lock = threading.Lock()
    
    def _rate_limit_wait(self):
        """Wait to respect rate limits."""
        with self._request_lock:
            now = time.time()
            elapsed = now - self._last_request_time
            
            if elapsed < self.config.MIN_DELAY_BETWEEN_REQUESTS:
                wait_time = self.config.MIN_DELAY_BETWEEN_REQUESTS - elapsed
                time.sleep(wait_time)
            
            self._last_request_time = time.time()
    
    def _handle_rate_limit_response(self, response: requests.Response) -> float:
        """
        Handle rate limit headers and return wait time.
        
        Returns
        -------
        float
            Time to wait before retrying (0 if no wait needed)
        """
        if response.status_code != 429:
            return 0
        
        # Check headers
        remaining = response.headers.get('X-RateLimit-Remaining')
        reset_time = response.headers.get('X-RateLimit-Reset')
        
        if reset_time:
            try:
                reset_timestamp = int(reset_time)
                wait_time = max(0, reset_timestamp - time.time()) + 1
                return min(wait_time, self.config.MAX_BACKOFF_SECONDS)
            except ValueError:
                pass
        
        # Default backoff
        return 60  # Wait 1 minute
    
    def _request_with_retry(
        self,
        method: str,
        url: str,
        max_retries: int = None,
        **kwargs
    ) -> requests.Response:
        """
        Make a request with retry logic and rate limiting.
        
        Parameters
        ----------
        method : str
            HTTP method
        url : str
            Request URL
        max_retries : int
            Maximum retry attempts
        **kwargs
            Additional arguments for requests
        
        Returns
        -------
        requests.Response
        """
        max_retries = max_retries or self.config.MAX_RETRIES
        backoff = self.config.INITIAL_BACKOFF_SECONDS
        
        for attempt in range(max_retries + 1):
            self._rate_limit_wait()
            
            try:
                response = self.session.request(method, url, **kwargs)
                
                # Handle rate limiting
                if response.status_code == 429:
                    wait_time = self._handle_rate_limit_response(response)
                    if attempt < max_retries:
                        print(f"  ⚠ Rate limited. Waiting {wait_time:.1f}s...")
                        time.sleep(wait_time)
                        continue
                
                return response
                
            except requests.exceptions.ConnectionError as e:
                if attempt < max_retries:
                    print(f"  ⚠ Connection error (attempt {attempt + 1}): {e}")
                    time.sleep(backoff)
                    backoff = min(backoff * 2, self.config.MAX_BACKOFF_SECONDS)
                    continue
                raise
            
            except requests.exceptions.Timeout as e:
                if attempt < max_retries:
                    print(f"  ⚠ Timeout (attempt {attempt + 1}): {e}")
                    time.sleep(backoff)
                    backoff = min(backoff * 2, self.config.MAX_BACKOFF_SECONDS)
                    continue
                raise
        
        raise Exception(f"Max retries exceeded for {url}")
    
    def get_deposit(self, deposit_id: str) -> Dict:
        """Get deposit information."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}"
        response = self._request_with_retry("GET", url)
        
        if response.status_code == 200:
            return response.json()
        
        raise Exception(f"Failed to get deposit: {response.status_code} - {response.text}")
    
    def get_record(self, record_id: str) -> Dict:
        """Get published record information."""
        url = f"{self.base_url}/records/{record_id}"
        response = self._request_with_retry("GET", url)
        
        if response.status_code == 200:
            return response.json()
        
        raise Exception(f"Failed to get record: {response.status_code} - {response.text}")
    
    def create_new_version(self, deposit_id: str) -> Dict:
        """Create a new version of a deposit."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}/actions/newversion"
        
        response = self._request_with_retry("POST", url)
        
        if response.status_code in [200, 201]:
            data = response.json()
            
            # Get the draft URL
            draft_url = data.get("links", {}).get("latest_draft")
            if draft_url:
                # Fetch the draft
                draft_response = self._request_with_retry("GET", draft_url)
                if draft_response.status_code == 200:
                    return draft_response.json()
            
            raise Exception("Could not find draft URL in response")
        
        raise Exception(f"Failed to create new version: {response.status_code} - {response.text}")
    
    def get_or_create_draft(self, record_id: str) -> Dict:
        """
        Get existing draft or create a new version.
        
        Parameters
        ----------
        record_id : str
            Published record ID
        
        Returns
        -------
        Dict
            Draft deposit data
        """
        # First try to get the record
        record = self.get_record(record_id)
        
        # Check for existing draft
        draft_url = record.get("links", {}).get("draft")
        if draft_url:
            response = self._request_with_retry("GET", draft_url)
            if response.status_code == 200:
                print("  ✓ Found existing draft")
                return response.json()
        
        # Create new version
        print("  Creating new version...")
        return self.create_new_version(record_id)
    
    def update_metadata(self, deposit_id: str, metadata: Dict) -> None:
        """Update deposit metadata."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}"
        
        response = self._request_with_retry(
            "PUT", url,
            json={"metadata": metadata},
            headers={"Content-Type": "application/json"}
        )
        
        if response.status_code not in [200, 201]:
            print(f"  ⚠ Metadata update warning: {response.status_code}")
    
    def list_files(self, deposit_id: str) -> List[Dict]:
        """List files in a deposit."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}"
        response = self._request_with_retry("GET", url)
        
        if response.status_code == 200:
            return response.json().get("files", [])
        
        return []
    
    def delete_file(self, deposit_id: str, file_id: str) -> bool:
        """Delete a file from a deposit.
        
        Returns False if deletion fails (including 500 errors from Zenodo).
        """
        url = f"{self.base_url}/deposit/depositions/{deposit_id}/files/{file_id}"
        try:
            response = self._request_with_retry("DELETE", url)
            return response.status_code in [200, 204]
        except Exception as e:
            # 500 errors are expected for new version drafts
            if "500" in str(e):
                return False
            raise
    
    def delete_all_files(self, deposit_id: str) -> int:
        """Delete all files from a deposit.
        
        Note: New version drafts typically start empty, so deletion may fail
        with 500 errors. This is expected behavior.
        """
        try:
            files = self.list_files(deposit_id)
        except Exception as e:
            print(f"    ⚠ Could not list files (this is OK for new versions): {e}")
            return 0
        
        if not files:
            print("    ℹ No files to delete (draft is clean)")
            return 0
        
        deleted = 0
        
        for file_info in files:
            file_id = file_info.get("id")
            filename = file_info.get("filename", "unknown")
            
            try:
                if self.delete_file(deposit_id, file_id):
                    print(f"    ✓ Deleted: {filename}")
                    deleted += 1
                else:
                    print(f"    ✗ Failed to delete: {filename}")
            except Exception as e:
                print(f"    ⚠ Error deleting {filename}: {e}")
                # Continue with other files
        
        return deleted
    
    @staticmethod
    def _is_curl_available() -> bool:
        """
        Check if curl is available on the system.
        
        Returns
        -------
        bool
            True if curl is available and executable
        """
        try:
            result = subprocess.run(
                ["curl", "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return False
    
    def upload_file_curl(
        self,
        bucket_url: str,
        file_path: Path,
        remote_name: str,
        progress_callback=None
    ) -> bool:
        """
        Upload a file using curl subprocess (robust for large files on Windows).
        
        This method is a PROVEN workaround for Windows socket issues (error 10053/10054)
        that plague Python's requests library with very large file uploads.
        
        Based on jhpoelen/zenodo-upload which successfully uploads 9GB+ files:
        https://github.com/jhpoelen/zenodo-upload
        
        curl uses native OS networking with optimized TCP buffer management,
        bypassing all Python/urllib3 networking issues.
        
        Parameters
        ----------
        bucket_url : str
            Bucket URL from deposit links
        file_path : Path
            Local file path
        remote_name : str
            Name for the file in Zenodo
        progress_callback : callable
            Optional callback for progress updates (limited support with curl)
        
        Returns
        -------
        bool
            True if upload successful
        """
        # Construct the upload URL exactly as Zenodo expects
        # Format: {bucket_url}/{filename}?access_token={token}
        #
        # CRITICAL: Zenodo does NOT support folder/directory structure at all.
        # Looking at existing files in deposit 18165020, they use FLAT FILENAMES:
        #   - data_processed_mechanism_graphs_TCF7L2_10q25_T2D.json
        #   - data_processed_benchmark_tier1_gold_standard_genes.tsv
        #
        # The paths were FLATTENED by replacing "/" with "_".
        # This is the approach we MUST use to match the existing convention.
        flat_name = remote_name.replace("/", "_")
        upload_url = f"{bucket_url}/{flat_name}?access_token={self.access_token}"
        
        file_size = file_path.stat().st_size
        file_size_str = format_size(file_size)
        
        print(f"    🌐 Using curl for upload (bypasses Python networking issues)")
        print(f"      File: {file_path.name}")
        print(f"      Size: {file_size_str}")
        print(f"      Original path: {remote_name}")
        print(f"      Flat name: {flat_name}")
        
        # Build curl command based on jhpoelen/zenodo-upload
        # Key flags:
        #   --upload-file: Stream file directly (efficient for large files)
        #   --progress-bar: Show progress
        #   --retry 5: Retry up to 5 times
        #   --retry-delay 5: Wait 5 seconds between retries
        #   --fail: Return error code on HTTP errors
        #   --silent + --show-error: Only show errors, not progress to stderr
        #   -H "Content-Type": Required for Zenodo
        curl_cmd = [
            "curl",
            "--fail",                   # Return error on HTTP errors (4xx, 5xx)
            "--progress-bar",           # Show progress bar
            "--retry", "5",             # Retry 5 times
            "--retry-delay", "5",       # Wait 5 seconds between retries
            "--retry-connrefused",      # Retry on connection refused
            "--connect-timeout", "30",  # 30 second connection timeout
            "--max-time", "0",          # No timeout for transfer (unlimited)
            "-H", "Content-Type: application/octet-stream",
            "--upload-file", str(file_path),
            upload_url
        ]
        
        try:
            # Run curl with real-time progress streaming
            # Use Popen to capture stderr (where curl outputs progress) and show it in real-time
            print(f"    ⏳ Starting curl upload...")
            print(f"    💡 This may take a while for large files. Progress will be shown below.")
            print(f"    📊 Estimated time: {file_size/1024/1024/4:.0f} minutes at 4 Mbps upload speed")
            print()
            
            # Use Popen with direct output to terminal for real-time progress
            # stderr=None means curl's progress goes directly to terminal
            process = subprocess.Popen(
                curl_cmd,
                stdout=subprocess.PIPE,  # Capture stdout (curl response)
                stderr=None,  # Let progress bar go directly to terminal
            )
            
            # Wait for completion
            stdout, _ = process.communicate()
            result_code = process.returncode
            
            print()  # Newline after progress bar
            
            if result_code == 0:
                print(f"    ✓ curl upload successful!")
                return True
            else:
                print(f"    ✗ curl upload failed with exit code {result_code}")
                        
                # Provide guidance for common curl errors
                if result_code == 22:  # HTTP error
                    print(f"      ℹ HTTP error occurred. Check Zenodo status or token validity.")
                elif result_code == 28:  # Timeout
                    print(f"      ℹ Connection timed out. Check network connectivity.")
                elif result_code == 7:   # Connection refused
                    print(f"      ℹ Connection refused. Zenodo may be temporarily unavailable.")
                elif result_code == 56:  # Receive error
                    print(f"      ℹ Network receive error. Check network stability.")
                
                return False
                
        except subprocess.TimeoutExpired:
            print(f"    ✗ curl upload timed out")
            return False
        except FileNotFoundError:
            print(f"    ✗ curl not found! Please install curl or add it to PATH.")
            print(f"      On Windows: winget install curl  OR  choco install curl")
            return False
        except subprocess.SubprocessError as e:
            print(f"    ✗ curl subprocess error: {e}")
            return False
        except Exception as e:
            print(f"    ✗ Unexpected error running curl: {type(e).__name__}: {e}")
            return False
    
    def upload_file_bucket_api(
        self,
        bucket_url: str,
        file_path: Path,
        remote_name: str,
        progress_callback=None
    ) -> bool:
        """
        Upload a file using the bucket API (supports large files up to 50GB).
        
        CRITICAL: This method opens the file fresh for EACH call to avoid urllib3 #459 bug
        where file generators are not reset on retry. The HTTPAdapter retry is disabled
        for PUT requests, so all retries are handled at the application level by the
        calling upload_file() method which calls this method multiple times.
        
        Uses optimized chunking and timeout configuration for large files:
        - Files > 1GB: Use 50MB chunks, no read timeout
        - Files < 1GB: Use 10MB chunks, 1-hour read timeout
        - Always use 30-second connect timeout
        - Keep-alive connection headers
        
        Parameters
        ----------
        bucket_url : str
            Bucket URL from deposit links
        file_path : Path
            Local file path
        remote_name : str
            Name for the file in Zenodo
        progress_callback : callable
            Optional callback for progress updates
        
        Returns
        -------
        bool
            True if upload successful
        """
        upload_url = f"{bucket_url}/{remote_name}"
        file_size = file_path.stat().st_size
        
        # Determine if this is a large file requiring special handling
        is_large_file = file_size > self.config.LARGE_FILE_THRESHOLD
        chunk_size = (
            self.config.LARGE_FILE_CHUNK_SIZE if is_large_file 
            else self.config.CHUNK_SIZE
        )
        
        # Configure timeout based on file size
        # Tuple: (connect_timeout, read_timeout)
        # For large files, use no read timeout (None) to prevent premature disconnection
        timeout = (
            self.config.CONNECT_TIMEOUT,
            self.config.LARGE_FILE_READ_TIMEOUT if is_large_file 
            else self.config.READ_TIMEOUT
        )
        
        if is_large_file:
            file_size_str = format_size(file_size)
            print(f"    ℹ Large file detected ({file_size_str}): Using optimized upload strategy")
            print(f"      - Chunk size: {format_size(chunk_size)}")
            print(f"      - Connect timeout: {timeout[0]}s")
            print(f"      - Read timeout: {'None (unlimited)' if timeout[1] is None else f'{timeout[1]}s'}")
        
        # CRITICAL: Generator function that will be called fresh for this attempt
        # This ensures the file is read from the beginning, not from a previous position
        def file_reader():
            """Generator that reads file in chunks. Called fresh each retry attempt."""
            with open(file_path, 'rb') as f:
                bytes_sent = 0
                while True:
                    chunk = f.read(chunk_size)
                    if not chunk:
                        break
                    bytes_sent += len(chunk)
                    if progress_callback:
                        progress_callback(len(chunk))
                    yield chunk
        
        self._rate_limit_wait()
        
        try:
            # CRITICAL: Must include Content-Length header to prevent hanging
            # Connection and Keep-Alive headers for stable long-duration uploads
            headers = {
                "Content-Type": "application/octet-stream",
                "Content-Length": str(file_size),
                "Connection": "keep-alive",
                "Keep-Alive": "timeout=3600, max=1"  # Keep connection alive for 1 hour
            }
            
            # Use streaming upload with optimized settings
            # Note: We pass the generator directly - requests will call it once
            # If this attempt fails, the calling method will call us again with a NEW generator
            response = self.session.put(
                upload_url,
                data=file_reader(),  # Fresh generator for this attempt
                headers=headers,
                timeout=timeout,
                stream=False  # Don't stream response, wait for complete confirmation
            )
            
            if response.status_code in [200, 201]:
                return True
            
            print(f"    ✗ Upload failed: {response.status_code}")
            print(f"    Response: {response.text[:500]}")
            return False
            
        except requests.exceptions.Timeout as e:
            print(f"    ✗ Timeout uploading file: {e}")
            print(f"      File size: {format_size(file_size)}")
            print(f"      Timeout config: connect={timeout[0]}s, read={timeout[1]}")
            return False
        except requests.exceptions.ConnectionError as e:
            error_str = str(e)
            print(f"    ✗ Connection error: {error_str[:200]}")
            
            # Provide specific guidance for Windows socket errors
            if "10054" in error_str or "10053" in error_str:
                print(f"      ⚠ Windows socket error detected:")
                print(f"        10053 = Connection aborted by local software")
                print(f"        10054 = Connection forcibly closed by remote host")
                print(f"        This often indicates network instability or server-side issues")
            
            return False
        except Exception as e:
            print(f"    ✗ Unexpected error: {type(e).__name__}: {e}")
            return False
    
    def upload_file(
        self,
        deposit: Dict,
        file_entry: FileEntry,
        progress: UploadProgress
    ) -> bool:
        """
        Upload a file to a deposit using curl for ALL files.
        
        CRITICAL FIX: Python requests on Windows has multiple issues:
        - Error 10053/10054 socket errors for medium/large files
        - URL path encoding issues causing 400 Bad Request
        
        curl uses native OS networking (WinHTTP on Windows) which:
        - Has proper TCP buffer management
        - Handles long-running uploads reliably
        - Works consistently for all file sizes
        
        Based on proven solution from jhpoelen/zenodo-upload:
        https://github.com/jhpoelen/zenodo-upload
        
        Parameters
        ----------
        deposit : Dict
            Deposit data (must include links.bucket)
        file_entry : FileEntry
            File to upload
        progress : UploadProgress
            Progress tracker
        
        Returns
        -------
        bool
            True if upload successful
        """
        # Check if already uploaded
        if file_entry.remote_name in progress.uploaded_files:
            print(f"  ⏭ Already uploaded: {file_entry.remote_name}")
            return True
        
        bucket_url = deposit.get("links", {}).get("bucket")
        if not bucket_url:
            print("  ✗ No bucket URL found in deposit")
            return False
        
        file_size = file_entry.size_bytes
        file_size_str = format_size(file_size)
        
        print(f"  📤 Uploading: {file_entry.remote_name} ({file_size_str})")
        
        # =================================================================
        # USE CURL FOR ALL FILES - PROVEN RELIABLE ON WINDOWS
        # =================================================================
        if self._is_curl_available():
            for attempt in range(self.config.MAX_RETRIES):
                attempt_start = time.time()
                
                success = self.upload_file_curl(
                    bucket_url,
                    file_entry.local_path,
                    file_entry.remote_name
                )
                
                attempt_duration = time.time() - attempt_start
                
                if success:
                    progress.uploaded_files.add(file_entry.remote_name)
                    progress.total_bytes_uploaded += file_size
                    if file_entry.remote_name in progress.failed_files:
                        progress.failed_files.remove(file_entry.remote_name)
                    
                    print(f"    ✓ Upload SUCCESS! ({attempt_duration:.1f}s)")
                    return True
                
                # curl failed - retry with backoff
                if attempt < self.config.MAX_RETRIES - 1:
                    base_wait = min(
                        self.config.INITIAL_BACKOFF_SECONDS * (2 ** attempt),
                        self.config.MAX_BACKOFF_SECONDS
                    )
                    jitter = random.uniform(0.5, 1.5)
                    wait_time = base_wait * jitter
                    
                    print(f"    ↻ Retry {attempt + 2}/{self.config.MAX_RETRIES} in {wait_time:.1f}s...")
                    time.sleep(wait_time)
            
            # All curl retries failed
            progress.failed_files.add(file_entry.remote_name)
            print(f"    ✗ Failed after {self.config.MAX_RETRIES} curl attempts")
            return False
        
        # =================================================================
        # FALLBACK: Python requests (only if curl not available)
        # =================================================================
        print(f"    ⚠ curl not available - using Python requests (less reliable on Windows)")
        print(f"      Tip: Install curl: winget install curl")
        
        # Create progress bar if available
        pbar = None
        if tqdm and file_size > 10 * 1024 * 1024:  # > 10MB
            pbar = tqdm(total=file_size, unit='B', unit_scale=True)
        
        def update_progress(bytes_sent):
            if pbar:
                pbar.update(bytes_sent)
        
        attempt_times = []
        
        try:
            for attempt in range(self.config.MAX_RETRIES):
                attempt_start = time.time()
                
                success = self.upload_file_bucket_api(
                    bucket_url,
                    file_entry.local_path,
                    file_entry.remote_name,
                    progress_callback=update_progress if pbar else None
                )
                
                attempt_duration = time.time() - attempt_start
                attempt_times.append(attempt_duration)
                
                if success:
                    progress.uploaded_files.add(file_entry.remote_name)
                    progress.total_bytes_uploaded += file_size
                    
                    if file_entry.remote_name in progress.failed_files:
                        progress.failed_files.remove(file_entry.remote_name)
                    
                    print(f"    ✓ Success ({attempt_duration:.1f}s)")
                    return True
                
                if attempt < self.config.MAX_RETRIES - 1:
                    base_wait = min(
                        self.config.INITIAL_BACKOFF_SECONDS * (2 ** attempt),
                        self.config.MAX_BACKOFF_SECONDS
                    )
                    jitter = random.uniform(0.5, 1.5)
                    wait_time = base_wait * jitter
                    
                    print(f"    ↻ Retry {attempt + 2}/{self.config.MAX_RETRIES} in {wait_time:.1f}s...")
                    time.sleep(wait_time)
                    
                    if pbar:
                        pbar.reset()
            
            progress.failed_files.add(file_entry.remote_name)
            print(f"    ✗ Failed after {self.config.MAX_RETRIES} attempts")
            return False
            
        finally:
            if pbar:
                pbar.close()
    
    def publish_deposit(self, deposit_id: str) -> Dict:
        """Publish a deposit."""
        url = f"{self.base_url}/deposit/depositions/{deposit_id}/actions/publish"
        
        response = self._request_with_retry("POST", url)
        
        if response.status_code in [200, 201, 202]:
            return response.json()
        
        raise Exception(f"Failed to publish: {response.status_code} - {response.text}")


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def format_size(size_bytes: int) -> str:
    """Format bytes as human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"


def get_folder_config() -> Dict[str, dict]:
    """
    Get configuration for which folders to ZIP.
    
    This is the key decision: which folders get zipped vs uploaded individually.
    
    Strategy:
    - "flat_only": Upload only files at root level, skip subdirectories
                   (Use when existing .tar.gz/.zip already contain subdirs)
    - "zip": True - Create ZIP archive of entire directory
    - "zip": False - Upload files individually (including recursive)
    - "skip_subdirs": List of subdirectory names to skip
    
    Returns
    -------
    Dict[str, dict]
        Folder configurations
    """
    return {
        # FLAMES - Has existing .tar.gz archives! Upload those, skip extracted dirs
        "data/external/flames": {
            "flat_only": True,  # Upload only root-level files (incl. tar.gz)
            "skip_subdirs": ["Annotation_data", "example_data", "model"],  # Skip extracted
        },
        
        # MODERATE FILE COLLECTIONS - Zip by subfolder
        "data/external/ot_25_12_l2g_prediction": {"zip": True},  # 404 files
        "data/external/cs2g": {"zip": True},  # 144 files
        "data/external/mpra_abell_2022": {"zip": True},  # 130 files
        "data/external/E2G_benchmarking": {"zip": True},  # 128 files
        "data/external/mpra_finemapped": {"zip": True},  # 121 files
        "data/external/open_targets": {"zip": True},  # 83 files
        "data/external/gold_standards": {"zip": True},  # 29 files
        "data/external/crispr_benchmark": {"zip": True},  # 29 files
        "data/external/pops": {"zip": True},  # 27 files
        "data/external/opentargets_l2g": {"zip": True},  # 22 files
        "data/external/mpra_schizophrenia": {"zip": True},  # 13 files
        
        # LARGE SINGLE FILES - Upload directly
        "data/external/drug_targets": {
            "flat_only": True,  # Skip chembl_36 extracted folder
            "skip_subdirs": ["chembl_36"],  # Keep only tar.gz
        },
        "data/external/crispr_validation": {"zip": False},  # 4 files, 3.62 GB
        "data/external/gtex": {"zip": False},  # 2 files, 1.46 GB
        "data/external/clinvar": {"zip": False},  # 1 file, 0.38 GB
        "data/external/abc": {"zip": False},  # 1 file, 0.32 GB
        "data/external/encode": {"zip": False},  # 1 file
        "data/external/ensembl": {"zip": False},  # 3 files
        "data/external/gnomad": {"zip": False},  # 1 file
        "data/external/benchmarks": {"zip": False},  # 4 files
        "data/external/benchmark_gold_standards": {"zip": False},  # 1 file
        "data/external/sting_seq": {"zip": False},  # 8 files
        
        # RAW DATA - Most are already compressed, upload individually
        "data/raw/gwas_sumstats": {"zip": False},  # 15 files, 25.40 GB - already compressed
        "data/raw/eqtl": {"zip": False},  # 1 file, 1.46 GB
        "data/raw/regulatory_annotations": {"zip": True},  # 3 files
        "data/raw/reference": {"zip": False},  # 1 file
        "data/raw/chromatin_interactions": {"zip": True},  # 4 files
        "data/raw/functional": {"zip": False},  # 1 file
        "data/raw/gene_annotations": {"zip": False},  # 2 files
        "data/raw/gwas_catalog": {"zip": False},  # 1 file
        "data/raw/ld_reference": {"zip": False},  # 1 file
        
        # PROCESSED DATA - Usually small, zip all
        "data/processed": {"zip": True},  # 93 files, 0.04 GB
        
        # BENCHMARKS
        "benchmarks": {"zip": True},  # Small
        
        # MANIFESTS
        "data/manifests": {"zip": True},  # Keep manifests accessible
    }


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Upload large dataset to Zenodo with proper handling of limits",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
IMPORTANT: Zenodo has strict limits!
  - 100 files maximum per record
  - 50 GB maximum per record (can request 200 GB quota)

This script creates ZIP archives to stay under the file limit.

WORKFLOW:
  1. First, create ZIP archives:
     python scripts/upload_zenodo_v3.py --create-zips
     
  2. Then upload:
     python scripts/upload_zenodo_v3.py --upload
     
  3. Or do both:
     python scripts/upload_zenodo_v3.py --create-zips --upload

RESUME:
  The script saves progress to zenodo_upload_progress.json
  If upload is interrupted, re-run and it will resume.

QUOTA INCREASE:
  If your data exceeds 50GB, request a quota increase:
  https://zenodo.org/support
        """
    )
    
    parser.add_argument(
        "--token",
        type=str,
        default=os.environ.get("ZENODO_TOKEN", "v0vwEqX8u9dw6MUFZqAQJSGjwcqA3JImFA5zQbPJx4MIJrhlfQgVp77jJz7p"),
        help="Zenodo access token (or set ZENODO_TOKEN env var)"
    )
    
    parser.add_argument(
        "--record-id",
        type=str,
        default="17880202",
        help="Record ID (from DOI 10.5281/zenodo.XXXXX). Use actual record ID, not concept ID."
    )
    
    parser.add_argument(
        "--deposit-id",
        type=str,
        default=None,
        help="Use existing draft deposit ID directly (skips version creation)"
    )
    
    parser.add_argument(
        "--version",
        type=str,
        default="2",
        help="Version number for the new version"
    )
    
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).parent.parent,
        help="Project root directory"
    )
    
    parser.add_argument(
        "--create-zips",
        action="store_true",
        help="Create ZIP archives for large directories"
    )
    
    parser.add_argument(
        "--upload",
        action="store_true",
        help="Upload files to Zenodo"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without doing it"
    )
    
    parser.add_argument(
        "--sandbox",
        action="store_true",
        help="Use Zenodo sandbox for testing"
    )
    
    parser.add_argument(
        "--skip-delete",
        action="store_true",
        default=True,
        help="Don't delete existing files in draft (default: True, since new versions start clean)"
    )
    
    parser.add_argument(
        "--force-delete",
        action="store_true",
        help="Force deletion of existing files (use if files were imported from previous version)"
    )
    
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous progress"
    )
    
    parser.add_argument(
        "--publish",
        action="store_true",
        help="Publish after successful upload"
    )
    
    args = parser.parse_args()
    
    if not args.create_zips and not args.upload:
        print("ERROR: Specify --create-zips and/or --upload")
        print()
        print("Examples:")
        print("  python scripts/upload_zenodo_v3.py --create-zips")
        print("  python scripts/upload_zenodo_v3.py --upload")
        print("  python scripts/upload_zenodo_v3.py --create-zips --upload")
        sys.exit(1)
    
    if args.upload and not args.token:
        print("ERROR: Zenodo access token required for upload!")
        print()
        print("Set environment variable:")
        print("  export ZENODO_TOKEN='your_token_here'")
        sys.exit(1)
    
    root_dir = args.root.resolve()
    config = ZenodoConfig()
    
    print("=" * 70)
    print("Zenodo Upload Script v3 - Large Dataset Handler")
    print("=" * 70)
    print(f"Project root: {root_dir}")
    print(f"Zenodo limits: {config.MAX_FILES_PER_RECORD} files, {config.MAX_SIZE_PER_RECORD_GB} GB")
    print()
    
    # Setup directories
    zip_output_dir = root_dir / "zenodo_archives"
    progress_file = root_dir / "zenodo_upload_progress.json"
    
    # Load progress
    if args.resume:
        progress = UploadProgress.load(progress_file)
        print(f"Resuming: {len(progress.uploaded_files)} files already uploaded")
    else:
        progress = UploadProgress()
    
    if not progress.start_time:
        progress.start_time = datetime.now()
    
    # Get folder configuration
    folder_config = get_folder_config()
    
    # ==========================================================================
    # STEP 1: Create ZIP archives
    # ==========================================================================
    
    files_to_upload = []
    
    if args.create_zips:
        print()
        print("=" * 70)
        print("STEP 1: Creating ZIP Archives")
        print("=" * 70)
        print(f"Output directory: {zip_output_dir}")
        
        if args.dry_run:
            print("DRY RUN - Would create archives for:")
            for folder, cfg in folder_config.items():
                if cfg.get('flat_only', False):
                    print(f"  📂 {folder} (root files only, skip subdirs)")
                elif cfg.get('zip', False):
                    print(f"  📦 {folder}")
                elif cfg.get('split_zip', False):
                    print(f"  📦 {folder} (split by subdir)")
                else:
                    print(f"  📄 {folder} (individual files)")
        else:
            zip_manager = ZipArchiveManager(root_dir, zip_output_dir, config)
            files_to_upload = zip_manager.create_all_archives(folder_config)
    
    # NOTE: create_all_archives already adds:
    # - ZIP archives for zip: True folders
    # - Root-level files for flat_only: True folders  
    # - Individual files for zip: False folders (with skip_subdirs support)
    # So we DON'T need to add more files here - that was causing double-counting!
    
    # Add non-zipped files from config (for upload)
    if args.upload:
        if not files_to_upload:
            # Load from existing zips if we didn't just create them
            print()
            print("Loading existing archives...")
            
            if zip_output_dir.exists():
                for zip_file in zip_output_dir.glob("*.zip"):
                    files_to_upload.append(FileEntry(
                        local_path=zip_file,
                        remote_name=zip_file.name,
                        size_bytes=zip_file.stat().st_size,
                        is_zip=True
                    ))
            
            # Add non-zipped files (respecting flat_only and skip_subdirs)
            for folder, cfg in folder_config.items():
                folder_path = root_dir / folder
                if not folder_path.exists():
                    continue
                
                # Handle flat_only - add only root-level files
                if cfg.get('flat_only', False):
                    for item in folder_path.iterdir():
                        if item.is_file():
                            files_to_upload.append(FileEntry(
                                local_path=item,
                                remote_name=str(item.relative_to(root_dir)).replace('\\', '/'),
                                size_bytes=item.stat().st_size,
                                is_zip=False
                            ))
                # Handle non-zip folders (recursive with skip support)
                elif not cfg.get('zip', False):
                    skip_subdirs = set(cfg.get('skip_subdirs', []))
                    for f in folder_path.rglob("*"):
                        if f.is_file():
                            # Check if file is in a skipped subdirectory
                            try:
                                rel_to_folder = f.relative_to(folder_path)
                                parts = rel_to_folder.parts
                                if parts and parts[0] in skip_subdirs:
                                    continue  # Skip this file
                            except ValueError:
                                pass
                            
                            files_to_upload.append(FileEntry(
                                local_path=f,
                                remote_name=str(f.relative_to(root_dir)).replace('\\', '/'),
                                size_bytes=f.stat().st_size,
                                is_zip=False
                            ))
    
    # Add manifest and other root files from data folder
    manifest_path = root_dir / "data" / "MANIFEST.json"
    if manifest_path.exists() and args.upload:
        files_to_upload.append(FileEntry(
            local_path=manifest_path,
            remote_name="data/MANIFEST.json",
            size_bytes=manifest_path.stat().st_size,
            is_zip=False
        ))
    
    # ==========================================================================
    # Summary
    # ==========================================================================
    
    print()
    print("=" * 70)
    print("UPLOAD SUMMARY")
    print("=" * 70)
    
    total_files = len(files_to_upload)
    total_size = sum(f.size_bytes for f in files_to_upload)
    zip_files = [f for f in files_to_upload if f.is_zip]
    regular_files = [f for f in files_to_upload if not f.is_zip]
    
    print(f"Total files to upload: {total_files}")
    print(f"  - ZIP archives: {len(zip_files)}")
    print(f"  - Regular files: {len(regular_files)}")
    print(f"Total size: {format_size(total_size)}")
    print()
    
    # Check against limits
    if total_files > config.MAX_FILES_SAFE:
        print(f"⚠ WARNING: {total_files} files exceeds safe limit ({config.MAX_FILES_SAFE})")
        print("  Consider zipping more directories.")
    else:
        print(f"✓ File count OK: {total_files} < {config.MAX_FILES_SAFE}")
    
    total_size_gb = total_size / (1024**3)
    if total_size_gb > config.MAX_SIZE_PER_RECORD_SAFE_GB:
        print(f"⚠ WARNING: {total_size_gb:.2f} GB exceeds safe limit ({config.MAX_SIZE_PER_RECORD_SAFE_GB} GB)")
        print("  Request quota increase from Zenodo: https://zenodo.org/support")
    else:
        print(f"✓ Size OK: {total_size_gb:.2f} GB < {config.MAX_SIZE_PER_RECORD_SAFE_GB} GB")
    
    print()
    
    if args.dry_run:
        print("DRY RUN - Files that would be uploaded:")
        for f in sorted(files_to_upload, key=lambda x: -x.size_bytes)[:30]:
            marker = "📦" if f.is_zip else "📄"
            print(f"  {marker} {f.remote_name} ({format_size(f.size_bytes)})")
        if len(files_to_upload) > 30:
            print(f"  ... and {len(files_to_upload) - 30} more")
        print()
        print("DRY RUN complete. Run without --dry-run to proceed.")
        sys.exit(0)
    
    if not args.upload:
        print("ZIP creation complete. Run with --upload to upload to Zenodo.")
        sys.exit(0)
    
    # ==========================================================================
    # STEP 2: Upload to Zenodo
    # ==========================================================================
    
    print()
    print("=" * 70)
    print("STEP 2: Upload to Zenodo")
    print("=" * 70)
    
    # Initialize client
    client = ZenodoClient(
        access_token=args.token,
        sandbox=args.sandbox,
        config=config
    )
    
    # Get or create draft
    # PRIORITY ORDER for deposit ID:
    # 1. --deposit-id command line argument (explicit override)
    # 2. Saved deposit_id from progress file (for resume)
    # 3. Create new version from record_id (fresh start)
    
    effective_deposit_id = args.deposit_id or (progress.deposit_id if args.resume else None)
    
    if effective_deposit_id:
        source = "command line" if args.deposit_id else "saved progress (resume)"
        print(f"Using existing draft deposit ID: {effective_deposit_id} (from {source})...")
        try:
            deposit = client.get_deposit(effective_deposit_id)
            deposit_id = effective_deposit_id
            state = deposit.get('state', 'unknown')
            existing_files = len(deposit.get('files', []))
            print(f"  ✓ Found deposit: {state} state, {existing_files} files already uploaded")
            
            # Store in progress for future resumes
            progress.deposit_id = deposit_id
            progress.bucket_url = deposit.get("links", {}).get("bucket")
        except Exception as e:
            print(f"  ✗ Error: {e}")
            sys.exit(1)
    else:
        print(f"Connecting to Zenodo record {args.record_id}...")
        try:
            deposit = client.get_or_create_draft(args.record_id)
            deposit_id = deposit.get("id")
            print(f"  ✓ Draft deposit ID: {deposit_id}")
            
            # Store in progress for future resumes
            progress.deposit_id = deposit_id
            progress.bucket_url = deposit.get("links", {}).get("bucket")
            progress.save(progress_file)  # Save immediately so we don't lose it
            print(f"  ✓ Deposit ID saved to progress file for resume capability")
        except Exception as e:
            print(f"  ✗ Error: {e}")
            sys.exit(1)
    
    # Update metadata
    print()
    print("Updating metadata...")
    try:
        current_metadata = deposit.get("metadata", {})
        current_metadata["version"] = args.version
        client.update_metadata(deposit_id, current_metadata)
        print(f"  ✓ Version set to {args.version}")
    except Exception as e:
        print(f"  ⚠ Warning: {e}")
    
    # Delete existing files if explicitly requested
    # Note: New version drafts start clean, so deletion is usually unnecessary
    if args.force_delete and not args.resume:
        print()
        print("Deleting existing files...")
        print("  ℹ Note: New versions typically start with no files")
        try:
            deleted = client.delete_all_files(deposit_id)
            print(f"  ✓ Deleted {deleted} files")
        except Exception as e:
            print(f"  ⚠ Deletion failed (continuing anyway): {e}")
            print("  This is expected for new version drafts")
    else:
        print()
        print("Skipping file deletion (new versions start clean)")
        print("  Use --force-delete if you imported files from previous version")
    
    # Refresh deposit to get bucket URL
    deposit = client.get_deposit(deposit_id)
    
    # Upload files
    print()
    print("=" * 70)
    print("Uploading files...")
    print("=" * 70)
    
    success_count = 0
    fail_count = 0
    
    # Sort files: larger files first (more important)
    files_to_upload.sort(key=lambda x: -x.size_bytes)
    
    for i, file_entry in enumerate(files_to_upload, 1):
        print(f"\n[{i}/{total_files}]", end=" ")
        
        try:
            if client.upload_file(deposit, file_entry, progress):
                success_count += 1
            else:
                fail_count += 1
        except Exception as e:
            print(f"  ✗ Exception: {e}")
            fail_count += 1
            progress.failed_files.add(file_entry.remote_name)
        
        # Save progress periodically
        if i % 5 == 0:
            progress.save(progress_file)
    
    # Final progress save
    progress.save(progress_file)
    
    # ==========================================================================
    # Results
    # ==========================================================================
    
    print()
    print("=" * 70)
    print("UPLOAD COMPLETE")
    print("=" * 70)
    print(f"Successful: {success_count}")
    print(f"Failed: {fail_count}")
    print(f"Total uploaded: {format_size(progress.total_bytes_uploaded)}")
    print()
    
    if fail_count > 0:
        print("⚠ Some uploads failed. Run again with --resume to retry.")
        print(f"   Progress saved to: {progress_file}")
    
    # Publish if requested
    if args.publish and fail_count == 0:
        print()
        print("Publishing deposit...")
        try:
            result = client.publish_deposit(deposit_id)
            doi = result.get("doi", "N/A")
            print(f"  ✓ Published!")
            print(f"  DOI: {doi}")
        except Exception as e:
            print(f"  ✗ Error: {e}")
    elif fail_count == 0:
        print(f"Review at: https://zenodo.org/deposit/{deposit_id}")
        print("Run with --publish to publish, or publish manually.")
    
    print()
    print("Done!")


if __name__ == "__main__":
    main()
