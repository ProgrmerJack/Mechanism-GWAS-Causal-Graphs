#!/usr/bin/env python3
"""
Manifest Verification Script

Verifies data integrity by checking file hashes against manifests.
Essential for Nature-level reproducibility requirements.

Usage:
    python scripts/verify_manifests.py --manifest data/manifests/
    python scripts/verify_manifests.py --manifest data/manifests/gwas_sumstats.yaml --data data/gwas/
"""

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Optional

import yaml


def compute_md5(filepath: Path, chunk_size: int = 8192) -> str:
    """Compute MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def compute_sha256(filepath: Path, chunk_size: int = 8192) -> str:
    """Compute SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def verify_file(
    filepath: Path, expected_md5: Optional[str] = None, expected_sha256: Optional[str] = None
) -> dict:
    """
    Verify a file against expected hashes.

    Returns:
        dict with verification results
    """
    result = {
        "file": str(filepath),
        "exists": filepath.exists(),
        "md5_match": None,
        "sha256_match": None,
        "computed_md5": None,
        "computed_sha256": None,
    }

    if not filepath.exists():
        return result

    # Compute hashes
    result["computed_md5"] = compute_md5(filepath)
    result["computed_sha256"] = compute_sha256(filepath)

    # Verify against expected
    if expected_md5 and not expected_md5.startswith("placeholder"):
        result["md5_match"] = result["computed_md5"] == expected_md5

    if expected_sha256 and not expected_sha256.startswith("placeholder"):
        result["sha256_match"] = result["computed_sha256"] == expected_sha256

    return result


def load_manifest(manifest_path: Path) -> dict:
    """Load a YAML manifest file."""
    with open(manifest_path) as f:
        return yaml.safe_load(f)


def extract_files_from_manifest(manifest: dict, data_dir: Path) -> list:
    """
    Extract file entries from manifest for verification.

    Returns list of dicts with file info.
    """
    files_to_verify = []

    def _extract_recursive(obj, path=""):
        if isinstance(obj, dict):
            # Check if this is a file entry
            if "filename" in obj or "md5" in obj or "sha256" in obj:
                entry = {
                    "filename": obj.get("filename"),
                    "md5": obj.get("md5"),
                    "sha256": obj.get("sha256"),
                    "manifest_path": path,
                }
                if entry["filename"]:
                    files_to_verify.append(entry)

            # Recurse into nested dicts
            for key, value in obj.items():
                _extract_recursive(value, f"{path}.{key}" if path else key)

        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                _extract_recursive(item, f"{path}[{i}]")

    _extract_recursive(manifest)
    return files_to_verify


def verify_manifest(manifest_path: Path, data_dir: Optional[Path] = None) -> dict:
    """
    Verify all files in a manifest.

    Args:
        manifest_path: Path to YAML manifest
        data_dir: Optional data directory root

    Returns:
        Verification report dict
    """
    manifest = load_manifest(manifest_path)
    files = extract_files_from_manifest(manifest, data_dir or Path("."))

    report = {
        "manifest": str(manifest_path),
        "manifest_version": manifest.get("manifest_version", "unknown"),
        "total_files": len(files),
        "verified": 0,
        "missing": 0,
        "hash_mismatch": 0,
        "placeholder_hashes": 0,
        "details": [],
    }

    for file_entry in files:
        filename = file_entry["filename"]
        if data_dir:
            filepath = data_dir / filename
        else:
            filepath = Path(filename)

        result = verify_file(
            filepath, expected_md5=file_entry.get("md5"), expected_sha256=file_entry.get("sha256")
        )

        result["manifest_path"] = file_entry["manifest_path"]

        if not result["exists"]:
            report["missing"] += 1
            result["status"] = "MISSING"
        elif result["md5_match"] is False or result["sha256_match"] is False:
            report["hash_mismatch"] += 1
            result["status"] = "HASH_MISMATCH"
        elif result["md5_match"] is None and result["sha256_match"] is None:
            report["placeholder_hashes"] += 1
            result["status"] = "PLACEHOLDER_HASH"
        else:
            report["verified"] += 1
            result["status"] = "VERIFIED"

        report["details"].append(result)

    return report


def print_report(report: dict, verbose: bool = False):
    """Print verification report."""
    print(f"\n{'=' * 60}")
    print(f"Manifest Verification Report")
    print(f"{'=' * 60}")
    print(f"Manifest: {report['manifest']}")
    print(f"Version: {report['manifest_version']}")
    print(f"\nSummary:")
    print(f"  Total files:        {report['total_files']}")
    print(f"  Verified:           {report['verified']}")
    print(f"  Missing:            {report['missing']}")
    print(f"  Hash mismatch:      {report['hash_mismatch']}")
    print(f"  Placeholder hashes: {report['placeholder_hashes']}")

    if verbose or report["missing"] > 0 or report["hash_mismatch"] > 0:
        print(f"\nDetails:")
        for detail in report["details"]:
            status_icon = {
                "VERIFIED": "✓",
                "MISSING": "✗",
                "HASH_MISMATCH": "!",
                "PLACEHOLDER_HASH": "?",
            }.get(detail["status"], "?")

            print(f"  [{status_icon}] {detail['file']}: {detail['status']}")

            if verbose and detail["computed_md5"]:
                print(f"      MD5: {detail['computed_md5']}")
            if verbose and detail["computed_sha256"]:
                print(f"      SHA256: {detail['computed_sha256']}")


def main():
    parser = argparse.ArgumentParser(description="Verify data manifests for reproducibility")
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="Path to manifest YAML file or directory containing manifests",
    )
    parser.add_argument("--data", type=Path, default=None, help="Data directory root")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--update-hashes", action="store_true", help="Update placeholder hashes in manifest")

    args = parser.parse_args()

    manifest_paths = []
    if args.manifest.is_dir():
        manifest_paths = list(args.manifest.glob("*.yaml"))
    else:
        manifest_paths = [args.manifest]

    all_verified = True
    total_missing = 0
    total_mismatch = 0

    for manifest_path in manifest_paths:
        if manifest_path.name == "README.md":
            continue

        report = verify_manifest(manifest_path, args.data)
        print_report(report, args.verbose)

        total_missing += report["missing"]
        total_mismatch += report["hash_mismatch"]

        if report["missing"] > 0 or report["hash_mismatch"] > 0:
            all_verified = False

    print(f"\n{'=' * 60}")
    print(f"Overall Result: {'PASS' if all_verified else 'FAIL'}")
    if total_missing > 0:
        print(f"  {total_missing} files missing - download required data")
    if total_mismatch > 0:
        print(f"  {total_mismatch} hash mismatches - data may be corrupted")
    print(f"{'=' * 60}")

    sys.exit(0 if all_verified else 1)


if __name__ == "__main__":
    main()
