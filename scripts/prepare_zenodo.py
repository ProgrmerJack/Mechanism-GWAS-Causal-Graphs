#!/usr/bin/env python3
"""
Zenodo Archive Preparation Script

Prepares the repository for Zenodo archival by:
1. Validating all required files exist
2. Computing checksums for key files
3. Creating a release manifest
4. Validating CITATION.cff and .zenodo.json

Usage:
    python scripts/prepare_zenodo.py --output release/
    python scripts/prepare_zenodo.py --validate-only
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml


def compute_sha256(filepath: Path) -> str:
    """Compute SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def get_git_info() -> dict:
    """Get current git information."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
        ).strip()
        dirty = (
            subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
            != ""
        )
        return {
            "commit": commit,
            "branch": branch,
            "dirty": dirty,
        }
    except subprocess.CalledProcessError:
        return {"commit": "unknown", "branch": "unknown", "dirty": True}


def validate_citation_cff(root_dir: Path) -> list:
    """Validate CITATION.cff file."""
    errors = []
    citation_path = root_dir / "CITATION.cff"

    if not citation_path.exists():
        errors.append("CITATION.cff file not found")
        return errors

    with open(citation_path) as f:
        try:
            citation = yaml.safe_load(f)
        except yaml.YAMLError as e:
            errors.append(f"CITATION.cff YAML parsing error: {e}")
            return errors

    required_fields = ["cff-version", "title", "authors", "license"]
    for field in required_fields:
        if field not in citation:
            errors.append(f"CITATION.cff missing required field: {field}")

    return errors


def validate_zenodo_json(root_dir: Path) -> list:
    """Validate .zenodo.json file."""
    errors = []
    zenodo_path = root_dir / ".zenodo.json"

    if not zenodo_path.exists():
        errors.append(".zenodo.json file not found")
        return errors

    with open(zenodo_path) as f:
        try:
            zenodo = json.load(f)
        except json.JSONDecodeError as e:
            errors.append(f".zenodo.json JSON parsing error: {e}")
            return errors

    if "metadata" not in zenodo:
        errors.append(".zenodo.json missing 'metadata' key")
        return errors

    required_metadata = ["title", "description", "creators", "upload_type", "license"]
    for field in required_metadata:
        if field not in zenodo["metadata"]:
            errors.append(f".zenodo.json missing metadata field: {field}")

    return errors


def validate_required_files(root_dir: Path) -> list:
    """Check that all required files exist."""
    errors = []

    required_files = [
        "README.md",
        "LICENSE",
        "pyproject.toml",
        "environment.yml",
        "Dockerfile",
        "CITATION.cff",
        ".zenodo.json",
        "workflow/Snakefile",
        "workflow/config.yaml",
        "manuscript/main.tex",
        "manuscript/references.bib",
    ]

    required_dirs = [
        "src",
        "scripts",
        "tests",
        "config",
        "data/manifests",
    ]

    for filepath in required_files:
        if not (root_dir / filepath).exists():
            errors.append(f"Required file missing: {filepath}")

    for dirpath in required_dirs:
        if not (root_dir / dirpath).is_dir():
            errors.append(f"Required directory missing: {dirpath}")

    return errors


def create_release_manifest(root_dir: Path, output_dir: Path) -> dict:
    """Create a manifest of all files to be archived."""
    manifest = {
        "version": "1.0.0",
        "created": datetime.now().isoformat(),
        "git": get_git_info(),
        "files": [],
    }

    # Patterns to include
    include_patterns = [
        "*.py",
        "*.R",
        "*.yaml",
        "*.yml",
        "*.json",
        "*.md",
        "*.tex",
        "*.bib",
        "*.toml",
        "*.cfg",
        "*.cff",
        "Dockerfile",
        "Snakefile",
        "LICENSE",
        ".zenodo.json",
    ]

    # Patterns to exclude
    exclude_patterns = [
        "__pycache__",
        "*.pyc",
        ".git",
        ".pytest_cache",
        "*.egg-info",
        ".mypy_cache",
        "node_modules",
        ".codacy",
    ]

    for pattern in include_patterns:
        for filepath in root_dir.rglob(pattern):
            # Skip excluded patterns
            if any(excl in str(filepath) for excl in exclude_patterns):
                continue

            rel_path = filepath.relative_to(root_dir)
            manifest["files"].append(
                {
                    "path": str(rel_path),
                    "size": filepath.stat().st_size,
                    "sha256": compute_sha256(filepath),
                }
            )

    # Sort files by path
    manifest["files"].sort(key=lambda x: x["path"])
    manifest["total_files"] = len(manifest["files"])

    return manifest


def prepare_archive(root_dir: Path, output_dir: Path) -> Path:
    """Prepare archive directory structure."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Files and directories to include in archive
    include = [
        "src",
        "scripts",
        "tests",
        "config",
        "workflow",
        "manuscript",
        "data/manifests",
        "README.md",
        "LICENSE",
        "pyproject.toml",
        "environment.yml",
        "Dockerfile",
        "CITATION.cff",
        ".zenodo.json",
    ]

    archive_dir = output_dir / "mechanism-gwas-causal-graphs-v1.0.0"
    archive_dir.mkdir(parents=True, exist_ok=True)

    for item in include:
        src_path = root_dir / item
        dst_path = archive_dir / item

        if src_path.is_file():
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_path, dst_path)
        elif src_path.is_dir():
            if dst_path.exists():
                shutil.rmtree(dst_path)
            shutil.copytree(
                src_path,
                dst_path,
                ignore=shutil.ignore_patterns(
                    "__pycache__", "*.pyc", ".git", ".pytest_cache", "*.egg-info"
                ),
            )

    return archive_dir


def main():
    parser = argparse.ArgumentParser(description="Prepare repository for Zenodo archival")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).parent.parent,
        help="Repository root directory",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("release"), help="Output directory for archive"
    )
    parser.add_argument(
        "--validate-only", action="store_true", help="Only validate, don't create archive"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")

    args = parser.parse_args()
    root_dir = args.root.resolve()

    print(f"Repository root: {root_dir}")
    print()

    # Run all validations
    print("Running validations...")
    all_errors = []

    errors = validate_required_files(root_dir)
    if errors:
        print(f"  ✗ Required files: {len(errors)} errors")
        all_errors.extend(errors)
    else:
        print("  ✓ Required files present")

    errors = validate_citation_cff(root_dir)
    if errors:
        print(f"  ✗ CITATION.cff: {len(errors)} errors")
        all_errors.extend(errors)
    else:
        print("  ✓ CITATION.cff valid")

    errors = validate_zenodo_json(root_dir)
    if errors:
        print(f"  ✗ .zenodo.json: {len(errors)} errors")
        all_errors.extend(errors)
    else:
        print("  ✓ .zenodo.json valid")

    # Check git status
    git_info = get_git_info()
    if git_info["dirty"]:
        print("  ⚠ Repository has uncommitted changes")
    else:
        print(f"  ✓ Git clean at commit {git_info['commit'][:8]}")

    print()

    if all_errors:
        print("Validation errors:")
        for error in all_errors:
            print(f"  - {error}")
        print()
        print("Please fix errors before archiving.")
        sys.exit(1)

    if args.validate_only:
        print("Validation passed!")
        sys.exit(0)

    # Create manifest
    print("Creating release manifest...")
    manifest = create_release_manifest(root_dir, args.output)

    manifest_path = args.output / "MANIFEST.json"
    args.output.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Manifest written to {manifest_path}")
    print(f"  Total files: {manifest['total_files']}")

    # Create archive
    print("\nPreparing archive...")
    archive_dir = prepare_archive(root_dir, args.output)
    print(f"  Archive prepared at {archive_dir}")

    # Create zip archive
    zip_path = shutil.make_archive(str(archive_dir), "zip", args.output, archive_dir.name)
    print(f"  ZIP archive created: {zip_path}")

    print()
    print("=" * 60)
    print("Archive preparation complete!")
    print()
    print("Next steps:")
    print("  1. Create a GitHub release with tag v1.0.0")
    print("  2. Zenodo will automatically archive the release")
    print("  3. Update CITATION.cff with the Zenodo DOI")
    print("  4. Update manuscript with Zenodo DOI reference")
    print("=" * 60)


if __name__ == "__main__":
    main()
