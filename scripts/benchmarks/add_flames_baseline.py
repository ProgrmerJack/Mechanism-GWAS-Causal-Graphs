#!/usr/bin/env python3
"""
FLAMES Baseline Evaluation for RegulatoryBench.

This script integrates FLAMES (Functional Linking and Mapping of Enhancers to Genes and 
Genes to Phenotypes) scores as a baseline for the RegulatoryBench benchmark.

FLAMES paper: https://www.nature.com/articles/s41588-024-02051-2
Data: Zenodo 12635505 (Annotation_data.tar.gz)

Usage:
    python add_flames_baseline.py --extract  # First: extract the archive
    python add_flames_baseline.py --evaluate # Then: run evaluation
"""

import argparse
import subprocess
import tarfile
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score

BASE_DIR = Path(__file__).parent.parent
FLAMES_DIR = BASE_DIR / "data" / "external" / "flames"
ARCHIVE_FILE = FLAMES_DIR / "Annotation_data.tar.gz"
CANDIDATES_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_cs2g.tsv"
OUTPUT_FILE = BASE_DIR / "data" / "processed" / "baselines" / "evaluation_candidates_with_all.tsv"


def extract_archive():
    """Extract the FLAMES annotation data archive."""
    print("=" * 70)
    print("Extracting FLAMES Annotation Data")
    print("=" * 70)
    
    if not ARCHIVE_FILE.exists():
        print(f"ERROR: Archive not found: {ARCHIVE_FILE}")
        print("Please download from: https://zenodo.org/records/12635505")
        return False
    
    print(f"Archive: {ARCHIVE_FILE}")
    print(f"Size: {ARCHIVE_FILE.stat().st_size / 1e9:.2f} GB")
    
    print("\nExtracting... (this may take a few minutes)")
    
    with tarfile.open(ARCHIVE_FILE, 'r:gz') as tar:
        tar.extractall(path=FLAMES_DIR)
    
    print("Extraction complete!")
    
    # List extracted contents
    print("\nExtracted files:")
    for item in FLAMES_DIR.iterdir():
        if item.is_dir():
            n_files = len(list(item.rglob('*')))
            print(f"  {item.name}/ ({n_files} files)")
        else:
            print(f"  {item.name}")
    
    return True


def find_flames_scores():
    """Locate FLAMES E2G score files."""
    
    # FLAMES provides enhancer-to-gene scores
    # Look for common file patterns
    potential_patterns = [
        "Annotation_data/**/E2G*.tsv*",
        "Annotation_data/**/enhancer*.tsv*",
        "Annotation_data/**/*scores*.tsv*",
        "**/E2G*.tsv*",
        "**/*E2G*.bed*",
    ]
    
    for pattern in potential_patterns:
        matches = list(FLAMES_DIR.glob(pattern))
        if matches:
            print(f"\nFound files matching '{pattern}':")
            for m in matches[:10]:  # Show first 10
                print(f"  {m}")
            return matches
    
    # If no patterns match, list directory structure
    print("\nNo E2G files found. Directory structure:")
    for item in FLAMES_DIR.rglob('*'):
        if item.is_file():
            rel_path = item.relative_to(FLAMES_DIR)
            print(f"  {rel_path}")
    
    return []


def load_flames_e2g_scores(score_files):
    """Load and merge FLAMES E2G scores."""
    
    print(f"\nLoading FLAMES scores from {len(score_files)} files...")
    
    all_scores = []
    for f in score_files:
        try:
            # Try different separators and headers
            df = pd.read_csv(f, sep='\t', comment='#')
            all_scores.append(df)
            print(f"  Loaded {f.name}: {len(df)} entries")
        except Exception as e:
            print(f"  Warning: Could not load {f.name}: {e}")
    
    if not all_scores:
        return None
    
    combined = pd.concat(all_scores, ignore_index=True)
    print(f"\nTotal FLAMES entries: {len(combined):,}")
    print(f"Columns: {combined.columns.tolist()}")
    
    return combined


def evaluate_flames_baseline():
    """Evaluate FLAMES as a baseline on RegulatoryBench."""
    
    print("=" * 70)
    print("Evaluating FLAMES Baseline")
    print("=" * 70)
    
    # Load candidates
    if not CANDIDATES_FILE.exists():
        print(f"ERROR: Candidates file not found: {CANDIDATES_FILE}")
        return
    
    print(f"\nLoading candidates...")
    candidates = pd.read_csv(CANDIDATES_FILE, sep='\t', low_memory=False)
    print(f"  Loaded {len(candidates):,} candidates")
    
    # Find FLAMES score files
    score_files = find_flames_scores()
    
    if not score_files:
        print("\nERROR: No FLAMES score files found.")
        print("Please first extract the archive with: python add_flames_baseline.py --extract")
        return
    
    # Load FLAMES scores
    flames_df = load_flames_e2g_scores(score_files)
    
    if flames_df is None:
        print("ERROR: Could not load FLAMES scores")
        return
    
    # The exact matching logic depends on FLAMES data structure
    # Common approaches:
    # 1. Match by enhancer position (chr, start, end) + gene
    # 2. Match by gene symbol/ENSEMBL + distance
    
    print("\nFLAMES data sample:")
    print(flames_df.head())
    
    # TODO: Implement matching logic based on actual data structure
    # This is a placeholder - actual implementation depends on FLAMES format
    
    print("\n" + "=" * 70)
    print("FLAMES integration requires examining the actual data structure")
    print("Once extracted, inspect files to determine matching strategy")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="FLAMES baseline for RegulatoryBench")
    parser.add_argument('--extract', action='store_true', help='Extract archive')
    parser.add_argument('--evaluate', action='store_true', help='Run evaluation')
    parser.add_argument('--list', action='store_true', help='List archive contents without extracting')
    
    args = parser.parse_args()
    
    if args.list:
        if ARCHIVE_FILE.exists():
            print("Archive contents (first 50 entries):")
            with tarfile.open(ARCHIVE_FILE, 'r:gz') as tar:
                for i, member in enumerate(tar.getmembers()):
                    if i >= 50:
                        print("  ... (truncated)")
                        break
                    print(f"  {member.name}")
        else:
            print(f"Archive not found: {ARCHIVE_FILE}")
        return
    
    if args.extract:
        extract_archive()
    elif args.evaluate:
        evaluate_flames_baseline()
    else:
        print("Usage:")
        print("  python add_flames_baseline.py --extract   # Extract archive")
        print("  python add_flames_baseline.py --evaluate  # Run evaluation")
        print("  python add_flames_baseline.py --list      # List archive contents")


if __name__ == "__main__":
    main()
