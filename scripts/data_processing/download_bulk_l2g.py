#!/usr/bin/env python3
"""
download_bulk_l2g.py
====================
Download and process bulk L2G (Locus-to-Gene) data from Open Targets.

Open Targets provides bulk data via FTP/HTTPS:
https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/

This script:
1. Downloads L2G parquet files from Open Targets
2. Processes and indexes them by variant position
3. Creates a lookup table for benchmark matching

Author: Mechanism-GWAS-Causal-Graphs
Date: 2025
"""

import os
import sys
import json
import requests
import pandas as pd
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import pyarrow.parquet as pq
import hashlib

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
L2G_DIR = DATA_DIR / "external" / "opentargets_l2g"
L2G_RAW_DIR = L2G_DIR / "raw_parquet"
L2G_PROCESSED_DIR = L2G_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

# Ensure directories exist
L2G_RAW_DIR.mkdir(parents=True, exist_ok=True)
L2G_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Open Targets FTP base URL
OT_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics"
# Latest release with L2G data (symlink to 22.09)
L2G_RELEASE = "latest"


def get_l2g_file_list():
    """Get list of L2G parquet files from Open Targets FTP."""
    # L2G files are typically in the l2g/ directory
    l2g_index_url = f"{OT_FTP_BASE}/{L2G_RELEASE}/l2g/"
    
    print(f"Fetching L2G file list from: {l2g_index_url}")
    
    try:
        response = requests.get(l2g_index_url, timeout=30)
        response.raise_for_status()
        
        # Parse HTML to find parquet files
        import re
        files = re.findall(r'href="([^"]+\.parquet)"', response.text)
        
        print(f"Found {len(files)} L2G parquet files")
        return [f"{l2g_index_url}{f}" for f in files]
    except Exception as e:
        print(f"Error fetching file list: {e}")
        # Fallback: try to construct file names
        return [f"{l2g_index_url}part-{i:05d}.parquet" for i in range(200)]


def download_file(url, dest_path, max_retries=3):
    """Download a single file with retry logic."""
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=120, stream=True)
            response.raise_for_status()
            
            with open(dest_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            return True
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"  Retry {attempt + 1} for {url}: {e}")
            else:
                print(f"  Failed to download {url}: {e}")
                return False
    return False


def download_l2g_data(max_files=None, force_redownload=False):
    """Download all L2G parquet files."""
    file_urls = get_l2g_file_list()
    
    if max_files:
        file_urls = file_urls[:max_files]
    
    print(f"\nDownloading {len(file_urls)} L2G files...")
    
    downloaded = 0
    skipped = 0
    failed = 0
    
    for i, url in enumerate(file_urls):
        filename = url.split('/')[-1]
        dest_path = L2G_RAW_DIR / filename
        
        if dest_path.exists() and not force_redownload:
            skipped += 1
            continue
        
        print(f"  [{i+1}/{len(file_urls)}] Downloading {filename}...")
        
        if download_file(url, dest_path):
            downloaded += 1
        else:
            failed += 1
    
    print(f"\nDownload complete: {downloaded} downloaded, {skipped} skipped, {failed} failed")
    return downloaded + skipped


def process_l2g_parquet_files():
    """Process all downloaded L2G parquet files into a unified dataset."""
    print("\n" + "=" * 70)
    print("PROCESSING L2G PARQUET FILES")
    print("=" * 70)
    
    parquet_files = list(L2G_RAW_DIR.glob("*.parquet"))
    print(f"Found {len(parquet_files)} parquet files to process")
    
    if not parquet_files:
        print("No parquet files found. Attempting download...")
        download_l2g_data(max_files=10)  # Download first 10 for testing
        parquet_files = list(L2G_RAW_DIR.glob("*.parquet"))
    
    # Read and concatenate all parquet files
    all_dfs = []
    total_rows = 0
    
    for i, pf in enumerate(parquet_files):
        try:
            df = pd.read_parquet(pf)
            all_dfs.append(df)
            total_rows += len(df)
            
            if (i + 1) % 20 == 0:
                print(f"  Processed {i + 1}/{len(parquet_files)} files ({total_rows:,} rows)")
        except Exception as e:
            print(f"  Error reading {pf.name}: {e}")
    
    if not all_dfs:
        print("No valid parquet files found!")
        return None
    
    # Concatenate all dataframes
    print(f"\nConcatenating {len(all_dfs)} dataframes ({total_rows:,} total rows)...")
    l2g_full = pd.concat(all_dfs, ignore_index=True)
    
    print(f"Combined L2G dataset: {len(l2g_full):,} rows")
    print(f"Columns: {list(l2g_full.columns)}")
    
    # Clean up and process
    print("\nProcessing L2G data...")
    
    # Standardize column names
    column_mapping = {
        'y_proba_full_model': 'l2g_score',
        'y_proba_loco_model': 'l2g_score_loco',
        'gene_id': 'ensembl_gene_id'
    }
    
    for old_name, new_name in column_mapping.items():
        if old_name in l2g_full.columns:
            l2g_full = l2g_full.rename(columns={old_name: new_name})
    
    # Create variant ID for matching
    if all(col in l2g_full.columns for col in ['chrom', 'pos', 'ref', 'alt']):
        l2g_full['variant_id'] = (
            l2g_full['chrom'].astype(str) + '_' + 
            l2g_full['pos'].astype(str) + '_' + 
            l2g_full['ref'].astype(str) + '_' + 
            l2g_full['alt'].astype(str)
        )
    
    # Get best L2G gene per study-variant combination
    print("\nExtracting top L2G gene per locus...")
    
    if 'study_id' in l2g_full.columns and 'l2g_score' in l2g_full.columns:
        # Sort by score and get top gene per study-variant
        l2g_full_sorted = l2g_full.sort_values('l2g_score', ascending=False)
        
        # Group by study + variant to get top gene
        grouping_cols = ['study_id', 'variant_id'] if 'variant_id' in l2g_full.columns else ['study_id', 'chrom', 'pos']
        
        l2g_top = l2g_full_sorted.groupby(grouping_cols).first().reset_index()
        print(f"Top L2G genes per locus: {len(l2g_top):,} entries")
    else:
        l2g_top = l2g_full
    
    # Save processed data
    output_path = L2G_PROCESSED_DIR / "l2g_processed.parquet"
    l2g_full.to_parquet(output_path, index=False)
    print(f"\nSaved full L2G data to: {output_path}")
    
    top_output_path = L2G_PROCESSED_DIR / "l2g_top_genes.parquet"
    l2g_top.to_parquet(top_output_path, index=False)
    print(f"Saved top L2G genes to: {top_output_path}")
    
    # Summary statistics
    print("\n" + "-" * 50)
    print("L2G DATA SUMMARY")
    print("-" * 50)
    
    if 'study_id' in l2g_full.columns:
        print(f"Unique studies: {l2g_full['study_id'].nunique():,}")
    if 'variant_id' in l2g_full.columns:
        print(f"Unique variants: {l2g_full['variant_id'].nunique():,}")
    if 'ensembl_gene_id' in l2g_full.columns:
        print(f"Unique genes: {l2g_full['ensembl_gene_id'].nunique():,}")
    if 'l2g_score' in l2g_full.columns:
        print(f"L2G score range: {l2g_full['l2g_score'].min():.3f} - {l2g_full['l2g_score'].max():.3f}")
        print(f"L2G score mean: {l2g_full['l2g_score'].mean():.3f}")
    
    return l2g_full


def match_l2g_to_benchmark():
    """Match L2G data to benchmark loci."""
    print("\n" + "=" * 70)
    print("MATCHING L2G TO BENCHMARK LOCI")
    print("=" * 70)
    
    # Load benchmark
    benchmark_path = DATA_DIR / "processed" / "baselines" / "mechanism_stratified_bench_v1.tsv"
    if not benchmark_path.exists():
        benchmark_path = DATA_DIR / "external" / "gold_standards" / "curated_gold_standards.tsv"
    
    if not benchmark_path.exists():
        print(f"Benchmark file not found: {benchmark_path}")
        return None
    
    benchmark = pd.read_csv(benchmark_path, sep='\t')
    print(f"Loaded benchmark with {len(benchmark)} loci")
    
    # Load L2G data
    l2g_path = L2G_PROCESSED_DIR / "l2g_processed.parquet"
    if not l2g_path.exists():
        print("L2G processed data not found. Processing...")
        process_l2g_parquet_files()
    
    if not l2g_path.exists():
        print("Could not create L2G processed data")
        return None
    
    l2g = pd.read_parquet(l2g_path)
    print(f"Loaded L2G data with {len(l2g):,} entries")
    
    # Create matching keys
    # Try multiple matching strategies
    
    matched_results = []
    
    for idx, row in benchmark.iterrows():
        locus_id = row.get('locus_id', f"locus_{idx}")
        gold_gene = row.get('gold_gene', row.get('gene_symbol', row.get('gene', None)))
        trait = row.get('trait', row.get('disease', None))
        
        # Try to find matching L2G entries
        # Strategy 1: Match by rsID if available
        rsid = row.get('sentinel_rsid', row.get('rsid', None))
        
        # Strategy 2: Match by chromosome position
        chrom = row.get('chromosome', row.get('chrom', None))
        pos = row.get('position', row.get('pos', None))
        
        l2g_match = None
        match_method = None
        
        # Try position-based matching
        if chrom is not None and pos is not None:
            chrom_str = str(chrom).replace('chr', '')
            mask = (l2g['chrom'].astype(str) == chrom_str) & (l2g['pos'] == int(pos))
            if mask.any():
                l2g_match = l2g[mask].nlargest(1, 'l2g_score').iloc[0] if 'l2g_score' in l2g.columns else l2g[mask].iloc[0]
                match_method = 'exact_position'
        
        # Try window-based matching if exact match failed
        if l2g_match is None and chrom is not None and pos is not None:
            window = 500000  # 500kb window
            chrom_str = str(chrom).replace('chr', '')
            mask = (
                (l2g['chrom'].astype(str) == chrom_str) & 
                (l2g['pos'] >= int(pos) - window) & 
                (l2g['pos'] <= int(pos) + window)
            )
            if mask.any():
                l2g_match = l2g[mask].nlargest(1, 'l2g_score').iloc[0] if 'l2g_score' in l2g.columns else l2g[mask].iloc[0]
                match_method = 'window_500kb'
        
        # Record result
        result = {
            'locus_id': locus_id,
            'gold_gene': gold_gene,
            'trait': trait,
            'chromosome': chrom,
            'position': pos,
            'rsid': rsid,
            'l2g_available': l2g_match is not None,
            'match_method': match_method
        }
        
        if l2g_match is not None:
            result['l2g_top_gene'] = l2g_match.get('ensembl_gene_id', l2g_match.get('gene_id', None))
            result['l2g_score'] = l2g_match.get('l2g_score', l2g_match.get('y_proba_full_model', None))
            result['l2g_study_id'] = l2g_match.get('study_id', None)
        
        matched_results.append(result)
    
    matched_df = pd.DataFrame(matched_results)
    
    # Summary
    coverage = matched_df['l2g_available'].sum() / len(matched_df) * 100
    print(f"\nL2G Coverage: {matched_df['l2g_available'].sum()}/{len(matched_df)} loci ({coverage:.1f}%)")
    
    # Save results
    output_path = RESULTS_DIR / "l2g_benchmark_matching.parquet"
    matched_df.to_parquet(output_path, index=False)
    print(f"Saved matching results to: {output_path}")
    
    # Also save as TSV for inspection
    tsv_path = RESULTS_DIR / "l2g_benchmark_matching.tsv"
    matched_df.to_csv(tsv_path, sep='\t', index=False)
    print(f"Saved TSV to: {tsv_path}")
    
    return matched_df


def main():
    """Main execution function."""
    print("=" * 70)
    print("OPEN TARGETS L2G BULK DATA DOWNLOAD & PROCESSING")
    print("=" * 70)
    print(f"Started at: {datetime.now().isoformat()}")
    print(f"Release: {L2G_RELEASE}")
    print(f"Output directory: {L2G_DIR}")
    
    # Step 1: Download L2G parquet files
    print("\n" + "=" * 70)
    print("STEP 1: DOWNLOADING L2G DATA")
    print("=" * 70)
    
    # First, try to download a few files to test
    download_l2g_data(max_files=20)
    
    # Step 2: Process parquet files
    l2g_data = process_l2g_parquet_files()
    
    # Step 3: Match to benchmark
    if l2g_data is not None:
        match_l2g_to_benchmark()
    
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"Finished at: {datetime.now().isoformat()}")


if __name__ == "__main__":
    main()
