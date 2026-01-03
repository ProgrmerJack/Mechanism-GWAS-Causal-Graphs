#!/usr/bin/env python3
"""
Download and Process Open Targets Genetics L2G Bulk Data

This script:
1. Downloads all L2G parquet files from Open Targets Genetics FTP
2. Filters to benchmark loci based on chromosome:position matching
3. Creates a unified L2G results file with proper gene predictions

Reference: Open Targets Genetics (deprecated July 2025)
FTP: https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g/
"""

import os
import sys
import json
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configuration
BASE_DIR = Path(__file__).parent.parent
L2G_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g"
L2G_BULK_DIR = BASE_DIR / "data" / "external" / "open_targets" / "l2g_bulk"
GOLD_STANDARDS_FILE = BASE_DIR / "data" / "external" / "open_targets" / "curated_gold_standards.tsv"
OUTPUT_DIR = BASE_DIR / "results" / "l2g_bulk_evaluation"

# File pattern for L2G parquet files (200 files: part-00000 to part-00199)
L2G_FILE_PATTERN = "part-{:05d}-38980086-a2d9-4855-88f8-4936870e64a0-c000.snappy.parquet"
N_FILES = 200


def download_file(file_idx: int) -> tuple[int, bool, str]:
    """Download a single L2G parquet file."""
    filename = L2G_FILE_PATTERN.format(file_idx)
    url = f"{L2G_FTP_BASE}/{filename}"
    local_path = L2G_BULK_DIR / f"part-{file_idx:05d}.snappy.parquet"
    
    if local_path.exists():
        return file_idx, True, "already exists"
    
    try:
        urllib.request.urlretrieve(url, local_path)
        return file_idx, True, "downloaded"
    except Exception as e:
        return file_idx, False, str(e)


def download_all_files(max_workers: int = 8, file_indices: list = None) -> int:
    """Download all L2G parquet files in parallel."""
    L2G_BULK_DIR.mkdir(parents=True, exist_ok=True)
    
    if file_indices is None:
        file_indices = list(range(N_FILES))
    
    downloaded = 0
    skipped = 0
    failed = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(download_file, i): i for i in file_indices}
        
        for future in as_completed(futures):
            idx, success, msg = future.result()
            if success:
                if msg == "downloaded":
                    downloaded += 1
                    logger.info(f"Downloaded part-{idx:05d} ({downloaded}/{len(file_indices)})")
                else:
                    skipped += 1
            else:
                failed += 1
                logger.warning(f"Failed part-{idx:05d}: {msg}")
    
    logger.info(f"Download complete: {downloaded} new, {skipped} existing, {failed} failed")
    return downloaded + skipped


def load_gold_standards() -> pd.DataFrame:
    """Load and process gold standard loci."""
    if not GOLD_STANDARDS_FILE.exists():
        raise FileNotFoundError(f"Gold standards file not found: {GOLD_STANDARDS_FILE}")
    
    gs = pd.read_csv(GOLD_STANDARDS_FILE, sep='\t')
    logger.info(f"Loaded {len(gs)} gold standard records")
    
    # Create locus key for matching (chromosome:position)
    gs['chrom'] = gs['chromosome_GRCh38'].astype(str)
    gs['pos'] = gs['position_GRCh38'].astype(int)
    gs['locus_key'] = gs['chrom'] + ':' + gs['pos'].astype(str)
    
    # Get unique loci
    unique_loci = gs.drop_duplicates(subset=['locus_key'])[['locus_key', 'chrom', 'pos', 'gene_id', 'gene_symbol', 'trait_name', 'evidence_class']]
    logger.info(f"Unique loci: {len(unique_loci)}")
    
    return gs, unique_loci


def process_l2g_files(unique_loci: pd.DataFrame, sample_n: int = None) -> pd.DataFrame:
    """Process L2G parquet files and extract predictions for benchmark loci."""
    # Create lookup set for fast matching
    loci_positions = set(zip(unique_loci['chrom'].astype(str), unique_loci['pos']))
    logger.info(f"Looking for {len(loci_positions)} unique loci")
    
    all_matches = []
    files_processed = 0
    
    parquet_files = sorted(L2G_BULK_DIR.glob("*.parquet"))
    if sample_n:
        parquet_files = parquet_files[:sample_n]
    
    logger.info(f"Processing {len(parquet_files)} parquet files...")
    
    for pf in parquet_files:
        try:
            df = pd.read_parquet(pf)
            df['chrom'] = df['chrom'].astype(str)
            
            # Filter to our loci
            mask = df.apply(lambda row: (str(row['chrom']), row['pos']) in loci_positions, axis=1)
            matches = df[mask].copy()
            
            if len(matches) > 0:
                all_matches.append(matches)
                logger.debug(f"{pf.name}: {len(matches)} matches")
            
            files_processed += 1
            if files_processed % 20 == 0:
                logger.info(f"Processed {files_processed}/{len(parquet_files)} files")
                
        except Exception as e:
            logger.warning(f"Error processing {pf}: {e}")
    
    if not all_matches:
        logger.warning("No matches found!")
        return pd.DataFrame()
    
    result = pd.concat(all_matches, ignore_index=True)
    logger.info(f"Total L2G predictions found: {len(result)}")
    
    return result


def evaluate_l2g_predictions(l2g_predictions: pd.DataFrame, gold_standards: pd.DataFrame) -> dict:
    """Evaluate L2G predictions against gold standards."""
    if l2g_predictions.empty:
        return {"coverage": 0, "top1_accuracy": None, "covered_loci": 0, "total_loci": len(gold_standards)}
    
    # Create locus key in predictions
    l2g_predictions['locus_key'] = l2g_predictions['chrom'].astype(str) + ':' + l2g_predictions['pos'].astype(str)
    
    # Get top gene per locus (highest L2G score)
    l2g_top = l2g_predictions.loc[l2g_predictions.groupby('locus_key')['y_proba_full_model'].idxmax()]
    l2g_top = l2g_top[['locus_key', 'gene_id', 'y_proba_full_model', 'study_id']].copy()
    l2g_top.columns = ['locus_key', 'l2g_gene_id', 'l2g_score', 'l2g_study_id']
    
    # Merge with gold standards
    gs_subset = gold_standards.drop_duplicates(subset=['locus_key'])[['locus_key', 'gene_id', 'gene_symbol', 'evidence_class']]
    gs_subset.columns = ['locus_key', 'true_gene_id', 'true_gene_symbol', 'evidence_class']
    
    merged = pd.merge(gs_subset, l2g_top, on='locus_key', how='left')
    
    # Calculate metrics
    covered = merged['l2g_gene_id'].notna()
    n_covered = covered.sum()
    n_total = len(merged)
    coverage = n_covered / n_total
    
    # Top-1 accuracy on covered loci
    if n_covered > 0:
        correct = merged.loc[covered, 'l2g_gene_id'] == merged.loc[covered, 'true_gene_id']
        top1_accuracy = correct.sum() / n_covered
    else:
        top1_accuracy = None
    
    results = {
        "coverage": coverage,
        "covered_loci": n_covered,
        "total_loci": n_total,
        "top1_accuracy": top1_accuracy,
        "correct_predictions": int(correct.sum()) if n_covered > 0 else 0
    }
    
    logger.info(f"L2G Evaluation Results:")
    logger.info(f"  Coverage: {n_covered}/{n_total} = {coverage:.1%}")
    if top1_accuracy is not None:
        logger.info(f"  Top-1 Accuracy: {results['correct_predictions']}/{n_covered} = {top1_accuracy:.1%}")
    
    return results, merged


def main():
    """Main execution flow."""
    import argparse
    parser = argparse.ArgumentParser(description='Download and process L2G bulk data')
    parser.add_argument('--download', action='store_true', help='Download L2G files')
    parser.add_argument('--download-n', type=int, default=None, help='Number of files to download (for testing)')
    parser.add_argument('--process', action='store_true', help='Process downloaded files')
    parser.add_argument('--workers', type=int, default=8, help='Download workers')
    args = parser.parse_args()
    
    if args.download:
        logger.info("=== Downloading L2G bulk data ===")
        file_indices = list(range(args.download_n)) if args.download_n else None
        download_all_files(max_workers=args.workers, file_indices=file_indices)
    
    if args.process:
        logger.info("=== Processing L2G data ===")
        
        # Load gold standards
        gs, unique_loci = load_gold_standards()
        
        # Process L2G files
        l2g_predictions = process_l2g_files(unique_loci)
        
        if not l2g_predictions.empty:
            # Evaluate
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            
            results, merged = evaluate_l2g_predictions(l2g_predictions, gs)
            
            # Save results
            l2g_predictions.to_parquet(OUTPUT_DIR / "l2g_bulk_predictions.parquet", index=False)
            merged.to_csv(OUTPUT_DIR / "l2g_evaluation_merged.csv", index=False)
            
            with open(OUTPUT_DIR / "l2g_bulk_evaluation.json", 'w') as f:
                json.dump(results, f, indent=2)
            
            logger.info(f"Results saved to {OUTPUT_DIR}")
    
    if not args.download and not args.process:
        parser.print_help()


if __name__ == "__main__":
    main()
