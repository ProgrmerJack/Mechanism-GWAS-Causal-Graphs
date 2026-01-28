#!/usr/bin/env python3
"""
Build cS2G Benchmark with Proper Coordinate Liftover
=====================================================

The cS2G scores from Gazal et al. 2022 use hg19/GRCh37 coordinates,
but our post-2021 benchmark uses hg38/GRCh38 coordinates.

This script:
1. Uses pyliftover to convert hg38 → hg19 coordinates
2. Matches benchmark loci to cS2G pre-computed scores
3. Evaluates cS2G accuracy on the benchmark

Dependencies:
    pip install pyliftover pandas tqdm

Author: Mechanism-GWAS-Causal-Graphs
"""

import gzip
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Check for pyliftover
try:
    from pyliftover import LiftOver
    LIFTOVER_AVAILABLE = True
except ImportError:
    LIFTOVER_AVAILABLE = False
    logger.warning("pyliftover not installed. Will install it now...")


def install_pyliftover():
    """Install pyliftover package."""
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyliftover"])
    global LiftOver, LIFTOVER_AVAILABLE
    from pyliftover import LiftOver
    LIFTOVER_AVAILABLE = True


class CS2GLiftoverMatcher:
    """Match benchmark loci to cS2G scores using coordinate liftover."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.data_dir = project_root / "data"
        
        # Initialize liftover
        logger.info("Initializing hg38 → hg19 liftover...")
        self.liftover = LiftOver('hg38', 'hg19')
        
        # Paths
        self.cs2g_dir = self.data_dir / "external" / "cS2G" / "cS2G_extracted" / "cS2G_UKBB"
        self.benchmark_path = self.data_dir / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
        self.output_dir = self.data_dir / "processed" / "baselines"
        
    def liftover_coordinate(self, chrom: str, pos: int) -> Optional[Tuple[str, int]]:
        """
        Convert hg38 coordinate to hg19.
        
        Args:
            chrom: Chromosome (e.g., '1', 'chr1')
            pos: Position in hg38 (1-based)
        
        Returns:
            Tuple of (chrom_hg19, pos_hg19) or None if unmappable
        """
        # Ensure chromosome has 'chr' prefix for liftover
        if not str(chrom).startswith('chr'):
            chrom_str = f"chr{chrom}"
        else:
            chrom_str = str(chrom)
        
        # pyliftover uses 0-based coordinates
        result = self.liftover.convert_coordinate(chrom_str, pos - 1)
        
        if result and len(result) > 0:
            # Take the first result (primary mapping)
            new_chrom, new_pos, _, _ = result[0]
            # Return without 'chr' prefix and 1-based position
            new_chrom_clean = new_chrom.replace('chr', '')
            return (new_chrom_clean, int(new_pos) + 1)
        
        return None
    
    def load_cs2g_scores(self) -> Tuple[pd.DataFrame, Dict[str, str]]:
        """
        Load cS2G scores and rsID mapping from extracted files.
        
        Returns:
            Tuple of (scores_df, rsid_to_hg19_position_map)
        """
        logger.info("Loading cS2G scores...")
        
        if not self.cs2g_dir.exists():
            raise FileNotFoundError(f"cS2G directory not found: {self.cs2g_dir}")
        
        # Load allsnps.txt.gz for rsID to position mapping
        allsnps_path = self.cs2g_dir / "allsnps.txt.gz"
        logger.info("Loading rsID → hg19 position mapping...")
        
        rsid_to_pos = {}  # rsID → "chr:pos"
        pos_to_rsid = {}  # "chr:pos" → rsID (for debugging)
        
        with gzip.open(allsnps_path, 'rt') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    rsid = parts[0]
                    chr_num = parts[1]
                    pos_hg19 = parts[2]
                    key = f"{chr_num}:{pos_hg19}"
                    rsid_to_pos[rsid] = key
                    pos_to_rsid[key] = rsid
        
        logger.info(f"Loaded {len(rsid_to_pos):,} rsID mappings")
        
        # Load scores for all chromosomes
        all_scores = []
        for chr_num in range(1, 23):
            score_file = self.cs2g_dir / f"cS2G.{chr_num}.SGscore.gz"
            if score_file.exists():
                with gzip.open(score_file, 'rt') as f:
                    header = f.readline()  # Skip header
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            snp_id = parts[0]  # Format: CHR:POS_REF_ALT
                            gene = parts[1]
                            score = float(parts[2])
                            
                            # Parse position from SNP ID
                            if ':' in snp_id and '_' in snp_id:
                                chr_pos = snp_id.split('_')[0]  # CHR:POS
                                chrom, pos = chr_pos.split(':')
                                
                                all_scores.append({
                                    'snp_id_full': snp_id,
                                    'chr': int(chrom),
                                    'pos_hg19': int(pos),
                                    'chr_pos_hg19': chr_pos,
                                    'gene': gene,
                                    'cs2g_score': score
                                })
        
        scores_df = pd.DataFrame(all_scores)
        logger.info(f"Loaded {len(scores_df):,} cS2G score entries across all chromosomes")
        
        return scores_df, rsid_to_pos
    
    def load_benchmark(self) -> pd.DataFrame:
        """Load post-2021 independent benchmark."""
        logger.info(f"Loading benchmark from {self.benchmark_path}")
        benchmark = pd.read_csv(self.benchmark_path, sep='\t')
        logger.info(f"Benchmark: {len(benchmark)} loci")
        return benchmark
    
    def match_benchmark_to_cs2g(
        self, 
        benchmark: pd.DataFrame, 
        cs2g_scores: pd.DataFrame,
        rsid_map: Dict[str, str],
        position_window: int = 10  # Allow small position differences
    ) -> pd.DataFrame:
        """
        Match benchmark loci to cS2G scores using liftover coordinates.
        
        Strategy:
        1. Lift over hg38 → hg19 position
        2. Look up cS2G score at lifted position for the truth gene
        3. Also try rsID matching as fallback
        """
        logger.info("Matching benchmark loci to cS2G scores...")
        
        results = []
        
        for idx, row in tqdm(benchmark.iterrows(), total=len(benchmark), desc="Matching loci"):
            locus_id = row['locus_id']
            gene = row['gene_symbol']
            chr_hg38 = row['chr']
            pos_hg38 = row['pos_hg38']
            rsid = row.get('lead_snp', '')
            
            cs2g_score = None
            match_method = None
            pos_hg19 = None
            
            # Strategy 1: Liftover hg38 → hg19
            try:
                lifted = self.liftover_coordinate(str(chr_hg38), int(pos_hg38))
            except Exception as e:
                logger.debug(f"Liftover failed for {chr_hg38}:{pos_hg38}: {e}")
                lifted = None
            
            if lifted:
                lifted_chr, lifted_pos = lifted
                pos_hg19 = lifted_pos
                
                # Skip non-autosomal chromosomes (cS2G only has 1-22)
                try:
                    lifted_chr_int = int(lifted_chr)
                except ValueError:
                    logger.debug(f"Skipping non-autosomal chromosome: {lifted_chr}")
                    lifted = None
                
            if lifted:
                # Search for cS2G score at lifted position
                # Allow small window due to potential liftover imprecision
                chr_matches = cs2g_scores[
                    (cs2g_scores['chr'] == lifted_chr_int) &
                    (cs2g_scores['gene'].str.upper() == gene.upper()) &
                    (abs(cs2g_scores['pos_hg19'] - lifted_pos) <= position_window)
                ]
                
                if len(chr_matches) > 0:
                    # Take closest position if multiple matches
                    chr_matches = chr_matches.copy()
                    chr_matches['pos_diff'] = abs(chr_matches['pos_hg19'] - lifted_pos)
                    closest = chr_matches.loc[chr_matches['pos_diff'].idxmin()]
                    cs2g_score = closest['cs2g_score']
                    match_method = 'liftover'
                    pos_hg19 = closest['pos_hg19']
            
            # Strategy 2: rsID matching as fallback
            if cs2g_score is None and rsid and rsid.startswith('rs'):
                if rsid in rsid_map:
                    chr_pos_hg19 = rsid_map[rsid]
                    chr_num, pos_str = chr_pos_hg19.split(':')
                    rsid_pos = int(pos_str)
                    
                    # Search for score at this position
                    rsid_matches = cs2g_scores[
                        (cs2g_scores['chr'] == int(chr_num)) &
                        (cs2g_scores['gene'].str.upper() == gene.upper()) &
                        (cs2g_scores['pos_hg19'] == rsid_pos)
                    ]
                    
                    if len(rsid_matches) > 0:
                        cs2g_score = rsid_matches['cs2g_score'].max()
                        match_method = 'rsid'
                        pos_hg19 = rsid_pos
            
            # Record result
            prediction = 1 if cs2g_score is not None and cs2g_score > 0.5 else 0
            correct = 1 if cs2g_score is not None and cs2g_score > 0.5 else None  # All benchmark genes are true positives
            
            results.append({
                'locus_id': locus_id,
                'gene_symbol': gene,
                'chr': chr_hg38,
                'pos_hg38': pos_hg38,
                'pos_hg19': pos_hg19,
                'lead_snp': rsid,
                'cs2g_score': cs2g_score,
                'cs2g_prediction': prediction,
                'match_method': match_method,
                'correct': correct
            })
        
        return pd.DataFrame(results)
    
    def calculate_accuracy(self, results: pd.DataFrame) -> Dict:
        """Calculate accuracy metrics for cS2G baseline."""
        # Coverage: loci with cS2G scores
        matched = results['cs2g_score'].notna().sum()
        total = len(results)
        coverage = matched / total
        
        # Accuracy: among matched loci, how many predict the truth gene?
        # (All benchmark genes are true positives by definition)
        matched_df = results[results['cs2g_score'].notna()]
        if len(matched_df) > 0:
            correct = (matched_df['cs2g_prediction'] == 1).sum()
            accuracy = correct / len(matched_df)
        else:
            correct = 0
            accuracy = 0.0
        
        # Match method breakdown
        method_counts = results['match_method'].value_counts().to_dict()
        
        summary = {
            'total_loci': total,
            'matched_loci': int(matched),
            'coverage': float(coverage),
            'correct_predictions': int(correct),
            'accuracy_on_matched': float(accuracy),
            'match_methods': method_counts
        }
        
        return summary
    
    def run(self) -> Tuple[pd.DataFrame, Dict]:
        """Run full cS2G matching pipeline."""
        logger.info("=" * 70)
        logger.info("cS2G BENCHMARK MATCHING WITH COORDINATE LIFTOVER")
        logger.info("=" * 70)
        
        # Load data
        cs2g_scores, rsid_map = self.load_cs2g_scores()
        benchmark = self.load_benchmark()
        
        # Match benchmark to cS2G
        results = self.match_benchmark_to_cs2g(benchmark, cs2g_scores, rsid_map)
        
        # Calculate accuracy
        summary = self.calculate_accuracy(results)
        
        # Report
        logger.info("=" * 70)
        logger.info("RESULTS:")
        logger.info(f"  Total loci: {summary['total_loci']}")
        logger.info(f"  Matched loci: {summary['matched_loci']} ({summary['coverage']:.1%} coverage)")
        logger.info(f"  Correct predictions: {summary['correct_predictions']}")
        logger.info(f"  Accuracy (on matched): {summary['accuracy_on_matched']:.1%}")
        logger.info(f"  Match methods: {summary['match_methods']}")
        logger.info("=" * 70)
        
        # Save results
        output_path = self.output_dir / "cs2g_benchmark_liftover.tsv"
        results.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved results to {output_path}")
        
        summary_path = self.output_dir / "cs2g_benchmark_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Saved summary to {summary_path}")
        
        return results, summary


def main():
    """Main entry point."""
    # Install pyliftover if needed
    if not LIFTOVER_AVAILABLE:
        install_pyliftover()
    
    # Find project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    # Run matching
    matcher = CS2GLiftoverMatcher(project_root)
    results, summary = matcher.run()
    
    print("\n" + "=" * 70)
    print("SUMMARY FOR MANUSCRIPT:")
    print("=" * 70)
    print(f"cS2G (Gazal et al. 2022):")
    print(f"  Coverage: {summary['matched_loci']}/{summary['total_loci']} ({summary['coverage']:.1%})")
    print(f"  Accuracy: {summary['accuracy_on_matched']:.1%}")
    print("=" * 70)


if __name__ == "__main__":
    main()
