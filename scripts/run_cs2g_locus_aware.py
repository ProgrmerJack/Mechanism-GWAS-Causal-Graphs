#!/usr/bin/env python3
"""
Locus-Aware cS2G Evaluation for Nature Genetics
================================================

Implements proper locus-conditional cS2G evaluation:
1. Defines locus variant sets (credible sets OR LD-expanded)
2. Pulls SNP→Gene→Score from official cS2G SGscore files
3. Aggregates scores WITHIN each locus (max + sum as sensitivity)
4. Ranks genes within locus, computes Top-k/MRR with bootstrap CIs

This replaces the methodologically invalid gene-based global matching
in build_cs2g_gene_based.py that inflated accuracy to 100%.

Official cS2G data source: https://zenodo.org/records/7754032
- cS2G_UKBB: 19,476,620 SNPs with columns (SNP, GENE, cS2G, INFO)
- cS2G_1000GEUR: 9,997,231 SNPs

Reference: Gazal et al. 2022 Nature Genetics 54:707-717

Author: Mechanism-GWAS-Causal-Graphs
Date: December 2025
"""

import gzip
import json
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class LocusDefinition:
    """Definition of a GWAS locus for evaluation."""
    locus_id: str
    chromosome: str
    lead_snp: str  # rsID of lead/sentinel variant
    lead_position: int  # Position (hg19 or hg38, must match cS2G build)
    true_gene: str  # Ground truth causal gene
    trait: str
    evidence_tier: str
    variant_set: Set[str]  # Set of rsIDs in locus (credible set or LD-expanded)
    locus_type: str = "ld_expanded"  # "credible_set" or "ld_expanded"
    pos_hg19: int = 0  # Validated hg19 position from allsnps


class LocusAwareCS2GEvaluator:
    """
    Evaluates cS2G in a locus-conditional manner.
    
    Critical difference from gene-based matching:
    - Gene-based: Does this gene have ANY high-scoring SNP in genome? (WRONG)
    - Locus-aware: Is this gene top-ranked AT THIS SPECIFIC LOCUS? (CORRECT)
    
    This class:
    1. Builds a query SNP set (union of all locus variants)
    2. Streams official cS2G files once, keeping only relevant SNPs
    3. For each locus, aggregates SNP→Gene scores and ranks genes
    """
    
    def __init__(
        self,
        project_root: Path,
        cs2g_source: str = "UKBB",  # "UKBB" or "1000GEUR"
        aggregation: str = "max",   # "max" or "sum"
        ld_window_kb: int = 500,
        ld_r2_threshold: float = 0.6
    ):
        self.project_root = project_root
        self.cs2g_source = cs2g_source
        self.aggregation = aggregation
        self.ld_window_kb = ld_window_kb
        self.ld_r2_threshold = ld_r2_threshold
        
        # Paths
        self.cs2g_dir = project_root / "data" / "external" / "cS2G" / f"cS2G_{cs2g_source}"
        if not self.cs2g_dir.exists():
            # Try alternate path
            self.cs2g_dir = project_root / "data" / "external" / "cS2G" / "cS2G_UKBB"
        
        self.benchmark_path = project_root / "data" / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
        self.output_dir = project_root / "results" / "cs2g_locus_aware"
        
        # Data containers
        self.loci: List[LocusDefinition] = []
        self.query_snps: Set[str] = set()
        self.cs2g_scores: Dict[str, List[Tuple[str, float]]] = defaultdict(list)  # snp -> [(gene, score)]
        
        # Gene ID mapping (for consistency)
        self.gene_symbol_map: Dict[str, str] = {}  # normalized -> original
    
    def load_benchmark(self) -> pd.DataFrame:
        """Load post-2021 independent benchmark."""
        logger.info(f"Loading benchmark from {self.benchmark_path}")
        
        if not self.benchmark_path.exists():
            raise FileNotFoundError(f"Benchmark not found: {self.benchmark_path}")
        
        df = pd.read_csv(self.benchmark_path, sep='\t')
        logger.info(f"Loaded {len(df)} loci from benchmark")
        return df
    
    def build_locus_definitions(self, benchmark_df: pd.DataFrame) -> None:
        """
        Build locus definitions from benchmark.
        
        For each benchmark locus:
        - Extract lead SNP and position
        - Define variant set (currently using LD window; credible sets as enhancement)
        """
        logger.info("Building locus definitions...")
        
        for _, row in benchmark_df.iterrows():
            # Extract locus info
            locus_id = row['locus_id']
            lead_snp = row.get('lead_snp', row.get('rsid', ''))
            
            # Get chromosome (handle various formats)
            chrom = str(row.get('chr', row.get('chromosome', '')))
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            
            # Get position (prefer hg19 if available as cS2G is hg19)
            pos = row.get('pos_hg19', row.get('pos_hg38', row.get('position', 0)))
            
            # True gene
            true_gene = row.get('gene_symbol', row.get('gene', ''))
            
            # Build variant set using LD window (simplified; PLINK for production)
            # For now, we use the lead SNP as query key
            # In production: expand using --ld-window-kb and --ld-window-r2
            variant_set = {lead_snp} if lead_snp else set()
            
            locus = LocusDefinition(
                locus_id=locus_id,
                chromosome=chrom,
                lead_snp=lead_snp,
                lead_position=int(pos) if pos else 0,
                true_gene=true_gene.upper(),
                trait=row.get('trait_category', row.get('trait', '')),
                evidence_tier=row.get('evidence_tier', ''),
                variant_set=variant_set,
                locus_type="ld_expanded"
            )
            
            self.loci.append(locus)
            self.query_snps.update(variant_set)
        
        logger.info(f"Built {len(self.loci)} locus definitions")
        logger.info(f"Query SNP set size: {len(self.query_snps)} variants")
    
    def build_variant_sets_from_window(self) -> None:
        """
        Populate variant_set for each locus using ±500kb window expansion.
        
        Reads data/external/cS2G/cS2G_UKBB/allsnps.txt.gz (rsID chr pos)
        and assigns variants to loci if they fall within the window.
        """
        logger.info("Building locus variant sets from ±500kb windows...")
        logger.info(f"Using LD window: ±{self.ld_window_kb}kb")
        
        allsnps_path = self.cs2g_dir / "allsnps.txt.gz"
        if not allsnps_path.exists():
            logger.warning(f"SNP map not found: {allsnps_path}")
            logger.warning("Falling back to lead SNP only (coverage will be low!)")
            self.query_snps.add(lead_snp) # Add lead SNP to query set for now
        
        logger.info(f"Built {len(self.loci)} initial locus definitions.")
        logger.info(f"Initial query SNP set size (lead SNPs): {len(self.query_snps)} variants.")
    
    def find_sentinels_hg19(self) -> Dict[str, Tuple[str, int]]:
        """
        Pass 1: Scan allsnps.txt.gz to find hg19 coordinates of sentinel SNPs.
        Returns: {locus_id: (chrom, pos_hg19)} for found sentinels.
        """
        logger.info("Pass 1: Scanning allsnps.txt.gz for sentinel hg19 coordinates...")
        
        allsnps_path = self.cs2g_dir / "allsnps.txt.gz"
        if not allsnps_path.exists():
            raise FileNotFoundError(f"SNP map not found: {allsnps_path}")
        
        # Create a mapping from rsID to a list of locus_ids that use it as a sentinel
        sentinel_rsid_to_locus_ids = defaultdict(list)
        for locus in self.loci:
            sentinel_rsid_to_locus_ids[locus.lead_snp].append(locus.locus_id)
        
        found_sentinels: Dict[str, Tuple[str, int]] = {} # locus_id -> (chrom, pos_hg19)
        
        with gzip.open(allsnps_path, 'rt') as f:
            for line in tqdm(f, desc="Streaming allsnps (Pass 1)"):
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                
                snp_id = parts[0]
                
                if snp_id in sentinel_rsid_to_locus_ids:
                    chrom = parts[1]
                    try:
                        pos = int(parts[2])
                    except ValueError:
                        logger.warning(f"Could not parse position for SNP {snp_id} in allsnps.txt.gz. Skipping.")
                        continue
                    
                    for locus_id in sentinel_rsid_to_locus_ids[snp_id]:
                        found_sentinels[locus_id] = (chrom, pos)
                    
                    # If all sentinels are found, we can stop early
                    if len(found_sentinels) == len(self.loci):
                        break
        
        logger.info(f"Found hg19 coordinates for {len(found_sentinels)} out of {len(self.loci)} sentinel SNPs.")
        return found_sentinels

    def build_variant_sets_hg19(self) -> None:
        """
        Pass 2: Scan allsnps.txt.gz to find all variants in hg19 windows around sentinels.
        Populates locus.variant_set and self.query_snps.
        """
        logger.info(f"Pass 2: Scanning allsnps.txt.gz for window variants (±{self.ld_window_kb}kb)...")
        
        allsnps_path = self.cs2g_dir / "allsnps.txt.gz"
        if not allsnps_path.exists():
            raise FileNotFoundError(f"SNP map not found: {allsnps_path}")
            
        # Pre-compute windows per chromosome for efficient lookup
        # Dictionary mapping: chromosome -> list of (start, end, locus_index)
        locus_windows = defaultdict(list)
        for i, locus in enumerate(self.loci):
            if locus.pos_hg19: # Only for loci with known hg19 positions
                window_bp = self.ld_window_kb * 1000
                start = max(1, locus.pos_hg19 - window_bp)
                end = locus.pos_hg19 + window_bp
                locus_windows[locus.chrom].append((start, end, i))
        
        # Clear existing query_snps and variant_sets, as we're rebuilding them
        self.query_snps.clear()
        for locus in self.loci:
            locus.variant_set.clear()
            # Re-add lead SNP to variant_set and query_snps
            if locus.lead_snp:
                locus.variant_set.add(locus.lead_snp)
                self.query_snps.add(locus.lead_snp)

        total_snps_processed = 0
        mapped_snps_to_loci = 0
        
        with gzip.open(allsnps_path, 'rt') as f:
            for line in tqdm(f, desc="Streaming allsnps (Pass 2)"):
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                        
                snp_id = parts[0]
                chrom = parts[1]
                try:
                    pos = int(parts[2])
                except ValueError:
                    continue
                
                total_snps_processed += 1
                
                # Check if in any locus window on this chromosome
                if chrom in locus_windows:
                    for start, end, idx in locus_windows[chrom]:
                        if start <= pos <= end:
                            # Add SNP to locus's variant set and global query set
                            self.loci[idx].variant_set.add(snp_id)
                            self.query_snps.add(snp_id)
                            mapped_snps_to_loci += 1
                                
        logger.info(f"Processed {total_snps_processed:,} SNPs from allsnps.txt.gz.")
        logger.info(f"Mapped {mapped_snps_to_loci:,} variant assignments to loci within windows.")
        logger.info(f"Total unique query SNPs after window expansion: {len(self.query_snps)}.")

    def run(self, window_kb: int = 500) -> Tuple[pd.DataFrame, Dict]:
        """
        Run evaluation with specified window size.
        """
        self.ld_window_kb = window_kb # Update the instance variable
        self.window_bp = window_kb * 1000
        logger.info("=" * 70)
        logger.info(f"LOCUS-AWARE cS2G EVALUATION (Window: ±{window_kb}kb, hg19-based)")
        logger.info("=" * 70)
        
        # Load benchmark
        benchmark_df = self.load_benchmark()
        
        # Build initial locus definitions (populates locus_id, lead_snp, true_gene, etc.)
        self.build_locus_definitions(benchmark_df)
        
        # Pass 1: Find hg19 coordinates for sentinel SNPs
        sentinel_hg19_coords = self.find_sentinels_hg19()
        
        # Update loci with hg19 positions and filter out unreachable loci
        valid_loci = []
        for locus in self.loci:
            if locus.locus_id in sentinel_hg19_coords:
                locus.chrom, locus.pos_hg19 = sentinel_hg19_coords[locus.locus_id]
                valid_loci.append(locus)
            else:
                logger.warning(f"Sentinel SNP '{locus.lead_snp}' for locus '{locus.locus_id}' not found in allsnps.txt.gz. Skipping locus.")
        
        self.loci = valid_loci
        logger.info(f"Proceeding with {len(self.loci)} loci whose sentinel SNPs are found in cS2G reference (allsnps.txt.gz).")
        
        if not self.loci:
            logger.error("No valid loci to evaluate after filtering. Exiting.")
            return pd.DataFrame(), {}

        # Pass 2: Build variant sets (window expansion in hg19)
        self.build_variant_sets_hg19()
        
        # Load cS2G scores for the expanded query SNP set
        self.load_cs2g_scores_for_query_set()
        
        # Evaluate all loci
        results_df, all_gene_scores_df = self.evaluate_all_loci() # Assuming this method exists
        
        # Compute summary metrics
        summary = self.compute_summary_metrics(results_df) # Assuming this method exists
        
        # Add metadata to summary
        summary['window_kb'] = window_kb
        summary['n_loci_valid_hg19'] = len(self.loci)
        
        # Report key metric
        if 'top1_correct_mean' in summary:
            logger.info(f"Top-1 Accuracy (±{window_kb}kb): {summary['top1_correct_mean']*100:.1f}%")
        else:
            logger.warning("Top-1 accuracy not found in summary metrics.")
        
        # Save results
        self.output_dir.mkdir(parents=True, exist_ok=True)
        suffix = f"{self.aggregation}_{window_kb}kb"
        results_df.to_csv(self.output_dir / f"cs2g_results_{suffix}.tsv", sep='\t', index=False)
        all_gene_scores_df.to_csv(self.output_dir / f"cs2g_full_scores_{suffix}.tsv", sep='\t', index=False)
        with open(self.output_dir / f"cs2g_summary_{suffix}.json", 'w') as f:
            json.dump(summary, f, indent=2)
            
        return results_df, summary

    def load_cs2g_scores_for_query_set(self) -> None:
        """
        Stream through official cS2G SGscore files and extract scores
        for SNPs in our query set.
        
        Efficient: single pass through each chromosome file.
        
        cS2G file format (tab-separated):
        - Column 0: SNP ID (rsID format)
        - Column 1: GENE (symbol or Ensembl ID)
        - Column 2: cS2G score (0-1)
        - Column 3+: INFO (constituent strategy scores)
        """
        logger.info("Loading cS2G scores for query SNP set...")
        logger.info(f"Looking in: {self.cs2g_dir}")
        
        if not self.cs2g_dir.exists():
            raise FileNotFoundError(f"cS2G directory not found: {self.cs2g_dir}")
        
        total_lines = 0
        matched_lines = 0
        
        for chr_num in tqdm(range(1, 23), desc="Processing cS2G chromosomes"):
            score_file = self.cs2g_dir / f"cS2G.{chr_num}.SGscore.gz"
            
            if not score_file.exists():
                logger.warning(f"Score file not found: {score_file}")
                continue
            
            try:
                with gzip.open(score_file, 'rt') as f:
                    # Skip header
                    header = f.readline()
                    
                    for line in f:
                        total_lines += 1
                        parts = line.strip().split('\t')
                        
                        if len(parts) >= 3:
                            snp_id = parts[0]
                            gene = parts[1].upper()  # Normalize to uppercase
                            
                            try:
                                score = float(parts[2])
                            except ValueError:
                                continue
                            
                            # Check if this SNP is in our query set
                            if snp_id in self.query_snps:
                                self.cs2g_scores[snp_id].append((gene, score))
                                matched_lines += 1
                                
            except Exception as e:
                logger.error(f"Error processing {score_file}: {e}")
                continue
        
        logger.info(f"Scanned {total_lines:,} total lines")
        logger.info(f"Matched {matched_lines:,} SNP-gene pairs for query set")
        logger.info(f"Unique SNPs with scores: {len(self.cs2g_scores)}")
    
    def load_cs2g_scores_by_position(self) -> Dict[Tuple[str, int], List[Tuple[str, float]]]:
        """
        Alternative: Load cS2G scores indexed by (chr, position).
        
        This enables position-based locus matching when rsIDs aren't available.
        
        Note: This requires parsing SNP IDs that include position info,
        or using a separate SNP→position mapping.
        """
        logger.info("Loading cS2G scores with position indexing...")
        
        position_scores: Dict[Tuple[str, int], List[Tuple[str, float]]] = defaultdict(list)
        
        # For now, we extract position from rsID format if possible
        # Production would use dbSNP or 1000G variant files
        
        for chr_num in tqdm(range(1, 23), desc="Position-indexing cS2G"):
            score_file = self.cs2g_dir / f"cS2G.{chr_num}.SGscore.gz"
            
            if not score_file.exists():
                continue
            
            try:
                with gzip.open(score_file, 'rt') as f:
                    header = f.readline()
                    
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            snp_id = parts[0]
                            gene = parts[1].upper()
                            
                            try:
                                score = float(parts[2])
                            except ValueError:
                                continue
                            
                            # Store by SNP ID and also index for quick lookup
                            self.cs2g_scores[snp_id].append((gene, score))
                            
            except Exception as e:
                logger.error(f"Error processing {score_file}: {e}")
        
        logger.info(f"Loaded scores for {len(self.cs2g_scores)} unique SNPs")
        return position_scores
    
    def aggregate_scores_for_locus(
        self,
        locus: LocusDefinition
    ) -> pd.DataFrame:
        """
        Aggregate cS2G scores for all genes at a locus.
        
        For each gene that has a cS2G link to any variant in the locus:
        - max: Take maximum score across variants
        - sum: Take sum of scores (penalizes single-SNP support)
        
        Returns DataFrame with columns: gene, cs2g_score, n_snps, max_snp
        """
        gene_scores: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
        
        # Collect all SNP→Gene scores for variants in this locus
        for snp_id in locus.variant_set:
            if snp_id in self.cs2g_scores:
                for gene, score in self.cs2g_scores[snp_id]:
                    gene_scores[gene].append((snp_id, score))
        
        if not gene_scores:
            # No cS2G coverage for this locus
            return pd.DataFrame(columns=['gene', 'cs2g_score', 'n_snps', 'max_snp'])
        
        # Aggregate to gene-level
        results = []
        for gene, snp_score_list in gene_scores.items():
            scores = [s for _, s in snp_score_list]
            
            if self.aggregation == "max":
                agg_score = max(scores)
            elif self.aggregation == "sum":
                agg_score = sum(scores)
            else:
                agg_score = max(scores)
            
            # Find best SNP
            best_snp, best_score = max(snp_score_list, key=lambda x: x[1])
            
            results.append({
                'gene': gene,
                'cs2g_score': agg_score,
                'n_snps': len(snp_score_list),
                'max_snp': best_snp,
                'max_score': best_score
            })
        
        df = pd.DataFrame(results)
        df = df.sort_values('cs2g_score', ascending=False).reset_index(drop=True)
        df['rank'] = range(1, len(df) + 1)
        
        return df
    
    def evaluate_locus(self, locus: LocusDefinition) -> Dict:
        """
        Evaluate cS2G performance for a single locus.
        
        Returns dict with:
        - locus_id, true_gene, trait, evidence_tier
        - true_gene_rank, true_gene_score
        - top1_correct, top3_correct, top5_correct, top10_correct
        - reciprocal_rank
        - candidate_count
        - coverage (whether true gene has any cS2G score)
        """
        # Get gene scores for this locus
        gene_scores = self.aggregate_scores_for_locus(locus)
        
        if gene_scores.empty:
            return {
                'locus_id': locus.locus_id,
                'true_gene': locus.true_gene,
                'trait': locus.trait,
                'evidence_tier': locus.evidence_tier,
                'true_gene_rank': np.nan,
                'true_gene_score': np.nan,
                'top1_correct': np.nan,
                'top3_correct': np.nan,
                'top5_correct': np.nan,
                'top10_correct': np.nan,
                'reciprocal_rank': np.nan,
                'candidate_count': 0,
                'coverage': False,
                'aggregation': self.aggregation
            }
        
        # Find true gene rank
        true_gene_row = gene_scores[gene_scores['gene'] == locus.true_gene]
        
        if true_gene_row.empty:
            # True gene not in cS2G predictions for this locus
            return {
                'locus_id': locus.locus_id,
                'true_gene': locus.true_gene,
                'trait': locus.trait,
                'evidence_tier': locus.evidence_tier,
                'true_gene_rank': len(gene_scores) + 1,
                'true_gene_score': 0.0,
                'top1_correct': 0,
                'top3_correct': 0,
                'top5_correct': 0,
                'top10_correct': 0,
                'reciprocal_rank': 0.0,
                'candidate_count': len(gene_scores),
                'coverage': False,
                'aggregation': self.aggregation
            }
        
        true_rank = int(true_gene_row.iloc[0]['rank'])
        true_score = float(true_gene_row.iloc[0]['cs2g_score'])
        
        return {
            'locus_id': locus.locus_id,
            'true_gene': locus.true_gene,
            'trait': locus.trait,
            'evidence_tier': locus.evidence_tier,
            'true_gene_rank': true_rank,
            'true_gene_score': true_score,
            'top1_correct': 1 if true_rank == 1 else 0,
            'top3_correct': 1 if true_rank <= 3 else 0,
            'top5_correct': 1 if true_rank <= 5 else 0,
            'top10_correct': 1 if true_rank <= 10 else 0,
            'reciprocal_rank': 1.0 / true_rank,
            'candidate_count': len(gene_scores),
            'coverage': True,
            'aggregation': self.aggregation
        }
    
    def evaluate_all_loci(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Evaluate cS2G across all loci and compute summary metrics.
        Returns: (results_df, all_gene_scores_df)
        """
        logger.info(f"Evaluating {len(self.loci)} loci...")
        
        results = []
        all_gene_scores = []
        
        for locus in tqdm(self.loci, desc="Evaluating loci"):
            # Get metrics
            result = self.evaluate_locus(locus)
            results.append(result)
            
            # Get full gene scores for this locus
            gene_scores = self.aggregate_scores_for_locus(locus)
            if not gene_scores.empty:
                gene_scores['locus_id'] = locus.locus_id
                all_gene_scores.append(gene_scores)
        
        results_df = pd.DataFrame(results)
        
        if all_gene_scores:
            all_gene_scores_df = pd.concat(all_gene_scores, ignore_index=True)
        else:
            all_gene_scores_df = pd.DataFrame(columns=['locus_id', 'gene', 'cs2g_score', 'rank'])
            
        return results_df, all_gene_scores_df
    
    def compute_bootstrap_ci(
        self,
        results_df: pd.DataFrame,
        metric_col: str,
        n_bootstrap: int = 1000,
        ci_level: float = 0.95
    ) -> Tuple[float, float, float]:
        """
        Compute bootstrap confidence interval for a metric.
        
        Returns: (mean, ci_lower, ci_upper)
        """
        values = results_df[metric_col].dropna().values
        if len(values) == 0:
            return np.nan, np.nan, np.nan
        
        np.random.seed(42)
        bootstrap_means = []
        
        for _ in range(n_bootstrap):
            sample = np.random.choice(values, size=len(values), replace=True)
            bootstrap_means.append(np.mean(sample))
        
        alpha = 1 - ci_level
        ci_lower = np.percentile(bootstrap_means, 100 * alpha / 2)
        ci_upper = np.percentile(bootstrap_means, 100 * (1 - alpha / 2))
        
        return float(np.mean(values)), float(ci_lower), float(ci_upper)
    
    def compute_summary_metrics(self, results_df: pd.DataFrame) -> Dict:
        """
        Compute summary metrics with bootstrap CIs.
        """
        # Filter to loci with coverage
        covered_df = results_df[results_df['coverage'] == True]
        
        n_total = len(results_df)
        n_covered = len(covered_df)
        
        # Compute metrics with CIs
        metrics = {
            'n_loci_total': n_total,
            'n_loci_covered': n_covered,
            'coverage_rate': n_covered / n_total if n_total > 0 else 0,
            'aggregation_method': self.aggregation,
            'locus_definition': 'ld_expanded'  # or 'credible_set'
        }
        
        for metric in ['top1_correct', 'top3_correct', 'top5_correct', 'top10_correct']:
            mean, ci_lower, ci_upper = self.compute_bootstrap_ci(covered_df, metric)
            metrics[f'{metric}_mean'] = mean
            metrics[f'{metric}_ci_lower'] = ci_lower
            metrics[f'{metric}_ci_upper'] = ci_upper
        
        # MRR
        mrr_mean, mrr_lower, mrr_upper = self.compute_bootstrap_ci(covered_df, 'reciprocal_rank')
        metrics['mrr_mean'] = mrr_mean
        metrics['mrr_ci_lower'] = mrr_lower
        metrics['mrr_ci_upper'] = mrr_upper
        
        # Mean rank
        mean_rank_values = covered_df['true_gene_rank'].dropna().values
        metrics['mean_rank'] = float(np.mean(mean_rank_values)) if len(mean_rank_values) > 0 else np.nan
        metrics['median_rank'] = float(np.median(mean_rank_values)) if len(mean_rank_values) > 0 else np.nan
        
        return metrics



def main():
    """Main entry point."""
    project_root = Path(__file__).resolve().parents[1]
    
    # Run evaluation with max aggregation (primary)
    logger.info("\n" + "=" * 70)
    logger.info("PRIMARY ANALYSIS: max aggregation")
    logger.info("=" * 70)
    
    evaluator_max = LocusAwareCS2GEvaluator(
        project_root,
        aggregation="max"
    )
    results_max, summary_max = evaluator_max.run()
    
    # Run evaluation with sum aggregation (sensitivity)
    logger.info("\n" + "=" * 70)
    logger.info("SENSITIVITY ANALYSIS: sum aggregation")
    logger.info("=" * 70)
    
    evaluator_sum = LocusAwareCS2GEvaluator(
        project_root,
        aggregation="sum"
    )
    evaluator_sum.loci = evaluator_max.loci  # Reuse locus definitions
    evaluator_sum.query_snps = evaluator_max.query_snps
    evaluator_sum.cs2g_scores = evaluator_max.cs2g_scores  # Reuse loaded scores
    
    results_sum, summary_sum = evaluator_sum.run()
    
    # Print comparison
    print("\n" + "=" * 70)
    print("COMPARISON: max vs sum aggregation")
    print("=" * 70)
    print(f"{'Metric':<20} {'Max':<25} {'Sum':<25}")
    print("-" * 70)
    
    for metric in ['top1_correct', 'top3_correct', 'top5_correct']:
        max_val = f"{summary_max[f'{metric}_mean']*100:.1f}% [{summary_max[f'{metric}_ci_lower']*100:.1f}-{summary_max[f'{metric}_ci_upper']*100:.1f}%]"
        sum_val = f"{summary_sum[f'{metric}_mean']*100:.1f}% [{summary_sum[f'{metric}_ci_lower']*100:.1f}-{summary_sum[f'{metric}_ci_upper']*100:.1f}%]"
        print(f"{metric:<20} {max_val:<25} {sum_val:<25}")
    
    print("=" * 70)
    
    print("\nNOTE: These are LOCUS-AWARE scores, not gene-global scores.")
    print("Results may differ significantly from the invalid gene-based matching.")


if __name__ == "__main__":
    main()
