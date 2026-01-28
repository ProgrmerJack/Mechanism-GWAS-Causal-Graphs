#!/usr/bin/env python3
"""
COMPREHENSIVE BASELINE COMPARISON
=================================

This script extracts baseline scores for all post-2021 benchmark loci:
1. Open Targets L2G (via Platform API v25.4.4)
2. cS2G (Combined SNP-to-Gene strategy)

Key considerations:
- Platform API uses GRCh38 coordinates
- cS2G uses hg19/GRCh37 coordinates
- Benchmark data is in GRCh38

For coordinate liftover, we use NCBI remapping service or local chain files.
"""

import requests
import pandas as pd
import numpy as np
import json
import gzip
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from io import StringIO
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# API endpoints
PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

class BaselineExtractor:
    """Extract baseline method scores for benchmark evaluation."""
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = Path(data_dir)
        self.output_dir = self.data_dir / "processed" / "baselines"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = []
        
        # Load benchmark
        benchmark_path = self.output_dir / "post2021_independent_benchmark_FINAL.tsv"
        self.benchmark = pd.read_csv(benchmark_path, sep='\t')
        logger.info(f"Loaded benchmark with {len(self.benchmark)} loci")
        
        # Cache for API queries
        self.l2g_cache = {}
        
    def query_l2g_for_variant(self, chr_num: str, pos: int) -> List[Dict]:
        """
        Query Open Targets L2G predictions for a genomic position.
        
        Since we don't have exact variant IDs, we search nearby credible sets.
        """
        cache_key = f"{chr_num}_{pos}"
        if cache_key in self.l2g_cache:
            return self.l2g_cache[cache_key]
        
        # Query credible sets that might contain this position
        # We use a region-based approach
        query = """
        query getCredibleSetsNearPosition($chr: String!) {
          credibleSets(
            chromosome: $chr,
            page: {size: 1000, index: 0}
          ) {
            count
            rows {
              studyLocusId
              chromosome
              position
              region
              study {
                id
                traitFromSource
              }
              locus {
                variantId
                posteriorProbability
              }
              l2GPredictions {
                rows {
                  score
                  target {
                    id
                    approvedSymbol
                  }
                }
              }
            }
          }
        }
        """
        
        try:
            r = requests.post(PLATFORM_API, json={
                "query": query,
                "variables": {"chr": str(chr_num)}
            }, timeout=120)
            data = r.json()
            
            if 'data' in data and data['data'].get('credibleSets'):
                rows = data['data']['credibleSets'].get('rows', [])
                
                # Filter to credible sets near our position (within 1Mb)
                nearby = []
                for row in rows:
                    cs_pos = row.get('position', 0)
                    if abs(cs_pos - pos) < 1_000_000:  # Within 1Mb
                        nearby.append(row)
                
                self.l2g_cache[cache_key] = nearby
                return nearby
        except Exception as e:
            logger.error(f"Error querying L2G for chr{chr_num}:{pos}: {e}")
        
        return []
    
    def query_l2g_by_rsid(self, rsid: str, gene_symbol: str) -> Optional[float]:
        """
        Query L2G score for a specific rsID and gene combination.
        """
        # Search for the variant using its rsID
        search_query = """
        query searchVariant($rsid: String!) {
          search(queryString: $rsid, entityNames: ["variant"], page: {size: 5, index: 0}) {
            hits {
              id
              name
              entity
            }
          }
        }
        """
        
        try:
            r = requests.post(PLATFORM_API, json={
                "query": search_query,
                "variables": {"rsid": rsid}
            }, timeout=30)
            data = r.json()
            
            if 'data' in data and data['data'].get('search'):
                hits = data['data']['search'].get('hits', [])
                for hit in hits:
                    if hit.get('entity') == 'variant':
                        variant_id = hit.get('id')
                        
                        # Get detailed variant info with credible sets
                        var_query = """
                        query getVariantL2G($variantId: String!) {
                          variant(variantId: $variantId) {
                            id
                            rsIds
                            credibleSets(page: {size: 100, index: 0}) {
                              rows {
                                l2GPredictions {
                                  rows {
                                    score
                                    target {
                                      approvedSymbol
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                        """
                        
                        r2 = requests.post(PLATFORM_API, json={
                            "query": var_query,
                            "variables": {"variantId": variant_id}
                        }, timeout=30)
                        data2 = r2.json()
                        
                        if 'data' in data2 and data2['data'].get('variant'):
                            v = data2['data']['variant']
                            for cs in v.get('credibleSets', {}).get('rows', []):
                                for l2g in cs.get('l2GPredictions', {}).get('rows', []):
                                    if l2g['target']['approvedSymbol'] == gene_symbol:
                                        return l2g['score']
        except Exception as e:
            logger.error(f"Error querying L2G for {rsid}/{gene_symbol}: {e}")
        
        return None
    
    def query_l2g_direct(self, chr_num: str, pos: int, gene_symbol: str) -> Optional[float]:
        """
        Direct L2G query using variant position and gene symbol.
        
        This uses a more efficient approach - querying by gene and filtering.
        """
        # First, get gene's Ensembl ID
        gene_query = """
        query getGene($gene: String!) {
          search(queryString: $gene, entityNames: ["target"], page: {size: 5, index: 0}) {
            hits {
              id
              name
              entity
            }
          }
        }
        """
        
        try:
            r = requests.post(PLATFORM_API, json={
                "query": gene_query,
                "variables": {"gene": gene_symbol}
            }, timeout=30)
            data = r.json()
            
            gene_id = None
            if 'data' in data and data['data'].get('search'):
                for hit in data['data']['search'].get('hits', []):
                    if hit.get('entity') == 'target' and hit.get('name', '').upper() == gene_symbol.upper():
                        gene_id = hit.get('id')
                        break
            
            if not gene_id:
                return None
            
            # Now query credible sets for this chromosome and look for our gene
            cs_query = """
            query getChromosomeCredibleSets($chr: String!, $geneId: String!) {
              credibleSets(
                chromosome: $chr,
                page: {size: 500, index: 0}
              ) {
                rows {
                  position
                  l2GPredictions {
                    rows {
                      score
                      target {
                        id
                        approvedSymbol
                      }
                    }
                  }
                }
              }
            }
            """
            
            r2 = requests.post(PLATFORM_API, json={
                "query": cs_query,
                "variables": {"chr": str(chr_num), "geneId": gene_id}
            }, timeout=120)
            data2 = r2.json()
            
            if 'data' in data2 and data2['data'].get('credibleSets'):
                best_score = None
                min_dist = float('inf')
                
                for cs in data2['data']['credibleSets'].get('rows', []):
                    cs_pos = cs.get('position', 0)
                    dist = abs(cs_pos - pos)
                    
                    if dist < 500_000:  # Within 500kb
                        for l2g in cs.get('l2GPredictions', {}).get('rows', []):
                            if l2g['target']['approvedSymbol'] == gene_symbol:
                                # Keep the score from the nearest credible set
                                if dist < min_dist:
                                    min_dist = dist
                                    best_score = l2g['score']
                
                return best_score
                
        except Exception as e:
            logger.error(f"Error in direct L2G query for chr{chr_num}:{pos}/{gene_symbol}: {e}")
        
        return None
    
    def extract_l2g_scores(self) -> pd.DataFrame:
        """Extract L2G scores for all benchmark loci."""
        logger.info("Extracting L2G scores for benchmark loci...")
        
        l2g_results = []
        
        for idx, row in self.benchmark.iterrows():
            locus_id = row['locus_id']
            chr_num = str(row['chr']).replace('chr', '')
            pos = row['pos_hg38']
            gene = row['gene_symbol']
            rsid = row.get('lead_snp', '')
            
            logger.info(f"[{idx+1}/{len(self.benchmark)}] {locus_id}: chr{chr_num}:{pos} - {gene}")
            
            # Try multiple approaches
            l2g_score = None
            
            # Approach 1: Direct position + gene query
            l2g_score = self.query_l2g_direct(chr_num, pos, gene)
            
            # Approach 2: Try with rsID if available
            if l2g_score is None and rsid and rsid.startswith('rs'):
                l2g_score = self.query_l2g_by_rsid(rsid, gene)
            
            l2g_results.append({
                'locus_id': locus_id,
                'chr': chr_num,
                'pos_hg38': pos,
                'gene_symbol': gene,
                'lead_snp': rsid,
                'l2g_score': l2g_score,
                'l2g_prediction': 1 if l2g_score and l2g_score > 0.5 else 0
            })
            
            # Rate limiting
            time.sleep(0.5)
        
        df = pd.DataFrame(l2g_results)
        
        # Save results
        output_path = self.output_dir / "l2g_benchmark_scores.tsv"
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved L2G scores to {output_path}")
        
        return df
    
    def load_cs2g_scores(self) -> pd.DataFrame:
        """
        Load cS2G scores from extracted files.
        
        Note: cS2G uses hg19 coordinates. We need to either:
        1. Lift over our benchmark to hg19
        2. Use rsID matching
        """
        logger.info("Loading cS2G scores...")
        
        cs2g_dir = self.data_dir / "external" / "cS2G" / "cS2G_extracted" / "cS2G_UKBB"
        
        if not cs2g_dir.exists():
            logger.warning(f"cS2G directory not found: {cs2g_dir}")
            return pd.DataFrame()
        
        # Load allsnps.txt.gz for rsID to position mapping
        allsnps_path = cs2g_dir / "allsnps.txt.gz"
        if not allsnps_path.exists():
            logger.warning("allsnps.txt.gz not found")
            return pd.DataFrame()
        
        logger.info("Loading cS2G rsID mapping...")
        # File format: rsid chr pos (space-separated)
        # e.g., "rs367896724 1 10177"
        with gzip.open(allsnps_path, 'rt') as f:
            rsid_map = {}
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    rsid = parts[0]
                    chr_num = parts[1]
                    pos = parts[2]
                    # Create SNP ID in cS2G format: CHR:POS_REF_ALT
                    # But we don't have ref/alt here, so store chr:pos for partial matching
                    rsid_map[rsid] = f"{chr_num}:{pos}"
        
        logger.info(f"Loaded {len(rsid_map)} rsID mappings")
        
        # Now load scores for each chromosome
        all_scores = []
        for chr_num in range(1, 23):
            score_file = cs2g_dir / f"cS2G.{chr_num}.SGscore.gz"
            if score_file.exists():
                with gzip.open(score_file, 'rt') as f:
                    header = f.readline().strip().split('\t')
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            all_scores.append({
                                'snp_hg19': parts[0],
                                'gene': parts[1],
                                'cs2g_score': float(parts[2]),
                                'chr': chr_num
                            })
        
        scores_df = pd.DataFrame(all_scores)
        logger.info(f"Loaded {len(scores_df)} cS2G scores")
        
        return scores_df, rsid_map
    
    def match_cs2g_to_benchmark(self, cs2g_scores: pd.DataFrame, rsid_map: Dict) -> pd.DataFrame:
        """Match cS2G scores to benchmark loci using rsID and position."""
        logger.info("Matching cS2G scores to benchmark...")
        
        cs2g_results = []
        
        for idx, row in self.benchmark.iterrows():
            locus_id = row['locus_id']
            gene = row['gene_symbol']
            rsid = row.get('lead_snp', '')
            
            cs2g_score = None
            matched_snp = None
            
            if rsid and rsid in rsid_map:
                # rsid_map value is "chr:pos" in hg19
                chr_pos_hg19 = rsid_map[rsid]
                chr_num, pos_hg19 = chr_pos_hg19.split(':')
                
                # Look for any SNP at this position (any allele combination)
                # cS2G format in scores file: CHR:POS_REF_ALT
                prefix = f"{chr_num}:{pos_hg19}_"
                
                # Find matching scores
                matches = cs2g_scores[
                    (cs2g_scores['snp_hg19'].str.startswith(prefix)) &
                    (cs2g_scores['gene'] == gene)
                ]
                
                if len(matches) > 0:
                    cs2g_score = matches['cs2g_score'].max()  # Take highest score if multiple alleles
                    matched_snp = matches['snp_hg19'].values[0]
            
            cs2g_results.append({
                'locus_id': locus_id,
                'gene_symbol': gene,
                'lead_snp': rsid,
                'cs2g_snp_hg19': matched_snp,
                'cs2g_score': cs2g_score,
                'cs2g_prediction': 1 if cs2g_score and cs2g_score > 0.5 else 0
            })
        
        df = pd.DataFrame(cs2g_results)
        
        # Save results
        output_path = self.output_dir / "cs2g_benchmark_scores.tsv"
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved cS2G scores to {output_path}")
        
        return df
    
    def run_extraction(self):
        """Run full baseline extraction pipeline."""
        logger.info("=" * 70)
        logger.info("BASELINE EXTRACTION FOR POST-2021 BENCHMARK")
        logger.info("=" * 70)
        
        # Extract L2G scores
        l2g_df = self.extract_l2g_scores()
        
        # Extract cS2G scores
        try:
            cs2g_scores, rsid_map = self.load_cs2g_scores()
            if len(cs2g_scores) > 0:
                cs2g_df = self.match_cs2g_to_benchmark(cs2g_scores, rsid_map)
            else:
                cs2g_df = pd.DataFrame()
        except Exception as e:
            logger.error(f"Error loading cS2G: {e}")
            cs2g_df = pd.DataFrame()
        
        # Merge results
        results = self.benchmark[['locus_id', 'chr', 'pos_hg38', 'gene_symbol', 'lead_snp', 
                                   'trait', 'evidence_tier', 'validation_type']].copy()
        
        if len(l2g_df) > 0:
            results = results.merge(l2g_df[['locus_id', 'l2g_score', 'l2g_prediction']], 
                                   on='locus_id', how='left')
        
        if len(cs2g_df) > 0:
            results = results.merge(cs2g_df[['locus_id', 'cs2g_score', 'cs2g_prediction']], 
                                   on='locus_id', how='left')
        
        # Save combined results
        output_path = self.output_dir / "benchmark_baseline_scores_combined.tsv"
        results.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved combined results to {output_path}")
        
        # Print summary
        logger.info("\n" + "=" * 70)
        logger.info("SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Total benchmark loci: {len(self.benchmark)}")
        
        if 'l2g_score' in results.columns:
            l2g_coverage = results['l2g_score'].notna().sum()
            l2g_correct = (results['l2g_prediction'] == 1).sum()
            logger.info(f"\nL2G Results:")
            logger.info(f"  Coverage: {l2g_coverage}/{len(results)} ({100*l2g_coverage/len(results):.1f}%)")
            logger.info(f"  Correct (>0.5): {l2g_correct}/{l2g_coverage} ({100*l2g_correct/max(1,l2g_coverage):.1f}%)")
        
        if 'cs2g_score' in results.columns:
            cs2g_coverage = results['cs2g_score'].notna().sum()
            cs2g_correct = (results['cs2g_prediction'] == 1).sum()
            logger.info(f"\ncS2G Results:")
            logger.info(f"  Coverage: {cs2g_coverage}/{len(results)} ({100*cs2g_coverage/len(results):.1f}%)")
            logger.info(f"  Correct (>0.5): {cs2g_correct}/{cs2g_coverage} ({100*cs2g_correct/max(1,cs2g_coverage):.1f}%)")
        
        return results


def main():
    import os
    os.chdir(r"C:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs")
    
    extractor = BaselineExtractor()
    results = extractor.run_extraction()
    
    print("\n\nFinal Results Preview:")
    print(results.head(10).to_string())

if __name__ == "__main__":
    main()
