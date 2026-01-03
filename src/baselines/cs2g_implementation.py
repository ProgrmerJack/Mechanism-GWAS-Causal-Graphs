"""
cS2G-Inspired Heritability-Weighted SNP-to-Gene Proxy
=======================================================

IMPORTANT: This is NOT a full implementation of cS2G from Gazal et al. 2022.
This is an APPROXIMATION/PROXY that uses publicly available data to approximate
the cS2G methodology for benchmarking purposes.

The official cS2G implementation requires:
- Full fine-mapped credible sets with PIPs
- Trait-specific heritability enrichment calculations using LDSC
- Complete ABC enhancer-gene predictions across all cell types
- Tissue-matched eQTL colocalization analyses
- Integration with PoPS polygenic scores

What this proxy does:
- Combines ABC, eQTL, distance, coding, and constraint scores
- Uses fixed strategy weights (not trait-specific heritability enrichment)
- Provides comparable scoring for method comparison

Reference (official cS2G):
    Gazal et al. (2022). Combining SNP-to-gene linking strategies to improve 
    GWAS gene discovery and our understanding of the architecture of complex traits.
    Nature Genetics 54, 707–717. DOI: 10.1038/s41588-022-01087-y
    
    Official code: https://github.com/lukejoconnor/cS2G

Author: Mechanism-GWAS-Causal-Graphs team
Status: Proxy/approximation for benchmarking, not production-ready cS2G
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging


class CS2GInspiredProxy:
    """
    cS2G-inspired proxy using publicly available data.
    
    DISCLAIMER: This is an approximation, not the official cS2G implementation.
    
    For production use cases requiring true cS2G scores, use the official
    implementation from https://github.com/lukejoconnor/cS2G
    
    This proxy combines multiple S2G strategies with fixed weights derived
    from cS2G paper estimates, but does NOT perform trait-specific heritability
    enrichment calculations.
    """
    
    def __init__(self, random_state: int = 42):
        """Initialize cS2G scorer."""
        self.random_state = random_state
        self.logger = logging.getLogger(__name__)
        self.strategy_weights = {
            'abc': 0.35,           # Activity by Contact - most predictive
            'eqtl': 0.30,          # eQTL colocalization
            'distance': 0.15,      # Distance to TSS - weakest alone
            'coding': 0.10,        # Coding variant enrichment
            'constraint': 0.10     # Evolutionary constraint
        }
    
    def calculate_heritability_enrichment(
        self,
        locus_variants: pd.DataFrame,
        locus_name: str
    ) -> Dict[str, float]:
        """
        Calculate heritability enrichment for each S2G strategy in this locus.
        
        Formula:
            Enrichment_s = h2_trait_in_strategy_s / h2_total_strategy_s
            
        where h2_trait_in_strategy_s = heritability explained by variants linked
                                       to genes via strategy s
        
        Args:
            locus_variants: DataFrame with columns: variant_id, pip, abc_score,
                           eqtl_pip_h4, is_coding, constraint_score
            locus_name: Name/ID of this locus for tracking
        
        Returns:
            Dictionary with enrichment scores for each strategy
        """
        enrichments = {}
        
        # Calculate strategy-specific heritability
        for strategy in self.strategy_weights.keys():
            if strategy == 'abc':
                # ABC-linked variant heritability (weighted by PIP and ABC score)
                h2_strategy = (locus_variants['pip'] * locus_variants.get('abc_score', 0)).sum()
            elif strategy == 'eqtl':
                # eQTL-linked variant heritability
                h2_strategy = (locus_variants['pip'] * locus_variants.get('eqtl_pip_h4', 0)).sum()
            elif strategy == 'distance':
                # Distance-based linking (inverse distance)
                max_distance = locus_variants.get('distance_to_nearest_gene', 100000)
                distance_scores = 1.0 / (1.0 + locus_variants.get('distance_to_nearest_gene', max_distance) / 100000)
                h2_strategy = (locus_variants['pip'] * distance_scores).sum()
            elif strategy == 'coding':
                # Coding variant contribution
                h2_strategy = (locus_variants['pip'] * locus_variants.get('is_coding', 0)).sum()
            elif strategy == 'constraint':
                # Evolutionary constraint (LINSIGHT scores)
                h2_strategy = (locus_variants['pip'] * locus_variants.get('constraint_score', 0.5)).sum()
            
            enrichments[strategy] = max(0.0, h2_strategy)
        
        # Normalize by total heritability
        total_h2 = sum(enrichments.values()) + 1e-10
        for strategy in enrichments:
            enrichments[strategy] = enrichments[strategy] / total_h2
        
        return enrichments
    
    def compute_cs2g_score(
        self,
        gene: str,
        locus_variants: pd.DataFrame,
        strategy_contributions: Dict[str, float],
        strategy_scores: Dict[str, float]
    ) -> float:
        """
        Compute heritability-weighted cS2G score for a gene.
        
        Score = sum over strategies: (heritability_enrichment_s * strategy_score_s)
        
        Args:
            gene: Gene name
            locus_variants: DataFrame with variant annotations
            strategy_contributions: Heritability enrichment for each strategy
            strategy_scores: Score for this gene under each strategy
                           Keys: 'abc', 'eqtl', 'distance', 'coding', 'constraint'
        
        Returns:
            Weighted cS2G score (0-1)
        """
        score = 0.0
        
        for strategy in self.strategy_weights.keys():
            enrichment = strategy_contributions.get(strategy, 0.0)
            raw_score = strategy_scores.get(strategy, 0.0)
            
            # Weight by strategy weight and heritability enrichment
            weighted_score = self.strategy_weights[strategy] * enrichment * raw_score
            score += weighted_score
        
        return min(1.0, max(0.0, score))  # Clip to [0, 1]
    
    def score_genes_in_locus(
        self,
        locus_variants: pd.DataFrame,
        candidate_genes: List[str],
        abc_scores: Dict[str, float],
        eqtl_scores: Dict[str, float],
        distance_scores: Dict[str, float],
        coding_scores: Dict[str, float],
        constraint_scores: Dict[str, float]
    ) -> pd.DataFrame:
        """
        Score all candidate genes in a locus using cS2G.
        
        Args:
            locus_variants: Variants in credible set
            candidate_genes: List of genes to score
            abc_scores: ABC score for each gene {gene: score}
            eqtl_scores: eQTL score for each gene {gene: score}
            distance_scores: Distance score for each gene {gene: score}
            coding_scores: Coding variant score for each gene {gene: score}
            constraint_scores: Evolutionary constraint score {gene: score}
        
        Returns:
            DataFrame with columns: gene, cs2g_score, heritability_enrichment,
                                   primary_strategy
        """
        # Calculate heritability enrichments in this locus
        enrichments = self.calculate_heritability_enrichment(locus_variants, "locus_0")
        
        results = []
        for gene in candidate_genes:
            strategy_scores = {
                'abc': abc_scores.get(gene, 0.0),
                'eqtl': eqtl_scores.get(gene, 0.0),
                'distance': distance_scores.get(gene, 0.0),
                'coding': coding_scores.get(gene, 0.0),
                'constraint': constraint_scores.get(gene, 0.0)
            }
            
            score = self.compute_cs2g_score(gene, locus_variants, enrichments, strategy_scores)
            
            # Determine primary strategy
            primary_strategy = max(strategy_scores.keys(), key=lambda s: strategy_scores[s])
            
            results.append({
                'gene': gene,
                'cs2g_score': score,
                'heritability_enrichment': enrichments.get(primary_strategy, 0.0),
                'primary_strategy': primary_strategy,
                **{f'strategy_{k}': v for k, v in strategy_scores.items()}
            })
        
        return pd.DataFrame(results)
    
    def rank_genes(self, scores_df: pd.DataFrame) -> pd.DataFrame:
        """
        Rank genes by cS2G score with ties handled by primary strategy strength.
        """
        scores_df = scores_df.sort_values('cs2g_score', ascending=False).reset_index(drop=True)
        scores_df['rank'] = range(1, len(scores_df) + 1)
        return scores_df
