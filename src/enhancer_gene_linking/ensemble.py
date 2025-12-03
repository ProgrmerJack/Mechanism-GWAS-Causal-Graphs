"""
Ensemble Enhancer-Gene Linking

Combines multiple evidence sources for enhancer-gene links:
- ABC model (Activity-by-Contact)
- PCHi-C (Promoter Capture Hi-C)
- Distance-based (fallback)

This module implements the "bridge ablation" analysis requested by reviewers:
- Distance-only baseline
- ABC-only
- PCHi-C-only  
- Ensemble (recommended)

The ensemble uses a principled combination that:
1. Prioritizes ABC when available (best precision)
2. Augments with PCHi-C for coverage
3. Falls back to distance only when no experimental evidence exists
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ..utils.logging import get_logger
from .abc_links import ABCLinks, ABCPrediction
from .pchic_links import PCHiCLinks, PCHiCLink


logger = get_logger("ensemble_links")


@dataclass
class EnsembleLinkEvidence:
    """
    Combined evidence for an enhancer-gene link.
    
    Attributes
    ----------
    enhancer_id : str
        Enhancer identifier (chr:start-end).
    gene : str
        Target gene symbol.
    ensemble_probability : float
        Combined probability (0-1).
    abc_probability : float
        ABC-derived probability.
    pchic_probability : float
        PCHi-C-derived probability.
    distance_probability : float
        Distance-based probability.
    evidence_sources : list
        Sources contributing to this link.
    cell_types : dict
        Cell types by source.
    abc_score : float
        Best ABC score.
    chicago_score : float
        Best CHiCAGO score.
    distance : int
        Distance to TSS.
    is_validated : bool
        CRISPRi validation status.
    """
    
    enhancer_id: str
    gene: str
    ensemble_probability: float = 0.0
    abc_probability: float = 0.0
    pchic_probability: float = 0.0
    distance_probability: float = 0.0
    evidence_sources: List[str] = field(default_factory=list)
    cell_types: Dict[str, List[str]] = field(default_factory=dict)
    abc_score: float = 0.0
    chicago_score: float = 0.0
    distance: int = 0
    is_validated: bool = False
    
    @property
    def has_experimental_evidence(self) -> bool:
        """Whether link has ABC or PCHi-C support."""
        return self.abc_probability > 0 or self.pchic_probability > 0
    
    @property
    def n_evidence_sources(self) -> int:
        """Number of evidence sources."""
        return len(self.evidence_sources)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "enhancer_id": self.enhancer_id,
            "gene": self.gene,
            "ensemble_probability": self.ensemble_probability,
            "abc_probability": self.abc_probability,
            "pchic_probability": self.pchic_probability,
            "distance_probability": self.distance_probability,
            "evidence_sources": self.evidence_sources,
            "n_sources": self.n_evidence_sources,
            "has_experimental": self.has_experimental_evidence,
            "cell_types": self.cell_types,
            "abc_score": self.abc_score,
            "chicago_score": self.chicago_score,
            "distance": self.distance,
            "is_validated": self.is_validated,
        }


class EnhancerGeneEnsemble:
    """
    Ensemble enhancer-gene linker combining multiple evidence sources.
    
    Implements the "bridge ablation" analysis:
    - method="distance": Distance-only baseline
    - method="abc": ABC-only
    - method="pchic": PCHi-C-only
    - method="ensemble": Combined (default, recommended)
    
    Example
    -------
    >>> abc = ABCLinks()
    >>> abc.load_nasser2021("abc_predictions.tsv.gz")
    >>> pchic = PCHiCLinks()
    >>> pchic.load_blueprint("pchic_calls.txt.gz")
    >>> 
    >>> ensemble = EnhancerGeneEnsemble(abc, pchic)
    >>> prob, evidence = ensemble.compute_edge_probability(
    ...     "chr1", 1000000, 1010000, "GENE1"
    ... )
    """
    
    # Ensemble weights (can be calibrated)
    DEFAULT_WEIGHTS = {
        "abc": 0.5,      # ABC has best precision
        "pchic": 0.3,    # PCHi-C provides coverage
        "distance": 0.2,  # Distance as fallback
    }
    
    # When experimental evidence exists, downweight distance
    EXPERIMENTAL_WEIGHTS = {
        "abc": 0.6,
        "pchic": 0.35,
        "distance": 0.05,  # Minimal contribution when we have experimental data
    }
    
    def __init__(
        self,
        abc_links: Optional[ABCLinks] = None,
        pchic_links: Optional[PCHiCLinks] = None,
        weights: Optional[Dict[str, float]] = None,
        distance_decay_kb: float = 100.0,
    ):
        """
        Initialize ensemble linker.
        
        Parameters
        ----------
        abc_links : ABCLinks, optional
            ABC predictions manager.
        pchic_links : PCHiCLinks, optional
            PCHi-C links manager.
        weights : dict, optional
            Custom weights for each source.
        distance_decay_kb : float
            Half-life for distance decay in kb.
        """
        self.abc = abc_links
        self.pchic = pchic_links
        self.weights = weights or self.DEFAULT_WEIGHTS
        self.distance_decay_kb = distance_decay_kb
    
    def compute_edge_probability(
        self,
        enhancer_chr: str,
        enhancer_start: int,
        enhancer_end: int,
        gene: str,
        gene_tss: Optional[int] = None,
        cell_types: Optional[List[str]] = None,
        method: str = "ensemble",
    ) -> Tuple[float, EnsembleLinkEvidence]:
        """
        Compute edge probability for enhancer→gene link.
        
        Parameters
        ----------
        enhancer_chr : str
            Enhancer chromosome.
        enhancer_start : int
            Enhancer start position.
        enhancer_end : int
            Enhancer end position.
        gene : str
            Target gene symbol.
        gene_tss : int, optional
            Gene TSS position (for distance calculation).
        cell_types : list, optional
            Filter to specific cell types.
        method : str
            Linking method: "ensemble", "abc", "pchic", "distance".
            
        Returns
        -------
        tuple
            (probability, EnsembleLinkEvidence)
        """
        enhancer_id = f"{enhancer_chr}:{enhancer_start}-{enhancer_end}"
        enhancer_mid = (enhancer_start + enhancer_end) // 2
        
        # Initialize evidence
        evidence = EnsembleLinkEvidence(
            enhancer_id=enhancer_id,
            gene=gene,
        )
        
        # Compute distance probability
        if gene_tss is not None:
            distance = abs(enhancer_mid - gene_tss)
            evidence.distance = distance
            evidence.distance_probability = self._distance_probability(distance)
        
        # ABC probability
        if self.abc is not None:
            abc_prob, abc_evidence = self.abc.compute_edge_probability(
                enhancer_chr, enhancer_start, enhancer_end, gene, cell_types
            )
            evidence.abc_probability = abc_prob
            if abc_prob > 0:
                evidence.evidence_sources.append("abc")
                evidence.abc_score = abc_evidence.get("abc_score", 0)
                evidence.cell_types["abc"] = abc_evidence.get("cell_types", [])
                evidence.is_validated = abc_evidence.get("validated", False)
        
        # PCHi-C probability
        if self.pchic is not None:
            pchic_prob, pchic_evidence = self.pchic.compute_edge_probability(
                enhancer_chr, enhancer_start, enhancer_end, gene, cell_types
            )
            evidence.pchic_probability = pchic_prob
            if pchic_prob > 0:
                evidence.evidence_sources.append("pchic")
                evidence.chicago_score = pchic_evidence.get("chicago_score", 0)
                evidence.cell_types["pchic"] = pchic_evidence.get("cell_types", [])
        
        # Add distance as source if no experimental evidence
        if evidence.distance_probability > 0:
            if not evidence.has_experimental_evidence:
                evidence.evidence_sources.append("distance")
        
        # Compute final probability based on method
        if method == "distance":
            evidence.ensemble_probability = evidence.distance_probability
        elif method == "abc":
            evidence.ensemble_probability = evidence.abc_probability
        elif method == "pchic":
            evidence.ensemble_probability = evidence.pchic_probability
        elif method == "ensemble":
            evidence.ensemble_probability = self._ensemble_probability(evidence)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return evidence.ensemble_probability, evidence
    
    def _distance_probability(self, distance: int) -> float:
        """
        Compute distance-based probability.
        
        Uses exponential decay with configurable half-life.
        
        Parameters
        ----------
        distance : int
            Distance in bp.
            
        Returns
        -------
        float
            Probability (0-1).
        """
        # Promoter proximity (< 2kb) gets high probability
        if distance < 2000:
            return 0.9
        
        # Exponential decay
        decay_bp = self.distance_decay_kb * 1000
        prob = np.exp(-distance / decay_bp) * 0.5
        
        return prob
    
    def _ensemble_probability(
        self,
        evidence: EnsembleLinkEvidence,
    ) -> float:
        """
        Compute ensemble probability from multiple sources.
        
        Uses weighted combination with adaptive weights based on
        whether experimental evidence is available.
        
        Parameters
        ----------
        evidence : EnsembleLinkEvidence
            Evidence from all sources.
            
        Returns
        -------
        float
            Ensemble probability (0-1).
        """
        # Use different weights based on experimental evidence availability
        if evidence.has_experimental_evidence:
            weights = self.EXPERIMENTAL_WEIGHTS
        else:
            weights = self.weights
        
        # Weighted combination
        prob = (
            weights.get("abc", 0) * evidence.abc_probability +
            weights.get("pchic", 0) * evidence.pchic_probability +
            weights.get("distance", 0) * evidence.distance_probability
        )
        
        # Normalize by sum of weights for available sources
        weight_sum = sum(
            weights.get(s, 0) for s in evidence.evidence_sources
        )
        
        if weight_sum > 0:
            prob = prob / weight_sum
        
        # Boost for multiple concordant sources
        if evidence.n_evidence_sources >= 2:
            # Concordance boost: if both ABC and PCHi-C support, increase confidence
            if evidence.abc_probability > 0.5 and evidence.pchic_probability > 0.5:
                prob = min(0.95, prob * 1.2)
        
        # Boost for validated links
        if evidence.is_validated:
            prob = min(0.98, prob + 0.15)
        
        return prob
    
    def get_links_for_gene(
        self,
        gene: str,
        cell_types: Optional[List[str]] = None,
        method: str = "ensemble",
    ) -> List[EnsembleLinkEvidence]:
        """
        Get all enhancer links for a gene.
        
        Parameters
        ----------
        gene : str
            Gene symbol.
        cell_types : list, optional
            Filter cell types.
        method : str
            Linking method.
            
        Returns
        -------
        list
            List of EnsembleLinkEvidence objects.
        """
        results = []
        seen_enhancers = set()
        
        # Get ABC links
        if self.abc is not None:
            for pred in self.abc.get_links_for_gene(gene, cell_types):
                enhancer_id = pred.enhancer_id
                if enhancer_id not in seen_enhancers:
                    seen_enhancers.add(enhancer_id)
                    _, evidence = self.compute_edge_probability(
                        pred.enhancer_chr, pred.enhancer_start, pred.enhancer_end,
                        gene, cell_types=cell_types, method=method
                    )
                    results.append(evidence)
        
        # Get PCHi-C links
        if self.pchic is not None:
            for link in self.pchic.get_links_for_gene(gene, cell_types):
                enhancer_id = link.oe_id
                if enhancer_id not in seen_enhancers:
                    seen_enhancers.add(enhancer_id)
                    _, evidence = self.compute_edge_probability(
                        link.oe_chr, link.oe_start, link.oe_end,
                        gene, cell_types=cell_types, method=method
                    )
                    results.append(evidence)
        
        return sorted(results, key=lambda x: x.ensemble_probability, reverse=True)
    
    def bridge_ablation_analysis(
        self,
        enhancer_chr: str,
        enhancer_start: int,
        enhancer_end: int,
        gene: str,
        gene_tss: Optional[int] = None,
        cell_types: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Run bridge ablation analysis comparing all methods.
        
        This directly addresses reviewer request to compare:
        - Distance-only vs ABC vs PCHi-C vs Ensemble
        
        Parameters
        ----------
        enhancer_chr, enhancer_start, enhancer_end : str, int, int
            Enhancer coordinates.
        gene : str
            Target gene.
        gene_tss : int, optional
            Gene TSS.
        cell_types : list, optional
            Filter cell types.
            
        Returns
        -------
        dict
            Comparison of all methods.
        """
        results = {}
        
        for method in ["distance", "abc", "pchic", "ensemble"]:
            prob, evidence = self.compute_edge_probability(
                enhancer_chr, enhancer_start, enhancer_end,
                gene, gene_tss, cell_types, method=method
            )
            results[method] = {
                "probability": prob,
                "evidence": evidence.to_dict(),
            }
        
        # Add comparison metrics
        results["comparison"] = {
            "distance_vs_abc": results["abc"]["probability"] - results["distance"]["probability"],
            "distance_vs_pchic": results["pchic"]["probability"] - results["distance"]["probability"],
            "distance_vs_ensemble": results["ensemble"]["probability"] - results["distance"]["probability"],
            "abc_vs_pchic": results["abc"]["probability"] - results["pchic"]["probability"],
            "best_method": max(results.keys() - {"comparison"}, 
                             key=lambda m: results[m]["probability"]),
        }
        
        return results
    
    def summary(self) -> Dict[str, Any]:
        """Get summary of available evidence sources."""
        summary = {
            "abc_available": self.abc is not None,
            "pchic_available": self.pchic is not None,
            "weights": self.weights,
            "distance_decay_kb": self.distance_decay_kb,
        }
        
        if self.abc is not None:
            summary["abc_summary"] = self.abc.summary()
        
        if self.pchic is not None:
            summary["pchic_summary"] = self.pchic.summary()
        
        return summary


def compute_ensemble_score(
    abc_score: float,
    chicago_score: float,
    distance: int,
    distance_decay_kb: float = 100.0,
) -> float:
    """
    Quick ensemble score computation.
    
    Convenience function for computing ensemble probability from raw scores.
    
    Parameters
    ----------
    abc_score : float
        ABC score (0-1).
    chicago_score : float
        CHiCAGO score.
    distance : int
        Distance to TSS in bp.
    distance_decay_kb : float
        Distance decay half-life in kb.
        
    Returns
    -------
    float
        Ensemble probability.
    """
    # Convert ABC score to probability
    abc_prob = 1 / (1 + np.exp(-50 * (abc_score - 0.015))) if abc_score > 0 else 0
    
    # Convert CHiCAGO score to probability
    pchic_prob = 1 / (1 + np.exp(-0.5 * (chicago_score - 5))) if chicago_score > 0 else 0
    
    # Distance probability
    if distance < 2000:
        dist_prob = 0.9
    else:
        dist_prob = np.exp(-distance / (distance_decay_kb * 1000)) * 0.5
    
    # Ensemble
    has_experimental = abc_prob > 0 or pchic_prob > 0
    
    if has_experimental:
        weights = (0.6, 0.35, 0.05)
    else:
        weights = (0.5, 0.3, 0.2)
    
    prob = weights[0] * abc_prob + weights[1] * pchic_prob + weights[2] * dist_prob
    
    # Normalize
    active_weights = sum(
        w for w, p in zip(weights, [abc_prob, pchic_prob, dist_prob]) if p > 0
    )
    
    if active_weights > 0:
        prob = prob / active_weights * sum(weights)
    
    return min(prob, 1.0)
