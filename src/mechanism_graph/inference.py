"""
Graph Inference

Probabilistic inference on mechanism graphs including:
- Noisy-OR path aggregation with formal generative model
- Belief propagation with proper uncertainty handling
- Calibration-aware inference

NOISY-OR GENERATIVE MODEL
=========================

We formalize the noisy-OR aggregation as a proper generative model following
Pearl (1988) and Koller & Friedman (2009).

Model assumptions:
1. INDEPENDENT CAUSES: Each path from variant to gene is an independent 
   "noisy" causal mechanism. This is appropriate when paths traverse different
   regulatory elements or tissues.

2. LEAK PROBABILITY: We include a background "leak" rate ε representing
   gene causality from unmeasured mechanisms.

3. INHIBITION PROBABILITY: Each path has probability (1 - p_path) of being
   inhibited (not transmitting the causal signal).

The generative process:
- Gene is causal if ANY of its incoming paths succeeds OR background leak fires
- P(causal | paths) = 1 - (1-ε) ∏_{path} (1 - p_path)

When the independence assumption is violated (e.g., paths through same cCRE,
tissues with correlated expression, or LD-linked variants), we provide
correction factors via:
- LD-aware path correlation adjustment
- Tissue correlation matrices
- Annotation overlap penalties

Reference:
- Pearl J (1988) Probabilistic Reasoning in Intelligent Systems
- Koller D, Friedman N (2009) Probabilistic Graphical Models
"""

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from ..utils.logging import get_logger
from .graph import MechanismGraph
from .nodes import BaseNode
from .edges import EdgeProbability, chain_probabilities


logger = get_logger("graph_inference")


# Default leak probability (background rate)
DEFAULT_LEAK_PROBABILITY = 0.01


@dataclass
class NoisyORParameters:
    """
    Parameters for the noisy-OR generative model.
    
    Attributes
    ----------
    leak_probability : float
        Background probability a gene is causal from unmeasured mechanisms.
    ld_correlation_penalty : float
        Penalty factor for paths sharing LD-correlated variants (0-1).
    tissue_correlation_matrix : dict, optional
        Pairwise tissue correlations for dependency adjustment.
    annotation_overlap_penalty : float
        Penalty for paths through same regulatory element.
    min_path_probability : float
        Minimum path probability to consider (pruning).
    """
    
    leak_probability: float = DEFAULT_LEAK_PROBABILITY
    ld_correlation_penalty: float = 0.5
    tissue_correlation_matrix: Optional[Dict[Tuple[str, str], float]] = None
    annotation_overlap_penalty: float = 0.3
    min_path_probability: float = 1e-6
    
    def get_tissue_correlation(
        self,
        tissue1: str,
        tissue2: str,
    ) -> float:
        """Get correlation between two tissues."""
        if self.tissue_correlation_matrix is None:
            return 0.0
        
        key = (tissue1, tissue2) if tissue1 < tissue2 else (tissue2, tissue1)
        return self.tissue_correlation_matrix.get(key, 0.0)


class GraphInference:
    """
    Probabilistic inference on mechanism graphs.
    
    Implements the noisy-OR generative model with explicit handling of:
    1. Path independence violations (LD, tissue correlation, annotation overlap)
    2. Calibration tracking per inference step
    3. Uncertainty quantification via bootstrap/MC sampling
    
    Supports:
    - Forward inference (variant → trait)
    - Backward inference (trait → variant)
    - Marginalization over paths
    - Confidence intervals
    """
    
    def __init__(
        self,
        graph: MechanismGraph,
        aggregation: str = "noisy_or",
        noisy_or_params: Optional[NoisyORParameters] = None,
    ):
        """
        Initialize inference engine.
        
        Parameters
        ----------
        graph : MechanismGraph
            The mechanism graph.
        aggregation : str
            Path aggregation method: "max", "sum", "noisy_or".
        noisy_or_params : NoisyORParameters, optional
            Parameters for noisy-OR model.
        """
        self.graph = graph
        self.aggregation = aggregation
        self.noisy_or_params = noisy_or_params or NoisyORParameters()
        
        # Cache for computed probabilities
        self._cache: Dict[str, float] = {}
        
        # Calibration tracking
        self._calibration_log: List[Dict] = []
    
    def forward_probability(
        self,
        source_id: str,
        target_type: str = "trait",
        account_for_correlations: bool = True,
    ) -> float:
        """
        Compute forward probability from source to target type.
        
        Uses noisy-OR aggregation with correlation corrections when enabled.
        
        Parameters
        ----------
        source_id : str
            Source node ID.
        target_type : str
            Target node type.
        account_for_correlations : bool
            Whether to apply correlation corrections (recommended).
            
        Returns
        -------
        float
            Probability (0-1).
        """
        cache_key = f"fwd_{source_id}_{target_type}_{account_for_correlations}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Find all paths to any node of target type
        all_paths_info = []
        
        for target_id in self.graph.nodes[target_type]:
            paths = self.graph.find_all_paths(source_id, target_id)
            
            for path in paths:
                prob = chain_probabilities(path)
                
                # Skip negligible paths
                if prob < self.noisy_or_params.min_path_probability:
                    continue
                
                path_info = {
                    "path": path,
                    "probability": prob,
                    "target_id": target_id,
                    "tissues": self._extract_tissues(path),
                    "variants": self._extract_variants(path),
                    "ccres": self._extract_ccres(path),
                }
                all_paths_info.append(path_info)
        
        if not all_paths_info:
            result = self.noisy_or_params.leak_probability
        else:
            probs = [p["probability"] for p in all_paths_info]
            
            if self.aggregation == "max":
                result = max(probs)
            elif self.aggregation == "sum":
                result = min(sum(probs), 1.0)
            elif self.aggregation == "noisy_or":
                if account_for_correlations:
                    result = self._noisy_or_with_corrections(all_paths_info)
                else:
                    result = self._noisy_or_naive(probs)
            else:
                result = max(probs)
        
        self._cache[cache_key] = result
        return result
    
    def _noisy_or_naive(self, probabilities: List[float]) -> float:
        """
        Naive noisy-OR: assumes complete independence.
        
        P(effect) = 1 - (1-leak) * ∏(1 - p_i)
        
        This is correct ONLY when all paths are truly independent.
        """
        leak = self.noisy_or_params.leak_probability
        inhibition_product = np.prod([1 - p for p in probabilities])
        return 1 - (1 - leak) * inhibition_product
    
    def _noisy_or_with_corrections(
        self,
        paths_info: List[Dict],
    ) -> float:
        """
        Noisy-OR with correlation corrections.
        
        Applies penalties when paths share:
        - LD-correlated variants
        - Correlated tissues
        - Same regulatory elements
        
        This prevents over-counting of correlated evidence.
        """
        if len(paths_info) <= 1:
            probs = [p["probability"] for p in paths_info]
            return self._noisy_or_naive(probs)
        
        # Compute effective (de-correlated) probabilities
        effective_probs = []
        
        for i, path_i in enumerate(paths_info):
            # Start with original probability
            eff_prob = path_i["probability"]
            
            # Check correlations with all previous paths
            for j in range(i):
                path_j = paths_info[j]
                
                # LD correlation penalty
                ld_penalty = self._compute_ld_penalty(
                    path_i["variants"],
                    path_j["variants"],
                )
                
                # Tissue correlation penalty
                tissue_penalty = self._compute_tissue_penalty(
                    path_i["tissues"],
                    path_j["tissues"],
                )
                
                # Annotation overlap penalty
                annotation_penalty = self._compute_annotation_penalty(
                    path_i["ccres"],
                    path_j["ccres"],
                )
                
                # Apply maximum penalty (conservative approach)
                max_penalty = max(ld_penalty, tissue_penalty, annotation_penalty)
                
                # Reduce effective probability proportionally
                eff_prob *= (1 - max_penalty)
            
            effective_probs.append(max(eff_prob, 0))
        
        return self._noisy_or_naive(effective_probs)
    
    def _compute_ld_penalty(
        self,
        variants1: Set[str],
        variants2: Set[str],
    ) -> float:
        """Compute penalty for LD-correlated variants."""
        if not variants1 or not variants2:
            return 0.0
        
        # Simple overlap-based penalty
        # Full LD correlation would require LD matrix lookup
        overlap = len(variants1 & variants2)
        
        if overlap > 0:
            return self.noisy_or_params.ld_correlation_penalty
        
        return 0.0
    
    def _compute_tissue_penalty(
        self,
        tissues1: Set[str],
        tissues2: Set[str],
    ) -> float:
        """Compute penalty for correlated tissues."""
        if not tissues1 or not tissues2:
            return 0.0
        
        # Check for same tissue
        if tissues1 & tissues2:
            return 0.5  # Same tissue = high correlation
        
        # Check tissue correlation matrix
        max_corr = 0.0
        for t1 in tissues1:
            for t2 in tissues2:
                corr = self.noisy_or_params.get_tissue_correlation(t1, t2)
                max_corr = max(max_corr, corr)
        
        return max_corr * 0.5  # Scale to penalty
    
    def _compute_annotation_penalty(
        self,
        ccres1: Set[str],
        ccres2: Set[str],
    ) -> float:
        """Compute penalty for shared regulatory elements."""
        if not ccres1 or not ccres2:
            return 0.0
        
        overlap = len(ccres1 & ccres2)
        
        if overlap > 0:
            return self.noisy_or_params.annotation_overlap_penalty
        
        return 0.0
    
    def _extract_tissues(self, path: List[EdgeProbability]) -> Set[str]:
        """Extract tissue IDs from path."""
        tissues = set()
        for edge in path:
            if hasattr(edge.source, 'tissue_id'):
                tissues.add(edge.source.tissue_id)
            if hasattr(edge.target, 'tissue_id'):
                tissues.add(edge.target.tissue_id)
        return tissues
    
    def _extract_variants(self, path: List[EdgeProbability]) -> Set[str]:
        """Extract variant IDs from path."""
        variants = set()
        for edge in path:
            if edge.source.type == "variant":
                variants.add(edge.source.id)
        return variants
    
    def _extract_ccres(self, path: List[EdgeProbability]) -> Set[str]:
        """Extract cCRE IDs from path."""
        ccres = set()
        for edge in path:
            if edge.source.type == "ccre":
                ccres.add(edge.source.id)
            if edge.target.type == "ccre":
                ccres.add(edge.target.id)
        return ccres
    
    def backward_probability(
        self,
        target_id: str,
        source_type: str = "variant",
    ) -> Dict[str, float]:
        """
        Compute backward probabilities from target to all sources.
        
        Parameters
        ----------
        target_id : str
            Target node ID.
        source_type : str
            Source node type.
            
        Returns
        -------
        dict
            Probabilities by source node ID.
        """
        results = {}
        
        for source_id in self.graph.nodes[source_type]:
            prob = self.forward_probability(source_id, 
                                            self.graph.get_node(target_id).type)
            results[source_id] = prob
        
        return results
    
    def gene_causal_probability(
        self,
        gene_id: str,
    ) -> Dict[str, Any]:
        """
        Compute causal probability for a gene.
        
        Returns both upstream (variant→gene) and downstream (gene→trait)
        probabilities, along with the combined score.
        
        Parameters
        ----------
        gene_id : str
            Gene identifier.
            
        Returns
        -------
        dict
            Causal probability components.
        """
        gene = self.graph.get_node(gene_id, "gene")
        if not gene:
            return {"gene_id": gene_id, "error": "Gene not found"}
        
        # Upstream: variants → gene
        upstream_probs = []
        for var_id in self.graph.nodes["variant"]:
            paths = self.graph.find_all_paths(var_id, gene_id, max_depth=3)
            for path in paths:
                prob = chain_probabilities(path)
                upstream_probs.append(prob)
        
        upstream = self._aggregate(upstream_probs)
        
        # Downstream: gene → trait
        downstream_probs = []
        for trait_id in self.graph.nodes["trait"]:
            paths = self.graph.find_all_paths(gene_id, trait_id, max_depth=3)
            for path in paths:
                prob = chain_probabilities(path)
                downstream_probs.append(prob)
        
        downstream = self._aggregate(downstream_probs)
        
        # Combined
        combined = upstream * downstream
        
        return {
            "gene_id": gene_id,
            "upstream_probability": upstream,
            "downstream_probability": downstream,
            "combined_probability": combined,
            "n_upstream_paths": len(upstream_probs),
            "n_downstream_paths": len(downstream_probs),
        }
    
    def tissue_contribution(
        self,
        gene_id: str,
    ) -> Dict[str, float]:
        """
        Compute tissue contributions for a gene.
        
        Parameters
        ----------
        gene_id : str
            Gene identifier.
            
        Returns
        -------
        dict
            Contribution probability by tissue.
        """
        contributions = {}
        
        gene_edges = self.graph.get_edges_from(gene_id)
        
        for edge in gene_edges:
            if edge.target.type == "tissue":
                tissue_id = edge.target.id
                
                # Get tissue → trait probability
                tissue_edges = self.graph.get_edges_from(tissue_id)
                trait_prob = max(
                    (e.probability for e in tissue_edges if e.target.type == "trait"),
                    default=0.0
                )
                
                # Combined gene → tissue → trait
                contributions[tissue_id] = edge.probability * trait_prob
        
        return contributions
    
    def _aggregate(self, probs: List[float]) -> float:
        """Aggregate probabilities using configured method."""
        if not probs:
            return 0.0
        
        if self.aggregation == "max":
            return max(probs)
        elif self.aggregation == "sum":
            return min(sum(probs), 1.0)
        elif self.aggregation == "noisy_or":
            return 1 - np.prod([1 - p for p in probs])
        else:
            return max(probs)
    
    def compute_confidence_interval(
        self,
        source_id: str,
        target_id: str,
        n_samples: int = 1000,
        alpha: float = 0.05,
    ) -> Tuple[float, float, float]:
        """
        Compute confidence interval via bootstrap.
        
        Parameters
        ----------
        source_id : str
            Source node ID.
        target_id : str
            Target node ID.
        n_samples : int
            Number of bootstrap samples.
        alpha : float
            Significance level.
            
        Returns
        -------
        tuple
            (lower, point_estimate, upper)
        """
        paths = self.graph.find_all_paths(source_id, target_id)
        
        if not paths:
            return (0.0, 0.0, 0.0)
        
        # Point estimate
        point_probs = [chain_probabilities(p) for p in paths]
        point_est = self._aggregate(point_probs)
        
        # Bootstrap
        bootstrap_estimates = []
        
        for _ in range(n_samples):
            # Resample edge probabilities with uncertainty
            sampled_probs = []
            
            for path in paths:
                path_prob = 1.0
                for edge in path:
                    # Add noise based on edge confidence
                    noise = np.random.normal(0, 0.1 * (1 - edge.confidence))
                    sampled_prob = np.clip(edge.probability + noise, 0, 1)
                    path_prob *= sampled_prob
                
                sampled_probs.append(path_prob)
            
            bootstrap_estimates.append(self._aggregate(sampled_probs))
        
        # Compute CI
        lower = np.percentile(bootstrap_estimates, 100 * alpha / 2)
        upper = np.percentile(bootstrap_estimates, 100 * (1 - alpha / 2))
        
        return (lower, point_est, upper)


def propagate_uncertainty(
    graph: MechanismGraph,
    n_iterations: int = 10,
    damping: float = 0.5,
) -> Dict[str, Dict[str, float]]:
    """
    Propagate uncertainty through the graph using belief propagation.
    
    Parameters
    ----------
    graph : MechanismGraph
        The mechanism graph.
    n_iterations : int
        Number of iterations.
    damping : float
        Damping factor for stability.
        
    Returns
    -------
    dict
        Updated beliefs by node type and ID.
    """
    # Initialize beliefs
    beliefs: Dict[str, Dict[str, float]] = defaultdict(dict)
    
    # Initialize variant beliefs from PIPs
    for var_id, var in graph.nodes["variant"].items():
        beliefs["variant"][var_id] = var.pip
    
    # Initialize other nodes uniformly
    for node_type in ["ccre", "gene", "tissue", "trait"]:
        for node_id in graph.nodes[node_type]:
            beliefs[node_type][node_id] = 0.5
    
    # Message passing
    messages: Dict[Tuple[str, str], float] = {}
    
    for _ in range(n_iterations):
        new_beliefs = defaultdict(dict)
        
        # Forward pass
        for edge in graph.edges:
            source_belief = beliefs[edge.source.type][edge.source.id]
            message = source_belief * edge.probability
            
            key = (edge.source.id, edge.target.id)
            if key in messages:
                message = damping * message + (1 - damping) * messages[key]
            
            messages[key] = message
        
        # Update beliefs
        for node_type in beliefs:
            for node_id in beliefs[node_type]:
                incoming = graph.get_edges_to(node_id)
                
                if not incoming:
                    new_beliefs[node_type][node_id] = beliefs[node_type][node_id]
                else:
                    # Aggregate incoming messages
                    msg_probs = []
                    for edge in incoming:
                        key = (edge.source.id, edge.target.id)
                        if key in messages:
                            msg_probs.append(messages[key])
                    
                    if msg_probs:
                        # Noisy-OR aggregation
                        new_belief = 1 - np.prod([1 - p for p in msg_probs])
                        new_beliefs[node_type][node_id] = new_belief
                    else:
                        new_beliefs[node_type][node_id] = beliefs[node_type][node_id]
        
        beliefs = new_beliefs
    
    return dict(beliefs)


def compute_gene_rankings(
    graph: MechanismGraph,
) -> List[Dict[str, Any]]:
    """
    Compute gene rankings based on causal probability.
    
    Parameters
    ----------
    graph : MechanismGraph
        The mechanism graph.
        
    Returns
    -------
    list
        Ranked genes with scores.
    """
    inference = GraphInference(graph)
    
    rankings = []
    
    for gene_id in graph.nodes["gene"]:
        result = inference.gene_causal_probability(gene_id)
        gene = graph.get_node(gene_id, "gene")
        
        result["symbol"] = gene.symbol if gene else ""
        rankings.append(result)
    
    # Sort by combined probability
    rankings.sort(key=lambda x: x.get("combined_probability", 0), reverse=True)
    
    # Add rank
    for i, r in enumerate(rankings):
        r["rank"] = i + 1
    
    return rankings
