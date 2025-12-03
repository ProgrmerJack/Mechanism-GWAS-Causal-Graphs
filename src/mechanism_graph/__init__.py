"""
Mechanism Graph Module

Constructs probabilistic mechanism graphs linking:
Variants → cCREs → Genes → Tissues → Traits

This is the core innovation - a structured graphical representation
with calibrated edge probabilities that captures the full causal chain.
"""

from .graph import MechanismGraph, MechanismGraphBuilder
from .nodes import VariantNode, CCRENode, GeneNode, TissueNode, TraitNode
from .edges import EdgeProbability, compute_edge_weights
from .inference import GraphInference, propagate_uncertainty

__all__ = [
    "MechanismGraph",
    "MechanismGraphBuilder",
    "VariantNode",
    "CCRENode",
    "GeneNode",
    "TissueNode",
    "TraitNode",
    "EdgeProbability",
    "compute_edge_weights",
    "GraphInference",
    "propagate_uncertainty",
]
