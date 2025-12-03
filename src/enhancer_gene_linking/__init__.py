"""
Enhancer-Gene Linking Module

Provides experimentally-grounded enhancer-gene links beyond distance-only heuristics.

Key components:
- abc_links: Activity-by-Contact model links (Fulco et al. 2019, Nasser et al. 2021)
- pchic_links: Promoter Capture Hi-C links (Jung et al. 2019)
- ensemble: Combined evidence from multiple sources

This module addresses the "weakest link in the causal chain" critique:
distance-only cCRE→gene edges have high false positive rates.

References:
- Fulco CP et al. (2019) Activity-by-contact model. Nat Genet.
- Nasser J et al. (2021) Genome-wide enhancer maps. Nature.
- Jung I et al. (2019) PCHi-C enhancer maps. Nat Genet.
"""

from .abc_links import ABCLinks, load_abc_predictions
from .pchic_links import PCHiCLinks, load_pchic_links
from .ensemble import EnhancerGeneEnsemble, compute_ensemble_score

__all__ = [
    "ABCLinks",
    "load_abc_predictions",
    "PCHiCLinks", 
    "load_pchic_links",
    "EnhancerGeneEnsemble",
    "compute_ensemble_score",
]
