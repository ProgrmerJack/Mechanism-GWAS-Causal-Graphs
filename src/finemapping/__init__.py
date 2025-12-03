"""
Fine-Mapping Module

Implements statistical fine-mapping using SuSiE and related methods
with support for tissue-specific regulatory priors.
"""

from .locus import LocusDefinition, define_loci
from .susie import SuSiEFinemapper, run_finemapping
from .priors import RegulatoryPriors, TissuePriors

__all__ = [
    "LocusDefinition",
    "define_loci",
    "SuSiEFinemapper",
    "run_finemapping",
    "RegulatoryPriors",
    "TissuePriors",
]
