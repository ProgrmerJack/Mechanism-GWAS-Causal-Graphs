"""
Colocalization Module

Implements Bayesian colocalization testing (COLOC) to assess
whether GWAS and QTL signals share a common causal variant.

Key components:
- ColocAnalysis: Standard COLOC (Giambartolomei 2014)
- run_susie_coloc: coloc.susie for multi-signal handling (Wallace 2021)
- EQTLCatalogueLoader: Cross-study replication via eQTL Catalogue
- CrossStudyReplication: Discovery in GTEx → Replication test
"""

from .coloc import (
    ColocAnalysis,
    run_coloc,
    run_susie_coloc,
    SuSiEColocResult,
)
from .eqtl import EQTLLoader, query_eqtl_catalogue
from .colocquake import ColocQuake
from .eqtl_catalogue import (
    EQTLCatalogueLoader,
    CrossStudyReplication,
    EQTLSignal,
    ReplicationResult,
    match_tissues,
)

__all__ = [
    # Core colocalization
    "ColocAnalysis",
    "run_coloc",
    "run_susie_coloc",
    "SuSiEColocResult",
    # eQTL loading
    "EQTLLoader",
    "query_eqtl_catalogue",
    # Multi-signal
    "ColocQuake",
    # Cross-study replication
    "EQTLCatalogueLoader",
    "CrossStudyReplication",
    "EQTLSignal",
    "ReplicationResult",
    "match_tissues",
]
