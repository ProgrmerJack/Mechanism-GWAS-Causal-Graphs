"""
Baseline Comparison Module
===========================

This module implements contemporary baseline methods for SNP-to-gene linking
to enable fair comparison with our path-probability mechanism graph approach.

Baselines implemented:
1. FLAMES: XGBoost + PoPS integration (Schipper et al. 2025)
2. cS2G: Heritability-based S2G combination (Gazal et al. 2022)
3. PoPS: Polygenic priority scores (Weeks et al. 2023)
4. Effector Index: Evidence-weighted scoring (Smemo et al. 2022)

Author: Mechanism-GWAS-Causal-Graphs team
Date: December 2025
"""

from .flames_approximation import FLAMESApproximation
from .cs2g_implementation import CS2GImplementation
from .pops_runner import POPSRunner
from .effector_index import EffectorIndex
from .baseline_runner import UnifiedBaselineRunner, ValidationMetrics

__all__ = [
    'FLAMESApproximation',
    'CS2GImplementation',
    'POPSRunner',
    'EffectorIndex',
    'UnifiedBaselineRunner',
    'ValidationMetrics',
]

