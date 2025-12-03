"""
Calibration Module

Implements calibration and benchmarking framework for mechanism graphs.
This is essential for a Nature paper - proving that probabilities are meaningful.

Key components:
- ProbabilityCalibration: ECE, MCE, reliability curves
- ModuleCalibrationSuite: Per-module calibration analysis
- ThreeTierBenchmark: Anti-leak benchmark with provenance
- PlattScaling: Recalibration via Platt scaling
"""

from .calibration import CalibrationAnalysis, reliability_curve
from .benchmarks import (
    BenchmarkSuite,
    GoldStandardGenes,
    BenchmarkTier,
    ThreeTierBenchmark,
)
from .metrics import compute_metrics, CalibrationMetrics
from .probability_calibration import (
    ProbabilityCalibration,
    CalibrationResult,
    CalibrationBin,
    ModuleCalibrationSuite,
    PlattScaling,
    compute_ece,
    compute_mce,
)

__all__ = [
    # Core calibration
    "CalibrationAnalysis",
    "reliability_curve",
    # Probability calibration
    "ProbabilityCalibration",
    "CalibrationResult",
    "CalibrationBin",
    "ModuleCalibrationSuite",
    "PlattScaling",
    "compute_ece",
    "compute_mce",
    # Benchmarks
    "BenchmarkSuite",
    "GoldStandardGenes",
    "BenchmarkTier",
    "ThreeTierBenchmark",
    # Metrics
    "compute_metrics",
    "CalibrationMetrics",
]
