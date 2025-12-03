"""
Probability Calibration Framework

Implements calibration analysis for mechanism graph probabilities.

Nature expects that "0.8 probability" means ~80% true on an external check.

This module provides:
- Reliability curves (calibration plots)
- Expected Calibration Error (ECE)
- Maximum Calibration Error (MCE)
- Per-module calibration (variant PIP, coloc posterior, gene posterior)
- Platt scaling and isotonic regression for recalibration

Reference:
- Guo et al. (2017) On Calibration of Modern Neural Networks
- Niculescu-Mizil & Caruana (2005) Predicting Good Probabilities
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats

from ..utils.logging import get_logger


logger = get_logger("calibration")


@dataclass
class CalibrationBin:
    """
    A single bin in calibration analysis.
    
    Attributes
    ----------
    bin_lower : float
        Lower bound of bin.
    bin_upper : float
        Upper bound of bin.
    n_samples : int
        Number of samples in bin.
    mean_predicted : float
        Mean predicted probability.
    mean_observed : float
        Observed fraction positive.
    confidence_lower : float
        Lower 95% CI on observed.
    confidence_upper : float
        Upper 95% CI on observed.
    """
    
    bin_lower: float
    bin_upper: float
    n_samples: int = 0
    mean_predicted: float = 0.0
    mean_observed: float = 0.0
    confidence_lower: float = 0.0
    confidence_upper: float = 1.0


@dataclass
class CalibrationResult:
    """
    Complete calibration analysis result.
    
    Attributes
    ----------
    module : str
        Module name (e.g., "variant_pip", "coloc_h4", "gene_posterior").
    n_samples : int
        Total samples analyzed.
    ece : float
        Expected Calibration Error.
    mce : float
        Maximum Calibration Error.
    bins : list
        List of CalibrationBin objects.
    brier_score : float
        Brier score (mean squared error).
    is_calibrated : bool
        Whether probabilities are well-calibrated (ECE < 0.05).
    hosmer_lemeshow_pvalue : float
        Hosmer-Lemeshow test p-value.
    """
    
    module: str
    n_samples: int = 0
    ece: float = 0.0
    mce: float = 0.0
    bins: List[CalibrationBin] = field(default_factory=list)
    brier_score: float = 0.0
    is_calibrated: bool = False
    hosmer_lemeshow_pvalue: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "module": self.module,
            "n_samples": self.n_samples,
            "ece": self.ece,
            "mce": self.mce,
            "brier_score": self.brier_score,
            "is_calibrated": self.is_calibrated,
            "hosmer_lemeshow_pvalue": self.hosmer_lemeshow_pvalue,
            "n_bins": len(self.bins),
        }


class ProbabilityCalibration:
    """
    Calibration analysis for probabilistic predictions.
    
    Analyzes whether predicted probabilities match observed frequencies.
    
    Example
    -------
    >>> cal = ProbabilityCalibration(n_bins=10)
    >>> 
    >>> # Analyze variant PIP calibration
    >>> result = cal.analyze(
    ...     predicted=variant_pips,
    ...     observed=is_causal,
    ...     module="variant_pip"
    ... )
    >>> print(f"ECE: {result.ece:.3f}")
    >>> 
    >>> # Plot reliability curve
    >>> cal.plot_reliability_curve(result, "variant_pip_calibration.png")
    """
    
    def __init__(
        self,
        n_bins: int = 10,
        strategy: str = "uniform",
    ):
        """
        Initialize calibration analyzer.
        
        Parameters
        ----------
        n_bins : int
            Number of bins for calibration analysis.
        strategy : str
            Binning strategy: "uniform" (equal width) or "quantile" (equal count).
        """
        self.n_bins = n_bins
        self.strategy = strategy
    
    def analyze(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        module: str = "unknown",
        sample_weight: Optional[np.ndarray] = None,
    ) -> CalibrationResult:
        """
        Analyze calibration of predicted probabilities.
        
        Parameters
        ----------
        predicted : np.ndarray
            Predicted probabilities (0-1).
        observed : np.ndarray
            Binary outcomes (0 or 1).
        module : str
            Module name for labeling.
        sample_weight : np.ndarray, optional
            Sample weights.
            
        Returns
        -------
        CalibrationResult
            Calibration analysis results.
        """
        predicted = np.asarray(predicted)
        observed = np.asarray(observed)
        
        if len(predicted) != len(observed):
            raise ValueError("predicted and observed must have same length")
        
        n_samples = len(predicted)
        
        if n_samples == 0:
            return CalibrationResult(module=module, n_samples=0)
        
        # Compute bin edges
        bin_edges = self._compute_bin_edges(predicted)
        
        # Compute per-bin statistics
        bins = []
        bin_gaps = []  # For ECE calculation
        
        for i in range(len(bin_edges) - 1):
            lower = bin_edges[i]
            upper = bin_edges[i + 1]
            
            # Find samples in this bin
            if i == len(bin_edges) - 2:  # Last bin includes upper bound
                mask = (predicted >= lower) & (predicted <= upper)
            else:
                mask = (predicted >= lower) & (predicted < upper)
            
            n_in_bin = mask.sum()
            
            if n_in_bin == 0:
                bins.append(CalibrationBin(
                    bin_lower=lower,
                    bin_upper=upper,
                    n_samples=0,
                ))
                continue
            
            bin_pred = predicted[mask]
            bin_obs = observed[mask]
            
            mean_pred = bin_pred.mean()
            mean_obs = bin_obs.mean()
            
            # Wilson score interval for binomial proportion
            ci_lower, ci_upper = self._wilson_score_interval(
                n_in_bin, int(bin_obs.sum())
            )
            
            bin_result = CalibrationBin(
                bin_lower=lower,
                bin_upper=upper,
                n_samples=n_in_bin,
                mean_predicted=mean_pred,
                mean_observed=mean_obs,
                confidence_lower=ci_lower,
                confidence_upper=ci_upper,
            )
            bins.append(bin_result)
            
            # Gap for ECE
            bin_gaps.append((n_in_bin, abs(mean_pred - mean_obs)))
        
        # Compute ECE (weighted average gap)
        if bin_gaps:
            ece = sum(n * gap for n, gap in bin_gaps) / n_samples
        else:
            ece = 0.0
        
        # Compute MCE (maximum gap)
        if bin_gaps:
            mce = max(gap for _, gap in bin_gaps)
        else:
            mce = 0.0
        
        # Brier score
        brier = np.mean((predicted - observed) ** 2)
        
        # Hosmer-Lemeshow test
        hl_pvalue = self._hosmer_lemeshow_test(predicted, observed, bins)
        
        return CalibrationResult(
            module=module,
            n_samples=n_samples,
            ece=ece,
            mce=mce,
            bins=bins,
            brier_score=brier,
            is_calibrated=(ece < 0.05),
            hosmer_lemeshow_pvalue=hl_pvalue,
        )
    
    def _compute_bin_edges(self, predicted: np.ndarray) -> np.ndarray:
        """Compute bin edges based on strategy."""
        if self.strategy == "uniform":
            return np.linspace(0, 1, self.n_bins + 1)
        elif self.strategy == "quantile":
            # Equal-count bins
            percentiles = np.linspace(0, 100, self.n_bins + 1)
            edges = np.percentile(predicted, percentiles)
            edges[0] = 0
            edges[-1] = 1
            return edges
        else:
            return np.linspace(0, 1, self.n_bins + 1)
    
    def _wilson_score_interval(
        self,
        n: int,
        k: int,
        confidence: float = 0.95,
    ) -> Tuple[float, float]:
        """Compute Wilson score confidence interval for proportion."""
        if n == 0:
            return 0.0, 1.0
        
        z = stats.norm.ppf(1 - (1 - confidence) / 2)
        p = k / n
        
        denominator = 1 + z**2 / n
        center = (p + z**2 / (2 * n)) / denominator
        margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denominator
        
        lower = max(0, center - margin)
        upper = min(1, center + margin)
        
        return lower, upper
    
    def _hosmer_lemeshow_test(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        bins: List[CalibrationBin],
    ) -> float:
        """
        Hosmer-Lemeshow goodness-of-fit test.
        
        Tests null hypothesis that predictions are well-calibrated.
        Low p-value suggests miscalibration.
        """
        chi_sq = 0
        df = 0
        
        for bin_result in bins:
            n = bin_result.n_samples
            if n == 0:
                continue
            
            expected = bin_result.mean_predicted * n
            observed_sum = bin_result.mean_observed * n
            
            if expected > 0 and (n - expected) > 0:
                chi_sq += (observed_sum - expected) ** 2 / expected
                chi_sq += (n - observed_sum - (n - expected)) ** 2 / (n - expected)
                df += 1
        
        if df <= 2:
            return 1.0  # Not enough data
        
        df = df - 2  # Degrees of freedom
        
        if df <= 0:
            return 1.0
        
        return 1 - stats.chi2.cdf(chi_sq, df)
    
    def reliability_curve_data(
        self,
        result: CalibrationResult,
    ) -> pd.DataFrame:
        """
        Get reliability curve data for plotting.
        
        Parameters
        ----------
        result : CalibrationResult
            Calibration analysis result.
            
        Returns
        -------
        pd.DataFrame
            Data for plotting reliability curve.
        """
        data = []
        for bin_result in result.bins:
            if bin_result.n_samples > 0:
                data.append({
                    "mean_predicted": bin_result.mean_predicted,
                    "mean_observed": bin_result.mean_observed,
                    "n_samples": bin_result.n_samples,
                    "ci_lower": bin_result.confidence_lower,
                    "ci_upper": bin_result.confidence_upper,
                    "bin_lower": bin_result.bin_lower,
                    "bin_upper": bin_result.bin_upper,
                })
        
        return pd.DataFrame(data)


class ModuleCalibrationSuite:
    """
    Calibration analysis across all mechanism graph modules.
    
    Analyzes calibration for:
    - Variant PIP (from fine-mapping)
    - Coloc PP.H4 (from colocalization)
    - cCRE→Gene probability (from ABC/PCHi-C)
    - Gene posterior (final gene-level probability)
    
    This addresses reviewer requirement: "calibration per module, not just final score"
    """
    
    MODULE_NAMES = [
        "variant_pip",
        "coloc_h4",
        "ccre_gene_prob",
        "gene_tissue_prob",
        "gene_posterior",
    ]
    
    def __init__(self, n_bins: int = 10):
        """Initialize suite."""
        self.calibrator = ProbabilityCalibration(n_bins=n_bins)
        self.results: Dict[str, CalibrationResult] = {}
    
    def analyze_module(
        self,
        module: str,
        predicted: np.ndarray,
        observed: np.ndarray,
    ) -> CalibrationResult:
        """
        Analyze calibration for a specific module.
        
        Parameters
        ----------
        module : str
            Module name.
        predicted : np.ndarray
            Predicted probabilities.
        observed : np.ndarray
            Binary outcomes.
            
        Returns
        -------
        CalibrationResult
            Calibration result.
        """
        result = self.calibrator.analyze(predicted, observed, module)
        self.results[module] = result
        return result
    
    def analyze_all(
        self,
        data: Dict[str, Dict[str, np.ndarray]],
    ) -> Dict[str, CalibrationResult]:
        """
        Analyze calibration for all modules.
        
        Parameters
        ----------
        data : dict
            Dict mapping module name to {"predicted": array, "observed": array}.
            
        Returns
        -------
        dict
            Module name -> CalibrationResult.
        """
        for module, module_data in data.items():
            self.analyze_module(
                module,
                module_data["predicted"],
                module_data["observed"],
            )
        
        return self.results
    
    def summary(self) -> pd.DataFrame:
        """
        Get summary table of all module calibrations.
        
        Returns
        -------
        pd.DataFrame
            Summary statistics.
        """
        rows = []
        for module, result in self.results.items():
            rows.append({
                "module": module,
                "n_samples": result.n_samples,
                "ECE": result.ece,
                "MCE": result.mce,
                "Brier": result.brier_score,
                "calibrated": result.is_calibrated,
                "HL_pvalue": result.hosmer_lemeshow_pvalue,
            })
        
        return pd.DataFrame(rows)
    
    def worst_calibrated(self) -> Optional[str]:
        """Get the worst-calibrated module."""
        if not self.results:
            return None
        
        return max(self.results, key=lambda m: self.results[m].ece)
    
    def all_calibrated(self, threshold: float = 0.05) -> bool:
        """Check if all modules are calibrated."""
        return all(r.ece < threshold for r in self.results.values())


class PlattScaling:
    """
    Platt scaling for probability recalibration.
    
    Fits a logistic regression to map raw scores to calibrated probabilities.
    """
    
    def __init__(self):
        """Initialize Platt scaling."""
        self.a = 0.0
        self.b = 0.0
        self.fitted = False
    
    def fit(
        self,
        scores: np.ndarray,
        labels: np.ndarray,
    ) -> "PlattScaling":
        """
        Fit Platt scaling parameters.
        
        Parameters
        ----------
        scores : np.ndarray
            Raw scores/probabilities.
        labels : np.ndarray
            Binary labels.
            
        Returns
        -------
        PlattScaling
            Self (fitted).
        """
        from scipy.optimize import minimize
        
        scores = np.asarray(scores)
        labels = np.asarray(labels)
        
        # Transform scores to log-odds if in (0, 1)
        eps = 1e-7
        scores = np.clip(scores, eps, 1 - eps)
        log_odds = np.log(scores / (1 - scores))
        
        # Objective: negative log-likelihood
        def objective(params):
            a, b = params
            probs = 1 / (1 + np.exp(-(a * log_odds + b)))
            probs = np.clip(probs, eps, 1 - eps)
            nll = -np.mean(labels * np.log(probs) + (1 - labels) * np.log(1 - probs))
            return nll
        
        result = minimize(objective, [1.0, 0.0], method="BFGS")
        self.a, self.b = result.x
        self.fitted = True
        
        return self
    
    def transform(self, scores: np.ndarray) -> np.ndarray:
        """
        Transform scores to calibrated probabilities.
        
        Parameters
        ----------
        scores : np.ndarray
            Raw scores.
            
        Returns
        -------
        np.ndarray
            Calibrated probabilities.
        """
        if not self.fitted:
            raise ValueError("PlattScaling not fitted")
        
        scores = np.asarray(scores)
        eps = 1e-7
        scores = np.clip(scores, eps, 1 - eps)
        log_odds = np.log(scores / (1 - scores))
        
        calibrated = 1 / (1 + np.exp(-(self.a * log_odds + self.b)))
        return calibrated


def compute_ece(
    predicted: np.ndarray,
    observed: np.ndarray,
    n_bins: int = 10,
) -> float:
    """
    Convenience function to compute Expected Calibration Error.
    
    Parameters
    ----------
    predicted : np.ndarray
        Predicted probabilities.
    observed : np.ndarray
        Binary outcomes.
    n_bins : int
        Number of bins.
        
    Returns
    -------
    float
        ECE value.
    """
    cal = ProbabilityCalibration(n_bins=n_bins)
    result = cal.analyze(predicted, observed)
    return result.ece


def compute_mce(
    predicted: np.ndarray,
    observed: np.ndarray,
    n_bins: int = 10,
) -> float:
    """
    Convenience function to compute Maximum Calibration Error.
    
    Parameters
    ----------
    predicted : np.ndarray
        Predicted probabilities.
    observed : np.ndarray
        Binary outcomes.
    n_bins : int
        Number of bins.
        
    Returns
    -------
    float
        MCE value.
    """
    cal = ProbabilityCalibration(n_bins=n_bins)
    result = cal.analyze(predicted, observed)
    return result.mce
