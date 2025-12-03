"""
Calibration Analysis

Evaluates whether predicted probabilities are well-calibrated
(i.e., a 70% probability should be correct ~70% of the time).
"""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from collections import defaultdict

from ..utils.logging import get_logger


logger = get_logger("calibration")


class CalibrationAnalysis:
    """
    Calibration analysis for probabilistic predictions.
    
    Computes:
    - Reliability diagrams
    - Expected Calibration Error (ECE)
    - Maximum Calibration Error (MCE)
    - Brier score decomposition
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
            Number of bins for reliability curve.
        strategy : str
            Binning strategy: "uniform" or "quantile".
        """
        self.n_bins = n_bins
        self.strategy = strategy
    
    def compute_reliability_curve(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """
        Compute reliability curve data.
        
        Parameters
        ----------
        y_true : np.ndarray
            True binary labels.
        y_prob : np.ndarray
            Predicted probabilities.
            
        Returns
        -------
        dict
            Reliability curve data.
        """
        y_true = np.asarray(y_true)
        y_prob = np.asarray(y_prob)
        
        if self.strategy == "uniform":
            bins = np.linspace(0, 1, self.n_bins + 1)
        else:  # quantile
            bins = np.percentile(y_prob, np.linspace(0, 100, self.n_bins + 1))
            bins = np.unique(bins)
        
        bin_indices = np.digitize(y_prob, bins[:-1]) - 1
        bin_indices = np.clip(bin_indices, 0, len(bins) - 2)
        
        bin_means = np.zeros(len(bins) - 1)
        bin_true_means = np.zeros(len(bins) - 1)
        bin_counts = np.zeros(len(bins) - 1)
        
        for i in range(len(bins) - 1):
            mask = bin_indices == i
            if mask.sum() > 0:
                bin_means[i] = y_prob[mask].mean()
                bin_true_means[i] = y_true[mask].mean()
                bin_counts[i] = mask.sum()
        
        return {
            "bin_edges": bins,
            "mean_predicted_prob": bin_means,
            "fraction_positives": bin_true_means,
            "bin_counts": bin_counts,
        }
    
    def expected_calibration_error(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
    ) -> float:
        """
        Compute Expected Calibration Error (ECE).
        
        ECE = Σ (n_i / N) * |acc_i - conf_i|
        
        Parameters
        ----------
        y_true : np.ndarray
            True labels.
        y_prob : np.ndarray
            Predicted probabilities.
            
        Returns
        -------
        float
            ECE value.
        """
        curve = self.compute_reliability_curve(y_true, y_prob)
        
        total = curve["bin_counts"].sum()
        if total == 0:
            return 0.0
        
        ece = 0.0
        for i in range(len(curve["bin_counts"])):
            if curve["bin_counts"][i] > 0:
                weight = curve["bin_counts"][i] / total
                gap = abs(curve["fraction_positives"][i] - curve["mean_predicted_prob"][i])
                ece += weight * gap
        
        return ece
    
    def maximum_calibration_error(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
    ) -> float:
        """
        Compute Maximum Calibration Error (MCE).
        
        Parameters
        ----------
        y_true : np.ndarray
            True labels.
        y_prob : np.ndarray
            Predicted probabilities.
            
        Returns
        -------
        float
            MCE value.
        """
        curve = self.compute_reliability_curve(y_true, y_prob)
        
        mce = 0.0
        for i in range(len(curve["bin_counts"])):
            if curve["bin_counts"][i] > 0:
                gap = abs(curve["fraction_positives"][i] - curve["mean_predicted_prob"][i])
                mce = max(mce, gap)
        
        return mce
    
    def brier_score_decomposition(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
    ) -> Dict[str, float]:
        """
        Decompose Brier score into reliability, resolution, uncertainty.
        
        Brier = Reliability - Resolution + Uncertainty
        
        Parameters
        ----------
        y_true : np.ndarray
            True labels.
        y_prob : np.ndarray
            Predicted probabilities.
            
        Returns
        -------
        dict
            Brier score components.
        """
        y_true = np.asarray(y_true)
        y_prob = np.asarray(y_prob)
        
        n = len(y_true)
        base_rate = y_true.mean()
        
        # Brier score
        brier = np.mean((y_prob - y_true) ** 2)
        
        # Uncertainty (irreducible)
        uncertainty = base_rate * (1 - base_rate)
        
        # Get reliability curve
        curve = self.compute_reliability_curve(y_true, y_prob)
        
        # Reliability (calibration)
        reliability = 0.0
        for i in range(len(curve["bin_counts"])):
            if curve["bin_counts"][i] > 0:
                weight = curve["bin_counts"][i] / n
                gap = (curve["mean_predicted_prob"][i] - curve["fraction_positives"][i]) ** 2
                reliability += weight * gap
        
        # Resolution (discrimination)
        resolution = 0.0
        for i in range(len(curve["bin_counts"])):
            if curve["bin_counts"][i] > 0:
                weight = curve["bin_counts"][i] / n
                gap = (curve["fraction_positives"][i] - base_rate) ** 2
                resolution += weight * gap
        
        return {
            "brier_score": brier,
            "reliability": reliability,
            "resolution": resolution,
            "uncertainty": uncertainty,
        }
    
    def analyze(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
    ) -> Dict[str, Any]:
        """
        Run full calibration analysis.
        
        Parameters
        ----------
        y_true : np.ndarray
            True labels.
        y_prob : np.ndarray
            Predicted probabilities.
            
        Returns
        -------
        dict
            Comprehensive calibration results.
        """
        return {
            "reliability_curve": self.compute_reliability_curve(y_true, y_prob),
            "ece": self.expected_calibration_error(y_true, y_prob),
            "mce": self.maximum_calibration_error(y_true, y_prob),
            "brier_decomposition": self.brier_score_decomposition(y_true, y_prob),
            "n_samples": len(y_true),
            "n_positives": int(y_true.sum()),
            "base_rate": float(y_true.mean()),
        }


def reliability_curve(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_bins: int = 10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convenience function for reliability curve.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels.
    y_prob : np.ndarray
        Predicted probabilities.
    n_bins : int
        Number of bins.
        
    Returns
    -------
    tuple
        (mean_predicted, fraction_positive)
    """
    analyzer = CalibrationAnalysis(n_bins=n_bins)
    curve = analyzer.compute_reliability_curve(y_true, y_prob)
    
    # Filter out empty bins
    mask = curve["bin_counts"] > 0
    
    return (
        curve["mean_predicted_prob"][mask],
        curve["fraction_positives"][mask],
    )


class PlattScaling:
    """
    Platt scaling for probability calibration.
    
    Fits a logistic regression to transform scores to probabilities.
    """
    
    def __init__(self):
        self.a = 0.0
        self.b = 1.0
        self._fitted = False
    
    def fit(
        self,
        y_true: np.ndarray,
        y_score: np.ndarray,
    ) -> "PlattScaling":
        """
        Fit Platt scaling parameters.
        
        Parameters
        ----------
        y_true : np.ndarray
            True labels.
        y_score : np.ndarray
            Raw scores or uncalibrated probabilities.
            
        Returns
        -------
        self
        """
        from scipy.optimize import minimize
        
        y_true = np.asarray(y_true, dtype=float)
        y_score = np.asarray(y_score, dtype=float)
        
        # Logit transform if scores are probabilities
        eps = 1e-10
        y_score = np.clip(y_score, eps, 1 - eps)
        
        def loss(params):
            a, b = params
            logits = a + b * y_score
            probs = 1 / (1 + np.exp(-logits))
            # Cross-entropy loss
            return -np.mean(y_true * np.log(probs + eps) + (1 - y_true) * np.log(1 - probs + eps))
        
        result = minimize(loss, [0.0, 1.0], method='BFGS')
        
        self.a, self.b = result.x
        self._fitted = True
        
        return self
    
    def transform(self, y_score: np.ndarray) -> np.ndarray:
        """
        Transform scores to calibrated probabilities.
        
        Parameters
        ----------
        y_score : np.ndarray
            Raw scores.
            
        Returns
        -------
        np.ndarray
            Calibrated probabilities.
        """
        if not self._fitted:
            raise ValueError("Must fit before transform")
        
        y_score = np.asarray(y_score, dtype=float)
        logits = self.a + self.b * y_score
        return 1 / (1 + np.exp(-logits))


class IsotonicCalibration:
    """
    Isotonic regression for probability calibration.
    
    Non-parametric calibration that preserves ordering.
    """
    
    def __init__(self):
        self._fitted = False
        self._ir = None
    
    def fit(
        self,
        y_true: np.ndarray,
        y_score: np.ndarray,
    ) -> "IsotonicCalibration":
        """
        Fit isotonic calibration.
        
        Parameters
        ----------
        y_true : np.ndarray
            True labels.
        y_score : np.ndarray
            Raw scores.
            
        Returns
        -------
        self
        """
        from scipy.interpolate import interp1d
        
        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)
        
        # Sort by score
        order = np.argsort(y_score)
        y_score_sorted = y_score[order]
        y_true_sorted = y_true[order]
        
        # Pool Adjacent Violators Algorithm (PAVA)
        calibrated = self._pava(y_true_sorted)
        
        # Create interpolation function
        self._ir = interp1d(
            y_score_sorted,
            calibrated,
            kind='linear',
            bounds_error=False,
            fill_value=(calibrated[0], calibrated[-1]),
        )
        
        self._fitted = True
        return self
    
    def _pava(self, y: np.ndarray) -> np.ndarray:
        """Pool Adjacent Violators Algorithm."""
        n = len(y)
        y = y.copy().astype(float)
        
        # Block-based PAVA
        blocks = list(range(n))
        
        while True:
            changed = False
            i = 0
            while i < len(blocks) - 1:
                # Get block means
                if i == 0:
                    start1 = 0
                else:
                    start1 = blocks[i-1] + 1 if i > 0 else 0
                end1 = blocks[i] + 1
                
                start2 = end1
                end2 = blocks[i+1] + 1 if i+1 < len(blocks) else n
                
                mean1 = y[start1:end1].mean()
                mean2 = y[start2:end2].mean()
                
                if mean1 > mean2:
                    # Merge blocks
                    combined_mean = y[start1:end2].mean()
                    y[start1:end2] = combined_mean
                    blocks.pop(i+1)
                    changed = True
                else:
                    i += 1
            
            if not changed:
                break
        
        return y
    
    def transform(self, y_score: np.ndarray) -> np.ndarray:
        """
        Transform scores to calibrated probabilities.
        
        Parameters
        ----------
        y_score : np.ndarray
            Raw scores.
            
        Returns
        -------
        np.ndarray
            Calibrated probabilities.
        """
        if not self._fitted:
            raise ValueError("Must fit before transform")
        
        y_score = np.asarray(y_score)
        return np.clip(self._ir(y_score), 0, 1)
