"""
Quality Control for GWAS Summary Statistics

Implements comprehensive QC checks including:
- Z-score consistency
- P-value recomputation
- Genomic inflation
- Effect size distribution
"""

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from ..utils.logging import get_logger


logger = get_logger("qc")


class QualityControl:
    """
    Quality control checks for GWAS summary statistics.
    """
    
    def __init__(
        self,
        z_threshold: float = 37.0,
        max_se: float = 10.0,
        min_n: int = 1000,
        max_pval_error: float = 0.1,
    ):
        """
        Initialize QC checks.
        
        Parameters
        ----------
        z_threshold : float
            Maximum plausible |Z| score.
        max_se : float
            Maximum plausible standard error.
        min_n : int
            Minimum sample size.
        max_pval_error : float
            Maximum relative error for p-value recomputation.
        """
        self.z_threshold = z_threshold
        self.max_se = max_se
        self.min_n = min_n
        self.max_pval_error = max_pval_error
        
    def check_z_score_consistency(
        self,
        df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Check Z = beta/SE consistency.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with beta, se, z columns.
            
        Returns
        -------
        dict
            QC results including flag counts and statistics.
        """
        results = {
            "check": "z_score_consistency",
            "passed": True,
            "n_variants": len(df),
        }
        
        # Calculate Z from beta/SE
        z_calc = df["beta"] / df["se"]
        
        # Compare to reported Z (if available)
        if "z" in df.columns:
            z_diff = np.abs(df["z"] - z_calc)
            results["z_diff_mean"] = z_diff.mean()
            results["z_diff_max"] = z_diff.max()
            results["n_inconsistent"] = (z_diff > 0.1).sum()
            
            if results["n_inconsistent"] > 0.01 * len(df):
                results["passed"] = False
                results["message"] = f"{results['n_inconsistent']} variants with Z inconsistency"
        
        # Check for extreme Z values
        z_extreme = np.abs(z_calc) > self.z_threshold
        results["n_extreme_z"] = z_extreme.sum()
        
        if results["n_extreme_z"] > 0:
            logger.warning(f"{results['n_extreme_z']} variants with |Z| > {self.z_threshold}")
        
        return results
    
    def check_pvalue_consistency(
        self,
        df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Check p-value consistency with Z-score.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with pval, beta, se columns.
            
        Returns
        -------
        dict
            QC results.
        """
        results = {
            "check": "pvalue_consistency",
            "passed": True,
            "n_variants": len(df),
        }
        
        # Calculate expected p-value from Z
        z = df["beta"] / df["se"]
        pval_expected = 2 * stats.norm.sf(np.abs(z))
        
        # Compare to reported p-value
        # Use log scale for very small p-values
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            log_pval_reported = np.log10(df["pval"].clip(lower=1e-300))
            log_pval_expected = np.log10(pval_expected.clip(lower=1e-300))
        
        log_diff = np.abs(log_pval_reported - log_pval_expected)
        
        # Allow larger differences for very small p-values
        relative_error = log_diff / np.abs(log_pval_expected).clip(lower=1)
        
        results["log_pval_diff_median"] = log_diff.median()
        results["log_pval_diff_max"] = log_diff.max()
        results["n_inconsistent"] = (relative_error > self.max_pval_error).sum()
        
        if results["n_inconsistent"] > 0.01 * len(df):
            results["passed"] = False
            results["message"] = f"{results['n_inconsistent']} variants with p-value inconsistency"
        
        return results
    
    def check_genomic_inflation(
        self,
        df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Calculate genomic inflation factor (lambda_GC).
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with pval column.
            
        Returns
        -------
        dict
            QC results including lambda_GC.
        """
        results = {
            "check": "genomic_inflation",
            "passed": True,
            "n_variants": len(df),
        }
        
        # Calculate chi-square statistics
        z = df["beta"] / df["se"]
        chi2 = z ** 2
        
        # Lambda GC is median chi2 / expected median
        lambda_gc = np.median(chi2) / stats.chi2.ppf(0.5, df=1)
        
        results["lambda_gc"] = lambda_gc
        
        # Flag if lambda is too high or too low
        if lambda_gc > 1.2:
            results["passed"] = False
            results["message"] = f"High genomic inflation (lambda = {lambda_gc:.3f})"
            logger.warning(results["message"])
        elif lambda_gc < 0.9:
            results["passed"] = False
            results["message"] = f"Deflation detected (lambda = {lambda_gc:.3f})"
            logger.warning(results["message"])
        
        return results
    
    def check_effect_distribution(
        self,
        df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Check effect size distribution.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with beta column.
            
        Returns
        -------
        dict
            QC results including distribution statistics.
        """
        results = {
            "check": "effect_distribution",
            "passed": True,
            "n_variants": len(df),
        }
        
        beta = df["beta"].dropna()
        
        results["beta_mean"] = beta.mean()
        results["beta_std"] = beta.std()
        results["beta_median"] = beta.median()
        results["beta_min"] = beta.min()
        results["beta_max"] = beta.max()
        
        # Check for extreme effects
        n_extreme = (np.abs(beta) > 5).sum()
        results["n_extreme_beta"] = n_extreme
        
        if n_extreme > 0.001 * len(beta):
            logger.warning(f"{n_extreme} variants with |beta| > 5")
        
        # Check for suspiciously uniform distribution
        # (could indicate synthetic data or major errors)
        _, pval_uniform = stats.kstest(beta, 'uniform', args=(beta.min(), beta.max() - beta.min()))
        results["ks_pval_uniform"] = pval_uniform
        
        if pval_uniform > 0.05:
            results["passed"] = False
            results["message"] = "Effect sizes suspiciously uniform"
        
        return results
    
    def check_standard_errors(
        self,
        df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Check standard error distribution.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with se column.
            
        Returns
        -------
        dict
            QC results.
        """
        results = {
            "check": "standard_errors",
            "passed": True,
            "n_variants": len(df),
        }
        
        se = df["se"].dropna()
        
        results["se_mean"] = se.mean()
        results["se_std"] = se.std()
        results["se_median"] = se.median()
        results["se_min"] = se.min()
        results["se_max"] = se.max()
        
        # Check for zero or negative SE
        n_invalid = (se <= 0).sum()
        results["n_invalid_se"] = n_invalid
        
        if n_invalid > 0:
            results["passed"] = False
            results["message"] = f"{n_invalid} variants with SE <= 0"
        
        # Check for extreme SE
        n_extreme = (se > self.max_se).sum()
        results["n_extreme_se"] = n_extreme
        
        return results
    
    def check_sample_size(
        self,
        df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Check sample size distribution.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with n column.
            
        Returns
        -------
        dict
            QC results.
        """
        results = {
            "check": "sample_size",
            "passed": True,
            "n_variants": len(df),
        }
        
        if "n" not in df.columns or df["n"].isna().all():
            results["passed"] = True
            results["message"] = "No sample size information"
            return results
        
        n = df["n"].dropna()
        
        results["n_mean"] = n.mean()
        results["n_min"] = n.min()
        results["n_max"] = n.max()
        results["n_unique"] = n.nunique()
        
        if n.min() < self.min_n:
            results["passed"] = False
            results["message"] = f"Minimum sample size ({n.min()}) below threshold ({self.min_n})"
        
        return results
    
    def run_all_checks(
        self,
        df: pd.DataFrame,
        verbose: bool = True,
    ) -> Dict[str, Dict[str, Any]]:
        """
        Run all QC checks.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics.
        verbose : bool
            Whether to print results.
            
        Returns
        -------
        dict
            All QC results.
        """
        results = {}
        
        checks = [
            ("z_score", self.check_z_score_consistency),
            ("pvalue", self.check_pvalue_consistency),
            ("genomic_inflation", self.check_genomic_inflation),
            ("effect_distribution", self.check_effect_distribution),
            ("standard_errors", self.check_standard_errors),
            ("sample_size", self.check_sample_size),
        ]
        
        all_passed = True
        
        for name, check_func in checks:
            try:
                result = check_func(df)
                results[name] = result
                
                if not result["passed"]:
                    all_passed = False
                    
                if verbose:
                    status = "✓" if result["passed"] else "✗"
                    msg = result.get("message", "OK")
                    logger.info(f"  {status} {name}: {msg}")
                    
            except Exception as e:
                results[name] = {
                    "check": name,
                    "passed": False,
                    "error": str(e),
                }
                all_passed = False
                if verbose:
                    logger.error(f"  ✗ {name}: Error - {e}")
        
        results["overall_passed"] = all_passed
        
        return results


def run_qc(
    df: pd.DataFrame,
    **kwargs,
) -> Dict[str, Dict[str, Any]]:
    """
    Convenience function to run all QC checks.
    
    Parameters
    ----------
    df : pd.DataFrame
        Summary statistics.
    **kwargs
        Arguments for QualityControl.
        
    Returns
    -------
    dict
        QC results.
    """
    qc = QualityControl(**kwargs)
    return qc.run_all_checks(df)
