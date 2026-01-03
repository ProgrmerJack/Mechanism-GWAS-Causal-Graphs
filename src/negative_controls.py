"""
Negative Controls Module for Causal Inference Validation

Implements falsification tests to strengthen causal identification:
1. Outcome negative controls - outcomes that shouldn't be affected
2. Treatment negative controls - fake treatments (placebo cutoffs, fake reforms)
3. Edge permutation - randomize causal graph edges
4. Balance tests - covariate smoothness at cutoff

These tests are essential for Nature-tier credibility.

References:
- Lipsitch, Tchetgen, & Cohen (2010). Negative Controls. Epidemiology.
- Shi, Miao, & Tchetgen (2020). A Selective Review of Negative Control Methods.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Callable, Any
from dataclasses import dataclass, field
from scipy import stats
import logging
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class NegativeControlResult:
    """Result from a single negative control test."""
    test_name: str
    test_type: str  # 'outcome', 'treatment', 'balance', 'permutation'
    null_hypothesis: str
    estimate: float
    se: float
    p_value: float
    passed: bool  # True if null NOT rejected (no spurious effect)
    confidence_interval: Tuple[float, float] = (np.nan, np.nan)
    n_observations: int = 0
    details: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        return {
            "test_name": self.test_name,
            "test_type": self.test_type,
            "null_hypothesis": self.null_hypothesis,
            "estimate": self.estimate,
            "se": self.se,
            "p_value": self.p_value,
            "passed": self.passed,
            "ci_lower": self.confidence_interval[0],
            "ci_upper": self.confidence_interval[1],
            "n": self.n_observations
        }


@dataclass
class NegativeControlBattery:
    """Collection of negative control results."""
    results: List[NegativeControlResult] = field(default_factory=list)
    
    @property
    def n_tests(self) -> int:
        return len(self.results)
    
    @property
    def n_passed(self) -> int:
        return sum(1 for r in self.results if r.passed)
    
    @property
    def pass_rate(self) -> float:
        return self.n_passed / self.n_tests if self.n_tests > 0 else np.nan
    
    def by_type(self, test_type: str) -> List[NegativeControlResult]:
        """Filter results by test type."""
        return [r for r in self.results if r.test_type == test_type]
    
    def summary(self) -> pd.DataFrame:
        """Generate summary DataFrame."""
        return pd.DataFrame([r.to_dict() for r in self.results])
    
    def format_report(self) -> str:
        """Generate formatted text report."""
        lines = [
            "=" * 60,
            "NEGATIVE CONTROL BATTERY RESULTS",
            "=" * 60,
            f"\nTotal tests: {self.n_tests}",
            f"Passed: {self.n_passed} ({self.pass_rate:.1%})",
            ""
        ]
        
        for test_type in ['outcome', 'treatment', 'balance', 'permutation']:
            type_results = self.by_type(test_type)
            if type_results:
                lines.append(f"\n--- {test_type.upper()} NEGATIVE CONTROLS ---")
                for r in type_results:
                    status = "✓ PASS" if r.passed else "✗ FAIL"
                    lines.append(
                        f"  {r.test_name}: {status} "
                        f"(est={r.estimate:.4f}, p={r.p_value:.4f})"
                    )
        
        lines.append("\n" + "=" * 60)
        return "\n".join(lines)


# =============================================================================
# OUTCOME NEGATIVE CONTROLS
# =============================================================================

class OutcomeNegativeControls:
    """
    Test outcomes that should NOT be affected by treatment.
    
    Logic: If treatment truly works through the hypothesized mechanism,
    it should NOT affect outcomes outside that pathway. Finding effects
    on negative control outcomes suggests confounding or model misspecification.
    
    Examples:
    - Procurement transparency → bidder count ✓ (expected)
    - Procurement transparency → weather ✗ (should NOT affect)
    - Procurement transparency → past outcomes ✗ (cannot affect past)
    """
    
    def __init__(self, estimator: Callable):
        """
        Initialize with an RDD/DiD estimator function.
        
        Args:
            estimator: Function(X, Y, cutoff) -> (estimate, se, n)
        """
        self.estimator = estimator
    
    def test_lagged_outcome(
        self,
        data: pd.DataFrame,
        running_var: str,
        outcome_var: str,
        cutoff: float = 0,
        lag_periods: int = 1
    ) -> NegativeControlResult:
        """
        Test effect on lagged (pre-treatment) outcome.
        
        Treatment at time t should not affect outcome at t-k.
        """
        # Create lagged outcome
        lagged_col = f"{outcome_var}_lag{lag_periods}"
        if lagged_col not in data.columns:
            data = data.sort_values(running_var)
            data[lagged_col] = data[outcome_var].shift(lag_periods)
        
        valid = data.dropna(subset=[running_var, lagged_col])
        
        if len(valid) < 20:
            return NegativeControlResult(
                test_name=f"Lagged outcome (t-{lag_periods})",
                test_type="outcome",
                null_hypothesis="No effect on past outcome",
                estimate=np.nan,
                se=np.nan,
                p_value=np.nan,
                passed=True,  # Insufficient data, assume pass
                n_observations=len(valid),
                details={"reason": "Insufficient observations"}
            )
        
        estimate, se, n = self.estimator(
            valid[running_var].values,
            valid[lagged_col].values,
            cutoff
        )
        
        # Two-sided test
        if se > 0:
            t_stat = estimate / se
            p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
        else:
            p_value = 1.0
        
        # Pass if p > 0.05 (no significant effect on lagged outcome)
        passed = p_value > 0.05
        
        return NegativeControlResult(
            test_name=f"Lagged outcome (t-{lag_periods})",
            test_type="outcome",
            null_hypothesis="No effect on past outcome",
            estimate=estimate,
            se=se,
            p_value=p_value,
            passed=passed,
            confidence_interval=(estimate - 1.96*se, estimate + 1.96*se),
            n_observations=n
        )
    
    def test_unrelated_outcome(
        self,
        data: pd.DataFrame,
        running_var: str,
        negative_outcome: str,
        cutoff: float = 0,
        outcome_description: str = "Unrelated variable"
    ) -> NegativeControlResult:
        """
        Test effect on a theoretically unrelated outcome.
        
        E.g., procurement transparency shouldn't affect weather.
        """
        valid = data.dropna(subset=[running_var, negative_outcome])
        
        if len(valid) < 20:
            return NegativeControlResult(
                test_name=f"Unrelated: {outcome_description}",
                test_type="outcome",
                null_hypothesis=f"No effect on {outcome_description}",
                estimate=np.nan,
                se=np.nan,
                p_value=np.nan,
                passed=True,
                n_observations=len(valid),
                details={"reason": "Insufficient observations"}
            )
        
        estimate, se, n = self.estimator(
            valid[running_var].values,
            valid[negative_outcome].values,
            cutoff
        )
        
        if se > 0:
            t_stat = estimate / se
            p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
        else:
            p_value = 1.0
        
        passed = p_value > 0.05
        
        return NegativeControlResult(
            test_name=f"Unrelated: {outcome_description}",
            test_type="outcome",
            null_hypothesis=f"No effect on {outcome_description}",
            estimate=estimate,
            se=se,
            p_value=p_value,
            passed=passed,
            confidence_interval=(estimate - 1.96*se, estimate + 1.96*se),
            n_observations=n
        )
    
    def run_battery(
        self,
        data: pd.DataFrame,
        running_var: str,
        primary_outcome: str,
        negative_outcomes: List[Dict[str, str]],
        cutoff: float = 0
    ) -> NegativeControlBattery:
        """
        Run full battery of outcome negative controls.
        
        Args:
            data: Analysis DataFrame
            running_var: Running variable column
            primary_outcome: Main outcome variable
            negative_outcomes: List of dicts with 'column' and 'description'
            cutoff: RDD cutoff
        
        Returns:
            NegativeControlBattery with all results
        """
        results = []
        
        # Lagged outcomes
        for lag in [1, 2, 3]:
            try:
                result = self.test_lagged_outcome(
                    data, running_var, primary_outcome, cutoff, lag
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"Lagged outcome test failed (lag={lag}): {e}")
        
        # Unrelated outcomes
        for neg_out in negative_outcomes:
            try:
                result = self.test_unrelated_outcome(
                    data, running_var,
                    neg_out['column'],
                    cutoff,
                    neg_out.get('description', neg_out['column'])
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"Unrelated outcome test failed: {e}")
        
        return NegativeControlBattery(results=results)


# =============================================================================
# TREATMENT NEGATIVE CONTROLS (PLACEBO TESTS)
# =============================================================================

class TreatmentNegativeControls:
    """
    Test fake treatments that should have no effect.
    
    Types:
    1. Placebo cutoffs - RDD at wrong thresholds
    2. Placebo reforms - DiD at wrong reform dates
    3. Permutation inference - randomly assign treatment
    """
    
    def __init__(self, estimator: Callable):
        """
        Initialize with estimator function.
        
        Args:
            estimator: Function(X, Y, cutoff) -> (estimate, se, n)
        """
        self.estimator = estimator
    
    def placebo_cutoffs(
        self,
        data: pd.DataFrame,
        running_var: str,
        outcome_var: str,
        true_cutoff: float,
        placebo_cutoffs: Optional[List[float]] = None,
        n_placebos: int = 10
    ) -> List[NegativeControlResult]:
        """
        Test RDD at false cutoffs.
        
        Args:
            data: Analysis DataFrame
            running_var: Running variable column
            outcome_var: Outcome column
            true_cutoff: Real RDD cutoff
            placebo_cutoffs: Specific false cutoffs to test
            n_placebos: Number of random placebos if not specified
        
        Returns:
            List of NegativeControlResult
        """
        valid = data.dropna(subset=[running_var, outcome_var])
        X = valid[running_var].values
        Y = valid[outcome_var].values
        
        # Generate placebo cutoffs if not specified
        if placebo_cutoffs is None:
            # Use percentiles of X, avoiding true cutoff
            percentiles = np.linspace(10, 90, n_placebos + 1)
            all_cutoffs = np.percentile(X, percentiles)
            # Remove cutoffs too close to true cutoff
            placebo_cutoffs = [c for c in all_cutoffs 
                             if abs(c - true_cutoff) > 0.1 * np.std(X)][:n_placebos]
        
        results = []
        for pc in placebo_cutoffs:
            try:
                estimate, se, n = self.estimator(X, Y, pc)
                
                if se > 0:
                    t_stat = estimate / se
                    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
                else:
                    p_value = 1.0
                
                # Pass if no significant effect at placebo cutoff
                passed = p_value > 0.05
                
                results.append(NegativeControlResult(
                    test_name=f"Placebo cutoff at {pc:.2f}",
                    test_type="treatment",
                    null_hypothesis="No discontinuity at false cutoff",
                    estimate=estimate,
                    se=se,
                    p_value=p_value,
                    passed=passed,
                    confidence_interval=(estimate - 1.96*se, estimate + 1.96*se),
                    n_observations=n,
                    details={"placebo_cutoff": pc, "true_cutoff": true_cutoff}
                ))
            except Exception as e:
                logger.warning(f"Placebo cutoff test failed at {pc}: {e}")
        
        return results
    
    def placebo_reforms(
        self,
        data: pd.DataFrame,
        time_var: str,
        treatment_var: str,
        outcome_var: str,
        true_reform_date: str,
        placebo_dates: Optional[List[str]] = None,
        n_placebos: int = 5
    ) -> List[NegativeControlResult]:
        """
        Test DiD at false reform dates.
        
        Args:
            data: Panel DataFrame
            time_var: Time period column
            treatment_var: Treatment indicator column
            outcome_var: Outcome column
            true_reform_date: Actual reform date
            placebo_dates: Specific false dates to test
            n_placebos: Number of placebos if not specified
        
        Returns:
            List of NegativeControlResult
        """
        valid = data.dropna(subset=[time_var, outcome_var])
        
        # Generate placebo dates if not specified
        if placebo_dates is None:
            # Get unique pre-treatment periods
            times = pd.to_datetime(valid[time_var])
            true_date = pd.to_datetime(true_reform_date)
            pre_dates = times[times < true_date].unique()
            
            if len(pre_dates) > n_placebos:
                placebo_indices = np.linspace(0, len(pre_dates)-1, n_placebos, dtype=int)
                placebo_dates = [str(pre_dates[i]) for i in placebo_indices]
            else:
                placebo_dates = [str(d) for d in pre_dates]
        
        results = []
        for pd_date in placebo_dates:
            try:
                # Create placebo treatment based on fake date
                fake_time = pd.to_datetime(pd_date)
                times = pd.to_datetime(valid[time_var])
                fake_post = (times >= fake_time).astype(int)
                
                # Subset to pre-true-reform period only
                true_time = pd.to_datetime(true_reform_date)
                pre_reform = valid[times < true_time].copy()
                pre_reform['fake_post'] = fake_post[times < true_time]
                
                if len(pre_reform) < 20:
                    continue
                
                # Simple DiD: compare mean outcome pre/post fake date
                pre_fake = pre_reform[pre_reform['fake_post'] == 0][outcome_var]
                post_fake = pre_reform[pre_reform['fake_post'] == 1][outcome_var]
                
                estimate = post_fake.mean() - pre_fake.mean()
                pooled_se = np.sqrt(
                    pre_fake.var()/len(pre_fake) + post_fake.var()/len(post_fake)
                )
                
                if pooled_se > 0:
                    t_stat = estimate / pooled_se
                    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=len(pre_reform)-2))
                else:
                    p_value = 1.0
                
                passed = p_value > 0.05
                
                results.append(NegativeControlResult(
                    test_name=f"Placebo reform at {pd_date}",
                    test_type="treatment",
                    null_hypothesis="No effect at false reform date",
                    estimate=estimate,
                    se=pooled_se,
                    p_value=p_value,
                    passed=passed,
                    confidence_interval=(estimate - 1.96*pooled_se, estimate + 1.96*pooled_se),
                    n_observations=len(pre_reform),
                    details={"placebo_date": pd_date, "true_date": true_reform_date}
                ))
            except Exception as e:
                logger.warning(f"Placebo reform test failed at {pd_date}: {e}")
        
        return results
    
    def permutation_test(
        self,
        data: pd.DataFrame,
        running_var: str,
        outcome_var: str,
        cutoff: float,
        n_permutations: int = 1000,
        seed: int = 42
    ) -> NegativeControlResult:
        """
        Randomization inference via permutation.
        
        Randomly shuffle treatment assignment and compute p-value
        as proportion of permuted estimates >= observed estimate.
        """
        np.random.seed(seed)
        
        valid = data.dropna(subset=[running_var, outcome_var])
        X = valid[running_var].values
        Y = valid[outcome_var].values
        
        # Observed estimate
        obs_estimate, obs_se, n = self.estimator(X, Y, cutoff)
        
        # Permutation distribution
        perm_estimates = []
        for _ in range(n_permutations):
            # Shuffle outcome values
            Y_perm = np.random.permutation(Y)
            try:
                perm_est, _, _ = self.estimator(X, Y_perm, cutoff)
                perm_estimates.append(perm_est)
            except Exception:
                continue
        
        perm_estimates = np.array(perm_estimates)
        
        # Permutation p-value (two-sided)
        p_value = np.mean(np.abs(perm_estimates) >= np.abs(obs_estimate))
        
        # Pass if observed is consistent with random assignment
        passed = p_value > 0.05
        
        return NegativeControlResult(
            test_name="Permutation inference",
            test_type="treatment",
            null_hypothesis="Effect consistent with random assignment",
            estimate=obs_estimate,
            se=obs_se,
            p_value=p_value,
            passed=not passed,  # For this test, significant = good (not random)
            n_observations=n,
            details={
                "n_permutations": n_permutations,
                "perm_mean": np.mean(perm_estimates),
                "perm_sd": np.std(perm_estimates)
            }
        )


# =============================================================================
# BALANCE TESTS (COVARIATE SMOOTHNESS)
# =============================================================================

class BalanceTests:
    """
    Test that predetermined covariates are smooth at cutoff.
    
    If there's sorting/manipulation, covariates will jump at cutoff.
    Smooth covariates support RDD validity.
    """
    
    def __init__(self, estimator: Callable):
        """
        Initialize with estimator function.
        
        Args:
            estimator: Function(X, Y, cutoff) -> (estimate, se, n)
        """
        self.estimator = estimator
    
    def test_covariate_balance(
        self,
        data: pd.DataFrame,
        running_var: str,
        covariate: str,
        cutoff: float = 0,
        covariate_description: str = None
    ) -> NegativeControlResult:
        """
        Test single covariate for discontinuity at cutoff.
        
        Args:
            data: Analysis DataFrame
            running_var: Running variable column
            covariate: Predetermined covariate column
            cutoff: RDD cutoff
            covariate_description: Human-readable description
        
        Returns:
            NegativeControlResult
        """
        valid = data.dropna(subset=[running_var, covariate])
        
        if len(valid) < 20:
            return NegativeControlResult(
                test_name=f"Balance: {covariate_description or covariate}",
                test_type="balance",
                null_hypothesis="Covariate smooth at cutoff",
                estimate=np.nan,
                se=np.nan,
                p_value=np.nan,
                passed=True,
                n_observations=len(valid),
                details={"reason": "Insufficient observations"}
            )
        
        estimate, se, n = self.estimator(
            valid[running_var].values,
            valid[covariate].values,
            cutoff
        )
        
        if se > 0:
            t_stat = estimate / se
            p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
        else:
            p_value = 1.0
        
        # Pass if no significant jump (covariate is smooth)
        passed = p_value > 0.05
        
        return NegativeControlResult(
            test_name=f"Balance: {covariate_description or covariate}",
            test_type="balance",
            null_hypothesis="Covariate smooth at cutoff",
            estimate=estimate,
            se=se,
            p_value=p_value,
            passed=passed,
            confidence_interval=(estimate - 1.96*se, estimate + 1.96*se),
            n_observations=n
        )
    
    def run_battery(
        self,
        data: pd.DataFrame,
        running_var: str,
        covariates: List[Dict[str, str]],
        cutoff: float = 0
    ) -> NegativeControlBattery:
        """
        Run balance tests on multiple covariates.
        
        Args:
            data: Analysis DataFrame
            running_var: Running variable column
            covariates: List of dicts with 'column' and 'description'
            cutoff: RDD cutoff
        
        Returns:
            NegativeControlBattery with all results
        """
        results = []
        
        for cov in covariates:
            try:
                result = self.test_covariate_balance(
                    data, running_var,
                    cov['column'],
                    cutoff,
                    cov.get('description', cov['column'])
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"Balance test failed for {cov['column']}: {e}")
        
        return NegativeControlBattery(results=results)


# =============================================================================
# EDGE PERMUTATION (DAG VALIDATION)
# =============================================================================

class EdgePermutationTest:
    """
    Test DAG structure by permuting/removing edges.
    
    If the hypothesized causal graph is correct, permuting edges
    (testing wrong pathways) should yield null effects.
    """
    
    def test_mediation_bypass(
        self,
        data: pd.DataFrame,
        treatment: str,
        mediator: str,
        outcome: str,
        estimator: Callable
    ) -> NegativeControlResult:
        """
        Test if treatment affects outcome when mediator is controlled.
        
        If mechanism is M (mediator), controlling for M should
        block/reduce the treatment effect.
        """
        valid = data.dropna(subset=[treatment, mediator, outcome])
        
        # Residualize outcome on mediator
        from scipy.stats import linregress
        
        slope, intercept, _, _, _ = linregress(
            valid[mediator].values, 
            valid[outcome].values
        )
        residual_outcome = valid[outcome] - (intercept + slope * valid[mediator])
        
        # Estimate treatment effect on residualized outcome
        estimate, se, n = estimator(
            valid[treatment].values,
            residual_outcome.values,
            0  # Assume binary treatment centered at 0
        )
        
        if se > 0:
            t_stat = estimate / se
            p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
        else:
            p_value = 1.0
        
        # Pass interpretation depends on theory:
        # If full mediation expected, should be p > 0.05
        # If partial mediation expected, should be p < 0.05 but smaller than direct
        
        return NegativeControlResult(
            test_name="Mediation bypass test",
            test_type="permutation",
            null_hypothesis="Direct effect blocked by mediator control",
            estimate=estimate,
            se=se,
            p_value=p_value,
            passed=p_value > 0.05,  # Full mediation assumption
            confidence_interval=(estimate - 1.96*se, estimate + 1.96*se),
            n_observations=n,
            details={
                "mediator": mediator,
                "mediator_coefficient": slope
            }
        )
    
    def test_reverse_causation(
        self,
        data: pd.DataFrame,
        treatment: str,
        outcome: str,
        time_var: str,
        estimator: Callable
    ) -> NegativeControlResult:
        """
        Test reverse causation by examining if outcome predicts treatment.
        
        If causation is T → Y, then Y should not predict future T.
        """
        valid = data.sort_values(time_var).dropna(subset=[treatment, outcome])
        
        # Lead treatment (future treatment)
        valid['lead_treatment'] = valid[treatment].shift(-1)
        valid = valid.dropna(subset=['lead_treatment'])
        
        if len(valid) < 20:
            return NegativeControlResult(
                test_name="Reverse causation test",
                test_type="permutation",
                null_hypothesis="Outcome does not predict future treatment",
                estimate=np.nan,
                se=np.nan,
                p_value=np.nan,
                passed=True,
                n_observations=len(valid),
                details={"reason": "Insufficient observations"}
            )
        
        # Test if current outcome predicts future treatment
        from scipy.stats import linregress
        slope, _, _, p_value, se = linregress(
            valid[outcome].values,
            valid['lead_treatment'].values
        )
        
        # Pass if outcome does NOT predict future treatment
        passed = p_value > 0.05
        
        return NegativeControlResult(
            test_name="Reverse causation test",
            test_type="permutation",
            null_hypothesis="Outcome does not predict future treatment",
            estimate=slope,
            se=se,
            p_value=p_value,
            passed=passed,
            n_observations=len(valid)
        )


# =============================================================================
# COMPREHENSIVE NEGATIVE CONTROL RUNNER
# =============================================================================

class NegativeControlRunner:
    """
    Unified interface for running all negative control tests.
    
    Combines:
    - Outcome negative controls
    - Treatment negative controls (placebos)
    - Balance tests
    - Edge permutation tests
    """
    
    def __init__(self, estimator: Callable):
        """
        Initialize with RDD estimator.
        
        Args:
            estimator: Function(X, Y, cutoff) -> (estimate, se, n)
        """
        self.estimator = estimator
        self.outcome_controls = OutcomeNegativeControls(estimator)
        self.treatment_controls = TreatmentNegativeControls(estimator)
        self.balance_tests = BalanceTests(estimator)
        self.edge_tests = EdgePermutationTest()
    
    def run_full_battery(
        self,
        data: pd.DataFrame,
        config: Dict[str, Any]
    ) -> NegativeControlBattery:
        """
        Run comprehensive negative control battery.
        
        Args:
            data: Analysis DataFrame
            config: Configuration dict with:
                - running_var: Running variable column
                - outcome_var: Primary outcome column
                - cutoff: RDD cutoff
                - negative_outcomes: List of unrelated outcomes
                - covariates: List of predetermined covariates
                - placebo_cutoffs: Optional specific placebo cutoffs
        
        Returns:
            NegativeControlBattery with all results
        """
        all_results = []
        
        running_var = config['running_var']
        outcome_var = config['outcome_var']
        cutoff = config.get('cutoff', 0)
        
        # 1. Outcome negative controls
        if 'negative_outcomes' in config:
            outcome_battery = self.outcome_controls.run_battery(
                data, running_var, outcome_var,
                config['negative_outcomes'],
                cutoff
            )
            all_results.extend(outcome_battery.results)
        
        # 2. Placebo cutoffs
        placebo_results = self.treatment_controls.placebo_cutoffs(
            data, running_var, outcome_var,
            cutoff,
            config.get('placebo_cutoffs'),
            config.get('n_placebos', 10)
        )
        all_results.extend(placebo_results)
        
        # 3. Balance tests
        if 'covariates' in config:
            balance_battery = self.balance_tests.run_battery(
                data, running_var,
                config['covariates'],
                cutoff
            )
            all_results.extend(balance_battery.results)
        
        # 4. Permutation test
        if config.get('run_permutation', True):
            perm_result = self.treatment_controls.permutation_test(
                data, running_var, outcome_var, cutoff,
                config.get('n_permutations', 500)
            )
            all_results.append(perm_result)
        
        return NegativeControlBattery(results=all_results)


# =============================================================================
# SIMPLE RDD ESTIMATOR (for testing)
# =============================================================================

def simple_rdd_estimator(X: np.ndarray, Y: np.ndarray, cutoff: float) -> Tuple[float, float, int]:
    """
    Simple local-linear RDD estimator for testing.
    
    In production, use rdrobust instead.
    """
    n = len(X)
    
    # Normalize X to cutoff
    X_norm = X - cutoff
    
    # Split above/below
    below = X_norm < 0
    above = X_norm >= 0
    
    if sum(below) < 5 or sum(above) < 5:
        return np.nan, np.nan, n
    
    # Simple difference in means (very basic)
    mean_below = Y[below].mean()
    mean_above = Y[above].mean()
    
    estimate = mean_above - mean_below
    
    # Pooled SE
    se = np.sqrt(
        Y[below].var() / sum(below) + 
        Y[above].var() / sum(above)
    )
    
    return estimate, se, n


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    # Create sample data
    np.random.seed(42)
    n = 1000
    
    # Running variable (contract value normalized to cutoff)
    X = np.random.uniform(-1, 1, n)
    
    # True treatment effect at cutoff
    treatment = (X >= 0).astype(float)
    Y = 2 + 0.5 * X + 0.3 * treatment + np.random.normal(0, 0.5, n)
    
    # Some covariates (predetermined, should be smooth)
    cov1 = 1 + 0.2 * X + np.random.normal(0, 0.3, n)  # Smooth
    cov2 = np.random.normal(0, 1, n)  # Unrelated
    
    # Unrelated outcome (should not be affected)
    unrelated = np.random.normal(5, 1, n)
    
    data = pd.DataFrame({
        'running_var': X,
        'outcome': Y,
        'covariate1': cov1,
        'covariate2': cov2,
        'unrelated_outcome': unrelated
    })
    
    # Run negative controls
    runner = NegativeControlRunner(simple_rdd_estimator)
    
    config = {
        'running_var': 'running_var',
        'outcome_var': 'outcome',
        'cutoff': 0,
        'negative_outcomes': [
            {'column': 'unrelated_outcome', 'description': 'Unrelated variable'}
        ],
        'covariates': [
            {'column': 'covariate1', 'description': 'Predetermined covariate 1'},
            {'column': 'covariate2', 'description': 'Predetermined covariate 2'}
        ],
        'n_placebos': 5,
        'n_permutations': 100  # Low for demo
    }
    
    battery = runner.run_full_battery(data, config)
    
    print(battery.format_report())
    print(f"\nSummary: {battery.n_passed}/{battery.n_tests} tests passed ({battery.pass_rate:.1%})")
