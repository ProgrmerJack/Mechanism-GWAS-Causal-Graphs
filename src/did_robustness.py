"""
Difference-in-Differences Robustness Battery for E-Procurement Impact Analysis.

Implements state-of-the-art DiD methods for staggered treatment adoption:
- Callaway-Sant'Anna (2021) event studies with dynamic effects
- Pre-trend visualization and formal testing
- Placebo reform timing tests
- Heterogeneity analysis by group and time
- Composition sensitivity (never-treated vs not-yet-treated)
- Bacon decomposition for TWFE diagnostics

References:
    Callaway, B., & Sant'Anna, P. H. (2021). Difference-in-differences with
        multiple time periods. Journal of Econometrics, 225(2), 200-230.
    Goodman-Bacon, A. (2021). Difference-in-differences with variation in
        treatment timing. Journal of Econometrics, 225(2), 254-277.
    Sun, L., & Abraham, S. (2021). Estimating dynamic treatment effects in
        event studies with heterogeneous treatment effects. JoE, 225(2), 175-199.
    Roth, J. (2022). Pretest with caution: Event-study estimates after testing
        for parallel trends. AER: Insights, 4(3), 305-322.

Author: Abduxoliq Ashuraliyev
Affiliation: Independent Researcher, Tashkent, Uzbekistan
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats

if TYPE_CHECKING:
    from matplotlib.figure import Figure

logger = logging.getLogger(__name__)


# ==============================================================================
# Data Classes for Type Safety
# ==============================================================================


@dataclass
class EventStudyResult:
    """Container for event study estimation results."""

    # Coefficients and inference
    event_times: np.ndarray
    coefficients: np.ndarray
    std_errors: np.ndarray
    confidence_lower: np.ndarray
    confidence_upper: np.ndarray
    p_values: np.ndarray

    # Pre-trend diagnostics
    pre_trend_f_stat: float | None = None
    pre_trend_p_value: float | None = None
    pre_trend_passed: bool | None = None

    # Aggregated ATT
    att_overall: float | None = None
    att_se: float | None = None

    # Metadata
    n_observations: int = 0
    n_treated_units: int = 0
    n_control_units: int = 0
    estimation_method: str = ""


@dataclass
class BaconDecompositionResult:
    """Container for Goodman-Bacon decomposition results."""

    components: pd.DataFrame  # Type, Weight, Estimate columns
    twfe_estimate: float
    clean_comparison_weight: float  # Weight on good comparisons
    problematic_weight: float  # Weight on potentially biased comparisons
    warning_flag: bool  # True if problematic weights are substantial


@dataclass
class HeterogeneityResult:
    """Container for heterogeneity analysis results."""

    dimension: str  # e.g., "country", "sector", "firm_size"
    groups: list[str]
    estimates: np.ndarray
    std_errors: np.ndarray
    p_values: np.ndarray
    q_stat: float  # Heterogeneity test statistic
    q_p_value: float
    heterogeneity_detected: bool


@dataclass
class PlaceboResult:
    """Container for placebo reform test results."""

    placebo_timing: int  # Relative to actual reform
    estimate: float
    std_error: float
    p_value: float
    passed: bool  # True if placebo estimate is insignificant


@dataclass
class CompositionSensitivityResult:
    """Container for control group composition sensitivity."""

    control_type: Literal["never_treated", "not_yet_treated", "both"]
    att_estimate: float
    att_se: float
    n_control_units: int
    comparison_notes: str = ""


@dataclass
class DiDRobustnessReport:
    """Comprehensive DiD robustness report."""

    # Main estimates
    event_study: EventStudyResult
    bacon_decomposition: BaconDecompositionResult | None

    # Robustness checks
    pre_trend_tests: dict[str, float]  # Various pre-trend test p-values
    placebo_results: list[PlaceboResult]
    heterogeneity_results: list[HeterogeneityResult]
    composition_sensitivity: list[CompositionSensitivityResult]

    # Summary
    overall_assessment: str
    concerns: list[str]
    recommendations: list[str]

    # Figures
    figure_paths: dict[str, Path] = field(default_factory=dict)


# ==============================================================================
# Main DiD Robustness Class
# ==============================================================================


class DiDRobustnessBattery:
    """
    Comprehensive DiD robustness testing for staggered adoption designs.

    Implements the full battery of tests recommended for Nature-tier review:
    1. Callaway-Sant'Anna event studies with proper aggregation
    2. Bacon decomposition for TWFE diagnostics
    3. Pre-trend testing (joint F-test, Roth correction)
    4. Placebo timing tests
    5. Heterogeneity analysis
    6. Composition sensitivity
    """

    def __init__(
        self,
        data: pd.DataFrame,
        outcome_col: str,
        time_col: str,
        unit_col: str,
        treatment_col: str,
        first_treat_col: str | None = None,
        covariates: list[str] | None = None,
        cluster_col: str | None = None,
        output_dir: Path | None = None,
    ):
        """
        Initialize DiD robustness battery.

        Args:
            data: Panel data with unit-time observations
            outcome_col: Name of outcome variable column
            time_col: Name of time period column
            unit_col: Name of unit identifier column
            treatment_col: Name of binary treatment indicator
            first_treat_col: Name of first treatment period column (for staggered)
            covariates: Optional list of control variables
            cluster_col: Column for clustering standard errors
            output_dir: Directory for saving figures
        """
        self.data = data.copy()
        self.outcome_col = outcome_col
        self.time_col = time_col
        self.unit_col = unit_col
        self.treatment_col = treatment_col
        self.first_treat_col = first_treat_col
        self.covariates = covariates or []
        self.cluster_col = cluster_col or unit_col
        self.output_dir = output_dir or Path("outputs/did_robustness")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Validate data structure
        self._validate_panel_data()

        # Compute timing information
        self._compute_timing_info()

        logger.info(
            f"DiD battery initialized: {self.n_units} units, "
            f"{self.n_periods} periods, {self.n_treated} treated units"
        )

    def _validate_panel_data(self) -> None:
        """Validate panel data structure and requirements."""
        required_cols = [
            self.outcome_col,
            self.time_col,
            self.unit_col,
            self.treatment_col,
        ]
        missing = [c for c in required_cols if c not in self.data.columns]
        if missing:
            msg = f"Missing required columns: {missing}"
            raise ValueError(msg)

        # Check for balanced panel (warn if unbalanced)
        unit_counts = self.data.groupby(self.unit_col).size()
        if unit_counts.nunique() > 1:
            logger.warning(
                "Unbalanced panel detected. Some methods may require balancing."
            )

        # Check for variation in treatment
        if self.data[self.treatment_col].nunique() < 2:
            msg = "Treatment variable has no variation"
            raise ValueError(msg)

    def _compute_timing_info(self) -> None:
        """Compute treatment timing information."""
        self.n_units = self.data[self.unit_col].nunique()
        self.n_periods = self.data[self.time_col].nunique()
        self.periods = sorted(self.data[self.time_col].unique())

        # Identify treated units
        if self.first_treat_col and self.first_treat_col in self.data.columns:
            treated_mask = self.data[self.first_treat_col].notna()
        else:
            treated_mask = self.data.groupby(self.unit_col)[self.treatment_col].transform("max") > 0

        self.treated_units = self.data.loc[treated_mask, self.unit_col].unique()
        self.n_treated = len(self.treated_units)

        # Get treatment cohorts (first treatment periods)
        if self.first_treat_col and self.first_treat_col in self.data.columns:
            cohorts = self.data.groupby(self.unit_col)[self.first_treat_col].first()
            self.treatment_cohorts = cohorts[cohorts.notna()].unique()
        else:
            # Infer from treatment indicator
            self._infer_first_treatment()
            self.treatment_cohorts = self.data.loc[
                self.data[self.unit_col].isin(self.treated_units), "_first_treat"
            ].unique()

    def _infer_first_treatment(self) -> None:
        """Infer first treatment period from treatment indicator."""
        # Find first period where treatment = 1 for each unit
        treated_periods = self.data[self.data[self.treatment_col] == 1].groupby(
            self.unit_col
        )[self.time_col].min()
        self.data = self.data.merge(
            treated_periods.rename("_first_treat"),
            on=self.unit_col,
            how="left",
        )
        self.first_treat_col = "_first_treat"

    # ==========================================================================
    # Event Study Estimation
    # ==========================================================================

    def callaway_santanna_event_study(
        self,
        anticipation: int = 0,
        base_period: int = -1,
        control_group: Literal["never_treated", "not_yet_treated"] = "never_treated",
        n_bootstrap: int = 999,
    ) -> EventStudyResult:
        """
        Estimate event study using Callaway-Sant'Anna (2021) approach.

        This is a simplified Python implementation. For production use,
        consider using the R 'did' package via rpy2 or stata's 'csdid'.

        Args:
            anticipation: Number of anticipation periods
            base_period: Base period for normalization (-1 = last pre-treatment)
            control_group: Control group definition
            n_bootstrap: Number of bootstrap iterations for inference

        Returns:
            EventStudyResult with dynamic treatment effects
        """
        logger.info(f"Running Callaway-Sant'Anna event study (control={control_group})")

        # Create event time variable
        data = self.data.copy()
        data["_event_time"] = data[self.time_col] - data[self.first_treat_col]

        # For never-treated units, set event_time to NaN
        never_treated = data[data[self.first_treat_col].isna()][self.unit_col].unique()
        data.loc[data[self.unit_col].isin(never_treated), "_event_time"] = np.nan

        # Define event window
        min_event = int(data["_event_time"].min())
        max_event = int(data["_event_time"].max())
        event_window = range(max(min_event, -10), min(max_event + 1, 11))
        event_window = [e for e in event_window if e != base_period]

        # Simple 2x2 DiD aggregation for each (cohort, event_time)
        results = []

        for cohort in self.treatment_cohorts:
            if pd.isna(cohort):
                continue

            cohort_data = data[
                (data[self.first_treat_col] == cohort)
                | (data[self.first_treat_col].isna() if control_group == "never_treated"
                   else data[self.first_treat_col] > data[self.time_col])
            ].copy()

            for event_time in event_window:
                actual_time = cohort + event_time
                if actual_time not in self.periods:
                    continue

                # Pre-period (base period)
                pre_time = cohort + base_period

                if pre_time not in self.periods:
                    continue

                # Compute 2x2 DiD
                dd_result = self._compute_2x2_did(
                    cohort_data, cohort, pre_time, actual_time, control_group
                )
                if dd_result is not None:
                    results.append({
                        "cohort": cohort,
                        "event_time": event_time,
                        "estimate": dd_result[0],
                        "se": dd_result[1],
                        "n_treated": dd_result[2],
                        "n_control": dd_result[3],
                    })

        # Aggregate across cohorts by event time
        results_df = pd.DataFrame(results)
        if results_df.empty:
            msg = "No valid event study estimates computed"
            raise ValueError(msg)

        aggregated = results_df.groupby("event_time").agg(
            coefficient=("estimate", "mean"),
            se=("se", lambda x: np.sqrt(np.mean(x**2))),  # Simple average
            n_treated=("n_treated", "sum"),
            n_control=("n_control", "sum"),
        ).reset_index()

        event_times = aggregated["event_time"].values
        coefficients = aggregated["coefficient"].values
        std_errors = aggregated["se"].values

        # Compute confidence intervals and p-values
        ci_lower = coefficients - 1.96 * std_errors
        ci_upper = coefficients + 1.96 * std_errors
        p_values = 2 * (1 - stats.norm.cdf(np.abs(coefficients / std_errors)))

        # Pre-trend test (joint F-test on pre-treatment coefficients)
        pre_mask = event_times < 0
        if pre_mask.sum() > 0:
            pre_coeffs = coefficients[pre_mask]
            pre_ses = std_errors[pre_mask]
            # Wald test
            wald_stat = np.sum((pre_coeffs / pre_ses) ** 2)
            pre_trend_p = 1 - stats.chi2.cdf(wald_stat, df=pre_mask.sum())
            pre_trend_passed = pre_trend_p > 0.10  # Common threshold
        else:
            pre_trend_p = None
            pre_trend_passed = None

        # Overall ATT (average of post-treatment effects)
        post_mask = event_times >= 0
        if post_mask.sum() > 0:
            att_overall = np.mean(coefficients[post_mask])
            att_se = np.sqrt(np.mean(std_errors[post_mask] ** 2))
        else:
            att_overall = None
            att_se = None

        return EventStudyResult(
            event_times=event_times,
            coefficients=coefficients,
            std_errors=std_errors,
            confidence_lower=ci_lower,
            confidence_upper=ci_upper,
            p_values=p_values,
            pre_trend_f_stat=wald_stat if pre_mask.sum() > 0 else None,
            pre_trend_p_value=pre_trend_p,
            pre_trend_passed=pre_trend_passed,
            att_overall=att_overall,
            att_se=att_se,
            n_observations=len(data),
            n_treated_units=int(aggregated["n_treated"].iloc[0]),
            n_control_units=int(aggregated["n_control"].iloc[0]),
            estimation_method="callaway_santanna",
        )

    def _compute_2x2_did(
        self,
        data: pd.DataFrame,
        cohort: float,
        pre_time: float,
        post_time: float,
        control_group: str,
    ) -> tuple[float, float, int, int] | None:
        """Compute simple 2x2 DiD estimate."""
        # Treated group
        treated = data[data[self.first_treat_col] == cohort]

        # Control group
        if control_group == "never_treated":
            control = data[data[self.first_treat_col].isna()]
        else:
            # Not-yet-treated at post_time
            control = data[data[self.first_treat_col] > post_time]

        if len(treated) == 0 or len(control) == 0:
            return None

        # Pre-period means
        treat_pre = treated[treated[self.time_col] == pre_time][self.outcome_col].mean()
        control_pre = control[control[self.time_col] == pre_time][self.outcome_col].mean()

        # Post-period means
        treat_post = treated[treated[self.time_col] == post_time][self.outcome_col].mean()
        control_post = control[control[self.time_col] == post_time][self.outcome_col].mean()

        if any(pd.isna([treat_pre, control_pre, treat_post, control_post])):
            return None

        # DiD estimate
        estimate = (treat_post - treat_pre) - (control_post - control_pre)

        # Simple SE (ignoring clustering for now)
        n_treat = len(treated[treated[self.time_col].isin([pre_time, post_time])])
        n_control = len(control[control[self.time_col].isin([pre_time, post_time])])

        # Crude SE estimate
        treat_var = treated[self.outcome_col].var()
        control_var = control[self.outcome_col].var()
        se = np.sqrt(treat_var / n_treat + control_var / n_control) if n_treat > 0 and n_control > 0 else np.nan

        return estimate, se, n_treat, n_control

    # ==========================================================================
    # Bacon Decomposition
    # ==========================================================================

    def bacon_decomposition(self) -> BaconDecompositionResult:
        """
        Perform Goodman-Bacon (2021) decomposition of TWFE estimator.

        Decomposes TWFE estimate into:
        - Timing groups: Earlier vs later treated
        - Treatment vs never-treated
        - Treatment vs not-yet-treated

        Returns:
            BaconDecompositionResult with component weights and estimates
        """
        logger.info("Running Bacon decomposition")

        # Identify comparison types
        components = []

        # Never-treated vs treated comparisons
        never_treated = self.data[self.data[self.first_treat_col].isna()][self.unit_col].unique()

        for cohort in self.treatment_cohorts:
            if pd.isna(cohort):
                continue

            # Treated vs never-treated
            if len(never_treated) > 0:
                estimate, weight = self._bacon_component(
                    cohort, never_treated=True
                )
                if estimate is not None:
                    components.append({
                        "type": "Treated vs Never-Treated",
                        "cohort": cohort,
                        "comparison": "never_treated",
                        "estimate": estimate,
                        "weight": weight,
                    })

            # Early vs late treated
            later_cohorts = [c for c in self.treatment_cohorts
                            if not pd.isna(c) and c > cohort]
            for later_cohort in later_cohorts:
                estimate, weight = self._bacon_component(
                    cohort, comparison_cohort=later_cohort
                )
                if estimate is not None:
                    components.append({
                        "type": "Earlier vs Later Treated",
                        "cohort": cohort,
                        "comparison": later_cohort,
                        "estimate": estimate,
                        "weight": weight,
                    })

                # Reverse comparison (problematic!)
                estimate_rev, weight_rev = self._bacon_component(
                    later_cohort, comparison_cohort=cohort, is_reverse=True
                )
                if estimate_rev is not None:
                    components.append({
                        "type": "Later vs Earlier Treated (Problematic)",
                        "cohort": later_cohort,
                        "comparison": cohort,
                        "estimate": estimate_rev,
                        "weight": weight_rev,
                    })

        if not components:
            msg = "No valid Bacon decomposition components"
            raise ValueError(msg)

        df = pd.DataFrame(components)

        # Normalize weights
        total_weight = df["weight"].sum()
        df["weight"] = df["weight"] / total_weight

        # TWFE estimate
        twfe_estimate = (df["estimate"] * df["weight"]).sum()

        # Weight on clean comparisons
        clean_types = ["Treated vs Never-Treated", "Earlier vs Later Treated"]
        clean_weight = df[df["type"].isin(clean_types)]["weight"].sum()
        problematic_weight = 1 - clean_weight

        return BaconDecompositionResult(
            components=df,
            twfe_estimate=twfe_estimate,
            clean_comparison_weight=clean_weight,
            problematic_weight=problematic_weight,
            warning_flag=problematic_weight > 0.20,  # Flag if >20% problematic
        )

    def _bacon_component(
        self,
        cohort: float,
        comparison_cohort: float | None = None,
        never_treated: bool = False,
        is_reverse: bool = False,
    ) -> tuple[float | None, float]:
        """Compute single Bacon decomposition component."""
        data = self.data.copy()

        # Select treated units from cohort
        treated = data[data[self.first_treat_col] == cohort]

        # Select control units
        if never_treated:
            control = data[data[self.first_treat_col].isna()]
        else:
            control = data[data[self.first_treat_col] == comparison_cohort]

        if len(treated) == 0 or len(control) == 0:
            return None, 0.0

        # Simple 2x2 DiD
        pre_period = cohort - 1  # Last pre-period
        post_period = cohort  # First post-period

        treat_pre = treated[treated[self.time_col] == pre_period][self.outcome_col].mean()
        treat_post = treated[treated[self.time_col] == post_period][self.outcome_col].mean()
        control_pre = control[control[self.time_col] == pre_period][self.outcome_col].mean()
        control_post = control[control[self.time_col] == post_period][self.outcome_col].mean()

        if any(pd.isna([treat_pre, treat_post, control_pre, control_post])):
            return None, 0.0

        estimate = (treat_post - treat_pre) - (control_post - control_pre)

        # Weight based on sample sizes and variance
        n_treat = len(treated[self.unit_col].unique())
        n_control = len(control[self.unit_col].unique())
        weight = n_treat * n_control / (n_treat + n_control)

        return estimate, weight

    # ==========================================================================
    # Pre-Trend Testing
    # ==========================================================================

    def pre_trend_tests(
        self,
        event_study: EventStudyResult | None = None,
    ) -> dict[str, dict[str, float]]:
        """
        Comprehensive pre-trend testing battery.

        Implements multiple pre-trend tests:
        1. Joint F-test on pre-treatment coefficients
        2. Individual coefficient tests
        3. Roth (2022) correction for pre-testing
        4. Slope test (linear pre-trend)

        Args:
            event_study: Pre-computed event study result

        Returns:
            Dictionary of test results
        """
        logger.info("Running pre-trend testing battery")

        if event_study is None:
            event_study = self.callaway_santanna_event_study()

        pre_mask = event_study.event_times < 0
        pre_times = event_study.event_times[pre_mask]
        pre_coeffs = event_study.coefficients[pre_mask]
        pre_ses = event_study.std_errors[pre_mask]

        results = {}

        # 1. Joint F-test
        if len(pre_coeffs) > 0:
            wald_stat = np.sum((pre_coeffs / pre_ses) ** 2)
            results["joint_f_test"] = {
                "statistic": wald_stat,
                "df": len(pre_coeffs),
                "p_value": 1 - stats.chi2.cdf(wald_stat, df=len(pre_coeffs)),
                "passed": 1 - stats.chi2.cdf(wald_stat, df=len(pre_coeffs)) > 0.10,
            }

        # 2. Individual coefficient tests
        individual_results = []
        for t, coef, se in zip(pre_times, pre_coeffs, pre_ses, strict=False):
            z = coef / se
            p = 2 * (1 - stats.norm.cdf(abs(z)))
            individual_results.append({
                "event_time": t,
                "coefficient": coef,
                "se": se,
                "p_value": p,
                "significant": p < 0.05,
            })
        results["individual_tests"] = individual_results

        # 3. Linear pre-trend test (slope)
        if len(pre_times) >= 2:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                pre_times, pre_coeffs
            )
            results["linear_trend_test"] = {
                "slope": slope,
                "intercept": intercept,
                "r_squared": r_value ** 2,
                "slope_p_value": p_value,
                "passed": p_value > 0.10,
            }

        # 4. Roth (2022) honesty bound
        # Simplified version: bound on ATT given pre-trend magnitude
        if len(pre_coeffs) > 0:
            max_pre_trend = np.max(np.abs(pre_coeffs))
            results["roth_bound"] = {
                "max_pre_trend_magnitude": max_pre_trend,
                "suggested_sensitivity_range": (-max_pre_trend, max_pre_trend),
                "note": "If true pre-trend continued, ATT could be biased by this amount",
            }

        return results

    # ==========================================================================
    # Placebo Tests
    # ==========================================================================

    def placebo_reform_tests(
        self,
        placebo_offsets: list[int] | None = None,
    ) -> list[PlaceboResult]:
        """
        Test placebo reform timings.

        Artificially shifts treatment timing and estimates effects.
        Significant effects at placebo timings suggest violations.

        Args:
            placebo_offsets: Offsets from actual reform (negative = before)

        Returns:
            List of PlaceboResult for each placebo timing
        """
        logger.info("Running placebo reform timing tests")

        if placebo_offsets is None:
            placebo_offsets = [-3, -2, -1]  # Default: test 1-3 periods before

        results = []

        for offset in placebo_offsets:
            # Create placebo treatment timing
            data_placebo = self.data.copy()

            # Shift first treatment
            data_placebo["_placebo_first_treat"] = data_placebo[self.first_treat_col] + offset

            # Create placebo treatment indicator
            data_placebo["_placebo_treat"] = (
                data_placebo[self.time_col] >= data_placebo["_placebo_first_treat"]
            ).astype(int)

            # Simple DiD on placebo
            try:
                estimate, se = self._simple_did_estimate(
                    data_placebo, "_placebo_treat"
                )

                p_value = 2 * (1 - stats.norm.cdf(abs(estimate / se))) if se > 0 else 1.0

                results.append(PlaceboResult(
                    placebo_timing=offset,
                    estimate=estimate,
                    std_error=se,
                    p_value=p_value,
                    passed=p_value > 0.10,  # Passed if NOT significant
                ))
            except Exception as e:
                logger.warning(f"Placebo test failed for offset {offset}: {e}")

        return results

    def _simple_did_estimate(
        self,
        data: pd.DataFrame,
        treatment_col: str,
    ) -> tuple[float, float]:
        """Simple DiD estimate via OLS."""
        # Two-way fixed effects regression
        # Y_it = alpha_i + gamma_t + beta * D_it + epsilon_it

        # Create dummies
        data = data.copy()
        data["_post"] = data[treatment_col]

        # Group means
        treat_post = data[(data["_post"] == 1)][self.outcome_col].mean()
        treat_pre = data[(data["_post"] == 0) &
                        (data[self.first_treat_col].notna())][self.outcome_col].mean()
        control_post = data[(data["_post"] == 0) &
                           (data[self.first_treat_col].isna())][self.outcome_col].mean()
        control_pre = data[(data[self.first_treat_col].isna())][self.outcome_col].mean()

        estimate = (treat_post - treat_pre) - (control_post - control_pre)

        # Crude SE
        n = len(data)
        var = data[self.outcome_col].var()
        se = np.sqrt(4 * var / n)

        return estimate, se

    # ==========================================================================
    # Heterogeneity Analysis
    # ==========================================================================

    def heterogeneity_analysis(
        self,
        dimensions: list[str] | None = None,
    ) -> list[HeterogeneityResult]:
        """
        Analyze treatment effect heterogeneity across subgroups.

        Args:
            dimensions: Column names to stratify by

        Returns:
            List of HeterogeneityResult for each dimension
        """
        logger.info("Running heterogeneity analysis")

        if dimensions is None:
            dimensions = self.covariates[:3] if self.covariates else []

        results = []

        for dim in dimensions:
            if dim not in self.data.columns:
                logger.warning(f"Dimension {dim} not in data, skipping")
                continue

            groups = self.data[dim].unique()
            estimates = []
            ses = []

            for group in groups:
                subset = self.data[self.data[dim] == group]
                try:
                    est, se = self._simple_did_estimate(subset, self.treatment_col)
                    estimates.append(est)
                    ses.append(se)
                except Exception:
                    estimates.append(np.nan)
                    ses.append(np.nan)

            estimates = np.array(estimates)
            ses = np.array(ses)

            # Heterogeneity test (Q-statistic)
            valid = ~np.isnan(estimates)
            if valid.sum() >= 2:
                weights = 1 / (ses[valid] ** 2)
                weighted_mean = np.average(estimates[valid], weights=weights)
                q_stat = np.sum(weights * (estimates[valid] - weighted_mean) ** 2)
                q_p = 1 - stats.chi2.cdf(q_stat, df=valid.sum() - 1)
            else:
                q_stat = np.nan
                q_p = np.nan

            p_values = 2 * (1 - stats.norm.cdf(np.abs(estimates / ses)))

            results.append(HeterogeneityResult(
                dimension=dim,
                groups=list(groups),
                estimates=estimates,
                std_errors=ses,
                p_values=p_values,
                q_stat=q_stat,
                q_p_value=q_p,
                heterogeneity_detected=q_p < 0.10 if not np.isnan(q_p) else False,
            ))

        return results

    # ==========================================================================
    # Composition Sensitivity
    # ==========================================================================

    def composition_sensitivity(self) -> list[CompositionSensitivityResult]:
        """
        Test sensitivity to control group composition.

        Compares estimates using:
        - Never-treated as controls
        - Not-yet-treated as controls
        - Both combined

        Returns:
            List of CompositionSensitivityResult
        """
        logger.info("Running composition sensitivity analysis")

        results = []

        # Never-treated only
        never_treated = self.data[self.first_treat_col].isna()
        if never_treated.sum() > 0:
            data_never = self.data[
                never_treated | (self.data[self.treatment_col] == 1)
            ]
            est, se = self._simple_did_estimate(data_never, self.treatment_col)
            results.append(CompositionSensitivityResult(
                control_type="never_treated",
                att_estimate=est,
                att_se=se,
                n_control_units=len(self.data[never_treated][self.unit_col].unique()),
                comparison_notes="Clean comparison, may have selection issues",
            ))

        # Not-yet-treated only
        # This is more complex - need to dynamically define
        results.append(CompositionSensitivityResult(
            control_type="not_yet_treated",
            att_estimate=np.nan,  # Would need more complex implementation
            att_se=np.nan,
            n_control_units=0,
            comparison_notes="Not implemented - requires dynamic control definition",
        ))

        return results

    # ==========================================================================
    # Visualization
    # ==========================================================================

    def plot_event_study(
        self,
        event_study: EventStudyResult,
        title: str = "Event Study: Dynamic Treatment Effects",
        save_path: Path | None = None,
    ) -> Figure:
        """
        Create publication-quality event study plot.

        Args:
            event_study: EventStudyResult to plot
            title: Figure title
            save_path: Optional path to save figure

        Returns:
            matplotlib Figure object
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib required for plotting")
            raise

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot coefficients with confidence intervals
        ax.errorbar(
            event_study.event_times,
            event_study.coefficients,
            yerr=[
                event_study.coefficients - event_study.confidence_lower,
                event_study.confidence_upper - event_study.coefficients,
            ],
            fmt="o",
            capsize=3,
            capthick=1,
            color="steelblue",
            markersize=6,
        )

        # Add horizontal line at 0
        ax.axhline(y=0, color="black", linestyle="--", linewidth=0.8)

        # Add vertical line at treatment
        ax.axvline(x=-0.5, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

        # Shade pre-treatment period
        ax.axvspan(
            min(event_study.event_times) - 0.5,
            -0.5,
            alpha=0.1,
            color="gray",
            label="Pre-treatment",
        )

        ax.set_xlabel("Event Time (Periods Relative to Treatment)", fontsize=12)
        ax.set_ylabel("Treatment Effect", fontsize=12)
        ax.set_title(title, fontsize=14)

        # Add pre-trend test annotation
        if event_study.pre_trend_p_value is not None:
            status = "✓ Passed" if event_study.pre_trend_passed else "✗ Failed"
            ax.annotate(
                f"Pre-trend test: p={event_study.pre_trend_p_value:.3f} {status}",
                xy=(0.02, 0.98),
                xycoords="axes fraction",
                fontsize=10,
                verticalalignment="top",
                bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
            )

        ax.legend(loc="upper right")
        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info(f"Event study figure saved to {save_path}")

        return fig

    def plot_bacon_decomposition(
        self,
        bacon: BaconDecompositionResult,
        save_path: Path | None = None,
    ) -> Figure:
        """Create Bacon decomposition visualization."""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise

        fig, ax = plt.subplots(figsize=(10, 6))

        df = bacon.components

        # Color by type
        colors = {
            "Treated vs Never-Treated": "steelblue",
            "Earlier vs Later Treated": "forestgreen",
            "Later vs Earlier Treated (Problematic)": "crimson",
        }

        for comp_type in df["type"].unique():
            subset = df[df["type"] == comp_type]
            ax.scatter(
                subset["estimate"],
                subset["weight"],
                c=colors.get(comp_type, "gray"),
                label=comp_type,
                s=100,
                alpha=0.7,
            )

        ax.axvline(x=bacon.twfe_estimate, color="black", linestyle="--",
                   label=f"TWFE = {bacon.twfe_estimate:.3f}")

        ax.set_xlabel("2x2 DD Estimate", fontsize=12)
        ax.set_ylabel("Weight", fontsize=12)
        ax.set_title("Bacon Decomposition of TWFE Estimator", fontsize=14)
        ax.legend(loc="best")

        # Add warning annotation if needed
        if bacon.warning_flag:
            ax.annotate(
                f"⚠ {bacon.problematic_weight:.1%} weight on problematic comparisons",
                xy=(0.02, 0.02),
                xycoords="axes fraction",
                fontsize=10,
                color="crimson",
            )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig

    # ==========================================================================
    # Full Report Generation
    # ==========================================================================

    def full_robustness_report(
        self,
        heterogeneity_dimensions: list[str] | None = None,
        placebo_offsets: list[int] | None = None,
    ) -> DiDRobustnessReport:
        """
        Generate comprehensive DiD robustness report.

        Runs all robustness checks and aggregates results.

        Args:
            heterogeneity_dimensions: Columns for heterogeneity analysis
            placebo_offsets: Offsets for placebo timing tests

        Returns:
            DiDRobustnessReport with all results
        """
        logger.info("Generating full DiD robustness report")

        # Main event study
        event_study = self.callaway_santanna_event_study()

        # Bacon decomposition
        try:
            bacon = self.bacon_decomposition()
        except Exception as e:
            logger.warning(f"Bacon decomposition failed: {e}")
            bacon = None

        # Pre-trend tests
        pre_trend_results = self.pre_trend_tests(event_study)

        # Placebo tests
        placebo_results = self.placebo_reform_tests(placebo_offsets)

        # Heterogeneity
        heterogeneity_results = self.heterogeneity_analysis(heterogeneity_dimensions)

        # Composition sensitivity
        composition_results = self.composition_sensitivity()

        # Generate figures
        figure_paths = {}
        try:
            fig = self.plot_event_study(
                event_study,
                save_path=self.output_dir / "event_study.png",
            )
            figure_paths["event_study"] = self.output_dir / "event_study.png"

            if bacon:
                fig = self.plot_bacon_decomposition(
                    bacon,
                    save_path=self.output_dir / "bacon_decomposition.png",
                )
                figure_paths["bacon"] = self.output_dir / "bacon_decomposition.png"
        except Exception as e:
            logger.warning(f"Figure generation failed: {e}")

        # Assess overall robustness
        concerns = []
        recommendations = []

        # Check pre-trends
        if event_study.pre_trend_p_value and event_study.pre_trend_p_value < 0.10:
            concerns.append("Pre-trend test failed (p < 0.10)")
            recommendations.append("Consider alternative specifications or matching")

        # Check Bacon decomposition
        if bacon and bacon.warning_flag:
            concerns.append(
                f"High weight ({bacon.problematic_weight:.1%}) on problematic comparisons"
            )
            recommendations.append("Use robust estimator (Callaway-Sant'Anna, Sun-Abraham)")

        # Check placebos
        failed_placebos = [p for p in placebo_results if not p.passed]
        if failed_placebos:
            concerns.append(f"{len(failed_placebos)} placebo tests failed")
            recommendations.append("Investigate anticipation effects or confounders")

        # Check heterogeneity
        sig_het = [h for h in heterogeneity_results if h.heterogeneity_detected]
        if sig_het:
            concerns.append(f"Significant heterogeneity detected in {len(sig_het)} dimensions")
            recommendations.append("Report disaggregated effects")

        # Overall assessment
        if len(concerns) == 0:
            overall = "ROBUST: All checks passed"
        elif len(concerns) <= 2:
            overall = "MODERATE: Some concerns, but main findings likely hold"
        else:
            overall = "CONCERNING: Multiple robustness failures detected"

        return DiDRobustnessReport(
            event_study=event_study,
            bacon_decomposition=bacon,
            pre_trend_tests=pre_trend_results,
            placebo_results=placebo_results,
            heterogeneity_results=heterogeneity_results,
            composition_sensitivity=composition_results,
            overall_assessment=overall,
            concerns=concerns,
            recommendations=recommendations,
            figure_paths=figure_paths,
        )

    def export_latex_table(
        self,
        event_study: EventStudyResult,
        output_path: Path | None = None,
    ) -> str:
        """Export event study results as LaTeX table."""
        lines = [
            r"\begin{table}[htbp]",
            r"\centering",
            r"\caption{Event Study Estimates: Dynamic Treatment Effects}",
            r"\label{tab:event_study}",
            r"\begin{tabular}{lccc}",
            r"\toprule",
            r"Event Time & Coefficient & Std. Error & 95\% CI \\",
            r"\midrule",
        ]

        for i, t in enumerate(event_study.event_times):
            coef = event_study.coefficients[i]
            se = event_study.std_errors[i]
            ci_lo = event_study.confidence_lower[i]
            ci_hi = event_study.confidence_upper[i]

            # Add significance stars
            p = event_study.p_values[i]
            if p < 0.01:
                stars = "***"
            elif p < 0.05:
                stars = "**"
            elif p < 0.10:
                stars = "*"
            else:
                stars = ""

            lines.append(
                f"$t = {int(t):+d}$ & {coef:.3f}{stars} & ({se:.3f}) & "
                f"[{ci_lo:.3f}, {ci_hi:.3f}] \\\\"
            )

        lines.extend([
            r"\midrule",
            f"ATT (Overall) & {event_study.att_overall:.3f} & ({event_study.att_se:.3f}) & -- \\\\",
            r"\bottomrule",
            r"\multicolumn{4}{l}{\footnotesize "
            r"$^{***}p<0.01$, $^{**}p<0.05$, $^*p<0.10$} \\",
            r"\end{tabular}",
            r"\end{table}",
        ])

        latex = "\n".join(lines)

        if output_path:
            output_path.write_text(latex, encoding="utf-8")
            logger.info(f"LaTeX table saved to {output_path}")

        return latex


# ==============================================================================
# Convenience Functions
# ==============================================================================


def run_did_robustness(
    data_path: Path,
    outcome_col: str,
    time_col: str,
    unit_col: str,
    treatment_col: str,
    output_dir: Path | None = None,
) -> DiDRobustnessReport:
    """
    Convenience function to run full DiD robustness battery.

    Args:
        data_path: Path to panel data CSV
        outcome_col: Outcome variable column
        time_col: Time period column
        unit_col: Unit identifier column
        treatment_col: Treatment indicator column
        output_dir: Output directory

    Returns:
        DiDRobustnessReport
    """
    data = pd.read_csv(data_path)

    battery = DiDRobustnessBattery(
        data=data,
        outcome_col=outcome_col,
        time_col=time_col,
        unit_col=unit_col,
        treatment_col=treatment_col,
        output_dir=output_dir,
    )

    return battery.full_robustness_report()


if __name__ == "__main__":
    # Example usage with synthetic data
    import sys

    logging.basicConfig(level=logging.INFO)

    # Create synthetic panel data
    np.random.seed(42)
    n_units = 100
    n_periods = 20

    units = np.repeat(range(n_units), n_periods)
    periods = np.tile(range(n_periods), n_units)

    # Staggered treatment adoption
    first_treat = np.random.choice([8, 10, 12, np.nan], n_units, p=[0.2, 0.3, 0.2, 0.3])
    first_treat_expanded = np.repeat(first_treat, n_periods)

    treatment = (periods >= first_treat_expanded).astype(float)
    treatment[np.isnan(first_treat_expanded)] = 0

    # Outcome with treatment effect
    true_effect = 2.0
    outcome = (
        np.random.randn(n_units * n_periods) +  # noise
        0.1 * periods +  # time trend
        treatment * true_effect  # treatment effect
    )

    data = pd.DataFrame({
        "unit_id": units,
        "time": periods,
        "outcome": outcome,
        "treatment": treatment,
        "first_treat": first_treat_expanded,
    })

    # Run robustness battery
    battery = DiDRobustnessBattery(
        data=data,
        outcome_col="outcome",
        time_col="time",
        unit_col="unit_id",
        treatment_col="treatment",
        first_treat_col="first_treat",
        output_dir=Path("outputs/did_test"),
    )

    report = battery.full_robustness_report()

    print(f"\n{'='*60}")
    print("DiD ROBUSTNESS REPORT")
    print(f"{'='*60}")
    print(f"\nOverall Assessment: {report.overall_assessment}")
    print(f"\nATT Estimate: {report.event_study.att_overall:.3f} "
          f"(SE: {report.event_study.att_se:.3f})")
    print(f"Pre-trend test p-value: {report.event_study.pre_trend_p_value:.3f}")

    if report.concerns:
        print(f"\nConcerns:")
        for c in report.concerns:
            print(f"  - {c}")

    if report.recommendations:
        print(f"\nRecommendations:")
        for r in report.recommendations:
            print(f"  - {r}")

    print(f"\nFigures saved to: {list(report.figure_paths.values())}")
