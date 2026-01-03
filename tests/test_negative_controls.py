"""
Tests for Negative Controls Module.

Tests the falsification test framework for causal inference validation.
"""

import pytest
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pathlib import Path
from negative_controls import (
    NegativeControlAnalyzer,
    NegativeControlResult,
    TemporalPlaceboResult,
    GeographicPlaceboResult,
)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def sample_procurement_data():
    """Create sample procurement data for testing."""
    np.random.seed(42)
    n = 500
    
    # Create data with clear treatment effect on main outcome
    # but no effect on negative control outcomes
    treatment = np.concatenate([np.zeros(250), np.ones(250)])
    
    # Main outcome: clear treatment effect
    main_outcome = 10 + 5 * treatment + np.random.normal(0, 2, n)
    
    # Negative control outcome: no treatment effect (pre-determined)
    negative_control_outcome = 8 + np.random.normal(0, 2, n)
    
    # Another outcome that shouldn't be affected
    unrelated_outcome = 15 + np.random.normal(0, 3, n)
    
    # Contract values (running variable for RDD)
    contract_value = np.random.lognormal(10, 1, n)
    
    # Time periods
    base_date = datetime(2018, 1, 1)
    dates = [base_date + timedelta(days=int(d)) for d in np.random.uniform(0, 730, n)]
    
    return pd.DataFrame({
        'contract_id': range(n),
        'treatment': treatment,
        'main_outcome': main_outcome,
        'negative_control_outcome': negative_control_outcome,
        'unrelated_outcome': unrelated_outcome,
        'contract_value': contract_value,
        'date': dates,
        'region': np.random.choice(['A', 'B', 'C'], n),
        'sector': np.random.choice(['health', 'infrastructure', 'defense'], n),
    })


@pytest.fixture
def sample_did_data():
    """Create sample data for DiD negative control tests."""
    np.random.seed(123)
    
    # Panel data: 50 units, 10 time periods
    n_units = 50
    n_periods = 10
    
    data = []
    for unit in range(n_units):
        # Treatment timing varies
        treatment_period = np.random.choice([5, 6, 7, None], p=[0.3, 0.3, 0.2, 0.2])
        
        for period in range(n_periods):
            treated = treatment_period is not None and period >= treatment_period
            
            # Main outcome: treatment effect of 3
            main_y = 10 + 3 * treated + 0.5 * period + np.random.normal(0, 1)
            
            # Negative control: no treatment effect
            control_y = 8 + 0.3 * period + np.random.normal(0, 1)
            
            data.append({
                'unit_id': unit,
                'period': period,
                'treated': int(treated),
                'treatment_period': treatment_period,
                'main_outcome': main_y,
                'negative_control_outcome': control_y,
                'region': 'treated' if treatment_period is not None else 'control',
            })
    
    return pd.DataFrame(data)


@pytest.fixture
def analyzer(sample_procurement_data):
    """Create analyzer instance with sample data."""
    return NegativeControlAnalyzer(sample_procurement_data)


# ============================================================================
# Test NegativeControlAnalyzer Initialization
# ============================================================================

class TestNegativeControlAnalyzerInit:
    """Tests for NegativeControlAnalyzer initialization."""
    
    def test_init_with_dataframe(self, sample_procurement_data):
        """Test initialization with DataFrame."""
        analyzer = NegativeControlAnalyzer(sample_procurement_data)
        assert analyzer.data is not None
        assert len(analyzer.data) == len(sample_procurement_data)
    
    def test_init_without_data(self):
        """Test initialization without data."""
        analyzer = NegativeControlAnalyzer()
        assert analyzer.data is None
    
    def test_set_data(self, sample_procurement_data):
        """Test setting data after initialization."""
        analyzer = NegativeControlAnalyzer()
        analyzer.set_data(sample_procurement_data)
        assert analyzer.data is not None


# ============================================================================
# Test Outcome Negative Controls
# ============================================================================

class TestOutcomeNegativeControls:
    """Tests for outcome negative control analysis."""
    
    def test_outcome_control_returns_result(self, analyzer):
        """Test that outcome control returns proper result structure."""
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='pre_determined_outcome'
        )
        
        assert isinstance(result, NegativeControlResult)
        assert result.control_type == 'outcome'
        assert result.control_name == 'pre_determined_outcome'
        assert hasattr(result, 'estimate')
        assert hasattr(result, 'std_error')
        assert hasattr(result, 'p_value')
        assert hasattr(result, 'passed')
    
    def test_negative_control_should_pass(self, analyzer):
        """Test that true negative control passes (effect ≈ 0)."""
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='pre_determined_outcome'
        )
        
        # Should pass: estimate should be close to 0
        assert result.passed, f"Expected pass, got estimate={result.estimate:.3f}, p={result.p_value:.3f}"
        assert abs(result.estimate) < 1.0  # Effect should be small
    
    def test_positive_control_should_fail(self, analyzer):
        """Test that affected outcome fails negative control test."""
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='main_outcome',  # This HAS treatment effect
            control_name='affected_outcome'
        )
        
        # Should fail: there IS a treatment effect
        assert not result.passed, f"Expected fail, got estimate={result.estimate:.3f}"
        assert abs(result.estimate) > 2.0  # Effect should be substantial
    
    def test_multiple_outcome_controls(self, analyzer):
        """Test running multiple outcome negative controls."""
        outcomes = ['negative_control_outcome', 'unrelated_outcome']
        results = []
        
        for outcome in outcomes:
            result = analyzer.outcome_negative_controls(
                treatment_col='treatment',
                outcome_col=outcome,
                control_name=outcome
            )
            results.append(result)
        
        assert len(results) == 2
        assert all(isinstance(r, NegativeControlResult) for r in results)
    
    def test_with_covariates(self, analyzer):
        """Test outcome control with covariate adjustment."""
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='with_covariates',
            covariates=['contract_value']
        )
        
        assert isinstance(result, NegativeControlResult)
        assert result.passed


# ============================================================================
# Test Exposure Negative Controls
# ============================================================================

class TestExposureNegativeControls:
    """Tests for exposure negative control analysis."""
    
    def test_exposure_control_returns_result(self, analyzer, sample_procurement_data):
        """Test that exposure control returns proper result structure."""
        # Define exempt group (e.g., defense sector)
        sample_procurement_data['exempt'] = sample_procurement_data['sector'] == 'defense'
        analyzer.set_data(sample_procurement_data)
        
        result = analyzer.exposure_negative_controls(
            treatment_col='treatment',
            outcome_col='main_outcome',
            exempt_col='exempt',
            control_name='defense_exempt'
        )
        
        assert isinstance(result, NegativeControlResult)
        assert result.control_type == 'exposure'
        assert result.control_name == 'defense_exempt'
    
    def test_exempt_group_effect(self, sample_procurement_data):
        """Test that exempt group shows no treatment effect."""
        # Create data where exempt group truly has no effect
        sample_procurement_data = sample_procurement_data.copy()
        sample_procurement_data['exempt'] = sample_procurement_data['sector'] == 'defense'
        
        # Remove treatment effect for defense sector
        defense_mask = sample_procurement_data['sector'] == 'defense'
        sample_procurement_data.loc[defense_mask, 'main_outcome'] = (
            10 + np.random.normal(0, 2, defense_mask.sum())
        )
        
        analyzer = NegativeControlAnalyzer(sample_procurement_data)
        result = analyzer.exposure_negative_controls(
            treatment_col='treatment',
            outcome_col='main_outcome',
            exempt_col='exempt',
            control_name='defense_exempt'
        )
        
        # Effect in exempt group should be near zero
        assert abs(result.estimate) < 2.0


# ============================================================================
# Test Temporal Placebos
# ============================================================================

class TestTemporalPlacebos:
    """Tests for temporal placebo analysis."""
    
    def test_temporal_placebo_returns_result(self, sample_did_data):
        """Test that temporal placebo returns proper result structure."""
        analyzer = NegativeControlAnalyzer(sample_did_data)
        
        result = analyzer.temporal_placebos(
            unit_col='unit_id',
            time_col='period',
            outcome_col='main_outcome',
            treatment_col='treated',
            true_treatment_time=6,
            placebo_times=[2, 3, 4]
        )
        
        assert isinstance(result, TemporalPlaceboResult)
        assert hasattr(result, 'placebo_estimates')
        assert hasattr(result, 'joint_test_pvalue')
        assert hasattr(result, 'passed')
    
    def test_placebo_times_before_treatment(self, sample_did_data):
        """Test that placebo effects at fake times are near zero."""
        analyzer = NegativeControlAnalyzer(sample_did_data)
        
        result = analyzer.temporal_placebos(
            unit_col='unit_id',
            time_col='period',
            outcome_col='negative_control_outcome',  # No true effect
            treatment_col='treated',
            true_treatment_time=6,
            placebo_times=[2, 3, 4]
        )
        
        # Placebo effects should be small
        for time, estimate in result.placebo_estimates.items():
            assert abs(estimate) < 2.0, f"Placebo at t={time} too large: {estimate}"
    
    def test_multiple_placebo_times(self, sample_did_data):
        """Test with multiple placebo times."""
        analyzer = NegativeControlAnalyzer(sample_did_data)
        
        placebo_times = [1, 2, 3, 4]
        result = analyzer.temporal_placebos(
            unit_col='unit_id',
            time_col='period',
            outcome_col='main_outcome',
            treatment_col='treated',
            true_treatment_time=6,
            placebo_times=placebo_times
        )
        
        assert len(result.placebo_estimates) == len(placebo_times)


# ============================================================================
# Test Geographic Placebos
# ============================================================================

class TestGeographicPlacebos:
    """Tests for geographic placebo analysis."""
    
    def test_geographic_placebo_returns_result(self, analyzer):
        """Test that geographic placebo returns proper result structure."""
        result = analyzer.geographic_placebos(
            treatment_col='treatment',
            outcome_col='main_outcome',
            region_col='region',
            treated_regions=['A'],
            placebo_regions=['B', 'C']
        )
        
        assert isinstance(result, GeographicPlaceboResult)
        assert hasattr(result, 'placebo_estimates')
        assert hasattr(result, 'passed')
    
    def test_placebo_regions_no_effect(self, sample_procurement_data):
        """Test that placebo regions show no spillover effect."""
        # Create data where only region A is truly treated
        sample_procurement_data = sample_procurement_data.copy()
        
        # Only apply treatment effect in region A
        region_a_mask = sample_procurement_data['region'] == 'A'
        sample_procurement_data['treatment'] = 0
        sample_procurement_data.loc[region_a_mask, 'treatment'] = np.random.binomial(1, 0.5, region_a_mask.sum())
        
        # Outcome only affected in treated region
        sample_procurement_data['outcome'] = 10 + np.random.normal(0, 2, len(sample_procurement_data))
        sample_procurement_data.loc[
            (sample_procurement_data['treatment'] == 1) & region_a_mask,
            'outcome'
        ] += 5
        
        analyzer = NegativeControlAnalyzer(sample_procurement_data)
        result = analyzer.geographic_placebos(
            treatment_col='treatment',
            outcome_col='outcome',
            region_col='region',
            treated_regions=['A'],
            placebo_regions=['B', 'C']
        )
        
        # Placebo regions should have smaller effects
        for region, estimate in result.placebo_estimates.items():
            assert abs(estimate) < 3.0, f"Placebo region {region} effect too large: {estimate}"


# ============================================================================
# Test Effect Comparison Methods
# ============================================================================

class TestEffectComparison:
    """Tests for effect size comparison utilities."""
    
    def test_compare_to_main_effect(self, analyzer):
        """Test comparison between main and control effects."""
        main_result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='main_outcome',
            control_name='main'
        )
        
        control_result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='control'
        )
        
        comparison = analyzer.compare_effects(main_result, control_result)
        
        assert 'difference' in comparison
        assert 'ratio' in comparison
        assert comparison['difference'] > 0  # Main effect > control effect
    
    def test_equivalence_test(self, analyzer):
        """Test TOST equivalence testing."""
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='test'
        )
        
        # Test equivalence to zero with margin
        equivalence = analyzer.equivalence_test(
            result,
            equivalence_margin=1.0
        )
        
        assert 'equivalent_to_zero' in equivalence
        assert 'tost_pvalue' in equivalence


# ============================================================================
# Test Summary Report Generation
# ============================================================================

class TestSummaryReport:
    """Tests for summary report generation."""
    
    def test_generate_summary_report(self, analyzer):
        """Test that summary report is generated correctly."""
        # Run several negative control tests
        results = []
        
        for outcome in ['negative_control_outcome', 'unrelated_outcome']:
            result = analyzer.outcome_negative_controls(
                treatment_col='treatment',
                outcome_col=outcome,
                control_name=outcome
            )
            results.append(result)
        
        summary = analyzer.generate_summary_report(results)
        
        assert 'total_tests' in summary
        assert 'passed' in summary
        assert 'failed' in summary
        assert 'pass_rate' in summary
        assert summary['total_tests'] == 2
    
    def test_summary_to_dataframe(self, analyzer):
        """Test converting results to DataFrame."""
        results = []
        
        for outcome in ['negative_control_outcome', 'unrelated_outcome']:
            result = analyzer.outcome_negative_controls(
                treatment_col='treatment',
                outcome_col=outcome,
                control_name=outcome
            )
            results.append(result)
        
        df = analyzer.results_to_dataframe(results)
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert 'control_name' in df.columns
        assert 'estimate' in df.columns
        assert 'passed' in df.columns


# ============================================================================
# Test Edge Cases
# ============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_missing_column_error(self, analyzer):
        """Test error when required column is missing."""
        with pytest.raises((KeyError, ValueError)):
            analyzer.outcome_negative_controls(
                treatment_col='nonexistent_column',
                outcome_col='main_outcome',
                control_name='test'
            )
    
    def test_empty_data_error(self):
        """Test error with empty data."""
        analyzer = NegativeControlAnalyzer(pd.DataFrame())
        
        with pytest.raises((ValueError, IndexError)):
            analyzer.outcome_negative_controls(
                treatment_col='treatment',
                outcome_col='outcome',
                control_name='test'
            )
    
    def test_small_sample_warning(self, sample_procurement_data):
        """Test warning with very small sample."""
        small_data = sample_procurement_data.head(20)
        analyzer = NegativeControlAnalyzer(small_data)
        
        # Should still run but may have warnings
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='small_sample'
        )
        
        assert isinstance(result, NegativeControlResult)
    
    def test_perfect_collinearity_handling(self, sample_procurement_data):
        """Test handling of perfectly collinear variables."""
        # Add perfectly collinear covariate
        sample_procurement_data['collinear'] = sample_procurement_data['treatment'] * 2
        
        analyzer = NegativeControlAnalyzer(sample_procurement_data)
        
        # Should handle gracefully
        try:
            result = analyzer.outcome_negative_controls(
                treatment_col='treatment',
                outcome_col='negative_control_outcome',
                control_name='collinear_test',
                covariates=['collinear']
            )
            # Either succeeds or raises appropriate error
            assert isinstance(result, NegativeControlResult)
        except (ValueError, np.linalg.LinAlgError):
            # Acceptable to raise error for perfect collinearity
            pass


# ============================================================================
# Test Integration with RDD/DiD Modules
# ============================================================================

class TestIntegration:
    """Tests for integration with RDD and DiD modules."""
    
    def test_rdd_negative_control(self, sample_procurement_data):
        """Test negative control in RDD context."""
        # Add RDD-specific columns
        threshold = np.median(sample_procurement_data['contract_value'])
        sample_procurement_data['above_threshold'] = (
            sample_procurement_data['contract_value'] > threshold
        ).astype(int)
        
        analyzer = NegativeControlAnalyzer(sample_procurement_data)
        
        result = analyzer.outcome_negative_controls(
            treatment_col='above_threshold',
            outcome_col='negative_control_outcome',
            control_name='rdd_control'
        )
        
        assert isinstance(result, NegativeControlResult)
        # No true effect at threshold for this outcome
        assert result.passed
    
    def test_did_negative_control(self, sample_did_data):
        """Test negative control in DiD context."""
        analyzer = NegativeControlAnalyzer(sample_did_data)
        
        result = analyzer.outcome_negative_controls(
            treatment_col='treated',
            outcome_col='negative_control_outcome',
            control_name='did_control'
        )
        
        assert isinstance(result, NegativeControlResult)
        assert result.passed


# ============================================================================
# Test Reporting Format
# ============================================================================

class TestReportingFormat:
    """Tests for Nature Sustainability reporting format."""
    
    def test_latex_table_output(self, analyzer):
        """Test LaTeX table generation for negative controls."""
        results = []
        for outcome in ['negative_control_outcome', 'unrelated_outcome']:
            result = analyzer.outcome_negative_controls(
                treatment_col='treatment',
                outcome_col=outcome,
                control_name=outcome
            )
            results.append(result)
        
        latex = analyzer.to_latex_table(results)
        
        assert isinstance(latex, str)
        assert 'tabular' in latex or 'begin' in latex
        assert 'negative_control_outcome' in latex
    
    def test_interpretation_text(self, analyzer):
        """Test interpretation text generation."""
        result = analyzer.outcome_negative_controls(
            treatment_col='treatment',
            outcome_col='negative_control_outcome',
            control_name='test_control'
        )
        
        interpretation = analyzer.interpret_result(result)
        
        assert isinstance(interpretation, str)
        assert len(interpretation) > 0
        # Should mention pass/fail status
        assert 'pass' in interpretation.lower() or 'fail' in interpretation.lower()


# ============================================================================
# Run tests
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
