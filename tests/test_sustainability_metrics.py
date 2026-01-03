"""
Tests for Sustainability Metrics Module

Tests cover:
- Carbon footprint calculation (CPV → sector → CO2e/$)
- Resilience metrics (HHI, single-bid rate, completion)
- Green procurement classification
- Composite sustainability index
- Pipeline integration
"""

import pytest
import pandas as pd
import numpy as np
from unittest.mock import Mock, patch
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sustainability_metrics import (
    SectorEmissionsIntensity,
    EmissionsIntensityDatabase,
    CPVToSectorMapper,
    CarbonFootprintCalculator,
    ProcurementFootprint,
    ResilienceMetrics,
    ResilienceCalculator,
    GreenProcurementClassifier,
    SustainabilityScore,
    SustainabilityIndexCalculator,
    SustainabilityPipeline,
)


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def sample_contracts():
    """Sample procurement contracts DataFrame."""
    np.random.seed(42)
    return pd.DataFrame({
        "ocid": [f"contract_{i}" for i in range(50)],
        "value_amount": [100000, 50000, 200000, 75000, 150000] * 10,
        "cpv_code": [
            "45000000",  # Construction
            "72000000",  # IT services
            "33000000",  # Medical
            "09310000",  # Solar (green)
            "60000000",  # Transport
        ] * 10,
        "supplier_id": [f"supplier_{i % 5}" for i in range(50)],
        "n_bidders": [3, 1, 5, 2, 4] * 10,
        "title": ["Construction project"] * 10 + 
                 ["Software development"] * 10 +
                 ["Medical equipment"] * 10 +
                 ["Solar panel installation"] * 10 +
                 ["Transport services"] * 10,
        "description": ["Standard description"] * 50,
        "status": ["completed"] * 45 + ["cancelled"] * 5,
    })


@pytest.fixture
def emissions_db():
    """EmissionsIntensityDatabase instance."""
    return EmissionsIntensityDatabase()


@pytest.fixture
def cpv_mapper(emissions_db):
    """CPVToSectorMapper instance."""
    return CPVToSectorMapper(emissions_db)


@pytest.fixture
def carbon_calculator(cpv_mapper):
    """CarbonFootprintCalculator instance."""
    return CarbonFootprintCalculator(cpv_mapper)


@pytest.fixture
def resilience_calculator():
    """ResilienceCalculator instance."""
    return ResilienceCalculator()


@pytest.fixture
def green_classifier():
    """GreenProcurementClassifier instance."""
    return GreenProcurementClassifier()


@pytest.fixture
def index_calculator():
    """SustainabilityIndexCalculator instance."""
    return SustainabilityIndexCalculator()


@pytest.fixture
def pipeline():
    """Full SustainabilityPipeline instance."""
    return SustainabilityPipeline()


# =============================================================================
# EMISSIONS INTENSITY DATABASE TESTS
# =============================================================================

class TestEmissionsIntensityDatabase:
    """Tests for EmissionsIntensityDatabase."""
    
    def test_default_intensities_exist(self, emissions_db):
        """Verify default sector intensities are populated."""
        assert len(emissions_db.intensities) > 0
        assert "construction" in emissions_db.intensities
        assert "electricity" in emissions_db.intensities
        assert "software_services" in emissions_db.intensities
    
    def test_get_intensity_exact_match(self, emissions_db):
        """Test exact sector name matching."""
        intensity = emissions_db.get_intensity("construction")
        assert isinstance(intensity, SectorEmissionsIntensity)
        assert intensity.intensity_kg_co2e_per_usd > 0
    
    def test_get_intensity_partial_match(self, emissions_db):
        """Test partial sector name matching."""
        intensity = emissions_db.get_intensity("electricity_grid")
        # Should match "electricity"
        assert intensity.sector_name == "Electricity, gas, steam"
    
    def test_get_intensity_fallback(self, emissions_db):
        """Test fallback to 'other' for unknown sectors."""
        intensity = emissions_db.get_intensity("nonexistent_sector_xyz")
        assert intensity.sector_code == "X"
        assert intensity.sector_name == "Other/Unknown"
    
    def test_custom_intensities(self):
        """Test custom intensity override."""
        custom = {
            "custom_sector": SectorEmissionsIntensity(
                sector_code="CUSTOM",
                sector_name="Custom Sector",
                intensity_kg_co2e_per_usd=99.0,
                uncertainty_lower=90.0,
                uncertainty_upper=110.0,
                source="Test",
                year=2024
            )
        }
        db = EmissionsIntensityDatabase(custom_intensities=custom)
        
        intensity = db.get_intensity("custom_sector")
        assert intensity.intensity_kg_co2e_per_usd == 99.0
    
    def test_uncertainty_bounds(self, emissions_db):
        """Verify uncertainty bounds are sensible."""
        for name, intensity in emissions_db.intensities.items():
            assert intensity.uncertainty_lower <= intensity.intensity_kg_co2e_per_usd
            assert intensity.intensity_kg_co2e_per_usd <= intensity.uncertainty_upper


# =============================================================================
# CPV MAPPER TESTS
# =============================================================================

class TestCPVToSectorMapper:
    """Tests for CPVToSectorMapper."""
    
    def test_map_construction(self, cpv_mapper):
        """Test mapping construction CPV."""
        sector = cpv_mapper.map_cpv_to_sector("45000000")
        assert sector == "construction"
    
    def test_map_it_services(self, cpv_mapper):
        """Test mapping IT services CPV."""
        sector = cpv_mapper.map_cpv_to_sector("72000000")
        assert sector == "software_services"
    
    def test_map_medical(self, cpv_mapper):
        """Test mapping medical equipment CPV."""
        sector = cpv_mapper.map_cpv_to_sector("33000000")
        assert sector == "pharmaceuticals"
    
    def test_map_unknown(self, cpv_mapper):
        """Test unknown CPV falls back to 'other'."""
        sector = cpv_mapper.map_cpv_to_sector("99999999")
        assert sector == "other"
    
    def test_map_empty_code(self, cpv_mapper):
        """Test empty CPV code handling."""
        sector = cpv_mapper.map_cpv_to_sector("")
        assert sector == "other"
        
        sector = cpv_mapper.map_cpv_to_sector(None)
        assert sector == "other"
    
    def test_get_intensity_for_cpv(self, cpv_mapper):
        """Test full CPV to intensity lookup."""
        intensity = cpv_mapper.get_intensity_for_cpv("45000000")
        assert intensity.sector_name == "Construction"
        assert intensity.intensity_kg_co2e_per_usd == 0.42


# =============================================================================
# CARBON FOOTPRINT CALCULATOR TESTS
# =============================================================================

class TestCarbonFootprintCalculator:
    """Tests for CarbonFootprintCalculator."""
    
    def test_calculate_single_footprint(self, carbon_calculator):
        """Test single contract footprint calculation."""
        footprint = carbon_calculator.calculate_footprint(
            contract_id="test_001",
            spend_usd=100000,
            cpv_code="45000000"  # Construction
        )
        
        assert isinstance(footprint, ProcurementFootprint)
        assert footprint.contract_id == "test_001"
        assert footprint.spend_usd == 100000
        assert footprint.sector == "construction"
        # Construction intensity is 0.42 kg CO2e/$
        assert footprint.footprint_kg_co2e == 100000 * 0.42
    
    def test_footprint_uncertainty(self, carbon_calculator):
        """Test uncertainty bounds propagation."""
        footprint = carbon_calculator.calculate_footprint(
            contract_id="test_002",
            spend_usd=100000,
            cpv_code="45000000"
        )
        
        assert footprint.footprint_lower < footprint.footprint_kg_co2e
        assert footprint.footprint_kg_co2e < footprint.footprint_upper
    
    def test_calculate_batch(self, carbon_calculator, sample_contracts):
        """Test batch footprint calculation."""
        footprints = carbon_calculator.calculate_batch(
            sample_contracts,
            id_col="ocid",
            value_col="value_amount",
            cpv_col="cpv_code"
        )
        
        assert len(footprints) == len(sample_contracts)
        assert "footprint_kg_co2e" in footprints.columns
        assert footprints["footprint_kg_co2e"].sum() > 0
    
    def test_aggregate_footprint_total(self, carbon_calculator, sample_contracts):
        """Test total footprint aggregation."""
        footprints = carbon_calculator.calculate_batch(sample_contracts)
        aggregated = carbon_calculator.aggregate_footprint(footprints)
        
        assert len(aggregated) == 1
        assert aggregated["total_spend_usd"].iloc[0] == sample_contracts["value_amount"].sum()
        assert aggregated["n_contracts"].iloc[0] == len(sample_contracts)
    
    def test_aggregate_footprint_by_sector(self, carbon_calculator, sample_contracts):
        """Test footprint aggregation by sector."""
        footprints = carbon_calculator.calculate_batch(sample_contracts)
        aggregated = carbon_calculator.aggregate_footprint(
            footprints,
            group_by=["sector"]
        )
        
        assert "sector" in aggregated.columns
        assert len(aggregated) > 1  # Multiple sectors


# =============================================================================
# RESILIENCE CALCULATOR TESTS
# =============================================================================

class TestResilienceCalculator:
    """Tests for ResilienceCalculator."""
    
    def test_supplier_concentration_uniform(self, resilience_calculator):
        """Test HHI for uniformly distributed suppliers."""
        contracts = pd.DataFrame({
            "supplier_id": [f"s_{i}" for i in range(10)],
            "value_amount": [100000] * 10
        })
        
        hhi, top3, n_suppliers = resilience_calculator.calculate_supplier_concentration(
            contracts
        )
        
        # 10 equal suppliers: HHI = 10 * (0.1)^2 * 10000 = 1000
        assert n_suppliers == 10
        assert abs(hhi - 1000) < 1  # Allow small float error
        assert abs(top3 - 0.3) < 0.01
    
    def test_supplier_concentration_monopoly(self, resilience_calculator):
        """Test HHI for monopoly situation."""
        contracts = pd.DataFrame({
            "supplier_id": ["monopolist"] * 10,
            "value_amount": [100000] * 10
        })
        
        hhi, top3, n_suppliers = resilience_calculator.calculate_supplier_concentration(
            contracts
        )
        
        # Monopoly: HHI = 10000
        assert n_suppliers == 1
        assert hhi == 10000
        assert top3 == 1.0
    
    def test_competition_metrics(self, resilience_calculator, sample_contracts):
        """Test competition metrics calculation."""
        avg_bidders, single_bid = resilience_calculator.calculate_competition_metrics(
            sample_contracts
        )
        
        assert avg_bidders > 0
        assert 0 <= single_bid <= 1
        # Our sample has 10 single-bid contracts (10/50 = 0.2)
        assert single_bid == 0.2
    
    def test_calculate_resilience_full(self, resilience_calculator, sample_contracts):
        """Test full resilience metrics calculation."""
        metrics = resilience_calculator.calculate_resilience(
            sample_contracts,
            entity_id="test_authority",
            period="2024-Q1"
        )
        
        assert isinstance(metrics, ResilienceMetrics)
        assert metrics.entity_id == "test_authority"
        assert metrics.period == "2024-Q1"
        assert metrics.supplier_hhi is not None
        assert metrics.avg_bidders is not None
    
    def test_empty_contracts(self, resilience_calculator):
        """Test handling of empty contracts."""
        empty = pd.DataFrame()
        
        hhi, top3, n = resilience_calculator.calculate_supplier_concentration(empty)
        assert np.isnan(hhi)
        assert n == 0


# =============================================================================
# GREEN PROCUREMENT CLASSIFIER TESTS
# =============================================================================

class TestGreenProcurementClassifier:
    """Tests for GreenProcurementClassifier."""
    
    def test_is_green_cpv_solar(self, green_classifier):
        """Test solar CPV classification."""
        assert green_classifier.is_green_cpv("09310000") is True
        assert green_classifier.is_green_cpv("09331200") is True  # PV modules
    
    def test_is_green_cpv_recycling(self, green_classifier):
        """Test recycling CPV classification."""
        assert green_classifier.is_green_cpv("90514000") is True  # Recycling
    
    def test_is_green_cpv_regular(self, green_classifier):
        """Test non-green CPV."""
        assert green_classifier.is_green_cpv("45000000") is False  # Construction
        assert green_classifier.is_green_cpv("72000000") is False  # IT
    
    def test_is_green_text_keywords(self, green_classifier):
        """Test keyword-based classification."""
        assert green_classifier.is_green_text("Solar panel installation") is True
        assert green_classifier.is_green_text("Energy efficient building") is True
        assert green_classifier.is_green_text("Renewable energy project") is True
        assert green_classifier.is_green_text("Electric vehicle fleet") is True
    
    def test_is_green_text_negative(self, green_classifier):
        """Test non-green text."""
        assert green_classifier.is_green_text("Standard office supplies") is False
        assert green_classifier.is_green_text("Road construction") is False
    
    def test_classify_combined(self, green_classifier):
        """Test combined CPV + text classification."""
        # Green CPV
        is_green, reason = green_classifier.classify(cpv_code="09310000")
        assert is_green is True
        assert "cpv:" in reason
        
        # Green title
        is_green, reason = green_classifier.classify(title="Renewable energy project")
        assert is_green is True
        assert "title_keyword" in reason
        
        # Both
        is_green, reason = green_classifier.classify(
            cpv_code="09310000",
            title="Solar installation"
        )
        assert is_green is True
        assert "cpv:" in reason
        assert "title_keyword" in reason
    
    def test_calculate_green_share(self, green_classifier, sample_contracts):
        """Test green share calculation."""
        result = green_classifier.calculate_green_share(sample_contracts)
        
        assert "green_count_share" in result
        assert "green_value_share" in result
        assert result["n_total"] == len(sample_contracts)
        # 10 solar contracts out of 50
        assert result["n_green"] >= 10
    
    def test_empty_text_handling(self, green_classifier):
        """Test empty/None text handling."""
        assert green_classifier.is_green_text("") is False
        assert green_classifier.is_green_text(None) is False


# =============================================================================
# SUSTAINABILITY INDEX CALCULATOR TESTS
# =============================================================================

class TestSustainabilityIndexCalculator:
    """Tests for SustainabilityIndexCalculator."""
    
    def test_score_carbon_excellent(self, index_calculator):
        """Test carbon scoring for low intensity."""
        score = index_calculator.score_carbon(0.10)  # Below excellent
        assert score == 100.0
    
    def test_score_carbon_poor(self, index_calculator):
        """Test carbon scoring for high intensity."""
        score = index_calculator.score_carbon(0.90)  # Above poor
        assert score == 0.0
    
    def test_score_carbon_middle(self, index_calculator):
        """Test carbon scoring for middle intensity."""
        score = index_calculator.score_carbon(0.475)  # Midpoint
        assert 40 < score < 60
    
    def test_score_resilience(self, index_calculator):
        """Test resilience scoring."""
        metrics = ResilienceMetrics(
            entity_id="test",
            period="2024",
            supplier_hhi=1500,  # Unconcentrated
            single_bid_rate=0.1,  # Low
            completion_rate=0.95,  # High
            avg_bidders=4.0  # Good
        )
        
        score = index_calculator.score_resilience(metrics)
        assert score > 70  # Should be high for good metrics
    
    def test_score_green(self, index_calculator):
        """Test green share scoring."""
        # 30% green = 100 score
        assert index_calculator.score_green(0.30) == 100.0
        # 15% green = 50 score
        assert abs(index_calculator.score_green(0.15) - 50.0) < 1
        # 0% green = 0 score
        assert index_calculator.score_green(0.0) == 0.0
    
    def test_calculate_composite(self, index_calculator):
        """Test composite score calculation."""
        metrics = ResilienceMetrics(
            entity_id="test",
            period="2024",
            supplier_hhi=1500,
            single_bid_rate=0.1,
            avg_bidders=4.0
        )
        
        score = index_calculator.calculate(
            entity_id="test",
            period="2024",
            carbon_intensity=0.20,
            resilience_metrics=metrics,
            green_share=0.20
        )
        
        assert isinstance(score, SustainabilityScore)
        assert score.carbon_score > 0
        assert score.resilience_score > 0
        assert score.green_score > 0
        assert 0 <= score.composite_score <= 100


# =============================================================================
# SUSTAINABILITY PIPELINE TESTS
# =============================================================================

class TestSustainabilityPipeline:
    """Tests for SustainabilityPipeline."""
    
    def test_process_contracts(self, pipeline, sample_contracts):
        """Test full pipeline processing."""
        results = pipeline.process(
            sample_contracts,
            entity_id="test_authority",
            period="2024-Q1"
        )
        
        assert results["entity_id"] == "test_authority"
        assert results["period"] == "2024-Q1"
        assert results["n_contracts"] == len(sample_contracts)
        
        # Carbon metrics
        assert "carbon" in results
        assert "intensity_kg_per_usd" in results["carbon"]
        
        # Resilience metrics
        assert "resilience" in results
        assert "supplier_hhi" in results["resilience"]
        
        # Green metrics
        assert "green" in results
        assert "green_value_share" in results["green"]
        
        # Sustainability score
        assert "sustainability_score" in results
        assert "composite_score" in results["sustainability_score"]
    
    def test_compare_periods(self, pipeline):
        """Test period comparison for DiD-style analysis."""
        # Create contracts with period labels
        contracts = pd.DataFrame({
            "ocid": [f"c_{i}" for i in range(100)],
            "value_amount": [100000] * 100,
            "cpv_code": ["45000000"] * 50 + ["09310000"] * 50,
            "supplier_id": [f"s_{i % 10}" for i in range(100)],
            "n_bidders": [2] * 50 + [4] * 50,  # More competition post
            "title": ["Construction"] * 50 + ["Solar energy"] * 50,
            "description": [""] * 100,
            "period": ["pre"] * 50 + ["post"] * 50
        })
        
        comparison = pipeline.compare_periods(
            contracts,
            entity_id="test",
            period_col="period",
            pre_period="pre",
            post_period="post"
        )
        
        assert "pre" in comparison
        assert "post" in comparison
        assert "difference" in comparison
        
        # Post should have better sustainability (more green, more competition)
        diff = comparison["difference"]
        assert "sustainability_score_change" in diff
    
    def test_column_mapping(self, pipeline):
        """Test custom column mapping."""
        contracts = pd.DataFrame({
            "contract_id": ["c1", "c2"],
            "amount": [100000, 200000],
            "cpv": ["45000000", "72000000"],
            "vendor": ["s1", "s2"],
            "bids": [3, 2],
            "name": ["Project A", "Project B"],
            "desc": ["", ""]
        })
        
        results = pipeline.process(
            contracts,
            entity_id="test",
            period="2024",
            column_mapping={
                "id": "contract_id",
                "value": "amount",
                "cpv": "cpv",
                "supplier": "vendor",
                "bidders": "bids",
                "title": "name",
                "description": "desc"
            }
        )
        
        assert results["n_contracts"] == 2


# =============================================================================
# EDGE CASES AND ERROR HANDLING
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_zero_value_contracts(self, carbon_calculator):
        """Test handling of zero-value contracts."""
        footprint = carbon_calculator.calculate_footprint(
            contract_id="zero",
            spend_usd=0,
            cpv_code="45000000"
        )
        assert footprint.footprint_kg_co2e == 0
    
    def test_missing_columns(self, pipeline):
        """Test handling of missing columns."""
        minimal = pd.DataFrame({
            "ocid": ["c1"],
            "value_amount": [100000],
            "cpv_code": ["45000000"]
        })
        
        # Should not raise, should handle missing columns gracefully
        results = pipeline.process(minimal, "test", "2024")
        assert results["n_contracts"] == 1
    
    def test_nan_handling(self, index_calculator):
        """Test NaN handling in scoring."""
        assert np.isnan(index_calculator.score_carbon(np.nan))
        assert np.isnan(index_calculator.score_green(np.nan))


# =============================================================================
# INTEGRATION TESTS
# =============================================================================

class TestIntegration:
    """Integration tests for full workflow."""
    
    def test_full_workflow(self):
        """Test complete workflow from raw data to sustainability score."""
        # Create realistic sample data
        np.random.seed(42)
        n_contracts = 200
        
        contracts = pd.DataFrame({
            "ocid": [f"ocid-{i:04d}" for i in range(n_contracts)],
            "value_amount": np.random.lognormal(11, 1, n_contracts),
            "cpv_code": np.random.choice([
                "45000000", "72000000", "33000000", 
                "09310000", "90514000", "60000000"
            ], n_contracts, p=[0.3, 0.2, 0.15, 0.15, 0.1, 0.1]),
            "supplier_id": np.random.choice(
                [f"supplier_{i}" for i in range(30)], n_contracts
            ),
            "n_bidders": np.random.poisson(3, n_contracts) + 1,
            "title": np.random.choice([
                "Construction works",
                "IT development",
                "Medical supplies",
                "Solar installation",
                "Waste recycling",
                "Transport services"
            ], n_contracts),
            "description": [""] * n_contracts,
            "status": np.random.choice(
                ["completed", "active", "cancelled"], n_contracts,
                p=[0.85, 0.10, 0.05]
            )
        })
        
        # Process through pipeline
        pipeline = SustainabilityPipeline()
        results = pipeline.process(
            contracts,
            entity_id="integration_test",
            period="2024"
        )
        
        # Verify all components present
        assert results["carbon"]["intensity_kg_per_usd"] > 0
        assert results["resilience"]["supplier_hhi"] > 0
        assert results["green"]["n_green"] > 0
        assert 0 <= results["sustainability_score"]["composite_score"] <= 100
        
        # Verify reasonable values
        assert results["carbon"]["intensity_kg_per_usd"] < 5  # Not astronomical
        assert results["resilience"]["supplier_hhi"] < 10000  # Not monopoly


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
