"""
Sustainability Metrics Module for Procurement Analysis

This module implements sustainability outcome metrics for Nature Sustainability submission:
1. Carbon/Environmental Footprint - CO2e intensity of procurement spend
2. Resilience/Reliability - Supply chain robustness metrics
3. Green Procurement Share - Sustainability-tagged contract proportion

These transform "rules → competition" into "rules → market structure → sustainability outcomes"
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from pathlib import Path
import json
import logging
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


# =============================================================================
# CARBON FOOTPRINT METRICS
# =============================================================================

@dataclass
class SectorEmissionsIntensity:
    """
    Sector-level CO2e emissions intensity (kg CO2e per USD spent).
    
    Based on environmentally-extended input-output (EEIO) analysis.
    Sources:
    - EXIOBASE 3 (multi-regional EEIO database)
    - US EPA USEEIO model
    - OECD Inter-Country Input-Output (ICIO) tables
    """
    sector_code: str
    sector_name: str
    intensity_kg_co2e_per_usd: float
    uncertainty_lower: float  # 95% CI lower bound
    uncertainty_upper: float  # 95% CI upper bound
    source: str
    year: int
    region: str = "global"


class EmissionsIntensityDatabase:
    """
    Database of sector-level emissions intensities for procurement footprinting.
    
    Implements the approach from:
    - Stadler et al. (2018) EXIOBASE 3
    - Yang et al. (2017) US EPA USEEIO
    """
    
    # Default intensities by broad sector (kg CO2e / USD, 2020 values)
    # Source: Aggregated from EXIOBASE 3 + USEEIO
    DEFAULT_INTENSITIES = {
        # Construction & Infrastructure
        "construction": SectorEmissionsIntensity(
            sector_code="F", sector_name="Construction",
            intensity_kg_co2e_per_usd=0.42,
            uncertainty_lower=0.35, uncertainty_upper=0.52,
            source="EXIOBASE3", year=2020
        ),
        # Energy & Utilities
        "electricity": SectorEmissionsIntensity(
            sector_code="D35", sector_name="Electricity, gas, steam",
            intensity_kg_co2e_per_usd=1.85,
            uncertainty_lower=1.20, uncertainty_upper=2.80,
            source="EXIOBASE3", year=2020
        ),
        "fossil_fuels": SectorEmissionsIntensity(
            sector_code="B", sector_name="Mining and quarrying",
            intensity_kg_co2e_per_usd=2.10,
            uncertainty_lower=1.50, uncertainty_upper=3.00,
            source="EXIOBASE3", year=2020
        ),
        # Manufacturing
        "manufacturing_heavy": SectorEmissionsIntensity(
            sector_code="C24", sector_name="Basic metals",
            intensity_kg_co2e_per_usd=1.45,
            uncertainty_lower=1.10, uncertainty_upper=1.90,
            source="EXIOBASE3", year=2020
        ),
        "manufacturing_light": SectorEmissionsIntensity(
            sector_code="C", sector_name="Manufacturing (other)",
            intensity_kg_co2e_per_usd=0.38,
            uncertainty_lower=0.28, uncertainty_upper=0.52,
            source="EXIOBASE3", year=2020
        ),
        "pharmaceuticals": SectorEmissionsIntensity(
            sector_code="C21", sector_name="Pharmaceuticals",
            intensity_kg_co2e_per_usd=0.22,
            uncertainty_lower=0.15, uncertainty_upper=0.32,
            source="EXIOBASE3", year=2020
        ),
        # IT & Electronics
        "it_equipment": SectorEmissionsIntensity(
            sector_code="C26", sector_name="Computer, electronic products",
            intensity_kg_co2e_per_usd=0.35,
            uncertainty_lower=0.25, uncertainty_upper=0.48,
            source="EXIOBASE3", year=2020
        ),
        "software_services": SectorEmissionsIntensity(
            sector_code="J62", sector_name="IT services",
            intensity_kg_co2e_per_usd=0.08,
            uncertainty_lower=0.05, uncertainty_upper=0.12,
            source="EXIOBASE3", year=2020
        ),
        # Transport
        "transport_road": SectorEmissionsIntensity(
            sector_code="H49", sector_name="Land transport",
            intensity_kg_co2e_per_usd=0.65,
            uncertainty_lower=0.50, uncertainty_upper=0.85,
            source="EXIOBASE3", year=2020
        ),
        "transport_air": SectorEmissionsIntensity(
            sector_code="H51", sector_name="Air transport",
            intensity_kg_co2e_per_usd=1.20,
            uncertainty_lower=0.90, uncertainty_upper=1.60,
            source="EXIOBASE3", year=2020
        ),
        # Health & Social
        "healthcare": SectorEmissionsIntensity(
            sector_code="Q", sector_name="Human health activities",
            intensity_kg_co2e_per_usd=0.18,
            uncertainty_lower=0.12, uncertainty_upper=0.26,
            source="EXIOBASE3", year=2020
        ),
        # Professional Services
        "professional_services": SectorEmissionsIntensity(
            sector_code="M", sector_name="Professional services",
            intensity_kg_co2e_per_usd=0.10,
            uncertainty_lower=0.06, uncertainty_upper=0.15,
            source="EXIOBASE3", year=2020
        ),
        # Food & Agriculture
        "food_products": SectorEmissionsIntensity(
            sector_code="C10", sector_name="Food products",
            intensity_kg_co2e_per_usd=0.55,
            uncertainty_lower=0.40, uncertainty_upper=0.75,
            source="EXIOBASE3", year=2020
        ),
        # Default fallback
        "other": SectorEmissionsIntensity(
            sector_code="X", sector_name="Other/Unknown",
            intensity_kg_co2e_per_usd=0.30,
            uncertainty_lower=0.15, uncertainty_upper=0.50,
            source="Average", year=2020
        )
    }
    
    def __init__(self, custom_intensities: Optional[Dict] = None):
        """
        Initialize emissions database.
        
        Args:
            custom_intensities: Optional dict of custom sector intensities
        """
        self.intensities = self.DEFAULT_INTENSITIES.copy()
        if custom_intensities:
            self.intensities.update(custom_intensities)
    
    def get_intensity(self, sector: str) -> SectorEmissionsIntensity:
        """Get emissions intensity for a sector."""
        sector_lower = sector.lower().replace(" ", "_").replace("-", "_")
        
        # Try exact match first
        if sector_lower in self.intensities:
            return self.intensities[sector_lower]
        
        # Try partial matching
        for key, intensity in self.intensities.items():
            if key in sector_lower or sector_lower in key:
                return intensity
        
        # Return default
        logger.warning(f"No intensity found for sector '{sector}', using default")
        return self.intensities["other"]


class CPVToSectorMapper:
    """
    Maps CPV (Common Procurement Vocabulary) codes to emission sectors.
    
    CPV is the EU standard classification for public procurement.
    This mapper enables carbon footprinting of procurement data.
    """
    
    # CPV Division (first 2 digits) to sector mapping
    CPV_DIVISION_MAP = {
        # Agricultural products
        "03": "food_products",
        # Petroleum, fuel, electricity
        "09": "fossil_fuels",
        # Mining, basic metals
        "14": "manufacturing_heavy",
        # Food, beverages, tobacco
        "15": "food_products",
        # Textiles, clothing
        "18": "manufacturing_light",
        # Leather, footwear
        "19": "manufacturing_light",
        # Wood, paper products
        "20": "manufacturing_light",
        "21": "manufacturing_light",
        "22": "manufacturing_light",
        # Chemicals, pharmaceuticals
        "24": "manufacturing_heavy",
        # Rubber, plastics
        "25": "manufacturing_light",
        # Glass, ceramics
        "26": "manufacturing_heavy",
        # Fabricated metals
        "27": "manufacturing_heavy",
        # Machinery, equipment
        "28": "manufacturing_light",
        "29": "manufacturing_light",
        # Office machinery, computers
        "30": "it_equipment",
        # Electrical equipment
        "31": "it_equipment",
        # Radio, TV, communications
        "32": "it_equipment",
        # Medical equipment
        "33": "pharmaceuticals",
        # Transport equipment
        "34": "manufacturing_light",
        # Security, defense
        "35": "manufacturing_light",
        # Musical instruments, sports
        "37": "manufacturing_light",
        # Miscellaneous products
        "38": "manufacturing_light",
        "39": "manufacturing_light",
        # Construction works
        "45": "construction",
        # Repair, maintenance
        "50": "professional_services",
        # Installation services
        "51": "construction",
        # Hotel, restaurant services
        "55": "professional_services",
        # Transport services
        "60": "transport_road",
        "61": "transport_air",
        "62": "transport_road",
        "63": "transport_road",
        # Postal, telecommunications
        "64": "software_services",
        # Utilities
        "65": "electricity",
        # Financial services
        "66": "professional_services",
        # Real estate
        "70": "professional_services",
        # Architectural, engineering
        "71": "professional_services",
        # IT services
        "72": "software_services",
        # R&D services
        "73": "professional_services",
        # Business services
        "75": "professional_services",
        "76": "professional_services",
        "77": "professional_services",
        # Admin, defense, social
        "78": "professional_services",
        "79": "professional_services",
        # Health, social services
        "80": "healthcare",
        "85": "healthcare",
        # Recreation, culture
        "92": "professional_services",
        # Other services
        "90": "professional_services",
        "98": "professional_services",
    }
    
    def __init__(self, emissions_db: Optional[EmissionsIntensityDatabase] = None):
        """
        Initialize mapper with emissions database.
        
        Args:
            emissions_db: Optional EmissionsIntensityDatabase instance
        """
        self.emissions_db = emissions_db or EmissionsIntensityDatabase()
    
    def map_cpv_to_sector(self, cpv_code: str) -> str:
        """
        Map CPV code to emissions sector.
        
        Args:
            cpv_code: CPV code (e.g., "45000000" for construction)
            
        Returns:
            Sector key for emissions lookup
        """
        if not cpv_code or len(cpv_code) < 2:
            return "other"
        
        division = cpv_code[:2]
        return self.CPV_DIVISION_MAP.get(division, "other")
    
    def get_intensity_for_cpv(self, cpv_code: str) -> SectorEmissionsIntensity:
        """
        Get emissions intensity for a CPV code.
        
        Args:
            cpv_code: CPV code
            
        Returns:
            SectorEmissionsIntensity with CO2e/$ values
        """
        sector = self.map_cpv_to_sector(cpv_code)
        return self.emissions_db.get_intensity(sector)


@dataclass
class ProcurementFootprint:
    """Carbon footprint calculation for a procurement contract."""
    contract_id: str
    spend_usd: float
    cpv_code: str
    sector: str
    intensity: SectorEmissionsIntensity
    footprint_kg_co2e: float
    footprint_lower: float  # 95% CI
    footprint_upper: float  # 95% CI
    
    def to_dict(self) -> Dict:
        return {
            "contract_id": self.contract_id,
            "spend_usd": self.spend_usd,
            "cpv_code": self.cpv_code,
            "sector": self.sector,
            "intensity_kg_per_usd": self.intensity.intensity_kg_co2e_per_usd,
            "footprint_kg_co2e": self.footprint_kg_co2e,
            "footprint_lower": self.footprint_lower,
            "footprint_upper": self.footprint_upper,
        }


class CarbonFootprintCalculator:
    """
    Calculate carbon footprint of procurement spending.
    
    Implements the approach for Nature Sustainability:
    - Map procurement categories to sectors
    - Apply sector-level emissions intensities
    - Propagate uncertainty through calculations
    - Enable causal analysis of rule changes on footprint
    """
    
    def __init__(self, mapper: Optional[CPVToSectorMapper] = None):
        """
        Initialize calculator.
        
        Args:
            mapper: Optional CPVToSectorMapper instance
        """
        self.mapper = mapper or CPVToSectorMapper()
    
    def calculate_footprint(
        self,
        contract_id: str,
        spend_usd: float,
        cpv_code: str
    ) -> ProcurementFootprint:
        """
        Calculate carbon footprint for a single contract.
        
        Args:
            contract_id: Unique contract identifier
            spend_usd: Contract value in USD
            cpv_code: CPV classification code
            
        Returns:
            ProcurementFootprint with CO2e estimates and uncertainty
        """
        sector = self.mapper.map_cpv_to_sector(cpv_code)
        intensity = self.mapper.get_intensity_for_cpv(cpv_code)
        
        footprint = spend_usd * intensity.intensity_kg_co2e_per_usd
        footprint_lower = spend_usd * intensity.uncertainty_lower
        footprint_upper = spend_usd * intensity.uncertainty_upper
        
        return ProcurementFootprint(
            contract_id=contract_id,
            spend_usd=spend_usd,
            cpv_code=cpv_code,
            sector=sector,
            intensity=intensity,
            footprint_kg_co2e=footprint,
            footprint_lower=footprint_lower,
            footprint_upper=footprint_upper
        )
    
    def calculate_batch(
        self,
        contracts: pd.DataFrame,
        id_col: str = "ocid",
        value_col: str = "value_amount",
        cpv_col: str = "cpv_code"
    ) -> pd.DataFrame:
        """
        Calculate footprints for a batch of contracts.
        
        Args:
            contracts: DataFrame with contract data
            id_col: Column name for contract ID
            value_col: Column name for contract value
            cpv_col: Column name for CPV code
            
        Returns:
            DataFrame with footprint calculations
        """
        results = []
        
        for _, row in contracts.iterrows():
            try:
                footprint = self.calculate_footprint(
                    contract_id=str(row.get(id_col, "")),
                    spend_usd=float(row.get(value_col, 0)),
                    cpv_code=str(row.get(cpv_col, ""))
                )
                results.append(footprint.to_dict())
            except Exception as e:
                logger.warning(f"Failed to calculate footprint for {row.get(id_col)}: {e}")
                continue
        
        return pd.DataFrame(results)
    
    def aggregate_footprint(
        self,
        footprints: pd.DataFrame,
        group_by: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Aggregate footprints by specified grouping.
        
        Args:
            footprints: DataFrame with footprint calculations
            group_by: Optional list of columns to group by
            
        Returns:
            Aggregated footprint statistics
        """
        if group_by is None:
            # Total aggregation
            return pd.DataFrame([{
                "total_spend_usd": footprints["spend_usd"].sum(),
                "total_footprint_kg_co2e": footprints["footprint_kg_co2e"].sum(),
                "total_footprint_lower": footprints["footprint_lower"].sum(),
                "total_footprint_upper": footprints["footprint_upper"].sum(),
                "mean_intensity": footprints["intensity_kg_per_usd"].mean(),
                "n_contracts": len(footprints)
            }])
        else:
            # Group aggregation
            return footprints.groupby(group_by).agg({
                "spend_usd": "sum",
                "footprint_kg_co2e": "sum",
                "footprint_lower": "sum",
                "footprint_upper": "sum",
                "intensity_kg_per_usd": "mean",
                "contract_id": "count"
            }).rename(columns={"contract_id": "n_contracts"}).reset_index()


# =============================================================================
# RESILIENCE METRICS
# =============================================================================

@dataclass
class ResilienceMetrics:
    """
    Supply chain resilience metrics for procurement.
    
    Captures dimensions of procurement system robustness:
    - Delivery reliability
    - Supplier diversification
    - Price stability
    - Market concentration
    """
    entity_id: str  # Contracting authority or aggregation unit
    period: str  # Time period (e.g., "2020-Q1")
    
    # Delivery reliability
    on_time_delivery_rate: Optional[float] = None  # % contracts delivered on time
    completion_rate: Optional[float] = None  # % contracts completed vs cancelled
    avg_delivery_delay_days: Optional[float] = None  # Mean delay from planned
    
    # Supplier diversification
    n_unique_suppliers: Optional[int] = None
    supplier_hhi: Optional[float] = None  # Herfindahl-Hirschman Index (0-10000)
    top3_supplier_share: Optional[float] = None  # Market share of top 3 suppliers
    new_supplier_rate: Optional[float] = None  # % contracts to new suppliers
    
    # Price stability
    price_cv: Optional[float] = None  # Coefficient of variation of unit prices
    price_trend: Optional[float] = None  # % change in prices over period
    
    # Competition
    avg_bidders: Optional[float] = None
    single_bid_rate: Optional[float] = None  # % single-bidder contracts
    
    def to_dict(self) -> Dict:
        return {
            "entity_id": self.entity_id,
            "period": self.period,
            "on_time_delivery_rate": self.on_time_delivery_rate,
            "completion_rate": self.completion_rate,
            "avg_delivery_delay_days": self.avg_delivery_delay_days,
            "n_unique_suppliers": self.n_unique_suppliers,
            "supplier_hhi": self.supplier_hhi,
            "top3_supplier_share": self.top3_supplier_share,
            "new_supplier_rate": self.new_supplier_rate,
            "price_cv": self.price_cv,
            "price_trend": self.price_trend,
            "avg_bidders": self.avg_bidders,
            "single_bid_rate": self.single_bid_rate,
        }


class ResilienceCalculator:
    """
    Calculate supply chain resilience metrics from procurement data.
    
    These metrics enable causal analysis of how procurement rules
    affect system resilience - a key sustainability outcome.
    """
    
    def calculate_supplier_concentration(
        self,
        contracts: pd.DataFrame,
        supplier_col: str = "supplier_id",
        value_col: str = "value_amount"
    ) -> Tuple[float, float, int]:
        """
        Calculate supplier concentration metrics.
        
        Args:
            contracts: DataFrame with contract data
            supplier_col: Column with supplier identifier
            value_col: Column with contract value
            
        Returns:
            Tuple of (HHI, top3_share, n_unique_suppliers)
        """
        if contracts.empty:
            return np.nan, np.nan, 0
        
        # Calculate market shares
        supplier_totals = contracts.groupby(supplier_col)[value_col].sum()
        total_value = supplier_totals.sum()
        
        if total_value == 0:
            return np.nan, np.nan, 0
        
        market_shares = supplier_totals / total_value
        
        # HHI (sum of squared market shares, scaled to 0-10000)
        hhi = (market_shares ** 2).sum() * 10000
        
        # Top 3 share
        top3_share = market_shares.nlargest(3).sum()
        
        # Unique suppliers
        n_suppliers = len(supplier_totals)
        
        return hhi, top3_share, n_suppliers
    
    def calculate_delivery_reliability(
        self,
        contracts: pd.DataFrame,
        planned_date_col: str = "planned_delivery_date",
        actual_date_col: str = "actual_delivery_date",
        status_col: str = "status"
    ) -> Tuple[float, float, float]:
        """
        Calculate delivery reliability metrics.
        
        Args:
            contracts: DataFrame with contract data
            planned_date_col: Column with planned delivery date
            actual_date_col: Column with actual delivery date
            status_col: Column with contract status
            
        Returns:
            Tuple of (on_time_rate, completion_rate, avg_delay_days)
        """
        if contracts.empty:
            return np.nan, np.nan, np.nan
        
        # Completion rate
        if status_col in contracts.columns:
            completed = contracts[status_col].isin(["complete", "completed", "closed"])
            completion_rate = completed.mean()
        else:
            completion_rate = np.nan
        
        # On-time delivery and delays
        if planned_date_col in contracts.columns and actual_date_col in contracts.columns:
            valid_dates = contracts.dropna(subset=[planned_date_col, actual_date_col])
            
            if not valid_dates.empty:
                planned = pd.to_datetime(valid_dates[planned_date_col])
                actual = pd.to_datetime(valid_dates[actual_date_col])
                delays = (actual - planned).dt.days
                
                on_time_rate = (delays <= 0).mean()
                avg_delay = delays[delays > 0].mean() if (delays > 0).any() else 0
            else:
                on_time_rate = np.nan
                avg_delay = np.nan
        else:
            on_time_rate = np.nan
            avg_delay = np.nan
        
        return on_time_rate, completion_rate, avg_delay
    
    def calculate_competition_metrics(
        self,
        contracts: pd.DataFrame,
        bidders_col: str = "n_bidders"
    ) -> Tuple[float, float]:
        """
        Calculate competition metrics.
        
        Args:
            contracts: DataFrame with contract data
            bidders_col: Column with number of bidders
            
        Returns:
            Tuple of (avg_bidders, single_bid_rate)
        """
        if contracts.empty or bidders_col not in contracts.columns:
            return np.nan, np.nan
        
        valid = contracts[bidders_col].dropna()
        
        if valid.empty:
            return np.nan, np.nan
        
        avg_bidders = valid.mean()
        single_bid_rate = (valid == 1).mean()
        
        return avg_bidders, single_bid_rate
    
    def calculate_price_stability(
        self,
        contracts: pd.DataFrame,
        price_col: str = "unit_price",
        date_col: str = "award_date"
    ) -> Tuple[float, float]:
        """
        Calculate price stability metrics.
        
        Args:
            contracts: DataFrame with contract data
            price_col: Column with unit price
            date_col: Column with date for trend calculation
            
        Returns:
            Tuple of (price_cv, price_trend)
        """
        if contracts.empty or price_col not in contracts.columns:
            return np.nan, np.nan
        
        prices = contracts[price_col].dropna()
        
        if len(prices) < 2:
            return np.nan, np.nan
        
        # Coefficient of variation
        price_cv = prices.std() / prices.mean() if prices.mean() != 0 else np.nan
        
        # Price trend (if dates available)
        if date_col in contracts.columns:
            dated = contracts.dropna(subset=[price_col, date_col]).copy()
            if len(dated) >= 2:
                dated = dated.sort_values(date_col)
                first_half = dated.iloc[:len(dated)//2][price_col].mean()
                second_half = dated.iloc[len(dated)//2:][price_col].mean()
                price_trend = (second_half - first_half) / first_half if first_half != 0 else np.nan
            else:
                price_trend = np.nan
        else:
            price_trend = np.nan
        
        return price_cv, price_trend
    
    def calculate_resilience(
        self,
        contracts: pd.DataFrame,
        entity_id: str,
        period: str,
        **kwargs
    ) -> ResilienceMetrics:
        """
        Calculate full resilience metrics for an entity-period.
        
        Args:
            contracts: DataFrame with contract data for the entity-period
            entity_id: Identifier for the contracting authority
            period: Time period label
            **kwargs: Column name overrides
            
        Returns:
            ResilienceMetrics dataclass
        """
        # Concentration
        hhi, top3_share, n_suppliers = self.calculate_supplier_concentration(
            contracts,
            supplier_col=kwargs.get("supplier_col", "supplier_id"),
            value_col=kwargs.get("value_col", "value_amount")
        )
        
        # Delivery
        on_time, completion, avg_delay = self.calculate_delivery_reliability(
            contracts,
            planned_date_col=kwargs.get("planned_date_col", "planned_delivery_date"),
            actual_date_col=kwargs.get("actual_date_col", "actual_delivery_date"),
            status_col=kwargs.get("status_col", "status")
        )
        
        # Competition
        avg_bidders, single_bid = self.calculate_competition_metrics(
            contracts,
            bidders_col=kwargs.get("bidders_col", "n_bidders")
        )
        
        # Price stability
        price_cv, price_trend = self.calculate_price_stability(
            contracts,
            price_col=kwargs.get("price_col", "unit_price"),
            date_col=kwargs.get("date_col", "award_date")
        )
        
        return ResilienceMetrics(
            entity_id=entity_id,
            period=period,
            on_time_delivery_rate=on_time,
            completion_rate=completion,
            avg_delivery_delay_days=avg_delay,
            n_unique_suppliers=n_suppliers,
            supplier_hhi=hhi,
            top3_supplier_share=top3_share,
            avg_bidders=avg_bidders,
            single_bid_rate=single_bid,
            price_cv=price_cv,
            price_trend=price_trend
        )


# =============================================================================
# GREEN PROCUREMENT METRICS
# =============================================================================

class GreenProcurementClassifier:
    """
    Classify procurement as 'green' based on CPV codes and keywords.
    
    Green public procurement (GPP) is a key sustainability indicator.
    This classifier identifies environmentally-preferable procurement.
    """
    
    # CPV codes associated with green/sustainable procurement
    GREEN_CPV_CODES = {
        # Renewable energy
        "09310000": "solar_energy",
        "09320000": "thermal_energy",
        "09330000": "solar_heat",
        "09331000": "solar_panels",
        "09331100": "solar_collectors",
        "09331200": "photovoltaic_modules",
        "09332000": "solar_installations",
        # Wind energy
        "42112500": "wind_turbines",
        # Electric vehicles
        "34144900": "electric_vehicles",
        "34144910": "electric_buses",
        # Energy efficiency
        "71314000": "energy_consultancy",
        "71314300": "energy_efficiency_consultancy",
        "45261215": "solar_panel_installation",
        # Recycling
        "90500000": "refuse_services",
        "90510000": "refuse_disposal",
        "90511000": "refuse_collection",
        "90512000": "refuse_transport",
        "90513000": "non_hazardous_refuse",
        "90514000": "refuse_recycling",
        # Water treatment
        "90400000": "sewage_services",
        "45232400": "sewer_construction",
        "45232420": "sewage_treatment",
        # Environmental services
        "90700000": "environmental_services",
        "90710000": "environmental_management",
        "90711000": "environmental_assessment",
        "90712000": "environmental_planning",
        "90720000": "environmental_protection",
        "90721000": "environmental_security",
        "90722000": "environment_rehabilitation",
    }
    
    # Keywords indicating green procurement (case-insensitive)
    GREEN_KEYWORDS = [
        "renewable", "solar", "wind", "photovoltaic", "biomass",
        "energy efficient", "energy efficiency", "low carbon",
        "zero emission", "electric vehicle", "hybrid",
        "recycled", "recyclable", "sustainable", "eco-friendly",
        "organic", "biodegradable", "green", "environmental",
        "climate", "carbon neutral", "net zero"
    ]
    
    def __init__(self):
        """Initialize classifier."""
        self.green_cpv_set = set(self.GREEN_CPV_CODES.keys())
        # Compile regex for keyword matching
        import re
        self.keyword_pattern = re.compile(
            r'\b(' + '|'.join(re.escape(kw) for kw in self.GREEN_KEYWORDS) + r')\b',
            re.IGNORECASE
        )
    
    def is_green_cpv(self, cpv_code: str) -> bool:
        """Check if CPV code is green-classified."""
        if not cpv_code:
            return False
        
        # Check exact match
        if cpv_code in self.green_cpv_set:
            return True
        
        # Check if any green CPV is a prefix
        for green_cpv in self.green_cpv_set:
            if cpv_code.startswith(green_cpv[:4]):  # Match at division level
                return True
        
        return False
    
    def is_green_text(self, text: str) -> bool:
        """Check if text contains green keywords."""
        if not text:
            return False
        return bool(self.keyword_pattern.search(text))
    
    def classify(
        self,
        cpv_code: Optional[str] = None,
        title: Optional[str] = None,
        description: Optional[str] = None
    ) -> Tuple[bool, str]:
        """
        Classify a procurement as green or not.
        
        Args:
            cpv_code: CPV classification code
            title: Contract title
            description: Contract description
            
        Returns:
            Tuple of (is_green, reason)
        """
        reasons = []
        
        if cpv_code and self.is_green_cpv(cpv_code):
            category = self.GREEN_CPV_CODES.get(cpv_code, "green_category")
            reasons.append(f"cpv:{category}")
        
        if title and self.is_green_text(title):
            reasons.append("title_keyword")
        
        if description and self.is_green_text(description):
            reasons.append("description_keyword")
        
        is_green = len(reasons) > 0
        reason = ";".join(reasons) if reasons else "not_green"
        
        return is_green, reason
    
    def calculate_green_share(
        self,
        contracts: pd.DataFrame,
        cpv_col: str = "cpv_code",
        title_col: str = "title",
        description_col: str = "description",
        value_col: str = "value_amount"
    ) -> Dict[str, float]:
        """
        Calculate green procurement share for a set of contracts.
        
        Args:
            contracts: DataFrame with contract data
            cpv_col: Column with CPV code
            title_col: Column with title
            description_col: Column with description
            value_col: Column with contract value
            
        Returns:
            Dict with share metrics
        """
        if contracts.empty:
            return {
                "green_count_share": np.nan,
                "green_value_share": np.nan,
                "n_green": 0,
                "n_total": 0
            }
        
        # Classify each contract
        green_flags = []
        for _, row in contracts.iterrows():
            is_green, _ = self.classify(
                cpv_code=row.get(cpv_col),
                title=row.get(title_col),
                description=row.get(description_col)
            )
            green_flags.append(is_green)
        
        contracts = contracts.copy()
        contracts["is_green"] = green_flags
        
        # Count share
        n_green = contracts["is_green"].sum()
        n_total = len(contracts)
        green_count_share = n_green / n_total if n_total > 0 else np.nan
        
        # Value share
        if value_col in contracts.columns:
            total_value = contracts[value_col].sum()
            green_value = contracts.loc[contracts["is_green"], value_col].sum()
            green_value_share = green_value / total_value if total_value > 0 else np.nan
        else:
            green_value_share = np.nan
        
        return {
            "green_count_share": green_count_share,
            "green_value_share": green_value_share,
            "n_green": int(n_green),
            "n_total": n_total
        }


# =============================================================================
# COMPOSITE SUSTAINABILITY INDEX
# =============================================================================

@dataclass
class SustainabilityScore:
    """
    Composite sustainability score for procurement entity/period.
    
    Combines carbon footprint, resilience, and green procurement metrics
    into interpretable sustainability indicators.
    """
    entity_id: str
    period: str
    
    # Component scores (0-100, higher = more sustainable)
    carbon_score: float  # Based on footprint intensity vs benchmark
    resilience_score: float  # Based on concentration, competition, reliability
    green_score: float  # Based on green procurement share
    
    # Composite
    composite_score: float  # Weighted average
    
    # Weights used
    weights: Dict[str, float] = field(default_factory=dict)
    
    # Raw metrics for transparency
    raw_metrics: Dict[str, float] = field(default_factory=dict)


class SustainabilityIndexCalculator:
    """
    Calculate composite sustainability index from component metrics.
    
    This index enables:
    1. Cross-entity comparison
    2. Time-series analysis
    3. Causal analysis of rule changes on sustainability
    """
    
    # Default weights
    DEFAULT_WEIGHTS = {
        "carbon": 0.40,
        "resilience": 0.35,
        "green": 0.25
    }
    
    # Benchmarks for scoring (higher intensity = lower score)
    CARBON_BENCHMARKS = {
        "excellent": 0.15,  # kg CO2e/$ below this = 100
        "poor": 0.80,  # kg CO2e/$ above this = 0
    }
    
    def __init__(
        self,
        weights: Optional[Dict[str, float]] = None,
        carbon_benchmarks: Optional[Dict[str, float]] = None
    ):
        """
        Initialize calculator.
        
        Args:
            weights: Optional custom weights for components
            carbon_benchmarks: Optional custom benchmarks for carbon scoring
        """
        self.weights = weights or self.DEFAULT_WEIGHTS.copy()
        self.carbon_benchmarks = carbon_benchmarks or self.CARBON_BENCHMARKS.copy()
        
        # Normalize weights
        total = sum(self.weights.values())
        self.weights = {k: v/total for k, v in self.weights.items()}
    
    def score_carbon(self, intensity: float) -> float:
        """
        Convert carbon intensity to 0-100 score.
        
        Lower intensity = higher score (more sustainable).
        """
        if np.isnan(intensity):
            return np.nan
        
        excellent = self.carbon_benchmarks["excellent"]
        poor = self.carbon_benchmarks["poor"]
        
        if intensity <= excellent:
            return 100.0
        elif intensity >= poor:
            return 0.0
        else:
            # Linear interpolation
            return 100 * (poor - intensity) / (poor - excellent)
    
    def score_resilience(self, metrics: ResilienceMetrics) -> float:
        """
        Convert resilience metrics to 0-100 score.
        
        Higher competition, lower concentration, better reliability = higher score.
        """
        scores = []
        
        # HHI scoring (lower = better, 0-10000 scale)
        if metrics.supplier_hhi is not None and not np.isnan(metrics.supplier_hhi):
            # HHI < 1500 = unconcentrated, > 2500 = highly concentrated
            hhi_score = max(0, min(100, 100 * (2500 - metrics.supplier_hhi) / 2500))
            scores.append(hhi_score)
        
        # Single-bid rate scoring (lower = better)
        if metrics.single_bid_rate is not None and not np.isnan(metrics.single_bid_rate):
            single_score = 100 * (1 - metrics.single_bid_rate)
            scores.append(single_score)
        
        # Completion rate scoring (higher = better)
        if metrics.completion_rate is not None and not np.isnan(metrics.completion_rate):
            completion_score = 100 * metrics.completion_rate
            scores.append(completion_score)
        
        # Average bidders scoring (higher = better, capped at 5)
        if metrics.avg_bidders is not None and not np.isnan(metrics.avg_bidders):
            bidder_score = min(100, 100 * metrics.avg_bidders / 5)
            scores.append(bidder_score)
        
        return np.mean(scores) if scores else np.nan
    
    def score_green(self, green_share: float) -> float:
        """
        Convert green share to 0-100 score.
        
        Higher green share = higher score.
        """
        if np.isnan(green_share):
            return np.nan
        return 100 * min(1.0, green_share / 0.30)  # 30% green = 100 score
    
    def calculate(
        self,
        entity_id: str,
        period: str,
        carbon_intensity: float,
        resilience_metrics: ResilienceMetrics,
        green_share: float
    ) -> SustainabilityScore:
        """
        Calculate composite sustainability score.
        
        Args:
            entity_id: Entity identifier
            period: Time period
            carbon_intensity: Average kg CO2e per USD
            resilience_metrics: ResilienceMetrics instance
            green_share: Proportion of green procurement
            
        Returns:
            SustainabilityScore with component and composite scores
        """
        carbon_score = self.score_carbon(carbon_intensity)
        resilience_score = self.score_resilience(resilience_metrics)
        green_score = self.score_green(green_share)
        
        # Composite (weighted average, handling missing)
        components = {
            "carbon": carbon_score,
            "resilience": resilience_score,
            "green": green_score
        }
        
        valid_components = {k: v for k, v in components.items() if not np.isnan(v)}
        
        if valid_components:
            # Renormalize weights for available components
            valid_weights = {k: self.weights[k] for k in valid_components}
            total_weight = sum(valid_weights.values())
            
            composite = sum(
                valid_components[k] * valid_weights[k] / total_weight
                for k in valid_components
            )
        else:
            composite = np.nan
        
        return SustainabilityScore(
            entity_id=entity_id,
            period=period,
            carbon_score=carbon_score,
            resilience_score=resilience_score,
            green_score=green_score,
            composite_score=composite,
            weights=self.weights.copy(),
            raw_metrics={
                "carbon_intensity": carbon_intensity,
                "supplier_hhi": resilience_metrics.supplier_hhi,
                "single_bid_rate": resilience_metrics.single_bid_rate,
                "green_share": green_share
            }
        )


# =============================================================================
# INTEGRATION: FULL SUSTAINABILITY PIPELINE
# =============================================================================

class SustainabilityPipeline:
    """
    End-to-end pipeline for computing sustainability metrics from procurement data.
    
    Integrates:
    - Carbon footprinting
    - Resilience metrics
    - Green procurement classification
    - Composite index calculation
    
    Designed for causal analysis: compute metrics before/after threshold
    or reform to estimate sustainability effects.
    """
    
    def __init__(
        self,
        carbon_calculator: Optional[CarbonFootprintCalculator] = None,
        resilience_calculator: Optional[ResilienceCalculator] = None,
        green_classifier: Optional[GreenProcurementClassifier] = None,
        index_calculator: Optional[SustainabilityIndexCalculator] = None
    ):
        """
        Initialize pipeline with component calculators.
        """
        self.carbon = carbon_calculator or CarbonFootprintCalculator()
        self.resilience = resilience_calculator or ResilienceCalculator()
        self.green = green_classifier or GreenProcurementClassifier()
        self.index = index_calculator or SustainabilityIndexCalculator()
    
    def process(
        self,
        contracts: pd.DataFrame,
        entity_id: str,
        period: str,
        column_mapping: Optional[Dict[str, str]] = None
    ) -> Dict:
        """
        Process contracts to compute all sustainability metrics.
        
        Args:
            contracts: DataFrame with procurement data
            entity_id: Entity identifier
            period: Time period label
            column_mapping: Optional mapping of standard names to actual columns
            
        Returns:
            Dict with all metrics and scores
        """
        cols = column_mapping or {}
        
        # Carbon footprint
        footprints = self.carbon.calculate_batch(
            contracts,
            id_col=cols.get("id", "ocid"),
            value_col=cols.get("value", "value_amount"),
            cpv_col=cols.get("cpv", "cpv_code")
        )
        
        if not footprints.empty:
            total_spend = footprints["spend_usd"].sum()
            total_footprint = footprints["footprint_kg_co2e"].sum()
            carbon_intensity = total_footprint / total_spend if total_spend > 0 else np.nan
        else:
            carbon_intensity = np.nan
        
        # Resilience
        resilience_metrics = self.resilience.calculate_resilience(
            contracts, entity_id, period,
            supplier_col=cols.get("supplier", "supplier_id"),
            value_col=cols.get("value", "value_amount"),
            bidders_col=cols.get("bidders", "n_bidders")
        )
        
        # Green share
        green_result = self.green.calculate_green_share(
            contracts,
            cpv_col=cols.get("cpv", "cpv_code"),
            title_col=cols.get("title", "title"),
            description_col=cols.get("description", "description"),
            value_col=cols.get("value", "value_amount")
        )
        
        # Composite index
        sustainability_score = self.index.calculate(
            entity_id=entity_id,
            period=period,
            carbon_intensity=carbon_intensity,
            resilience_metrics=resilience_metrics,
            green_share=green_result["green_value_share"]
        )
        
        return {
            "entity_id": entity_id,
            "period": period,
            "carbon": {
                "intensity_kg_per_usd": carbon_intensity,
                "total_footprint_kg": total_footprint if not footprints.empty else np.nan,
                "total_spend_usd": total_spend if not footprints.empty else np.nan,
            },
            "resilience": resilience_metrics.to_dict(),
            "green": green_result,
            "sustainability_score": {
                "carbon_score": sustainability_score.carbon_score,
                "resilience_score": sustainability_score.resilience_score,
                "green_score": sustainability_score.green_score,
                "composite_score": sustainability_score.composite_score,
            },
            "n_contracts": len(contracts)
        }
    
    def compare_periods(
        self,
        contracts: pd.DataFrame,
        entity_id: str,
        period_col: str,
        pre_period: str,
        post_period: str,
        column_mapping: Optional[Dict[str, str]] = None
    ) -> Dict:
        """
        Compare sustainability metrics between two periods.
        
        Useful for DiD-style analysis of reform effects.
        
        Args:
            contracts: DataFrame with procurement data
            entity_id: Entity identifier
            period_col: Column with period labels
            pre_period: Pre-treatment period label
            post_period: Post-treatment period label
            column_mapping: Optional column mapping
            
        Returns:
            Dict with pre, post, and difference metrics
        """
        pre_contracts = contracts[contracts[period_col] == pre_period]
        post_contracts = contracts[contracts[period_col] == post_period]
        
        pre_metrics = self.process(pre_contracts, entity_id, pre_period, column_mapping)
        post_metrics = self.process(post_contracts, entity_id, post_period, column_mapping)
        
        # Calculate differences
        def safe_diff(post_val, pre_val):
            if pd.isna(post_val) or pd.isna(pre_val):
                return np.nan
            return post_val - pre_val
        
        diff = {
            "carbon_intensity_change": safe_diff(
                post_metrics["carbon"]["intensity_kg_per_usd"],
                pre_metrics["carbon"]["intensity_kg_per_usd"]
            ),
            "sustainability_score_change": safe_diff(
                post_metrics["sustainability_score"]["composite_score"],
                pre_metrics["sustainability_score"]["composite_score"]
            ),
            "single_bid_rate_change": safe_diff(
                post_metrics["resilience"]["single_bid_rate"],
                pre_metrics["resilience"]["single_bid_rate"]
            ),
            "green_share_change": safe_diff(
                post_metrics["green"]["green_value_share"],
                pre_metrics["green"]["green_value_share"]
            ),
        }
        
        return {
            "entity_id": entity_id,
            "pre_period": pre_period,
            "post_period": post_period,
            "pre": pre_metrics,
            "post": post_metrics,
            "difference": diff
        }


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    # Example usage with synthetic data
    
    # Create sample contracts
    sample_contracts = pd.DataFrame({
        "ocid": [f"contract_{i}" for i in range(100)],
        "value_amount": np.random.lognormal(10, 1, 100),
        "cpv_code": np.random.choice(
            ["45000000", "72000000", "33000000", "09310000", "90514000"],
            100
        ),
        "supplier_id": np.random.choice([f"supplier_{i}" for i in range(20)], 100),
        "n_bidders": np.random.poisson(3, 100) + 1,
        "title": [f"Contract for goods/services {i}" for i in range(100)],
        "description": ["Standard procurement description"] * 100,
        "status": np.random.choice(["completed", "active", "cancelled"], 100, p=[0.8, 0.15, 0.05]),
        "period": np.random.choice(["pre", "post"], 100)
    })
    
    # Add some green keywords to random contracts
    for i in np.random.choice(100, 15, replace=False):
        sample_contracts.loc[i, "title"] = "Solar panel installation and energy efficiency"
    
    # Initialize pipeline
    pipeline = SustainabilityPipeline()
    
    # Process contracts
    results = pipeline.process(
        sample_contracts,
        entity_id="test_authority",
        period="2024"
    )
    
    print("\n=== Sustainability Metrics ===")
    print(f"Entity: {results['entity_id']}")
    print(f"Period: {results['period']}")
    print(f"Contracts: {results['n_contracts']}")
    print(f"\nCarbon Intensity: {results['carbon']['intensity_kg_per_usd']:.3f} kg CO2e/$")
    print(f"Single-Bid Rate: {results['resilience']['single_bid_rate']:.1%}")
    print(f"Green Share: {results['green']['green_value_share']:.1%}")
    print(f"\nComposite Sustainability Score: {results['sustainability_score']['composite_score']:.1f}/100")
    
    # Compare periods
    comparison = pipeline.compare_periods(
        sample_contracts,
        entity_id="test_authority",
        period_col="period",
        pre_period="pre",
        post_period="post"
    )
    
    print("\n=== Period Comparison ===")
    print(f"Sustainability Score Change: {comparison['difference']['sustainability_score_change']:+.1f}")
    print(f"Carbon Intensity Change: {comparison['difference']['carbon_intensity_change']:+.3f} kg CO2e/$")
