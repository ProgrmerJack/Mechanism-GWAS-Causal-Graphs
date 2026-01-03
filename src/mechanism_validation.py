"""
Mechanism Index Validation Framework for E-Procurement Transparency Analysis.

This module implements a rigorous validation framework for the Mechanism Index,
ensuring it reliably measures procurement transparency mechanisms across:
- Multiple countries (Ukraine, Colombia, UK)
- Multiple languages (Ukrainian, Spanish, English)
- Multiple procurement stages

Key Components:
1. Labeled dataset schema and management
2. Cross-lingual validation strategy
3. Inter-annotator agreement metrics
4. Construct validity testing
5. External validation against known benchmarks

References:
    Krippendorff, K. (2018). Content analysis: An introduction to its methodology.
    Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters.
    Cohen, J. (1960). A coefficient of agreement for nominal scales.
    Cronbach, L.J. (1951). Coefficient alpha and the internal structure of tests.

Author: Abduxoliq Ashuraliyev
Affiliation: Independent Researcher, Tashkent, Uzbekistan
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


# ==============================================================================
# Enumerations and Constants
# ==============================================================================


class Language(str, Enum):
    """Supported languages for procurement documents."""

    ENGLISH = "en"
    UKRAINIAN = "uk"
    SPANISH = "es"


class MechanismDimension(str, Enum):
    """Dimensions of procurement transparency measured by the index."""

    DISCLOSURE = "disclosure"  # Information availability
    PARTICIPATION = "participation"  # Bidder access and engagement
    OVERSIGHT = "oversight"  # Monitoring and accountability
    INTEGRITY = "integrity"  # Anti-corruption safeguards
    EFFICIENCY = "efficiency"  # Process optimization


class AnnotationStatus(str, Enum):
    """Status of annotation task."""

    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    REVIEWED = "reviewed"
    DISPUTED = "disputed"


# Component structure for Mechanism Index
MECHANISM_COMPONENTS = {
    MechanismDimension.DISCLOSURE: [
        "tender_notice_published",
        "award_notice_published",
        "contract_published",
        "bidder_info_disclosed",
        "evaluation_criteria_published",
        "price_breakdown_available",
    ],
    MechanismDimension.PARTICIPATION: [
        "electronic_submission",
        "minimum_bidder_threshold",
        "complaint_mechanism",
        "pre_bid_clarification",
        "extension_allowed",
    ],
    MechanismDimension.OVERSIGHT: [
        "audit_trail_complete",
        "amendment_justified",
        "price_reasonableness_check",
        "conflict_of_interest_check",
    ],
    MechanismDimension.INTEGRITY: [
        "debarment_check",
        "beneficial_ownership_disclosed",
        "red_flag_monitoring",
        "whistleblower_protection",
    ],
    MechanismDimension.EFFICIENCY: [
        "e_procurement_used",
        "framework_agreement",
        "centralized_purchasing",
        "timeline_compliance",
    ],
}


# ==============================================================================
# Data Classes
# ==============================================================================


@dataclass
class AnnotationItem:
    """Single item to be annotated."""

    item_id: str
    contract_id: str
    country: str
    language: Language
    text_content: str
    component: str
    dimension: MechanismDimension

    # Metadata
    source_url: str = ""
    extraction_date: str = ""

    # Ground truth (if available)
    gold_label: bool | None = None
    gold_confidence: float | None = None


@dataclass
class AnnotatorResponse:
    """Response from a single annotator."""

    annotator_id: str
    item_id: str
    timestamp: datetime
    label: bool  # True = mechanism present, False = absent
    confidence: float  # 0-1 scale
    reasoning: str = ""
    time_spent_seconds: int = 0


@dataclass
class ValidationDataset:
    """Container for validation dataset."""

    name: str
    version: str
    created_date: datetime
    items: list[AnnotationItem]
    responses: list[AnnotatorResponse] = field(default_factory=list)

    # Metadata
    n_annotators: int = 0
    countries: list[str] = field(default_factory=list)
    languages: list[Language] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "name": self.name,
            "version": self.version,
            "created_date": self.created_date.isoformat(),
            "n_items": len(self.items),
            "n_responses": len(self.responses),
            "n_annotators": self.n_annotators,
            "countries": self.countries,
            "languages": [lang.value for lang in self.languages],
        }


@dataclass
class AgreementMetrics:
    """Inter-annotator agreement metrics."""

    cohens_kappa: float | None = None  # Pairwise agreement
    fleiss_kappa: float | None = None  # Multi-rater agreement
    krippendorff_alpha: float | None = None  # Reliability coefficient
    percent_agreement: float | None = None  # Raw agreement
    gwet_ac1: float | None = None  # Chance-corrected agreement

    # Confidence intervals
    kappa_ci_lower: float | None = None
    kappa_ci_upper: float | None = None

    # Sample info
    n_items: int = 0
    n_annotators: int = 0


@dataclass
class ValidityMetrics:
    """Construct and criterion validity metrics."""

    # Internal consistency
    cronbachs_alpha: float | None = None
    split_half_reliability: float | None = None

    # Construct validity
    convergent_correlation: float | None = None  # With related measures
    discriminant_correlation: float | None = None  # With unrelated measures

    # Criterion validity
    predictive_r2: float | None = None  # Against outcome
    concurrent_correlation: float | None = None  # With existing measure

    # Cross-lingual
    cross_lingual_correlation: float | None = None
    translation_equivalence: float | None = None


@dataclass
class ValidationReport:
    """Complete validation report."""

    dataset: ValidationDataset
    agreement_metrics: AgreementMetrics
    validity_metrics: ValidityMetrics
    dimension_results: dict[str, AgreementMetrics]
    country_results: dict[str, AgreementMetrics]
    language_results: dict[str, AgreementMetrics]

    # Quality assessment
    overall_quality: str  # "excellent", "good", "acceptable", "poor"
    concerns: list[str]
    recommendations: list[str]


# ==============================================================================
# Validation Dataset Manager
# ==============================================================================


class ValidationDatasetManager:
    """
    Manages the creation and maintenance of validation datasets.

    Implements stratified sampling across:
    - Countries
    - Languages
    - Mechanism dimensions
    - Procurement stages
    """

    def __init__(
        self,
        output_dir: Path,
        seed: int = 42,
    ):
        """
        Initialize dataset manager.

        Args:
            output_dir: Directory for dataset storage
            seed: Random seed for reproducibility
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.rng = np.random.default_rng(seed)

        logger.info(f"Validation dataset manager initialized: {output_dir}")

    def create_stratified_sample(
        self,
        source_data: pd.DataFrame,
        n_per_stratum: int = 50,
        strata_cols: list[str] | None = None,
    ) -> list[AnnotationItem]:
        """
        Create stratified sample from source data.

        Args:
            source_data: DataFrame with procurement records
            n_per_stratum: Number of items per stratum
            strata_cols: Columns to stratify by

        Returns:
            List of AnnotationItem for annotation
        """
        if strata_cols is None:
            strata_cols = ["country", "dimension"]

        items = []

        for _, stratum_data in source_data.groupby(strata_cols):
            n_sample = min(n_per_stratum, len(stratum_data))
            sampled = stratum_data.sample(n=n_sample, random_state=self.rng)

            for _, row in sampled.iterrows():
                item_id = self._generate_item_id(row)
                items.append(self._row_to_item(row, item_id))

        logger.info(f"Created stratified sample: {len(items)} items")
        return items

    def _generate_item_id(self, row: pd.Series) -> str:
        """Generate unique item ID from row content."""
        content = f"{row.get('contract_id', '')}{row.get('text_content', '')}"
        return hashlib.md5(content.encode()).hexdigest()[:12]

    def _row_to_item(self, row: pd.Series, item_id: str) -> AnnotationItem:
        """Convert DataFrame row to AnnotationItem."""
        return AnnotationItem(
            item_id=item_id,
            contract_id=str(row.get("contract_id", "")),
            country=str(row.get("country", "")),
            language=Language(row.get("language", "en")),
            text_content=str(row.get("text_content", "")),
            component=str(row.get("component", "")),
            dimension=MechanismDimension(row.get("dimension", "disclosure")),
            source_url=str(row.get("source_url", "")),
            extraction_date=str(row.get("extraction_date", "")),
        )

    def create_gold_standard_set(
        self,
        items: list[AnnotationItem],
        expert_labels: dict[str, tuple[bool, float]],
    ) -> list[AnnotationItem]:
        """
        Create gold standard set with expert labels.

        Args:
            items: Items to label
            expert_labels: Dict mapping item_id to (label, confidence)

        Returns:
            Items with gold labels
        """
        for item in items:
            if item.item_id in expert_labels:
                label, conf = expert_labels[item.item_id]
                item.gold_label = label
                item.gold_confidence = conf

        labeled = [i for i in items if i.gold_label is not None]
        logger.info(f"Created gold standard: {len(labeled)} labeled items")
        return items

    def export_annotation_format(
        self,
        items: list[AnnotationItem],
        output_path: Path,
        format_type: str = "csv",
    ) -> None:
        """
        Export items in format suitable for annotation tools.

        Args:
            items: Items to export
            output_path: Output file path
            format_type: "csv", "json", or "labelstudio"
        """
        if format_type == "csv":
            records = []
            for item in items:
                records.append({
                    "item_id": item.item_id,
                    "contract_id": item.contract_id,
                    "country": item.country,
                    "language": item.language.value,
                    "component": item.component,
                    "dimension": item.dimension.value,
                    "text_content": item.text_content[:500],  # Truncate for display
                    "source_url": item.source_url,
                })
            pd.DataFrame(records).to_csv(output_path, index=False)

        elif format_type == "json":
            records = [
                {
                    "id": item.item_id,
                    "data": {
                        "text": item.text_content,
                        "contract_id": item.contract_id,
                        "country": item.country,
                        "language": item.language.value,
                        "component": item.component,
                        "dimension": item.dimension.value,
                    },
                }
                for item in items
            ]
            output_path.write_text(json.dumps(records, indent=2, ensure_ascii=False))

        elif format_type == "labelstudio":
            # Label Studio import format
            tasks = []
            for item in items:
                tasks.append({
                    "id": item.item_id,
                    "data": {
                        "text": item.text_content,
                        "meta": {
                            "contract_id": item.contract_id,
                            "country": item.country,
                            "language": item.language.value,
                        },
                    },
                    "predictions": [],
                })
            output_path.write_text(json.dumps(tasks, indent=2, ensure_ascii=False))

        logger.info(f"Exported {len(items)} items to {output_path}")

    def import_annotations(
        self,
        annotation_path: Path,
        format_type: str = "csv",
    ) -> list[AnnotatorResponse]:
        """
        Import annotations from annotation tool export.

        Args:
            annotation_path: Path to annotation file
            format_type: "csv" or "json"

        Returns:
            List of AnnotatorResponse
        """
        responses = []

        if format_type == "csv":
            df = pd.read_csv(annotation_path)
            for _, row in df.iterrows():
                responses.append(AnnotatorResponse(
                    annotator_id=str(row["annotator_id"]),
                    item_id=str(row["item_id"]),
                    timestamp=pd.to_datetime(row.get("timestamp", datetime.now())),
                    label=bool(row["label"]),
                    confidence=float(row.get("confidence", 1.0)),
                    reasoning=str(row.get("reasoning", "")),
                    time_spent_seconds=int(row.get("time_spent", 0)),
                ))

        elif format_type == "json":
            data = json.loads(annotation_path.read_text())
            for item in data:
                for annotation in item.get("annotations", []):
                    responses.append(AnnotatorResponse(
                        annotator_id=str(annotation["completed_by"]),
                        item_id=str(item["id"]),
                        timestamp=datetime.fromisoformat(annotation["created_at"]),
                        label=self._parse_label(annotation["result"]),
                        confidence=float(annotation.get("confidence", 1.0)),
                    ))

        logger.info(f"Imported {len(responses)} annotations")
        return responses

    def _parse_label(self, result: list[dict]) -> bool:
        """Parse label from annotation result."""
        for r in result:
            if r.get("type") == "choices":
                choices = r.get("value", {}).get("choices", [])
                return "present" in [c.lower() for c in choices]
        return False


# ==============================================================================
# Agreement Calculator
# ==============================================================================


class InterAnnotatorAgreement:
    """
    Calculates inter-annotator agreement metrics.

    Implements:
    - Cohen's Kappa (pairwise)
    - Fleiss' Kappa (multi-rater)
    - Krippendorff's Alpha (any scale)
    - Gwet's AC1 (paradox-resistant)
    """

    @staticmethod
    def cohens_kappa(
        rater1: np.ndarray,
        rater2: np.ndarray,
    ) -> tuple[float, float, float]:
        """
        Calculate Cohen's Kappa between two raters.

        Args:
            rater1: Binary ratings from rater 1
            rater2: Binary ratings from rater 2

        Returns:
            Tuple of (kappa, ci_lower, ci_upper)
        """
        # Confusion matrix
        n = len(rater1)
        a = np.sum((rater1 == 1) & (rater2 == 1))  # Both agree positive
        b = np.sum((rater1 == 1) & (rater2 == 0))  # Disagree
        c = np.sum((rater1 == 0) & (rater2 == 1))  # Disagree
        d = np.sum((rater1 == 0) & (rater2 == 0))  # Both agree negative

        # Observed agreement
        po = (a + d) / n

        # Expected agreement
        pe = ((a + b) * (a + c) + (c + d) * (b + d)) / (n * n)

        # Kappa
        if pe == 1:
            kappa = 1.0
        else:
            kappa = (po - pe) / (1 - pe)

        # Standard error (Fleiss et al., 1969)
        se = np.sqrt(po * (1 - po) / (n * (1 - pe) ** 2))

        ci_lower = kappa - 1.96 * se
        ci_upper = kappa + 1.96 * se

        return kappa, ci_lower, ci_upper

    @staticmethod
    def fleiss_kappa(
        ratings_matrix: np.ndarray,
    ) -> float:
        """
        Calculate Fleiss' Kappa for multiple raters.

        Args:
            ratings_matrix: n_items x n_categories matrix
                           Each row sums to n_raters

        Returns:
            Fleiss' Kappa value
        """
        n_items, n_categories = ratings_matrix.shape
        n_raters = int(ratings_matrix.sum(axis=1).mean())

        # Proportion of ratings in each category
        p = ratings_matrix.sum(axis=0) / (n_items * n_raters)

        # P_i: agreement for each item
        p_i = (ratings_matrix ** 2).sum(axis=1) - n_raters
        p_i = p_i / (n_raters * (n_raters - 1))

        # P-bar: mean of P_i
        p_bar = p_i.mean()

        # P-bar_e: expected agreement by chance
        p_bar_e = (p ** 2).sum()

        # Fleiss' Kappa
        if p_bar_e == 1:
            return 1.0
        return (p_bar - p_bar_e) / (1 - p_bar_e)

    @staticmethod
    def krippendorff_alpha(
        reliability_data: np.ndarray,
        level: str = "nominal",
    ) -> float:
        """
        Calculate Krippendorff's Alpha.

        Args:
            reliability_data: n_raters x n_items matrix (NaN for missing)
            level: "nominal", "ordinal", "interval", or "ratio"

        Returns:
            Alpha coefficient
        """
        # Remove items with fewer than 2 ratings
        n_ratings = (~np.isnan(reliability_data)).sum(axis=0)
        valid_items = n_ratings >= 2
        data = reliability_data[:, valid_items]

        n_raters, n_items = data.shape

        # Get all pairs of values
        values = data[~np.isnan(data)].flatten()
        unique_values = np.unique(values)

        # Observed disagreement
        d_o = 0.0
        n_pairs = 0

        for i in range(n_items):
            item_values = data[:, i][~np.isnan(data[:, i])]
            m = len(item_values)
            if m < 2:
                continue

            for v1 in item_values:
                for v2 in item_values:
                    if v1 != v2:
                        d_o += InterAnnotatorAgreement._difference(v1, v2, level)
                        n_pairs += 1

        d_o = d_o / n_pairs if n_pairs > 0 else 0

        # Expected disagreement
        d_e = 0.0
        n_total = len(values)

        for v1 in unique_values:
            for v2 in unique_values:
                n1 = np.sum(values == v1)
                n2 = np.sum(values == v2)
                d_e += n1 * n2 * InterAnnotatorAgreement._difference(v1, v2, level)

        d_e = d_e / (n_total * (n_total - 1))

        # Alpha
        if d_e == 0:
            return 1.0
        return 1 - d_o / d_e

    @staticmethod
    def _difference(v1: float, v2: float, level: str) -> float:
        """Calculate difference metric based on measurement level."""
        if level == "nominal":
            return 0 if v1 == v2 else 1
        elif level == "ordinal":
            return abs(v1 - v2)  # Simplified
        elif level == "interval":
            return (v1 - v2) ** 2
        elif level == "ratio":
            return ((v1 - v2) / (v1 + v2)) ** 2 if (v1 + v2) != 0 else 0
        return 0

    @staticmethod
    def gwet_ac1(
        rater1: np.ndarray,
        rater2: np.ndarray,
    ) -> float:
        """
        Calculate Gwet's AC1 (first-order agreement coefficient).

        More robust than Kappa when marginal distributions are unbalanced.

        Args:
            rater1: Binary ratings from rater 1
            rater2: Binary ratings from rater 2

        Returns:
            AC1 coefficient
        """
        n = len(rater1)

        # Observed agreement
        po = np.mean(rater1 == rater2)

        # Marginal probabilities
        p1 = (np.mean(rater1 == 1) + np.mean(rater2 == 1)) / 2
        p0 = 1 - p1

        # Expected agreement by chance (AC1)
        pe = 2 * p1 * p0

        # AC1
        if pe == 1:
            return 1.0
        return (po - pe) / (1 - pe)

    @classmethod
    def compute_all_metrics(
        cls,
        responses: list[AnnotatorResponse],
        items: list[AnnotationItem],
    ) -> AgreementMetrics:
        """
        Compute all agreement metrics from responses.

        Args:
            responses: List of annotator responses
            items: List of annotation items

        Returns:
            AgreementMetrics with all coefficients
        """
        # Build response matrix
        item_ids = [item.item_id for item in items]
        annotator_ids = list({r.annotator_id for r in responses})

        n_items = len(item_ids)
        n_annotators = len(annotator_ids)

        if n_annotators < 2:
            logger.warning("Need at least 2 annotators for agreement metrics")
            return AgreementMetrics(n_items=n_items, n_annotators=n_annotators)

        # Create ratings matrix
        ratings = np.full((n_annotators, n_items), np.nan)

        for resp in responses:
            try:
                i = annotator_ids.index(resp.annotator_id)
                j = item_ids.index(resp.item_id)
                ratings[i, j] = 1.0 if resp.label else 0.0
            except ValueError:
                continue

        # Pairwise kappa (average across pairs)
        kappas = []
        for i in range(n_annotators):
            for j in range(i + 1, n_annotators):
                mask = ~np.isnan(ratings[i]) & ~np.isnan(ratings[j])
                if mask.sum() >= 10:
                    k, _, _ = cls.cohens_kappa(
                        ratings[i, mask].astype(int),
                        ratings[j, mask].astype(int),
                    )
                    kappas.append(k)

        avg_kappa = np.mean(kappas) if kappas else None

        # Fleiss Kappa (convert to category counts)
        # Count how many raters said 0 vs 1 for each item
        fleiss_matrix = np.zeros((n_items, 2))
        for j in range(n_items):
            item_ratings = ratings[:, j]
            valid = ~np.isnan(item_ratings)
            fleiss_matrix[j, 0] = np.sum(item_ratings[valid] == 0)
            fleiss_matrix[j, 1] = np.sum(item_ratings[valid] == 1)

        # Only include items with at least 2 ratings
        valid_items = fleiss_matrix.sum(axis=1) >= 2
        if valid_items.sum() >= 10:
            fleiss = cls.fleiss_kappa(fleiss_matrix[valid_items])
        else:
            fleiss = None

        # Krippendorff's Alpha
        if n_items >= 10:
            alpha = cls.krippendorff_alpha(ratings, level="nominal")
        else:
            alpha = None

        # Percent agreement
        n_agree = 0
        n_total = 0
        for j in range(n_items):
            item_ratings = ratings[:, j]
            valid = ~np.isnan(item_ratings)
            if valid.sum() >= 2:
                unique_vals = np.unique(item_ratings[valid])
                if len(unique_vals) == 1:
                    n_agree += 1
                n_total += 1

        pct_agree = n_agree / n_total if n_total > 0 else None

        # Gwet's AC1 (pairwise average)
        ac1s = []
        for i in range(n_annotators):
            for j in range(i + 1, n_annotators):
                mask = ~np.isnan(ratings[i]) & ~np.isnan(ratings[j])
                if mask.sum() >= 10:
                    ac1 = cls.gwet_ac1(
                        ratings[i, mask].astype(int),
                        ratings[j, mask].astype(int),
                    )
                    ac1s.append(ac1)

        avg_ac1 = np.mean(ac1s) if ac1s else None

        return AgreementMetrics(
            cohens_kappa=avg_kappa,
            fleiss_kappa=fleiss,
            krippendorff_alpha=alpha,
            percent_agreement=pct_agree,
            gwet_ac1=avg_ac1,
            n_items=n_items,
            n_annotators=n_annotators,
        )


# ==============================================================================
# Validity Assessment
# ==============================================================================


class ValidityAssessment:
    """
    Assesses construct and criterion validity of the Mechanism Index.
    """

    @staticmethod
    def cronbachs_alpha(items_matrix: np.ndarray) -> float:
        """
        Calculate Cronbach's Alpha for internal consistency.

        Args:
            items_matrix: n_observations x n_components matrix

        Returns:
            Alpha coefficient
        """
        n_items = items_matrix.shape[1]

        # Item variances
        item_vars = items_matrix.var(axis=0, ddof=1)

        # Total score variance
        total_scores = items_matrix.sum(axis=1)
        total_var = total_scores.var(ddof=1)

        if total_var == 0:
            return np.nan

        alpha = (n_items / (n_items - 1)) * (1 - item_vars.sum() / total_var)
        return alpha

    @staticmethod
    def split_half_reliability(
        items_matrix: np.ndarray,
        n_iterations: int = 100,
    ) -> tuple[float, float]:
        """
        Calculate split-half reliability with Spearman-Brown correction.

        Args:
            items_matrix: n_observations x n_components matrix
            n_iterations: Number of random splits

        Returns:
            Tuple of (mean reliability, std)
        """
        n_items = items_matrix.shape[1]
        reliabilities = []

        rng = np.random.default_rng(42)

        for _ in range(n_iterations):
            # Random split
            idx = rng.permutation(n_items)
            half1 = idx[: n_items // 2]
            half2 = idx[n_items // 2:]

            # Scores for each half
            score1 = items_matrix[:, half1].sum(axis=1)
            score2 = items_matrix[:, half2].sum(axis=1)

            # Correlation
            if score1.std() > 0 and score2.std() > 0:
                r = np.corrcoef(score1, score2)[0, 1]

                # Spearman-Brown correction
                reliability = 2 * r / (1 + r)
                reliabilities.append(reliability)

        return np.mean(reliabilities), np.std(reliabilities)

    @staticmethod
    def convergent_validity(
        index_scores: np.ndarray,
        related_measure: np.ndarray,
    ) -> tuple[float, float]:
        """
        Assess convergent validity with related measure.

        Args:
            index_scores: Mechanism Index scores
            related_measure: Theoretically related measure (e.g., TI CPI)

        Returns:
            Tuple of (correlation, p-value)
        """
        mask = ~np.isnan(index_scores) & ~np.isnan(related_measure)
        if mask.sum() < 10:
            return np.nan, np.nan

        r, p = stats.pearsonr(index_scores[mask], related_measure[mask])
        return r, p

    @staticmethod
    def discriminant_validity(
        index_scores: np.ndarray,
        unrelated_measure: np.ndarray,
    ) -> tuple[float, float]:
        """
        Assess discriminant validity with unrelated measure.

        Args:
            index_scores: Mechanism Index scores
            unrelated_measure: Theoretically unrelated measure

        Returns:
            Tuple of (correlation, p-value)
        """
        # Same computation, different interpretation
        return ValidityAssessment.convergent_validity(index_scores, unrelated_measure)

    @staticmethod
    def cross_lingual_equivalence(
        scores_lang1: np.ndarray,
        scores_lang2: np.ndarray,
        matched_items: np.ndarray,  # Boolean mask for matched items
    ) -> tuple[float, float, float]:
        """
        Assess cross-lingual measurement equivalence.

        Args:
            scores_lang1: Scores from language 1
            scores_lang2: Scores from language 2 (same items translated)
            matched_items: Boolean mask for paired items

        Returns:
            Tuple of (correlation, mean_diff, max_diff)
        """
        s1 = scores_lang1[matched_items]
        s2 = scores_lang2[matched_items]

        if len(s1) < 10:
            return np.nan, np.nan, np.nan

        correlation = np.corrcoef(s1, s2)[0, 1]
        mean_diff = np.abs(s1 - s2).mean()
        max_diff = np.abs(s1 - s2).max()

        return correlation, mean_diff, max_diff


# ==============================================================================
# Cross-Lingual Validation Strategy
# ==============================================================================


@dataclass
class TranslationValidationPlan:
    """Plan for cross-lingual validation."""

    source_language: Language
    target_languages: list[Language]

    # Translation approach
    method: str  # "back_translation", "parallel_annotation", "expert_review"

    # Sampling
    n_items_per_dimension: int
    stratify_by: list[str]

    # Quality thresholds
    min_agreement: float = 0.80
    min_equivalence: float = 0.90


class CrossLingualValidator:
    """
    Validates measurement equivalence across languages.

    Key methods:
    1. Back-translation validation
    2. Parallel annotation comparison
    3. DIF (Differential Item Functioning) analysis
    """

    def __init__(
        self,
        source_lang: Language,
        target_langs: list[Language],
    ):
        """Initialize cross-lingual validator."""
        self.source_lang = source_lang
        self.target_langs = target_langs

    def create_translation_sample(
        self,
        items: list[AnnotationItem],
        n_per_dimension: int = 30,
    ) -> dict[Language, list[AnnotationItem]]:
        """
        Create stratified sample for translation validation.

        Args:
            items: Source language items
            n_per_dimension: Items per mechanism dimension

        Returns:
            Dict mapping language to items
        """
        # Filter source language items
        source_items = [i for i in items if i.language == self.source_lang]

        # Stratified sample by dimension
        selected = []
        for dim in MechanismDimension:
            dim_items = [i for i in source_items if i.dimension == dim]
            n_select = min(n_per_dimension, len(dim_items))
            selected.extend(np.random.choice(dim_items, n_select, replace=False))

        # Create translation target sets
        result = {self.source_lang: selected}
        for target_lang in self.target_langs:
            result[target_lang] = []  # Placeholder for translations

        return result

    def compute_equivalence_metrics(
        self,
        source_responses: list[AnnotatorResponse],
        target_responses: list[AnnotatorResponse],
    ) -> dict[str, float]:
        """
        Compute cross-lingual equivalence metrics.

        Args:
            source_responses: Responses in source language
            target_responses: Responses in target language (same items)

        Returns:
            Dict of equivalence metrics
        """
        # Match items by ID
        source_by_item = {}
        for r in source_responses:
            if r.item_id not in source_by_item:
                source_by_item[r.item_id] = []
            source_by_item[r.item_id].append(r.label)

        target_by_item = {}
        for r in target_responses:
            if r.item_id not in target_by_item:
                target_by_item[r.item_id] = []
            target_by_item[r.item_id].append(r.label)

        # Find common items
        common_items = set(source_by_item.keys()) & set(target_by_item.keys())

        if len(common_items) < 10:
            logger.warning("Insufficient common items for equivalence analysis")
            return {}

        # Compute agreement on same items
        agreements = []
        for item_id in common_items:
            source_majority = np.mean(source_by_item[item_id]) > 0.5
            target_majority = np.mean(target_by_item[item_id]) > 0.5
            agreements.append(source_majority == target_majority)

        return {
            "cross_lingual_agreement": np.mean(agreements),
            "n_common_items": len(common_items),
            "equivalence_achieved": np.mean(agreements) >= 0.85,
        }


# ==============================================================================
# Main Validation Controller
# ==============================================================================


class MechanismIndexValidator:
    """
    Main controller for Mechanism Index validation.

    Orchestrates:
    - Dataset creation and management
    - Annotation collection
    - Agreement computation
    - Validity assessment
    - Cross-lingual validation
    - Report generation
    """

    def __init__(
        self,
        output_dir: Path,
        config: dict[str, Any] | None = None,
    ):
        """
        Initialize validator.

        Args:
            output_dir: Directory for all validation outputs
            config: Optional configuration dict
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.config = config or {}

        # Sub-components
        self.dataset_manager = ValidationDatasetManager(
            output_dir / "datasets",
            seed=self.config.get("seed", 42),
        )

        logger.info(f"Mechanism Index validator initialized: {output_dir}")

    def run_full_validation(
        self,
        source_data: pd.DataFrame,
        responses_path: Path | None = None,
    ) -> ValidationReport:
        """
        Run complete validation workflow.

        Args:
            source_data: Source procurement data
            responses_path: Path to annotation responses (if available)

        Returns:
            Complete ValidationReport
        """
        logger.info("Starting full validation workflow")

        # 1. Create validation dataset
        items = self.dataset_manager.create_stratified_sample(
            source_data,
            n_per_stratum=self.config.get("n_per_stratum", 50),
        )

        dataset = ValidationDataset(
            name="mechanism_index_validation",
            version="1.0",
            created_date=datetime.now(),
            items=items,
            countries=list(source_data["country"].unique()),
            languages=[Language(l) for l in source_data["language"].unique()],
        )

        # 2. Load or collect responses
        if responses_path and responses_path.exists():
            responses = self.dataset_manager.import_annotations(responses_path)
            dataset.responses = responses
            dataset.n_annotators = len({r.annotator_id for r in responses})
        else:
            logger.warning("No responses provided - exporting annotation template")
            self.dataset_manager.export_annotation_format(
                items,
                self.output_dir / "annotation_template.csv",
            )
            # Return partial report
            return self._create_partial_report(dataset)

        # 3. Compute agreement metrics
        agreement = InterAnnotatorAgreement.compute_all_metrics(
            dataset.responses,
            dataset.items,
        )

        # 4. Compute by dimension
        dimension_results = {}
        for dim in MechanismDimension:
            dim_items = [i for i in dataset.items if i.dimension == dim]
            dim_responses = [r for r in dataset.responses
                           if any(i.item_id == r.item_id for i in dim_items)]
            if dim_responses:
                dimension_results[dim.value] = InterAnnotatorAgreement.compute_all_metrics(
                    dim_responses, dim_items
                )

        # 5. Compute by country
        country_results = {}
        for country in dataset.countries:
            country_items = [i for i in dataset.items if i.country == country]
            country_responses = [r for r in dataset.responses
                                if any(i.item_id == r.item_id for i in country_items)]
            if country_responses:
                country_results[country] = InterAnnotatorAgreement.compute_all_metrics(
                    country_responses, country_items
                )

        # 6. Compute by language
        language_results = {}
        for lang in dataset.languages:
            lang_items = [i for i in dataset.items if i.language == lang]
            lang_responses = [r for r in dataset.responses
                            if any(i.item_id == r.item_id for i in lang_items)]
            if lang_responses:
                language_results[lang.value] = InterAnnotatorAgreement.compute_all_metrics(
                    lang_responses, lang_items
                )

        # 7. Validity assessment (if we have gold standard)
        validity = self._compute_validity_metrics(dataset)

        # 8. Generate report
        report = self._generate_report(
            dataset, agreement, validity,
            dimension_results, country_results, language_results
        )

        # 9. Save report
        self._save_report(report)

        return report

    def _compute_validity_metrics(
        self,
        dataset: ValidationDataset,
    ) -> ValidityMetrics:
        """Compute validity metrics if possible."""
        # This would require additional data (gold standard, external measures)
        # Placeholder implementation
        return ValidityMetrics()

    def _generate_report(
        self,
        dataset: ValidationDataset,
        agreement: AgreementMetrics,
        validity: ValidityMetrics,
        dimension_results: dict,
        country_results: dict,
        language_results: dict,
    ) -> ValidationReport:
        """Generate validation report."""
        concerns = []
        recommendations = []

        # Assess overall quality
        if agreement.fleiss_kappa is not None:
            if agreement.fleiss_kappa >= 0.80:
                quality = "excellent"
            elif agreement.fleiss_kappa >= 0.60:
                quality = "good"
            elif agreement.fleiss_kappa >= 0.40:
                quality = "acceptable"
                concerns.append(f"Moderate agreement (κ={agreement.fleiss_kappa:.2f})")
                recommendations.append("Consider additional annotator training")
            else:
                quality = "poor"
                concerns.append(f"Low agreement (κ={agreement.fleiss_kappa:.2f})")
                recommendations.append("Review annotation guidelines and retrain")
        else:
            quality = "unknown"
            concerns.append("Could not compute agreement metrics")

        # Check dimension consistency
        for dim, metrics in dimension_results.items():
            if metrics.fleiss_kappa and metrics.fleiss_kappa < 0.50:
                concerns.append(f"Low agreement for {dim} dimension")

        # Check cross-country consistency
        kappas = [m.fleiss_kappa for m in country_results.values() if m.fleiss_kappa]
        if kappas and np.std(kappas) > 0.15:
            concerns.append("Substantial variation in agreement across countries")
            recommendations.append("Consider country-specific annotation guidelines")

        return ValidationReport(
            dataset=dataset,
            agreement_metrics=agreement,
            validity_metrics=validity,
            dimension_results=dimension_results,
            country_results=country_results,
            language_results=language_results,
            overall_quality=quality,
            concerns=concerns,
            recommendations=recommendations,
        )

    def _create_partial_report(self, dataset: ValidationDataset) -> ValidationReport:
        """Create partial report when responses are missing."""
        return ValidationReport(
            dataset=dataset,
            agreement_metrics=AgreementMetrics(),
            validity_metrics=ValidityMetrics(),
            dimension_results={},
            country_results={},
            language_results={},
            overall_quality="pending",
            concerns=["No annotations available"],
            recommendations=["Complete annotation task using exported template"],
        )

    def _save_report(self, report: ValidationReport) -> None:
        """Save validation report to files."""
        # Summary JSON
        summary = {
            "dataset": report.dataset.to_dict(),
            "agreement": {
                "cohens_kappa": report.agreement_metrics.cohens_kappa,
                "fleiss_kappa": report.agreement_metrics.fleiss_kappa,
                "krippendorff_alpha": report.agreement_metrics.krippendorff_alpha,
                "percent_agreement": report.agreement_metrics.percent_agreement,
            },
            "quality": report.overall_quality,
            "concerns": report.concerns,
            "recommendations": report.recommendations,
        }

        summary_path = self.output_dir / "validation_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, default=str))

        # Detailed markdown report
        md_lines = [
            "# Mechanism Index Validation Report",
            "",
            f"**Generated:** {datetime.now().isoformat()}",
            f"**Quality Assessment:** {report.overall_quality.upper()}",
            "",
            "## Dataset Summary",
            f"- Items: {len(report.dataset.items)}",
            f"- Annotators: {report.dataset.n_annotators}",
            f"- Countries: {', '.join(report.dataset.countries)}",
            f"- Languages: {', '.join(l.value for l in report.dataset.languages)}",
            "",
            "## Agreement Metrics",
        ]

        if report.agreement_metrics.fleiss_kappa:
            md_lines.append(f"- Fleiss' Kappa: {report.agreement_metrics.fleiss_kappa:.3f}")
        if report.agreement_metrics.krippendorff_alpha:
            md_lines.append(f"- Krippendorff's Alpha: {report.agreement_metrics.krippendorff_alpha:.3f}")
        if report.agreement_metrics.percent_agreement:
            md_lines.append(f"- Percent Agreement: {report.agreement_metrics.percent_agreement:.1%}")

        md_lines.extend([
            "",
            "## Concerns",
        ])
        for concern in report.concerns:
            md_lines.append(f"- ⚠️ {concern}")

        md_lines.extend([
            "",
            "## Recommendations",
        ])
        for rec in report.recommendations:
            md_lines.append(f"- 📋 {rec}")

        md_path = self.output_dir / "validation_report.md"
        md_path.write_text("\n".join(md_lines))

        logger.info(f"Validation report saved to {self.output_dir}")


# ==============================================================================
# Convenience Functions
# ==============================================================================


def create_annotation_template(
    output_dir: Path,
    n_items: int = 500,
) -> Path:
    """
    Create empty annotation template with proper structure.

    Args:
        output_dir: Output directory
        n_items: Number of placeholder items

    Returns:
        Path to template file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create template DataFrame
    template = pd.DataFrame({
        "item_id": [f"item_{i:04d}" for i in range(n_items)],
        "contract_id": [""] * n_items,
        "country": [""] * n_items,
        "language": [""] * n_items,
        "component": [""] * n_items,
        "dimension": [""] * n_items,
        "text_content": ["[Text to annotate]"] * n_items,
        "annotator_id": [""] * n_items,
        "label": [""] * n_items,  # 1 = present, 0 = absent
        "confidence": [""] * n_items,  # 0-1
        "reasoning": [""] * n_items,
        "timestamp": [""] * n_items,
    })

    output_path = output_dir / "annotation_template.csv"
    template.to_csv(output_path, index=False)

    logger.info(f"Annotation template created: {output_path}")
    return output_path


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO)

    # Example usage
    output_dir = Path("outputs/validation")

    # Create validator
    validator = MechanismIndexValidator(output_dir)

    # Create sample data
    np.random.seed(42)
    n_items = 100

    sample_data = pd.DataFrame({
        "contract_id": [f"contract_{i}" for i in range(n_items)],
        "country": np.random.choice(["Ukraine", "Colombia", "UK"], n_items),
        "language": np.random.choice(["en", "uk", "es"], n_items),
        "text_content": [f"Sample procurement text {i}" for i in range(n_items)],
        "component": np.random.choice(
            ["tender_notice_published", "award_notice_published", "electronic_submission"],
            n_items
        ),
        "dimension": np.random.choice(["disclosure", "participation", "oversight"], n_items),
    })

    # Run validation (will export template since no responses)
    report = validator.run_full_validation(sample_data)

    print(f"\nValidation Quality: {report.overall_quality}")
    print(f"Concerns: {report.concerns}")
    print(f"Recommendations: {report.recommendations}")
