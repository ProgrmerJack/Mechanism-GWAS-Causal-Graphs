"""
FLAMES Approximation Implementation
====================================

Approximate FLAMES (Fine-mapped Locus-Associated cis-eQTL Mechanism Scoring)
based on the methodology from Schipper et al. 2025.

FLAMES combines:
1. Fine-mapped credible set annotations (ABC, eQTL, distance, coding)
2. XGBoost classifier trained on curated locus-gene pairs
3. Integration with PoPS convergence scores (0.725 * XGBoost + 0.275 * PoPS)

Reference:
    Schipper et al. (2025). FLAMES: A flexible gene prioritization framework for
    GWAS leveraging functional annotations, fine-mapping, and polygenic methods.
    Nature Genetics (in press, preprint: bioRxiv).

Author: Mechanism-GWAS-Causal-Graphs team
Date: December 2025
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

try:
    import xgboost as xgb
    HAS_XGBOOST = True
except ImportError:
    HAS_XGBOOST = False
    logging.warning("XGBoost not available - FLAMES approximation will be limited")


class FLAMESApproximation:
    """
    Approximate FLAMES scoring using publicly available methodology.
    
    Strategy:
    1. Annotate fine-mapped credible sets with ABC/eQTL/distance scores
    2. Multiply scores by variant PIPs to get locus-level gene scores
    3. Train XGBoost on curated gene-locus pairs
    4. Integrate with PoPS scores using 0.725:0.275 weighting
    
    Attributes:
        xgb_model: XGBoost classifier for locus-gene scoring
        pops_weight: Weight for PoPS integration (default 0.275)
        xgb_weight: Weight for XGBoost integration (default 0.725)
    """
    
    def __init__(
        self,
        pops_weight: float = 0.275,
        xgb_weight: float = 0.725,
        n_estimators: int = 100,
        max_depth: int = 6,
        random_state: int = 42
    ):
        """
        Initialize FLAMES approximation.
        
        Args:
            pops_weight: Weight for PoPS scores in final integration
            xgb_weight: Weight for XGBoost scores in final integration
            n_estimators: Number of XGBoost trees
            max_depth: Maximum tree depth
            random_state: Random seed for reproducibility
        """
        if not HAS_XGBOOST:
            raise ImportError(
                "XGBoost required for FLAMES. Install: pip install xgboost"
            )
        
        self.pops_weight = pops_weight
        self.xgb_weight = xgb_weight
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.random_state = random_state
        
        self.xgb_model = None
        self.feature_names = None
        
        self.logger = logging.getLogger(__name__)
    
    def create_feature_matrix(
        self,
        locus_data: pd.DataFrame,
        gene: str
    ) -> np.ndarray:
        """
        Create feature matrix for a locus-gene pair.
        
        Features (following FLAMES methodology):
        1. Max ABC score (max across variants in credible set)
        2. PIP-weighted ABC score (sum of variant_pip * abc_score)
        3. Max eQTL colocalization PP.H4
        4. PIP-weighted eQTL PP.H4
        5. Distance to TSS (closest variant in credible set)
        6. Coding variant indicator (any coding variant in credible set)
        7. Max variant PIP
        8. Sum of PIPs (total posterior probability in credible set)
        
        Args:
            locus_data: DataFrame with credible set variants
                Required columns: variant_id, pip, abc_score, eqtl_pp4,
                                 distance_to_tss, is_coding
            gene: Gene symbol for this locus-gene pair
        
        Returns:
            Feature vector (8 features)
        """
        # Filter to variants near this gene
        gene_data = locus_data[locus_data['target_gene'] == gene].copy()
        
        if len(gene_data) == 0:
            # No variants linked to this gene - return zero features
            return np.zeros(8)
        
        features = np.array([
            # ABC features
            gene_data['abc_score'].max() if 'abc_score' in gene_data else 0,
            (gene_data['pip'] * gene_data.get('abc_score', 0)).sum(),
            
            # eQTL features
            gene_data.get('eqtl_pp4', 0).max(),
            (gene_data['pip'] * gene_data.get('eqtl_pp4', 0)).sum(),
            
            # Distance feature (log-transform)
            np.log10(gene_data['distance_to_tss'].min() + 1) if 'distance_to_tss' in gene_data else 0,
            
            # Coding variant indicator
            float(gene_data.get('is_coding', False).any()),
            
            # Fine-mapping features
            gene_data['pip'].max(),
            gene_data['pip'].sum()
        ])
        
        return features
    
    def prepare_training_data(
        self,
        credible_sets: Dict[str, pd.DataFrame],
        gold_standard_genes: Dict[str, str]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepare training data from credible sets and curated genes.
        
        Args:
            credible_sets: Dict mapping locus_id → credible set DataFrame
            gold_standard_genes: Dict mapping locus_id → true gene symbol
        
        Returns:
            X: Feature matrix (n_samples, 8)
            y: Binary labels (1 = true gene, 0 = negative example)
        """
        X_list = []
        y_list = []
        
        for locus_id, true_gene in gold_standard_genes.items():
            if locus_id not in credible_sets:
                self.logger.warning(f"Locus {locus_id} missing from credible sets")
                continue
            
            locus_data = credible_sets[locus_id]
            
            # Get all genes at this locus
            candidate_genes = locus_data['target_gene'].unique()
            
            for gene in candidate_genes:
                # Create features
                features = self.create_feature_matrix(locus_data, gene)
                X_list.append(features)
                
                # Label: 1 if true gene, 0 otherwise
                y_list.append(1 if gene == true_gene else 0)
        
        X = np.vstack(X_list)
        y = np.array(y_list)
        
        self.logger.info(
            f"Prepared training data: {len(y)} examples, "
            f"{y.sum()} positives, {(~y.astype(bool)).sum()} negatives"
        )
        
        return X, y
    
    def train(
        self,
        X: np.ndarray,
        y: np.ndarray,
        sample_weight: Optional[np.ndarray] = None
    ) -> Dict[str, float]:
        """
        Train XGBoost classifier on locus-gene features.
        
        Args:
            X: Feature matrix (n_samples, 8)
            y: Binary labels
            sample_weight: Optional sample weights (e.g., for class imbalance)
        
        Returns:
            Training metrics (accuracy, AUC, etc.)
        """
        # Handle class imbalance with scale_pos_weight
        scale_pos_weight = (y == 0).sum() / (y == 1).sum()
        
        self.xgb_model = xgb.XGBClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=0.1,
            scale_pos_weight=scale_pos_weight,
            random_state=self.random_state,
            objective='binary:logistic',
            eval_metric='auc'
        )
        
        # Train model
        self.xgb_model.fit(
            X, y,
            sample_weight=sample_weight,
            verbose=False
        )
        
        # Store feature names
        self.feature_names = [
            'max_abc', 'pip_weighted_abc',
            'max_eqtl_pp4', 'pip_weighted_eqtl',
            'log_min_distance', 'has_coding_variant',
            'max_pip', 'sum_pip'
        ]
        
        # Compute training metrics
        y_pred_proba = self.xgb_model.predict_proba(X)[:, 1]
        
        from sklearn.metrics import roc_auc_score, average_precision_score
        
        metrics = {
            'train_auc': roc_auc_score(y, y_pred_proba),
            'train_auprc': average_precision_score(y, y_pred_proba),
        }
        
        self.logger.info(f"Training complete: AUC={metrics['train_auc']:.3f}")
        
        return metrics
    
    def predict_xgboost_score(
        self,
        locus_data: pd.DataFrame,
        gene: str
    ) -> float:
        """
        Predict XGBoost score for a locus-gene pair.
        
        Args:
            locus_data: Credible set data for locus
            gene: Gene symbol
        
        Returns:
            XGBoost probability score [0, 1]
        """
        if self.xgb_model is None:
            raise ValueError("Model not trained - call train() first")
        
        # Create features
        features = self.create_feature_matrix(locus_data, gene)
        
        # Predict probability
        score = self.xgb_model.predict_proba(features.reshape(1, -1))[0, 1]
        
        return float(score)
    
    def integrate_with_pops(
        self,
        xgb_score: float,
        pops_score: float
    ) -> float:
        """
        Integrate XGBoost and PoPS scores using FLAMES weighting.
        
        FLAMES final score = 0.725 * XGBoost + 0.275 * PoPS
        
        Args:
            xgb_score: XGBoost prediction [0, 1]
            pops_score: PoPS convergence score [0, 1]
        
        Returns:
            Integrated FLAMES score [0, 1]
        """
        flames_score = self.xgb_weight * xgb_score + self.pops_weight * pops_score
        return flames_score
    
    def score_locus(
        self,
        locus_data: pd.DataFrame,
        pops_scores: Optional[Dict[str, float]] = None
    ) -> pd.DataFrame:
        """
        Score all genes at a locus.
        
        Args:
            locus_data: Credible set data with variants and gene links
            pops_scores: Optional dict mapping gene → PoPS score
        
        Returns:
            DataFrame with columns: gene, xgb_score, pops_score, flames_score
        """
        genes = locus_data['target_gene'].unique()
        
        results = []
        for gene in genes:
            xgb_score = self.predict_xgboost_score(locus_data, gene)
            
            # Get PoPS score if available
            pops_score = pops_scores.get(gene, 0.0) if pops_scores else 0.0
            
            # Integrate
            flames_score = self.integrate_with_pops(xgb_score, pops_score)
            
            results.append({
                'gene': gene,
                'xgb_score': xgb_score,
                'pops_score': pops_score,
                'flames_score': flames_score
            })
        
        df = pd.DataFrame(results)
        df = df.sort_values('flames_score', ascending=False)
        
        return df
    
    def get_feature_importance(self) -> pd.DataFrame:
        """
        Get feature importance from trained XGBoost model.
        
        Returns:
            DataFrame with feature names and importance scores
        """
        if self.xgb_model is None:
            raise ValueError("Model not trained")
        
        importance = self.xgb_model.feature_importances_
        
        df = pd.DataFrame({
            'feature': self.feature_names,
            'importance': importance
        })
        df = df.sort_values('importance', ascending=False)
        
        return df
    
    def save_model(self, filepath: Path):
        """Save trained XGBoost model."""
        if self.xgb_model is None:
            raise ValueError("No model to save")
        
        self.xgb_model.save_model(str(filepath))
        self.logger.info(f"Model saved to {filepath}")
    
    def load_model(self, filepath: Path):
        """Load trained XGBoost model."""
        self.xgb_model = xgb.XGBClassifier()
        self.xgb_model.load_model(str(filepath))
        self.logger.info(f"Model loaded from {filepath}")


def main():
    """Example usage."""
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Example: Train FLAMES approximation
    flames = FLAMESApproximation()
    
    # Prepare dummy data (replace with real credible sets)
    credible_sets = {
        'locus_1': pd.DataFrame({
            'variant_id': ['rs1', 'rs2'],
            'pip': [0.8, 0.2],
            'abc_score': [0.9, 0.1],
            'eqtl_pp4': [0.95, 0.05],
            'distance_to_tss': [1000, 5000],
            'is_coding': [True, False],
            'target_gene': ['SORT1', 'PSRC1']
        })
    }
    
    gold_standard = {'locus_1': 'SORT1'}
    
    # Prepare and train
    X, y = flames.prepare_training_data(credible_sets, gold_standard)
    metrics = flames.train(X, y)
    
    print(f"Training metrics: {metrics}")
    print("\nFeature importance:")
    print(flames.get_feature_importance())


if __name__ == "__main__":
    main()
