"""
CRISPR Perturbation Validation Framework
=========================================

Validates baseline methods against independent CRISPR perturbation data.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score
import logging


class CRISPRValidator:
    """Validates predictions against CRISPR perturbation screens."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def load_fulco_data(self, filepath: str) -> pd.DataFrame:
        """Load Fulco 2019 CRISPRi-FlowFISH data."""
        try:
            df = pd.read_excel(filepath, sheet_name='Table S6a')
            self.logger.info(f"Loaded Fulco data: {len(df)} entries")
            return df
        except Exception as e:
            self.logger.error(f"Failed to load Fulco data: {e}")
            return pd.DataFrame()
    
    def load_gasperini_data(self, filepath: str) -> pd.DataFrame:
        """Load Gasperini 2019 screen results."""
        try:
            import gzip
            df = pd.read_csv(filepath, sep='\t', compression='gzip')
            self.logger.info(f"Loaded Gasperini data: {len(df)} entries")
            return df
        except Exception as e:
            self.logger.error(f"Failed to load Gasperini data: {e}")
            return pd.DataFrame()
    
    def validate_predictions(
        self,
        predictions: pd.DataFrame,
        validation_data: pd.DataFrame,
        score_column: str = 'score'
    ) -> Dict[str, float]:
        """
        Validate predictions against validation data.
        
        Returns dictionary with metrics.
        """
        # Create labels
        validated_genes = set(validation_data['target_gene'].unique())
        predictions_copy = predictions.copy()
        predictions_copy['validated'] = predictions_copy['gene'].isin(validated_genes).astype(int)
        
        # Sort by score
        predictions_copy = predictions_copy.sort_values(score_column, ascending=False)
        
        y_true = predictions_copy['validated'].values
        y_score = predictions_copy[score_column].values
        
        # Compute metrics
        metrics = {}
        
        for k in [10, 20, 50, 100]:
            if k <= len(y_true):
                metrics[f'precision@{k}'] = y_true[:k].sum() / k
                if y_true.sum() > 0:
                    metrics[f'recall@{k}'] = y_true[:k].sum() / y_true.sum()
        
        try:
            precision_vals, recall_vals, _ = precision_recall_curve(y_true, y_score)
            metrics['auprc'] = auc(recall_vals, precision_vals)
        except:
            metrics['auprc'] = 0.0
        
        try:
            if len(set(y_true)) > 1:
                metrics['auroc'] = roc_auc_score(y_true, y_score)
            else:
                metrics['auroc'] = 0.0
        except:
            metrics['auroc'] = 0.0
        
        return metrics
