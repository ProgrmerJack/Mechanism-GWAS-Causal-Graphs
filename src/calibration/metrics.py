"""
Calibration Metrics

Comprehensive metrics for evaluating mechanism graph predictions.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..utils.logging import get_logger


logger = get_logger("metrics")


@dataclass
class CalibrationMetrics:
    """Container for calibration metrics."""
    
    ece: float  # Expected Calibration Error
    mce: float  # Maximum Calibration Error
    brier: float  # Brier Score
    auroc: float  # Area Under ROC
    auprc: float  # Area Under Precision-Recall
    
    def to_dict(self) -> Dict[str, float]:
        return {
            "ece": self.ece,
            "mce": self.mce,
            "brier": self.brier,
            "auroc": self.auroc,
            "auprc": self.auprc,
        }


def compute_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """
    Compute comprehensive evaluation metrics.
    
    Parameters
    ----------
    y_true : np.ndarray
        True binary labels.
    y_prob : np.ndarray
        Predicted probabilities.
    threshold : float
        Classification threshold.
        
    Returns
    -------
    dict
        Dictionary of metrics.
    """
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    y_pred = (y_prob >= threshold).astype(int)
    
    n = len(y_true)
    n_pos = y_true.sum()
    n_neg = n - n_pos
    
    # Classification metrics
    tp = ((y_pred == 1) & (y_true == 1)).sum()
    fp = ((y_pred == 1) & (y_true == 0)).sum()
    tn = ((y_pred == 0) & (y_true == 0)).sum()
    fn = ((y_pred == 0) & (y_true == 1)).sum()
    
    accuracy = (tp + tn) / n if n > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    # Probabilistic metrics
    brier = np.mean((y_prob - y_true) ** 2)
    log_loss = -np.mean(
        y_true * np.log(np.clip(y_prob, 1e-10, 1)) +
        (1 - y_true) * np.log(np.clip(1 - y_prob, 1e-10, 1))
    )
    
    # Discrimination metrics
    auroc = compute_auroc(y_true, y_prob)
    auprc = compute_auprc(y_true, y_prob)
    
    # Calibration metrics
    ece, mce = compute_calibration_errors(y_true, y_prob)
    
    return {
        # Classification
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        
        # Probabilistic
        "brier_score": brier,
        "log_loss": log_loss,
        
        # Discrimination
        "auroc": auroc,
        "auprc": auprc,
        
        # Calibration
        "ece": ece,
        "mce": mce,
        
        # Metadata
        "n_samples": n,
        "n_positives": int(n_pos),
        "base_rate": n_pos / n if n > 0 else 0,
        "threshold": threshold,
    }


def compute_auroc(
    y_true: np.ndarray,
    y_prob: np.ndarray,
) -> float:
    """
    Compute Area Under ROC Curve.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels.
    y_prob : np.ndarray
        Predicted probabilities.
        
    Returns
    -------
    float
        AUROC.
    """
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    
    n_pos = y_true.sum()
    n_neg = len(y_true) - n_pos
    
    if n_pos == 0 or n_neg == 0:
        return 0.5
    
    # Mann-Whitney U statistic
    pairs = list(zip(y_prob, y_true))
    pairs.sort(key=lambda x: -x[0])
    
    auc = 0.0
    pos_above = 0
    
    for prob, label in pairs:
        if label == 0:
            auc += pos_above
        else:
            pos_above += 1
    
    return auc / (n_pos * n_neg)


def compute_auprc(
    y_true: np.ndarray,
    y_prob: np.ndarray,
) -> float:
    """
    Compute Area Under Precision-Recall Curve.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels.
    y_prob : np.ndarray
        Predicted probabilities.
        
    Returns
    -------
    float
        AUPRC.
    """
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    
    n_pos = y_true.sum()
    
    if n_pos == 0:
        return 0.0
    
    # Sort by probability
    order = np.argsort(-y_prob)
    y_true_sorted = y_true[order]
    
    # Compute precision at each recall level
    precisions = []
    recalls = []
    
    tp = 0
    for i in range(len(y_true_sorted)):
        if y_true_sorted[i] == 1:
            tp += 1
            precision = tp / (i + 1)
            recall = tp / n_pos
            precisions.append(precision)
            recalls.append(recall)
    
    if not precisions:
        return 0.0
    
    # Average precision
    return np.mean(precisions)


def compute_calibration_errors(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_bins: int = 10,
) -> Tuple[float, float]:
    """
    Compute ECE and MCE.
    
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
        (ECE, MCE)
    """
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    
    bins = np.linspace(0, 1, n_bins + 1)
    bin_indices = np.digitize(y_prob, bins[:-1]) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    ece = 0.0
    mce = 0.0
    n = len(y_true)
    
    for i in range(n_bins):
        mask = bin_indices == i
        if mask.sum() > 0:
            bin_prob = y_prob[mask].mean()
            bin_true = y_true[mask].mean()
            gap = abs(bin_prob - bin_true)
            
            weight = mask.sum() / n
            ece += weight * gap
            mce = max(mce, gap)
    
    return ece, mce


def compute_ranking_metrics(
    rankings: List[Tuple[str, float]],
    gold_set: set,
) -> Dict[str, float]:
    """
    Compute ranking-based metrics.
    
    Parameters
    ----------
    rankings : list
        List of (item, score) tuples, sorted by score.
    gold_set : set
        Set of gold standard items.
        
    Returns
    -------
    dict
        Ranking metrics.
    """
    n = len(rankings)
    n_gold = len(gold_set)
    
    if n_gold == 0:
        return {"error": "empty_gold_set"}
    
    metrics = {}
    
    # Recall@k
    for k in [10, 20, 50, 100, 200]:
        if k > n:
            continue
        
        top_k = {item for item, _ in rankings[:k]}
        recall = len(top_k & gold_set) / n_gold
        metrics[f"recall@{k}"] = recall
    
    # Precision@k
    for k in [10, 20, 50]:
        if k > n:
            continue
        
        top_k = {item for item, _ in rankings[:k]}
        precision = len(top_k & gold_set) / k
        metrics[f"precision@{k}"] = precision
    
    # Mean Reciprocal Rank
    mrr = 0.0
    for i, (item, _) in enumerate(rankings):
        if item in gold_set:
            mrr = 1 / (i + 1)
            break
    metrics["mrr"] = mrr
    
    # Mean rank of gold standard items
    gold_ranks = []
    for i, (item, _) in enumerate(rankings):
        if item in gold_set:
            gold_ranks.append(i + 1)
    
    if gold_ranks:
        metrics["mean_gold_rank"] = np.mean(gold_ranks)
        metrics["median_gold_rank"] = np.median(gold_ranks)
        metrics["min_gold_rank"] = min(gold_ranks)
        metrics["max_gold_rank"] = max(gold_ranks)
    
    # Normalized Discounted Cumulative Gain
    dcg = 0.0
    for i, (item, _) in enumerate(rankings):
        if item in gold_set:
            dcg += 1 / np.log2(i + 2)
    
    # Ideal DCG
    idcg = sum(1 / np.log2(i + 2) for i in range(min(n_gold, n)))
    
    metrics["ndcg"] = dcg / idcg if idcg > 0 else 0
    
    return metrics


def enrichment_at_threshold(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    threshold: float = 0.5,
) -> float:
    """
    Compute fold enrichment at probability threshold.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels.
    y_prob : np.ndarray
        Predicted probabilities.
    threshold : float
        Probability threshold.
        
    Returns
    -------
    float
        Fold enrichment.
    """
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    
    base_rate = y_true.mean()
    
    high_prob_mask = y_prob >= threshold
    
    if high_prob_mask.sum() == 0:
        return 0.0
    
    observed_rate = y_true[high_prob_mask].mean()
    
    return observed_rate / base_rate if base_rate > 0 else 0.0
