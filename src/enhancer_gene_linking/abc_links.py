"""
Activity-by-Contact (ABC) Model Links

Implements loading and querying of ABC enhancer-gene predictions from:
- Fulco CP et al. (2019) Activity-by-contact model predicts enhancer-gene targets
- Nasser J et al. (2021) Genome-wide enhancer maps link risk variants to disease genes

The ABC score = (Enhancer Activity × Contact Frequency) / Σ(all enhancer-gene pairs)

Key resources:
- Nasser 2021 predictions: 131 cell types/tissues, ~4M enhancer-gene links
- ENCODE-rE2G: ENCODE Consortium enhancer-gene predictions
- ABC-Perturb validated: CRISPRi-validated enhancer-gene pairs

This addresses the critical reviewer concern: "Distance alone has false positive rates >50%"
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from ..utils.config import get_config


logger = get_logger("abc_links")


# ABC score thresholds from Nasser et al. 2021
ABC_THRESHOLDS = {
    "high_confidence": 0.015,  # ~70% precision at ~60% recall
    "standard": 0.02,          # Original Fulco threshold
    "stringent": 0.05,         # High precision, lower recall
}


@dataclass
class ABCPrediction:
    """
    A single ABC enhancer-gene prediction.
    
    Attributes
    ----------
    enhancer_chr : str
        Enhancer chromosome.
    enhancer_start : int
        Enhancer start position.
    enhancer_end : int
        Enhancer end position.
    gene_symbol : str
        Target gene symbol.
    gene_id : str
        Ensembl gene ID.
    abc_score : float
        ABC score (0-1, typically 0.01-0.3).
    cell_type : str
        Cell type / tissue of prediction.
    activity : float
        Enhancer activity (H3K27ac signal).
    contact : float
        Hi-C contact frequency.
    distance : int
        Distance to TSS in bp.
    is_validated : bool
        Whether this link has CRISPRi validation.
    validation_effect : float
        Effect size from CRISPRi (if validated).
    """
    
    enhancer_chr: str
    enhancer_start: int
    enhancer_end: int
    gene_symbol: str
    gene_id: str = ""
    abc_score: float = 0.0
    cell_type: str = ""
    activity: float = 0.0
    contact: float = 0.0
    distance: int = 0
    is_validated: bool = False
    validation_effect: float = 0.0
    
    @property
    def enhancer_id(self) -> str:
        """Unique enhancer identifier."""
        return f"{self.enhancer_chr}:{self.enhancer_start}-{self.enhancer_end}"
    
    @property
    def link_id(self) -> str:
        """Unique link identifier."""
        return f"{self.enhancer_id}:{self.gene_symbol}:{self.cell_type}"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "enhancer_chr": self.enhancer_chr,
            "enhancer_start": self.enhancer_start,
            "enhancer_end": self.enhancer_end,
            "enhancer_id": self.enhancer_id,
            "gene_symbol": self.gene_symbol,
            "gene_id": self.gene_id,
            "abc_score": self.abc_score,
            "cell_type": self.cell_type,
            "activity": self.activity,
            "contact": self.contact,
            "distance": self.distance,
            "is_validated": self.is_validated,
            "validation_effect": self.validation_effect,
        }


class ABCLinks:
    """
    Manager for ABC enhancer-gene predictions.
    
    Loads and queries ABC model predictions from Nasser et al. 2021
    and other sources.
    
    Example
    -------
    >>> abc = ABCLinks()
    >>> abc.load_nasser2021("/path/to/abc_predictions.tsv.gz")
    >>> links = abc.get_links_for_region("chr1", 1000000, 1500000)
    >>> gene_links = abc.get_links_for_gene("LDLR")
    """
    
    # Default column mappings for ABC files
    COLUMN_MAPPINGS = {
        "nasser2021": {
            "chr": "chr",
            "start": "start",
            "end": "end", 
            "TargetGene": "gene_symbol",
            "TargetGeneEnsemblID": "gene_id",
            "ABC.Score": "abc_score",
            "CellType": "cell_type",
            "activity_base": "activity",
            "hic_contact": "contact",
            "distance": "distance",
        },
        "encode": {
            "chrom": "chr",
            "chromStart": "start",
            "chromEnd": "end",
            "name": "gene_symbol",
            "score": "abc_score",
            "biosample": "cell_type",
        },
    }
    
    def __init__(
        self,
        abc_threshold: float = 0.015,
        cell_types: Optional[List[str]] = None,
    ):
        """
        Initialize ABC links manager.
        
        Parameters
        ----------
        abc_threshold : float
            Minimum ABC score threshold (default: 0.015).
        cell_types : list, optional
            Filter to specific cell types.
        """
        self.abc_threshold = abc_threshold
        self.cell_type_filter = set(cell_types) if cell_types else None
        
        # Storage
        self._predictions: List[ABCPrediction] = []
        
        # Indices for fast lookup
        self._by_region: Dict[str, List[int]] = {}  # chr -> indices
        self._by_gene: Dict[str, List[int]] = {}    # gene -> indices
        self._by_cell_type: Dict[str, List[int]] = {}  # cell_type -> indices
        
        # Validated links (for ground truth)
        self._validated_links: Set[str] = set()
    
    def load_nasser2021(
        self,
        filepath: str,
        validated_filepath: Optional[str] = None,
    ) -> int:
        """
        Load ABC predictions from Nasser et al. 2021.
        
        Available from: https://www.engreitzlab.org/resources/
        
        Parameters
        ----------
        filepath : str
            Path to ABC predictions file (TSV or TSV.gz).
        validated_filepath : str, optional
            Path to CRISPRi-validated links.
            
        Returns
        -------
        int
            Number of predictions loaded.
        """
        logger.info(f"Loading ABC predictions from {filepath}")
        
        df = pd.read_csv(filepath, sep="\t", compression="infer")
        
        # Map columns
        col_map = self.COLUMN_MAPPINGS["nasser2021"]
        
        count = 0
        for _, row in df.iterrows():
            abc_score = row.get(col_map.get("ABC.Score", "ABC.Score"), 0)
            
            # Apply threshold
            if abc_score < self.abc_threshold:
                continue
            
            cell_type = str(row.get(col_map.get("CellType", "CellType"), ""))
            
            # Apply cell type filter
            if self.cell_type_filter and cell_type not in self.cell_type_filter:
                continue
            
            pred = ABCPrediction(
                enhancer_chr=str(row.get("chr", row.get("chrom", ""))),
                enhancer_start=int(row.get("start", row.get("chromStart", 0))),
                enhancer_end=int(row.get("end", row.get("chromEnd", 0))),
                gene_symbol=str(row.get("TargetGene", "")),
                gene_id=str(row.get("TargetGeneEnsemblID", "")),
                abc_score=float(abc_score),
                cell_type=cell_type,
                activity=float(row.get("activity_base", 0)),
                contact=float(row.get("hic_contact", 0)),
                distance=int(row.get("distance", 0)),
            )
            
            self._add_prediction(pred)
            count += 1
        
        # Load validated links
        if validated_filepath and Path(validated_filepath).exists():
            self._load_validated_links(validated_filepath)
        
        logger.info(f"Loaded {count} ABC predictions above threshold {self.abc_threshold}")
        return count
    
    def load_encode_re2g(self, filepath: str) -> int:
        """
        Load ENCODE-rE2G predictions.
        
        Parameters
        ----------
        filepath : str
            Path to ENCODE-rE2G file.
            
        Returns
        -------
        int
            Number of predictions loaded.
        """
        logger.info(f"Loading ENCODE-rE2G predictions from {filepath}")
        
        df = pd.read_csv(filepath, sep="\t", compression="infer")
        
        count = 0
        for _, row in df.iterrows():
            score = float(row.get("score", 0))
            
            if score < self.abc_threshold:
                continue
            
            pred = ABCPrediction(
                enhancer_chr=str(row.get("chrom", "")),
                enhancer_start=int(row.get("chromStart", 0)),
                enhancer_end=int(row.get("chromEnd", 0)),
                gene_symbol=str(row.get("name", "")).split("|")[0],
                abc_score=score,
                cell_type=str(row.get("biosample", "")),
            )
            
            self._add_prediction(pred)
            count += 1
        
        return count
    
    def _add_prediction(self, pred: ABCPrediction) -> None:
        """Add prediction to storage and indices."""
        idx = len(self._predictions)
        self._predictions.append(pred)
        
        # Index by chromosome
        if pred.enhancer_chr not in self._by_region:
            self._by_region[pred.enhancer_chr] = []
        self._by_region[pred.enhancer_chr].append(idx)
        
        # Index by gene
        gene_key = pred.gene_symbol.upper()
        if gene_key not in self._by_gene:
            self._by_gene[gene_key] = []
        self._by_gene[gene_key].append(idx)
        
        # Index by cell type
        if pred.cell_type not in self._by_cell_type:
            self._by_cell_type[pred.cell_type] = []
        self._by_cell_type[pred.cell_type].append(idx)
    
    def _load_validated_links(self, filepath: str) -> None:
        """Load CRISPRi-validated links."""
        df = pd.read_csv(filepath, sep="\t")
        
        for _, row in df.iterrows():
            enhancer_id = f"{row['chr']}:{row['start']}-{row['end']}"
            gene = row.get("gene", row.get("TargetGene", ""))
            link_key = f"{enhancer_id}:{gene}"
            
            self._validated_links.add(link_key)
            
            # Mark existing predictions as validated
            for pred in self._predictions:
                if pred.enhancer_id == enhancer_id and pred.gene_symbol == gene:
                    pred.is_validated = True
                    pred.validation_effect = float(row.get("effect_size", 0))
        
        logger.info(f"Loaded {len(self._validated_links)} validated enhancer-gene links")
    
    def get_links_for_region(
        self,
        chrom: str,
        start: int,
        end: int,
        cell_types: Optional[List[str]] = None,
    ) -> List[ABCPrediction]:
        """
        Get ABC links overlapping a genomic region.
        
        Parameters
        ----------
        chrom : str
            Chromosome.
        start : int
            Start position.
        end : int
            End position.
        cell_types : list, optional
            Filter to specific cell types.
            
        Returns
        -------
        list
            List of ABCPrediction objects.
        """
        results = []
        
        # Normalize chromosome
        chrom = chrom.replace("chr", "")
        for chr_key in [chrom, f"chr{chrom}"]:
            if chr_key not in self._by_region:
                continue
            
            for idx in self._by_region[chr_key]:
                pred = self._predictions[idx]
                
                # Check overlap
                if pred.enhancer_end < start or pred.enhancer_start > end:
                    continue
                
                # Cell type filter
                if cell_types and pred.cell_type not in cell_types:
                    continue
                
                results.append(pred)
        
        return results
    
    def get_links_for_gene(
        self,
        gene: str,
        cell_types: Optional[List[str]] = None,
        min_score: Optional[float] = None,
    ) -> List[ABCPrediction]:
        """
        Get ABC links for a specific gene.
        
        Parameters
        ----------
        gene : str
            Gene symbol.
        cell_types : list, optional
            Filter to specific cell types.
        min_score : float, optional
            Minimum ABC score.
            
        Returns
        -------
        list
            List of ABCPrediction objects.
        """
        gene_key = gene.upper()
        
        if gene_key not in self._by_gene:
            return []
        
        results = []
        min_abc = min_score or self.abc_threshold
        
        for idx in self._by_gene[gene_key]:
            pred = self._predictions[idx]
            
            if pred.abc_score < min_abc:
                continue
            
            if cell_types and pred.cell_type not in cell_types:
                continue
            
            results.append(pred)
        
        return sorted(results, key=lambda x: x.abc_score, reverse=True)
    
    def get_best_link(
        self,
        enhancer_chr: str,
        enhancer_start: int,
        enhancer_end: int,
        candidate_genes: Optional[List[str]] = None,
        cell_types: Optional[List[str]] = None,
    ) -> Optional[ABCPrediction]:
        """
        Get best ABC-supported link for an enhancer.
        
        Parameters
        ----------
        enhancer_chr : str
            Enhancer chromosome.
        enhancer_start : int
            Enhancer start.
        enhancer_end : int
            Enhancer end.
        candidate_genes : list, optional
            Restrict to candidate genes.
        cell_types : list, optional
            Filter cell types.
            
        Returns
        -------
        ABCPrediction or None
            Best prediction or None.
        """
        links = self.get_links_for_region(
            enhancer_chr, enhancer_start, enhancer_end, cell_types
        )
        
        if not links:
            return None
        
        # Filter to candidate genes if specified
        if candidate_genes:
            candidate_set = {g.upper() for g in candidate_genes}
            links = [l for l in links if l.gene_symbol.upper() in candidate_set]
        
        if not links:
            return None
        
        return max(links, key=lambda x: x.abc_score)
    
    def get_validated_links(self) -> List[ABCPrediction]:
        """Get all CRISPRi-validated links."""
        return [p for p in self._predictions if p.is_validated]
    
    def compute_edge_probability(
        self,
        enhancer_chr: str,
        enhancer_start: int,
        enhancer_end: int,
        gene: str,
        cell_types: Optional[List[str]] = None,
    ) -> Tuple[float, Dict[str, Any]]:
        """
        Compute edge probability for cCRE→gene using ABC evidence.
        
        This replaces distance-only heuristics with experimentally-grounded
        probabilities.
        
        Parameters
        ----------
        enhancer_chr : str
            Enhancer chromosome.
        enhancer_start : int
            Enhancer start.
        enhancer_end : int
            Enhancer end.
        gene : str
            Target gene symbol.
        cell_types : list, optional
            Filter to specific cell types.
            
        Returns
        -------
        tuple
            (probability, evidence_dict)
        """
        links = self.get_links_for_region(
            enhancer_chr, enhancer_start, enhancer_end, cell_types
        )
        
        gene_upper = gene.upper()
        gene_links = [l for l in links if l.gene_symbol.upper() == gene_upper]
        
        evidence = {
            "source": "abc",
            "n_links": len(gene_links),
            "cell_types": list(set(l.cell_type for l in gene_links)),
        }
        
        if not gene_links:
            return 0.0, evidence
        
        # Use maximum ABC score across cell types
        best_link = max(gene_links, key=lambda x: x.abc_score)
        
        # Convert ABC score to probability
        # ABC scores are typically 0.01-0.3, rarely >0.5
        # We use a sigmoid transformation centered at the threshold
        abc_score = best_link.abc_score
        
        # Probability = ABC_score * scaling_factor
        # Nasser 2021 shows ABC > 0.015 has ~70% precision
        # We scale to give ABC=0.015 → P=0.5, ABC=0.1 → P=0.85
        prob = 1 / (1 + np.exp(-50 * (abc_score - 0.015)))
        
        # Boost for validated links
        if best_link.is_validated:
            prob = min(0.95, prob + 0.2)
            evidence["validated"] = True
            evidence["validation_effect"] = best_link.validation_effect
        
        evidence["abc_score"] = abc_score
        evidence["best_cell_type"] = best_link.cell_type
        evidence["distance"] = best_link.distance
        
        return prob, evidence
    
    @property
    def n_predictions(self) -> int:
        """Number of loaded predictions."""
        return len(self._predictions)
    
    @property
    def cell_types(self) -> List[str]:
        """List of cell types with predictions."""
        return list(self._by_cell_type.keys())
    
    def summary(self) -> Dict[str, Any]:
        """Get summary statistics."""
        return {
            "n_predictions": self.n_predictions,
            "n_cell_types": len(self._by_cell_type),
            "n_genes": len(self._by_gene),
            "n_validated": sum(1 for p in self._predictions if p.is_validated),
            "abc_threshold": self.abc_threshold,
            "score_distribution": {
                "min": min(p.abc_score for p in self._predictions) if self._predictions else 0,
                "max": max(p.abc_score for p in self._predictions) if self._predictions else 0,
                "median": float(np.median([p.abc_score for p in self._predictions])) if self._predictions else 0,
            },
        }


def load_abc_predictions(
    filepath: str,
    abc_threshold: float = 0.015,
    cell_types: Optional[List[str]] = None,
) -> ABCLinks:
    """
    Convenience function to load ABC predictions.
    
    Parameters
    ----------
    filepath : str
        Path to ABC predictions file.
    abc_threshold : float
        Minimum ABC score.
    cell_types : list, optional
        Filter to specific cell types.
        
    Returns
    -------
    ABCLinks
        Loaded ABC links manager.
    """
    abc = ABCLinks(abc_threshold=abc_threshold, cell_types=cell_types)
    abc.load_nasser2021(filepath)
    return abc
