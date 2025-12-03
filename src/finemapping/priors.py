"""
Prior Construction for Fine-Mapping

Implements regulatory and tissue-specific priors for fine-mapping.
This is a key "Nature lever" - using biological annotation to improve resolution.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from ..utils.config import get_config
from ..utils.logging import get_logger


logger = get_logger("priors")


class RegulatoryPriors:
    """
    Constructs prior probabilities based on regulatory annotations.
    
    Variants in regulatory elements (promoters, enhancers) get higher
    prior probabilities than those in non-functional regions.
    """
    
    def __init__(
        self,
        ccre_path: Optional[str] = None,
        weights: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize regulatory priors.
        
        Parameters
        ----------
        ccre_path : str, optional
            Path to ENCODE cCRE BED file.
        weights : dict, optional
            Prior weights by element type.
        """
        self.ccre_path = ccre_path
        
        # Default weights from configuration
        config = get_config("config")
        default_weights = config.get("finemapping", {}).get("priors", {})
        
        self.weights = weights or {
            "promoter": default_weights.get("prior_weight_promoter", 10.0),
            "enhancer": default_weights.get("prior_weight_enhancer", 5.0),
            "ctcf": default_weights.get("prior_weight_ctcf", 2.0),
            "other": default_weights.get("prior_weight_other", 1.0),
        }
        
        self.ccre_data = None
        
        if ccre_path:
            self._load_ccre(ccre_path)
    
    def _load_ccre(self, path: str) -> None:
        """Load ENCODE cCRE annotations."""
        logger.info(f"Loading cCRE annotations from {path}")
        
        self.ccre_data = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "id", "score", "strand",
                   "thickStart", "thickEnd", "rgb", "accession", "type"],
            usecols=["chr", "start", "end", "id", "type"],
        )
        
        # Standardize chromosome names
        self.ccre_data["chr"] = self.ccre_data["chr"].str.replace("chr", "")
        
        logger.info(f"Loaded {len(self.ccre_data):,} cCRE elements")
    
    def annotate_variant(
        self,
        chrom: str,
        pos: int,
    ) -> Dict[str, Any]:
        """
        Annotate a variant with regulatory information.
        
        Parameters
        ----------
        chrom : str
            Chromosome.
        pos : int
            Position.
            
        Returns
        -------
        dict
            Annotation including element type and prior weight.
        """
        if self.ccre_data is None:
            return {"element_type": "none", "prior_weight": self.weights["other"]}
        
        # Find overlapping cCREs
        chrom = str(chrom).replace("chr", "")
        
        overlaps = self.ccre_data[
            (self.ccre_data["chr"] == chrom) &
            (self.ccre_data["start"] <= pos) &
            (self.ccre_data["end"] >= pos)
        ]
        
        if len(overlaps) == 0:
            return {"element_type": "none", "prior_weight": self.weights["other"]}
        
        # Get highest-weight annotation
        best_weight = self.weights["other"]
        best_type = "none"
        ccre_ids = []
        
        for _, ccre in overlaps.iterrows():
            ccre_type = ccre["type"]
            ccre_ids.append(ccre["id"])
            
            # Map cCRE types to weights
            if "PLS" in ccre_type or "promoter" in ccre_type.lower():
                weight = self.weights["promoter"]
                elem_type = "promoter"
            elif "ELS" in ccre_type or "enhancer" in ccre_type.lower():
                weight = self.weights["enhancer"]
                elem_type = "enhancer"
            elif "CTCF" in ccre_type:
                weight = self.weights["ctcf"]
                elem_type = "ctcf"
            else:
                weight = self.weights["other"]
                elem_type = "other_regulatory"
            
            if weight > best_weight:
                best_weight = weight
                best_type = elem_type
        
        return {
            "element_type": best_type,
            "prior_weight": best_weight,
            "ccre_ids": ccre_ids,
        }
    
    def compute_priors(
        self,
        variants: List[Dict],
    ) -> Dict[str, float]:
        """
        Compute prior probabilities for a set of variants.
        
        Parameters
        ----------
        variants : list
            List of variant dictionaries with chr, pos, rsid.
            
        Returns
        -------
        dict
            Prior probabilities by rsid (normalized).
        """
        weights = {}
        
        for variant in variants:
            annotation = self.annotate_variant(
                variant["chr"],
                variant["pos"],
            )
            weights[variant["rsid"]] = annotation["prior_weight"]
        
        # Normalize to probabilities
        total = sum(weights.values())
        priors = {k: v / total for k, v in weights.items()}
        
        return priors


class TissuePriors:
    """
    Constructs tissue-specific priors based on trait-tissue mappings.
    
    Variants active in relevant tissues get higher priors.
    """
    
    def __init__(
        self,
        trait_name: str,
        ccre_activity_path: Optional[str] = None,
    ):
        """
        Initialize tissue priors.
        
        Parameters
        ----------
        trait_name : str
            Name of the trait.
        ccre_activity_path : str, optional
            Path to tissue-specific cCRE activity data.
        """
        self.trait_name = trait_name
        
        # Load trait-tissue configuration
        traits_config = get_config("traits")
        self.trait_config = traits_config.get("trait_tissue_priors", {}).get(
            trait_name, {}
        )
        self.tissue_ontology = traits_config.get("tissue_ontology", {})
        
        # Get relevant tissues and their priors
        self.tissue_priors = self._build_tissue_priors()
        
        # Load activity data if provided
        self.activity_data = None
        if ccre_activity_path:
            self._load_activity_data(ccre_activity_path)
    
    def _build_tissue_priors(self) -> Dict[str, float]:
        """Build tissue prior weights from configuration."""
        priors = {}
        
        # Primary tissues
        for tissue_info in self.trait_config.get("primary_tissues", []):
            tissue = tissue_info["tissue"]
            prior = tissue_info["prior"]
            priors[tissue] = prior
        
        # Secondary tissues
        for tissue_info in self.trait_config.get("secondary_tissues", []):
            tissue = tissue_info["tissue"]
            prior = tissue_info.get("prior", 0.3)
            priors[tissue] = prior
        
        return priors
    
    def _load_activity_data(self, path: str) -> None:
        """Load tissue-specific cCRE activity matrix."""
        logger.info(f"Loading tissue activity data from {path}")
        
        # Format: cCRE_id, tissue1_activity, tissue2_activity, ...
        self.activity_data = pd.read_csv(path, sep="\t", index_col=0)
        
        logger.info(f"Loaded activity for {len(self.activity_data):,} elements")
    
    def compute_tissue_weight(
        self,
        ccre_ids: List[str],
    ) -> float:
        """
        Compute tissue-weighted score for cCRE(s).
        
        Parameters
        ----------
        ccre_ids : list
            List of cCRE identifiers.
            
        Returns
        -------
        float
            Tissue-weighted score.
        """
        if self.activity_data is None or not ccre_ids:
            return 1.0
        
        weights = []
        
        for ccre_id in ccre_ids:
            if ccre_id not in self.activity_data.index:
                continue
            
            activity = self.activity_data.loc[ccre_id]
            
            # Compute weighted score across relevant tissues
            tissue_score = 0.0
            for tissue, prior in self.tissue_priors.items():
                tissue_col = self.tissue_ontology.get(tissue, {}).get("encode_biosample")
                if tissue_col and tissue_col in activity.index:
                    tissue_score += prior * activity[tissue_col]
            
            weights.append(tissue_score)
        
        return max(weights) if weights else 1.0
    
    def compute_combined_priors(
        self,
        variants: List[Dict],
        regulatory_priors: Dict[str, float],
        ccre_annotations: Dict[str, List[str]],
    ) -> Dict[str, float]:
        """
        Compute combined regulatory + tissue priors.
        
        Parameters
        ----------
        variants : list
            List of variant dictionaries.
        regulatory_priors : dict
            Regulatory priors by rsid.
        ccre_annotations : dict
            cCRE IDs by rsid.
            
        Returns
        -------
        dict
            Combined prior probabilities.
        """
        combined = {}
        
        for variant in variants:
            rsid = variant["rsid"]
            
            # Get regulatory prior
            reg_prior = regulatory_priors.get(rsid, 1.0)
            
            # Get tissue weight
            ccre_ids = ccre_annotations.get(rsid, [])
            tissue_weight = self.compute_tissue_weight(ccre_ids)
            
            # Combine multiplicatively
            combined[rsid] = reg_prior * tissue_weight
        
        # Normalize
        total = sum(combined.values())
        if total > 0:
            combined = {k: v / total for k, v in combined.items()}
        
        return combined


def build_priors(
    variants: List[Dict],
    trait_name: str,
    ccre_path: Optional[str] = None,
    ccre_activity_path: Optional[str] = None,
) -> Dict[str, float]:
    """
    Build combined priors for fine-mapping.
    
    Parameters
    ----------
    variants : list
        List of variant dictionaries.
    trait_name : str
        Name of the trait.
    ccre_path : str, optional
        Path to cCRE BED file.
    ccre_activity_path : str, optional
        Path to tissue activity data.
        
    Returns
    -------
    dict
        Prior probabilities by rsid.
    """
    # Build regulatory priors
    reg_priors = RegulatoryPriors(ccre_path)
    regulatory = reg_priors.compute_priors(variants)
    
    # Build tissue priors
    tissue_priors = TissuePriors(trait_name, ccre_activity_path)
    
    # Get cCRE annotations
    ccre_annotations = {}
    for variant in variants:
        annotation = reg_priors.annotate_variant(variant["chr"], variant["pos"])
        ccre_annotations[variant["rsid"]] = annotation.get("ccre_ids", [])
    
    # Combine
    combined = tissue_priors.compute_combined_priors(
        variants, regulatory, ccre_annotations
    )
    
    return combined
