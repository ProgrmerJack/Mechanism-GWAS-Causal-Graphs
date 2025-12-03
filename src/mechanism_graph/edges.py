"""
Edge Probability Computation

Computes calibrated edge probabilities between nodes in the mechanism graph.
Each edge represents a step in the causal chain with associated uncertainty.

Key improvements over distance-only linking:
- ABC Model scores (Nasser 2021) for activity-contact based linking
- PCHi-C contacts (Jung 2019, Javierre 2016) for physical chromatin links
- Ensemble approach combining multiple evidence types
- Bridge ablation to demonstrate ABC/PCHi-C contribute beyond distance
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from .nodes import BaseNode, VariantNode, CCRENode, GeneNode, TissueNode, TraitNode


@dataclass
class EdgeProbability:
    """
    Represents an edge in the mechanism graph.
    
    Attributes
    ----------
    source : BaseNode
        Source node.
    target : BaseNode
        Target node.
    probability : float
        Edge probability (0-1).
    evidence : dict
        Supporting evidence for this edge.
    edge_type : str
        Type of relationship.
    confidence : float
        Confidence in the probability estimate.
    method : str
        Method used to compute probability.
    """
    
    source: BaseNode
    target: BaseNode
    probability: float
    evidence: Dict[str, Any] = field(default_factory=dict)
    edge_type: str = ""
    confidence: float = 1.0
    method: str = "default"
    
    def __post_init__(self):
        if not self.edge_type:
            self.edge_type = f"{self.source.type}_to_{self.target.type}"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "source": self.source.id,
            "target": self.target.id,
            "source_type": self.source.type,
            "target_type": self.target.type,
            "probability": self.probability,
            "confidence": self.confidence,
            "edge_type": self.edge_type,
            "method": self.method,
            "evidence": self.evidence,
        }
    
    def __hash__(self):
        return hash((self.source.id, self.target.id))


class VariantToCCREEdge:
    """
    Computes probability that a variant affects a cCRE.
    
    Evidence sources:
    - Physical overlap (variant in cCRE)
    - Variant annotation (coding, splice, regulatory)
    - ChromHMM/functional state
    """
    
    @staticmethod
    def compute(
        variant: VariantNode,
        ccre: CCRENode,
        annotations: Optional[Dict] = None,
    ) -> EdgeProbability:
        """
        Compute variant → cCRE edge probability.
        
        Parameters
        ----------
        variant : VariantNode
            Source variant.
        ccre : CCRENode
            Target cCRE.
        annotations : dict, optional
            Additional variant annotations.
            
        Returns
        -------
        EdgeProbability
            Edge with computed probability.
        """
        evidence = {}
        prob = 0.0
        
        # Check physical overlap
        if ccre.contains_position(variant.pos):
            overlap_prob = 0.8  # High probability if overlapping
            evidence["overlap"] = True
            prob = overlap_prob
        else:
            # Distance-based decay
            distance = min(
                abs(variant.pos - ccre.start),
                abs(variant.pos - ccre.end),
            )
            
            # Exponential decay with 10kb half-life
            distance_prob = np.exp(-distance / 10000) * 0.5
            evidence["distance_bp"] = distance
            evidence["overlap"] = False
            prob = distance_prob
        
        # Boost for promoter elements
        if ccre.element_type in ["PLS", "promoter"]:
            prob *= 1.2
            evidence["element_boost"] = "promoter"
        elif ccre.element_type in ["ELS", "enhancer"]:
            prob *= 1.1
            evidence["element_boost"] = "enhancer"
        
        # Cap at 1.0
        prob = min(prob, 1.0)
        
        # Weight by variant PIP
        weighted_prob = prob * variant.pip
        evidence["variant_pip"] = variant.pip
        
        return EdgeProbability(
            source=variant,
            target=ccre,
            probability=weighted_prob,
            evidence=evidence,
            edge_type="variant_affects_ccre",
        )


class CCREToGeneEdge:
    """
    Computes probability that a cCRE regulates a gene.
    
    Evidence sources (hierarchical):
    1. ABC Model score - Activity-by-Contact (Nasser 2021)
    2. PCHi-C contact - Promoter Capture Hi-C (Jung 2019)
    3. Distance-based fallback
    
    This addresses reviewer concern about distance-only linking.
    """
    
    def __init__(
        self,
        abc_linker: Optional[Any] = None,
        pchic_linker: Optional[Any] = None,
        ensemble_linker: Optional[Any] = None,
        method: str = "ensemble",
    ):
        """
        Initialize edge computer.
        
        Parameters
        ----------
        abc_linker : ABCLinks, optional
            ABC model linker.
        pchic_linker : PCHiCLinks, optional
            PCHi-C linker.
        ensemble_linker : EnsembleLinker, optional
            Ensemble linker.
        method : str
            Linking method: "abc", "pchic", "distance", "ensemble".
        """
        self.abc_linker = abc_linker
        self.pchic_linker = pchic_linker
        self.ensemble_linker = ensemble_linker
        self.method = method
    
    def compute(
        self,
        ccre: CCRENode,
        gene: GeneNode,
        cell_type: Optional[str] = None,
        abc_score: Optional[float] = None,
        pchic_score: Optional[float] = None,
        hic_contact: Optional[float] = None,
    ) -> EdgeProbability:
        """
        Compute cCRE → gene edge probability.
        
        Uses hierarchical evidence:
        1. Ensemble (ABC + PCHi-C + distance) if available
        2. ABC score if available
        3. PCHi-C score if available
        4. Distance-based fallback
        
        Parameters
        ----------
        ccre : CCRENode
            Source cCRE.
        gene : GeneNode
            Target gene.
        cell_type : str, optional
            Cell type for ABC lookup.
        abc_score : float, optional
            Pre-computed ABC score.
        pchic_score : float, optional
            Pre-computed PCHi-C score.
        hic_contact : float, optional
            Raw Hi-C contact frequency.
            
        Returns
        -------
        EdgeProbability
            Edge with computed probability.
        """
        evidence = {}
        method_used = "distance"
        
        # Check if same chromosome
        if str(ccre.chr) != str(gene.chr):
            return EdgeProbability(
                source=ccre,
                target=gene,
                probability=0.0,
                evidence={"error": "different_chromosomes"},
                method="none",
            )
        
        # Distance to TSS
        distance = abs((ccre.start + ccre.end) // 2 - gene.tss)
        evidence["distance_to_tss"] = distance
        
        # Try to get ABC score from linker if not provided
        if abc_score is None and self.abc_linker is not None:
            try:
                abc_score = self.abc_linker.get_link_score(
                    ccre.chr, ccre.start, ccre.end, gene.id, cell_type
                )
            except Exception:
                pass
        
        # Try to get PCHi-C score from linker if not provided
        if pchic_score is None and self.pchic_linker is not None:
            try:
                pchic_score = self.pchic_linker.get_link_score(
                    ccre.chr, ccre.start, ccre.end, gene.id
                )
            except Exception:
                pass
        
        # Use ensemble linker if available
        if self.method == "ensemble" and self.ensemble_linker is not None:
            try:
                prob, link_evidence = self.ensemble_linker.get_ensemble_score(
                    ccre.chr, ccre.start, ccre.end, gene.id, cell_type
                )
                evidence.update(link_evidence)
                evidence["abc_score"] = abc_score
                evidence["pchic_score"] = pchic_score
                method_used = "ensemble"
                
                return EdgeProbability(
                    source=ccre,
                    target=gene,
                    probability=min(prob, 1.0),
                    evidence=evidence,
                    edge_type="ccre_regulates_gene",
                    method=method_used,
                )
            except Exception:
                pass  # Fall through to individual methods
        
        # Individual method handling
        prob = 0.0
        
        # ABC score (strongest evidence)
        if abc_score is not None and abc_score > 0:
            # ABC score is already calibrated to represent regulatory probability
            prob = abc_score
            evidence["abc_score"] = abc_score
            method_used = "abc"
        # PCHi-C contact
        elif pchic_score is not None and pchic_score > 0:
            # PCHi-C CHiCAGO score, convert to probability
            # CHiCAGO >= 5 is significant contact
            prob = 1 - np.exp(-pchic_score / 5)
            evidence["pchic_score"] = pchic_score
            method_used = "pchic"
        # Distance-based fallback
        else:
            # Promoter cCRE very close to TSS
            if ccre.element_type in ["PLS", "promoter"] and distance < 2000:
                prob = 0.9
                evidence["promoter_proximal"] = True
            else:
                # Distance-based prior
                # ~100kb effective range for enhancers
                distance_prob = np.exp(-distance / 100000) * 0.5
                prob = distance_prob
                evidence["promoter_proximal"] = False
            
            method_used = "distance"
        
        # Use Hi-C contact to boost if available
        if hic_contact is not None and method_used == "distance":
            prob = prob * (0.5 + 0.5 * min(hic_contact / 10, 1.0))
            evidence["hic_contact"] = hic_contact
        
        return EdgeProbability(
            source=ccre,
            target=gene,
            probability=min(prob, 1.0),
            evidence=evidence,
            edge_type="ccre_regulates_gene",
            method=method_used,
        )
    
    @staticmethod
    def compute_distance_only(
        ccre: CCRENode,
        gene: GeneNode,
    ) -> EdgeProbability:
        """
        Compute cCRE → gene edge using distance only.
        
        Used for baseline/ablation comparison.
        
        Parameters
        ----------
        ccre : CCRENode
            Source cCRE.
        gene : GeneNode
            Target gene.
            
        Returns
        -------
        EdgeProbability
            Edge with distance-only probability.
        """
        evidence = {}
        
        # Check if same chromosome
        if str(ccre.chr) != str(gene.chr):
            return EdgeProbability(
                source=ccre,
                target=gene,
                probability=0.0,
                evidence={"error": "different_chromosomes"},
                method="distance_only",
            )
        
        # Distance to TSS
        distance = abs((ccre.start + ccre.end) // 2 - gene.tss)
        evidence["distance_to_tss"] = distance
        
        # Promoter cCRE very close to TSS
        if ccre.element_type in ["PLS", "promoter"] and distance < 2000:
            prob = 0.9
            evidence["promoter_proximal"] = True
        else:
            # Distance-based prior (~100kb effective range)
            prob = np.exp(-distance / 100000) * 0.5
            evidence["promoter_proximal"] = False
        
        return EdgeProbability(
            source=ccre,
            target=gene,
            probability=min(prob, 1.0),
            evidence=evidence,
            edge_type="ccre_regulates_gene",
            method="distance_only",
        )


class GeneToTissueEdge:
    """
    Computes probability that gene-trait association is mediated by a tissue.
    
    Evidence sources:
    - Gene expression in tissue
    - eQTL colocalization (PP.H4)
    - sQTL colocalization
    - Tissue relevance prior
    - Cross-study replication via eQTL Catalogue
    """
    
    def __init__(
        self,
        replication_tester: Optional[Any] = None,
    ):
        """
        Initialize edge computer.
        
        Parameters
        ----------
        replication_tester : CrossStudyReplication, optional
            For testing GTEx eQTL replication in eQTL Catalogue.
        """
        self.replication_tester = replication_tester
    
    def compute(
        self,
        gene: GeneNode,
        tissue: TissueNode,
        coloc_h4: Optional[float] = None,
        expression_tpm: Optional[float] = None,
        replicated_in_eqtl_catalogue: Optional[bool] = None,
    ) -> EdgeProbability:
        """
        Compute gene → tissue edge probability.
        
        Parameters
        ----------
        gene : GeneNode
            Source gene.
        tissue : TissueNode
            Target tissue.
        coloc_h4 : float, optional
            PP.H4 from colocalization.
        expression_tpm : float, optional
            Gene expression (TPM) in tissue.
        replicated_in_eqtl_catalogue : bool, optional
            Whether eQTL replicates in eQTL Catalogue.
            
        Returns
        -------
        EdgeProbability
            Edge with computed probability.
        """
        evidence = {}
        
        # Start with tissue relevance prior
        prob = tissue.relevance_prior * 0.3
        evidence["relevance_prior"] = tissue.relevance_prior
        
        # Gene expression contribution
        if expression_tpm is not None:
            # Log-transform and normalize
            expr_score = np.log1p(expression_tpm) / 10
            expr_score = min(expr_score, 1.0)
            prob += expr_score * 0.2
            evidence["expression_tpm"] = expression_tpm
            evidence["expression_score"] = expr_score
        elif gene.is_expressed(tissue.tissue_id):
            prob += 0.15
            evidence["is_expressed"] = True
        
        # Colocalization is strongest evidence
        if coloc_h4 is not None:
            # Weight coloc more heavily
            prob += coloc_h4 * 0.4
            evidence["coloc_h4"] = coloc_h4
            
            # Boost if replicated in eQTL Catalogue
            if replicated_in_eqtl_catalogue is True:
                prob *= 1.2  # 20% boost for cross-study replication
                evidence["eqtl_catalogue_replicated"] = True
            elif replicated_in_eqtl_catalogue is False:
                # Slight penalty for non-replication (but not too harsh)
                prob *= 0.9
                evidence["eqtl_catalogue_replicated"] = False
        
        return EdgeProbability(
            source=gene,
            target=tissue,
            probability=min(prob, 1.0),
            evidence=evidence,
            edge_type="gene_active_in_tissue",
            method="coloc_expression" if coloc_h4 else "expression_only",
        )
    
    @staticmethod
    def compute_static(
        gene: GeneNode,
        tissue: TissueNode,
        coloc_h4: Optional[float] = None,
        expression_tpm: Optional[float] = None,
    ) -> EdgeProbability:
        """
        Static method for backward compatibility.
        
        Parameters
        ----------
        gene : GeneNode
            Source gene.
        tissue : TissueNode
            Target tissue.
        coloc_h4 : float, optional
            PP.H4 from colocalization.
        expression_tpm : float, optional
            Gene expression (TPM) in tissue.
            
        Returns
        -------
        EdgeProbability
            Edge with computed probability.
        """
        edge_computer = GeneToTissueEdge()
        return edge_computer.compute(gene, tissue, coloc_h4, expression_tpm)


class TissueToTraitEdge:
    """
    Computes probability that tissue contributes to trait.
    
    Evidence sources:
    - Prior tissue-trait relevance
    - Enrichment of GWAS signal in tissue-active elements
    - LD score regression tissue enrichment
    """
    
    @staticmethod
    def compute(
        tissue: TissueNode,
        trait: TraitNode,
        enrichment_score: Optional[float] = None,
        ldsc_coefficient: Optional[float] = None,
    ) -> EdgeProbability:
        """
        Compute tissue → trait edge probability.
        
        Parameters
        ----------
        tissue : TissueNode
            Source tissue.
        trait : TraitNode
            Target trait.
        enrichment_score : float, optional
            GWAS enrichment in tissue elements.
        ldsc_coefficient : float, optional
            LD score regression coefficient.
            
        Returns
        -------
        EdgeProbability
            Edge with computed probability.
        """
        evidence = {}
        
        # Base on tissue relevance prior
        prob = tissue.relevance_prior
        evidence["relevance_prior"] = tissue.relevance_prior
        
        # Enrichment evidence
        if enrichment_score is not None:
            # Convert to probability
            enrich_prob = 1 - np.exp(-enrichment_score)
            prob = prob * 0.5 + enrich_prob * 0.5
            evidence["enrichment_score"] = enrichment_score
        
        # LDSC evidence
        if ldsc_coefficient is not None and ldsc_coefficient > 0:
            ldsc_prob = min(ldsc_coefficient / 2, 1.0)
            prob = prob * 0.6 + ldsc_prob * 0.4
            evidence["ldsc_coefficient"] = ldsc_coefficient
        
        return EdgeProbability(
            source=tissue,
            target=trait,
            probability=min(prob, 1.0),
            evidence=evidence,
            edge_type="tissue_contributes_to_trait",
        )


def compute_edge_weights(
    edge_type: str,
    source: BaseNode,
    target: BaseNode,
    **kwargs,
) -> EdgeProbability:
    """
    Factory function to compute edge weights.
    
    Parameters
    ----------
    edge_type : str
        Type of edge.
    source : BaseNode
        Source node.
    target : BaseNode
        Target node.
    **kwargs
        Additional evidence.
        
    Returns
    -------
    EdgeProbability
        Computed edge.
    """
    edge_computers = {
        "variant_to_ccre": VariantToCCREEdge.compute,
        "ccre_to_gene": CCREToGeneEdge.compute_distance_only,  # Default to distance
        "gene_to_tissue": GeneToTissueEdge.compute_static,
        "tissue_to_trait": TissueToTraitEdge.compute,
    }
    
    if edge_type not in edge_computers:
        raise ValueError(f"Unknown edge type: {edge_type}")
    
    return edge_computers[edge_type](source, target, **kwargs)


def chain_probabilities(
    edges: List[EdgeProbability],
) -> float:
    """
    Compute probability along a chain of edges.
    
    Uses product of probabilities (independence assumption).
    
    Parameters
    ----------
    edges : list
        List of edges in the chain.
        
    Returns
    -------
    float
        Chain probability.
    """
    if not edges:
        return 0.0
    
    prob = 1.0
    for edge in edges:
        prob *= edge.probability
    
    return prob


def max_path_probability(
    all_paths: List[List[EdgeProbability]],
) -> Tuple[float, List[EdgeProbability]]:
    """
    Find the maximum probability path.
    
    Parameters
    ----------
    all_paths : list
        List of paths, each path is a list of edges.
        
    Returns
    -------
    tuple
        (max_probability, best_path)
    """
    if not all_paths:
        return 0.0, []
    
    max_prob = 0.0
    best_path = []
    
    for path in all_paths:
        prob = chain_probabilities(path)
        if prob > max_prob:
            max_prob = prob
            best_path = path
    
    return max_prob, best_path


def create_edge_computer_with_linkers(
    abc_linker: Optional[Any] = None,
    pchic_linker: Optional[Any] = None,
    ensemble_linker: Optional[Any] = None,
    replication_tester: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Create edge computers with functional linkers.
    
    Parameters
    ----------
    abc_linker : ABCLinks, optional
        ABC model linker from enhancer_gene_linking module.
    pchic_linker : PCHiCLinks, optional  
        PCHi-C linker from enhancer_gene_linking module.
    ensemble_linker : EnsembleLinker, optional
        Ensemble linker from enhancer_gene_linking module.
    replication_tester : CrossStudyReplication, optional
        For eQTL cross-study replication.
        
    Returns
    -------
    dict
        Edge computers keyed by edge type.
        
    Example
    -------
    >>> from ..enhancer_gene_linking import ABCLinks, PCHiCLinks, EnsembleLinker
    >>> 
    >>> abc = ABCLinks()
    >>> abc.load_abc_predictions("abc_predictions.tsv.gz")
    >>> 
    >>> pchic = PCHiCLinks()
    >>> pchic.load_pchic_data("pchic_contacts.tsv.gz")
    >>> 
    >>> ensemble = EnsembleLinker(abc_linker=abc, pchic_linker=pchic)
    >>> 
    >>> edge_computers = create_edge_computer_with_linkers(
    ...     abc_linker=abc,
    ...     pchic_linker=pchic,
    ...     ensemble_linker=ensemble,
    ... )
    >>> 
    >>> # Use ensemble linker for cCRE→gene edges
    >>> ccre_gene_edge = edge_computers["ccre_to_gene"].compute(ccre, gene, cell_type="K562")
    """
    ccre_gene = CCREToGeneEdge(
        abc_linker=abc_linker,
        pchic_linker=pchic_linker,
        ensemble_linker=ensemble_linker,
        method="ensemble" if ensemble_linker else ("abc" if abc_linker else "distance"),
    )
    
    gene_tissue = GeneToTissueEdge(
        replication_tester=replication_tester,
    )
    
    return {
        "variant_to_ccre": VariantToCCREEdge(),
        "ccre_to_gene": ccre_gene,
        "gene_to_tissue": gene_tissue,
        "tissue_to_trait": TissueToTraitEdge(),
    }
