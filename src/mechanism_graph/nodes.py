"""
Node Definitions for Mechanism Graphs

Defines the five node types in the causal chain:
1. Variant - genetic variant with fine-mapping PIP
2. cCRE - candidate cis-regulatory element
3. Gene - target gene with evidence
4. Tissue - relevant tissue context
5. Trait - the phenotype of interest
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class BaseNode:
    """Base class for all nodes."""
    
    id: str
    type: str
    attributes: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "id": self.id,
            "type": self.type,
            "attributes": self.attributes,
        }
    
    def __hash__(self):
        return hash(self.id)
    
    def __eq__(self, other):
        if isinstance(other, BaseNode):
            return self.id == other.id
        return False


@dataclass
class VariantNode(BaseNode):
    """
    Genetic variant node.
    
    Attributes
    ----------
    rsid : str
        Variant identifier.
    chr : str
        Chromosome.
    pos : int
        Genomic position (GRCh38).
    ref : str
        Reference allele.
    alt : str
        Alternative allele.
    pip : float
        Posterior inclusion probability from fine-mapping.
    cs_id : int, optional
        Credible set ID.
    z_score : float, optional
        Z-score from GWAS.
    """
    
    type: str = "variant"
    rsid: str = ""
    chr: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    pip: float = 0.0
    cs_id: Optional[int] = None
    z_score: Optional[float] = None
    
    def __post_init__(self):
        if not self.id:
            self.id = self.rsid or f"chr{self.chr}_{self.pos}_{self.ref}_{self.alt}"
        
        self.attributes.update({
            "rsid": self.rsid,
            "chr": self.chr,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "pip": self.pip,
            "cs_id": self.cs_id,
            "z_score": self.z_score,
        })
    
    @classmethod
    def from_finemapping(cls, fm_result: Dict) -> "VariantNode":
        """Create from fine-mapping result."""
        return cls(
            id=fm_result.get("rsid", ""),
            rsid=fm_result.get("rsid", ""),
            chr=str(fm_result.get("chr", "")),
            pos=fm_result.get("pos", 0),
            ref=fm_result.get("ref", ""),
            alt=fm_result.get("alt", ""),
            pip=fm_result.get("pip", 0.0),
            cs_id=fm_result.get("cs_id"),
            z_score=fm_result.get("z_score"),
        )


@dataclass
class CCRENode(BaseNode):
    """
    Candidate cis-Regulatory Element node.
    
    Attributes
    ----------
    ccre_id : str
        ENCODE cCRE accession.
    chr : str
        Chromosome.
    start : int
        Start position.
    end : int
        End position.
    element_type : str
        Element type (PLS, ELS, CTCF-only, etc.).
    activity_score : float
        Activity score in relevant context.
    """
    
    type: str = "ccre"
    ccre_id: str = ""
    chr: str = ""
    start: int = 0
    end: int = 0
    element_type: str = ""
    activity_score: float = 0.0
    tissue_activities: Dict[str, float] = field(default_factory=dict)
    
    def __post_init__(self):
        if not self.id:
            self.id = self.ccre_id or f"ccre_{self.chr}_{self.start}_{self.end}"
        
        self.attributes.update({
            "ccre_id": self.ccre_id,
            "chr": self.chr,
            "start": self.start,
            "end": self.end,
            "element_type": self.element_type,
            "activity_score": self.activity_score,
            "tissue_activities": self.tissue_activities,
        })
    
    def contains_position(self, pos: int) -> bool:
        """Check if position falls within this cCRE."""
        return self.start <= pos <= self.end
    
    def get_tissue_activity(self, tissue: str) -> float:
        """Get activity score for a specific tissue."""
        return self.tissue_activities.get(tissue, 0.0)


@dataclass
class GeneNode(BaseNode):
    """
    Gene node.
    
    Attributes
    ----------
    gene_id : str
        Ensembl gene ID (without version).
    symbol : str
        Gene symbol.
    chr : str
        Chromosome.
    tss : int
        Transcription start site.
    strand : str
        Strand (+ or -).
    biotype : str
        Gene biotype (protein_coding, etc.).
    """
    
    type: str = "gene"
    gene_id: str = ""
    symbol: str = ""
    chr: str = ""
    tss: int = 0
    strand: str = "+"
    biotype: str = "protein_coding"
    expression: Dict[str, float] = field(default_factory=dict)
    
    def __post_init__(self):
        if not self.id:
            self.id = self.gene_id or self.symbol
        
        self.attributes.update({
            "gene_id": self.gene_id,
            "symbol": self.symbol,
            "chr": self.chr,
            "tss": self.tss,
            "strand": self.strand,
            "biotype": self.biotype,
            "expression": self.expression,
        })
    
    def get_expression(self, tissue: str) -> float:
        """Get expression in a tissue."""
        return self.expression.get(tissue, 0.0)
    
    def is_expressed(self, tissue: str, threshold: float = 1.0) -> bool:
        """Check if gene is expressed in a tissue."""
        return self.get_expression(tissue) >= threshold


@dataclass
class TissueNode(BaseNode):
    """
    Tissue context node.
    
    Attributes
    ----------
    tissue_id : str
        Tissue identifier.
    name : str
        Human-readable name.
    ontology_id : str
        UBERON or CL ontology ID.
    relevance_prior : float
        Prior probability of relevance to the trait.
    """
    
    type: str = "tissue"
    tissue_id: str = ""
    name: str = ""
    ontology_id: str = ""
    relevance_prior: float = 0.5
    
    def __post_init__(self):
        if not self.id:
            self.id = self.tissue_id or self.name
        
        self.attributes.update({
            "tissue_id": self.tissue_id,
            "name": self.name,
            "ontology_id": self.ontology_id,
            "relevance_prior": self.relevance_prior,
        })


@dataclass
class TraitNode(BaseNode):
    """
    Trait/phenotype node.
    
    Attributes
    ----------
    trait_id : str
        Trait identifier.
    name : str
        Trait name.
    category : str
        Trait category (cardiometabolic, etc.).
    ontology_id : str
        EFO or other ontology ID.
    """
    
    type: str = "trait"
    trait_id: str = ""
    name: str = ""
    category: str = ""
    ontology_id: str = ""
    
    def __post_init__(self):
        if not self.id:
            self.id = self.trait_id or self.name
        
        self.attributes.update({
            "trait_id": self.trait_id,
            "name": self.name,
            "category": self.category,
            "ontology_id": self.ontology_id,
        })


def create_node(node_type: str, **kwargs) -> BaseNode:
    """
    Factory function to create nodes.
    
    Parameters
    ----------
    node_type : str
        Type of node to create.
    **kwargs
        Node attributes.
        
    Returns
    -------
    BaseNode
        Created node.
    """
    node_classes = {
        "variant": VariantNode,
        "ccre": CCRENode,
        "gene": GeneNode,
        "tissue": TissueNode,
        "trait": TraitNode,
    }
    
    if node_type not in node_classes:
        raise ValueError(f"Unknown node type: {node_type}")
    
    return node_classes[node_type](**kwargs)
