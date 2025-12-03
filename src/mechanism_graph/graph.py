"""
Mechanism Graph Core

Main graph structure and builder for mechanism graphs.
"""

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Set, Tuple

import json
import numpy as np

from ..utils.logging import get_logger
from .nodes import (
    BaseNode, VariantNode, CCRENode, GeneNode, TissueNode, TraitNode,
    create_node
)
from .edges import (
    EdgeProbability, VariantToCCREEdge, CCREToGeneEdge,
    GeneToTissueEdge, TissueToTraitEdge, chain_probabilities
)


logger = get_logger("mechanism_graph")


class MechanismGraph:
    """
    Probabilistic mechanism graph linking variants to traits.
    
    Structure:
    Variants → cCREs → Genes → Tissues → Traits
    
    Each edge has a calibrated probability representing
    the confidence that the relationship is causal.
    """
    
    def __init__(self, trait_id: str):
        """
        Initialize mechanism graph.
        
        Parameters
        ----------
        trait_id : str
            Identifier for the trait.
        """
        self.trait_id = trait_id
        
        # Node storage by type
        self.nodes: Dict[str, Dict[str, BaseNode]] = {
            "variant": {},
            "ccre": {},
            "gene": {},
            "tissue": {},
            "trait": {},
        }
        
        # Edge storage
        self.edges: List[EdgeProbability] = []
        
        # Adjacency lists for traversal
        self._outgoing: Dict[str, List[EdgeProbability]] = defaultdict(list)
        self._incoming: Dict[str, List[EdgeProbability]] = defaultdict(list)
    
    def add_node(self, node: BaseNode) -> None:
        """Add a node to the graph."""
        self.nodes[node.type][node.id] = node
    
    def add_edge(self, edge: EdgeProbability) -> None:
        """Add an edge to the graph."""
        # Ensure nodes exist
        if edge.source.id not in self.nodes[edge.source.type]:
            self.add_node(edge.source)
        if edge.target.id not in self.nodes[edge.target.type]:
            self.add_node(edge.target)
        
        self.edges.append(edge)
        self._outgoing[edge.source.id].append(edge)
        self._incoming[edge.target.id].append(edge)
    
    def get_node(self, node_id: str, node_type: Optional[str] = None) -> Optional[BaseNode]:
        """Get a node by ID."""
        if node_type:
            return self.nodes[node_type].get(node_id)
        
        for type_nodes in self.nodes.values():
            if node_id in type_nodes:
                return type_nodes[node_id]
        return None
    
    def get_edges_from(self, node_id: str) -> List[EdgeProbability]:
        """Get all outgoing edges from a node."""
        return self._outgoing[node_id]
    
    def get_edges_to(self, node_id: str) -> List[EdgeProbability]:
        """Get all incoming edges to a node."""
        return self._incoming[node_id]
    
    def find_all_paths(
        self,
        source_id: str,
        target_id: str,
        max_depth: int = 5,
    ) -> List[List[EdgeProbability]]:
        """
        Find all paths between two nodes.
        
        Parameters
        ----------
        source_id : str
            Source node ID.
        target_id : str
            Target node ID.
        max_depth : int
            Maximum path length.
            
        Returns
        -------
        list
            List of paths (each path is a list of edges).
        """
        paths = []
        
        def dfs(current: str, path: List[EdgeProbability], visited: Set[str]):
            if current == target_id:
                paths.append(path.copy())
                return
            
            if len(path) >= max_depth:
                return
            
            for edge in self._outgoing[current]:
                if edge.target.id not in visited:
                    visited.add(edge.target.id)
                    path.append(edge)
                    dfs(edge.target.id, path, visited)
                    path.pop()
                    visited.discard(edge.target.id)
        
        dfs(source_id, [], {source_id})
        return paths
    
    def compute_variant_to_trait_probability(
        self,
        variant_id: str,
        aggregation: str = "max",
    ) -> float:
        """
        Compute probability that a variant affects the trait.
        
        Parameters
        ----------
        variant_id : str
            Variant ID.
        aggregation : str
            How to aggregate paths: "max", "sum", or "noisy_or".
            
        Returns
        -------
        float
            Probability (0-1).
        """
        trait_nodes = list(self.nodes["trait"].keys())
        
        if not trait_nodes:
            return 0.0
        
        all_probs = []
        
        for trait_id in trait_nodes:
            paths = self.find_all_paths(variant_id, trait_id)
            
            for path in paths:
                prob = chain_probabilities(path)
                all_probs.append(prob)
        
        if not all_probs:
            return 0.0
        
        if aggregation == "max":
            return max(all_probs)
        elif aggregation == "sum":
            return min(sum(all_probs), 1.0)
        elif aggregation == "noisy_or":
            # P(at least one) = 1 - P(none)
            return 1 - np.prod([1 - p for p in all_probs])
        else:
            return max(all_probs)
    
    def get_gene_scores(self) -> Dict[str, Dict[str, Any]]:
        """
        Compute gene prioritization scores.
        
        Returns
        -------
        dict
            Gene scores with evidence.
        """
        gene_scores = {}
        
        for gene_id, gene in self.nodes["gene"].items():
            # Get all paths through this gene
            incoming = self.get_edges_to(gene_id)
            outgoing = self.get_edges_from(gene_id)
            
            # Compute upstream probability (variant → gene)
            upstream_prob = 0.0
            variant_evidence = []
            
            for edge in incoming:
                if edge.source.type == "ccre":
                    # Find variants affecting this cCRE
                    ccre_edges = self.get_edges_to(edge.source.id)
                    for var_edge in ccre_edges:
                        if var_edge.source.type == "variant":
                            path_prob = var_edge.probability * edge.probability
                            upstream_prob = max(upstream_prob, path_prob)
                            variant_evidence.append({
                                "variant": var_edge.source.id,
                                "ccre": edge.source.id,
                                "probability": path_prob,
                            })
            
            # Compute downstream probability (gene → trait)
            downstream_prob = 0.0
            tissue_evidence = []
            
            for edge in outgoing:
                if edge.source.type == "gene" and edge.target.type == "tissue":
                    tissue_edges = self.get_edges_from(edge.target.id)
                    for trait_edge in tissue_edges:
                        if trait_edge.target.type == "trait":
                            path_prob = edge.probability * trait_edge.probability
                            downstream_prob = max(downstream_prob, path_prob)
                            tissue_evidence.append({
                                "tissue": edge.target.id,
                                "trait": trait_edge.target.id,
                                "probability": path_prob,
                            })
            
            # Combined score
            combined = upstream_prob * downstream_prob
            
            gene_scores[gene_id] = {
                "gene_id": gene_id,
                "symbol": gene.symbol,
                "combined_score": combined,
                "upstream_probability": upstream_prob,
                "downstream_probability": downstream_prob,
                "variant_evidence": variant_evidence,
                "tissue_evidence": tissue_evidence,
            }
        
        return gene_scores
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert graph to dictionary for JSON serialization."""
        return {
            "trait_id": self.trait_id,
            "nodes": {
                node_type: {
                    node_id: node.to_dict()
                    for node_id, node in type_nodes.items()
                }
                for node_type, type_nodes in self.nodes.items()
            },
            "edges": [edge.to_dict() for edge in self.edges],
            "statistics": {
                "n_variants": len(self.nodes["variant"]),
                "n_ccres": len(self.nodes["ccre"]),
                "n_genes": len(self.nodes["gene"]),
                "n_tissues": len(self.nodes["tissue"]),
                "n_traits": len(self.nodes["trait"]),
                "n_edges": len(self.edges),
            },
        }
    
    def save(self, path: str) -> None:
        """Save graph to JSON file."""
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
    
    @classmethod
    def load(cls, path: str) -> "MechanismGraph":
        """Load graph from JSON file."""
        with open(path) as f:
            data = json.load(f)
        
        graph = cls(data["trait_id"])
        
        # Reconstruct nodes
        for node_type, type_nodes in data["nodes"].items():
            for node_id, node_data in type_nodes.items():
                node = create_node(node_type, id=node_id, **node_data.get("attributes", {}))
                graph.add_node(node)
        
        # Reconstruct edges
        for edge_data in data["edges"]:
            source = graph.get_node(edge_data["source"], edge_data["source_type"])
            target = graph.get_node(edge_data["target"], edge_data["target_type"])
            
            if source and target:
                edge = EdgeProbability(
                    source=source,
                    target=target,
                    probability=edge_data["probability"],
                    evidence=edge_data.get("evidence", {}),
                    edge_type=edge_data["edge_type"],
                    confidence=edge_data.get("confidence", 1.0),
                )
                graph.edges.append(edge)
                graph._outgoing[source.id].append(edge)
                graph._incoming[target.id].append(edge)
        
        return graph


class MechanismGraphBuilder:
    """
    Builder for constructing mechanism graphs from analysis results.
    """
    
    def __init__(
        self,
        trait_name: str,
        trait_id: str,
    ):
        """
        Initialize builder.
        
        Parameters
        ----------
        trait_name : str
            Human-readable trait name.
        trait_id : str
            Trait identifier.
        """
        self.trait_name = trait_name
        self.trait_id = trait_id
        self.graph = MechanismGraph(trait_id)
        
        # Add trait node
        trait = TraitNode(
            trait_id=trait_id,
            name=trait_name,
            category="cardiometabolic",
        )
        self.graph.add_node(trait)
    
    def add_variants_from_finemapping(
        self,
        finemapping_results: List[Dict],
        pip_threshold: float = 0.01,
    ) -> None:
        """
        Add variant nodes from fine-mapping results.
        
        Parameters
        ----------
        finemapping_results : list
            Fine-mapping results with PIPs.
        pip_threshold : float
            Minimum PIP to include.
        """
        for result in finemapping_results:
            pip_dict = result.get("pip", {})
            
            for rsid, pip in pip_dict.items():
                if pip >= pip_threshold:
                    # Find variant info
                    variant_info = next(
                        (v for v in result.get("variants", []) if v.get("rsid") == rsid),
                        {"rsid": rsid}
                    )
                    
                    variant = VariantNode(
                        rsid=rsid,
                        chr=str(variant_info.get("chr", "")),
                        pos=variant_info.get("pos", 0),
                        ref=variant_info.get("ref", ""),
                        alt=variant_info.get("alt", ""),
                        pip=pip,
                    )
                    
                    self.graph.add_node(variant)
        
        logger.info(f"Added {len(self.graph.nodes['variant'])} variants")
    
    def add_ccres_from_annotations(
        self,
        ccre_annotations: List[Dict],
    ) -> None:
        """
        Add cCRE nodes from annotations.
        
        Parameters
        ----------
        ccre_annotations : list
            cCRE annotations with coordinates.
        """
        for ccre_data in ccre_annotations:
            ccre = CCRENode(
                ccre_id=ccre_data.get("id", ""),
                chr=str(ccre_data.get("chr", "")),
                start=ccre_data.get("start", 0),
                end=ccre_data.get("end", 0),
                element_type=ccre_data.get("type", ""),
                activity_score=ccre_data.get("activity", 0.0),
            )
            
            self.graph.add_node(ccre)
        
        logger.info(f"Added {len(self.graph.nodes['ccre'])} cCREs")
    
    def add_genes(
        self,
        genes: List[Dict],
    ) -> None:
        """
        Add gene nodes.
        
        Parameters
        ----------
        genes : list
            Gene information.
        """
        for gene_data in genes:
            gene = GeneNode(
                gene_id=gene_data.get("gene_id", ""),
                symbol=gene_data.get("symbol", ""),
                chr=str(gene_data.get("chr", "")),
                tss=gene_data.get("tss", 0),
                strand=gene_data.get("strand", "+"),
                biotype=gene_data.get("biotype", "protein_coding"),
            )
            
            self.graph.add_node(gene)
        
        logger.info(f"Added {len(self.graph.nodes['gene'])} genes")
    
    def add_tissues(
        self,
        tissues: List[Dict],
    ) -> None:
        """
        Add tissue nodes with relevance priors.
        
        Parameters
        ----------
        tissues : list
            Tissue information with priors.
        """
        for tissue_data in tissues:
            tissue = TissueNode(
                tissue_id=tissue_data.get("tissue_id", ""),
                name=tissue_data.get("name", ""),
                ontology_id=tissue_data.get("ontology_id", ""),
                relevance_prior=tissue_data.get("prior", 0.5),
            )
            
            self.graph.add_node(tissue)
        
        logger.info(f"Added {len(self.graph.nodes['tissue'])} tissues")
    
    def connect_variants_to_ccres(self) -> None:
        """Connect variants to overlapping/nearby cCREs."""
        for variant in self.graph.nodes["variant"].values():
            for ccre in self.graph.nodes["ccre"].values():
                if str(variant.chr) != str(ccre.chr):
                    continue
                
                edge = VariantToCCREEdge.compute(variant, ccre)
                
                if edge.probability > 0.01:
                    self.graph.add_edge(edge)
        
        logger.info(f"Added {len([e for e in self.graph.edges if e.edge_type == 'variant_affects_ccre'])} variant→cCRE edges")
    
    def connect_ccres_to_genes(
        self,
        abc_scores: Optional[Dict[Tuple[str, str], float]] = None,
    ) -> None:
        """
        Connect cCREs to target genes.
        
        Parameters
        ----------
        abc_scores : dict, optional
            ABC scores keyed by (ccre_id, gene_id).
        """
        for ccre in self.graph.nodes["ccre"].values():
            for gene in self.graph.nodes["gene"].values():
                if str(ccre.chr) != str(gene.chr):
                    continue
                
                abc = None
                if abc_scores and (ccre.ccre_id, gene.gene_id) in abc_scores:
                    abc = abc_scores[(ccre.ccre_id, gene.gene_id)]
                
                edge = CCREToGeneEdge.compute(ccre, gene, abc_score=abc)
                
                if edge.probability > 0.01:
                    self.graph.add_edge(edge)
        
        logger.info(f"Added {len([e for e in self.graph.edges if e.edge_type == 'ccre_regulates_gene'])} cCRE→gene edges")
    
    def connect_genes_to_tissues_from_coloc(
        self,
        coloc_results: List[Dict],
    ) -> None:
        """
        Connect genes to tissues based on colocalization.
        
        Parameters
        ----------
        coloc_results : list
            Colocalization results.
        """
        for result in coloc_results:
            for coloc in result.get("colocalizations", []):
                gene_id = coloc.get("gene")
                tissue_id = coloc.get("tissue")
                h4 = coloc.get("pp_h4", 0)
                
                gene = self.graph.get_node(gene_id, "gene")
                tissue = self.graph.get_node(tissue_id, "tissue")
                
                if gene and tissue:
                    edge = GeneToTissueEdge.compute(gene, tissue, coloc_h4=h4)
                    self.graph.add_edge(edge)
        
        logger.info(f"Added {len([e for e in self.graph.edges if e.edge_type == 'gene_active_in_tissue'])} gene→tissue edges")
    
    def connect_tissues_to_trait(self) -> None:
        """Connect all tissues to the trait node."""
        trait = list(self.graph.nodes["trait"].values())[0]
        
        for tissue in self.graph.nodes["tissue"].values():
            edge = TissueToTraitEdge.compute(tissue, trait)
            self.graph.add_edge(edge)
        
        logger.info(f"Added {len([e for e in self.graph.edges if e.edge_type == 'tissue_contributes_to_trait'])} tissue→trait edges")
    
    def build(self) -> MechanismGraph:
        """
        Build and return the complete graph.
        
        Returns
        -------
        MechanismGraph
            Constructed graph.
        """
        logger.info(
            f"Built mechanism graph: "
            f"{len(self.graph.nodes['variant'])} variants, "
            f"{len(self.graph.nodes['ccre'])} cCREs, "
            f"{len(self.graph.nodes['gene'])} genes, "
            f"{len(self.graph.nodes['tissue'])} tissues, "
            f"{len(self.graph.edges)} edges"
        )
        
        return self.graph
