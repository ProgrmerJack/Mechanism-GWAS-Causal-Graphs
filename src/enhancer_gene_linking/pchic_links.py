"""
Promoter Capture Hi-C (PCHi-C) Links

Loads and queries PCHi-C enhancer-gene predictions from:
- Jung et al. 2019 (17 human primary blood cell types)
- Javierre et al. 2016 (17 human primary blood cell types)
- EpiMap/ENCODE PCHi-C predictions

PCHi-C directly measures chromatin contacts between promoters and distal
elements, providing orthogonal evidence to ABC model.

Key data sources:
- CHiCAGO-processed PCHi-C from Blueprint Epigenome
- 3DIV database (3D Interaction Viewer)
- ENCODE CRISPR-validated chromatin contacts
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from ..utils.logging import get_logger


logger = get_logger("pchic_links")


# CHiCAGO score thresholds
CHICAGO_THRESHOLDS = {
    "significant": 5.0,     # Standard CHiCAGO threshold
    "stringent": 10.0,      # High confidence
    "permissive": 3.0,      # For discovery
}


@dataclass
class PCHiCLink:
    """
    A single PCHi-C interaction.
    
    Attributes
    ----------
    bait_chr : str
        Bait (promoter) chromosome.
    bait_start : int
        Bait start position.
    bait_end : int
        Bait end position.
    oe_chr : str
        Other end (enhancer) chromosome.
    oe_start : int
        Other end start position.
    oe_end : int  
        Other end end position.
    gene_symbol : str
        Target gene at bait.
    gene_id : str
        Ensembl gene ID.
    chicago_score : float
        CHiCAGO interaction score.
    cell_type : str
        Cell type of interaction.
    n_reads : int
        Number of supporting read pairs.
    distance : int
        Interaction distance in bp.
    is_validated : bool
        CRISPR validation status.
    """
    
    bait_chr: str
    bait_start: int
    bait_end: int
    oe_chr: str
    oe_start: int
    oe_end: int
    gene_symbol: str
    gene_id: str = ""
    chicago_score: float = 0.0
    cell_type: str = ""
    n_reads: int = 0
    distance: int = 0
    is_validated: bool = False
    
    @property
    def bait_id(self) -> str:
        """Bait (promoter) identifier."""
        return f"{self.bait_chr}:{self.bait_start}-{self.bait_end}"
    
    @property
    def oe_id(self) -> str:
        """Other end (enhancer) identifier."""
        return f"{self.oe_chr}:{self.oe_start}-{self.oe_end}"
    
    @property
    def link_id(self) -> str:
        """Unique interaction identifier."""
        return f"{self.bait_id}:{self.oe_id}:{self.cell_type}"
    
    def overlaps_region(
        self,
        chrom: str,
        start: int,
        end: int,
        side: str = "oe",
    ) -> bool:
        """Check if interaction overlaps a region."""
        if side == "oe":
            if chrom.replace("chr", "") != self.oe_chr.replace("chr", ""):
                return False
            return not (end < self.oe_start or start > self.oe_end)
        else:
            if chrom.replace("chr", "") != self.bait_chr.replace("chr", ""):
                return False
            return not (end < self.bait_start or start > self.bait_end)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "bait_chr": self.bait_chr,
            "bait_start": self.bait_start,
            "bait_end": self.bait_end,
            "bait_id": self.bait_id,
            "oe_chr": self.oe_chr,
            "oe_start": self.oe_start,
            "oe_end": self.oe_end,
            "oe_id": self.oe_id,
            "gene_symbol": self.gene_symbol,
            "gene_id": self.gene_id,
            "chicago_score": self.chicago_score,
            "cell_type": self.cell_type,
            "n_reads": self.n_reads,
            "distance": self.distance,
            "is_validated": self.is_validated,
        }


class PCHiCLinks:
    """
    Manager for PCHi-C enhancer-gene interactions.
    
    Loads and queries chromatin contact data from CHiCAGO-processed
    PCHi-C experiments.
    
    Example
    -------
    >>> pchic = PCHiCLinks()
    >>> pchic.load_blueprint("/path/to/pchic_calls.txt.gz")
    >>> links = pchic.get_links_for_region("chr1", 1000000, 1500000)
    >>> gene_links = pchic.get_links_for_gene("LDLR")
    """
    
    # Column mappings for different file formats
    COLUMN_MAPPINGS = {
        "chicago": {
            "baitChr": "bait_chr",
            "baitStart": "bait_start",
            "baitEnd": "bait_end",
            "oeChr": "oe_chr",
            "oeStart": "oe_start",
            "oeEnd": "oe_end",
            "baitName": "gene_symbol",
            "score": "chicago_score",
            "N_reads": "n_reads",
        },
        "3div": {
            "chrA": "bait_chr",
            "startA": "bait_start",
            "endA": "bait_end",
            "chrB": "oe_chr",
            "startB": "oe_start",
            "endB": "oe_end",
            "gene": "gene_symbol",
            "score": "chicago_score",
        },
    }
    
    def __init__(
        self,
        chicago_threshold: float = 5.0,
        cell_types: Optional[List[str]] = None,
    ):
        """
        Initialize PCHi-C links manager.
        
        Parameters
        ----------
        chicago_threshold : float
            Minimum CHiCAGO score threshold (default: 5.0).
        cell_types : list, optional
            Filter to specific cell types.
        """
        self.chicago_threshold = chicago_threshold
        self.cell_type_filter = set(cell_types) if cell_types else None
        
        # Storage
        self._links: List[PCHiCLink] = []
        
        # Indices
        self._by_oe_region: Dict[str, List[int]] = {}   # chr -> indices
        self._by_bait_region: Dict[str, List[int]] = {}  # chr -> indices  
        self._by_gene: Dict[str, List[int]] = {}         # gene -> indices
        self._by_cell_type: Dict[str, List[int]] = {}    # cell_type -> indices
    
    def load_chicago_file(
        self,
        filepath: str,
        cell_type: str = "",
        file_format: str = "chicago",
    ) -> int:
        """
        Load CHiCAGO-processed PCHi-C interactions.
        
        Parameters
        ----------
        filepath : str
            Path to CHiCAGO output file.
        cell_type : str
            Cell type name (if not in file).
        file_format : str
            File format: "chicago" or "3div".
            
        Returns
        -------
        int
            Number of links loaded.
        """
        logger.info(f"Loading PCHi-C from {filepath}")
        
        df = pd.read_csv(filepath, sep="\t", compression="infer")
        col_map = self.COLUMN_MAPPINGS.get(file_format, self.COLUMN_MAPPINGS["chicago"])
        
        count = 0
        for _, row in df.iterrows():
            score = float(row.get(col_map.get("score", "score"), 0))
            
            if score < self.chicago_threshold:
                continue
            
            # Get cell type from file or parameter
            ct = cell_type or str(row.get("cellType", row.get("cell_type", "")))
            
            if self.cell_type_filter and ct not in self.cell_type_filter:
                continue
            
            # Parse gene from bait name (may be "GENE:ENS..." format)
            bait_name = str(row.get(col_map.get("baitName", "baitName"), ""))
            gene = bait_name.split(":")[0] if ":" in bait_name else bait_name
            gene_id = bait_name.split(":")[1] if ":" in bait_name else ""
            
            # Calculate distance
            bait_mid = (int(row.get(col_map.get("bait_start", "baitStart"), 0)) + 
                       int(row.get(col_map.get("bait_end", "baitEnd"), 0))) // 2
            oe_mid = (int(row.get(col_map.get("oe_start", "oeStart"), 0)) + 
                     int(row.get(col_map.get("oe_end", "oeEnd"), 0))) // 2
            
            link = PCHiCLink(
                bait_chr=str(row.get(col_map.get("bait_chr", "baitChr"), "")),
                bait_start=int(row.get(col_map.get("bait_start", "baitStart"), 0)),
                bait_end=int(row.get(col_map.get("bait_end", "baitEnd"), 0)),
                oe_chr=str(row.get(col_map.get("oe_chr", "oeChr"), "")),
                oe_start=int(row.get(col_map.get("oe_start", "oeStart"), 0)),
                oe_end=int(row.get(col_map.get("oe_end", "oeEnd"), 0)),
                gene_symbol=gene,
                gene_id=gene_id,
                chicago_score=score,
                cell_type=ct,
                n_reads=int(row.get(col_map.get("n_reads", "N_reads"), 0)),
                distance=abs(bait_mid - oe_mid),
            )
            
            self._add_link(link)
            count += 1
        
        logger.info(f"Loaded {count} PCHi-C links above threshold {self.chicago_threshold}")
        return count
    
    def load_blueprint(
        self,
        filepath: str,
        cell_type: str = "",
    ) -> int:
        """
        Load Blueprint Epigenome PCHi-C data.
        
        Parameters
        ----------
        filepath : str
            Path to Blueprint PCHi-C file.
        cell_type : str
            Cell type name.
            
        Returns
        -------
        int
            Number of links loaded.
        """
        return self.load_chicago_file(filepath, cell_type, "chicago")
    
    def load_3div(self, filepath: str) -> int:
        """
        Load 3DIV database PCHi-C data.
        
        Parameters
        ----------
        filepath : str
            Path to 3DIV export file.
            
        Returns
        -------
        int
            Number of links loaded.
        """
        return self.load_chicago_file(filepath, file_format="3div")
    
    def _add_link(self, link: PCHiCLink) -> None:
        """Add link to storage and indices."""
        idx = len(self._links)
        self._links.append(link)
        
        # Index by OE (enhancer) chromosome
        oe_chr = link.oe_chr.replace("chr", "")
        if oe_chr not in self._by_oe_region:
            self._by_oe_region[oe_chr] = []
        self._by_oe_region[oe_chr].append(idx)
        
        # Index by bait (promoter) chromosome  
        bait_chr = link.bait_chr.replace("chr", "")
        if bait_chr not in self._by_bait_region:
            self._by_bait_region[bait_chr] = []
        self._by_bait_region[bait_chr].append(idx)
        
        # Index by gene
        gene_key = link.gene_symbol.upper()
        if gene_key not in self._by_gene:
            self._by_gene[gene_key] = []
        self._by_gene[gene_key].append(idx)
        
        # Index by cell type
        if link.cell_type not in self._by_cell_type:
            self._by_cell_type[link.cell_type] = []
        self._by_cell_type[link.cell_type].append(idx)
    
    def get_links_for_region(
        self,
        chrom: str,
        start: int,
        end: int,
        cell_types: Optional[List[str]] = None,
        side: str = "oe",
    ) -> List[PCHiCLink]:
        """
        Get PCHi-C links overlapping a region.
        
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
        side : str
            Which side to match: "oe" (enhancer) or "bait" (promoter).
            
        Returns
        -------
        list
            List of PCHiCLink objects.
        """
        results = []
        
        chrom_norm = chrom.replace("chr", "")
        index = self._by_oe_region if side == "oe" else self._by_bait_region
        
        if chrom_norm not in index:
            return []
        
        for idx in index[chrom_norm]:
            link = self._links[idx]
            
            if not link.overlaps_region(chrom, start, end, side):
                continue
            
            if cell_types and link.cell_type not in cell_types:
                continue
            
            results.append(link)
        
        return results
    
    def get_links_for_gene(
        self,
        gene: str,
        cell_types: Optional[List[str]] = None,
        min_score: Optional[float] = None,
    ) -> List[PCHiCLink]:
        """
        Get PCHi-C links for a specific gene.
        
        Parameters
        ----------
        gene : str
            Gene symbol.
        cell_types : list, optional
            Filter to specific cell types.
        min_score : float, optional
            Minimum CHiCAGO score.
            
        Returns
        -------
        list
            List of PCHiCLink objects.
        """
        gene_key = gene.upper()
        
        if gene_key not in self._by_gene:
            return []
        
        results = []
        min_chicago = min_score or self.chicago_threshold
        
        for idx in self._by_gene[gene_key]:
            link = self._links[idx]
            
            if link.chicago_score < min_chicago:
                continue
            
            if cell_types and link.cell_type not in cell_types:
                continue
            
            results.append(link)
        
        return sorted(results, key=lambda x: x.chicago_score, reverse=True)
    
    def get_gene_for_enhancer(
        self,
        enhancer_chr: str,
        enhancer_start: int,
        enhancer_end: int,
        cell_types: Optional[List[str]] = None,
    ) -> Optional[str]:
        """
        Get best gene target for an enhancer based on PCHi-C.
        
        Parameters
        ----------
        enhancer_chr : str
            Enhancer chromosome.
        enhancer_start : int
            Enhancer start.
        enhancer_end : int
            Enhancer end.
        cell_types : list, optional
            Filter cell types.
            
        Returns
        -------
        str or None
            Best target gene or None.
        """
        links = self.get_links_for_region(
            enhancer_chr, enhancer_start, enhancer_end,
            cell_types=cell_types, side="oe"
        )
        
        if not links:
            return None
        
        best_link = max(links, key=lambda x: x.chicago_score)
        return best_link.gene_symbol
    
    def compute_edge_probability(
        self,
        enhancer_chr: str,
        enhancer_start: int,
        enhancer_end: int,
        gene: str,
        cell_types: Optional[List[str]] = None,
    ) -> Tuple[float, Dict[str, Any]]:
        """
        Compute edge probability for cCRE→gene using PCHi-C evidence.
        
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
            enhancer_chr, enhancer_start, enhancer_end,
            cell_types=cell_types, side="oe"
        )
        
        gene_upper = gene.upper()
        gene_links = [l for l in links if l.gene_symbol.upper() == gene_upper]
        
        evidence = {
            "source": "pchic",
            "n_links": len(gene_links),
            "cell_types": list(set(l.cell_type for l in gene_links)),
        }
        
        if not gene_links:
            return 0.0, evidence
        
        # Use maximum CHiCAGO score
        best_link = max(gene_links, key=lambda x: x.chicago_score)
        chicago_score = best_link.chicago_score
        
        # Convert CHiCAGO score to probability
        # CHiCAGO 5 = ~95% specificity threshold
        # CHiCAGO 10 = very high confidence
        # Sigmoid transformation
        prob = 1 / (1 + np.exp(-0.5 * (chicago_score - 5)))
        
        evidence["chicago_score"] = chicago_score
        evidence["best_cell_type"] = best_link.cell_type
        evidence["distance"] = best_link.distance
        evidence["n_reads"] = best_link.n_reads
        
        return prob, evidence
    
    @property
    def n_links(self) -> int:
        """Number of loaded links."""
        return len(self._links)
    
    @property
    def cell_types(self) -> List[str]:
        """List of cell types with links."""
        return list(self._by_cell_type.keys())
    
    def summary(self) -> Dict[str, Any]:
        """Get summary statistics."""
        return {
            "n_links": self.n_links,
            "n_cell_types": len(self._by_cell_type),
            "n_genes": len(self._by_gene),
            "chicago_threshold": self.chicago_threshold,
            "score_distribution": {
                "min": min(l.chicago_score for l in self._links) if self._links else 0,
                "max": max(l.chicago_score for l in self._links) if self._links else 0,
                "median": float(np.median([l.chicago_score for l in self._links])) if self._links else 0,
            },
        }


def load_pchic_links(
    filepath: str,
    chicago_threshold: float = 5.0,
    cell_type: str = "",
) -> PCHiCLinks:
    """
    Convenience function to load PCHi-C links.
    
    Parameters
    ----------
    filepath : str
        Path to PCHi-C file.
    chicago_threshold : float
        Minimum CHiCAGO score.
    cell_type : str
        Cell type name.
        
    Returns
    -------
    PCHiCLinks
        Loaded PCHi-C links manager.
    """
    pchic = PCHiCLinks(chicago_threshold=chicago_threshold)
    pchic.load_chicago_file(filepath, cell_type=cell_type)
    return pchic
