"""
ColocQuake: Multi-Tissue Colocalization Pipeline

Orchestrates colocalization analysis across multiple tissues,
genes, and QTL types (eQTL, sQTL).
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..utils.config import get_config
from ..utils.logging import get_logger
from .coloc import ColocAnalysis
from .eqtl import EQTLLoader, MultiTissueEQTL


logger = get_logger("colocquake")


class ColocQuake:
    """
    Large-scale colocalization pipeline.
    
    For each GWAS locus, tests colocalization with:
    - eQTLs across relevant tissues
    - sQTLs across relevant tissues
    - Other molecular QTLs as available
    
    Named after the earth-shaking scale of analysis required
    for comprehensive multi-trait colocalization.
    """
    
    def __init__(
        self,
        eqtl_dir: str,
        trait_name: str,
        h4_threshold: float = 0.75,
        h3_h4_threshold: float = 0.9,
    ):
        """
        Initialize ColocQuake.
        
        Parameters
        ----------
        eqtl_dir : str
            Path to eQTL data directory.
        trait_name : str
            Name of the GWAS trait.
        h4_threshold : float
            PP.H4 threshold for significant colocalization.
        h3_h4_threshold : float
            Combined PP.H3+PP.H4 threshold.
        """
        self.eqtl_dir = Path(eqtl_dir)
        self.trait_name = trait_name
        self.h4_threshold = h4_threshold
        self.h3_h4_threshold = h3_h4_threshold
        
        # Load trait configuration
        traits_config = get_config("traits")
        self.trait_config = traits_config.get("trait_tissue_priors", {}).get(
            trait_name, {}
        )
        
        # Get relevant tissues
        self.tissues = self._get_relevant_tissues()
        
        # Initialize analyzers
        self.coloc = ColocAnalysis()
        self.multi_eqtl = MultiTissueEQTL(eqtl_dir, self.tissues)
    
    def _get_relevant_tissues(self) -> List[str]:
        """Get tissues relevant to the trait."""
        tissues = []
        
        # Primary tissues
        for t in self.trait_config.get("primary_tissues", []):
            tissues.append(t["tissue"])
        
        # Secondary tissues
        for t in self.trait_config.get("secondary_tissues", []):
            tissues.append(t["tissue"])
        
        return tissues
    
    def run_locus(
        self,
        locus: Dict[str, Any],
        gwas_df: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Run colocalization for a single locus.
        
        Parameters
        ----------
        locus : dict
            Locus definition with chr, start, end.
        gwas_df : pd.DataFrame
            GWAS summary statistics for the locus.
            
        Returns
        -------
        dict
            Colocalization results for all tissue-gene pairs.
        """
        logger.info(f"Running colocalization for {locus['locus_id']}")
        
        chrom = locus["chr"]
        start = locus["start"]
        end = locus["end"]
        
        # Get eQTLs across tissues
        tissue_eqtls = self.multi_eqtl.query_locus_all_tissues(chrom, start, end)
        
        results = {
            "locus_id": locus["locus_id"],
            "chr": chrom,
            "start": start,
            "end": end,
            "colocalizations": [],
            "summary": {
                "n_tissues_tested": 0,
                "n_genes_tested": 0,
                "n_significant_h4": 0,
            },
        }
        
        genes_tested = set()
        
        for tissue, eqtl_df in tissue_eqtls.items():
            logger.debug(f"Testing {tissue}: {len(eqtl_df)} eQTL associations")
            
            # Get unique genes
            if "gene_ensembl" in eqtl_df.columns:
                genes = eqtl_df["gene_ensembl"].unique()
            else:
                genes = eqtl_df["gene_id"].str.split(".").str[0].unique()
            
            for gene in genes:
                genes_tested.add(gene)
                
                # Get eQTLs for this gene
                gene_eqtl = eqtl_df[
                    eqtl_df["gene_id"].str.startswith(gene) |
                    (eqtl_df.get("gene_ensembl") == gene)
                ]
                
                if len(gene_eqtl) < 10:
                    continue
                
                # Run COLOC
                coloc_result = self.coloc.run_coloc(
                    gwas_df=gwas_df,
                    qtl_df=gene_eqtl,
                    gwas_type="cc",
                    qtl_type="quant",
                )
                
                # Store result
                result = {
                    "tissue": tissue,
                    "gene": gene,
                    "pp_h4": coloc_result["pp_h4"],
                    "pp_h3": coloc_result["pp_h3"],
                    "pp_h0": coloc_result["pp_h0"],
                    "pp_h1": coloc_result["pp_h1"],
                    "pp_h2": coloc_result["pp_h2"],
                    "n_snps": coloc_result["n_snps"],
                    "colocalized": coloc_result["pp_h4"] >= self.h4_threshold,
                }
                
                results["colocalizations"].append(result)
                
                if result["colocalized"]:
                    results["summary"]["n_significant_h4"] += 1
        
        results["summary"]["n_tissues_tested"] = len(tissue_eqtls)
        results["summary"]["n_genes_tested"] = len(genes_tested)
        
        logger.info(
            f"Tested {len(genes_tested)} genes across {len(tissue_eqtls)} tissues, "
            f"{results['summary']['n_significant_h4']} significant"
        )
        
        return results
    
    def run_all_loci(
        self,
        loci: List[Dict],
        gwas_df: pd.DataFrame,
    ) -> List[Dict]:
        """
        Run colocalization for all loci.
        
        Parameters
        ----------
        loci : list
            List of locus dictionaries.
        gwas_df : pd.DataFrame
            Full GWAS summary statistics.
            
        Returns
        -------
        list
            Colocalization results for all loci.
        """
        all_results = []
        
        for i, locus in enumerate(loci):
            logger.info(f"Processing locus {i+1}/{len(loci)}")
            
            # Extract GWAS data for this locus
            locus_gwas = gwas_df[
                (gwas_df["chr"].astype(str) == str(locus["chr"])) &
                (gwas_df["pos"] >= locus["start"]) &
                (gwas_df["pos"] <= locus["end"])
            ]
            
            if len(locus_gwas) < 10:
                logger.warning(f"Insufficient GWAS variants for {locus['locus_id']}")
                continue
            
            results = self.run_locus(locus, locus_gwas)
            all_results.append(results)
        
        return all_results
    
    def prioritize_genes(
        self,
        coloc_results: List[Dict],
        min_h4: float = 0.5,
    ) -> pd.DataFrame:
        """
        Prioritize genes based on colocalization evidence.
        
        Parameters
        ----------
        coloc_results : list
            Colocalization results from run_all_loci.
        min_h4 : float
            Minimum PP.H4 to include.
            
        Returns
        -------
        pd.DataFrame
            Prioritized genes with evidence.
        """
        records = []
        
        for locus_result in coloc_results:
            locus_id = locus_result["locus_id"]
            
            for coloc in locus_result["colocalizations"]:
                if coloc["pp_h4"] >= min_h4:
                    records.append({
                        "locus_id": locus_id,
                        "gene": coloc["gene"],
                        "tissue": coloc["tissue"],
                        "pp_h4": coloc["pp_h4"],
                        "pp_h3_h4": coloc["pp_h3"] + coloc["pp_h4"],
                        "n_snps": coloc["n_snps"],
                    })
        
        if not records:
            return pd.DataFrame()
        
        df = pd.DataFrame(records)
        
        # Aggregate by gene
        gene_summary = df.groupby("gene").agg({
            "pp_h4": ["max", "mean", "count"],
            "tissue": lambda x: list(x.unique()),
            "locus_id": lambda x: list(x.unique()),
        }).reset_index()
        
        gene_summary.columns = [
            "gene", "max_h4", "mean_h4", "n_coloc",
            "tissues", "loci",
        ]
        
        gene_summary = gene_summary.sort_values("max_h4", ascending=False)
        
        return gene_summary
    
    def get_tissue_activity(
        self,
        coloc_results: List[Dict],
    ) -> pd.DataFrame:
        """
        Summarize tissue activity patterns.
        
        Parameters
        ----------
        coloc_results : list
            Colocalization results.
            
        Returns
        -------
        pd.DataFrame
            Tissue activity summary.
        """
        tissue_stats = {}
        
        for locus_result in coloc_results:
            for coloc in locus_result["colocalizations"]:
                tissue = coloc["tissue"]
                
                if tissue not in tissue_stats:
                    tissue_stats[tissue] = {
                        "n_tested": 0,
                        "n_colocalized": 0,
                        "genes_colocalized": set(),
                    }
                
                tissue_stats[tissue]["n_tested"] += 1
                
                if coloc["colocalized"]:
                    tissue_stats[tissue]["n_colocalized"] += 1
                    tissue_stats[tissue]["genes_colocalized"].add(coloc["gene"])
        
        records = []
        for tissue, stats in tissue_stats.items():
            records.append({
                "tissue": tissue,
                "n_tested": stats["n_tested"],
                "n_colocalized": stats["n_colocalized"],
                "coloc_rate": stats["n_colocalized"] / max(stats["n_tested"], 1),
                "n_unique_genes": len(stats["genes_colocalized"]),
            })
        
        df = pd.DataFrame(records)
        df = df.sort_values("coloc_rate", ascending=False)
        
        return df
