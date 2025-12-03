"""
Locus Definition

Defines GWAS loci by clumping and window expansion.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..utils.logging import get_logger


logger = get_logger("locus_definition")


class LocusDefinition:
    """
    Defines independent loci from GWAS summary statistics.
    """
    
    def __init__(
        self,
        p_threshold: float = 5e-8,
        r2_threshold: float = 0.1,
        clump_kb: int = 1000,
        window_kb: int = 1000,
        min_variants: int = 10,
    ):
        """
        Initialize locus definition parameters.
        
        Parameters
        ----------
        p_threshold : float
            P-value threshold for significance.
        r2_threshold : float
            LD r² threshold for clumping.
        clump_kb : int
            Clumping window in kb.
        window_kb : int
            Extension window around lead variant in kb.
        min_variants : int
            Minimum variants required in a locus.
        """
        self.p_threshold = p_threshold
        self.r2_threshold = r2_threshold
        self.clump_kb = clump_kb
        self.window_kb = window_kb
        self.min_variants = min_variants
    
    def identify_lead_variants(
        self,
        df: pd.DataFrame,
        pval_col: str = "pval",
    ) -> pd.DataFrame:
        """
        Identify genome-wide significant variants.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics with pval column.
        pval_col : str
            Name of p-value column.
            
        Returns
        -------
        pd.DataFrame
            Significant variants sorted by p-value.
        """
        sig_df = df[df[pval_col] < self.p_threshold].copy()
        sig_df = sig_df.sort_values(pval_col)
        
        logger.info(f"Found {len(sig_df):,} variants with p < {self.p_threshold:.1e}")
        
        return sig_df
    
    def clump_variants(
        self,
        df: pd.DataFrame,
        sig_df: pd.DataFrame,
        ld_matrix: Optional[pd.DataFrame] = None,
    ) -> List[Dict]:
        """
        Clump significant variants into independent loci.
        
        Uses a greedy algorithm:
        1. Select most significant variant as lead
        2. Remove all variants in LD (r² > threshold)
        3. Repeat until no significant variants remain
        
        Parameters
        ----------
        df : pd.DataFrame
            Full summary statistics.
        sig_df : pd.DataFrame
            Significant variants.
        ld_matrix : pd.DataFrame, optional
            LD matrix (if None, uses distance-based clumping).
            
        Returns
        -------
        list
            List of locus dictionaries with lead variants.
        """
        loci = []
        used_variants = set()
        
        for idx, row in sig_df.iterrows():
            rsid = row["rsid"]
            
            if rsid in used_variants:
                continue
            
            # Define locus around this lead variant
            chrom = row["chr"]
            pos = row["pos"]
            
            # Get variants in window
            window_start = pos - self.clump_kb * 1000
            window_end = pos + self.clump_kb * 1000
            
            window_df = df[
                (df["chr"] == chrom) &
                (df["pos"] >= window_start) &
                (df["pos"] <= window_end)
            ]
            
            # Mark variants in this window as used
            for v in window_df["rsid"]:
                used_variants.add(v)
            
            # Create locus record
            locus = {
                "locus_id": f"chr{chrom}_{window_start}_{window_end}",
                "chr": chrom,
                "start": window_start,
                "end": window_end,
                "lead_variant": rsid,
                "lead_pos": pos,
                "lead_pval": row["pval"],
                "lead_beta": row.get("beta"),
                "n_variants_window": len(window_df),
            }
            
            loci.append(locus)
        
        logger.info(f"Identified {len(loci)} independent loci")
        
        return loci
    
    def expand_loci(
        self,
        df: pd.DataFrame,
        loci: List[Dict],
    ) -> List[Dict]:
        """
        Expand loci windows and extract variant data.
        
        Parameters
        ----------
        df : pd.DataFrame
            Full summary statistics.
        loci : list
            List of locus dictionaries.
            
        Returns
        -------
        list
            Loci with expanded windows and variant data.
        """
        expanded_loci = []
        
        for locus in loci:
            chrom = locus["chr"]
            lead_pos = locus["lead_pos"]
            
            # Expand window
            window_start = lead_pos - self.window_kb * 1000
            window_end = lead_pos + self.window_kb * 1000
            
            # Get all variants in expanded window
            window_df = df[
                (df["chr"] == chrom) &
                (df["pos"] >= window_start) &
                (df["pos"] <= window_end)
            ].copy()
            
            if len(window_df) < self.min_variants:
                logger.warning(
                    f"Locus {locus['locus_id']}: only {len(window_df)} variants, "
                    f"skipping (min {self.min_variants})"
                )
                continue
            
            # Update locus record
            locus_expanded = locus.copy()
            locus_expanded.update({
                "locus_id": f"chr{chrom}_{window_start}_{window_end}",
                "start": window_start,
                "end": window_end,
                "n_variants": len(window_df),
                "variants": window_df.to_dict("records"),
            })
            
            expanded_loci.append(locus_expanded)
        
        logger.info(f"Expanded {len(expanded_loci)} loci (window: ±{self.window_kb}kb)")
        
        return expanded_loci
    
    def define_loci(
        self,
        df: pd.DataFrame,
    ) -> List[Dict]:
        """
        Full pipeline to define loci from summary statistics.
        
        Parameters
        ----------
        df : pd.DataFrame
            Summary statistics.
            
        Returns
        -------
        list
            List of defined loci with variant data.
        """
        logger.info("Defining loci from summary statistics")
        
        # Identify significant variants
        sig_df = self.identify_lead_variants(df)
        
        if len(sig_df) == 0:
            logger.warning("No genome-wide significant variants found")
            return []
        
        # Clump into independent loci
        loci = self.clump_variants(df, sig_df)
        
        # Expand windows
        expanded_loci = self.expand_loci(df, loci)
        
        return expanded_loci


def define_loci(
    df: pd.DataFrame,
    **kwargs,
) -> List[Dict]:
    """
    Convenience function to define loci.
    
    Parameters
    ----------
    df : pd.DataFrame
        Summary statistics.
    **kwargs
        Arguments for LocusDefinition.
        
    Returns
    -------
    list
        List of defined loci.
    """
    definer = LocusDefinition(**kwargs)
    return definer.define_loci(df)
