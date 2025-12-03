"""
LiftOver Pipeline for Genomic Coordinates

Handles coordinate conversion between genome builds (e.g., hg19 → hg38).
"""

from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from ..utils.logging import get_logger


logger = get_logger("liftover")


class LiftOverPipeline:
    """
    Pipeline for lifting over genomic coordinates between genome builds.
    """
    
    def __init__(
        self,
        source_build: str = "hg19",
        target_build: str = "hg38",
    ):
        """
        Initialize the liftover pipeline.
        
        Parameters
        ----------
        source_build : str
            Source genome build ('hg19' or 'hg38').
        target_build : str
            Target genome build ('hg19' or 'hg38').
        """
        self.source_build = source_build
        self.target_build = target_build
        self._lifter = None
        
    @property
    def lifter(self):
        """Lazy load the lifter."""
        if self._lifter is None:
            try:
                from liftover import get_lifter
                self._lifter = get_lifter(self.source_build, self.target_build)
            except ImportError:
                raise ImportError(
                    "liftover package required. Install with: pip install liftover"
                )
        return self._lifter
    
    def lift_position(
        self,
        chrom: str,
        pos: int,
    ) -> Optional[Tuple[str, int]]:
        """
        Lift a single position to the target build.
        
        Parameters
        ----------
        chrom : str
            Chromosome (e.g., "chr1" or "1").
        pos : int
            Position (1-based).
            
        Returns
        -------
        tuple or None
            (chromosome, position) in target build, or None if failed.
        """
        # Ensure chr prefix
        if not str(chrom).startswith("chr"):
            chrom = f"chr{chrom}"
        
        try:
            result = self.lifter.convert_coordinate(chrom, int(pos))
            if result and len(result) > 0:
                new_chrom = result[0][0]
                new_pos = result[0][1]
                # Remove chr prefix for consistency
                new_chrom = new_chrom.replace("chr", "")
                return (new_chrom, new_pos)
        except Exception:
            pass
        
        return None
    
    def lift_dataframe(
        self,
        df: pd.DataFrame,
        chrom_col: str = "chr",
        pos_col: str = "pos",
        drop_unmapped: bool = True,
    ) -> pd.DataFrame:
        """
        Lift over coordinates in a dataframe.
        
        Parameters
        ----------
        df : pd.DataFrame
            Dataframe with genomic coordinates.
        chrom_col : str
            Name of chromosome column.
        pos_col : str
            Name of position column.
        drop_unmapped : bool
            Whether to drop variants that fail liftover.
            
        Returns
        -------
        pd.DataFrame
            Dataframe with lifted coordinates.
        """
        logger.info(f"Lifting {len(df):,} variants from {self.source_build} to {self.target_build}")
        
        df = df.copy()
        
        # Store original positions
        df[f"{pos_col}_original"] = df[pos_col]
        df[f"{chrom_col}_original"] = df[chrom_col]
        
        # Apply liftover
        n_success = 0
        n_failed = 0
        
        for idx in df.index:
            chrom = df.loc[idx, chrom_col]
            pos = df.loc[idx, pos_col]
            
            result = self.lift_position(chrom, pos)
            
            if result is not None:
                df.loc[idx, chrom_col] = result[0]
                df.loc[idx, pos_col] = result[1]
                n_success += 1
            else:
                df.loc[idx, pos_col] = np.nan
                n_failed += 1
        
        logger.info(f"LiftOver: {n_success:,} success, {n_failed:,} failed "
                   f"({100*n_success/(n_success+n_failed):.2f}%)")
        
        if drop_unmapped:
            df = df.dropna(subset=[pos_col])
            df[pos_col] = df[pos_col].astype("Int64")
            logger.info(f"After dropping unmapped: {len(df):,} variants")
        
        return df


def liftover(
    df: pd.DataFrame,
    source_build: str = "hg19",
    target_build: str = "hg38",
    **kwargs,
) -> pd.DataFrame:
    """
    Convenience function for lifting over coordinates.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with coordinates.
    source_build : str
        Source genome build.
    target_build : str
        Target genome build.
    **kwargs
        Additional arguments for lift_dataframe.
        
    Returns
    -------
    pd.DataFrame
        Dataframe with lifted coordinates.
    """
    pipeline = LiftOverPipeline(source_build, target_build)
    return pipeline.lift_dataframe(df, **kwargs)
