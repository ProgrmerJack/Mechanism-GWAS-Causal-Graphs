"""
GWAS Summary Statistics Harmonizer

Standardizes column names, allele coding, and coordinates across different
GWAS summary statistics formats.
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..utils.config import get_config
from ..utils.io import read_sumstats, write_sumstats
from ..utils.genomics import is_ambiguous_snp, harmonize_alleles
from ..utils.logging import get_logger, setup_logger


logger = setup_logger("harmonization")


class GWASHarmonizer:
    """
    Harmonizes GWAS summary statistics to a standard format.
    
    The standard format includes:
    - chr: Chromosome (1-22, X, Y, MT)
    - pos: Position (1-based, GRCh38)
    - rsid: RS identifier
    - effect_allele: Effect allele (the allele that beta refers to)
    - other_allele: Non-effect allele
    - beta: Effect size (log-OR for binary traits)
    - se: Standard error of beta
    - pval: P-value
    - eaf: Effect allele frequency
    - n: Sample size
    - n_cases: Number of cases (for binary traits)
    - n_controls: Number of controls (for binary traits)
    """
    
    # Column name mappings for common formats
    COLUMN_MAPPINGS = {
        # Chromosome
        "chr": ["chr", "chrom", "chromosome", "#chr", "#chrom", "CHR", "CHROM"],
        # Position
        "pos": ["pos", "position", "bp", "BP", "POS", "base_pair_location"],
        # RS ID
        "rsid": ["rsid", "snp", "variant_id", "SNP", "RSID", "rs_id", "MarkerName", "ID"],
        # Effect allele
        "effect_allele": ["effect_allele", "a1", "A1", "alt", "ALT", "allele1", "ea", "EA", 
                         "effect_allele_frequency", "Allele1", "tested_allele"],
        # Other allele
        "other_allele": ["other_allele", "a2", "A2", "ref", "REF", "allele2", "oa", "OA",
                        "Allele2", "non_effect_allele", "nea", "NEA"],
        # Beta/effect size
        "beta": ["beta", "effect", "b", "BETA", "Effect", "effect_size", "log_odds"],
        # Standard error
        "se": ["se", "stderr", "standard_error", "SE", "StdErr"],
        # P-value
        "pval": ["pval", "p", "pvalue", "p_value", "P", "Pvalue", "p-value"],
        # Effect allele frequency
        "eaf": ["eaf", "freq", "af", "maf", "MAF", "effect_allele_frequency", "EAF", 
               "frequency", "Freq1"],
        # Sample size
        "n": ["n", "N", "n_total", "sample_size", "TotalSampleSize"],
        # Number of cases
        "n_cases": ["n_cases", "Ncases", "N_cases", "cases", "n_case"],
        # Number of controls  
        "n_controls": ["n_controls", "Ncontrols", "N_controls", "controls", "n_control"],
        # INFO score
        "info": ["info", "INFO", "imputation_quality", "info_score"],
    }
    
    def __init__(
        self,
        target_build: str = "GRCh38",
        remove_ambiguous: bool = True,
        max_ambiguous_maf: float = 0.40,
        min_maf: float = 0.001,
    ):
        """
        Initialize the harmonizer.
        
        Parameters
        ----------
        target_build : str
            Target genome build.
        remove_ambiguous : bool
            Whether to remove ambiguous (A/T, C/G) SNPs.
        max_ambiguous_maf : float
            MAF threshold above which ambiguous SNPs are removed.
        min_maf : float
            Minimum MAF filter.
        """
        self.target_build = target_build
        self.remove_ambiguous = remove_ambiguous
        self.max_ambiguous_maf = max_ambiguous_maf
        self.min_maf = min_maf
        
    def _map_columns(self, df: pd.DataFrame) -> Dict[str, str]:
        """
        Map input column names to standard names.
        
        Parameters
        ----------
        df : pd.DataFrame
            Input dataframe.
            
        Returns
        -------
        dict
            Mapping from standard name to actual column name.
        """
        column_map = {}
        input_cols_lower = {c.lower(): c for c in df.columns}
        
        for standard_name, variants in self.COLUMN_MAPPINGS.items():
            for variant in variants:
                if variant.lower() in input_cols_lower:
                    column_map[standard_name] = input_cols_lower[variant.lower()]
                    break
        
        return column_map
    
    def _standardize_chromosome(self, chrom: Union[str, int]) -> str:
        """
        Standardize chromosome notation.
        
        Parameters
        ----------
        chrom : str or int
            Input chromosome.
            
        Returns
        -------
        str
            Standardized chromosome (e.g., "1", "X", "MT").
        """
        chrom = str(chrom).strip()
        
        # Remove 'chr' prefix
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        
        # Handle special chromosomes
        chrom_upper = chrom.upper()
        if chrom_upper in ["23", "X"]:
            return "X"
        elif chrom_upper in ["24", "Y"]:
            return "Y"
        elif chrom_upper in ["25", "MT", "M", "MITO"]:
            return "MT"
        
        return chrom
    
    def _standardize_allele(self, allele: str) -> str:
        """
        Standardize allele notation.
        
        Parameters
        ----------
        allele : str
            Input allele.
            
        Returns
        -------
        str
            Uppercase allele.
        """
        if pd.isna(allele):
            return ""
        return str(allele).strip().upper()
    
    def harmonize(
        self,
        df: pd.DataFrame,
        trait_name: Optional[str] = None,
        sample_size: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Harmonize GWAS summary statistics.
        
        Parameters
        ----------
        df : pd.DataFrame
            Input summary statistics.
        trait_name : str, optional
            Name of the trait (for logging).
        sample_size : int, optional
            Override sample size if not in data.
            
        Returns
        -------
        pd.DataFrame
            Harmonized summary statistics.
        """
        logger.info(f"Harmonizing summary statistics{f' for {trait_name}' if trait_name else ''}")
        logger.info(f"Input: {len(df):,} variants")
        
        # Map columns
        col_map = self._map_columns(df)
        logger.info(f"Mapped columns: {col_map}")
        
        # Check required columns
        required = ["chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval"]
        missing = [c for c in required if c not in col_map]
        
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        # Create standardized dataframe
        harmonized = pd.DataFrame()
        
        # Chromosome
        harmonized["chr"] = df[col_map["chr"]].apply(self._standardize_chromosome)
        
        # Position
        harmonized["pos"] = pd.to_numeric(df[col_map["pos"]], errors="coerce").astype("Int64")
        
        # RS ID
        if "rsid" in col_map:
            harmonized["rsid"] = df[col_map["rsid"]].astype(str)
        else:
            # Generate variant ID
            harmonized["rsid"] = harmonized.apply(
                lambda r: f"chr{r['chr']}:{r['pos']}", axis=1
            )
        
        # Alleles
        harmonized["effect_allele"] = df[col_map["effect_allele"]].apply(self._standardize_allele)
        harmonized["other_allele"] = df[col_map["other_allele"]].apply(self._standardize_allele)
        
        # Effect size
        harmonized["beta"] = pd.to_numeric(df[col_map["beta"]], errors="coerce")
        harmonized["se"] = pd.to_numeric(df[col_map["se"]], errors="coerce")
        
        # P-value
        harmonized["pval"] = pd.to_numeric(df[col_map["pval"]], errors="coerce")
        
        # Optional columns
        if "eaf" in col_map:
            harmonized["eaf"] = pd.to_numeric(df[col_map["eaf"]], errors="coerce")
        else:
            harmonized["eaf"] = np.nan
        
        if "n" in col_map:
            harmonized["n"] = pd.to_numeric(df[col_map["n"]], errors="coerce").astype("Int64")
        elif sample_size:
            harmonized["n"] = sample_size
        else:
            harmonized["n"] = pd.NA
        
        if "n_cases" in col_map:
            harmonized["n_cases"] = pd.to_numeric(df[col_map["n_cases"]], errors="coerce").astype("Int64")
        
        if "n_controls" in col_map:
            harmonized["n_controls"] = pd.to_numeric(df[col_map["n_controls"]], errors="coerce").astype("Int64")
        
        if "info" in col_map:
            harmonized["info"] = pd.to_numeric(df[col_map["info"]], errors="coerce")
        
        # Calculate Z-score
        harmonized["z"] = harmonized["beta"] / harmonized["se"]
        
        # Log initial count
        n_initial = len(harmonized)
        
        # Filter by chromosome
        valid_chroms = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
        harmonized = harmonized[harmonized["chr"].isin(valid_chroms)]
        logger.info(f"After chromosome filter: {len(harmonized):,} ({n_initial - len(harmonized):,} removed)")
        
        # Filter missing values
        harmonized = harmonized.dropna(subset=["chr", "pos", "beta", "se", "pval"])
        logger.info(f"After missing value filter: {len(harmonized):,}")
        
        # Filter ambiguous SNPs
        if self.remove_ambiguous:
            is_ambig = harmonized.apply(
                lambda r: is_ambiguous_snp(r["effect_allele"], r["other_allele"]),
                axis=1
            )
            
            # Keep ambiguous SNPs with low MAF
            if "eaf" in harmonized.columns:
                maf = harmonized["eaf"].apply(lambda x: min(x, 1-x) if pd.notna(x) else 0)
                keep_ambig = ~is_ambig | (maf < self.max_ambiguous_maf)
                n_ambig = (~keep_ambig).sum()
                harmonized = harmonized[keep_ambig]
                logger.info(f"Removed {n_ambig:,} ambiguous SNPs with MAF > {self.max_ambiguous_maf}")
            else:
                n_ambig = is_ambig.sum()
                harmonized = harmonized[~is_ambig]
                logger.info(f"Removed {n_ambig:,} ambiguous SNPs (no MAF for filtering)")
        
        # Filter by MAF
        if self.min_maf > 0 and "eaf" in harmonized.columns:
            maf = harmonized["eaf"].apply(lambda x: min(x, 1-x) if pd.notna(x) else 0.5)
            n_before = len(harmonized)
            harmonized = harmonized[maf >= self.min_maf]
            logger.info(f"After MAF filter (>{self.min_maf}): {len(harmonized):,} ({n_before - len(harmonized):,} removed)")
        
        # Sort by chromosome and position
        harmonized = harmonized.sort_values(["chr", "pos"]).reset_index(drop=True)
        
        # Add trait name
        if trait_name:
            harmonized["trait"] = trait_name
        
        logger.info(f"Final: {len(harmonized):,} variants")
        
        return harmonized
    
    def harmonize_file(
        self,
        input_path: str | Path,
        output_path: str | Path,
        trait_name: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Harmonize a summary statistics file.
        
        Parameters
        ----------
        input_path : str or Path
            Input file path.
        output_path : str or Path
            Output file path.
        trait_name : str, optional
            Name of the trait.
        **kwargs
            Additional arguments for read_sumstats.
            
        Returns
        -------
        pd.DataFrame
            Harmonized summary statistics.
        """
        # Read input
        df = read_sumstats(input_path, **kwargs)
        
        # Harmonize
        harmonized = self.harmonize(df, trait_name=trait_name)
        
        # Write output
        write_sumstats(harmonized, output_path)
        
        return harmonized


def harmonize_gwas(
    input_file: str | Path,
    output_file: str | Path,
    trait_name: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Convenience function for harmonizing GWAS summary statistics.
    
    Parameters
    ----------
    input_file : str or Path
        Path to input file.
    output_file : str or Path
        Path to output file.
    trait_name : str, optional
        Name of the trait.
    **kwargs
        Additional arguments for GWASHarmonizer.
        
    Returns
    -------
    pd.DataFrame
        Harmonized summary statistics.
    """
    harmonizer = GWASHarmonizer(**kwargs)
    return harmonizer.harmonize_file(input_file, output_file, trait_name=trait_name)
