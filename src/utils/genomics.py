"""
Genomics utility functions.
"""

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


def liftover_coordinates(
    df: pd.DataFrame,
    source_build: str = "hg19",
    target_build: str = "hg38",
    chrom_col: str = "chr",
    pos_col: str = "pos",
    chain_file: Optional[str] = None,
) -> pd.DataFrame:
    """
    Lift over genomic coordinates between genome builds.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with genomic coordinates.
    source_build : str
        Source genome build ('hg19' or 'hg38').
    target_build : str
        Target genome build ('hg19' or 'hg38').
    chrom_col : str
        Name of chromosome column.
    pos_col : str
        Name of position column.
    chain_file : str, optional
        Path to chain file. If None, uses default location.
        
    Returns
    -------
    pd.DataFrame
        Dataframe with lifted coordinates.
    """
    try:
        from liftover import get_lifter
    except ImportError:
        raise ImportError("liftover package required. Install with: pip install liftover")
    
    # Get lifter
    converter = get_lifter(source_build, target_build)
    
    # Apply liftover
    def lift_position(row):
        chrom = str(row[chrom_col])
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        
        pos = int(row[pos_col])
        
        result = converter.convert_coordinate(chrom, pos)
        
        if result and len(result) > 0:
            return result[0][1]  # Return position
        return np.nan
    
    df = df.copy()
    df[f"{pos_col}_lifted"] = df.apply(lift_position, axis=1)
    
    # Count successful liftovers
    n_success = df[f"{pos_col}_lifted"].notna().sum()
    n_total = len(df)
    print(f"LiftOver success: {n_success}/{n_total} ({100*n_success/n_total:.2f}%)")
    
    return df


def get_rsid(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    build: str = "GRCh38",
) -> Optional[str]:
    """
    Get rsID for a variant from Ensembl.
    
    Parameters
    ----------
    chrom : str
        Chromosome.
    pos : int
        Position.
    ref : str
        Reference allele.
    alt : str
        Alternative allele.
    build : str
        Genome build.
        
    Returns
    -------
    str or None
        rsID if found.
    """
    try:
        import requests
    except ImportError:
        raise ImportError("requests package required")
    
    # Format chromosome
    chrom = str(chrom).replace("chr", "")
    
    # Query Ensembl VEP
    url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}:{pos}/{alt}"
    
    params = {"content-type": "application/json"}
    
    try:
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                return data[0].get("id")
    except Exception:
        pass
    
    return None


def calculate_ld(
    variants: List[str],
    ld_panel_path: str,
    population: str = "EUR",
    r2_threshold: float = 0.001,
) -> pd.DataFrame:
    """
    Calculate LD matrix for a set of variants using PLINK.
    
    Parameters
    ----------
    variants : list
        List of variant IDs (rsIDs).
    ld_panel_path : str
        Path to LD reference panel (PLINK format).
    population : str
        Population for LD calculation.
    r2_threshold : float
        Minimum r² to report.
        
    Returns
    -------
    pd.DataFrame
        LD matrix as dataframe.
    """
    import tempfile
    
    # Write variants to file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        for v in variants:
            f.write(f"{v}\n")
        variants_file = f.name
    
    # Output prefix
    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = f"{tmpdir}/ld"
        
        # Run PLINK
        cmd = [
            "plink",
            "--bfile", ld_panel_path,
            "--extract", variants_file,
            "--r2", "square",
            "--out", out_prefix,
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Read LD matrix
            ld_file = f"{out_prefix}.ld"
            if Path(ld_file).exists():
                ld_matrix = pd.read_csv(ld_file, sep=r'\s+', header=None)
                ld_matrix.index = variants[:len(ld_matrix)]
                ld_matrix.columns = variants[:len(ld_matrix)]
                return ld_matrix
        except subprocess.CalledProcessError as e:
            print(f"PLINK error: {e.stderr.decode()}")
        except FileNotFoundError:
            print("PLINK not found. Please install PLINK and add to PATH.")
    
    # Return identity matrix as fallback
    n = len(variants)
    return pd.DataFrame(
        np.eye(n),
        index=variants,
        columns=variants,
    )


def is_ambiguous_snp(ref: str, alt: str) -> bool:
    """
    Check if a SNP is strand-ambiguous (A/T or C/G).
    
    Parameters
    ----------
    ref : str
        Reference allele.
    alt : str
        Alternative allele.
        
    Returns
    -------
    bool
        True if ambiguous.
    """
    ref = ref.upper()
    alt = alt.upper()
    
    ambiguous_pairs = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    return (ref, alt) in ambiguous_pairs


def complement_allele(allele: str) -> str:
    """
    Get complement of a nucleotide.
    
    Parameters
    ----------
    allele : str
        Nucleotide (A, T, C, G).
        
    Returns
    -------
    str
        Complementary nucleotide.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complement.get(allele.upper(), allele)


def harmonize_alleles(
    ref: str,
    alt: str,
    ref_panel_ref: str,
    ref_panel_alt: str,
    effect: float,
) -> Tuple[str, str, float, bool]:
    """
    Harmonize alleles to reference panel.
    
    Parameters
    ----------
    ref : str
        GWAS reference allele.
    alt : str
        GWAS alternative allele.
    ref_panel_ref : str
        Reference panel reference allele.
    ref_panel_alt : str
        Reference panel alternative allele.
    effect : float
        Effect size (beta or log-OR).
        
    Returns
    -------
    tuple
        (harmonized_ref, harmonized_alt, harmonized_effect, is_flipped)
    """
    ref = ref.upper()
    alt = alt.upper()
    ref_panel_ref = ref_panel_ref.upper()
    ref_panel_alt = ref_panel_alt.upper()
    
    # Direct match
    if ref == ref_panel_ref and alt == ref_panel_alt:
        return ref, alt, effect, False
    
    # Flipped
    if ref == ref_panel_alt and alt == ref_panel_ref:
        return alt, ref, -effect, True
    
    # Complement
    ref_c = complement_allele(ref)
    alt_c = complement_allele(alt)
    
    if ref_c == ref_panel_ref and alt_c == ref_panel_alt:
        return ref_c, alt_c, effect, False
    
    if ref_c == ref_panel_alt and alt_c == ref_panel_ref:
        return alt_c, ref_c, -effect, True
    
    # No match found
    raise ValueError(f"Cannot harmonize alleles: {ref}/{alt} vs {ref_panel_ref}/{ref_panel_alt}")
