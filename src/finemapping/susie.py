"""
SuSiE Fine-Mapping Implementation

Wrapper for SuSiE (Sum of Single Effects) fine-mapping using summary statistics.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..utils.logging import get_logger


logger = get_logger("susie")


class SuSiEFinemapper:
    """
    Fine-mapping using SuSiE (Sum of Single Effects).
    
    Supports both:
    - SuSiE-RSS (summary statistics + LD reference)
    - SuSiE with individual-level data (if available)
    """
    
    def __init__(
        self,
        max_causal: int = 10,
        coverage: float = 0.95,
        min_abs_corr: float = 0.5,
        prior_variance: float = 0.2,
        estimate_prior_variance: bool = True,
    ):
        """
        Initialize SuSiE parameters.
        
        Parameters
        ----------
        max_causal : int
            Maximum number of causal variants (L parameter).
        coverage : float
            Target coverage for credible sets.
        min_abs_corr : float
            Minimum absolute correlation for credible set purity.
        prior_variance : float
            Prior variance for effect sizes.
        estimate_prior_variance : bool
            Whether to estimate prior variance from data.
        """
        self.max_causal = max_causal
        self.coverage = coverage
        self.min_abs_corr = min_abs_corr
        self.prior_variance = prior_variance
        self.estimate_prior_variance = estimate_prior_variance
        
        # Check R/SuSiE availability
        self._check_susie()
    
    def _check_susie(self):
        """Check if R and susieR are available."""
        try:
            result = subprocess.run(
                ["Rscript", "-e", "library(susieR); cat('OK')"],
                capture_output=True,
                text=True,
            )
            if "OK" not in result.stdout:
                logger.warning("susieR not found. Install with: install.packages('susieR')")
        except FileNotFoundError:
            logger.warning("R not found. Please install R and susieR.")
    
    def _prepare_input(
        self,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n: int,
        variant_ids: List[str],
    ) -> Tuple[str, str, str]:
        """
        Prepare input files for SuSiE-RSS.
        
        Parameters
        ----------
        z_scores : np.ndarray
            Z-scores for each variant.
        ld_matrix : np.ndarray
            LD correlation matrix.
        n : int
            Sample size.
        variant_ids : list
            Variant identifiers.
            
        Returns
        -------
        tuple
            Paths to (z_file, ld_file, variants_file).
        """
        tmpdir = tempfile.mkdtemp()
        
        # Save Z-scores
        z_file = f"{tmpdir}/z_scores.txt"
        np.savetxt(z_file, z_scores, fmt="%.6f")
        
        # Save LD matrix
        ld_file = f"{tmpdir}/ld_matrix.txt"
        np.savetxt(ld_file, ld_matrix, fmt="%.6f")
        
        # Save variant IDs
        var_file = f"{tmpdir}/variants.txt"
        with open(var_file, 'w') as f:
            for v in variant_ids:
                f.write(f"{v}\n")
        
        return z_file, ld_file, var_file
    
    def _run_susie_rss(
        self,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n: int,
        variant_ids: List[str],
        priors: Optional[np.ndarray] = None,
    ) -> Dict[str, Any]:
        """
        Run SuSiE-RSS via R.
        
        Parameters
        ----------
        z_scores : np.ndarray
            Z-scores.
        ld_matrix : np.ndarray
            LD matrix.
        n : int
            Sample size.
        variant_ids : list
            Variant IDs.
        priors : np.ndarray, optional
            Prior probabilities for each variant.
            
        Returns
        -------
        dict
            SuSiE results including PIPs and credible sets.
        """
        z_file, ld_file, var_file = self._prepare_input(
            z_scores, ld_matrix, n, variant_ids
        )
        
        # Build R script
        prior_arg = ""
        if priors is not None:
            prior_file = f"{Path(z_file).parent}/priors.txt"
            np.savetxt(prior_file, priors, fmt="%.6f")
            prior_arg = f'prior_weights <- scan("{prior_file}", what=numeric())'
            prior_use = ", prior_weights = prior_weights"
        else:
            prior_use = ""
        
        r_script = f'''
library(susieR)

# Load data
z <- scan("{z_file}", what=numeric())
R <- as.matrix(read.table("{ld_file}"))
variants <- scan("{var_file}", what=character())
n <- {n}
L <- {self.max_causal}

{prior_arg}

# Run SuSiE-RSS
fit <- susie_rss(
    z = z,
    R = R,
    n = n,
    L = L,
    coverage = {self.coverage},
    min_abs_corr = {self.min_abs_corr},
    estimate_prior_variance = {"TRUE" if self.estimate_prior_variance else "FALSE"}
    {prior_use}
)

# Extract PIPs
pip <- susie_get_pip(fit)

# Extract credible sets
cs <- susie_get_cs(fit, coverage = {self.coverage})

# Save results
write.table(
    data.frame(variant = variants, pip = pip),
    "{Path(z_file).parent}/pip.txt",
    row.names = FALSE,
    quote = FALSE
)

# Save credible sets
cs_df <- data.frame()
if (length(cs$cs) > 0) {{
    for (i in seq_along(cs$cs)) {{
        cs_df <- rbind(cs_df, data.frame(
            cs_id = i,
            variant_idx = cs$cs[[i]],
            coverage = cs$coverage[i],
            purity = cs$purity[i]
        ))
    }}
}}
write.table(
    cs_df,
    "{Path(z_file).parent}/cs.txt",
    row.names = FALSE,
    quote = FALSE
)

# Save convergence info
cat(fit$converged, file = "{Path(z_file).parent}/converged.txt")
'''
        
        # Write and run R script
        script_file = f"{Path(z_file).parent}/susie.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        result = subprocess.run(
            ["Rscript", script_file],
            capture_output=True,
            text=True,
        )
        
        if result.returncode != 0:
            logger.error(f"SuSiE error: {result.stderr}")
            raise RuntimeError(f"SuSiE failed: {result.stderr}")
        
        # Read results
        pip_df = pd.read_csv(f"{Path(z_file).parent}/pip.txt", sep=r'\s+')
        cs_df = pd.read_csv(f"{Path(z_file).parent}/cs.txt", sep=r'\s+')
        
        with open(f"{Path(z_file).parent}/converged.txt") as f:
            converged = f.read().strip() == "TRUE"
        
        # Build credible sets
        credible_sets = []
        if len(cs_df) > 0:
            for cs_id in cs_df["cs_id"].unique():
                cs_data = cs_df[cs_df["cs_id"] == cs_id]
                indices = cs_data["variant_idx"].values - 1  # R is 1-indexed
                
                credible_sets.append({
                    "cs_id": int(cs_id),
                    "variants": [variant_ids[i] for i in indices],
                    "variant_indices": indices.tolist(),
                    "coverage": cs_data["coverage"].iloc[0],
                    "purity": cs_data["purity"].iloc[0],
                })
        
        return {
            "pip": dict(zip(pip_df["variant"], pip_df["pip"])),
            "credible_sets": credible_sets,
            "converged": converged,
            "n_variants": len(variant_ids),
        }
    
    def finemap_locus(
        self,
        locus: Dict[str, Any],
        ld_matrix: np.ndarray,
        priors: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        """
        Fine-map a single locus.
        
        Parameters
        ----------
        locus : dict
            Locus data including variants.
        ld_matrix : np.ndarray
            LD matrix for the locus.
        priors : dict, optional
            Prior probabilities by variant ID.
            
        Returns
        -------
        dict
            Fine-mapping results.
        """
        variants = locus["variants"]
        
        # Extract data
        variant_ids = [v["rsid"] for v in variants]
        z_scores = np.array([v["beta"] / v["se"] for v in variants])
        n = variants[0].get("n", 100000)  # Default sample size
        
        # Prepare priors
        prior_array = None
        if priors is not None:
            prior_array = np.array([priors.get(v, 1.0) for v in variant_ids])
            prior_array = prior_array / prior_array.sum()  # Normalize
        
        # Run SuSiE
        try:
            results = self._run_susie_rss(
                z_scores=z_scores,
                ld_matrix=ld_matrix,
                n=n,
                variant_ids=variant_ids,
                priors=prior_array,
            )
            
            # Add locus info
            results["locus_id"] = locus["locus_id"]
            results["chr"] = locus["chr"]
            results["start"] = locus["start"]
            results["end"] = locus["end"]
            results["lead_variant"] = locus["lead_variant"]
            
            return results
            
        except Exception as e:
            logger.error(f"Fine-mapping failed for {locus['locus_id']}: {e}")
            return {
                "locus_id": locus["locus_id"],
                "error": str(e),
                "pip": {},
                "credible_sets": [],
            }
    
    def finemap_all_loci(
        self,
        loci: List[Dict],
        ld_reference_path: str,
        priors: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> List[Dict]:
        """
        Fine-map all loci.
        
        Parameters
        ----------
        loci : list
            List of locus dictionaries.
        ld_reference_path : str
            Path to LD reference panel.
        priors : dict, optional
            Prior probabilities by locus_id and variant_id.
            
        Returns
        -------
        list
            Fine-mapping results for all loci.
        """
        results = []
        
        for i, locus in enumerate(loci):
            logger.info(f"Fine-mapping locus {i+1}/{len(loci)}: {locus['locus_id']}")
            
            # Get LD matrix for this locus
            # (In practice, this would query the LD reference)
            ld_matrix = self._get_ld_matrix(locus, ld_reference_path)
            
            # Get priors if available
            locus_priors = None
            if priors and locus["locus_id"] in priors:
                locus_priors = priors[locus["locus_id"]]
            
            # Fine-map
            result = self.finemap_locus(locus, ld_matrix, locus_priors)
            results.append(result)
        
        return results
    
    def _get_ld_matrix(
        self,
        locus: Dict,
        ld_reference_path: str,
    ) -> np.ndarray:
        """
        Get LD matrix for a locus from reference panel.
        
        Parameters
        ----------
        locus : dict
            Locus data.
        ld_reference_path : str
            Path to LD reference.
            
        Returns
        -------
        np.ndarray
            LD correlation matrix.
        """
        # Placeholder - returns identity matrix
        # In practice, would use PLINK or LDstore to compute LD
        n_variants = len(locus["variants"])
        
        logger.warning(f"Using identity LD matrix for {locus['locus_id']} (placeholder)")
        
        return np.eye(n_variants)


def run_finemapping(
    loci: List[Dict],
    ld_reference_path: str,
    **kwargs,
) -> List[Dict]:
    """
    Convenience function to run fine-mapping.
    
    Parameters
    ----------
    loci : list
        List of loci to fine-map.
    ld_reference_path : str
        Path to LD reference panel.
    **kwargs
        Arguments for SuSiEFinemapper.
        
    Returns
    -------
    list
        Fine-mapping results.
    """
    finemapper = SuSiEFinemapper(**kwargs)
    return finemapper.finemap_all_loci(loci, ld_reference_path)
