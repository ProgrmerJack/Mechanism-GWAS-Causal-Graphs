"""
COLOC Bayesian Colocalization

Implements COLOC methods for testing colocalization between
GWAS and molecular QTL signals.

Supports:
- coloc.abf: Original single-causal-variant method (Giambartolomei et al. 2014)
- coloc.susie: Multi-signal colocalization using SuSiE credible sets (Wallace 2021)

Reference:
- Wallace C (2021). A more accurate method for colocalisation analysis 
  allowing for multiple causal variants. PLOS Genetics.
  https://doi.org/10.1371/journal.pgen.1009440
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..utils.logging import get_logger


logger = get_logger("coloc")


# Colocalization interpretation thresholds
COLOC_THRESHOLDS = {
    "strong": 0.8,       # Strong evidence of colocalization
    "moderate": 0.5,     # Moderate evidence
    "weak": 0.25,        # Weak evidence
}


class ColocAnalysis:
    """
    Bayesian colocalization analysis using COLOC.
    
    Tests five hypotheses:
    - H0: No association in either trait
    - H1: Association with trait 1 only
    - H2: Association with trait 2 only
    - H3: Both traits associated, different causal variants
    - H4: Both traits associated, shared causal variant (colocalization)
    """
    
    def __init__(
        self,
        p1: float = 1e-4,
        p2: float = 1e-4,
        p12: float = 1e-5,
    ):
        """
        Initialize COLOC parameters.
        
        Parameters
        ----------
        p1 : float
            Prior probability a variant is associated with trait 1.
        p2 : float
            Prior probability a variant is associated with trait 2.
        p12 : float
            Prior probability a variant is associated with both traits.
        """
        self.p1 = p1
        self.p2 = p2
        self.p12 = p12
        
        # Check R/coloc availability
        self._check_coloc()
    
    def _check_coloc(self):
        """Check if R and coloc are available."""
        try:
            result = subprocess.run(
                ["Rscript", "-e", "library(coloc); cat('OK')"],
                capture_output=True,
                text=True,
            )
            if "OK" not in result.stdout:
                logger.warning("coloc not found. Install with: install.packages('coloc')")
        except FileNotFoundError:
            logger.warning("R not found. Please install R and coloc.")
    
    def _prepare_dataset(
        self,
        sumstats: pd.DataFrame,
        dataset_type: str = "quant",
        sample_size: int = 100000,
        sdY: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Prepare a dataset for COLOC.
        
        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics with columns: rsid, beta, se, pval, maf.
        dataset_type : str
            Type: "quant" for quantitative, "cc" for case-control.
        sample_size : int
            Sample size.
        sdY : float, optional
            Standard deviation of Y (for quantitative traits).
            
        Returns
        -------
        dict
            Dataset dictionary for COLOC.
        """
        required_cols = ["rsid", "beta", "se", "pval"]
        for col in required_cols:
            if col not in sumstats.columns:
                raise ValueError(f"Missing required column: {col}")
        
        dataset = {
            "type": dataset_type,
            "beta": sumstats["beta"].values,
            "varbeta": (sumstats["se"] ** 2).values,
            "pvalues": sumstats["pval"].values,
            "snp": sumstats["rsid"].values,
            "N": sample_size,
        }
        
        if "maf" in sumstats.columns:
            dataset["MAF"] = sumstats["maf"].values
        
        if dataset_type == "quant" and sdY is not None:
            dataset["sdY"] = sdY
        
        return dataset
    
    def run_coloc(
        self,
        gwas_df: pd.DataFrame,
        qtl_df: pd.DataFrame,
        gwas_type: str = "cc",
        qtl_type: str = "quant",
        gwas_n: int = 100000,
        qtl_n: int = 500,
    ) -> Dict[str, Any]:
        """
        Run COLOC analysis between GWAS and QTL.
        
        Parameters
        ----------
        gwas_df : pd.DataFrame
            GWAS summary statistics.
        qtl_df : pd.DataFrame
            QTL summary statistics.
        gwas_type : str
            GWAS trait type ("cc" or "quant").
        qtl_type : str
            QTL trait type (usually "quant").
        gwas_n : int
            GWAS sample size.
        qtl_n : int
            QTL sample size.
            
        Returns
        -------
        dict
            COLOC results including posterior probabilities.
        """
        # Find common variants
        common_snps = set(gwas_df["rsid"]) & set(qtl_df["rsid"])
        
        if len(common_snps) < 10:
            logger.warning(f"Only {len(common_snps)} common variants, skipping")
            return {
                "pp_h0": np.nan,
                "pp_h1": np.nan,
                "pp_h2": np.nan,
                "pp_h3": np.nan,
                "pp_h4": np.nan,
                "n_snps": len(common_snps),
                "error": "Insufficient variants",
            }
        
        # Filter to common variants
        gwas_common = gwas_df[gwas_df["rsid"].isin(common_snps)].copy()
        qtl_common = qtl_df[qtl_df["rsid"].isin(common_snps)].copy()
        
        # Align by rsid
        gwas_common = gwas_common.set_index("rsid").loc[list(common_snps)].reset_index()
        qtl_common = qtl_common.set_index("rsid").loc[list(common_snps)].reset_index()
        
        # Run via R
        return self._run_coloc_r(
            gwas_common, qtl_common, gwas_type, qtl_type, gwas_n, qtl_n
        )
    
    def _run_coloc_r(
        self,
        gwas_df: pd.DataFrame,
        qtl_df: pd.DataFrame,
        gwas_type: str,
        qtl_type: str,
        gwas_n: int,
        qtl_n: int,
    ) -> Dict[str, Any]:
        """Run COLOC via R."""
        tmpdir = tempfile.mkdtemp()
        
        # Save data
        gwas_file = f"{tmpdir}/gwas.txt"
        qtl_file = f"{tmpdir}/qtl.txt"
        
        gwas_df.to_csv(gwas_file, sep="\t", index=False)
        qtl_df.to_csv(qtl_file, sep="\t", index=False)
        
        # Build R script
        r_script = f'''
library(coloc)

# Load data
gwas <- read.table("{gwas_file}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)
qtl <- read.table("{qtl_file}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

# Build datasets
d1 <- list(
    type = "{gwas_type}",
    beta = gwas$beta,
    varbeta = gwas$se^2,
    pvalues = gwas$pval,
    snp = gwas$rsid,
    N = {gwas_n}
)

d2 <- list(
    type = "{qtl_type}",
    beta = qtl$beta,
    varbeta = qtl$se^2,
    pvalues = qtl$pval,
    snp = qtl$rsid,
    N = {qtl_n}
)

# Add MAF if available
if ("maf" %in% names(gwas)) d1$MAF <- gwas$maf
if ("maf" %in% names(qtl)) d2$MAF <- qtl$maf

# Run COLOC
result <- coloc.abf(
    dataset1 = d1,
    dataset2 = d2,
    p1 = {self.p1},
    p2 = {self.p2},
    p12 = {self.p12}
)

# Save results
pp <- result$summary
write.table(
    t(pp),
    "{tmpdir}/coloc_result.txt",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

# Save per-SNP posteriors
snp_results <- result$results
write.table(
    snp_results,
    "{tmpdir}/coloc_snps.txt",
    row.names = FALSE,
    quote = FALSE,
    sep = "\\t"
)
'''
        
        script_file = f"{tmpdir}/coloc.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        result = subprocess.run(
            ["Rscript", script_file],
            capture_output=True,
            text=True,
        )
        
        if result.returncode != 0:
            logger.error(f"COLOC error: {result.stderr}")
            return {
                "pp_h0": np.nan,
                "pp_h1": np.nan,
                "pp_h2": np.nan,
                "pp_h3": np.nan,
                "pp_h4": np.nan,
                "error": result.stderr,
            }
        
        # Read results
        pp_df = pd.read_csv(f"{tmpdir}/coloc_result.txt", sep=r'\s+')
        snp_df = pd.read_csv(f"{tmpdir}/coloc_snps.txt", sep="\t")
        
        return {
            "pp_h0": pp_df.iloc[0, 0] if len(pp_df.columns) > 0 else np.nan,
            "pp_h1": pp_df.iloc[0, 1] if len(pp_df.columns) > 1 else np.nan,
            "pp_h2": pp_df.iloc[0, 2] if len(pp_df.columns) > 2 else np.nan,
            "pp_h3": pp_df.iloc[0, 3] if len(pp_df.columns) > 3 else np.nan,
            "pp_h4": pp_df.iloc[0, 4] if len(pp_df.columns) > 4 else np.nan,
            "n_snps": len(gwas_df),
            "snp_posteriors": snp_df.to_dict("records"),
        }
    
    def run_susie_coloc(
        self,
        gwas_df: pd.DataFrame,
        qtl_df: pd.DataFrame,
        gwas_ld: np.ndarray,
        qtl_ld: np.ndarray,
        gwas_type: str = "cc",
        qtl_type: str = "quant",
        gwas_n: int = 100000,
        qtl_n: int = 500,
        gwas_s: Optional[float] = None,  # case proportion for cc
        susie_args: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        """
        Run SuSiE-based colocalization (coloc.susie) per Wallace 2021.
        
        This method properly handles multiple causal variants by:
        1. Running SuSiE on each dataset to identify credible sets
        2. Computing pairwise colocalization between all CS pairs
        3. Returning sensitivity analysis (coloc.susie vs coloc.abf)
        
        Key difference from coloc.abf: Does NOT assume single causal variant.
        Internal consistency with SuSiE fine-mapping is maintained.
        
        Reference: Wallace C (2021) PLOS Genetics. doi:10.1371/journal.pgen.1009440
        
        Parameters
        ----------
        gwas_df : pd.DataFrame
            GWAS summary statistics with columns: rsid, beta, se, maf.
        qtl_df : pd.DataFrame
            QTL summary statistics with columns: rsid, beta, se, maf.
        gwas_ld : np.ndarray
            LD matrix for GWAS variants (same order as gwas_df).
        qtl_ld : np.ndarray
            LD matrix for QTL variants (same order as qtl_df).
        gwas_type : str
            "cc" for case-control, "quant" for quantitative.
        qtl_type : str
            Usually "quant" for expression QTLs.
        gwas_n : int
            GWAS sample size.
        qtl_n : int
            QTL sample size.
        gwas_s : float, optional
            Case proportion for case-control GWAS.
        susie_args : dict, optional
            Additional arguments for SuSiE (L, coverage, min_abs_corr).
            
        Returns
        -------
        dict
            Comprehensive results including:
            - summary: Overall colocalization summary
            - credible_sets: All CS pairs with PP.H4
            - sensitivity: Comparison with coloc.abf
            - best_coloc: Best colocalized CS pair
        """
        susie_args = susie_args or {}
        
        # Find common variants
        common_snps = list(set(gwas_df["rsid"]) & set(qtl_df["rsid"]))
        
        if len(common_snps) < 50:
            logger.warning(f"Only {len(common_snps)} common variants, insufficient for SuSiE")
            return self._error_result("Insufficient common variants (<50)")
        
        # Align dataframes
        gwas_aligned, qtl_aligned, ld1, ld2 = self._align_datasets(
            gwas_df, qtl_df, gwas_ld, qtl_ld, common_snps
        )
        
        # Run via R coloc.susie
        return self._run_susie_coloc_r(
            gwas_aligned, qtl_aligned, ld1, ld2,
            gwas_type, qtl_type, gwas_n, qtl_n, gwas_s,
            susie_args
        )
    
    def _align_datasets(
        self,
        gwas_df: pd.DataFrame,
        qtl_df: pd.DataFrame,
        gwas_ld: np.ndarray,
        qtl_ld: np.ndarray,
        common_snps: List[str],
    ) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
        """Align datasets to common variants."""
        # Index by rsid
        gwas_idx = {r: i for i, r in enumerate(gwas_df["rsid"])}
        qtl_idx = {r: i for i, r in enumerate(qtl_df["rsid"])}
        
        # Get common indices
        gwas_common_idx = [gwas_idx[s] for s in common_snps]
        qtl_common_idx = [qtl_idx[s] for s in common_snps]
        
        # Subset dataframes
        gwas_aligned = gwas_df.iloc[gwas_common_idx].reset_index(drop=True)
        qtl_aligned = qtl_df.iloc[qtl_common_idx].reset_index(drop=True)
        
        # Subset LD matrices
        ld1 = gwas_ld[np.ix_(gwas_common_idx, gwas_common_idx)]
        ld2 = qtl_ld[np.ix_(qtl_common_idx, qtl_common_idx)]
        
        return gwas_aligned, qtl_aligned, ld1, ld2
    
    def _run_susie_coloc_r(
        self,
        gwas_df: pd.DataFrame,
        qtl_df: pd.DataFrame,
        ld1: np.ndarray,
        ld2: np.ndarray,
        gwas_type: str,
        qtl_type: str,
        gwas_n: int,
        qtl_n: int,
        gwas_s: Optional[float],
        susie_args: Dict,
    ) -> Dict[str, Any]:
        """Run coloc.susie via R."""
        tmpdir = tempfile.mkdtemp()
        
        # Save data files
        gwas_file = f"{tmpdir}/gwas.txt"
        qtl_file = f"{tmpdir}/qtl.txt"
        ld1_file = f"{tmpdir}/ld1.txt"
        ld2_file = f"{tmpdir}/ld2.txt"
        
        gwas_df.to_csv(gwas_file, sep="\t", index=False)
        qtl_df.to_csv(qtl_file, sep="\t", index=False)
        np.savetxt(ld1_file, ld1, delimiter="\t")
        np.savetxt(ld2_file, ld2, delimiter="\t")
        
        # SuSiE parameters
        L = susie_args.get("L", 10)
        coverage = susie_args.get("coverage", 0.95)
        min_abs_corr = susie_args.get("min_abs_corr", 0.5)
        
        # Build R script for coloc.susie
        s_param = f"s = {gwas_s}," if gwas_s and gwas_type == "cc" else ""
        
        r_script = f'''
library(coloc)
library(susieR)

# Load data
gwas <- read.table("{gwas_file}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)
qtl <- read.table("{qtl_file}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)
ld1 <- as.matrix(read.table("{ld1_file}", sep="\\t"))
ld2 <- as.matrix(read.table("{ld2_file}", sep="\\t"))

# Ensure LD matrices are positive definite
ld1 <- ld1 + diag(1e-6, nrow(ld1))
ld2 <- ld2 + diag(1e-6, nrow(ld2))

# Build dataset 1 (GWAS)
d1 <- list(
    type = "{gwas_type}",
    beta = gwas$beta,
    varbeta = gwas$se^2,
    snp = gwas$rsid,
    N = {gwas_n},
    {s_param}
    LD = ld1
)
if ("maf" %in% names(gwas)) d1$MAF <- gwas$maf

# Build dataset 2 (QTL)
d2 <- list(
    type = "{qtl_type}",
    beta = qtl$beta,
    varbeta = qtl$se^2,
    snp = qtl$rsid,
    N = {qtl_n},
    LD = ld2
)
if ("maf" %in% names(qtl)) d2$MAF <- qtl$maf

# Run SuSiE on each dataset
tryCatch({{
    s1 <- runsusie(d1, L = {L}, coverage = {coverage}, min_abs_corr = {min_abs_corr})
    s2 <- runsusie(d2, L = {L}, coverage = {coverage}, min_abs_corr = {min_abs_corr})
    
    # Run coloc.susie
    result <- coloc.susie(s1, s2)
    
    # Also run standard coloc.abf for sensitivity comparison
    abf_result <- coloc.abf(
        dataset1 = d1,
        dataset2 = d2,
        p1 = {self.p1},
        p2 = {self.p2},
        p12 = {self.p12}
    )
    
    # Save comprehensive results
    # 1. Summary statistics
    summary_df <- result$summary
    write.table(summary_df, "{tmpdir}/coloc_susie_summary.txt", 
                row.names = FALSE, quote = FALSE, sep = "\\t")
    
    # 2. Credible set information
    if (!is.null(s1$sets$cs)) {{
        cs1_info <- data.frame(
            cs_id = names(s1$sets$cs),
            n_snps = sapply(s1$sets$cs, length),
            coverage = s1$sets$coverage
        )
        write.table(cs1_info, "{tmpdir}/gwas_cs.txt", 
                    row.names = FALSE, quote = FALSE, sep = "\\t")
    }}
    
    if (!is.null(s2$sets$cs)) {{
        cs2_info <- data.frame(
            cs_id = names(s2$sets$cs),
            n_snps = sapply(s2$sets$cs, length),
            coverage = s2$sets$coverage
        )
        write.table(cs2_info, "{tmpdir}/qtl_cs.txt", 
                    row.names = FALSE, quote = FALSE, sep = "\\t")
    }}
    
    # 3. Standard coloc for comparison
    abf_summary <- data.frame(
        pp_h0 = abf_result$summary[1],
        pp_h1 = abf_result$summary[2],
        pp_h2 = abf_result$summary[3],
        pp_h3 = abf_result$summary[4],
        pp_h4 = abf_result$summary[5]
    )
    write.table(abf_summary, "{tmpdir}/coloc_abf_summary.txt", 
                row.names = FALSE, quote = FALSE, sep = "\\t")
    
    # 4. Number of credible sets found
    meta <- data.frame(
        n_gwas_cs = length(s1$sets$cs),
        n_qtl_cs = length(s2$sets$cs),
        status = "success"
    )
    write.table(meta, "{tmpdir}/meta.txt", 
                row.names = FALSE, quote = FALSE, sep = "\\t")
                
}}, error = function(e) {{
    meta <- data.frame(
        n_gwas_cs = 0,
        n_qtl_cs = 0,
        status = paste("error:", e$message)
    )
    write.table(meta, "{tmpdir}/meta.txt", 
                row.names = FALSE, quote = FALSE, sep = "\\t")
}})
'''
        
        script_file = f"{tmpdir}/coloc_susie.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        result = subprocess.run(
            ["Rscript", script_file],
            capture_output=True,
            text=True,
        )
        
        if result.returncode != 0:
            logger.error(f"coloc.susie error: {result.stderr}")
            return self._error_result(f"R error: {result.stderr}")
        
        # Parse results
        return self._parse_susie_coloc_results(tmpdir)
    
    def _parse_susie_coloc_results(self, tmpdir: str) -> Dict[str, Any]:
        """Parse coloc.susie results from R output files."""
        results = {
            "method": "coloc.susie",
            "credible_set_pairs": [],
            "best_coloc": None,
            "sensitivity": {},
        }
        
        # Read meta info
        meta_file = f"{tmpdir}/meta.txt"
        if Path(meta_file).exists():
            meta = pd.read_csv(meta_file, sep="\t")
            results["n_gwas_cs"] = int(meta["n_gwas_cs"].iloc[0])
            results["n_qtl_cs"] = int(meta["n_qtl_cs"].iloc[0])
            
            if "error" in str(meta["status"].iloc[0]):
                return self._error_result(str(meta["status"].iloc[0]))
        
        # Read coloc.susie summary (CS-pair colocalization)
        summary_file = f"{tmpdir}/coloc_susie_summary.txt"
        if Path(summary_file).exists():
            summary = pd.read_csv(summary_file, sep="\t")
            
            for _, row in summary.iterrows():
                pair = {
                    "idx1": int(row.get("idx1", 0)),
                    "idx2": int(row.get("idx2", 0)),
                    "pp_h0": row.get("PP.H0.abf", np.nan),
                    "pp_h1": row.get("PP.H1.abf", np.nan),
                    "pp_h2": row.get("PP.H2.abf", np.nan),
                    "pp_h3": row.get("PP.H3.abf", np.nan),
                    "pp_h4": row.get("PP.H4.abf", np.nan),
                }
                results["credible_set_pairs"].append(pair)
            
            # Find best colocalized pair
            if results["credible_set_pairs"]:
                best = max(results["credible_set_pairs"], 
                          key=lambda x: x.get("pp_h4", 0))
                results["best_coloc"] = best
                results["max_pp_h4"] = best.get("pp_h4", 0)
        
        # Read standard coloc.abf for sensitivity
        abf_file = f"{tmpdir}/coloc_abf_summary.txt"
        if Path(abf_file).exists():
            abf = pd.read_csv(abf_file, sep="\t")
            results["sensitivity"]["coloc_abf"] = {
                "pp_h0": float(abf["pp_h0"].iloc[0]),
                "pp_h1": float(abf["pp_h1"].iloc[0]),
                "pp_h2": float(abf["pp_h2"].iloc[0]),
                "pp_h3": float(abf["pp_h3"].iloc[0]),
                "pp_h4": float(abf["pp_h4"].iloc[0]),
            }
            
            # Flag potential discordance
            if results.get("max_pp_h4", 0) > 0:
                abf_h4 = results["sensitivity"]["coloc_abf"]["pp_h4"]
                susie_h4 = results["max_pp_h4"]
                
                results["sensitivity"]["concordant"] = (
                    (abf_h4 > 0.5 and susie_h4 > 0.5) or
                    (abf_h4 <= 0.5 and susie_h4 <= 0.5)
                )
                results["sensitivity"]["h4_difference"] = susie_h4 - abf_h4
        
        # Interpretation
        max_h4 = results.get("max_pp_h4", 0)
        if max_h4 >= COLOC_THRESHOLDS["strong"]:
            results["interpretation"] = "strong_colocalization"
        elif max_h4 >= COLOC_THRESHOLDS["moderate"]:
            results["interpretation"] = "moderate_colocalization"
        elif max_h4 >= COLOC_THRESHOLDS["weak"]:
            results["interpretation"] = "weak_colocalization"
        else:
            results["interpretation"] = "no_colocalization"
        
        return results
    
    def _error_result(self, message: str) -> Dict[str, Any]:
        """Return error result structure."""
        return {
            "method": "coloc.susie",
            "error": message,
            "max_pp_h4": np.nan,
            "credible_set_pairs": [],
            "interpretation": "error",
        }
    
    def run_coloc_sensitivity(
        self,
        gwas_df: pd.DataFrame,
        qtl_df: pd.DataFrame,
        gwas_ld: np.ndarray,
        qtl_ld: np.ndarray,
        **kwargs,
    ) -> Dict[str, Any]:
        """
        Run both coloc.abf and coloc.susie for sensitivity analysis.
        
        This is the recommended approach for publication: report both methods
        and highlight any discordance which may indicate multiple causal signals.
        
        Parameters
        ----------
        gwas_df, qtl_df : pd.DataFrame
            Summary statistics.
        gwas_ld, qtl_ld : np.ndarray
            LD matrices.
        **kwargs
            Additional arguments for both methods.
            
        Returns
        -------
        dict
            Combined results with sensitivity comparison.
        """
        # Run coloc.abf (original method)
        abf_result = self.run_coloc(
            gwas_df, qtl_df,
            gwas_type=kwargs.get("gwas_type", "cc"),
            qtl_type=kwargs.get("qtl_type", "quant"),
            gwas_n=kwargs.get("gwas_n", 100000),
            qtl_n=kwargs.get("qtl_n", 500),
        )
        
        # Run coloc.susie (multi-signal method)
        susie_result = self.run_susie_coloc(
            gwas_df, qtl_df, gwas_ld, qtl_ld, **kwargs
        )
        
        # Combine results
        return {
            "coloc_abf": abf_result,
            "coloc_susie": susie_result,
            "recommended_method": "coloc_susie",
            "sensitivity_analysis": {
                "abf_h4": abf_result.get("pp_h4", np.nan),
                "susie_max_h4": susie_result.get("max_pp_h4", np.nan),
                "n_gwas_signals": susie_result.get("n_gwas_cs", 0),
                "n_qtl_signals": susie_result.get("n_qtl_cs", 0),
                "concordant": susie_result.get("sensitivity", {}).get("concordant", None),
                "note": (
                    "coloc.susie is recommended when SuSiE finds >1 credible set, "
                    "as coloc.abf assumes single causal variant."
                ),
            },
        }


def run_coloc(
    gwas_df: pd.DataFrame,
    qtl_df: pd.DataFrame,
    **kwargs,
) -> Dict[str, Any]:
    """
    Convenience function to run COLOC.
    
    Parameters
    ----------
    gwas_df : pd.DataFrame
        GWAS summary statistics.
    qtl_df : pd.DataFrame
        QTL summary statistics.
    **kwargs
        Arguments for ColocAnalysis.
        
    Returns
    -------
    dict
        COLOC results.
    """
    analyzer = ColocAnalysis()
    return analyzer.run_coloc(gwas_df, qtl_df, **kwargs)
