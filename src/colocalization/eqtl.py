"""
eQTL Data Loading and Processing

Handles loading and querying eQTL/sQTL data from GTEx and eQTL Catalogue.
"""

import gzip
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import pandas as pd
import requests

from ..utils.config import get_config
from ..utils.logging import get_logger


logger = get_logger("eqtl")


class EQTLLoader:
    """
    Loads and queries eQTL data from various sources.
    
    Supports:
    - GTEx v8 eQTLs and sQTLs
    - eQTL Catalogue fine-mapped data
    """
    
    def __init__(
        self,
        data_dir: str,
        source: str = "gtex",
    ):
        """
        Initialize eQTL loader.
        
        Parameters
        ----------
        data_dir : str
            Path to eQTL data directory.
        source : str
            Data source: "gtex" or "eqtl_catalogue".
        """
        self.data_dir = Path(data_dir)
        self.source = source
        
        # Cache loaded data
        self._cache = {}
    
    def list_tissues(self) -> List[str]:
        """
        List available tissues.
        
        Returns
        -------
        list
            Available tissue names.
        """
        if self.source == "gtex":
            # GTEx tissue directory structure
            tissues = []
            tissue_dir = self.data_dir / "GTEx_Analysis_v8_eQTL"
            
            if tissue_dir.exists():
                for f in tissue_dir.glob("*.v8.egenes.txt.gz"):
                    tissue = f.name.replace(".v8.egenes.txt.gz", "")
                    tissues.append(tissue)
            
            return sorted(tissues)
        
        elif self.source == "eqtl_catalogue":
            # eQTL Catalogue uses dataset IDs
            datasets = []
            for f in self.data_dir.glob("*.credible_sets.tsv.gz"):
                dataset = f.name.replace(".credible_sets.tsv.gz", "")
                datasets.append(dataset)
            
            return sorted(datasets)
        
        return []
    
    def load_tissue(
        self,
        tissue: str,
        genes: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Load eQTL data for a tissue.
        
        Parameters
        ----------
        tissue : str
            Tissue name.
        genes : list, optional
            Filter to specific genes.
            
        Returns
        -------
        pd.DataFrame
            eQTL summary statistics.
        """
        cache_key = f"{tissue}_{hash(tuple(genes or []))}"
        
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        if self.source == "gtex":
            df = self._load_gtex_tissue(tissue, genes)
        elif self.source == "eqtl_catalogue":
            df = self._load_eqtl_catalogue(tissue, genes)
        else:
            raise ValueError(f"Unknown source: {self.source}")
        
        self._cache[cache_key] = df
        return df
    
    def _load_gtex_tissue(
        self,
        tissue: str,
        genes: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """Load GTEx eQTL data for a tissue."""
        # Path to allpairs file (full summary statistics)
        allpairs_file = (
            self.data_dir / "GTEx_Analysis_v8_eQTL" /
            f"{tissue}.v8.signif_variant_gene_pairs.txt.gz"
        )
        
        if not allpairs_file.exists():
            logger.warning(f"GTEx file not found: {allpairs_file}")
            return pd.DataFrame()
        
        logger.info(f"Loading GTEx eQTLs for {tissue}")
        
        # Read with chunking for large files
        chunks = []
        for chunk in pd.read_csv(
            allpairs_file,
            sep="\t",
            chunksize=100000,
        ):
            # Filter by genes if specified
            if genes:
                chunk = chunk[chunk["gene_id"].str.split(".").str[0].isin(genes)]
            
            chunks.append(chunk)
        
        df = pd.concat(chunks, ignore_index=True)
        
        # Standardize columns
        df = self._standardize_gtex(df)
        
        logger.info(f"Loaded {len(df):,} eQTL associations")
        
        return df
    
    def _standardize_gtex(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize GTEx column names."""
        # GTEx uses format: chr1_123456_A_G_b38
        if "variant_id" in df.columns:
            variant_parts = df["variant_id"].str.split("_", expand=True)
            df["chr"] = variant_parts[0].str.replace("chr", "")
            df["pos"] = variant_parts[1].astype(int)
            df["ref"] = variant_parts[2]
            df["alt"] = variant_parts[3]
        
        # Rename columns
        rename_map = {
            "slope": "beta",
            "slope_se": "se",
            "pval_nominal": "pval",
            "maf": "maf",
            "gene_id": "gene_id",
        }
        
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
        
        # Create rsid if not present
        if "rsid" not in df.columns and "variant_id" in df.columns:
            df["rsid"] = df["variant_id"]
        
        # Strip version from gene IDs
        if "gene_id" in df.columns:
            df["gene_ensembl"] = df["gene_id"].str.split(".").str[0]
        
        return df
    
    def _load_eqtl_catalogue(
        self,
        dataset: str,
        genes: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """Load eQTL Catalogue credible sets."""
        cs_file = self.data_dir / f"{dataset}.credible_sets.tsv.gz"
        
        if not cs_file.exists():
            logger.warning(f"eQTL Catalogue file not found: {cs_file}")
            return pd.DataFrame()
        
        logger.info(f"Loading eQTL Catalogue data for {dataset}")
        
        df = pd.read_csv(cs_file, sep="\t")
        
        # Filter by genes if specified
        if genes and "molecular_trait_id" in df.columns:
            df = df[df["molecular_trait_id"].str.split(".").str[0].isin(genes)]
        
        # Standardize columns
        rename_map = {
            "chromosome": "chr",
            "position": "pos",
            "variant": "rsid",
            "molecular_trait_id": "gene_id",
            "pip": "pip",
            "cs_id": "cs_id",
        }
        
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
        
        return df
    
    def query_region(
        self,
        tissue: str,
        chrom: str,
        start: int,
        end: int,
    ) -> pd.DataFrame:
        """
        Query eQTLs in a genomic region.
        
        Parameters
        ----------
        tissue : str
            Tissue name.
        chrom : str
            Chromosome.
        start : int
            Start position.
        end : int
            End position.
            
        Returns
        -------
        pd.DataFrame
            eQTLs in the region.
        """
        df = self.load_tissue(tissue)
        
        chrom = str(chrom).replace("chr", "")
        
        region_df = df[
            (df["chr"].astype(str) == chrom) &
            (df["pos"] >= start) &
            (df["pos"] <= end)
        ]
        
        return region_df
    
    def get_gene_eqtls(
        self,
        tissue: str,
        gene: str,
    ) -> pd.DataFrame:
        """
        Get all eQTLs for a gene.
        
        Parameters
        ----------
        tissue : str
            Tissue name.
        gene : str
            Gene ID (Ensembl ID without version).
            
        Returns
        -------
        pd.DataFrame
            eQTL data for the gene.
        """
        df = self.load_tissue(tissue, genes=[gene])
        return df


def query_eqtl_catalogue(
    chrom: str,
    start: int,
    end: int,
    dataset: str = "QTS000001",
    qtl_type: str = "eQTL",
) -> pd.DataFrame:
    """
    Query eQTL Catalogue API for a region.
    
    Parameters
    ----------
    chrom : str
        Chromosome.
    start : int
        Start position.
    end : int
        End position.
    dataset : str
        Dataset ID (e.g., "QTS000001").
    qtl_type : str
        QTL type: "eQTL" or "sQTL".
        
    Returns
    -------
    pd.DataFrame
        eQTL associations in the region.
    """
    base_url = "https://www.ebi.ac.uk/eqtl/api"
    
    # Build query
    chrom = str(chrom).replace("chr", "")
    
    url = (
        f"{base_url}/associations"
        f"?study={dataset}"
        f"&qtl_type={qtl_type}"
        f"&chromosome={chrom}"
        f"&position_start={start}"
        f"&position_end={end}"
    )
    
    logger.info(f"Querying eQTL Catalogue: {url}")
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        
        if not data:
            return pd.DataFrame()
        
        df = pd.DataFrame(data)
        
        # Standardize columns
        rename_map = {
            "chromosome": "chr",
            "position": "pos",
            "variant_id": "rsid",
            "gene_id": "gene_id",
            "beta": "beta",
            "se": "se",
            "pvalue": "pval",
        }
        
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
        
        return df
        
    except Exception as e:
        logger.error(f"eQTL Catalogue query failed: {e}")
        return pd.DataFrame()


class MultiTissueEQTL:
    """
    Handles multi-tissue eQTL analysis for a locus.
    """
    
    def __init__(
        self,
        data_dir: str,
        tissues: List[str],
    ):
        """
        Initialize multi-tissue loader.
        
        Parameters
        ----------
        data_dir : str
            Path to eQTL data.
        tissues : list
            List of tissues to include.
        """
        self.loader = EQTLLoader(data_dir, source="gtex")
        self.tissues = tissues
    
    def query_locus_all_tissues(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Dict[str, pd.DataFrame]:
        """
        Query eQTLs across all tissues for a locus.
        
        Parameters
        ----------
        chrom : str
            Chromosome.
        start : int
            Start position.
        end : int
            End position.
            
        Returns
        -------
        dict
            Dictionary of tissue -> eQTL DataFrame.
        """
        results = {}
        
        for tissue in self.tissues:
            df = self.loader.query_region(tissue, chrom, start, end)
            if len(df) > 0:
                results[tissue] = df
        
        logger.info(f"Found eQTLs in {len(results)}/{len(self.tissues)} tissues")
        
        return results
    
    def get_tissue_specific_genes(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Dict[str, List[str]]:
        """
        Get genes with eQTLs in each tissue.
        
        Returns
        -------
        dict
            Dictionary of tissue -> list of genes.
        """
        tissue_eqtls = self.query_locus_all_tissues(chrom, start, end)
        
        tissue_genes = {}
        for tissue, df in tissue_eqtls.items():
            if "gene_ensembl" in df.columns:
                genes = df["gene_ensembl"].unique().tolist()
            elif "gene_id" in df.columns:
                genes = df["gene_id"].str.split(".").str[0].unique().tolist()
            else:
                genes = []
            
            tissue_genes[tissue] = genes
        
        return tissue_genes
