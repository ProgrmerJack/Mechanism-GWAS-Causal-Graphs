"""
eQTL Catalogue Integration

Cross-study replication of eQTL signals using the eQTL Catalogue.

Key principle: Discover in GTEx, replicate in eQTL Catalogue.
- Single-study QTLs (GTEx only) raise concern about generalizability
- Cross-dataset replication dramatically increases credibility

References:
- Kerimov et al. (2021) Nature Genetics: eQTL Catalogue
- GTEx Consortium (2020) Science: GTEx v8
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple, Union
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from ..utils.genomics import liftover_coordinates


logger = get_logger("eqtl_catalogue")


# eQTL Catalogue tissue/cell type mapping to GTEx
EQTL_CATALOGUE_TISSUE_MAP = {
    # Blood
    "Whole_Blood": ["BLUEPRINT_SE", "BLUEPRINT_PE", "Alasoo_2018_macrophage_naive"],
    # Adipose
    "Adipose_Subcutaneous": ["FUSION_adipose"],
    "Adipose_Visceral_Omentum": ["FUSION_adipose"],
    # Heart
    "Heart_Left_Ventricle": ["GTEx_v8_heart_left_ventricle"],
    "Heart_Atrial_Appendage": [],
    # Liver
    "Liver": ["GTEx_v8_liver"],
    # Muscle
    "Muscle_Skeletal": ["FUSION_muscle"],
    # Pancreas
    "Pancreas": ["FUSION_pancreas"],
    # Artery
    "Artery_Aorta": [],
    "Artery_Coronary": [],
    "Artery_Tibial": [],
}


@dataclass
class EQTLSignal:
    """
    An eQTL signal from any source.
    
    Attributes
    ----------
    gene_id : str
        ENSG gene ID.
    gene_symbol : str
        Gene symbol.
    variant_id : str
        Variant ID (chr_pos_ref_alt).
    tissue : str
        Tissue/cell type.
    beta : float
        Effect size.
    se : float
        Standard error.
    pvalue : float
        P-value.
    source : str
        Data source (e.g., "GTEx_v8", "eQTL_Catalogue").
    study : str
        Specific study within source.
    sample_size : int
        Sample size.
    maf : float
        Minor allele frequency.
    """
    
    gene_id: str
    gene_symbol: str = ""
    variant_id: str = ""
    tissue: str = ""
    beta: float = 0.0
    se: float = 0.0
    pvalue: float = 1.0
    source: str = ""
    study: str = ""
    sample_size: int = 0
    maf: float = 0.0
    
    @property
    def zscore(self) -> float:
        """Compute Z-score."""
        if self.se > 0:
            return self.beta / self.se
        return 0.0


@dataclass
class ReplicationResult:
    """
    Result of cross-study replication test.
    
    Attributes
    ----------
    gene_id : str
        Gene ID.
    tissue : str
        Tissue for discovery.
    discovery_source : str
        Discovery study source.
    replication_source : str
        Replication study source.
    discovery_pvalue : float
        Discovery P-value.
    replication_pvalue : float
        Replication P-value.
    concordant_direction : bool
        Whether effect directions match.
    replicated : bool
        Whether signal replicated.
    replication_r2 : float
        Correlation of effect sizes.
    n_variants_tested : int
        Number of variants tested.
    """
    
    gene_id: str
    tissue: str = ""
    discovery_source: str = ""
    replication_source: str = ""
    discovery_pvalue: float = 1.0
    replication_pvalue: float = 1.0
    concordant_direction: bool = False
    replicated: bool = False
    replication_r2: float = 0.0
    n_variants_tested: int = 0


class EQTLCatalogueLoader:
    """
    Load eQTL data from the eQTL Catalogue.
    
    The eQTL Catalogue provides uniformly processed eQTL data from
    multiple studies, enabling cross-study replication.
    
    Example
    -------
    >>> loader = EQTLCatalogueLoader(cache_dir="data/eqtl_catalogue")
    >>> 
    >>> # Load blood eQTLs
    >>> blood_eqtls = loader.load_dataset("BLUEPRINT_SE")
    >>> 
    >>> # Query specific gene
    >>> gene_eqtls = loader.query_gene("ENSG00000012048")  # BRCA1
    """
    
    BASE_URL = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master"
    
    AVAILABLE_DATASETS = [
        "BLUEPRINT_SE",
        "BLUEPRINT_PE",
        "Alasoo_2018",
        "FUSION",
        "GEUVADIS",
        "HipSci",
        "Lepik_2017",
        "Nedelec_2016",
        "Quach_2016",
        "Schmiedel_2018",
        "van_de_Bunt_2015",
    ]
    
    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        use_tabix: bool = True,
    ):
        """
        Initialize loader.
        
        Parameters
        ----------
        cache_dir : str or Path, optional
            Directory for caching downloaded data.
        use_tabix : bool
            Whether to use tabix for indexed queries.
        """
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.use_tabix = use_tabix
        self._loaded_datasets: Dict[str, pd.DataFrame] = {}
    
    def load_dataset(
        self,
        dataset: str,
        tissue: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Load an eQTL Catalogue dataset.
        
        Parameters
        ----------
        dataset : str
            Dataset name (e.g., "BLUEPRINT_SE").
        tissue : str, optional
            Specific tissue/cell type within dataset.
            
        Returns
        -------
        pd.DataFrame
            eQTL summary statistics.
        """
        cache_key = f"{dataset}_{tissue}" if tissue else dataset
        
        if cache_key in self._loaded_datasets:
            return self._loaded_datasets[cache_key]
        
        # In production, this would download from eQTL Catalogue
        # For now, return empty DataFrame with expected schema
        logger.info(f"Loading eQTL Catalogue dataset: {dataset}")
        
        df = pd.DataFrame(columns=[
            "gene_id",
            "gene_symbol",
            "variant_id",
            "chromosome",
            "position",
            "ref",
            "alt",
            "beta",
            "se",
            "pvalue",
            "maf",
            "sample_size",
        ])
        
        self._loaded_datasets[cache_key] = df
        return df
    
    def query_gene(
        self,
        gene_id: str,
        datasets: Optional[List[str]] = None,
    ) -> List[EQTLSignal]:
        """
        Query eQTLs for a specific gene across datasets.
        
        Parameters
        ----------
        gene_id : str
            ENSG gene ID.
        datasets : list, optional
            Datasets to query. If None, queries all.
            
        Returns
        -------
        list
            List of EQTLSignal objects.
        """
        if datasets is None:
            datasets = self.AVAILABLE_DATASETS
        
        signals = []
        
        for dataset in datasets:
            df = self.load_dataset(dataset)
            
            gene_df = df[df["gene_id"] == gene_id]
            
            for _, row in gene_df.iterrows():
                signals.append(EQTLSignal(
                    gene_id=row["gene_id"],
                    gene_symbol=row.get("gene_symbol", ""),
                    variant_id=row["variant_id"],
                    tissue=dataset,
                    beta=row["beta"],
                    se=row["se"],
                    pvalue=row["pvalue"],
                    source="eQTL_Catalogue",
                    study=dataset,
                    sample_size=row.get("sample_size", 0),
                    maf=row.get("maf", 0.0),
                ))
        
        return signals
    
    def query_region(
        self,
        chromosome: str,
        start: int,
        end: int,
        dataset: str,
    ) -> pd.DataFrame:
        """
        Query eQTLs in a genomic region.
        
        Parameters
        ----------
        chromosome : str
            Chromosome.
        start : int
            Start position.
        end : int
            End position.
        dataset : str
            Dataset to query.
            
        Returns
        -------
        pd.DataFrame
            eQTLs in region.
        """
        df = self.load_dataset(dataset)
        
        mask = (
            (df["chromosome"] == chromosome) &
            (df["position"] >= start) &
            (df["position"] <= end)
        )
        
        return df[mask].copy()


class CrossStudyReplication:
    """
    Cross-study replication analysis for eQTL signals.
    
    Discovery in GTEx → Replication in eQTL Catalogue.
    
    This addresses reviewer concern about single-study effects:
    "GTEx-only effects" may not replicate in independent cohorts.
    
    Example
    -------
    >>> replicator = CrossStudyReplication()
    >>> 
    >>> # Test if a GTEx eQTL replicates
    >>> result = replicator.test_replication(
    ...     gene_id="ENSG00000012048",
    ...     tissue="Whole_Blood",
    ...     discovery_variants=gtex_variants,
    ...     discovery_betas=gtex_betas,
    ... )
    >>> 
    >>> if result.replicated:
    ...     print(f"Gene replicates with r²={result.replication_r2:.2f}")
    """
    
    def __init__(
        self,
        eqtl_catalogue: Optional[EQTLCatalogueLoader] = None,
        replication_pvalue_threshold: float = 0.05,
        require_direction_concordance: bool = True,
    ):
        """
        Initialize replication tester.
        
        Parameters
        ----------
        eqtl_catalogue : EQTLCatalogueLoader, optional
            eQTL Catalogue loader.
        replication_pvalue_threshold : float
            P-value threshold for replication.
        require_direction_concordance : bool
            Whether to require concordant effect directions.
        """
        self.eqtl_catalogue = eqtl_catalogue or EQTLCatalogueLoader()
        self.replication_pvalue_threshold = replication_pvalue_threshold
        self.require_direction_concordance = require_direction_concordance
    
    def test_replication(
        self,
        gene_id: str,
        tissue: str,
        discovery_variants: List[str],
        discovery_betas: List[float],
        discovery_pvalues: Optional[List[float]] = None,
    ) -> ReplicationResult:
        """
        Test if discovery eQTL signal replicates.
        
        Parameters
        ----------
        gene_id : str
            Gene ID.
        tissue : str
            GTEx tissue name.
        discovery_variants : list
            Variant IDs from discovery.
        discovery_betas : list
            Effect sizes from discovery.
        discovery_pvalues : list, optional
            P-values from discovery.
            
        Returns
        -------
        ReplicationResult
            Replication test result.
        """
        # Get matching datasets from eQTL Catalogue
        matching_datasets = EQTL_CATALOGUE_TISSUE_MAP.get(tissue, [])
        
        if not matching_datasets:
            logger.warning(f"No matching eQTL Catalogue datasets for tissue: {tissue}")
            return ReplicationResult(
                gene_id=gene_id,
                tissue=tissue,
                discovery_source="GTEx_v8",
                replication_source="none",
                replicated=False,
            )
        
        # Get replication signals
        replication_signals = []
        for dataset in matching_datasets:
            signals = self.eqtl_catalogue.query_gene(gene_id, [dataset])
            replication_signals.extend(signals)
        
        if not replication_signals:
            return ReplicationResult(
                gene_id=gene_id,
                tissue=tissue,
                discovery_source="GTEx_v8",
                replication_source=",".join(matching_datasets),
                replicated=False,
            )
        
        # Match variants between discovery and replication
        discovery_betas_dict = dict(zip(discovery_variants, discovery_betas))
        
        matched_discovery = []
        matched_replication = []
        
        for signal in replication_signals:
            if signal.variant_id in discovery_betas_dict:
                matched_discovery.append(discovery_betas_dict[signal.variant_id])
                matched_replication.append(signal.beta)
        
        if len(matched_discovery) < 2:
            return ReplicationResult(
                gene_id=gene_id,
                tissue=tissue,
                discovery_source="GTEx_v8",
                replication_source=",".join(matching_datasets),
                replicated=False,
                n_variants_tested=len(matched_discovery),
            )
        
        # Compute correlation
        from scipy.stats import pearsonr, spearmanr
        
        matched_discovery = np.array(matched_discovery)
        matched_replication = np.array(matched_replication)
        
        r, pvalue = pearsonr(matched_discovery, matched_replication)
        r2 = r ** 2
        
        # Check direction concordance
        concordant = np.sign(matched_discovery) == np.sign(matched_replication)
        concordance_rate = concordant.mean()
        
        # Determine if replicated
        replicated = (
            pvalue < self.replication_pvalue_threshold and
            (not self.require_direction_concordance or concordance_rate > 0.5)
        )
        
        return ReplicationResult(
            gene_id=gene_id,
            tissue=tissue,
            discovery_source="GTEx_v8",
            replication_source=",".join(matching_datasets),
            discovery_pvalue=discovery_pvalues[0] if discovery_pvalues else 1.0,
            replication_pvalue=pvalue,
            concordant_direction=(concordance_rate > 0.5),
            replicated=replicated,
            replication_r2=r2,
            n_variants_tested=len(matched_discovery),
        )
    
    def batch_replication_test(
        self,
        genes: List[str],
        tissue: str,
        discovery_data: Dict[str, Dict[str, Any]],
    ) -> pd.DataFrame:
        """
        Test replication for multiple genes.
        
        Parameters
        ----------
        genes : list
            Gene IDs.
        tissue : str
            GTEx tissue.
        discovery_data : dict
            Gene ID -> {"variants": [...], "betas": [...], ...}.
            
        Returns
        -------
        pd.DataFrame
            Replication results.
        """
        results = []
        
        for gene_id in genes:
            if gene_id not in discovery_data:
                continue
            
            data = discovery_data[gene_id]
            result = self.test_replication(
                gene_id=gene_id,
                tissue=tissue,
                discovery_variants=data.get("variants", []),
                discovery_betas=data.get("betas", []),
                discovery_pvalues=data.get("pvalues"),
            )
            
            results.append({
                "gene_id": result.gene_id,
                "tissue": result.tissue,
                "replicated": result.replicated,
                "replication_r2": result.replication_r2,
                "concordant_direction": result.concordant_direction,
                "n_variants": result.n_variants_tested,
                "replication_pvalue": result.replication_pvalue,
            })
        
        return pd.DataFrame(results)
    
    def compute_replication_rate(
        self,
        results: List[ReplicationResult],
    ) -> Dict[str, float]:
        """
        Compute replication statistics.
        
        Parameters
        ----------
        results : list
            ReplicationResult objects.
            
        Returns
        -------
        dict
            Replication statistics.
        """
        if not results:
            return {"replication_rate": 0.0, "n_tested": 0}
        
        n_tested = len(results)
        n_replicated = sum(1 for r in results if r.replicated)
        n_concordant = sum(1 for r in results if r.concordant_direction)
        
        return {
            "replication_rate": n_replicated / n_tested,
            "concordance_rate": n_concordant / n_tested,
            "n_tested": n_tested,
            "n_replicated": n_replicated,
            "mean_r2": np.mean([r.replication_r2 for r in results if r.replication_r2 > 0]),
        }


def match_tissues(
    gtex_tissue: str,
) -> List[str]:
    """
    Get matching eQTL Catalogue datasets for a GTEx tissue.
    
    Parameters
    ----------
    gtex_tissue : str
        GTEx tissue name.
        
    Returns
    -------
    list
        Matching eQTL Catalogue dataset names.
    """
    return EQTL_CATALOGUE_TISSUE_MAP.get(gtex_tissue, [])
