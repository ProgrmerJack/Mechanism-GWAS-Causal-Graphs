"""
Benchmark Suite

Three-tier benchmark structure with explicit provenance and anti-leak provisions:

Tier 1: Holdout Curated Genes (50-100)
    - Mendelian disease genes from OMIM
    - Drug targets with approved indications
    - Anti-leak: Excluded from any training data (L2G, ABC, etc.)

Tier 2: Drug Target Validation (100-200)
    - ChEMBL/OpenTargets drug targets
    - Phase 2+ clinical trial evidence
    - Genetic support from independent studies

Tier 3: CRISPR Perturbation Ground Truth (30-50)
    - CRISPRi/CRISPRa validated enhancer-gene pairs
    - Perturbation sequencing validation
    - Highest confidence, smallest set

Reference for benchmark design:
- Weeks et al. (2023) Leveraging polygenic enrichments
- Mountjoy et al. (2021) Open Targets Genetics
- Fulco et al. (2019) CRISPRi validation

Critical: This module implements anti-leak provisions to ensure benchmark genes
are NOT in training sets for:
- Open Targets L2G model
- ABC model training
- Any method we compare against
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from enum import Enum

import numpy as np
import pandas as pd

from ..utils.config import get_config
from ..utils.logging import get_logger


logger = get_logger("benchmarks")


class BenchmarkTier(Enum):
    """Benchmark tier classification."""
    TIER1_HOLDOUT = "tier1_holdout"      # Curated, anti-leak holdout
    TIER2_DRUG = "tier2_drug"             # Drug target validation
    TIER3_CRISPR = "tier3_crispr"         # CRISPR ground truth


@dataclass
class BenchmarkGene:
    """
    A benchmark gene with full provenance.
    
    Attributes
    ----------
    symbol : str
        Gene symbol.
    ensembl_id : str
        Ensembl gene ID.
    tier : BenchmarkTier
        Benchmark tier.
    trait : str
        Associated trait.
    evidence_type : str
        Type of evidence (mendelian, drug, crispr, etc.).
    evidence_source : str
        Source database/publication.
    evidence_id : str
        Source-specific identifier (OMIM ID, ChEMBL ID, etc.).
    confidence : float
        Confidence score (0-1).
    is_in_l2g_training : bool
        Whether gene was in L2G training (LEAK FLAG).
    is_in_abc_training : bool
        Whether gene was in ABC training (LEAK FLAG).
    provenance : dict
        Full provenance information.
    """
    
    symbol: str
    ensembl_id: str = ""
    tier: BenchmarkTier = BenchmarkTier.TIER1_HOLDOUT
    trait: str = ""
    evidence_type: str = ""
    evidence_source: str = ""
    evidence_id: str = ""
    confidence: float = 1.0
    is_in_l2g_training: bool = False
    is_in_abc_training: bool = False
    provenance: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def has_leak_flag(self) -> bool:
        """Check if gene has potential training data leak."""
        return self.is_in_l2g_training or self.is_in_abc_training
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "symbol": self.symbol,
            "ensembl_id": self.ensembl_id,
            "tier": self.tier.value,
            "trait": self.trait,
            "evidence_type": self.evidence_type,
            "evidence_source": self.evidence_source,
            "evidence_id": self.evidence_id,
            "confidence": self.confidence,
            "is_in_l2g_training": self.is_in_l2g_training,
            "is_in_abc_training": self.is_in_abc_training,
            "has_leak_flag": self.has_leak_flag,
            "provenance": self.provenance,
        }


class GoldStandardGenes:
    """
    Manages gold standard gene sets for benchmarking.
    
    Sources:
    - Mendelian disease genes (OMIM)
    - Drug targets (ChEMBL, OpenTargets)
    - Mouse knockout phenotypes (MGI)
    - Literature-curated sets
    """
    
    # Known lipid genes from Mendelian disorders
    LIPID_GENES = {
        "LDLR": {"mechanism": "receptor", "evidence": "mendelian", "score": 1.0},
        "APOB": {"mechanism": "ligand", "evidence": "mendelian", "score": 1.0},
        "PCSK9": {"mechanism": "degradation", "evidence": "mendelian", "score": 1.0},
        "ABCA1": {"mechanism": "efflux", "evidence": "mendelian", "score": 1.0},
        "ABCG5": {"mechanism": "efflux", "evidence": "mendelian", "score": 1.0},
        "ABCG8": {"mechanism": "efflux", "evidence": "mendelian", "score": 1.0},
        "APOE": {"mechanism": "ligand", "evidence": "mendelian", "score": 1.0},
        "LIPC": {"mechanism": "lipase", "evidence": "drug", "score": 0.9},
        "LPL": {"mechanism": "lipase", "evidence": "mendelian", "score": 1.0},
        "CETP": {"mechanism": "transfer", "evidence": "drug", "score": 0.9},
        "ANGPTL3": {"mechanism": "regulation", "evidence": "drug", "score": 0.9},
        "ANGPTL4": {"mechanism": "regulation", "evidence": "gwas", "score": 0.8},
        "NPC1L1": {"mechanism": "absorption", "evidence": "drug", "score": 0.95},
        "HMGCR": {"mechanism": "synthesis", "evidence": "drug", "score": 1.0},
        "APOA1": {"mechanism": "structural", "evidence": "mendelian", "score": 1.0},
        "APOC3": {"mechanism": "regulation", "evidence": "drug", "score": 0.9},
        "SORT1": {"mechanism": "trafficking", "evidence": "gwas", "score": 0.85},
    }
    
    # Known T2D genes
    T2D_GENES = {
        "TCF7L2": {"mechanism": "transcription", "evidence": "gwas", "score": 0.95},
        "KCNJ11": {"mechanism": "channel", "evidence": "mendelian", "score": 1.0},
        "ABCC8": {"mechanism": "channel", "evidence": "mendelian", "score": 1.0},
        "HNF1A": {"mechanism": "transcription", "evidence": "mendelian", "score": 1.0},
        "HNF4A": {"mechanism": "transcription", "evidence": "mendelian", "score": 1.0},
        "HNF1B": {"mechanism": "transcription", "evidence": "mendelian", "score": 1.0},
        "GCK": {"mechanism": "metabolism", "evidence": "mendelian", "score": 1.0},
        "INS": {"mechanism": "hormone", "evidence": "mendelian", "score": 1.0},
        "PPARG": {"mechanism": "transcription", "evidence": "drug", "score": 0.95},
        "SLC30A8": {"mechanism": "transport", "evidence": "gwas", "score": 0.85},
        "GCKR": {"mechanism": "regulation", "evidence": "gwas", "score": 0.8},
        "MTNR1B": {"mechanism": "signaling", "evidence": "gwas", "score": 0.8},
    }
    
    # Known CAD genes
    CAD_GENES = {
        "LDLR": {"mechanism": "receptor", "evidence": "mendelian", "score": 1.0},
        "APOB": {"mechanism": "ligand", "evidence": "mendelian", "score": 1.0},
        "PCSK9": {"mechanism": "degradation", "evidence": "drug", "score": 1.0},
        "LPA": {"mechanism": "lipoprotein", "evidence": "mendelian", "score": 1.0},
        "PHACTR1": {"mechanism": "signaling", "evidence": "gwas", "score": 0.8},
        "ADAMTS7": {"mechanism": "matrix", "evidence": "gwas", "score": 0.8},
        "IL6R": {"mechanism": "inflammation", "evidence": "drug", "score": 0.85},
        "GUCY1A1": {"mechanism": "signaling", "evidence": "gwas", "score": 0.8},
        "NOS3": {"mechanism": "signaling", "evidence": "gwas", "score": 0.85},
    }
    
    # Known blood pressure genes
    BP_GENES = {
        "AGT": {"mechanism": "angiotensin", "evidence": "drug", "score": 1.0},
        "ACE": {"mechanism": "enzyme", "evidence": "drug", "score": 1.0},
        "AGTR1": {"mechanism": "receptor", "evidence": "drug", "score": 1.0},
        "NR3C2": {"mechanism": "transcription", "evidence": "drug", "score": 0.95},
        "SCNN1A": {"mechanism": "channel", "evidence": "mendelian", "score": 1.0},
        "SCNN1B": {"mechanism": "channel", "evidence": "mendelian", "score": 1.0},
        "SCNN1G": {"mechanism": "channel", "evidence": "mendelian", "score": 1.0},
        "WNK1": {"mechanism": "kinase", "evidence": "mendelian", "score": 1.0},
        "WNK4": {"mechanism": "kinase", "evidence": "mendelian", "score": 1.0},
        "CYP11B2": {"mechanism": "enzyme", "evidence": "mendelian", "score": 1.0},
        "NPPA": {"mechanism": "hormone", "evidence": "gwas", "score": 0.85},
        "NPPB": {"mechanism": "hormone", "evidence": "gwas", "score": 0.85},
        "GUCY1A1": {"mechanism": "signaling", "evidence": "drug", "score": 0.9},
    }
    
    TRAIT_GENE_SETS = {
        "ldl_cholesterol": LIPID_GENES,
        "hdl_cholesterol": LIPID_GENES,
        "triglycerides": LIPID_GENES,
        "total_cholesterol": LIPID_GENES,
        "type_2_diabetes": T2D_GENES,
        "coronary_artery_disease": CAD_GENES,
        "systolic_blood_pressure": BP_GENES,
        "diastolic_blood_pressure": BP_GENES,
    }
    
    def __init__(
        self,
        trait: str,
        external_file: Optional[str] = None,
    ):
        """
        Initialize gold standard for a trait.
        
        Parameters
        ----------
        trait : str
            Trait name.
        external_file : str, optional
            Path to external gene list.
        """
        self.trait = trait
        
        # Get built-in genes
        self.genes = self.TRAIT_GENE_SETS.get(trait, {}).copy()
        
        # Load from config
        traits_config = get_config("traits")
        config_genes = (
            traits_config.get("trait_tissue_priors", {})
            .get(trait, {})
            .get("benchmark_genes", [])
        )
        
        for gene in config_genes:
            if gene not in self.genes:
                self.genes[gene] = {"evidence": "config", "score": 0.7}
        
        # Load external file if provided
        if external_file and Path(external_file).exists():
            self._load_external(external_file)
        
        logger.info(f"Loaded {len(self.genes)} gold standard genes for {trait}")
    
    def _load_external(self, path: str) -> None:
        """Load genes from external file."""
        df = pd.read_csv(path, sep="\t")
        
        for _, row in df.iterrows():
            gene = row.get("gene", row.get("symbol", ""))
            if gene and gene not in self.genes:
                self.genes[gene] = {
                    "evidence": "external",
                    "score": row.get("score", 0.5),
                }
    
    def is_gold_standard(self, gene: str) -> bool:
        """Check if gene is in gold standard."""
        return gene in self.genes or gene.upper() in self.genes
    
    def get_score(self, gene: str) -> float:
        """Get gold standard score for a gene."""
        if gene in self.genes:
            return self.genes[gene].get("score", 0.5)
        if gene.upper() in self.genes:
            return self.genes[gene.upper()].get("score", 0.5)
        return 0.0
    
    def get_gene_set(self, min_score: float = 0.0) -> Set[str]:
        """Get set of gold standard genes above threshold."""
        return {
            gene for gene, info in self.genes.items()
            if info.get("score", 0) >= min_score
        }


class BenchmarkSuite:
    """
    Comprehensive benchmarking for mechanism graphs.
    """
    
    def __init__(
        self,
        trait: str,
        gold_standard: Optional[GoldStandardGenes] = None,
    ):
        """
        Initialize benchmark suite.
        
        Parameters
        ----------
        trait : str
            Trait name.
        gold_standard : GoldStandardGenes, optional
            Gold standard genes.
        """
        self.trait = trait
        self.gold_standard = gold_standard or GoldStandardGenes(trait)
    
    def evaluate_ranking(
        self,
        ranked_genes: List[Dict],
        k_values: Optional[List[int]] = None,
    ) -> Dict[str, Any]:
        """
        Evaluate gene ranking against gold standard.
        
        Parameters
        ----------
        ranked_genes : list
            List of dicts with 'gene' and 'score' keys, sorted by score.
        k_values : list
            Values of k for recall@k.
            
        Returns
        -------
        dict
            Evaluation metrics.
        """
        if k_values is None:
            k_values = [10, 20, 50, 100]
            
        gold_set = self.gold_standard.get_gene_set(min_score=0.5)
        
        results = {
            "n_ranked": len(ranked_genes),
            "n_gold_standard": len(gold_set),
        }
        
        # Recall@k
        for k in k_values:
            top_k_genes = {
                g.get("gene", g.get("symbol", "")).upper()
                for g in ranked_genes[:k]
            }
            
            recall = len(top_k_genes & gold_set) / max(len(gold_set), 1)
            results[f"recall@{k}"] = recall
        
        # Precision@k
        for k in k_values:
            top_k_genes = {
                g.get("gene", g.get("symbol", "")).upper()
                for g in ranked_genes[:k]
            }
            
            precision = len(top_k_genes & gold_set) / max(k, 1)
            results[f"precision@{k}"] = precision
        
        # Mean rank of gold standard genes
        gold_ranks = []
        for i, g in enumerate(ranked_genes):
            gene = g.get("gene", g.get("symbol", "")).upper()
            if gene in gold_set:
                gold_ranks.append(i + 1)
        
        if gold_ranks:
            results["mean_gold_rank"] = np.mean(gold_ranks)
            results["median_gold_rank"] = np.median(gold_ranks)
        
        # AUROC approximation
        results["auroc"] = self._compute_auroc(ranked_genes, gold_set)
        
        # Average precision
        results["average_precision"] = self._compute_ap(ranked_genes, gold_set)
        
        return results
    
    def _compute_auroc(
        self,
        ranked_genes: List[Dict],
        gold_set: Set[str],
    ) -> float:
        """Compute AUROC."""
        y_true = []
        y_score = []
        
        for g in ranked_genes:
            gene = g.get("gene", g.get("symbol", "")).upper()
            score = g.get("score", g.get("combined_probability", 0))
            
            y_true.append(1 if gene in gold_set else 0)
            y_score.append(score)
        
        if sum(y_true) == 0 or sum(y_true) == len(y_true):
            return 0.5
        
        # Simple AUROC calculation
        n_pos = sum(y_true)
        n_neg = len(y_true) - n_pos
        
        # Sort by score descending
        pairs = sorted(zip(y_score, y_true), reverse=True)
        
        auc = 0.0
        pos_above = 0
        
        for score, label in pairs:
            if label == 0:
                auc += pos_above
            else:
                pos_above += 1
        
        return auc / (n_pos * n_neg) if (n_pos * n_neg) > 0 else 0.5
    
    def _compute_ap(
        self,
        ranked_genes: List[Dict],
        gold_set: Set[str],
    ) -> float:
        """Compute average precision."""
        precisions = []
        n_relevant = 0
        
        for i, g in enumerate(ranked_genes):
            gene = g.get("gene", g.get("symbol", "")).upper()
            
            if gene in gold_set:
                n_relevant += 1
                precision = n_relevant / (i + 1)
                precisions.append(precision)
        
        if not precisions:
            return 0.0
        
        return np.mean(precisions)
    
    def compare_methods(
        self,
        method_rankings: Dict[str, List[Dict]],
    ) -> pd.DataFrame:
        """
        Compare multiple methods.
        
        Parameters
        ----------
        method_rankings : dict
            Method name -> ranked gene list.
            
        Returns
        -------
        pd.DataFrame
            Comparison table.
        """
        results = []
        
        for method, ranking in method_rankings.items():
            metrics = self.evaluate_ranking(ranking)
            metrics["method"] = method
            results.append(metrics)
        
        return pd.DataFrame(results)
    
    def negative_control_test(
        self,
        ranked_genes: List[Dict],
        n_permutations: int = 1000,
    ) -> Dict[str, Any]:
        """
        Test against permuted gold standard.
        
        Parameters
        ----------
        ranked_genes : list
            Ranked gene list.
        n_permutations : int
            Number of permutations.
            
        Returns
        -------
        dict
            Permutation test results.
        """
        all_genes = [
            g.get("gene", g.get("symbol", "")).upper()
            for g in ranked_genes
        ]
        
        gold_set = self.gold_standard.get_gene_set()
        
        # Observed recall@20
        observed = self.evaluate_ranking(ranked_genes)["recall@20"]
        
        # Permuted recalls
        permuted_recalls = []
        
        for _ in range(n_permutations):
            # Random sample of same size as gold standard
            fake_gold = set(np.random.choice(
                all_genes,
                size=len(gold_set),
                replace=False,
            ))
            
            top_20 = {
                g.get("gene", g.get("symbol", "")).upper()
                for g in ranked_genes[:20]
            }
            
            recall = len(top_20 & fake_gold) / max(len(fake_gold), 1)
            permuted_recalls.append(recall)
        
        p_value = (np.array(permuted_recalls) >= observed).mean()
        
        return {
            "observed_recall@20": observed,
            "mean_permuted_recall@20": np.mean(permuted_recalls),
            "std_permuted_recall@20": np.std(permuted_recalls),
            "p_value": p_value,
            "significant": p_value < 0.05,
        }


class ThreeTierBenchmark:
    """
    Three-tier benchmark with anti-leak provisions.
    
    Implements the reviewer-requested benchmark structure:
    
    Tier 1: Holdout Curated (anti-leak)
        - Mendelian genes not in L2G/ABC training
        - ~50-100 genes per trait
        
    Tier 2: Drug Target Validation
        - Approved drugs with genetic support
        - ~100-200 genes per trait
        
    Tier 3: CRISPR Ground Truth
        - CRISPRi validated E-G pairs
        - ~30-50 genes (highest confidence)
    
    Example
    -------
    >>> benchmark = ThreeTierBenchmark("ldl_cholesterol")
    >>> benchmark.load_omim_genes("omim_genes.txt")
    >>> benchmark.load_drug_targets("chembl_targets.txt")
    >>> benchmark.load_crispr_validation("crispr_validated.txt")
    >>> benchmark.flag_training_leaks("l2g_training_genes.txt")
    >>> 
    >>> # Evaluate with anti-leak holdout only
    >>> results = benchmark.evaluate_ranking(ranked_genes, tier="tier1_holdout")
    """
    
    # L2G training genes to flag (from Open Targets)
    # These should NOT be used in Tier 1 evaluation
    L2G_TRAINING_GENES = {
        # Known training set genes from Mountjoy et al. 2021
        "LDLR", "PCSK9", "APOE", "SORT1", "HMGCR",  # Lipids
        "TCF7L2", "KCNJ11", "SLC30A8", "PPARG",     # T2D
        # Add more from published L2G training set
    }
    
    # ABC training genes (from Fulco et al. 2019)
    ABC_TRAINING_GENES = {
        # CRISPRi validation set used to train ABC
        "MYC", "GATA1", "BCL11A", "HBG1", "HBG2",
        # Add more from Fulco supplementary
    }
    
    def __init__(
        self,
        trait: str,
        l2g_training_file: Optional[str] = None,
        abc_training_file: Optional[str] = None,
    ):
        """
        Initialize three-tier benchmark.
        
        Parameters
        ----------
        trait : str
            Trait name.
        l2g_training_file : str, optional
            Path to L2G training genes (for leak flagging).
        abc_training_file : str, optional
            Path to ABC training genes (for leak flagging).
        """
        self.trait = trait
        
        # Gene storage by tier
        self.genes: Dict[BenchmarkTier, List[BenchmarkGene]] = {
            BenchmarkTier.TIER1_HOLDOUT: [],
            BenchmarkTier.TIER2_DRUG: [],
            BenchmarkTier.TIER3_CRISPR: [],
        }
        
        # Load training gene lists for leak detection
        self.l2g_training = self.L2G_TRAINING_GENES.copy()
        self.abc_training = self.ABC_TRAINING_GENES.copy()
        
        if l2g_training_file and Path(l2g_training_file).exists():
            self._load_training_genes(l2g_training_file, "l2g")
        
        if abc_training_file and Path(abc_training_file).exists():
            self._load_training_genes(abc_training_file, "abc")
        
        logger.info(f"Initialized ThreeTierBenchmark for {trait}")
    
    def _load_training_genes(self, filepath: str, source: str) -> None:
        """Load training genes for leak detection."""
        df = pd.read_csv(filepath, sep="\t")
        gene_col = "gene" if "gene" in df.columns else df.columns[0]
        genes = set(df[gene_col].str.upper())
        
        if source == "l2g":
            self.l2g_training.update(genes)
        elif source == "abc":
            self.abc_training.update(genes)
        
        logger.info(f"Loaded {len(genes)} {source} training genes for leak detection")
    
    def add_gene(
        self,
        symbol: str,
        tier: BenchmarkTier,
        evidence_type: str,
        evidence_source: str,
        evidence_id: str = "",
        confidence: float = 1.0,
        ensembl_id: str = "",
        provenance: Optional[Dict] = None,
    ) -> BenchmarkGene:
        """
        Add a benchmark gene.
        
        Parameters
        ----------
        symbol : str
            Gene symbol.
        tier : BenchmarkTier
            Benchmark tier.
        evidence_type : str
            Evidence type (mendelian, drug, crispr).
        evidence_source : str
            Source database.
        evidence_id : str
            Source identifier.
        confidence : float
            Confidence score.
        ensembl_id : str
            Ensembl ID.
        provenance : dict, optional
            Additional provenance.
            
        Returns
        -------
        BenchmarkGene
            Added gene.
        """
        # Check for training leaks
        symbol_upper = symbol.upper()
        is_in_l2g = symbol_upper in self.l2g_training
        is_in_abc = symbol_upper in self.abc_training
        
        gene = BenchmarkGene(
            symbol=symbol,
            ensembl_id=ensembl_id,
            tier=tier,
            trait=self.trait,
            evidence_type=evidence_type,
            evidence_source=evidence_source,
            evidence_id=evidence_id,
            confidence=confidence,
            is_in_l2g_training=is_in_l2g,
            is_in_abc_training=is_in_abc,
            provenance=provenance or {},
        )
        
        self.genes[tier].append(gene)
        
        if gene.has_leak_flag:
            logger.warning(
                f"Gene {symbol} has training leak flag "
                f"(L2G={is_in_l2g}, ABC={is_in_abc})"
            )
        
        return gene
    
    def load_omim_genes(
        self,
        filepath: str,
        confidence_threshold: float = 0.9,
    ) -> int:
        """
        Load Tier 1 genes from OMIM Mendelian disorders.
        
        Parameters
        ----------
        filepath : str
            Path to OMIM gene file.
        confidence_threshold : float
            Minimum confidence for inclusion.
            
        Returns
        -------
        int
            Number of genes loaded.
        """
        df = pd.read_csv(filepath, sep="\t")
        
        count = 0
        for _, row in df.iterrows():
            confidence = float(row.get("confidence", 1.0))
            if confidence < confidence_threshold:
                continue
            
            self.add_gene(
                symbol=str(row.get("gene", row.get("symbol", ""))),
                tier=BenchmarkTier.TIER1_HOLDOUT,
                evidence_type="mendelian",
                evidence_source="OMIM",
                evidence_id=str(row.get("omim_id", "")),
                confidence=confidence,
                ensembl_id=str(row.get("ensembl_id", "")),
                provenance={
                    "disorder": row.get("disorder", ""),
                    "inheritance": row.get("inheritance", ""),
                },
            )
            count += 1
        
        logger.info(f"Loaded {count} OMIM genes for Tier 1")
        return count
    
    def load_drug_targets(
        self,
        filepath: str,
        min_phase: int = 2,
    ) -> int:
        """
        Load Tier 2 genes from drug target databases.
        
        Parameters
        ----------
        filepath : str
            Path to drug target file (ChEMBL/OpenTargets format).
        min_phase : int
            Minimum clinical phase.
            
        Returns
        -------
        int
            Number of genes loaded.
        """
        df = pd.read_csv(filepath, sep="\t")
        
        count = 0
        for _, row in df.iterrows():
            phase = int(row.get("max_phase", row.get("phase", 0)))
            if phase < min_phase:
                continue
            
            # Confidence based on clinical phase
            if phase >= 4:
                confidence = 1.0  # Approved drug
            elif phase >= 3:
                confidence = 0.9
            else:
                confidence = 0.8
            
            self.add_gene(
                symbol=str(row.get("target_gene", row.get("gene", ""))),
                tier=BenchmarkTier.TIER2_DRUG,
                evidence_type="drug_target",
                evidence_source=str(row.get("source", "ChEMBL")),
                evidence_id=str(row.get("chembl_id", row.get("drug_id", ""))),
                confidence=confidence,
                provenance={
                    "drug_name": row.get("drug_name", ""),
                    "indication": row.get("indication", ""),
                    "mechanism": row.get("mechanism", ""),
                    "phase": phase,
                },
            )
            count += 1
        
        logger.info(f"Loaded {count} drug target genes for Tier 2")
        return count
    
    def load_crispr_validation(
        self,
        filepath: str,
        min_effect_size: float = 0.2,
    ) -> int:
        """
        Load Tier 3 genes from CRISPR validation experiments.
        
        Parameters
        ----------
        filepath : str
            Path to CRISPR validation file.
        min_effect_size : float
            Minimum absolute effect size for inclusion.
            
        Returns
        -------
        int
            Number of genes loaded.
        """
        df = pd.read_csv(filepath, sep="\t")
        
        count = 0
        for _, row in df.iterrows():
            effect = abs(float(row.get("effect_size", row.get("effect", 0))))
            if effect < min_effect_size:
                continue
            
            # Confidence based on effect size and significance
            pval = float(row.get("pvalue", row.get("p_value", 1.0)))
            if pval < 0.001 and effect > 0.5:
                confidence = 1.0
            elif pval < 0.01:
                confidence = 0.95
            else:
                confidence = 0.9
            
            self.add_gene(
                symbol=str(row.get("gene", row.get("target_gene", ""))),
                tier=BenchmarkTier.TIER3_CRISPR,
                evidence_type="crispr",
                evidence_source=str(row.get("source", row.get("study", ""))),
                evidence_id=str(row.get("experiment_id", "")),
                confidence=confidence,
                provenance={
                    "enhancer": row.get("enhancer", ""),
                    "effect_size": effect,
                    "pvalue": pval,
                    "cell_type": row.get("cell_type", ""),
                    "method": row.get("method", "CRISPRi"),
                },
            )
            count += 1
        
        logger.info(f"Loaded {count} CRISPR-validated genes for Tier 3")
        return count
    
    def get_gene_set(
        self,
        tier: Optional[BenchmarkTier] = None,
        exclude_leaks: bool = True,
        min_confidence: float = 0.0,
    ) -> Set[str]:
        """
        Get benchmark gene set.
        
        Parameters
        ----------
        tier : BenchmarkTier, optional
            Specific tier (None for all).
        exclude_leaks : bool
            Exclude genes with training leak flags.
        min_confidence : float
            Minimum confidence threshold.
            
        Returns
        -------
        set
            Set of gene symbols.
        """
        genes = set()
        
        tiers = [tier] if tier else list(BenchmarkTier)
        
        for t in tiers:
            for gene in self.genes[t]:
                if gene.confidence < min_confidence:
                    continue
                
                if exclude_leaks and gene.has_leak_flag:
                    continue
                
                genes.add(gene.symbol.upper())
        
        return genes
    
    def get_genes_with_provenance(
        self,
        tier: Optional[BenchmarkTier] = None,
    ) -> List[BenchmarkGene]:
        """
        Get benchmark genes with full provenance.
        
        Parameters
        ----------
        tier : BenchmarkTier, optional
            Specific tier (None for all).
            
        Returns
        -------
        list
            List of BenchmarkGene objects.
        """
        if tier:
            return self.genes[tier].copy()
        
        all_genes = []
        for t in BenchmarkTier:
            all_genes.extend(self.genes[t])
        return all_genes
    
    def evaluate_ranking(
        self,
        ranked_genes: List[Dict],
        tier: Optional[BenchmarkTier] = None,
        exclude_leaks: bool = True,
        k_values: Optional[List[int]] = None,
    ) -> Dict[str, Any]:
        """
        Evaluate gene ranking against benchmark tier.
        
        Parameters
        ----------
        ranked_genes : list
            List of dicts with 'gene' and 'score' keys.
        tier : BenchmarkTier, optional
            Specific tier to evaluate (None for all).
        exclude_leaks : bool
            Exclude genes with training leak flags.
        k_values : list, optional
            Values for recall@k and precision@k.
            
        Returns
        -------
        dict
            Evaluation metrics.
        """
        if k_values is None:
            k_values = [10, 20, 50, 100]
        
        gold_set = self.get_gene_set(
            tier=tier,
            exclude_leaks=exclude_leaks,
        )
        
        results = {
            "tier": tier.value if tier else "all",
            "exclude_leaks": exclude_leaks,
            "n_ranked": len(ranked_genes),
            "n_gold_standard": len(gold_set),
        }
        
        # Recall@k and Precision@k
        for k in k_values:
            top_k = {
                g.get("gene", g.get("symbol", "")).upper()
                for g in ranked_genes[:k]
            }
            
            n_hit = len(top_k & gold_set)
            results[f"recall@{k}"] = n_hit / max(len(gold_set), 1)
            results[f"precision@{k}"] = n_hit / k
        
        # AUROC and Average Precision
        results["auroc"] = self._compute_auroc(ranked_genes, gold_set)
        results["average_precision"] = self._compute_ap(ranked_genes, gold_set)
        
        # Mean/median rank of gold standard genes
        gold_ranks = []
        for i, g in enumerate(ranked_genes):
            gene = g.get("gene", g.get("symbol", "")).upper()
            if gene in gold_set:
                gold_ranks.append(i + 1)
        
        if gold_ranks:
            results["mean_gold_rank"] = np.mean(gold_ranks)
            results["median_gold_rank"] = np.median(gold_ranks)
            results["n_gold_in_ranking"] = len(gold_ranks)
        
        return results
    
    def _compute_auroc(
        self,
        ranked_genes: List[Dict],
        gold_set: Set[str],
    ) -> float:
        """Compute AUROC."""
        if not gold_set:
            return np.nan
        
        y_true = []
        for g in ranked_genes:
            gene = g.get("gene", g.get("symbol", "")).upper()
            y_true.append(1 if gene in gold_set else 0)
        
        n_pos = sum(y_true)
        n_neg = len(y_true) - n_pos
        
        if n_pos == 0 or n_neg == 0:
            return 0.5
        
        auc = 0.0
        pos_above = 0
        
        for label in y_true:
            if label == 0:
                auc += pos_above
            else:
                pos_above += 1
        
        return auc / (n_pos * n_neg)
    
    def _compute_ap(
        self,
        ranked_genes: List[Dict],
        gold_set: Set[str],
    ) -> float:
        """Compute average precision."""
        if not gold_set:
            return np.nan
        
        precisions = []
        n_relevant = 0
        
        for i, g in enumerate(ranked_genes):
            gene = g.get("gene", g.get("symbol", "")).upper()
            
            if gene in gold_set:
                n_relevant += 1
                precisions.append(n_relevant / (i + 1))
        
        return np.mean(precisions) if precisions else 0.0
    
    def tier_comparison(
        self,
        ranked_genes: List[Dict],
        exclude_leaks: bool = True,
    ) -> pd.DataFrame:
        """
        Compare performance across all tiers.
        
        Parameters
        ----------
        ranked_genes : list
            Ranked gene list.
        exclude_leaks : bool
            Exclude leaked genes.
            
        Returns
        -------
        pd.DataFrame
            Comparison table.
        """
        results = []
        
        for tier in BenchmarkTier:
            metrics = self.evaluate_ranking(
                ranked_genes,
                tier=tier,
                exclude_leaks=exclude_leaks,
            )
            metrics["tier_name"] = tier.name
            results.append(metrics)
        
        # Add combined (all tiers)
        all_metrics = self.evaluate_ranking(
            ranked_genes,
            tier=None,
            exclude_leaks=exclude_leaks,
        )
        all_metrics["tier_name"] = "ALL_TIERS"
        results.append(all_metrics)
        
        return pd.DataFrame(results)
    
    def leak_analysis(self) -> Dict[str, Any]:
        """
        Analyze potential training data leaks.
        
        Returns
        -------
        dict
            Leak analysis summary.
        """
        leak_summary = {
            "total_genes": 0,
            "genes_with_leaks": 0,
            "l2g_leaks": 0,
            "abc_leaks": 0,
            "by_tier": {},
        }
        
        for tier in BenchmarkTier:
            tier_total = len(self.genes[tier])
            tier_leaks = sum(1 for g in self.genes[tier] if g.has_leak_flag)
            
            leak_summary["by_tier"][tier.value] = {
                "total": tier_total,
                "leaked": tier_leaks,
                "clean": tier_total - tier_leaks,
                "leak_rate": tier_leaks / max(tier_total, 1),
            }
            
            leak_summary["total_genes"] += tier_total
            leak_summary["genes_with_leaks"] += tier_leaks
            leak_summary["l2g_leaks"] += sum(
                1 for g in self.genes[tier] if g.is_in_l2g_training
            )
            leak_summary["abc_leaks"] += sum(
                1 for g in self.genes[tier] if g.is_in_abc_training
            )
        
        leak_summary["overall_leak_rate"] = (
            leak_summary["genes_with_leaks"] / 
            max(leak_summary["total_genes"], 1)
        )
        
        return leak_summary
    
    def summary(self) -> Dict[str, Any]:
        """Get benchmark summary."""
        return {
            "trait": self.trait,
            "tier1_genes": len(self.genes[BenchmarkTier.TIER1_HOLDOUT]),
            "tier2_genes": len(self.genes[BenchmarkTier.TIER2_DRUG]),
            "tier3_genes": len(self.genes[BenchmarkTier.TIER3_CRISPR]),
            "total_genes": sum(len(g) for g in self.genes.values()),
            "leak_analysis": self.leak_analysis(),
            "l2g_training_size": len(self.l2g_training),
            "abc_training_size": len(self.abc_training),
        }
    
    def export_provenance(
        self,
        filepath: str,
        tier: Optional[BenchmarkTier] = None,
    ) -> None:
        """
        Export benchmark genes with full provenance.
        
        Parameters
        ----------
        filepath : str
            Output file path.
        tier : BenchmarkTier, optional
            Specific tier (None for all).
        """
        genes = self.get_genes_with_provenance(tier)
        
        rows = [g.to_dict() for g in genes]
        df = pd.DataFrame(rows)
        
        df.to_csv(filepath, sep="\t", index=False)
        logger.info(f"Exported {len(rows)} benchmark genes to {filepath}")
