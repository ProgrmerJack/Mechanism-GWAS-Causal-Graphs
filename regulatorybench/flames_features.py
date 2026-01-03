#!/usr/bin/env python3
"""
FLAMES Annotation Features for Task B Evaluation

The FLAMES method (Marijn-Schipper/FLAMES) is designed for Task A (GWAS→Gene)
and requires GWAS-derived inputs (credible sets, PoPS, MAGMA).

However, FLAMES uses annotation features that can be directly applied to
Task B (Enhancer→Gene) evaluation:

1. ABC_CRISPR: Activity-by-Contact scores for enhancer-gene pairs
2. EpiMap: Enhancer-gene correlations from epigenomic data  
3. GeneHancer: Curated enhancer-gene links
4. PCHiC: Promoter capture Hi-C contacts

This script extracts these features and evaluates them on the Task B benchmark.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
from sklearn.metrics import roc_auc_score, average_precision_score
from dataclasses import dataclass

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class MethodResult:
    """Results for a single method evaluation."""
    method_name: str
    task_type: str
    auroc: float
    auprc: float
    coverage: float
    n_pairs_scored: int
    n_pairs_total: int


class FLAMESFeatureExtractor:
    """Extract and evaluate FLAMES annotation features on Task B."""
    
    def __init__(self, flames_dir: Path = Path("data/external/flames/Annotation_data")):
        self.flames_dir = flames_dir
        self.features_available = self._check_available_features()
        
    def _check_available_features(self) -> Dict[str, bool]:
        """Check which FLAMES annotation features are available."""
        features = {
            'ABC_CRISPR': (self.flames_dir / 'ABC_CRISPR').exists(),
            'ABC_EP': (self.flames_dir / 'ABC_EP').exists(),
            'EpiMap': (self.flames_dir / 'EpiMap').exists(),
            'GeneHancer': (self.flames_dir / 'GeneHancer').exists(),
            'RoadMapEpi': (self.flames_dir / 'RoadMapEpi').exists(),
            'Jung_PCHiC': (self.flames_dir / 'Jung_PCHiC').exists(),
            'Javierre_PCHiC': (self.flames_dir / 'Javierre_PCHiC').exists(),
            'HACER': (self.flames_dir / 'HACER').exists(),
            'Cicero': (self.flames_dir / 'Cicero_whole_blood').exists(),
        }
        available = [k for k, v in features.items() if v]
        logger.info(f"Available FLAMES features: {available}")
        return features
    
    def load_abc_crispr_for_gene(self, ensembl_id: str) -> Optional[pd.DataFrame]:
        """Load ABC_CRISPR scores for a specific gene."""
        gene_file = self.flames_dir / 'ABC_CRISPR' / f'{ensembl_id}.parquet'
        if gene_file.exists():
            return pd.read_parquet(gene_file)
        return None
    
    def load_abc_ep_for_gene(self, ensembl_id: str) -> Optional[pd.DataFrame]:
        """Load ABC_EP scores for a specific gene."""
        gene_file = self.flames_dir / 'ABC_EP' / f'{ensembl_id}.parquet'
        if gene_file.exists():
            return pd.read_parquet(gene_file)
        return None
    
    def get_abc_score_for_pair(self, chrom: str, start: int, end: int, 
                               ensembl_id: str, source: str = 'ABC_CRISPR') -> Optional[float]:
        """
        Get ABC score for an enhancer-gene pair.
        
        The enhancer is defined by genomic coordinates, and we look up
        the maximum ABC score for any overlapping element.
        """
        if source == 'ABC_CRISPR':
            df = self.load_abc_crispr_for_gene(ensembl_id)
        elif source == 'ABC_EP':
            df = self.load_abc_ep_for_gene(ensembl_id)
        else:
            return None
            
        if df is None or len(df) == 0:
            return None
            
        # Find overlapping enhancers
        # Normalize chromosome format
        chrom_normalized = chrom.replace('chr', '')
        df_chrom = df[df['chr'].astype(str).str.replace('chr', '') == chrom_normalized]
        
        if len(df_chrom) == 0:
            return None
            
        # Find overlapping elements (allowing some slack)
        slack = 500  # 500bp overlap buffer
        overlapping = df_chrom[
            (df_chrom['start'] <= end + slack) & 
            (df_chrom['end'] >= start - slack)
        ]
        
        if len(overlapping) == 0:
            return None
            
        # Return max ABC score
        score_col = 'ABC_Score' if 'ABC_Score' in overlapping.columns else 'ABC_max'
        if score_col in overlapping.columns:
            return overlapping[score_col].max()
        return None
    
    def score_task_b_benchmark(self, benchmark_path: Path) -> Dict[str, MethodResult]:
        """
        Score the Task B benchmark using FLAMES annotation features.
        
        Returns results for each available feature source.
        """
        # Load benchmark
        logger.info(f"Loading Task B benchmark: {benchmark_path}")
        bench = pd.read_parquet(benchmark_path)
        
        results = {}
        
        # Extract gene IDs from benchmark
        # Need to map gene symbols to Ensembl IDs
        if 'measuredGeneEnsemblId' in bench.columns:
            bench['ensembl_id'] = bench['measuredGeneEnsemblId']
        elif 'ensembl_id' not in bench.columns:
            logger.warning("No Ensembl ID column found - attempting gene symbol lookup")
            # Would need to add gene symbol mapping here
            bench['ensembl_id'] = None
        
        # Score using ABC_CRISPR
        if self.features_available.get('ABC_CRISPR', False):
            logger.info("Scoring with FLAMES ABC_CRISPR features...")
            results['FLAMES_ABC'] = self._score_with_abc(bench, 'ABC_CRISPR')
            
        # Score using ABC_EP  
        if self.features_available.get('ABC_EP', False):
            logger.info("Scoring with FLAMES ABC_EP features...")
            results['FLAMES_ABC_EP'] = self._score_with_abc(bench, 'ABC_EP')
            
        return results
    
    def _score_with_abc(self, bench: pd.DataFrame, source: str) -> MethodResult:
        """Score benchmark pairs using ABC features."""
        
        scores = []
        labels = []
        
        for idx, row in bench.iterrows():
            if idx % 1000 == 0:
                logger.info(f"  Processing {idx}/{len(bench)}...")
                
            ensembl_id = row.get('ensembl_id')
            if pd.isna(ensembl_id):
                continue
                
            score = self.get_abc_score_for_pair(
                chrom=str(row.get('chrom', row.get('chr', ''))),
                start=int(row.get('chromStart', row.get('start', 0))),
                end=int(row.get('chromEnd', row.get('end', 0))),
                ensembl_id=ensembl_id,
                source=source
            )
            
            if score is not None:
                scores.append(score)
                labels.append(1 if row.get('is_positive', row.get('Regulated', False)) else 0)
        
        n_scored = len(scores)
        n_total = len(bench)
        coverage = n_scored / n_total if n_total > 0 else 0
        
        if n_scored < 10 or len(set(labels)) < 2:
            logger.warning(f"Insufficient data for {source}: {n_scored} scored, {len(set(labels))} unique labels")
            return MethodResult(
                method_name=f'FLAMES_{source}',
                task_type='B',
                auroc=np.nan,
                auprc=np.nan,
                coverage=coverage,
                n_pairs_scored=n_scored,
                n_pairs_total=n_total
            )
        
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
        
        logger.info(f"{source}: AUROC={auroc:.3f}, AUPRC={auprc:.3f}, coverage={coverage:.1%}")
        
        return MethodResult(
            method_name=f'FLAMES_{source}',
            task_type='B',
            auroc=auroc,
            auprc=auprc,
            coverage=coverage,
            n_pairs_scored=n_scored,
            n_pairs_total=n_total
        )


def main():
    """Run FLAMES feature evaluation on Task B benchmark."""
    
    extractor = FLAMESFeatureExtractor()
    
    benchmark_path = Path("regulatorybench/benchmarks/task_b_enhancer_to_gene.parquet")
    if not benchmark_path.exists():
        logger.error(f"Benchmark not found: {benchmark_path}")
        return
    
    results = extractor.score_task_b_benchmark(benchmark_path)
    
    print("\n" + "=" * 60)
    print("FLAMES FEATURE EVALUATION ON TASK B")
    print("=" * 60)
    
    for name, result in results.items():
        print(f"\n{result.method_name}:")
        print(f"  AUROC: {result.auroc:.3f}" if not np.isnan(result.auroc) else "  AUROC: N/A")
        print(f"  AUPRC: {result.auprc:.3f}" if not np.isnan(result.auprc) else "  AUPRC: N/A")
        print(f"  Coverage: {result.coverage:.1%}")
        print(f"  Pairs scored: {result.n_pairs_scored:,} / {result.n_pairs_total:,}")


if __name__ == "__main__":
    main()
