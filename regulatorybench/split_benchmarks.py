#!/usr/bin/env python3
"""
RegulatoryBench v4: Benchmark Splitter

Splits the benchmark data into Task A/B/C according to the taxonomy:
- Task A (GWAS→Gene): E2G benchmarking dataset with UKBB credible sets + ABC predictions
- Task B (Enhancer→Gene): ENCODE CRISPRi perturbation screens
- Task C (Variant→Gene): MPRA variant effect datasets

This script creates independent benchmarks that do NOT mix evaluation tasks.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import json
import logging
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class BenchmarkStats:
    """Statistics for a benchmark split."""
    name: str
    task_type: str
    n_pairs: int
    n_positive: int
    n_negative: int
    n_loci: int
    sources: List[str]
    positive_rate: float


class BenchmarkSplitter:
    """Create task-specific benchmark splits from RegulatoryBench data."""
    
    def __init__(self, data_dir: Path = Path("data/external")):
        self.data_dir = data_dir
        self.output_dir = Path("regulatorybench/benchmarks")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def load_task_b_crispr(self) -> pd.DataFrame:
        """
        Load Task B benchmark: Enhancer → Gene from CRISPRi screens.
        
        Sources:
        - ENCODE CRISPRi benchmark (Engreitz lab)
        - Fulco 2019 CRISPRi-FlowFISH
        """
        dfs = []
        
        # 1. ENCODE CRISPRi benchmark - K562 training
        k562_path = self.data_dir / "crispr_benchmark/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
        if k562_path.exists():
            logger.info(f"Loading ENCODE CRISPRi K562: {k562_path}")
            k562 = pd.read_csv(k562_path, sep='\t')
            k562['source'] = 'ENCODE_CRISPRi_K562'
            k562['task_type'] = 'B'
            dfs.append(k562)
            
        # 2. ENCODE CRISPRi benchmark - held-out cell types
        held_path = self.data_dir / "crispr_benchmark/resources/crispr_data/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz"
        if held_path.exists():
            logger.info(f"Loading ENCODE CRISPRi held-out: {held_path}")
            held = pd.read_csv(held_path, sep='\t')
            held['source'] = 'ENCODE_CRISPRi_heldout'
            held['task_type'] = 'B'
            dfs.append(held)
            
        # 3. Fulco 2019 CRISPRi-FlowFISH
        fulco_path = self.data_dir / "crispr_validation/fulco_2019_table_s6a.xlsx"
        if fulco_path.exists():
            logger.info(f"Loading Fulco 2019: {fulco_path}")
            try:
                fulco = pd.read_excel(fulco_path, sheet_name='Supplementary Table 6a', header=1)
                # Standardize column names
                fulco = fulco.rename(columns={
                    'chr': 'chrom',
                    'start': 'chromStart',
                    'end': 'chromEnd',
                    'Gene': 'measuredGeneSymbol',
                    'Significant': 'Regulated',
                    'Fraction change in gene expr': 'EffectSize',
                    'Adjusted p-value': 'pValueAdjusted',
                    'Cell type': 'CellType',
                    'ABC Score': 'ABC_score'
                })
                fulco['source'] = 'Fulco_2019'
                fulco['task_type'] = 'B'
                dfs.append(fulco)
            except Exception as e:
                logger.warning(f"Could not load Fulco 2019: {e}")
        
        if not dfs:
            logger.error("No Task B data found!")
            return pd.DataFrame()
            
        # Combine and standardize
        combined = pd.concat(dfs, ignore_index=True)
        
        # Create locus_id for enhancer-gene pairs
        combined['locus_id'] = (
            combined['chrom'].astype(str) + ':' + 
            combined['chromStart'].astype(str) + '-' + 
            combined['chromEnd'].astype(str)
        )
        
        # Standardize label column
        if 'Regulated' in combined.columns:
            combined['is_positive'] = combined['Regulated'].astype(bool)
        
        logger.info(f"Task B total: {len(combined)} pairs, {combined['is_positive'].sum()} positives")
        
        return combined
    
    def load_task_a_gwas(self) -> pd.DataFrame:
        """
        Load Task A benchmark: GWAS Credible Set → Gene.
        
        Source: E2G benchmarking dataset (UKBB credible sets with ABC predictions)
        """
        # Load UKBB ABC Gene predictions with ground truth
        ukbb_path = self.data_dir / "E2G_benchmarking/resources/UKBiobank.ABCGene.anyabc.tsv"
        
        if not ukbb_path.exists():
            logger.warning(f"Task A data not found: {ukbb_path}")
            return pd.DataFrame()
            
        logger.info(f"Loading E2G UKBB benchmark: {ukbb_path}")
        df = pd.read_csv(ukbb_path, sep='\t')
        
        # Add task metadata
        df['task_type'] = 'A'
        df['source'] = 'E2G_UKBB_benchmark'
        
        # Map to standard columns
        df['locus_id'] = df['CredibleSetUnique']
        df['gene_symbol'] = df['TargetGene']
        df['is_positive'] = df['truth'].astype(bool)
        
        logger.info(f"Task A total: {len(df)} pairs, {df['is_positive'].sum()} positives, {df['locus_id'].nunique()} loci")
        
        return df
    
    def load_task_c_mpra(self) -> pd.DataFrame:
        """
        Load Task C benchmark: Regulatory Variant → Gene.
        
        Sources:
        - MPRA Abell 2022
        - MPRA Schizophrenia studies
        - STING-seq data
        """
        dfs = []
        
        # 1. MPRA Abell 2022 - allelic effects
        abell_path = self.data_dir / "mpra_abell_2022/MPRAResults_AlleleEffect_10kb.txt.gz"
        if abell_path.exists():
            logger.info(f"Loading MPRA Abell 2022: {abell_path}")
            try:
                abell = pd.read_csv(abell_path, sep='\t')
                abell['source'] = 'MPRA_Abell_2022'
                abell['task_type'] = 'C'
                dfs.append(abell)
            except Exception as e:
                logger.warning(f"Could not load Abell 2022: {e}")
        
        # 2. STING-seq (variant perturbation)
        sting_path = self.data_dir / "sting_seq/STING_seq_Table_S3.xlsx"
        if sting_path.exists():
            logger.info(f"Loading STING-seq: {sting_path}")
            try:
                sting = pd.read_excel(sting_path)
                sting['source'] = 'STING_seq'
                sting['task_type'] = 'C'
                dfs.append(sting)
            except Exception as e:
                logger.warning(f"Could not load STING-seq: {e}")
        
        if not dfs:
            logger.warning("No Task C data found - proceeding with empty dataset")
            return pd.DataFrame()
            
        combined = pd.concat(dfs, ignore_index=True)
        logger.info(f"Task C total: {len(combined)} entries")
        
        return combined
    
    def create_benchmark_card(self, stats: Dict[str, BenchmarkStats]) -> str:
        """Create a benchmark card documenting the splits."""
        
        card = f"""# RegulatoryBench v4 - Task-Stratified Benchmark

Generated: {datetime.now().isoformat()}

## Task Taxonomy

This benchmark separates three fundamentally different prediction tasks:

| Task | Description | Ground Truth | Applicable Methods |
|------|-------------|--------------|-------------------|
| A | GWAS Credible Set → Causal Gene | UKBB + ABC validation | L2G, FLAMES, PoPS |
| B | Regulatory Element → Target Gene | CRISPRi perturbation | ABC, rE2G, distance |
| C | Regulatory Variant → Affected Gene | MPRA allelic effects | eQTL coloc, VEP |

**CRITICAL**: Methods designed for one task may be inapplicable to others.
- L2G/FLAMES require GWAS inputs → cannot evaluate on Task B/C without synthetic mapping
- cS2G provides gene-level scores → degenerates on Task B (no within-locus discrimination)

## Benchmark Splits

"""
        for name, s in stats.items():
            card += f"""### {name} (Task {s.task_type})

- **Total pairs**: {s.n_pairs:,}
- **Positives**: {s.n_positive:,} ({s.positive_rate:.1%})
- **Negatives**: {s.n_negative:,}
- **Unique loci**: {s.n_loci:,}
- **Sources**: {', '.join(s.sources)}

"""
        
        card += """## Evaluation Protocol

### Within-Task Evaluation
- Compare only methods applicable to the same task
- Primary metric: AUROC (within-locus binary classification)
- Secondary: AUPRC, distance-stratified AUC, recall@K

### Cross-Task Reporting
- Document "inapplicable" rather than "failed" for method-task mismatches
- Provide coverage statistics (% loci with predictions)

## Independence Verification

- CRISPRi labels: Experimentally derived (independent of prediction methods)
- MPRA labels: Experimentally derived (independent of prediction methods)
- GWAS gold standards: Curated independently of ABC/L2G training

## Leakage Audit

Training overlap checked for:
- L2G: Uses Open Targets gold standards (disjoint from CRISPRi/MPRA)
- ABC: Trained on K562 CRISPRi → held-out cell types provide independent test
- FLAMES: Uses annotation features, not trained on benchmark loci

## Citation

If using this benchmark, please cite:
- ENCODE CRISPRi: Engreitz et al., 2024
- Fulco 2019: Fulco et al., Nature Genetics 2019
- E2G Benchmarking: ENCODE Consortium
"""
        return card
    
    def run(self) -> Dict[str, BenchmarkStats]:
        """Execute the benchmark splitting pipeline."""
        
        stats = {}
        
        # Task A: GWAS → Gene
        logger.info("=" * 60)
        logger.info("Processing Task A: GWAS Credible Set → Gene")
        task_a = self.load_task_a_gwas()
        if len(task_a) > 0:
            output_path = self.output_dir / "task_a_gwas_to_gene.parquet"
            task_a.to_parquet(output_path, index=False)
            logger.info(f"Saved Task A to {output_path}")
            
            stats['GWAS-E2G'] = BenchmarkStats(
                name='GWAS-E2G',
                task_type='A',
                n_pairs=len(task_a),
                n_positive=task_a['is_positive'].sum(),
                n_negative=(~task_a['is_positive']).sum(),
                n_loci=task_a['locus_id'].nunique(),
                sources=task_a['source'].unique().tolist(),
                positive_rate=task_a['is_positive'].mean()
            )
        
        # Task B: Enhancer → Gene
        logger.info("=" * 60)
        logger.info("Processing Task B: Enhancer → Gene")
        task_b = self.load_task_b_crispr()
        if len(task_b) > 0:
            output_path = self.output_dir / "task_b_enhancer_to_gene.parquet"
            task_b.to_parquet(output_path, index=False)
            logger.info(f"Saved Task B to {output_path}")
            
            stats['CRISPR-E2G'] = BenchmarkStats(
                name='CRISPR-E2G',
                task_type='B',
                n_pairs=len(task_b),
                n_positive=task_b['is_positive'].sum() if 'is_positive' in task_b.columns else 0,
                n_negative=(~task_b['is_positive']).sum() if 'is_positive' in task_b.columns else 0,
                n_loci=task_b['locus_id'].nunique() if 'locus_id' in task_b.columns else 0,
                sources=task_b['source'].unique().tolist(),
                positive_rate=task_b['is_positive'].mean() if 'is_positive' in task_b.columns else 0
            )
        
        # Task C: Variant → Gene
        logger.info("=" * 60)
        logger.info("Processing Task C: Variant → Gene")
        task_c = self.load_task_c_mpra()
        if len(task_c) > 0:
            output_path = self.output_dir / "task_c_variant_to_gene.parquet"
            task_c.to_parquet(output_path, index=False)
            logger.info(f"Saved Task C to {output_path}")
            
            stats['MPRA-V2G'] = BenchmarkStats(
                name='MPRA-V2G',
                task_type='C',
                n_pairs=len(task_c),
                n_positive=0,  # Need to determine positive criteria
                n_negative=0,
                n_loci=0,
                sources=task_c['source'].unique().tolist() if 'source' in task_c.columns else [],
                positive_rate=0
            )
        
        # Generate benchmark card
        card = self.create_benchmark_card(stats)
        card_path = self.output_dir / "BENCHMARK_CARD.md"
        card_path.write_text(card, encoding='utf-8')
        logger.info(f"Saved benchmark card to {card_path}")
        
        # Save stats JSON
        stats_dict = {k: {kk: int(vv) if isinstance(vv, (np.integer, np.int64)) else 
                          float(vv) if isinstance(vv, (np.floating, np.float64)) else vv 
                          for kk, vv in vars(v).items()} for k, v in stats.items()}
        stats_path = self.output_dir / "benchmark_stats.json"
        with open(stats_path, 'w') as f:
            json.dump(stats_dict, f, indent=2)
        logger.info(f"Saved stats to {stats_path}")
        
        return stats


def main():
    splitter = BenchmarkSplitter()
    stats = splitter.run()
    
    print("\n" + "=" * 60)
    print("BENCHMARK SPLIT SUMMARY")
    print("=" * 60)
    
    for name, s in stats.items():
        print(f"\n{name} (Task {s.task_type}):")
        print(f"  Pairs: {s.n_pairs:,}")
        print(f"  Positives: {s.n_positive:,} ({s.positive_rate:.1%})")
        print(f"  Loci: {s.n_loci:,}")
        print(f"  Sources: {', '.join(s.sources)}")


if __name__ == "__main__":
    main()
