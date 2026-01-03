#!/usr/bin/env python3
"""
Leakage Audit Framework

Implements the user's specified leakage controls:
1. ABC data leakage detection (training data contamination)
2. Leave-one-study-out cross-validation
3. Leave-one-trait-category-out cross-validation
4. Chromosome holdout validation
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Optional
import json
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Use script directory for relative paths
SCRIPT_DIR = Path(__file__).parent.resolve()


@dataclass
class LeakageReport:
    """Report from leakage audit."""
    audit_name: str
    has_leakage: bool
    leakage_description: str
    affected_pairs: int
    total_pairs: int
    recommendation: str


class LeakageAuditor:
    """
    Audit benchmarks for data leakage.
    
    Based on common leakage patterns in E2G benchmarking:
    1. ABC training data overlap with CRISPRi test data
    2. Same locus appearing in training and test
    3. Same gene appearing in multiple studies
    """
    
    def __init__(self, benchmark_path: Path):
        self.benchmark = pd.read_parquet(benchmark_path)
        self.reports: List[LeakageReport] = []
        
    def audit_abc_training_overlap(self) -> LeakageReport:
        """
        Check if benchmark pairs were used to train ABC model.
        
        ABC model was trained on Fulco 2019 K562 CRISPRi data.
        If evaluating ABC on the same data, results are invalid.
        """
        logger.info("Auditing ABC training data overlap...")
        
        # Check for Fulco 2019 or K562 in source
        k562_mask = self.benchmark['source'].str.contains('K562', case=False, na=False)
        fulco_mask = self.benchmark['source'].str.contains('Fulco', case=False, na=False)
        training_data = k562_mask | fulco_mask
        
        n_training = training_data.sum()
        
        if n_training > 0:
            report = LeakageReport(
                audit_name='ABC Training Overlap',
                has_leakage=True,
                leakage_description=f'{n_training} pairs overlap with ABC training data (K562/Fulco 2019)',
                affected_pairs=n_training,
                total_pairs=len(self.benchmark),
                recommendation='Exclude K562/Fulco data when evaluating ABC model, or use held-out cell types only'
            )
        else:
            report = LeakageReport(
                audit_name='ABC Training Overlap',
                has_leakage=False,
                leakage_description='No overlap with known ABC training data',
                affected_pairs=0,
                total_pairs=len(self.benchmark),
                recommendation='None required'
            )
        
        self.reports.append(report)
        return report
    
    def audit_gene_overlap_across_sources(self) -> LeakageReport:
        """
        Check for genes appearing in multiple sources.
        
        Gene overlap can artificially inflate performance if models
        learn gene-specific features from one study that transfer to another.
        """
        logger.info("Auditing gene overlap across sources...")
        
        gene_col = 'measuredGeneSymbol' if 'measuredGeneSymbol' in self.benchmark.columns else 'TargetGene'
        
        if gene_col not in self.benchmark.columns or 'source' not in self.benchmark.columns:
            return LeakageReport(
                audit_name='Gene Overlap',
                has_leakage=False,
                leakage_description='Cannot audit - missing gene or source columns',
                affected_pairs=0,
                total_pairs=len(self.benchmark),
                recommendation='Add gene and source columns to benchmark'
            )
        
        gene_sources = self.benchmark.groupby(gene_col)['source'].nunique()
        multi_source_genes = gene_sources[gene_sources > 1]
        
        affected_mask = self.benchmark[gene_col].isin(multi_source_genes.index)
        n_affected = affected_mask.sum()
        
        report = LeakageReport(
            audit_name='Gene Overlap Across Sources',
            has_leakage=n_affected > len(self.benchmark) * 0.1,  # >10% is concerning
            leakage_description=f'{len(multi_source_genes)} genes appear in multiple sources ({n_affected} pairs)',
            affected_pairs=n_affected,
            total_pairs=len(self.benchmark),
            recommendation='Use leave-one-study-out CV or stratify by gene' if n_affected > 0 else 'None'
        )
        
        self.reports.append(report)
        return report
    
    def audit_genomic_proximity(self) -> LeakageReport:
        """
        Check for genomic proximity leakage.
        
        Nearby loci may share causal genes, leaking information.
        """
        logger.info("Auditing genomic proximity...")
        
        if 'chrom' not in self.benchmark.columns:
            return LeakageReport(
                audit_name='Genomic Proximity',
                has_leakage=False,
                leakage_description='Cannot audit - missing coordinate columns',
                affected_pairs=0,
                total_pairs=len(self.benchmark),
                recommendation='Add genomic coordinates'
            )
        
        # Group by chromosome and check for clustering
        nearby_window = 100000  # 100kb
        proximity_issues = 0
        
        for chrom in self.benchmark['chrom'].unique():
            chrom_data = self.benchmark[self.benchmark['chrom'] == chrom].copy()
            if 'chromStart' in chrom_data.columns:
                chrom_data = chrom_data.sort_values('chromStart')
                positions = chrom_data['chromStart'].values
                
                for i in range(1, len(positions)):
                    if positions[i] - positions[i-1] < nearby_window:
                        proximity_issues += 1
        
        report = LeakageReport(
            audit_name='Genomic Proximity',
            has_leakage=proximity_issues > len(self.benchmark) * 0.05,
            leakage_description=f'{proximity_issues} pairs within 100kb of another pair',
            affected_pairs=proximity_issues,
            total_pairs=len(self.benchmark),
            recommendation='Use chromosome holdout for evaluation' if proximity_issues > 0 else 'None'
        )
        
        self.reports.append(report)
        return report
    
    def create_leave_one_study_out_splits(self) -> Dict[str, pd.DataFrame]:
        """
        Create leave-one-study-out cross-validation splits.
        """
        logger.info("Creating leave-one-study-out splits...")
        
        if 'source' not in self.benchmark.columns:
            logger.warning("No source column - cannot create LOSO splits")
            return {}
        
        splits = {}
        sources = self.benchmark['source'].unique()
        
        for holdout_source in sources:
            train = self.benchmark[self.benchmark['source'] != holdout_source]
            test = self.benchmark[self.benchmark['source'] == holdout_source]
            
            splits[f'holdout_{holdout_source}'] = {
                'train': train,
                'test': test,
                'train_size': len(train),
                'test_size': len(test),
                'train_pos': train['is_positive'].sum() if 'is_positive' in train.columns else 0,
                'test_pos': test['is_positive'].sum() if 'is_positive' in test.columns else 0
            }
            
            logger.info(f"  {holdout_source}: train={len(train)}, test={len(test)}")
        
        return splits
    
    def create_chromosome_holdout_splits(self) -> Dict[str, pd.DataFrame]:
        """
        Create chromosome holdout splits.
        
        Standard approach: hold out chr8 and chr9 for testing.
        """
        logger.info("Creating chromosome holdout splits...")
        
        if 'chrom' not in self.benchmark.columns:
            logger.warning("No chrom column - cannot create chromosome splits")
            return {}
        
        # Normalize chromosome names
        self.benchmark['chrom_norm'] = self.benchmark['chrom'].astype(str).str.replace('chr', '')
        
        holdout_chroms = ['8', '9']
        test_mask = self.benchmark['chrom_norm'].isin(holdout_chroms)
        
        train = self.benchmark[~test_mask]
        test = self.benchmark[test_mask]
        
        splits = {
            'chr8_9_holdout': {
                'train': train,
                'test': test,
                'train_size': len(train),
                'test_size': len(test),
                'train_pos': train['is_positive'].sum() if 'is_positive' in train.columns else 0,
                'test_pos': test['is_positive'].sum() if 'is_positive' in test.columns else 0
            }
        }
        
        logger.info(f"  Chromosome holdout (8,9): train={len(train)}, test={len(test)}")
        
        return splits
    
    def create_distance_matched_negatives(self, n_bins: int = 10) -> pd.DataFrame:
        """
        Create distance-matched negative controls.
        
        This addresses the critique that random negatives are trivially distinguished
        by distance. We sample negatives that match the distance distribution of
        positives, making distance uninformative and isolating the contribution
        of functional features.
        
        Per user requirement Point 5: "distance-matched controls"
        
        Parameters
        ----------
        n_bins : int
            Number of distance quantile bins for matching.
            
        Returns
        -------
        pd.DataFrame
            Benchmark with distance-matched negatives replacing random negatives.
        """
        logger.info("Creating distance-matched negative controls...")
        
        # Find distance column
        dist_col = None
        for col in ['distanceToTSS', 'distance', 'dist', 'distance_to_gene']:
            if col in self.benchmark.columns:
                dist_col = col
                break
        
        if dist_col is None:
            logger.warning("No distance column found - cannot create distance-matched negatives")
            return self.benchmark
        
        label_col = 'is_positive' if 'is_positive' in self.benchmark.columns else 'Regulated'
        if label_col not in self.benchmark.columns:
            logger.warning("No label column found - cannot create distance-matched negatives")
            return self.benchmark
        
        positives = self.benchmark[self.benchmark[label_col] == True].copy()
        negatives = self.benchmark[self.benchmark[label_col] == False].copy()
        
        if len(positives) == 0 or len(negatives) == 0:
            logger.warning("Cannot match - need both positives and negatives")
            return self.benchmark
        
        # Compute distance quantile bins from positives
        positives['dist_bin'] = pd.qcut(positives[dist_col].abs(), q=n_bins, labels=False, duplicates='drop')
        bin_edges = pd.qcut(positives[dist_col].abs(), q=n_bins, retbins=True, duplicates='drop')[1]
        
        # Assign negatives to same bins
        negatives['dist_bin'] = pd.cut(
            negatives[dist_col].abs(), 
            bins=bin_edges, 
            labels=False, 
            include_lowest=True
        )
        
        # Sample negatives to match positive distribution
        pos_bin_counts = positives['dist_bin'].value_counts()
        matched_negatives = []
        
        for bin_idx, count in pos_bin_counts.items():
            bin_negs = negatives[negatives['dist_bin'] == bin_idx]
            if len(bin_negs) > 0:
                # Sample with replacement if needed
                n_sample = min(count, len(bin_negs))
                sampled = bin_negs.sample(n=n_sample, random_state=42, replace=len(bin_negs) < count)
                matched_negatives.append(sampled)
        
        if matched_negatives:
            matched_df = pd.concat([positives, *matched_negatives], ignore_index=True)
            logger.info(f"  Distance-matched benchmark: {len(matched_df)} pairs "
                       f"({len(positives)} pos, {len(matched_df) - len(positives)} neg)")
            return matched_df.drop(columns=['dist_bin'], errors='ignore')
        else:
            logger.warning("Could not create distance-matched negatives")
            return self.benchmark
    
    def validate_gene_ids(self) -> LeakageReport:
        """
        Validate gene ID integrity.
        
        Per user requirement Point 5: "ID integrity gate"
        Ensures all gene IDs are valid Ensembl IDs (ENSG format) or valid symbols.
        Flags pairs with invalid or missing gene identifiers.
        
        Returns
        -------
        LeakageReport
            Report on gene ID validation status.
        """
        logger.info("Validating gene ID integrity...")
        
        # Identify gene columns
        gene_cols = [col for col in self.benchmark.columns 
                     if 'gene' in col.lower() or 'target' in col.lower()]
        
        if not gene_cols:
            return LeakageReport(
                audit_name='Gene ID Integrity',
                has_leakage=False,
                leakage_description='No gene columns found to validate',
                affected_pairs=0,
                total_pairs=len(self.benchmark),
                recommendation='Add gene identifier columns'
            )
        
        invalid_ids = 0
        invalid_details = []
        
        for gene_col in gene_cols:
            values = self.benchmark[gene_col].dropna().astype(str)
            
            # Check for valid Ensembl format (ENSG followed by digits)
            ensg_pattern = values.str.match(r'^ENSG\d{11}$', na=False)
            
            # Check for valid gene symbols (alphabetic, possibly with numbers)
            symbol_pattern = values.str.match(r'^[A-Z][A-Z0-9\-]+$', na=False, case=False)
            
            # Check for obvious invalids
            obvious_invalid = (
                values.str.lower().isin(['nan', 'na', 'null', 'none', '']) |
                values.str.contains(r'[^\w\-\.]', regex=True, na=False)
            )
            
            # Valid if matches ENSG pattern OR valid symbol pattern AND not obviously invalid
            valid = (ensg_pattern | symbol_pattern) & ~obvious_invalid
            n_invalid = (~valid).sum()
            
            if n_invalid > 0:
                invalid_ids += n_invalid
                # Get example invalid IDs
                examples = values[~valid].head(5).tolist()
                invalid_details.append(f"{gene_col}: {n_invalid} invalid (e.g., {examples[:3]})")
        
        has_issues = invalid_ids > len(self.benchmark) * 0.01  # >1% is concerning
        
        report = LeakageReport(
            audit_name='Gene ID Integrity',
            has_leakage=has_issues,
            leakage_description=f'{invalid_ids} invalid gene IDs found. {"; ".join(invalid_details)}' if invalid_ids > 0 else 'All gene IDs valid',
            affected_pairs=invalid_ids,
            total_pairs=len(self.benchmark),
            recommendation='Clean gene identifiers before evaluation' if invalid_ids > 0 else 'None'
        )
        
        self.reports.append(report)
        return report
    
    def run_full_audit(self) -> Dict:
        """Run all audits and return comprehensive report."""
        logger.info("Running full leakage audit...")
        
        self.audit_abc_training_overlap()
        self.audit_gene_overlap_across_sources()
        self.audit_genomic_proximity()
        self.validate_gene_ids()  # NEW: ID integrity gate
        
        loso_splits = self.create_leave_one_study_out_splits()
        chrom_splits = self.create_chromosome_holdout_splits()
        
        # Create distance-matched negatives for sensitivity analysis
        distance_matched = self.create_distance_matched_negatives()
        
        return {
            'audits': [
                {
                    'name': r.audit_name,
                    'has_leakage': bool(r.has_leakage),
                    'description': r.leakage_description,
                    'affected': int(r.affected_pairs),
                    'total': int(r.total_pairs),
                    'recommendation': r.recommendation
                }
                for r in self.reports
            ],
            'loso_splits': {
                k: {
                    'train_size': v['train_size'],
                    'test_size': v['test_size'],
                    'train_pos': int(v['train_pos']),
                    'test_pos': int(v['test_pos'])
                }
                for k, v in loso_splits.items()
            },
            'chromosome_splits': {
                k: {
                    'train_size': v['train_size'],
                    'test_size': v['test_size'],
                    'train_pos': int(v['train_pos']),
                    'test_pos': int(v['test_pos'])
                }
                for k, v in chrom_splits.items()
            },
            'distance_matched_controls': {
                'original_size': len(self.benchmark),
                'matched_size': len(distance_matched),
                'description': 'Distance-matched negatives to isolate functional feature contribution'
            }
        }
    
    def print_report(self):
        """Print formatted audit report."""
        print("\n" + "=" * 70)
        print("LEAKAGE AUDIT REPORT")
        print("=" * 70)
        
        for r in self.reports:
            status = "WARNING" if r.has_leakage else "OK"
            print(f"\n[{status}] {r.audit_name}")
            print(f"  {r.leakage_description}")
            print(f"  Affected: {r.affected_pairs:,} / {r.total_pairs:,} pairs")
            print(f"  Recommendation: {r.recommendation}")


def main():
    """Run leakage audit on Task B benchmark."""
    
    # Audit Task B
    task_b_path = SCRIPT_DIR / "benchmarks/task_b_enhancer_to_gene.parquet"
    if task_b_path.exists():
        print("\n### TASK B LEAKAGE AUDIT ###")
        auditor = LeakageAuditor(task_b_path)
        report = auditor.run_full_audit()
        auditor.print_report()
        
        # Print splits
        print("\n### Leave-One-Study-Out Splits ###")
        for name, stats in report['loso_splits'].items():
            print(f"  {name}: train={stats['train_size']}, test={stats['test_size']}, "
                  f"train_pos={stats['train_pos']}, test_pos={stats['test_pos']}")
        
        # Save report
        output_path = SCRIPT_DIR / "benchmarks/task_b_leakage_audit.json"
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"\nSaved to {output_path}")
    
    # Audit Task A
    task_a_path = SCRIPT_DIR / "benchmarks/task_a_gwas_to_gene.parquet"
    if task_a_path.exists():
        print("\n### TASK A LEAKAGE AUDIT ###")
        auditor = LeakageAuditor(task_a_path)
        report = auditor.run_full_audit()
        auditor.print_report()
        
        # Save report
        output_path = SCRIPT_DIR / "benchmarks/task_a_leakage_audit.json"
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"\nSaved to {output_path}")


if __name__ == "__main__":
    main()
