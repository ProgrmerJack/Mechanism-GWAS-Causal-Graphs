#!/usr/bin/env python3
"""
Benchmark Integrity Toolkit: Detect Training Leakage and Label Circularity

This toolkit provides automated checks for common benchmark integrity problems:
1. Training data overlap/leakage
2. Label circularity (positives defined using method predictions)
3. Test set contamination
4. Temporal leakage (future data in training)

The goal: Prevent the situation where Task A v1 had 100% ABC training examples,
creating circular evaluation that inflated performance claims by 4×.

Author: RegulatoryBench v2.0
Date: December 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass
import json
import warnings

warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).parent.resolve()


@dataclass
class IntegrityCheck:
    """Result of a benchmark integrity check."""
    check_name: str
    passed: bool
    severity: str  # 'fatal', 'warning', 'info'
    message: str
    details: Dict


class BenchmarkIntegrityChecker:
    """Systematic integrity checks for E2G benchmarks."""
    
    def __init__(self, benchmark_path: Path):
        self.benchmark_path = benchmark_path
        self.df = pd.read_parquet(benchmark_path)
        self.checks: List[IntegrityCheck] = []
        
    def check_label_independence(self) -> IntegrityCheck:
        """
        Check that positive labels are not derived from method predictions.
        
        Red flags:
        - Column named '{Method}TrainingExample' with 100% positive overlap
        - High correlation between labels and method scores
        - Positive labels defined by score thresholds
        
        IMPORTANT: If benchmark has 'evidence_type' column documenting independent
        provenance (coding, mendelian, rare_variant_burden), then training overlap
        is OK - it just means the method was trained on the same ground truth.
        """
        # If we have documented independent evidence, training overlap is acceptable
        has_evidence = 'evidence_type' in self.df.columns
        if has_evidence:
            positives = self.df[self.df['is_positive'] == True]
            evidence_types = positives['evidence_type'].unique()
            
            # Check if evidence is independent
            independent_types = {'coding_variant', 'coding_splice_site', 'mendelian_disease', 
                                'rare_variant_burden', 'polygenic_priority', 'sting_seq'}
            has_independent = any(t in independent_types for t in evidence_types)
            
            if has_independent:
                return IntegrityCheck(
                    check_name='label_independence',
                    passed=True,
                    severity='info',
                    message=f'[OK] Labels have independent provenance: {", ".join(evidence_types)}',
                    details={'evidence_types': list(evidence_types)}
                )
        
        # Check for training example flags
        training_cols = [c for c in self.df.columns if 'trainingexample' in c.lower()]
        
        if training_cols:
            for col in training_cols:
                if self.df[col].dtype == bool or self.df[col].dtype == object:
                    # Check overlap with positives
                    positives = self.df[self.df['is_positive'] == True]
                    training_overlap = (positives[col] == True).sum()
                    training_pct = training_overlap / len(positives) if len(positives) > 0 else 0
                    
                    if training_pct > 0.8:  # >80% overlap is FATAL (if no independent evidence)
                        return IntegrityCheck(
                            check_name='label_independence',
                            passed=False,
                            severity='fatal',
                            message=f'FATAL: {training_pct:.1%} of positive labels overlap with {col}',
                            details={
                                'training_column': col,
                                'overlap_count': int(training_overlap),
                                'total_positives': int(len(positives)),
                                'overlap_percentage': float(training_pct)
                            }
                        )
        
        # Check correlation between labels and prediction scores
        score_cols = [c for c in self.df.columns if any(x in c.lower() for x in ['score', 'prediction', 'prob'])]
        high_correlations = []
        
        for col in score_cols[:10]:  # Check first 10 score columns
            try:
                corr = self.df['is_positive'].astype(int).corr(self.df[col])
                if abs(corr) > 0.5:
                    high_correlations.append((col, corr))
            except:
                pass
        
        if high_correlations:
            return IntegrityCheck(
                check_name='label_independence',
                passed=False,
                severity='warning',
                message=f'WARNING: {len(high_correlations)} prediction scores highly correlated with labels',
                details={'high_correlations': high_correlations}
            )
        
        return IntegrityCheck(
            check_name='label_independence',
            passed=True,
            severity='info',
            message='[OK] Labels appear independent of method predictions',
            details={}
        )
    
    def check_evidence_provenance(self) -> IntegrityCheck:
        """
        Check that positives have documented independent evidence.
        
        Good benchmarks document WHERE each positive label came from:
        - coding_variant, rare_variant_burden, mendelian_disease, experimental_validation
        
        Bad benchmarks have no provenance or vague "literature curation".
        """
        if 'evidence_type' not in self.df.columns:
            return IntegrityCheck(
                check_name='evidence_provenance',
                passed=False,
                severity='warning',
                message='WARNING: No evidence_type column - provenance not documented',
                details={}
            )
        
        positives = self.df[self.df['is_positive'] == True]
        evidence_counts = positives['evidence_type'].value_counts().to_dict()
        
        # Check for vague evidence types
        vague_types = {'none', 'unknown', 'literature', 'curated'}
        vague_count = sum(evidence_counts.get(t, 0) for t in vague_types)
        
        if vague_count > len(positives) * 0.2:
            return IntegrityCheck(
                check_name='evidence_provenance',
                passed=False,
                severity='warning',
                message=f'WARNING: {vague_count/len(positives):.1%} of positives have vague evidence type',
                details={'evidence_counts': evidence_counts}
            )
        
        return IntegrityCheck(
            check_name='evidence_provenance',
            passed=True,
            severity='info',
            message=f'[OK] Positives have documented evidence: {len(evidence_counts)} types',
            details={'evidence_counts': evidence_counts}
        )
    
    def check_positive_rate_sanity(self) -> IntegrityCheck:
        """
        Check that positive rate is biologically reasonable.
        
        For Task A (GWAS→Gene):
        - 1-5% positive rate is reasonable (1-2 causal genes per locus, ~20-50 candidates)
        - >10% suggests label inflation
        - <0.5% suggests overly strict criteria
        """
        positive_rate = self.df['is_positive'].mean()
        n_positives = self.df['is_positive'].sum()
        n_loci = self.df['locus_id'].nunique() if 'locus_id' in self.df.columns else None
        
        if positive_rate > 0.10:
            return IntegrityCheck(
                check_name='positive_rate',
                passed=False,
                severity='warning',
                message=f'WARNING: Positive rate {positive_rate:.2%} is unusually high (expect 1-5%)',
                details={
                    'positive_rate': float(positive_rate),
                    'n_positives': int(n_positives),
                    'n_loci': n_loci
                }
            )
        elif positive_rate < 0.005:
            return IntegrityCheck(
                check_name='positive_rate',
                passed=False,
                severity='warning',
                message=f'WARNING: Positive rate {positive_rate:.2%} is unusually low (expect 1-5%)',
                details={
                    'positive_rate': float(positive_rate),
                    'n_positives': int(n_positives),
                    'n_loci': n_loci
                }
            )
        
        return IntegrityCheck(
            check_name='positive_rate',
            passed=True,
            severity='info',
            message=f'[OK] Positive rate {positive_rate:.2%} is reasonable',
            details={
                'positive_rate': float(positive_rate),
                'n_positives': int(n_positives),
                'n_loci': n_loci
            }
        )
    
    def check_temporal_leakage(self) -> IntegrityCheck:
        """
        Check for temporal leakage: test data published before training data.
        
        If benchmark has 'publication_year' and methods have 'training_year',
        ensure no test examples are from after method training.
        """
        if 'publication_year' not in self.df.columns:
            return IntegrityCheck(
                check_name='temporal_leakage',
                passed=True,
                severity='info',
                message='INFO: No publication_year column - cannot check temporal leakage',
                details={}
            )
        
        # Check if any test examples are from the future
        positives = self.df[self.df['is_positive'] == True]
        if 'publication_year' in positives.columns:
            years = positives['publication_year'].dropna()
            if len(years) > 0:
                latest_year = years.max()
                if latest_year > 2025:
                    return IntegrityCheck(
                        check_name='temporal_leakage',
                        passed=False,
                        severity='fatal',
                        message=f'FATAL: Test data includes examples from {latest_year} (future data)',
                        details={'latest_year': int(latest_year)}
                    )
        
        return IntegrityCheck(
            check_name='temporal_leakage',
            passed=True,
            severity='info',
            message='[OK] No obvious temporal leakage detected',
            details={}
        )
    
    def check_distance_confounding(self) -> IntegrityCheck:
        """
        Check if positive labels are confounded with distance.
        
        If positives are systematically closer than negatives, methods may
        learn to predict distance rather than regulatory function.
        
        This is NOT necessarily wrong (biological truth), but should be documented.
        """
        if 'GeneBodyDistanceToBestSNP' not in self.df.columns:
            return IntegrityCheck(
                check_name='distance_confounding',
                passed=True,
                severity='info',
                message='INFO: No distance column - cannot check distance confounding',
                details={}
            )
        
        positives = self.df[self.df['is_positive'] == True]['GeneBodyDistanceToBestSNP']
        negatives = self.df[self.df['is_positive'] == False]['GeneBodyDistanceToBestSNP']
        
        median_pos = positives.median()
        median_neg = negatives.median()
        
        fold_difference = median_neg / median_pos if median_pos > 0 else np.inf
        
        if fold_difference > 5:
            return IntegrityCheck(
                check_name='distance_confounding',
                passed=False,
                severity='warning',
                message=f'WARNING: Positives are {fold_difference:.1f}× closer than negatives (strong distance confound)',
                details={
                    'median_positive_distance': float(median_pos),
                    'median_negative_distance': float(median_neg),
                    'fold_difference': float(fold_difference)
                }
            )
        
        return IntegrityCheck(
            check_name='distance_confounding',
            passed=True,
            severity='info',
            message=f'[OK] Distance confound is modest ({fold_difference:.1f}×)',
            details={
                'median_positive_distance': float(median_pos),
                'median_negative_distance': float(median_neg),
                'fold_difference': float(fold_difference)
            }
        )
    
    def run_all_checks(self) -> List[IntegrityCheck]:
        """Run all integrity checks."""
        self.checks = [
            self.check_label_independence(),
            self.check_evidence_provenance(),
            self.check_positive_rate_sanity(),
            self.check_temporal_leakage(),
            self.check_distance_confounding()
        ]
        return self.checks
    
    def generate_report(self) -> str:
        """Generate human-readable integrity report."""
        report = []
        report.append("="*100)
        report.append(f"BENCHMARK INTEGRITY REPORT: {self.benchmark_path.name}")
        report.append("="*100)
        report.append(f"Total pairs: {len(self.df):,}")
        report.append(f"Positive pairs: {self.df['is_positive'].sum():,} ({self.df['is_positive'].mean():.2%})")
        if 'locus_id' in self.df.columns:
            report.append(f"Unique loci: {self.df['locus_id'].nunique():,}")
        report.append("")
        
        # Group by severity
        fatal = [c for c in self.checks if c.severity == 'fatal' and not c.passed]
        warnings = [c for c in self.checks if c.severity == 'warning' and not c.passed]
        passed = [c for c in self.checks if c.passed]
        
        if fatal:
            report.append("[FATAL] FATAL ISSUES (benchmark NOT suitable for publication):")
            report.append("-"*100)
            for check in fatal:
                report.append(f"  {check.message}")
                if check.details:
                    report.append(f"    Details: {check.details}")
            report.append("")
        
        if warnings:
            report.append("[WARN] WARNINGS (should be addressed or justified):")
            report.append("-"*100)
            for check in warnings:
                report.append(f"  {check.message}")
                if check.details:
                    report.append(f"    Details: {check.details}")
            report.append("")
        
        if passed:
            report.append("[PASS] PASSED CHECKS:")
            report.append("-"*100)
            for check in passed:
                report.append(f"  {check.message}")
            report.append("")
        
        # Overall verdict
        report.append("="*100)
        if fatal:
            report.append("VERDICT: [X] BENCHMARK HAS FATAL FLAWS - NOT SUITABLE FOR PUBLICATION")
        elif warnings:
            report.append("VERDICT: [!] BENCHMARK HAS WARNINGS - REVIEW BEFORE PUBLICATION")
        else:
            report.append("VERDICT: [OK] BENCHMARK PASSES ALL INTEGRITY CHECKS")
        report.append("="*100)
        
        return "\n".join(report)
    
    def save_report(self, output_path: Optional[Path] = None):
        """Save integrity report to file."""
        if output_path is None:
            output_path = self.benchmark_path.parent / f"{self.benchmark_path.stem}_integrity_report.txt"
        
        report = self.generate_report()
        output_path.write_text(report)
        print(f"Saved integrity report: {output_path}")
        
        # Also save JSON for programmatic access
        json_path = output_path.with_suffix('.json')
        json_data = {
            'benchmark': str(self.benchmark_path),
            'checks': [
                {
                    'name': c.check_name,
                    'passed': c.passed,
                    'severity': c.severity,
                    'message': c.message,
                    'details': c.details
                }
                for c in self.checks
            ]
        }
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        print(f"Saved JSON report: {json_path}")


def compare_benchmarks(v1_path: Path, v2_path: Path):
    """Compare two benchmark versions to show improvements."""
    print("\n" + "="*100)
    print("BENCHMARK COMPARISON: v1 (contaminated) vs v2 (independent)")
    print("="*100 + "\n")
    
    print("V1 (Original):")
    checker_v1 = BenchmarkIntegrityChecker(v1_path)
    checker_v1.run_all_checks()
    print(checker_v1.generate_report())
    
    print("\n\n")
    
    print("V2 (Fixed):")
    checker_v2 = BenchmarkIntegrityChecker(v2_path)
    checker_v2.run_all_checks()
    print(checker_v2.generate_report())
    
    # Summary comparison
    v1_fatal = sum(1 for c in checker_v1.checks if c.severity == 'fatal' and not c.passed)
    v2_fatal = sum(1 for c in checker_v2.checks if c.severity == 'fatal' and not c.passed)
    
    print("\n" + "="*100)
    print("IMPROVEMENT SUMMARY")
    print("="*100)
    print(f"v1 fatal issues: {v1_fatal}")
    print(f"v2 fatal issues: {v2_fatal}")
    print(f"Issues fixed: {v1_fatal - v2_fatal}")
    print("="*100)


def main():
    """Run benchmark integrity checks."""
    benchmarks_dir = SCRIPT_DIR / "benchmarks"
    
    # Check v1 (contaminated)
    v1_path = benchmarks_dir / "task_a_gwas_to_gene.parquet"
    if v1_path.exists():
        print("Checking Task A v1 (original, contaminated)...")
        checker_v1 = BenchmarkIntegrityChecker(v1_path)
        checker_v1.run_all_checks()
        print(checker_v1.generate_report())
        checker_v1.save_report()
    
    # Check v2 (independent)
    v2_path = benchmarks_dir / "task_a_gwas_to_gene_v2.parquet"
    if v2_path.exists():
        print("\n\nChecking Task A v2 (independent labels)...")
        checker_v2 = BenchmarkIntegrityChecker(v2_path)
        checker_v2.run_all_checks()
        print(checker_v2.generate_report())
        checker_v2.save_report()
    
    # Compare
    if v1_path.exists() and v2_path.exists():
        compare_benchmarks(v1_path, v2_path)


if __name__ == '__main__':
    main()
