#!/usr/bin/env python3
"""
RegulatoryBench Task Taxonomy Framework

This module defines the three distinct evaluation tasks for regulatory element
to gene assignment benchmarking, following Nature Genetics Analysis standards.

Task A: GWAS Credible Set → Causal Gene (L2G, FLAMES domain)
Task B: Regulatory Element → Target Gene (CRISPRi ground truth)
Task C: Regulatory Variant → Affected Gene (MPRA ground truth)

These tasks are NOT interchangeable. Methods designed for one task may be
inapplicable or degenerate on others.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set
from enum import Enum
from pathlib import Path
import pandas as pd
import numpy as np


class TaskType(Enum):
    """The three fundamental task types in regulatory-to-gene prediction."""
    
    TASK_A = "GWAS_credible_set_to_gene"
    TASK_B = "regulatory_element_to_gene"  
    TASK_C = "regulatory_variant_to_gene"


@dataclass
class TaskDefinition:
    """Formal definition of an evaluation task."""
    
    task_type: TaskType
    description: str
    input_type: str
    output_type: str
    ground_truth_source: str
    applicable_methods: List[str]
    inapplicable_methods: List[str]
    required_inputs: List[str]
    primary_metric: str
    secondary_metrics: List[str]
    

# Task A: GWAS Credible Set → Causal Gene
TASK_A = TaskDefinition(
    task_type=TaskType.TASK_A,
    description="Given a GWAS locus with credible set variants and fine-mapping "
                "posterior probabilities, predict the causal effector gene(s).",
    input_type="Credible set variants + PIPs + LD structure + GWAS summary stats",
    output_type="Ranked genes per locus with causal scores",
    ground_truth_source="Known causal genes from literature curation, drug targets, "
                        "or functional validation at GWAS loci",
    applicable_methods=[
        "Open Targets L2G",
        "FLAMES",
        "PoPS",
        "DEPICT",
        "MAGMA-based gene scores",
        "Colocalization-based rankers (coloc, enloc)",
        "Distance-to-TSS baseline"
    ],
    inapplicable_methods=[
        "ABC model (no GWAS input)",
        "rE2G (no GWAS input)",
        "Pure enhancer-gene linkers"
    ],
    required_inputs=[
        "GWAS summary statistics",
        "Fine-mapping results (PIPs)",
        "Credible set definitions",
        "LD reference panel"
    ],
    primary_metric="Mean Reciprocal Rank (MRR)",
    secondary_metrics=["Recall@K", "AUROC (binary)", "Precision@1"]
)


# Task B: Regulatory Element → Target Gene  
TASK_B = TaskDefinition(
    task_type=TaskType.TASK_B,
    description="Given a regulatory element (enhancer) defined by genomic coordinates "
                "and chromatin context, predict its target gene(s).",
    input_type="Enhancer coordinates + chromatin accessibility + histone marks",
    output_type="Target gene(s) with linking scores",
    ground_truth_source="CRISPRi perturbation screens (ENCODE/Engreitz, Gasperini, etc.)",
    applicable_methods=[
        "ABC model",
        "rE2G",
        "EpiMap enhancer-gene links",
        "GeneHancer",
        "HACER",
        "Distance-to-TSS baseline",
        "Promoter capture Hi-C"
    ],
    inapplicable_methods=[
        "L2G (requires GWAS input)",
        "FLAMES (requires GWAS input)",
        "cS2G (SNP-specific, gene-level output)"
    ],
    required_inputs=[
        "Enhancer coordinates (chr, start, end)",
        "Chromatin accessibility (ATAC-seq)",
        "Gene annotations (TSS positions)"
    ],
    primary_metric="AUROC (within-locus discrimination)",
    secondary_metrics=["AUPRC", "Recall@K", "Distance-stratified AUC"]
)


# Task C: Regulatory Variant → Affected Gene
TASK_C = TaskDefinition(
    task_type=TaskType.TASK_C,
    description="Given a regulatory variant with measured allelic effect (MPRA), "
                "predict which gene's expression is affected.",
    input_type="Variant coordinates + MPRA effect size + allele information",
    output_type="Target gene(s) with effect predictions",
    ground_truth_source="MPRA experiments with variant-to-gene mappings "
                        "(MPRAVarDB, Tewhey et al., etc.)",
    applicable_methods=[
        "Distance-to-TSS baseline",
        "eQTL colocalization",
        "Variant effect predictors (VEP, CADD)",
        "ABC model (if variant overlaps enhancer)",
        "cS2G (for common variants in eQTL panels)"
    ],
    inapplicable_methods=[
        "L2G (no GWAS context)",
        "FLAMES (no GWAS context)"
    ],
    required_inputs=[
        "Variant coordinates (chr, pos)",
        "MPRA activity evidence",
        "Gene annotations"
    ],
    primary_metric="AUROC (variant-gene pair classification)",
    secondary_metrics=["AUPRC", "Distance-stratified AUC", "Effect size correlation"]
)


def assign_task_type(evidence_source: str, has_gwas_context: bool = False) -> TaskType:
    """
    Determine the appropriate task type for a given evidence source.
    
    Args:
        evidence_source: Source of ground truth (CRISPRi, MPRA, GWAS, etc.)
        has_gwas_context: Whether GWAS summary statistics are available
        
    Returns:
        TaskType appropriate for this evaluation
    """
    if has_gwas_context and evidence_source in ['GWAS_gold', 'drug_target', 'known_causal']:
        return TaskType.TASK_A
    elif evidence_source in ['CRISPRi', 'ENCODE_CRISPRi', 'Gasperini', 'Fulco']:
        return TaskType.TASK_B
    elif evidence_source in ['MPRA', 'MPRAVarDB', 'Tewhey', 'Ulirsch']:
        return TaskType.TASK_C
    else:
        raise ValueError(f"Unknown evidence source: {evidence_source}")


def check_method_applicability(method: str, task_type: TaskType) -> Dict[str, any]:
    """
    Check if a method is applicable to a given task type.
    
    Returns dict with:
        - applicable: bool
        - reason: str explanation
        - warning: optional limitation notes
    """
    task_def = {
        TaskType.TASK_A: TASK_A,
        TaskType.TASK_B: TASK_B,
        TaskType.TASK_C: TASK_C
    }[task_type]
    
    if method in task_def.applicable_methods:
        return {
            "applicable": True,
            "reason": f"{method} is designed for {task_type.value}",
            "warning": None
        }
    elif method in task_def.inapplicable_methods:
        return {
            "applicable": False,
            "reason": f"{method} requires inputs not available for {task_type.value}",
            "warning": "Applying this method will yield degenerate or meaningless results"
        }
    else:
        return {
            "applicable": None,
            "reason": f"{method} applicability to {task_type.value} is uncertain",
            "warning": "Requires manual review of method inputs and outputs"
        }


def get_evaluation_contract() -> str:
    """
    Return the formal Evaluation Contract for RegulatoryBench.
    This document defines what can and cannot be compared.
    """
    contract = """
================================================================================
REGULATORYBENCH EVALUATION CONTRACT
================================================================================

This document defines the evaluation protocol for RegulatoryBench, establishing
clear boundaries between task types to prevent invalid method comparisons.

1. TASK DEFINITIONS
-------------------

TASK A: GWAS Credible Set → Causal Gene
  - Input: GWAS locus with credible set variants, PIPs, LD structure
  - Output: Ranked genes per locus
  - Methods: L2G, FLAMES, PoPS, DEPICT, colocalization
  - Evaluation: MRR, Recall@K, Precision@1

TASK B: Regulatory Element → Target Gene  
  - Input: Enhancer coordinates, chromatin context
  - Output: Target gene predictions
  - Methods: ABC, rE2G, distance baselines, Hi-C
  - Evaluation: AUROC, AUPRC (within-locus)

TASK C: Regulatory Variant → Affected Gene
  - Input: Variant coordinates, MPRA evidence
  - Output: Gene effect predictions
  - Methods: Distance baselines, eQTL coloc, VEP
  - Evaluation: AUROC, AUPRC (variant-gene pairs)

2. CROSS-TASK COMPARISONS
-------------------------

PROHIBITED comparisons:
  - L2G/FLAMES on CRISPRi data (no GWAS context → inapplicable)
  - cS2G on arbitrary enhancers (gene-level max → no within-locus discrimination)
  - ABC on GWAS loci without enhancer mapping (input mismatch)

PERMITTED comparisons:
  - Distance baselines across all tasks (universally applicable)
  - Within-task method comparisons (same input requirements)
  - Sensitivity analyses with explicit caveats

3. GROUND TRUTH BENCHMARKS
--------------------------

CRISPR-E2G (Task B primary):
  - Source: ENCODE/Engreitz CRISPRi perturbation screens
  - N entries: 661 enhancer-gene pairs
  - Ground truth: Significant expression changes upon enhancer knockdown

MPRA-V2G (Task C primary):
  - Source: MPRA variant effect databases
  - N entries: 8,601 variant-gene pairs  
  - Ground truth: Significant allelic activity differences

4. CANDIDATE GENERATION
-----------------------

Per-locus candidate sets:
  - Window: ±500kb from regulatory element/variant
  - Gene set: GENCODE v43 protein-coding genes
  - Positive rate: ~4% (true pairs / total candidates)

5. METRICS
----------

Primary: AUROC (within-locus binary classification)
  - For each locus: rank candidates by method score
  - Compute AUC treating true target as positive, all others as negative

Secondary:
  - AUPRC (precision-recall, sensitive to class imbalance)
  - Distance-stratified AUC (0-10kb, 10-100kb, 100-500kb)
  - Evidence-type stratified AUC (CRISPRi vs MPRA)

6. INDEPENDENCE REQUIREMENTS
----------------------------

Training overlap audit:
  - Methods must not have been trained on benchmark loci
  - L2G: 0% overlap verified

Evidence layer audit:
  - Labels must not derive from evaluated method's outputs
  - CRISPRi/MPRA labels are experimentally derived (independent)

Publication cutoff audit:
  - For generalization claims, test on post-cutoff data

7. REPORTING STANDARDS
----------------------

Required for all method evaluations:
  1. Task type (A/B/C) and justification
  2. Applicability statement for method on task
  3. Primary and secondary metrics
  4. Stratified results (distance, evidence type)
  5. Coverage rate (% of loci with method scores)
  6. Independence verification

================================================================================
"""
    return contract


def print_task_summary():
    """Print a summary of all task definitions."""
    for task_def in [TASK_A, TASK_B, TASK_C]:
        print(f"\n{'='*70}")
        print(f"{task_def.task_type.value}")
        print(f"{'='*70}")
        print(f"\nDescription: {task_def.description}")
        print(f"\nInput: {task_def.input_type}")
        print(f"Output: {task_def.output_type}")
        print(f"Ground Truth: {task_def.ground_truth_source}")
        print(f"\nApplicable Methods:")
        for m in task_def.applicable_methods:
            print(f"  ✓ {m}")
        print(f"\nInapplicable Methods:")
        for m in task_def.inapplicable_methods:
            print(f"  ✗ {m}")
        print(f"\nPrimary Metric: {task_def.primary_metric}")


if __name__ == "__main__":
    print_task_summary()
    print("\n")
    print(get_evaluation_contract())
