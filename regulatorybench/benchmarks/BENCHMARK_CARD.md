# RegulatoryBench v4 - Task-Stratified Benchmark

Generated: 2025-12-20T18:18:21.869490

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

### GWAS-E2G (Task A)

- **Total pairs**: 14,016
- **Positives**: 569 (4.1%)
- **Negatives**: 13,447
- **Unique loci**: 569
- **Sources**: E2G_UKBB_benchmark

### CRISPR-E2G (Task B)

- **Total pairs**: 19,825
- **Positives**: 863 (4.4%)
- **Negatives**: 18,962
- **Unique loci**: 6,927
- **Sources**: ENCODE_CRISPRi_K562, ENCODE_CRISPRi_heldout, Fulco_2019

## Evaluation Protocol

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
