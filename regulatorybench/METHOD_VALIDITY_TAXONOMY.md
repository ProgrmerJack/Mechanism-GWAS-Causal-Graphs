# Method Validity Taxonomy Table

**Supplementary Table for RegulatoryBench Analysis**

This table provides a comprehensive overview of 23 methods evaluated in RegulatoryBench, 
including their task applicability, training data sources, leakage status, and recommended 
evaluation splits.

---

## Table 1: Method Validity Taxonomy

| Method | Task A (GWAS→Gene) | Task B (E2G) | Task C (V2G) | Training Data | Leakage Risk | Recommended Split |
|--------|-------------------|--------------|--------------|---------------|--------------|-------------------|
| **Distance-Based** |
| Distance (Rank) | ✓ Primary | ✓ Secondary | ✓ | None | ○ Low | All |
| Distance to Gene Body | ✓ | ✓ | ✓ | None | ○ Low | All |
| Distance to TSS | ✓ | ✓ | ✓ | None | ○ Low | All |
| Distance to Promoter | ✓ | ✓ | ✓ | None | ○ Low | All |
| **ABC/Enhancer-Based** |
| ABC (Max Score) | ⚠ Secondary | ✓ Primary | ✗ | Fulco 2019 K562 | ● High | Held-out cell types |
| ABC (Prediction) | ⚠ Secondary | ✓ | ✗ | Fulco 2019 K562 | ● High | Held-out cell types |
| ABC (Loops) | ⚠ | ✓ | ✗ | Fulco 2019 K562 | ● High | Held-out cell types |
| rE2G | ⚠ | ✓ | ✗ | Multiple CRISPRi | ◐ Medium | LOSO |
| **GWAS-Derived** |
| PoPS | ✓ Primary | ✗ | ⚠ | MAGMA z-scores | ○ Low | All GWAS loci |
| FLAMES | ✓ | ⚠ | ⚠ | Multi-source | ◐ Medium | ExWAS benchmark |
| CALDERA | ✓ | ✗ | ⚠ | Open Targets | ◐ Medium | Held-out traits |
| L2G | ✓ | ⚠ | ⚠ | Open Targets | ◐ Medium | Held-out traits |
| **eQTL-Based** |
| eQTL CTS | ✗ | ✗ | ✓ Primary | GTEx | ○ Low | All tissues |
| eQTL Max | ✗ | ✗ | ✓ | GTEx | ○ Low | All tissues |
| Coloc | ⚠ | ⚠ | ✓ | GTEx | ◐ Medium | Held-out tissues |
| SMR | ⚠ | ⚠ | ✓ | GTEx | ◐ Medium | Held-out tissues |
| **Functional/Constraint** |
| Fishilevich 2017 | ⚠ | ✓ | ⚠ | GeneHancer | ◐ Medium | LOSO |
| pLI | ✓ | ⚠ | ✓ | gnomAD | ○ Low | All |
| LOEUF | ✓ | ⚠ | ✓ | gnomAD | ○ Low | All |
| EDS | ⚠ | ✓ | ⚠ | ENCODE | ◐ Medium | Held-out cell types |
| Coding Variants | ✓ | ✗ | ⚠ | None | ○ Low | All |

---

## Legend

**Task Applicability:**
- ✓ Primary = Method was designed for this task; recommended evaluation
- ✓ Secondary = Method can be applied but wasn't designed for this task
- ⚠ = Use with caution; may give misleading results
- ✗ = Inappropriate for this task; do not evaluate

**Leakage Risk:**
- ○ Low = No training data overlap with common benchmarks
- ◐ Medium = Partial overlap or gene-level leakage possible
- ● High = Direct training data overlap; use held-out splits only

**Tasks:**
- Task A (GWAS→Gene): Predicting causal genes from GWAS credible sets
- Task B (E2G): Predicting enhancer-gene regulatory connections
- Task C (V2G): Predicting variant effects on specific genes

---

## Recommended Evaluation Protocols

### Task A: GWAS-to-Gene Mapping

**Recommended methods:** Distance, PoPS, FLAMES, CALDERA, L2G

**Benchmark:** UK Biobank credible sets with experimental validation (14,016 pairs)

**Metrics:**
1. AUPRC (primary) - accounts for class imbalance
2. AUROC (secondary) - overall discrimination
3. Top-1 accuracy - per-locus ranking
4. Drug target OR - external validation per Ji et al. (2025)

**Recommended splits:**
- Full dataset for methods without training overlap
- Leave-one-trait-out for PoPS, L2G
- Leave-one-study-out for multi-source methods

### Task B: Enhancer-to-Gene Linking

**Recommended methods:** ABC, rE2G, Distance, EDS

**Benchmark:** CRISPRi validation data (19,825 pairs)

**CRITICAL:** ABC was trained on Fulco 2019 K562 data. Evaluate only on:
- ENCODE held-out cell types (4,378 pairs)
- Non-K562 sources

**Metrics:**
1. AUPRC (primary)
2. AUROC (secondary)
3. Distance-matched AUROC (controls for proximity confound)

### Task C: Variant-to-Gene Assignment

**Recommended methods:** eQTL-based, SMR, Coloc

**Benchmark:** MPRA validated variants (future work)

---

## References

1. Fulco CP, et al. Activity-by-contact model. Nat Genet 51, 1664–1669 (2019).
2. Schipper M, et al. FLAMES. Nat Genet 57, 323–333 (2025).
3. Mountjoy E, et al. L2G. Nat Genet 53, 1527–1533 (2021).
4. Weeks EM, et al. PoPS. Nat Genet 55, 1267–1276 (2023).
5. Schipper M, et al. CALDERA. medRxiv 10.1101/2024.07.26.24311057 (2024).
6. Ji C, et al. Drug target benchmarking. medRxiv 10.1101/2025.09.23.25336370 (2025).

---

## Data Availability

All processed benchmarks available at Zenodo: [DOI to be assigned]

Method implementations available at GitHub: [URL to be added]
