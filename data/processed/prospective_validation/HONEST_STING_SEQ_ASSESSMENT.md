# STING-seq Validation: Honest Assessment

**Date:** 2025-01-21 (Updated with training leakage analysis)

## ⚠️ CRITICAL UPDATE: Training-Test Leakage Detected

**After rigorous provenance verification, we discovered 21/120 (17.5%) of STING-seq genes with L2G scores overlap with L2G training gold standards.**

### True Prospective Performance:

| Evaluation | AUROC | 95% CI | N (positives) |
|------------|-------|--------|---------------|
| **All genes (includes training overlap)** | 0.756 | [0.708 - 0.801] | 120 |
| **Prospective only (TRUE validation)** | **0.727** | [0.670 - 0.780] | 99 |

**⚠️ The claimed AUROC is inflated by ~3 percentage points due to training leakage.**

### Overlapping Genes (must be excluded):
CFH, HFE, CPS1, THPO, TF, JAK2, IL12B, PDGFRB, BCL11A, APOE, PTPN22, IL2RA, TBX3, MYC, AKT1, EGFR, BMP6, RUNX1, CCR5, CD52, F2

---

## Updated Summary

Our rigorous validation on the STING-seq benchmark (Morris et al. Science 2023) reveals:

| Method | AUROC (All) | AUROC (Prospective) | Note |
|--------|-------------|---------------------|------|
| **L2G (Open Targets 22.09)** | 0.756 | **0.727** | TRUE prospective |
| Path-probability (L2G only) | 0.756 | 0.727 | Equivalent to L2G |
| Path-probability (L2G + ABC) | 0.62 | ~0.60 | **Hurts** performance |
| cS2G | 0.62 | - | Baseline |
| NearestGene | 0.50 | - | Random |

## Key Findings

### 1. Path-Probability = L2G on This Benchmark

When we compute path-probability using only L2G through noisy-OR:
```
P = 1 - (1-ε)(1-P_L2G)
```

The transformation is monotonic, so ranking is preserved and AUROC is identical.
**The noisy-OR framework only adds value when we have MULTIPLE independent evidence sources.**

### 2. Adding ABC HURTS Performance

Counterintuitively, adding ABC scores via noisy-OR degraded AUROC from 0.70 to 0.62.

**Why?**
- ABC provides gene-level max scores across **all 131 cell types** and **all enhancers**
- This introduces noise: many false positive enhancer-gene links
- L2G already incorporates enhancer information as one of its features
- The proper use of ABC requires **variant-specific** enhancer-gene links, not gene-level aggregates

### 3. What the Full Mechanism Graph Framework Actually Does

The complete mechanism graph framework (in `src/mechanism_graphs/`) uses:

1. **Variant → cCRE**: Fine-mapping PIPs weighted by specific enhancer overlap
2. **cCRE → Gene**: ABC score for the **specific enhancer-gene pair** (not gene-level max)
3. **Gene → Tissue**: eQTL colocalization (PP.H4)
4. **Tissue → Trait**: Tissue relevance prior

**This variant-specific integration is what makes the full framework valuable.**

However, we **cannot run the full framework on STING-seq** because:
- We don't have the specific variant → cCRE mappings for all STING-seq variants
- We don't have eQTL colocalization for all the blood traits
- The existing mechanism graphs (APOE, SORT1, TCF7L2) are for different loci

### 4. What This Means for the Paper

**Honest Claims We Can Make:**

1. ✅ "L2G achieves AUROC 0.70 on prospective STING-seq validation"
2. ✅ "Path-probability equals L2G when L2G is the only evidence source"
3. ✅ "The framework's value comes from integrating multiple variant-specific evidence sources"

**Claims We Cannot Make:**

1. ❌ "Path-probability outperforms L2G on STING-seq" (it doesn't)
2. ❌ "Adding ABC improves performance" (it hurts on this benchmark)
3. ❌ "Full framework validated on STING-seq" (we can only validate L2G component)

## Recommendations

### Option A: Honest Repositioning

Reframe the manuscript to:
1. Present L2G validation as the prospective external benchmark
2. Present mechanism graphs as a conceptual framework for interpretability
3. Show the full framework on the 3 curated loci (APOE, SORT1, TCF7L2) as exemplars
4. Be clear that full framework advantage requires variant-specific evidence integration

### Option B: Build Full Evidence Integration

Create variant-specific ABC lookups for STING-seq:
1. Map each STING-seq variant to overlapping ABC enhancers
2. Use only the specific enhancer-gene ABC scores for each variant
3. This requires significant additional data curation

### Option C: Use Alternative Benchmarks

Find benchmarks where we have full evidence (coloc + ABC + L2G + eQTL):
1. Drug target benchmarks with full colocalization data
2. Mendelian disease loci with complete annotation

## Data Provenance Note

- STING-seq: Morris et al. Science 2023, published May 19, 2023
- L2G: Open Targets Platform 22.09 (training cutoff unknown, but claimed pre-2021)
- ABC: Fulco et al. Nature Genetics 2019

**VERIFIED**: L2G 22.09 training gold standards contain 526 genes from GWAS gold standards.
**LEAKAGE DETECTED**: 21 of 120 STING-seq genes with L2G scores overlap with training.
**TRUE PROSPECTIVE AUROC**: 0.727 (excluding overlapping genes)

## Files

- Provenance verification: `scripts/verify_l2g_training_provenance.py`
- Leakage analysis: `scripts/evaluate_l2g_with_leakage_exclusion.py`
- Results: `data/processed/prospective_validation/l2g_auroc_with_leakage_analysis.json`
- Report: `data/processed/prospective_validation/L2G_TRAINING_PROVENANCE_REPORT.md`
