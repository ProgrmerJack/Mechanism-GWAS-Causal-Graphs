# RegulatoryBench: Nature Genetics Submission Summary

**Date:** January 2025  
**Status:** Ready for Submission

---

## Executive Summary

RegulatoryBench has been transformed into a Nature Genetics-grade Analysis submission implementing all 8 points of the comprehensive enhancement plan.

### Core Positioning

> "Benchmarks for GWAS gene prioritization are routinely distorted by task conflation and train–test leakage; when evaluated on leakage-audited experimental ground truth, proximity baselines remain hard to beat"

---

## Completed Enhancements

### 1. ✅ Format Compliance (Point 0)

- **Word count:** ~2,600 words (within 3,000 limit)
- **Display items:** 6 (5 figures + 1 table)
- **Structure:** Abstract → Introduction → Results → Discussion → Methods

### 2. ✅ Novelty Clarification (Point 1)

Three-part novelty claim clearly articulated:

| Novelty Component | Evidence |
|-------------------|----------|
| Task Separation | Taxonomy distinguishes GWAS→Gene vs Enhancer→Gene vs Variant→Gene |
| Leakage as Failure Mode | 78% ABC training overlap in Task B; held-out AUROC drops from 0.885 to 0.826 |
| Experimental Ground Truth | CRISPRi validation from ENCODE Phase 4 + drug targets from Open Targets |

### 3. ✅ SOTA Baselines Added (Point 2)

| Method | Source | Status |
|--------|--------|--------|
| FLAMES | PMID 39930082 (Jan 2025) | Mentioned as 2025 SOTA (AUPRC=55.8%) |
| CALDERA | PMC11312651 (Jul 2024) | Mentioned as 2024 baseline (AUPRC=84.4%) |
| cS2G | PMID 36646662 | Referenced in taxonomy |
| Open Targets L2G | platform.opentargets.org | Included via Ji et al. comparison |

### 4. ✅ Drug-Target Validation Aligned (Point 3)

Protocol standardized per Ji et al. 2025 medRxiv preprint:

| Method | Our OR | Ji et al. OR (reference) |
|--------|--------|--------------------------|
| Nearest gene baseline | 8.62 | 3.08 |
| L2G score | — | 3.14 (no improvement) |
| eQTL colocalization | — | 1.61 (not significant) |
| PoPS | 28.28*** | — |
| ABC Prediction | 0.64 (ns) | — |

**Key finding:** L2G ≤2% improvement over nearest gene in Ji et al., consistent with our Distance dominance.

### 5. ✅ Desk-Rejection Preemption (Point 4)

| Potential Critique | Preemption |
|-------------------|------------|
| "Just distance" | Distance reflects evolutionary constraint; explained mechanistically |
| "Limited novelty" | First systematic task separation + leakage audit framework |
| "Small benchmark" | 14,016 pairs, 560 loci, externally validated with drug targets |
| "Cherry-picked methods" | 23 methods evaluated across 3 task types |

### 6. ✅ Reviewer Trapdoor Analyses (Point 5)

Four pre-emptive analyses implemented:

| Analysis | Implementation | File |
|----------|---------------|------|
| LOSO splits | `create_loso_splits()` | leakage_audit.py |
| Distance-matched controls | `create_distance_matched_negatives()` | leakage_audit.py |
| AUPRC as primary | `compute_auprc_primary_metrics()` | enhanced_evaluation.py |
| Gene ID integrity gate | `validate_gene_ids()` | leakage_audit.py |

### 7. ✅ 6 Display Items (Point 6)

| Item | Content | File |
|------|---------|------|
| Fig 1 | Task taxonomy + leakage schematic | figure1_task_taxonomy.pdf |
| Fig 2 | Task A performance with 95% CIs | figure2_task_a_enhanced.pdf |
| Fig 3 | External validation (Ji et al. 2025) | figure3_external_validation.pdf |
| Fig 4 | Task B leakage sensitivity | figure4_task_b_leakage.pdf |
| Fig 5 | Benchmark integrity checklist | figure5_benchmark_integrity.pdf |
| Table 1 | Method validity taxonomy (23×3) | METHOD_VALIDITY_TAXONOMY.md |

### 8. ✅ NG Policy Compliance (Point 7)

| Requirement | Status |
|-------------|--------|
| Zenodo deposition | Placeholder DOI in manuscript |
| Container/environment | requirements.txt with pinned versions |
| Data availability | Complete statement in manuscript |
| Code availability | GitHub URL placeholder |
| Reproducibility card | REPRODUCIBILITY_CARD.md created |
| CITATION.cff | Present and complete |

---

## Key Results

### Task A: GWAS Credible Set → Causal Gene

| Method | AUROC [95% CI] | AUPRC | Top-1 Accuracy |
|--------|---------------|-------|----------------|
| Distance (Rank) | 0.930 [0.922-0.938] | 0.367 | 50.2% |
| Distance to Body | 0.873 [0.858-0.888] | 0.386 | 54.3% |
| PoPS | 0.786 [0.767-0.806] | 0.259 | 36.6% |
| ABC (Max) | 0.599 [0.581-0.617] | 0.128 | 26.4% |
| ABC Prediction | 0.474 [0.456-0.492] | 0.039 | 6.2% |

### Drug Target Enrichment

| Method | Odds Ratio | p-value |
|--------|-----------|---------|
| PoPS | 28.28 | <0.001 |
| Distance | 8.62 | <0.001 |
| ABC Prediction | 0.64 | 0.768 (ns) |

### Task B: Leakage Impact

| Condition | ABC AUROC | Distance AUROC |
|-----------|-----------|----------------|
| Full data | 0.885 | 0.877 |
| Held-out only | 0.826 | 0.826 |
| Leakage inflation | +7% | +6% |

---

## Files Modified/Created

### Core Manuscript
- `MANUSCRIPT_DRAFT.md` - Complete Nature Genetics Analysis manuscript

### Code Enhancements
- `leakage_audit.py` - Added distance-matched controls, gene ID validation
- `enhanced_evaluation.py` - Added Ji protocol validation, AUPRC metrics
- `generate_enhanced_figures.py` - Added Ji et al. external validation figure
- `generate_figures.py` - Fixed figure generation bug
- `requirements.txt` - Updated with scipy dependency

### New Documentation
- `REPRODUCIBILITY_CARD.md` - Complete reproducibility documentation
- `METHOD_VALIDITY_TAXONOMY.md` - 23 methods × 3 tasks validity matrix
- `SUBMISSION_SUMMARY.md` - This file

---

## Pipeline Verification

```
============================================================
PIPELINE SUMMARY
============================================================
  task_a: ✓ PASSED
  task_b: ✓ PASSED
  leakage: ✓ PASSED
  enhanced: ✓ PASSED
  figures: ✓ PASSED
  enhanced_figures: ✓ PASSED

✓ All steps completed successfully!
```

---

## Next Steps for Author

1. **Zenodo:** Upload benchmarks and obtain DOI
2. **GitHub:** Make repository public and add URL
3. **Author info:** Complete author names, affiliations, contributions
4. **ORCID:** Add author ORCID identifiers
5. **Cover letter:** Draft for Nature Genetics editors
6. **Formatting:** Convert markdown to journal submission format

---

## External Validation References

### Ji et al. 2025 (Critical External Support)
- **Title:** "Odds of Drug Approval by Gene Prioritization Methods"
- **DOI:** 10.1101/2025.09.23.25336370
- **Key finding:** Nearest gene (OR=3.08) ≈ L2G (OR=3.14), eQTL not significant

### FLAMES 2025
- **PMID:** 39930082
- **Key finding:** AUPRC=55.8% on enhancer-gene task

### CALDERA 2024
- **PMC:** PMC11312651
- **Key finding:** AUPRC=84.4% with functional annotations

---

## Submission Checklist

- [x] Word count under 3,000
- [x] 6 or fewer display items
- [x] All figures generated as PDF
- [x] Data availability statement
- [x] Code availability statement
- [x] Reproducibility documentation
- [x] Pipeline runs successfully
- [x] External validation cited
- [ ] Zenodo DOI obtained
- [ ] GitHub repository public
- [ ] Author information complete
- [ ] Cover letter drafted
