# Critical Peer-Review Analysis: PoPS Circularity Problem

**Date**: 2025-01-XX  
**Issue**: Task A v2 gold standard contains fatal circular evaluation for PoPS  
**Severity**: FATAL - Will cause immediate desk rejection at Nature Genetics

---

## Executive Summary

**The Problem**: v2 benchmark fixed ABC circularity but introduced PoPS circularity:
- 159 of 203 (78.3%) positive labels come from PoPS>0.95 computational predictions
- Then PoPS is evaluated on this benchmark, achieving 94.2% Top-1 accuracy
- This is **self-fulfilling evaluation**: PoPS "wins" because it defined the positives

**Manuscript Status**:
- Abstract claims "independent gold standards" ❌ FALSE for PoPS
- Methods line 268: Still says "ABC model predictions >0.015" ❌ OBSOLETE
- Results reports PoPS 94.2% without disclosing circularity ❌ MISLEADING
- Word count: ~4,300 (Nature Genetics limit ~4,000) ❌ OVER LIMIT

**Only 44 of 203 (21.7%) positives are truly independent** (Mendelian disease genes)

---

## 1. Evidence Breakdown (Current v2 Benchmark)

```
Total pairs: 14,016
Positives: 203 (1.45%)

Evidence sources:
- PoPS (polygenic_priority):  159 (78.3%) ← CIRCULAR ❌
- Mendelian disease:           44 (21.7%) ← CLEAN ✓
- Coding/LoF variants:          0 (0.0%)  ← NEED TO FIND
- Rare variant burden:          0 (0.0%)  ← NEED TO SEARCH
```

**Evidence manifest shows**:
- ALL 159 PoPS labels have: `evidence_source: "PoPS_high_confidence"`, `confidence: "medium"`
- PoPS scores range 1.03 to 2.36 (above 0.95 threshold)
- Then PoPS is evaluated → achieves 94.2% Top-1 accuracy = CIRCULAR

---

## 2. Performance Results (CONTAMINATED)

```
Method       Top-1   Top-5   Top-10  Top-20
PoPS         94.2%   96.8%   97.3%   98.1%  ← INFLATED (78% of labels from PoPS)
ABC          10.5%   18.9%   24.6%   32.5%
Distance     13.6%   24.1%   30.5%   38.9%
```

**Why PoPS performance is artificially inflated**:
1. Gold standard uses PoPS>0.95 as evidence source
2. Evaluation ranks genes by PoPS score
3. Obviously high overlap → high accuracy
4. **This is NOT generalization performance, it's self-grading**

---

## 3. Manuscript Internal Inconsistencies

### 3.1 Methods Section (Line 268-270)

**Current text**:
> "Ground truth labels derive from the intersection of: (1) ABC model predictions with score >0.015..."

**Problems**:
- ❌ References ABC predictions (we removed ABC contamination in v2)
- ❌ No mention of PoPS being used as evidence source
- ❌ Contradicts v2 fix that supposedly removed method-derived evidence
- ❌ Doesn't match actual evidence breakdown (159 PoPS, 44 Mendelian, 0 coding)

### 3.2 Abstract (Lines 7-15)

**Current text**:
> "We establish independent gold standard benchmarks..."

**Problems**:
- ❌ Claims "independent" but 78% of labels come from PoPS (an evaluated method)
- ❌ Reports PoPS 94.2% accuracy without disclosing circular evaluation
- ❌ Misleading representation of benchmark independence

### 3.3 Results Section (Lines 60-75)

**Current text**:
> "Task A gold standard comprises 203 positive (locus, gene) pairs with independent evidence..."
> "PoPS achieves 94.2% Top-1 accuracy..."

**Problems**:
- ❌ Emphasizes "independent evidence" but doesn't disclose PoPS circularity
- ❌ Reports PoPS performance as if it's fair comparison
- ❌ No mention that PoPS labels should be excluded from PoPS evaluation

### 3.4 Word Count

- Current: ~4,300 words (main text)
- Nature Genetics Article limit: ~4,000 words
- **Over by ~300 words** ❌

---

## 4. Why This Is a Peer-Review Disaster

### 4.1 Violates Core Benchmark Principles

From Weeks et al. (2023) and Mountjoy et al. (2021):
> "Gold standard labels must be independent of evaluated methods. Using method outputs as ground truth creates circular validation."

**We violated this principle**: PoPS is both:
1. Evidence source (159 genes with PoPS>0.95)
2. Evaluated method (ranked by PoPS score)

### 4.2 Parallels with Original ABC Problem

| Aspect | v1 (ABC circular) | v2 (PoPS circular) | What should be |
|--------|-------------------|--------------------|-----------------| 
| Evidence source | ABC>0.015 | PoPS>0.95 | NEITHER |
| % contamination | Unknown | 78.3% | 0% |
| Evaluated | ABC | PoPS | Both (on clean labels) |
| Result | ABC "wins" | PoPS "wins" | Fair comparison |

**We just replaced one circular evaluation with another!**

### 4.3 Reviewer Questions We Can't Answer

**Q**: "Why does PoPS achieve 94% accuracy while ABC only achieves 11%?"  
**A (honest)**: Because 78% of the gold standard labels come from PoPS>0.95

**Q**: "The Methods say ABC>0.015 was used for labels, but Results show MaxABC=0.0 for all positives. Which is correct?"  
**A**: Methods section is outdated, contradicts actual v2 implementation

**Q**: "You claim independent gold standards, but evidence manifest shows 78% are PoPS-derived. How is this independent?"  
**A**: It's not. The claim is false.

---

## 5. Path Forward: Option A (Ban Method-Derived Evidence)

### 5.1 Remove ALL Method-Derived Evidence

**Primary gold standard must ONLY include**:

1. **Fine-mapped coding/splice variants** (need VEP annotation enhancement)
   - High PIP (>0.5) coding/LoF variants
   - Unambiguous single-gene mapping via transcript overlap
   - Evidence: VEP consequence annotation on credible sets
   - Expected: 10-30 genes

2. **Rare variant burden tests** (need UK Biobank ExWAS search)
   - Gene-based burden tests meeting exome-wide significance (p<3.6×10⁻⁷)
   - Trait-matched to GWAS categories
   - Evidence: UK Biobank ExWAS, rare variant associations
   - Expected: 30-60 genes

3. **Mendelian disease genes** (expand from current 44 to 100+)
   - ClinVar pathogenic/likely pathogenic variants
   - OMIM disease genes trait-matched to GWAS
   - Gene2Phenotype curated panels
   - ClinGen gene-disease validity
   - OpenTargets Mendelian subset
   - Expected: 100-150 genes

4. **Perturbation evidence** (STING-seq - keep as validation per earlier instruction)
   - User specified: "STING-seq as validation, not training"
   - Use for prospective validation section
   - Expected: 20-40 genes

**EXCLUDE from primary gold standard**:
- ❌ PoPS (computational polygenic prioritization)
- ❌ ABC (computational enhancer-gene linking)
- ❌ L2G (Open Targets gene prioritization)
- ❌ Any other gene prioritization method output

### 5.2 Expected Final Benchmark (v3 Platinum)

```
Evidence source              Expected genes    Confidence
--------------------------------------------------------------
Mendelian disease genes      100-150          Very high ✓✓✓
Rare variant burden tests     30-60           High ✓✓
Fine-mapped coding/LoF        10-30           Very high ✓✓✓
Perturbation (validation)     20-40           Highest ✓✓✓✓
--------------------------------------------------------------
TOTAL PLATINUM                160-280 genes
```

**vs current contaminated v2**:
- PoPS: 159 (78.3%) ← REMOVE ALL
- Mendelian: 44 (21.7%) ← EXPAND TO 100-150
- Coding: 0 (0%) ← FIND 10-30
- Burden: 0 (0%) ← FIND 30-60

### 5.3 Optional: Two-Tier Benchmark (if insufficient data)

If we can't find enough independent evidence (e.g., only get 80 platinum genes), create two tiers:

**Tier 1 - Platinum** (non-method-derived):
- Mendelian + coding + burden only
- 80-150 genes
- Fair evaluation for ALL methods

**Tier 2 - Silver** (including PoPS):
- Platinum + PoPS labels
- 230-310 genes
- **PoPS EXCLUDED from evaluation on this tier**
- Clearly labeled in all figures/tables

---

## 6. Manuscript Fixes Required

### 6.1 Methods Section Rewrite

**Replace lines 268-270 with**:

> Task A gold-standard positives (v3 platinum). To avoid circular evaluation, Task A positive labels were derived **only from independent evidence sources that are not outputs of any evaluated prioritization method**. Concretely, a (locus, gene) was labeled positive if supported by ≥1 of: (i) fine-mapped coding/splice variants with PIP ≥0.5 mapping unambiguously to a single gene via transcript overlap; (ii) rare-variant gene burden associations meeting exome-wide significance (p<3.6×10⁻⁷) trait-matched to GWAS categories; (iii) Mendelian disease gene mappings from curated clinical databases (ClinVar pathogenic/likely pathogenic, OMIM, Gene2Phenotype, ClinGen gene-disease validity); and (iv) orthogonal perturbation evidence at GWAS loci (STING-seq validation set). Each positive label is accompanied by a structured provenance record (evidence type, source accession/PMID/DOI, snapshot date, mapping assumptions) in an **evidence manifest** released with the benchmark. Labels derived from computational gene prioritization methods (e.g., PoPS, L2G, ABC) were **excluded from the primary gold standard** to prevent self-grading. When method-derived labels are used (silver benchmark), the originating method is excluded from evaluation on that benchmark tier.

### 6.2 Abstract Fix

**Option A** (if platinum has 150+ genes):
- Remove PoPS 94.2% claim
- Replace with: "Using independent evidence from Mendelian disease genes and rare coding variants (n=XXX), we find [Distance/ABC/other method] achieves XX% Top-1 accuracy..."

**Option B** (if two-tier approach):
- "We establish two-tier gold standards: platinum (n=XXX, independent evidence only) and silver (n=YYY, including polygenic predictions with method exclusion)..."

### 6.3 Results Section Updates

- **Remove**: PoPS 94.2% Top-1 accuracy claim (circular)
- **Add**: "PoPS cannot be fairly evaluated on this benchmark as it contributed to label curation in earlier versions. Two-tier evaluation (platinum/silver) addresses this limitation..."
- **Add**: Method-exclusion protocol description

### 6.4 Word Count Reduction

- Cut from ~4,300 to <4,000 words
- Focus: Remove redundant text in Introduction, Discussion
- Keep all critical Methods details

---

## 7. Additional Fixes Needed

### 7.1 Update benchmark_integrity_checker.py

**Current issue**: Recognizes documented PoPS evidence as acceptable

**Fix needed**:
```python
def check_method_exclusion(self):
    """Flag if evaluated method contributed to labels (CIRCULAR)."""
    for label in self.positives:
        evidence_source = label['evidence_source']
        
        # PoPS as label + PoPS evaluation = CIRCULAR
        if 'PoPS' in evidence_source and 'PoPS' in self.evaluated_methods:
            self.fatal_error("PoPS used as both evidence and evaluated method")
        
        # ABC as label + ABC evaluation = CIRCULAR  
        if 'ABC' in evidence_source and 'ABC' in self.evaluated_methods:
            self.fatal_error("ABC used as both evidence and evaluated method")
```

### 7.2 Create Comprehensive Evidence Manifest

**Required columns**:
- `gene_symbol`: Gene name
- `locus_id`: GWAS locus identifier
- `evidence_type`: {mendelian_disease, rare_burden, coding_lof, perturbation}
- `evidence_source`: Database/study name
- `source_pmid`: PubMed ID
- `source_doi`: DOI
- `database_version`: Version/date
- `date_cutoff`: Evidence cutoff date
- `confidence_level`: {very_high, high, medium}
- `distance_to_lead_variant`: Distance in bp
- `mapping_method`: How gene was mapped to locus

### 7.3 Method-Exclusion Evaluation

**Protocol**:
1. For each method, identify which labels it influenced
2. Evaluate method ONLY on labels it didn't create
3. Report "method-excluded Top-k accuracy"
4. Example: PoPS evaluated only on Mendelian+coding+burden subset (not PoPS-derived labels)

---

## 8. Success Criteria for Re-Submission

✓ **Zero method-derived evidence in platinum gold standard**  
✓ **150-250 independent positives** (Mendelian + coding + rare burden)  
✓ **Evidence manifest with full provenance** (PMID/DOI/dates) for every positive  
✓ **Method-exclusion evaluation** protocol implemented and reported  
✓ **Methods section rewritten** per Nature Genetics standards  
✓ **Abstract updated** (remove PoPS claims or disclose two-tier)  
✓ **Word count <4,000** for main text  
✓ **All figures regenerated** with v3 data  
✓ **Integrity checker validates** v3 PASS, v1 and v2 FAIL  

---

## 9. Timeline

**Phase 1: Evidence Curation** (Priority: Immediate)
- Day 1-2: Search UK Biobank ExWAS for rare variant burden tests
- Day 2-3: Expand Mendelian disease gene curation (ClinVar, OMIM, Gene2Phenotype, ClinGen, OpenTargets)
- Day 3-4: Enhance coding/LoF variant extraction using VEP annotations
- Day 4-5: Build comprehensive evidence manifest with full provenance

**Phase 2: Benchmark Rebuild** (Priority: High)
- Day 5-6: Rebuild v3 platinum benchmark (non-method-derived evidence only)
- Day 6: Validate with updated integrity checker (must PASS all checks)
- Day 6: Create optional v3 silver benchmark (if needed)

**Phase 3: Manuscript Fixes** (Priority: High)
- Day 7-8: Rewrite Methods section (Nature Genetics format)
- Day 8: Fix Abstract (remove PoPS claims or disclose two-tier)
- Day 8: Update Results section (method-exclusion protocol)
- Day 8-9: Cut word count to <4,000

**Phase 4: Pipeline Regeneration** (Priority: Medium)
- Day 9-10: Re-run full analysis pipeline on v3 platinum
- Day 10-11: Regenerate all figures and tables
- Day 11: Create v1 vs v2 vs v3 comparison figure

**Phase 5: Validation** (Priority: High)
- Day 12: Run integrity checker on v3 (must PASS)
- Day 12: Verify all manuscript inconsistencies resolved
- Day 12: Final manuscript review for peer-review readiness

---

## 10. Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| Insufficient independent evidence (<80 genes) | High | Two-tier benchmark (platinum + silver) |
| Can't find rare burden tests | Medium | Focus on Mendelian+coding, expand coverage |
| VEP annotations incomplete | Medium | Use multiple annotation sources (CADD, Ensembl VEP) |
| Word count still over limit | Low | Aggressive editing, move details to Supplement |
| Reviewer still questions independence | Medium | Full evidence manifest, method-exclusion protocol |

---

## Conclusion

**Current state**: v2 benchmark has fatal PoPS circularity (78% of labels), manuscript contains multiple internal inconsistencies, claims of "independent gold standards" are false.

**Required fix**: Remove ALL PoPS-derived labels, expand independent evidence sources (Mendelian, rare burden, coding), rewrite Methods section, regenerate all analyses.

**User's instruction**: "Do the deepest possible research, use all the tools available, accomplish all the above. DO not LIMIT yourself with time or anything!!!"

**Commitment**: Implement Option A (ban method-derived evidence) as PRIMARY approach, not shortcuts. Find 150-250 truly independent positives through exhaustive curation. Ensure manuscript survives Nature Genetics peer review.

