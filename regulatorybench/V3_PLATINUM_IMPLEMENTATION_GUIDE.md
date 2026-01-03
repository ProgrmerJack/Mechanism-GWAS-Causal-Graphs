# RegulatoryBench v3 Platinum Gold Standard - Implementation Guide

**Status**: Ready for execution  
**Date**: 2025-01  
**Purpose**: Fix PoPS circularity by removing ALL method-derived evidence

---

## Problem Summary

The v2 benchmark replaced ABC circularity with **PoPS circularity**:
- **Evidence contamination**: 159/203 positives (78%) derived from PoPS>0.95 computational predictions
- **Circular evaluation**: PoPS evaluated on benchmark using PoPS-derived labels → reports 94.2% Top-1 accuracy (self-fulfilling)
- **Manuscript inconsistencies**: 
  - Methods (line 268) still references ABC predictions (obsolete)
  - Abstract claims "independent gold standards" but 78% from PoPS
  - Results reports PoPS 94.2% without disclosing circularity
- **Evidence manifest**: 100% of documented evidence is 'PoPS_high_confidence' (218 rows)
- **Peer-review risk**: FATAL - primary result invalidated by circularity

---

## Solution: v3 Platinum Gold Standard

**Principle**: Ban ALL computational prioritization methods from labels

**Evidence sources allowed** (truly independent):
1. **Mendelian disease genes** (ClinVar, OMIM, Gene2Phenotype, ClinGen, OpenTargets)
2. **Fine-mapped coding/LoF variants** (VEP annotations, PIP>0.5)
3. **Rare variant burden tests** (UK Biobank ExWAS, p<3.6×10⁻⁷)
4. **Perturbation data** (STING-seq - validation set only)

**Evidence sources BANNED**:
- ❌ PoPS scores (computational gene prioritization)
- ❌ ABC scores (computational regulatory activity)
- ❌ L2G scores (computational locus-to-gene mapping)
- ❌ Any other machine learning-based prioritization

**Target**: 150-250 independent positives
- Mendelian: 100-150 genes (expand from current 44)
- Rare burden: 30-60 genes
- Coding variants: 10-30 genes

---

## Implementation Steps

### Phase 1: Evidence Curation (Scripts Created ✓)

#### 1.1 Mendelian Disease Genes
**Script**: `expand_mendelian_curation.py` ✓ CREATED

**Data sources**:
- ClinVar pathogenic/likely pathogenic variants
- OMIM gene-disease associations
- Gene2Phenotype curated panels (DDG2P, cardiac, cancer, etc.)
- ClinGen gene-disease validity curations
- OpenTargets Genetics L2G Mendelian subset

**Curation criteria**:
- Disease must match GWAS trait by ontology
- Gene within locus boundaries (<500kb preferred, <1Mb max)
- Confidence scoring based on multi-source agreement
- Full provenance: PMID, DOI, database version, date

**Expected output**: `mendelian_genes_curated.tsv` (~100-150 genes)

**Execution**:
```bash
cd Mechanism-GWAS-Causal-Graphs/regulatorybench
python expand_mendelian_curation.py --benchmark-dir benchmarks --output mendelian_genes_curated.tsv
```

**Data acquisition needed**:
1. ClinVar: Download from NCBI FTP (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/)
2. OMIM: Requires license - download genemap2.txt
3. Gene2Phenotype: https://www.ebi.ac.uk/gene2phenotype/downloads
4. ClinGen: https://search.clinicalgenome.org/kb/gene-validity/download
5. OpenTargets: https://genetics.opentargets.org/ API

---

#### 1.2 Coding/LoF Variants
**Script**: `extract_coding_vep_enhanced.py` ✓ CREATED

**Data source**: FLAMES Annotation_data (321,501+ VEP-annotated parquet files)

**Extraction criteria**:
- Fine-mapped variants in credible sets (PIP>0.5 for LoF, PIP>0.8 for missense)
- VEP consequences: stop_gained, frameshift, splice_donor/acceptor, missense
- Unambiguous gene mapping (single gene only)
- CADD score filtering for missense variants

**Expected output**: `coding_lof_genes.tsv` (~10-30 genes)

**Execution**:
```bash
cd Mechanism-GWAS-Causal-Graphs/regulatorybench
python extract_coding_vep_enhanced.py \
    --annotations-dir data/external/flames/Annotation_data \
    --benchmark-dir benchmarks \
    --pip-threshold-lof 0.5 \
    --pip-threshold-missense 0.8 \
    --output coding_lof_genes.tsv
```

**Prerequisites**: Fine-mapping results (SuSiE credible sets)

---

#### 1.3 Rare Variant Burden Tests
**Script**: `search_exwas_burden.py` ✓ CREATED

**Data sources**:
- UK Biobank ExWAS gene-based burden tests
- GeneBass database (https://app.genebass.org/)
- Published ExWAS catalogs (DeBoever 2018, Backman 2021, Karczewski 2020)

**Search criteria**:
- Exome-wide significance (p<3.6×10⁻⁷)
- Phenotype matches GWAS trait by ontology
- Gene within locus boundaries

**Expected output**: `rare_burden_genes.tsv` (~30-60 genes)

**Execution**:
```bash
cd Mechanism-GWAS-Causal-Graphs/regulatorybench
python search_exwas_burden.py \
    --benchmark-dir benchmarks \
    --p-threshold 3.6e-7 \
    --sources uk_biobank_exwas genebass published_catalogs \
    --output rare_burden_genes.tsv
```

**Data acquisition needed**:
1. UK Biobank ExWAS: Download from published sources or request from UKBB
2. GeneBass: Download from https://app.genebass.org/downloads
3. Manual curation of published ExWAS results (PMIDs: 29785011, 34012112, 32461654)

---

### Phase 2: Benchmark Reconstruction

#### 2.1 Update rebuild_task_a_gold_standard.py

**Changes required**:
1. **REMOVE**: `extract_pops_high_confidence()` method entirely
2. **UPDATE**: `load_mendelian_disease_genes()` to use `mendelian_genes_curated.tsv`
3. **ADD**: `load_rare_burden_genes()` to load `rare_burden_genes.tsv`
4. **ADD**: `load_coding_variant_genes()` to load `coding_lof_genes.tsv`
5. **UPDATE**: Evidence manifest generation with full provenance

**Expected output**: 
- `task_a_gwas_to_gene_v3_platinum.parquet` (150-250 positives, ZERO PoPS/ABC/L2G)
- `task_a_v3_evidence_manifest.tsv` (full provenance for all positives)

**Validation**: Run `benchmark_integrity_checker.py` → should PASS all checks

---

#### 2.2 Create Optional v3 Silver Benchmark

**Purpose**: Keep PoPS labels for exploratory analysis, but clearly marked

**Output**: `task_a_gwas_to_gene_v3_silver.parquet`

**Warnings in metadata**:
- "Silver labels include PoPS computational predictions"
- "PoPS should NOT be evaluated on this benchmark"
- "Use for exploratory analysis only, NOT primary publication"

**Use cases**: Secondary validation, method development, exploratory analyses

---

### Phase 3: Evaluation & Validation

#### 3.1 Update benchmark_integrity_checker.py

**New checks**:
1. `check_method_exclusion()`: Detect if evaluated method appears in evidence_source
2. Update `check_label_independence()`: Flag computational prioritization in labels
3. Validate provenance: Every positive must have PMID/DOI/database version

**Validation targets**:
- v3 platinum: Should PASS all checks
- v2: Should FAIL with "PoPS circularity detected"

---

#### 3.2 Create method_exclusion_evaluator.py

**Purpose**: Evaluate each method only on labels it didn't create

**Protocol**:
- ABC: Exclude labels where MaxABC>0.015
- PoPS: Exclude labels where POPS.Score>0.95 (or use Mendelian+coding subset only)
- L2G: Exclude labels from Gene2Phenotype (L2G training data)
- Distance: No exclusions (non-parametric baseline)

**Output**: Performance metrics on truly independent subset

---

#### 3.3 Re-run Analysis Pipeline

**Steps**:
1. Update `evaluate_task_a.py` to use v3 platinum benchmark
2. Run `locus_level_evaluation.py` for standard metrics
3. Run `method_exclusion_evaluator.py` for method-excluded metrics
4. Compare v1 vs v2 vs v3 performance

**Expected changes**:
- Distance/ABC may improve relative to PoPS (no longer circular)
- Smaller sample size (150-250 vs 203) but cleaner evaluation
- Method-excluded results show true generalization ability

---

### Phase 4: Manuscript Rewrite

#### 4.1 Rewrite Methods Section

**Template** (from user's instruction):
```
Task A gold-standard positives (v3). To avoid circular evaluation, Task A positive 
labels were derived **only from independent evidence sources that are not outputs 
of any evaluated prioritization method**. Concretely, a (locus, gene) was labeled 
positive if supported by ≥1 of: (i) fine-mapped coding/splice variants with PIP 
≥0.5 mapping unambiguously to a single gene via VEP annotation; (ii) rare-variant 
gene burden associations meeting exome-wide significance (p<3.6×10⁻⁷) in UK Biobank 
ExWAS; (iii) Mendelian disease gene mappings from curated clinical databases 
(ClinVar pathogenic variants, OMIM, Gene2Phenotype, ClinGen gene-disease validity) 
matched to GWAS trait ontology; and (iv) orthogonal perturbation evidence at GWAS 
loci (STING-seq validation set). Each positive label is accompanied by a structured 
provenance record (evidence type, source accession/PMID/DOI, database version, 
snapshot date, mapping assumptions) in an **evidence manifest** released with the 
benchmark. Labels derived from computational gene prioritization methods (e.g., 
PoPS, L2G, ABC) were **excluded from the primary gold standard** to prevent 
self-grading.
```

**Actions**:
- Remove line 268 reference to "ABC model predictions >0.015" (obsolete)
- Add evidence manifest description
- Document method-exclusion evaluation protocol
- Explain v1→v2→v3 evolution

---

#### 4.2 Update Abstract

**Options**:
- **A (recommended)**: Focus on method-excluded results on v3 platinum
  - Remove PoPS 94.2% claim entirely
  - Report: "Using truly independent evidence (Mendelian, coding, rare burden, n=150-250), distance achieves X% recall, ABC achieves Y%..."
  
- **B (fallback)**: Disclose two-tier approach
  - "We present platinum (independent evidence, n=150-250) and silver (including PoPS, n=203) benchmarks..."
  - Clearly state PoPS excluded from primary evaluation

**Constraint**: Maintain ~148 words (current length)

---

#### 4.3 Update Results Section

**Changes**:
1. Remove/qualify PoPS 94.2% Top-1 accuracy (circular result)
2. Add explanation: "PoPS cannot be fairly evaluated on benchmarks using PoPS predictions as evidence. We therefore evaluate PoPS on Mendelian+coding subset only, or exclude it from primary leaderboard."
3. Present v3 platinum results (Distance vs ABC comparison)
4. Add method-exclusion results
5. Cut word count to reach <4,000 for Nature Genetics

---

#### 4.4 Regenerate Figures & Tables

**Updates required**:
- **Figure 2a**: Method comparison leaderboard (v3 platinum results)
- **Figure 2b**: Regime-stratified analysis (re-run on v3)
- **Table 1**: Update positive counts (v1: 569, v2: 203, v3: 150-250)
- **Supplementary Table S1**: Evidence manifest with full provenance
- **NEW Figure**: Venn diagram showing v1 vs v2 vs v3 overlap, evidence sources

---

## Success Criteria for Nature Genetics

- [x] Evidence curation scripts created (mendelian, coding, burden)
- [ ] Execute all curation scripts → generate evidence files
- [ ] Rebuild v3 platinum: 150-250 independent positives, ZERO method-derived
- [ ] Evidence manifest: Full provenance (PMID/DOI/version/date) for ALL positives
- [ ] Integrity checker: v3 PASS, v2 FAIL (PoPS circularity detected)
- [ ] Method-exclusion evaluation implemented
- [ ] Methods section rewritten (user template, ABC language removed)
- [ ] Abstract updated (PoPS circular claim removed)
- [ ] Results section fixed (PoPS excluded or evaluated on independent subset)
- [ ] Word count reduced to <4,000
- [ ] All figures/tables regenerated with v3 data
- [ ] Comprehensive v1→v2→v3 documentation

---

## Timeline Estimate

**Phase 1** (Evidence Curation): 2-3 days
- Data acquisition (ClinVar, OMIM, ExWAS): 1 day
- Script execution and validation: 1 day
- Quality control and deduplication: 0.5 day

**Phase 2** (Benchmark Reconstruction): 1 day
- Update rebuild script: 0.5 day
- Generate v3 platinum and silver: 0.25 day
- Run integrity checks: 0.25 day

**Phase 3** (Evaluation): 1 day
- Update evaluation scripts: 0.5 day
- Re-run full analysis pipeline: 0.5 day

**Phase 4** (Manuscript): 1-2 days
- Rewrite Methods/Abstract/Results: 1 day
- Regenerate figures/tables: 0.5 day
- Final validation: 0.5 day

**Total**: 5-7 days for complete v3 platinum implementation

---

## Data Acquisition Checklist

### Required Downloads

1. **ClinVar**
   - Source: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/
   - File: `variant_summary.txt.gz`
   - Size: ~2GB
   - Update: Monthly

2. **OMIM** (requires license)
   - Source: https://omim.org/downloads/
   - File: `genemap2.txt`
   - Requires institutional access

3. **Gene2Phenotype**
   - Source: https://www.ebi.ac.uk/gene2phenotype/downloads
   - Files: DDG2P, cardiac, cancer, skin, eye panels
   - Format: CSV

4. **ClinGen**
   - Source: https://search.clinicalgenome.org/kb/gene-validity/download
   - File: Gene-disease validity curations
   - Format: TSV

5. **OpenTargets Genetics**
   - Source: https://genetics.opentargets.org/
   - API or bulk download
   - Filter: L2G Mendelian associations

6. **UK Biobank ExWAS**
   - Source: Published papers or UKBB application
   - Papers: DeBoever 2018 (PMID: 29785011), Backman 2021 (PMID: 34012112)
   - Requires data access agreement

7. **GeneBass**
   - Source: https://app.genebass.org/downloads
   - File: Gene-based burden test results
   - Free download

---

## Next Immediate Actions

1. **Download required data sources** (see checklist above)
2. **Execute expand_mendelian_curation.py** to generate Mendelian gene list
3. **Execute extract_coding_vep_enhanced.py** to extract coding/LoF genes
4. **Execute search_exwas_burden.py** to find rare burden associations
5. **Update rebuild_task_a_gold_standard.py** with new evidence loaders
6. **Generate v3 platinum benchmark** and validate with integrity checker
7. **Re-run evaluation pipeline** on clean benchmark
8. **Rewrite manuscript** per templates above
9. **Regenerate all figures and tables** with v3 results

---

## Contact & Questions

For questions about implementation:
- Review PEER_REVIEW_CIRCULARITY_ANALYSIS.md for detailed problem analysis
- Check individual script docstrings for usage examples
- Validate outputs with benchmark_integrity_checker.py

**Remember**: The goal is "deepest possible research, no limits" - exhaustive multi-source curation to achieve truly independent gold standard for peer-review survival.
