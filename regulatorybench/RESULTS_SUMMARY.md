# RegulatoryBench v4: Task-Stratified Evaluation Results

## Key Findings Summary

### Task A: GWAS Credible Set → Causal Gene
**Benchmark**: 14,016 pairs across 560 loci (4.1% positive rate)

| Rank | Method | Category | AUROC | Coverage |
|------|--------|----------|-------|----------|
| 1 | Distance (Rank) | Distance | **0.930** | 100% |
| 2 | Distance to Gene Body | Distance | 0.873 | 100% |
| 3 | Distance to Promoter | Distance | 0.855 | 100% |
| 4 | PoPS | GWAS-derived | 0.786 | 95.3% |
| 5 | Fishilevich 2017 | Enhancer | 0.694 | 100% |
| 6 | Connection Strength | Regulatory | 0.604 | 100% |
| 7 | ABC (Max Score) | ABC | 0.599 | 100% |

**Critical Insight**: On GWAS→Gene prediction, **simple distance ranking achieves AUROC=0.930**, vastly outperforming all enhancer-gene linking methods. This suggests that for locus-to-gene mapping at GWAS scale, proximity remains the dominant signal.

---

### Task B: Regulatory Element → Target Gene
**Benchmark**: 19,825 pairs from 3 sources (4.4% positive rate)

| Rank | Method | AUROC | AUPRC | Coverage |
|------|--------|-------|-------|----------|
| 1 | ABC (Fulco 2019) | **0.885** | 0.443 | 25.7% |
| 2 | Distance to TSS | 0.877 | 0.403 | 74.3% |

**Distance-Stratified Analysis** (Distance baseline):
| Distance Range | N Pairs | N Positive | AUROC |
|---------------|---------|------------|-------|
| 0-10 kb | 281 | 163 | 0.672 |
| 10-100 kb | 2,267 | 343 | 0.707 |
| 100-500 kb | 6,635 | 121 | 0.652 |

**Critical Insight**: ABC outperforms distance slightly (0.885 vs 0.877), but this is evaluated on its own training data (Fulco 2019), creating severe leakage.

---

## Leakage Audit Results

### Task B Leakage Issues

| Audit | Status | Affected Pairs | Recommendation |
|-------|--------|----------------|----------------|
| ABC Training Overlap | ⚠️ WARNING | 15,447 / 19,825 (78%) | Exclude K562/Fulco when evaluating ABC |
| Gene Overlap Across Sources | ⚠️ WARNING | 12,759 / 19,825 (64%) | Use leave-one-study-out CV |
| Genomic Proximity | ⚠️ WARNING | 18,287 / 19,825 (92%) | Use chromosome holdout |

### Leave-One-Study-Out Splits (Task B)

| Holdout Study | Train Size | Test Size | Train Pos | Test Pos |
|---------------|------------|-----------|-----------|----------|
| ENCODE_CRISPRi_K562 | 9,469 | 10,356 | 392 | 471 |
| ENCODE_CRISPRi_heldout | 15,447 | 4,378 | 673 | 190 |
| Fulco_2019 | 14,734 | 5,091 | 661 | 202 |

---

## Task Taxonomy Framework

### Task A: GWAS Credible Set → Causal Gene
- **Input**: GWAS fine-mapping credible set
- **Output**: Causal gene at locus
- **Applicable Methods**: L2G, PoPS, FLAMES, proximity, coding variant
- **Inapplicable Methods**: ABC (not designed for GWAS interpretation)

### Task B: Regulatory Element → Target Gene
- **Input**: Defined enhancer region
- **Output**: Target gene(s)
- **Applicable Methods**: ABC, rE2G, Hi-C, eQTL colocalization
- **Inapplicable Methods**: L2G, FLAMES (require GWAS inputs)

### Task C: Regulatory Variant → Affected Gene
- **Input**: Specific variant position
- **Output**: Genes whose expression is affected
- **Applicable Methods**: eQTL, VEP, DeepSEA
- **Note**: Task C benchmark (MPRA) currently not loaded

---

## Implications for Nature Genetics Analysis

### Main Finding
The field has been **silently conflating Task A and Task B** in benchmarking comparisons. Methods designed for enhancer-gene linking (ABC) are being evaluated on GWAS→gene tasks where simple distance dominates, and vice versa.

### Proposed Display Items (6 for NG Analysis)

1. **Figure 1**: Task taxonomy schematic (Task A/B/C definitions)
2. **Figure 2**: Task A results - Distance dominance in GWAS→Gene
3. **Figure 3**: Task B results - ABC vs Distance in E2G
4. **Figure 4**: Leakage audit heatmap
5. **Table 1**: Method comparison by task type
6. **Table 2**: Leakage-controlled evaluation results

---

## Files Generated

- `regulatorybench/task_taxonomy.py` - Task definitions and evaluation contract
- `regulatorybench/split_benchmarks.py` - Benchmark splitter
- `regulatorybench/evaluate_task_a.py` - Task A evaluation
- `regulatorybench/evaluate_task_b.py` - Task B evaluation
- `regulatorybench/leakage_audit.py` - Leakage detection framework
- `regulatorybench/benchmarks/task_a_gwas_to_gene.parquet` - 14,016 pairs
- `regulatorybench/benchmarks/task_b_enhancer_to_gene.parquet` - 19,825 pairs
- `regulatorybench/benchmarks/task_a_evaluation_results.json`
- `regulatorybench/benchmarks/task_b_baseline_results.json`
- `regulatorybench/benchmarks/task_a_leakage_audit.json`
- `regulatorybench/benchmarks/task_b_leakage_audit.json`

---

## Next Steps

1. ✅ Task taxonomy created
2. ✅ Benchmarks split by task type
3. ✅ Baselines evaluated
4. ✅ Leakage audit completed
5. ⬜ Run leave-one-study-out evaluation (ABC on held-out data)
6. ⬜ Create publication figures
7. ⬜ Write Nature Genetics Analysis manuscript (~3,000 words)
8. ⬜ Prepare Zenodo package with DOI
