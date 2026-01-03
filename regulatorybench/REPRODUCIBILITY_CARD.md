# RegulatoryBench Reproducibility Card

**Nature Genetics Submission Compliance**

---

## Overview

This reproducibility card documents all information required for reproducing the 
RegulatoryBench analysis as per Nature Genetics' reporting requirements.

---

## 1. Data Availability

### Primary Benchmarks

| Dataset | Source | Pairs | Positives | DOI/Access |
|---------|--------|-------|-----------|------------|
| Task A: GWAS→Gene | UK Biobank + ABC validation | 14,016 | 569 (4.1%) | [Zenodo DOI] |
| Task B: Enhancer→Gene | ENCODE + Fulco CRISPRi | 19,825 | 863 (4.4%) | [Zenodo DOI] |

### External Validation

| Resource | Version | Access | Usage |
|----------|---------|--------|-------|
| Open Targets | v25.06 | https://platform.opentargets.org | Drug target genes |
| gnomAD | v4.0 | https://gnomad.broadinstitute.org | Constraint scores |
| GTEx | v8 | https://gtexportal.org | eQTL data |
| ENCODE | Phase 4 | https://encodeproject.org | CRISPRi data |

### Deposited Data

All processed benchmarks and evaluation results are deposited at:
- **Zenodo:** [DOI to be assigned upon acceptance]
- **GitHub:** https://github.com/[username]/regulatorybench

---

## 2. Code Availability

### Repository Structure

```
regulatorybench/
├── run_all.py                    # One-command reproduction
├── enhanced_evaluation.py        # Comprehensive evaluation
├── leakage_audit.py             # Leakage detection + LOSO
├── generate_enhanced_figures.py  # Publication figures
├── benchmarks/                   # Processed data
│   ├── task_a_gwas_to_gene.parquet
│   ├── task_b_enhancer_to_gene.parquet
│   └── *.json (results)
└── figures/                      # Generated figures
```

### Software Requirements

```
Python >= 3.9
pandas >= 2.0
numpy >= 1.24
scikit-learn >= 1.3
scipy >= 1.11
matplotlib >= 3.8
seaborn >= 0.13
```

### Reproduction Commands

```bash
# Clone repository
git clone https://github.com/[username]/regulatorybench.git
cd regulatorybench

# Install dependencies
pip install -r requirements.txt

# Run complete analysis (6 steps)
python run_all.py

# Expected output:
# ✓ Step 1: Task A evaluation
# ✓ Step 2: Task B evaluation  
# ✓ Step 3: Leakage audit
# ✓ Step 4: Enhanced evaluation
# ✓ Step 5: Generate figures
# ✓ Step 6: Enhanced figures
```

---

## 3. Computational Environment

### Hardware Used

- CPU: [To be specified]
- RAM: 32 GB minimum recommended
- Storage: ~5 GB for benchmarks + results

### Expected Runtime

| Step | Time (approx) |
|------|---------------|
| Data loading | 30 seconds |
| Task A evaluation | 2 minutes |
| Task B evaluation | 3 minutes |
| Bootstrap CIs (1000 iter) | 10 minutes |
| Figure generation | 2 minutes |
| **Total** | **~20 minutes** |

---

## 4. Statistical Analysis

### Primary Metrics

| Metric | Implementation | Bootstrap |
|--------|---------------|-----------|
| AUROC | sklearn.roc_auc_score | 1000 iterations |
| AUPRC | sklearn.average_precision_score | 1000 iterations |
| Top-1 Accuracy | Custom per-locus | N/A |
| Odds Ratio | Fisher's exact | Woolf CI |

### Confidence Intervals

- Method: Percentile bootstrap
- Iterations: 1000
- Confidence level: 95%
- Seed: 42 (reproducible)

### Multiple Testing

- No formal correction applied to main results
- Exploratory subgroup analyses clearly labeled
- Pre-specified primary outcomes: AUROC, AUPRC

---

## 5. Key Results Summary

### Task A: GWAS→Gene (n=14,016 pairs)

| Method | AUROC [95% CI] | AUPRC | p vs Distance |
|--------|---------------|-------|---------------|
| Distance (Rank) | 0.930 [0.918-0.941] | 0.52 | — |
| Distance to Body | 0.873 [0.858-0.888] | 0.38 | <0.001 |
| PoPS | 0.786 [0.768-0.804] | 0.29 | <0.001 |
| ABC (Max) | 0.599 [0.574-0.624] | 0.09 | <0.001 |

### Drug Target Validation

| Method | Odds Ratio [95% CI] | p-value |
|--------|-------------------|---------|
| PoPS | 28.28 [12.4-64.3] | <0.001 |
| Distance | 8.62 [4.2-17.6] | <0.001 |
| ABC Prediction | 0.64 [0.2-1.9] | 0.42 (ns) |

### External Validation (Ji et al. 2025)

| Method | OR (Drug Approval) | 95% CI | Source |
|--------|-------------------|--------|--------|
| Nearest gene | 3.08 | 2.25-4.11 | medRxiv 2025.09.23.25336370 |
| L2G | 3.14 | 2.31-4.28 | medRxiv 2025.09.23.25336370 |
| eQTL coloc | 1.61 | 0.92-2.83 (ns) | medRxiv 2025.09.23.25336370 |

---

## 6. Leakage Audit Results

### Task B Contamination

| Leakage Type | Affected Pairs | Percentage |
|--------------|----------------|------------|
| ABC training overlap | 15,447 | 78% |
| Gene overlap across sources | 12,759 | 64% |
| Genomic proximity (<100kb) | 18,287 | 92% |

### Recommended Held-Out Splits

| Split Name | Train Size | Test Size | Use For |
|------------|------------|-----------|---------|
| ENCODE K562 holdout | 9,469 | 10,356 | ABC evaluation |
| Fulco excluded | 14,734 | 5,091 | Clean baseline |
| Chromosome 8,9 holdout | ~16,000 | ~3,800 | Independence check |

---

## 7. Limitations and Caveats

1. **Ground truth source:** Task A uses ABC validation labels, which may favor distance
2. **Cell type coverage:** Task B limited to available CRISPRi cell types
3. **Rare variant effects:** Not captured in common disease GWAS
4. **Tissue specificity:** eQTL methods evaluated on GTEx v8 only

---

## 8. Version Control

| Component | Version/Date | Hash/DOI |
|-----------|--------------|----------|
| Code repository | v1.0 | [commit hash] |
| Task A benchmark | 2025-01 | [checksum] |
| Task B benchmark | 2025-01 | [checksum] |
| Zenodo deposit | [Date] | [DOI] |

---

## 9. Contact

For questions or issues reproducing results:
- GitHub Issues: [repository URL]
- Email: [corresponding author email]

---

## 10. Checklist for Reviewers

- [ ] Data accessible at Zenodo
- [ ] Code runs with documented dependencies
- [ ] Results match reported values
- [ ] Figures reproducible from scripts
- [ ] Statistical methods appropriate
- [ ] Leakage controls adequate
- [ ] External validation confirms findings
