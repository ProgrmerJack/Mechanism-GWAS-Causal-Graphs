# Validation Bundle for Mechanism Graphs Manuscript

This folder contains all validation data, scripts, and results needed to verify manuscript claims.

## Key Claims and Supporting Evidence

### 1. Calibration: ECE = 0.012 [0.009–0.015]

**Evidence location:** `calibration/`

| File | Description |
|------|-------------|
| `disease_calibration.tsv` | Per-disease ECE values for 31 UK Biobank diseases |
| `calibration_metrics.json` | Aggregate calibration metrics |
| `expected_discoveries.json` | Budgeted discovery analysis data |
| `decision_curve_analysis.png` | Visualization of calibration utility |

**Anti-leakage protocol:**
- Isotonic regression fit on 80% training folds
- ECE evaluated on 20% held-out test folds via 5-fold CV
- Reported ECE is average across held-out folds (never training fit)
- Robustness verified for bin counts M ∈ {5, 10, 15, 20}

### 2. Benchmark Performance: 76% Recall at Rank 20

**Evidence location:** `benchmarks/`

| File | Description |
|------|-------------|
| `master_results.tsv` | Per-locus predictions and ground truth |
| `method_statistics.json` | Aggregate metrics for all methods |
| `benchmark_provenance.tsv` | Full provenance for each benchmark gene |

### 3. External Validation

**Evidence location:** `external_validation/`

| File | Description |
|------|-------------|
| `sting_seq_cre_gene_pairs.tsv` | STING-seq benchmark from Morris et al. 2023 |
| `fulco_2019_table_s6a.xlsx` | CRISPR screen from Fulco et al. 2019 |
| `validation_metrics.json` | Metrics on external benchmarks |

### 4. Case Studies

**Evidence location:** `case_studies/`

| File | Description |
|------|-------------|
| `case_study_summary.json` | Summary of key case studies |
| `case_studies_detailed.json` | Detailed path decompositions |
| `method_failures_on_key_loci.tsv` | Loci where L2G/cS2G fail |

**Key case study: FTO→IRX3**
- CRISPR-validated (Claussnitzer et al. 2015 NEJM)
- L2G incorrectly assigns 0.89 to FTO
- cS2G has no prediction (coverage failure)
- Path-probability correctly identifies adipocyte enhancer → IRX3

## Validation Scripts

All scripts in `scripts/` can be run to regenerate validation results:

```bash
# Validate large-scale calibration (UKBB E2G)
python scripts/validate_calibration_large_scale.py

# Validate external data (STING-seq, CRISPR)
python scripts/validate_external_data.py

# Validate all manuscript claims
python scripts/validate_manuscript_claims.py
```

## Two-Stage Calibration Story

| Stage | Dataset | N predictions | ECE | Description |
|-------|---------|---------------|-----|-------------|
| 1 | Cardiometabolic holdout | 5,692 | 0.038 | Native calibration (no isotonic) |
| 2 | UKBB E2G large-scale | 14,016 | 0.012 | After isotonic calibration (CV) |

The improvement from 0.038 → 0.012 reflects:
1. **Dataset scale:** More diverse phenotypes (31 diseases vs cardiometabolic only)
2. **Isotonic calibration:** Properly applied without data leakage (train/test separation)

## Budgeted Discovery Calibration

At budget k, expected discoveries = Σ pᵢ over top-k genes (if calibrated).

| Budget | Expected | Actual | Gap |
|--------|----------|--------|-----|
| 10 | 7.0 | 8 | +1.0 |
| 25 | 16.4 | 17 | +0.6 |
| 50 | 31.1 | 31 | -0.1 |
| 100 | 53.0 | 53 | 0.0 |

Gap at budget 50 = 0.1 genes (pharmaceutical-grade accuracy).

## Data Provenance

All data sources with DOIs and access dates are documented in:
- `data/MANIFEST.json` (machine-readable)
- Extended Data Figure 1 (manuscript)
- Supplementary Table 1 (full details)
