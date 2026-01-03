# Nature Genetics Professional Figures - Documentation

**Generated:** Professional publication-quality figures for Nature Genetics submission  
**Location:** `figures/nature_genetics_v2/`  
**Script:** `scripts/generate_nature_genetics_professional_figures.py`

---

## Technical Specifications

All figures comply with Nature Genetics guidelines:

| Parameter | Specification |
|-----------|---------------|
| **Resolution** | 600 DPI (raster) |
| **Formats** | PDF (vector), PNG, TIFF |
| **Single Column** | 89 mm (~3.5 in) |
| **Double Column** | 183 mm (~7.2 in) |
| **Fonts** | Sans-serif (Arial), 7-12 pt |
| **Color Palette** | Okabe-Ito (colorblind-safe) |
| **Font Embedding** | Type 42 (PDF) |

### Colorblind-Safe Palette (Okabe-Ito)

| Color | Hex | Usage |
|-------|-----|-------|
| Blue | `#0072B2` | Primary - Mechanism Graphs |
| Vermilion | `#D55E00` | Secondary - L2G/competitors |
| Orange | `#E69F00` | Tertiary data |
| Sky Blue | `#56B4E9` | Light accents |
| Bluish Green | `#009E73` | Success/positive |
| Yellow | `#F0E442` | Highlights |
| Reddish Purple | `#CC79A7` | Alternative method |
| Gray | `#999999` | Neutral/baseline |

---

## Figure Inventory

### Main Figures (4)

#### Figure 1: Calibration Overview
**File:** `Figure_1_Calibration_Overview.{pdf,png,tiff}`  
**Size:** Double column (183 mm)

**Panels:**
- **a)** Reliability diagram comparing Mechanism Graphs (ECE=0.012) vs L2G (ECE=0.18)
- **b)** ECE comparison bar chart across all methods with 15× improvement annotation
- **c)** Expected vs Actual discoveries (budget validation)
- **d)** Per-disease calibration sorted waterfall (31 diseases)

**Key Message:** Decision-grade calibration enables rational resource allocation.

---

#### Figure 2: Stress Test Validation
**File:** `Figure_2_Stress_Test.{pdf,png,tiff}`  
**Size:** Double column

**Panels:**
- **a)** Leave-one-family-out ECE comparison (training vs held-out)
- **b)** Transfer ratio visualization (8/8 families pass)

**Key Message:** Robust generalization across disease domains.

---

#### Figure 3: Flagship Case Studies
**File:** `Figure_3_Case_Studies.{pdf,png,tiff}`  
**Size:** Double column

**Panels:**
- **a)** FTO→IRX3: Resolving a decade of misdirection
- **b)** APOE: Tissue-specific pathway decomposition (AD vs CAD)
- **c)** TCF7L2: Pleiotropic mechanism resolution

**Key Message:** Mechanism graphs provide interpretable, actionable insights.

---

#### Figure 4: Benchmark Performance Comparison
**File:** `Figure_4_Benchmark_Comparison.{pdf,png,tiff}`  
**Size:** Double column

**Panels:**
- **a)** Recall@K curves for all methods
- **b)** Top-1 accuracy ranking (horizontal bars, sorted)
- **c)** Performance matrix heatmap
- **d)** Mean Reciprocal Rank (MRR) comparison

**Key Message:** Superior performance across multiple metrics.

---

### Extended Data Figures (5)

#### Extended Data Figure 1: Detailed Reliability Analysis
**File:** `Extended_Data_Figure_1_Reliability.{pdf,png,tiff}`

**Panels:**
- **a)** ECE distribution across all diseases
- **b)** ECE vs base rate (colored by sample size)
- **c)** Confidence intervals for each disease
- **d)** ECE vs sample size

---

#### Extended Data Figure 2: Decision Curve Analysis
**File:** `Extended_Data_Figure_2_Decision_Curve.{pdf,png,tiff}`

**Panels:**
- **a)** Calibration error by budget
- **b)** Expected vs actual discoveries (line plot)

---

#### Extended Data Figure 3: Component Ablation
**File:** `Extended_Data_Figure_3_Ablation.{pdf,png,tiff}`

**Panels:**
- **a)** ECE by ablation condition (Full, -L2G, -ABC, -Coloc, Distance)
- **b)** ECE improvement attribution (pie chart)

---

#### Extended Data Figure 4: CRISPRi Validation
**File:** `Extended_Data_Figure_4_CRISPR.{pdf,png,tiff}`

**Panels:**
- **a)** Precision-recall curves on 863 CRISPRi pairs
- **b)** Performance by enhancer-gene distance

---

#### Extended Data Figure 5: Benchmark Statistical Details
**File:** `Extended_Data_Figure_5_Benchmark_Details.{pdf,png,tiff}`

**Panels:**
- **a)** Top-1 accuracy with 95% confidence intervals
- **b)** Performance stratified by gene evidence tier
- **c)** Statistical significance vs other methods

---

### Supplementary Figures (2)

#### Conceptual Framework
**File:** `Supplementary_Conceptual_Framework.{pdf,png,tiff}`

Visual representation of the mechanism graph framework:
- SORT1 example with variant → enhancer → gene pathway
- Probability propagation formulas
- Five-stage pipeline overview
- Comparison to traditional approaches

---

#### Graphical Abstract
**File:** `Graphical_Abstract.{pdf,png,tiff}`

Summary visualization for journal landing page:
- Key result: ECE = 0.012 (15× better)
- Three key features
- FTO case study highlight

---

## Data Sources

All figures are data-driven from validated results files:

| Figure | Data Source |
|--------|-------------|
| Fig 1a-d | `results/calibration_validation/cv_ece_results.json`, `disease_calibration.tsv` |
| Fig 1c | `results/decision_curve/expected_discoveries.json` |
| Fig 2 | `results/stress_test/leave_family_out_results.json` |
| Fig 3 | `results/case_studies/case_studies_detailed.json` |
| Fig 4 | `results/baselines/post2021_comparison_metrics.tsv` |

---

## Regeneration Instructions

To regenerate all figures:

```bash
python scripts/generate_nature_genetics_professional_figures.py
```

Output will be saved to `figures/nature_genetics_v2/` in PDF, PNG, and TIFF formats.

---

## Quality Checklist

- [x] 600 DPI resolution
- [x] Colorblind-safe palette (Okabe-Ito)
- [x] Sans-serif fonts (Arial)
- [x] Proper panel labels (lowercase bold)
- [x] Vector PDF for print
- [x] Raster PNG/TIFF for review
- [x] Nature Genetics dimensions
- [x] Data-driven from validated results
- [x] Minimal spines (top/right removed)
- [x] Consistent styling across figures

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 2.0 | 2025 | Complete rewrite with professional styling |

