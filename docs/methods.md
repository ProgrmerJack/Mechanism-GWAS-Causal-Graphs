# Methods Documentation

## Overview

This document describes the econometric methods used in the e-procurement transparency analysis. We employ two main identification strategies:

1. **Regression Discontinuity Design (RDD)** - exploiting procurement value thresholds
2. **Difference-in-Differences (DiD)** - exploiting staggered reform adoption

All methods follow best practices from the econometrics literature on causal inference.

---

## 1. Data Sources

### 1.1 Procurement Data

| Source | Country | Coverage | Fields |
|--------|---------|----------|--------|
| ProZorro API | Ukraine | 2016-present | OCDS v1.1 compliant |
| SECOP II | Colombia | 2017-present | OCDS v1.1 compliant |
| Contracts Finder | UK | 2015-present | Partial OCDS |

### 1.2 Innovation Data

| Source | Coverage | Use Case |
|--------|----------|----------|
| Lens API | Global | Primary cross-country innovation linkage |
| EPO OPS | European | Patent family resolution |
| PatentsView | US only | Triangulation for US patents |

See `src/patent_linkage.py` for implementation details.

### 1.3 Auxiliary Data

- **TI Corruption Perceptions Index**: Country-level corruption measure
- **World Bank WGI**: Governance quality indicators
- **OECD**: Public procurement statistics

---

## 2. Regression Discontinuity Design

### 2.1 Institutional Background

Procurement regulations impose transparency requirements that vary with contract value. Thresholds create quasi-experimental variation:

| Jurisdiction | Threshold | Requirement |
|--------------|-----------|-------------|
| EU | €139,000 (supplies) | Full OJEU publication |
| EU | €5.35M (works) | Full OJEU publication |
| Ukraine | UAH 200,000 | ProZorro publication |
| Colombia | Various SMMLV | SECOP II requirements |

### 2.2 Estimation Strategy

We estimate the local average treatment effect at the threshold using local polynomial regression (Imbens & Kalyanaraman, 2012):

$$\tau_{RDD} = \lim_{x \downarrow c} \mathbb{E}[Y_i | X_i = x] - \lim_{x \uparrow c} \mathbb{E}[Y_i | X_i = x]$$

**Implementation:**

```python
# Local linear regression specification
Y_i = α + τ·D_i + β₁(X_i - c) + β₂·D_i·(X_i - c) + ε_i
```

where:
- $Y_i$: Outcome (competition, price, quality)
- $D_i = \mathbf{1}(X_i \geq c)$: Treatment indicator
- $X_i$: Contract value (running variable)
- $c$: Threshold cutoff

### 2.3 Bandwidth Selection

We use three bandwidth selection methods:

1. **IK Optimal** (Imbens & Kalyanaraman, 2012):
   - Minimizes asymptotic mean squared error
   - Default for main specification

2. **CCT Robust** (Calonico, Cattaneo & Titiunik, 2014):
   - Bias-corrected with robust confidence intervals
   - Used for inference

3. **CV Cross-Validation**:
   - Leave-one-out cross-validation
   - Used for sensitivity analysis

**Bandwidth Sweep:**
We report estimates across $h \in [0.5h^*, 2h^*]$ where $h^*$ is the optimal bandwidth.

### 2.4 Robustness Checks

We implement a comprehensive robustness battery (see `src/rdd_robustness.py`):

| Test | Purpose | Reference |
|------|---------|-----------|
| McCrary density | Manipulation check | McCrary (2008) |
| Donut-hole RDD | Exclude observations near cutoff | Barreca et al. (2011) |
| Placebo cutoffs | False positive check | Imbens & Lemieux (2008) |
| Covariate smoothness | Balance verification | Lee & Lemieux (2010) |
| Bandwidth sweep | Sensitivity to bandwidth choice | IK (2012) |

**McCrary Density Test:**

Tests for discontinuity in the density of the running variable:
$$H_0: \lim_{x \downarrow c} f(x) = \lim_{x \uparrow c} f(x)$$

Using local polynomial density estimation (Cattaneo, Jansson & Ma, 2020).

**Donut-Hole Estimation:**

Excludes observations within $\pm \delta$ of the cutoff to address potential manipulation:
$$\hat{\tau}_{donut}(\delta) = \hat{\tau}_{RDD} \text{ using } |X_i - c| > \delta$$

### 2.5 Inference

Standard errors are computed using:
1. Heteroskedasticity-robust (HC1) standard errors
2. Clustering at the contracting authority level
3. Bootstrap (1,000 iterations) for small samples

---

## 3. Difference-in-Differences

### 3.1 Setting

We exploit staggered adoption of e-procurement reforms:

| Country | Reform | Timing |
|---------|--------|--------|
| Ukraine | ProZorro launch | 2016 |
| Colombia | SECOP II | 2017 |
| UK | Contracts Finder overhaul | 2015 |

### 3.2 Estimation Strategy

**Two-Way Fixed Effects (TWFE):**

$$Y_{it} = \alpha_i + \gamma_t + \tau D_{it} + X_{it}'\beta + \epsilon_{it}$$

**Problem:** TWFE can be biased with heterogeneous treatment effects and staggered timing (Goodman-Bacon, 2021).

**Solution:** Callaway-Sant'Anna (2021) estimator:

$$ATT(g, t) = \mathbb{E}[Y_{it}(g) - Y_{it}(\infty) | G_i = g]$$

where $G_i$ is the first treatment period for unit $i$.

### 3.3 Event Study Specification

We estimate dynamic treatment effects:

$$Y_{it} = \alpha_i + \gamma_t + \sum_{k \neq -1} \tau_k \cdot \mathbf{1}(t - G_i = k) + \epsilon_{it}$$

**Pre-trend Test:**
$$H_0: \tau_{-K} = \tau_{-K+1} = ... = \tau_{-2} = 0$$

Tested via joint F-test with Roth (2022) correction for pre-testing.

### 3.4 Robustness Checks

See `src/did_robustness.py` for implementation:

| Test | Purpose | Reference |
|------|---------|-----------|
| Bacon decomposition | Identify TWFE bias sources | Goodman-Bacon (2021) |
| Placebo reforms | False timing test | Standard |
| Heterogeneity analysis | Effect variation | Callaway-Sant'Anna (2021) |
| Composition sensitivity | Never vs not-yet-treated | Sun & Abraham (2021) |

**Bacon Decomposition:**

Decomposes TWFE estimate into:
- Clean comparisons (treated vs. never-treated)
- Potentially biased comparisons (earlier vs. later treated)

Report weight on each component type.

---

## 4. Mechanism Index

### 4.1 Construction

The Mechanism Index captures five dimensions of procurement transparency:

| Dimension | Components | Weight |
|-----------|------------|--------|
| Disclosure | Tender/award notices, prices, criteria | 0.25 |
| Participation | E-submission, complaint mechanisms | 0.20 |
| Oversight | Audit trails, amendments | 0.20 |
| Integrity | Debarment checks, ownership disclosure | 0.20 |
| Efficiency | E-procurement, framework agreements | 0.15 |

### 4.2 Validation

See `src/mechanism_validation.py` for validation framework:

1. **Inter-annotator agreement**: Krippendorff's α ≥ 0.80
2. **Internal consistency**: Cronbach's α ≥ 0.70
3. **Cross-lingual equivalence**: r ≥ 0.85 across languages
4. **Construct validity**: Convergent/discriminant correlations

### 4.3 Cross-Lingual Strategy

| Language | Country | Approach |
|----------|---------|----------|
| English | UK | Native annotation |
| Ukrainian | Ukraine | Back-translation |
| Spanish | Colombia | Parallel annotation |

---

## 5. Innovation Linkage

### 5.1 Patent Data Pipeline

We link procurement winners to patent activity using:

```
Supplier Name → Fuzzy Match → Patent Applicant → Patent Family
```

**Matching Algorithm:**
1. Text normalization (case, suffixes, punctuation)
2. TF-IDF vectorization
3. Cosine similarity with threshold 0.85
4. Manual review for borderline cases

See `src/patent_linkage.py` for implementation.

### 5.2 Data Sources

| Source | Coverage | API |
|--------|----------|-----|
| Lens | Global, 130M+ patents | REST API (500 req/day) |
| EPO OPS | European patents | REST API |
| PatentsView | US patents only | REST API |

**Triangulation:**
- Primary: Lens API
- Cross-validation: EPO OPS for European
- US validation: PatentsView

---

## 6. Statistical Inference

### 6.1 Multiple Hypothesis Correction

For multiple outcomes, we apply:
- **Bonferroni**: Conservative, used for primary outcomes
- **Benjamini-Hochberg**: FDR control for secondary outcomes
- **Romano-Wolf**: Stepdown, accounts for correlation

### 6.2 Confidence Intervals

All confidence intervals are 95% unless noted:
- RDD: CCT robust confidence intervals
- DiD: Cluster-robust at agency level

### 6.3 Effect Size Interpretation

| Cohen's d | Interpretation |
|-----------|----------------|
| 0.2 | Small |
| 0.5 | Medium |
| 0.8 | Large |

---

## 7. Sustainability Metrics

**Rationale:** For Nature Sustainability, sustainability must be an outcome, not just motivation. We quantify three core sustainability endpoints.

### 7.1 Carbon Footprint Intensity

**Approach:** Combine procurement data with EXIOBASE Multi-Regional Input-Output (MRIO) tables to estimate embodied carbon.

**Concordance Pipeline:**
```
CPV Code → ISIC Rev. 4 → EXIOBASE Sector → CO₂e/$ Intensity
```

**Data Sources:**
- EXIOBASE 3.8: Environmental extensions by sector/region
- UN SIEC: ISIC ↔ CPC correspondence tables
- CPV 2008: EU procurement classification

**Formula:**
$$\text{Carbon Intensity}_{it} = \sum_{s} w_{is} \cdot \text{EXIO}_s^{\text{CO}_2\text{e}/\$}$$

where $w_{is}$ is the share of contract value in EXIOBASE sector $s$.

**Estimation:**
- Unit: kg CO₂-equivalent per $ contract value
- Regional specificity: Use EXIOBASE country-specific coefficients where available
- Uncertainty: Monte Carlo propagation of EXIOBASE uncertainty ranges

### 7.2 Supply Chain Resilience

We construct a multi-dimensional resilience index:

| Component | Measure | Formula | Weight |
|-----------|---------|---------|--------|
| Concentration | HHI | $\sum s_i^2$ | 0.30 |
| Delivery | On-time % | % within contracted deadline | 0.25 |
| Diversification | Shannon entropy | $-\sum p_i \ln(p_i)$ | 0.25 |
| Geographic | Geographic spread | Std. dev. of supplier locations | 0.20 |

**HHI Calculation:**
$$\text{HHI} = \sum_{i=1}^{N} \left(\frac{\text{value}_i}{\text{total value}}\right)^2$$

Interpretation:
- HHI < 0.01: Competitive (low concentration)
- HHI 0.01-0.18: Moderate concentration
- HHI > 0.18: High concentration

**Composite Index:**
$$\text{Resilience} = 0.30 \times (1 - \text{HHI}) + 0.25 \times \text{Delivery} + 0.25 \times \text{Diversity} + 0.20 \times \text{Geographic}$$

### 7.3 Green Procurement Share

**Definition:** Fraction of procurement value meeting environmental criteria.

**Three Identification Methods:**

1. **Keyword Matching:**
   - Environmental: sustainable, green, eco-friendly, renewable, recyclable
   - Exclusions: "greenwashing" indicators
   
2. **Certified Suppliers:**
   - ISO 14001 (Environmental Management)
   - EMAS (EU Eco-Management)
   - EU Ecolabel suppliers

3. **Green CPV Categories:**
   - 09310000: Solar energy
   - 09320000: Steam/hot water
   - 34100000: Electric vehicles
   - 90500000: Refuse/waste services

**Formula:**
$$\text{Green Share} = \frac{\sum \text{Green Contract Values}}{\sum \text{All Contract Values}}$$

### 7.4 Cross-Country Harmonization

| Challenge | Solution |
|-----------|----------|
| Currency | Convert to USD PPP (World Bank ICP) |
| CPV versions | Map all to CPV 2008 |
| Sector definitions | Use EXIOBASE cross-walks |
| Temporal coverage | Align to common observation window |

See `src/sustainability_metrics.py` for implementation.

---

## 8. Negative Controls (Falsification Tests)

**Purpose:** Demonstrate causal identification by showing treatment does NOT affect outcomes it theoretically should NOT affect.

### 8.1 Outcome Negative Controls

**Logic:** Procurement reforms should NOT affect:
- Pre-determined outcomes (fixed before reform)
- Theoretically unrelated outcomes

**Implementation:**
| Control | Rationale |
|---------|-----------|
| Pre-period procurement | Fixed before treatment |
| Established firm characteristics | Not plausibly affected |
| Lagged outcomes (t-2, t-3) | Pre-determined |

**Test:** RDD/DiD on negative control outcomes should yield:
- Point estimate ≈ 0
- Confidence interval containing 0
- |t-stat| < 1.96

### 8.2 Exposure Negative Controls

**Logic:** Groups that should NOT be affected by treatment should show no effect.

**Implementation:**
| Control Group | Rationale |
|---------------|-----------|
| Defense procurement | Often exempt from rules |
| Non-CPV sectors | Below radar |
| Pre-existing contracts | Grandfathered |

**Estimation:**
$$\hat{\tau}_{\text{placebo}} = \hat{\tau}_{\text{control}} - \hat{\tau}_{\text{treated}}$$

Expected: $\hat{\tau}_{\text{placebo}} \approx 0$

### 8.3 Temporal Placebos

**Logic:** Treatment effects should NOT appear at fake reform dates.

**Implementation:**
```python
for fake_date in [true_date - 12m, true_date - 24m, true_date - 36m]:
    estimate DiD with fake_date as treatment
    test if effect ≈ 0
```

**Joint Test:**
$$H_0: \tau_{-12} = \tau_{-24} = \tau_{-36} = 0$$

Bonferroni-corrected p-values reported.

### 8.4 Geographic Placebos

**Logic:** Effects should NOT appear at fake geographic boundaries.

**Implementation:**
- Use adjacent untreated regions as placebo
- Estimate "treatment" effects at non-threshold boundaries
- Compare to true threshold effects

**Expected:** Placebo effects statistically indistinguishable from zero.

### 8.5 Interpretation Guidelines

| Result | Interpretation |
|--------|----------------|
| All controls ≈ 0 | Strong causal credibility |
| Some controls ≠ 0 | Investigate confounding |
| Most controls ≠ 0 | Identification strategy compromised |

**Reporting Standard:** Report ALL negative control results, even if unfavorable.

See `src/negative_controls.py` for implementation.

---

## 9. Software and Reproducibility

### 9.1 Software Stack

| Task | Package | Version |
|------|---------|---------|
| RDD estimation | rddensity (R) | 2.4 |
| DiD estimation | did (R) | 2.1 |
| Data processing | pandas | 2.0+ |
| Visualization | matplotlib | 3.7+ |
| Sustainability metrics | src/sustainability_metrics.py | Custom |
| Negative controls | src/negative_controls.py | Custom |

### 9.2 Reproducibility

All analyses are reproducible via:
```bash
snakemake --cores all
```

See `workflow/Snakefile` for complete pipeline.

---

## 10. References

### RDD Methods

- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica*, 82(6), 2295-2326.

- Cattaneo, M. D., Jansson, M., & Ma, X. (2020). Simple local polynomial density estimators. *Journal of the American Statistical Association*, 115(531), 1449-1455.

- Imbens, G., & Kalyanaraman, K. (2012). Optimal bandwidth choice for the regression discontinuity estimator. *Review of Economic Studies*, 79(3), 933-959.

- Lee, D. S., & Lemieux, T. (2010). Regression discontinuity designs in economics. *Journal of Economic Literature*, 48(2), 281-355.

- McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. *Journal of Econometrics*, 142(2), 698-714.

### DiD Methods

- Callaway, B., & Sant'Anna, P. H. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.

- Goodman-Bacon, A. (2021). Difference-in-differences with variation in treatment timing. *Journal of Econometrics*, 225(2), 254-277.

- Roth, J. (2022). Pretest with caution: Event-study estimates after testing for parallel trends. *AER: Insights*, 4(3), 305-322.

- Sun, L., & Abraham, S. (2021). Estimating dynamic treatment effects in event studies with heterogeneous treatment effects. *Journal of Econometrics*, 225(2), 175-199.

### Validation Methods

- Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. *Psychological Bulletin*, 76(5), 378-382.

- Krippendorff, K. (2018). *Content analysis: An introduction to its methodology* (4th ed.). Sage.

### Sustainability & Environmental Accounting

- Stadler, K., Wood, R., Buber, T., et al. (2018). EXIOBASE 3: Developing a time series of detailed environmentally extended multi-regional input-output tables. *Journal of Industrial Ecology*, 22(3), 502-515.

- Wiedmann, T. O., & Lenzen, M. (2018). Environmental and social footprints of international trade. *Nature Geoscience*, 11(5), 314-321.

- Södersten, C. J., Wood, R., & Hertwich, E. G. (2018). Environmental impacts of capital formation. *Journal of Industrial Ecology*, 22(1), 55-67.

### Negative Controls

- Lipsitch, M., Tchetgen Tchetgen, E., & Cohen, T. (2010). Negative controls: A tool for detecting confounding and bias in observational studies. *Epidemiology*, 21(3), 383-388.

- Shi, X., Miao, W., & Tchetgen Tchetgen, E. (2020). A selective review of negative control methods in epidemiology. *Current Epidemiology Reports*, 7(4), 190-202.
