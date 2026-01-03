# Baseline Expansion Comparisons

**Purpose:** Document the progression from simple to sophisticated models, demonstrating added value of each component for Nature Sustainability reviewers.

---

## Overview

Reviewers will ask: "Why do you need all this complexity?" This document provides systematic evidence that each methodological component adds explanatory/predictive power.

**Baseline Ladder:**
1. Naive cross-sectional comparison
2. Simple DiD (no heterogeneity)
3. Full causal inference (RDD + Callaway-Sant'Anna)
4. + Mechanism Index
5. + Sustainability metrics
6. + Cross-country pooling

---

## 1. Baseline Models

### 1.1 Model 0: Naive Cross-Sectional

**Specification:**
$$Y_i = \alpha + \beta \cdot \text{Transparent}_i + \epsilon_i$$

**What it captures:** Simple treated vs. untreated comparison

**What it misses:**
- Selection bias (better agencies adopt transparency)
- Time trends (outcomes improving anyway)
- Confounders (economic conditions, political will)

**Expected bias:** Likely **overestimates** treatment effect due to positive selection

### 1.2 Model 1: Simple DiD (TWFE)

**Specification:**
$$Y_{it} = \alpha_i + \gamma_t + \beta \cdot D_{it} + \epsilon_{it}$$

**What it captures:** 
- Unit fixed effects (time-invariant confounders)
- Time fixed effects (common shocks)

**What it misses:**
- Heterogeneous treatment effects
- Staggered adoption bias (Goodman-Bacon 2021)
- Treatment effect dynamics

**Known problem:** With staggered timing, TWFE weights can be negative → biased estimates

### 1.3 Model 2: Callaway-Sant'Anna DiD

**Specification:**
$$ATT(g,t) = \mathbb{E}[Y_t(g) - Y_t(\infty) | G = g]$$

**What it adds:**
- Group-time specific treatment effects
- Proper aggregation without negative weights
- Dynamic treatment effect estimation

**Comparison to TWFE:**

| Metric | TWFE | C-S DiD |
|--------|------|---------|
| Handles heterogeneity | ❌ | ✅ |
| Event study valid | ⚠️ | ✅ |
| Aggregation unbiased | ❌ | ✅ |

### 1.4 Model 3: RDD at Thresholds

**Specification:**
$$Y_i = \alpha + \tau \cdot D_i + \beta_1(X_i - c) + \beta_2 D_i(X_i - c) + \epsilon_i$$

**What it adds:**
- Local randomization at threshold
- Cleaner identification than DiD
- Free of parallel trends assumption

**Limitation:** 
- Local Average Treatment Effect (LATE) only
- External validity to infra-marginal contracts unclear

---

## 2. Component Value-Add Analysis

### 2.1 Mechanism Index Addition

**Model without Mechanism Index:**
$$\text{Sustainability}_i = f(\text{Transparency}_{binary})$$

**Model with Mechanism Index:**
$$\text{Sustainability}_i = f(\text{Mechanism Index}_{continuous})$$

**Value-Added:**

| Metric | Binary Treatment | Mechanism Index |
|--------|-----------------|-----------------|
| R² | [baseline] | [+X%] |
| RMSE | [baseline] | [-Y%] |
| AIC | [baseline] | [improvement] |
| Causal precision | Low | High |

**Interpretation:** The Mechanism Index captures *intensity* of transparency, not just presence/absence. This allows:
- Dose-response relationships
- Cross-country comparability
- Component-specific effects (disclosure vs. oversight)

### 2.2 Sustainability Metrics Addition

**Model with competition outcomes only:**
$$\text{Competition}_i = f(\text{Treatment}_i)$$

**Model with sustainability outcomes:**
$$\text{Sustainability}_i = f(\text{Treatment}_i)$$

**Nature Sustainability Value:**

| Outcome | Policy Relevance | Journal Fit |
|---------|-----------------|-------------|
| Competition (bids) | Medium | Economics journals |
| Price savings | Medium | Economics journals |
| **Carbon intensity** | **High** | **Nature Sustainability** |
| **Resilience** | **High** | **Nature Sustainability** |
| **Green share** | **High** | **Nature Sustainability** |

**Key insight:** Without sustainability metrics, this is an economics paper. With them, it's a Nature Sustainability paper.

### 2.3 Cross-Country Pooling

**Single-country model:**
$$Y_{it}^{Ukraine} = f(\text{ProZorro}_t)$$

**Pooled model:**
$$Y_{ijt} = f(\text{Reform}_{jt}) + \theta_j + \text{Country} \times \text{Time FE}$$

**Value-Added:**

| Metric | Single Country | Pooled 3 Countries |
|--------|---------------|-------------------|
| Statistical power | Low | High |
| External validity | Limited | Global |
| Mechanism generalizability | Uncertain | Demonstrated |

---

## 3. Model Comparison Table

### 3.1 Competition Outcome (Bids per Tender)

| Model | Estimate | SE | R² | AIC |
|-------|----------|----|----|-----|
| 0. Cross-sectional | [β₀] | [se₀] | [r₀] | [aic₀] |
| 1. Simple TWFE | [β₁] | [se₁] | [r₁] | [aic₁] |
| 2. Callaway-Sant'Anna | [β₂] | [se₂] | [r₂] | [aic₂] |
| 3. RDD | [β₃] | [se₃] | [r₃] | [aic₃] |
| 4. + Mechanism Index | [β₄] | [se₄] | [r₄] | [aic₄] |
| 5. + Cross-country | [β₅] | [se₅] | [r₅] | [aic₅] |

### 3.2 Carbon Footprint Intensity

| Model | Estimate | SE | R² | AIC |
|-------|----------|----|----|-----|
| 0. Cross-sectional | [β₀] | [se₀] | [r₀] | [aic₀] |
| 1. Simple TWFE | [β₁] | [se₁] | [r₁] | [aic₁] |
| 2. Callaway-Sant'Anna | [β₂] | [se₂] | [r₂] | [aic₂] |
| 3. + Sustainability | [β₃] | [se₃] | [r₃] | [aic₃] |
| 4. + Mechanism Index | [β₄] | [se₄] | [r₄] | [aic₄] |
| 5. + Cross-country | [β₅] | [se₅] | [r₅] | [aic₅] |

---

## 4. Ablation Studies

### 4.1 Component Removal Impact

| Component Removed | Δ Estimate | Δ SE | Interpretation |
|-------------------|------------|------|----------------|
| Remove RDD | [Δ] | [Δ] | Threshold design essential |
| Remove C-S DiD | [Δ] | [Δ] | Heterogeneity matters |
| Remove Mechanism Index | [Δ] | [Δ] | Intensity vs. binary |
| Remove cross-country | [Δ] | [Δ] | Generalizability |

### 4.2 Feature Importance (Mechanism Index Components)

| Component | Marginal R² | Importance Rank |
|-----------|-------------|-----------------|
| Disclosure | [r²] | [1-5] |
| Participation | [r²] | [1-5] |
| Oversight | [r²] | [1-5] |
| Integrity | [r²] | [1-5] |
| Efficiency | [r²] | [1-5] |

---

## 5. Alternative Specifications

### 5.1 Procurement-Only Models

**Specification:** Use only procurement data (no NLP, no mechanism index)

$$Y_i = f(\text{Contract Features Only})$$

**Variables:**
- Contract value
- Procurement method
- Number of bidders
- Award date
- Agency fixed effects

**Limitation:** Cannot capture *quality* of transparency, only presence

### 5.2 NLP-Only Models

**Specification:** Use text-derived features without causal design

$$Y_i = f(\text{NLP Features})$$

**Variables:**
- Tender description embeddings
- Keyword frequencies
- Document length
- Language complexity

**Limitation:** Purely correlational, no causal identification

### 5.3 Combined vs. Separate

| Approach | Advantages | Disadvantages |
|----------|------------|---------------|
| Procurement-only | Transparent, reproducible | Misses mechanism quality |
| NLP-only | Rich features | No causal identification |
| **Combined** | **Best of both** | **Complexity** |

---

## 6. Prediction Performance

### 6.1 Out-of-Sample R²

Train on 80% of data, test on 20%:

| Model | Train R² | Test R² | Overfitting? |
|-------|----------|---------|--------------|
| 0. Cross-sectional | [r] | [r] | [Yes/No] |
| 1. TWFE | [r] | [r] | [Yes/No] |
| 2. C-S DiD | [r] | [r] | [Yes/No] |
| Full model | [r] | [r] | [Yes/No] |

### 6.2 Cross-Country Prediction

Train on 2 countries, predict 3rd:

| Train | Test | R² |
|-------|------|-----|
| Ukraine + Colombia | UK | [r] |
| Ukraine + UK | Colombia | [r] |
| Colombia + UK | Ukraine | [r] |

---

## 7. Computational Considerations

| Model | Estimation Time | Memory | Scalability |
|-------|-----------------|--------|-------------|
| Cross-sectional | <1 min | Low | Excellent |
| TWFE | <5 min | Medium | Good |
| Callaway-Sant'Anna | 10-30 min | High | Moderate |
| Full pipeline | 1-2 hours | High | Requires HPC |

---

## 8. Reviewer Response Template

For "Why so complex?" questions:

> "Table S[X] shows systematic model comparisons. Simple cross-sectional models overestimate treatment effects by [X]% due to selection bias. Standard TWFE is biased by [X]% with staggered adoption. Our full specification—combining RDD at thresholds, Callaway-Sant'Anna for dynamics, Mechanism Index for intensity, and cross-country pooling for generalizability—provides the most credible estimates. Ablation studies confirm each component contributes meaningfully: removing [component] changes estimates by [Y]%. The added complexity is necessary for credible causal claims about sustainability outcomes."

---

## 9. Code References

| Analysis | File | Function |
|----------|------|----------|
| Baseline comparisons | `scripts/run_baselines.py` | `compare_models()` |
| Ablation | `scripts/run_ablation.py` | `ablation_study()` |
| Cross-validation | `src/validation.py` | `cross_validate()` |
| Model diagnostics | `src/diagnostics.py` | `model_comparison()` |

---

## 10. Summary Statistics for Baselines

### 10.1 Data Coverage

| Country | N (contracts) | T (periods) | Unique Agencies |
|---------|---------------|-------------|-----------------|
| Ukraine | [n] | [t] | [k] |
| Colombia | [n] | [t] | [k] |
| UK | [n] | [t] | [k] |
| **Pooled** | **[N]** | **[T]** | **[K]** |

### 10.2 Outcome Distributions

| Outcome | Mean | SD | Min | Max |
|---------|------|-----|-----|-----|
| Bids per tender | [μ] | [σ] | [min] | [max] |
| Carbon intensity | [μ] | [σ] | [min] | [max] |
| Resilience index | [μ] | [σ] | [min] | [max] |
| Green share | [μ] | [σ] | [min] | [max] |
