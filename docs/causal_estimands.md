# Causal Estimand Definitions

## Overview

This document provides formal definitions of causal estimands for the e-procurement transparency analysis. We employ two complementary identification strategies:

1. **Regression Discontinuity Design (RDD)** for threshold-based policies
2. **Difference-in-Differences (DiD)** for staggered reform adoption

All definitions follow the potential outcomes framework (Rubin, 1974; Imbens & Angrist, 1994).

---

## 1. Potential Outcomes Framework

### 1.1 Basic Setup

Let $i$ index procurement contracts and $t$ index time periods.

**Potential Outcomes:**
- $Y_i(1)$: Outcome if contract $i$ is subject to enhanced transparency (treated)
- $Y_i(0)$: Outcome if contract $i$ is subject to standard procedures (control)

**Observed Outcome:**
$$Y_i = D_i \cdot Y_i(1) + (1 - D_i) \cdot Y_i(0)$$

where $D_i \in \{0, 1\}$ is the treatment indicator.

**Fundamental Problem:** We observe only one potential outcome per unit.

### 1.2 Treatment Effects

**Individual Treatment Effect (ITE):**
$$\tau_i = Y_i(1) - Y_i(0)$$

**Average Treatment Effect (ATE):**
$$\tau_{ATE} = \mathbb{E}[Y_i(1) - Y_i(0)]$$

**Average Treatment Effect on the Treated (ATT):**
$$\tau_{ATT} = \mathbb{E}[Y_i(1) - Y_i(0) | D_i = 1]$$

**Local Average Treatment Effect (LATE):**
$$\tau_{LATE} = \mathbb{E}[Y_i(1) - Y_i(0) | \text{Compliers}]$$

---

## 2. Regression Discontinuity Design (RDD)

### 2.1 Setting

We exploit procurement value thresholds that trigger enhanced transparency requirements.

**Running Variable:** $X_i$ = Contract value
**Cutoff:** $c$ = Threshold value (e.g., €1M for EU directives)
**Treatment:** $D_i = \mathbf{1}(X_i \geq c)$

### 2.2 Sharp RDD Estimand

**Definition:** The Sharp RDD estimand is the local average treatment effect at the threshold:

$$\tau_{RDD} = \lim_{x \downarrow c} \mathbb{E}[Y_i | X_i = x] - \lim_{x \uparrow c} \mathbb{E}[Y_i | X_i = x]$$

**In potential outcomes notation:**
$$\tau_{RDD} = \mathbb{E}[Y_i(1) - Y_i(0) | X_i = c]$$

### 2.3 Identification Assumptions

**Assumption 1 (Continuity):**
$$\mathbb{E}[Y_i(0) | X_i = x] \text{ and } \mathbb{E}[Y_i(1) | X_i = x] \text{ are continuous in } x \text{ at } c$$

This requires:
- No manipulation of the running variable at the cutoff
- No other policies change discontinuously at $c$

**Assumption 2 (No Precise Manipulation):**
The density $f(x)$ of the running variable is continuous at $c$:
$$\lim_{x \downarrow c} f(x) = \lim_{x \uparrow c} f(x)$$

Testable via McCrary (2008) density test.

**Assumption 3 (SUTVA):**
- No interference: $Y_i$ depends only on $D_i$, not on $D_j$ for $j \neq i$
- No hidden treatment variations

### 2.4 Estimator Specification

**Local Linear Regression (Imbens & Kalyanaraman, 2012):**

$$Y_i = \alpha + \tau \cdot D_i + \beta_1 (X_i - c) + \beta_2 D_i (X_i - c) + \epsilon_i$$

for $X_i \in [c - h, c + h]$, where $h$ is the optimal bandwidth.

**Bandwidth Selection:**
- IK optimal bandwidth (Imbens & Kalyanaraman, 2012)
- CCT robust bandwidth (Calonico, Cattaneo & Titiunik, 2014)

### 2.5 Primary Outcomes

| Outcome | Definition | Hypothesis |
|---------|------------|------------|
| Competition | Number of qualifying bidders | $\tau > 0$: Transparency increases competition |
| Price | Log contract price / estimated value | $\tau < 0$: Transparency reduces prices |
| Time | Days from tender to contract | $\tau < 0$: Transparency reduces delays |
| Quality | Post-contract amendment rate | $\tau < 0$: Transparency reduces renegotiation |

### 2.6 Heterogeneous Effects

We estimate conditional average treatment effects:
$$\tau(z) = \mathbb{E}[Y_i(1) - Y_i(0) | X_i = c, Z_i = z]$$

where $Z_i$ includes:
- Sector (construction, IT, services)
- Contracting authority type (central, local)
- Prior corruption risk

---

## 3. Difference-in-Differences (DiD)

### 3.1 Setting

We exploit staggered adoption of e-procurement reforms across countries/agencies.

**Panel Structure:** $(i, t)$ where $i$ = procurement agency, $t$ = time period
**Treatment Timing:** $G_i$ = First period of reform adoption (∞ if never-treated)
**Treatment:** $D_{it} = \mathbf{1}(t \geq G_i)$

### 3.2 Standard DiD Estimand (2x2)

For a single treatment timing $G$:

$$\tau_{DiD} = \mathbb{E}[Y_{it}(1) - Y_{it}(0) | G_i = G, t \geq G] - \mathbb{E}[Y_{it}(1) - Y_{it}(0) | G_i = G, t < G]$$

Under parallel trends, this simplifies to:

$$\tau_{DiD} = (\bar{Y}_{treated, post} - \bar{Y}_{treated, pre}) - (\bar{Y}_{control, post} - \bar{Y}_{control, pre})$$

### 3.3 Staggered DiD Estimand (Callaway-Sant'Anna)

With staggered adoption, we define group-time average treatment effects:

**Definition (ATT(g,t)):**
$$ATT(g, t) = \mathbb{E}[Y_{it}(g) - Y_{it}(\infty) | G_i = g]$$

for groups $g$ (cohorts by first treatment period) and times $t \geq g$.

**Aggregated ATT:**
$$ATT = \sum_{g} \sum_{t \geq g} w_{g,t} \cdot ATT(g, t)$$

where weights $w_{g,t}$ depend on sample sizes.

### 3.4 Identification Assumptions

**Assumption 1 (Parallel Trends):**
$$\mathbb{E}[Y_{it}(0) - Y_{i,t-1}(0) | G_i = g] = \mathbb{E}[Y_{it}(0) - Y_{i,t-1}(0) | G_i = g']$$

for all $g, g'$ and $t < \min(g, g')$.

*Interpretation:* In the absence of treatment, treated and control groups would have followed the same trajectory.

**Assumption 2 (No Anticipation):**
$$Y_{it}(g) = Y_{it}(\infty) \text{ for all } t < g$$

*Interpretation:* No effect of treatment before it is implemented.

**Assumption 3 (SUTVA):**
- No interference across units
- Well-defined treatment

### 3.5 Event Study Specification

**Model:**
$$Y_{it} = \alpha_i + \gamma_t + \sum_{k \neq -1} \tau_k \cdot \mathbf{1}(t - G_i = k) + X_{it}'\beta + \epsilon_{it}$$

where:
- $\alpha_i$: Unit fixed effects
- $\gamma_t$: Time fixed effects
- $\tau_k$: Event-time coefficient (effect $k$ periods from treatment)
- $k = -1$ normalized to zero (base period)

**Pre-trend Test:**
$H_0: \tau_{-2} = \tau_{-3} = ... = \tau_{-K} = 0$

### 3.6 Primary Outcomes

| Outcome | Definition | Hypothesis |
|---------|------------|------------|
| Savings | % difference from estimated cost | $ATT > 0$: Reforms increase savings |
| SME Share | % contracts to small/medium enterprises | $ATT > 0$: Reforms improve SME access |
| Processing Time | Days from tender to contract | $ATT < 0$: Reforms reduce delays |
| Bid Count | Average bidders per tender | $ATT > 0$: Reforms increase competition |

---

## 4. Causal Mechanisms

### 4.1 Mechanism Index

We construct a Mechanism Index $M_i$ capturing transparency dimensions:

$$M_i = \sum_{d=1}^{D} w_d \cdot M_{id}$$

where:
- $M_{id}$: Score for dimension $d$ (disclosure, participation, oversight, integrity, efficiency)
- $w_d$: Dimension weight (equal or empirically derived)

### 4.2 Mediation Analysis

**Direct Effect:**
$$\tau_{direct} = \mathbb{E}[Y_i(1, M_i(0)) - Y_i(0, M_i(0))]$$

**Indirect Effect (via mechanism):**
$$\tau_{indirect} = \mathbb{E}[Y_i(1, M_i(1)) - Y_i(1, M_i(0))]$$

**Total Effect:**
$$\tau_{total} = \tau_{direct} + \tau_{indirect}$$

---

## 5. Threats to Identification

### 5.1 RDD Threats

| Threat | Diagnostic | Mitigation |
|--------|-----------|------------|
| Manipulation | McCrary density test | Report and bound estimates |
| Covariate discontinuity | Balance tests at cutoff | Include covariates |
| Bandwidth sensitivity | Sweep $h \in [0.5h^*, 2h^*]$ | Report range of estimates |
| Functional form | Vary polynomial order | Use local linear (preferred) |

### 5.2 DiD Threats

| Threat | Diagnostic | Mitigation |
|--------|-----------|------------|
| Pre-trend violation | Event study visualization | Bound using Roth (2022) |
| Negative weighting | Bacon decomposition | Use Callaway-Sant'Anna |
| Staggered timing | Check heterogeneity by cohort | Report by treatment group |
| Composition effects | Test never vs not-yet-treated | Report both control groups |

---

## 6. Power and Sample Size

### 6.1 Minimum Detectable Effect (MDE)

For RDD with effective sample $n_h$ in bandwidth:

$$MDE = 2.8 \sqrt{\frac{\sigma^2}{n_h}} \cdot \frac{1}{\sqrt{1 - R^2}}$$

where $\sigma^2$ is outcome variance and $R^2$ is from covariates.

### 6.2 Target Sample Sizes

| Analysis | Required N | Available N | Power |
|----------|-----------|-------------|-------|
| RDD (Ukraine) | 2,000 | ~50,000 | >0.99 |
| RDD (Colombia) | 2,000 | ~30,000 | >0.99 |
| DiD (cross-country) | 500 agencies | 1,000+ | >0.95 |

---

## 7. Pre-Registration Elements

### 7.1 Primary Hypotheses

1. **H1 (Competition):** Enhanced transparency requirements increase the number of bidders
   - Estimand: $\tau_{RDD}$ for bid count at threshold
   - Direction: Positive
   - Significance: $\alpha = 0.05$, two-tailed

2. **H2 (Prices):** Enhanced transparency reduces contract prices
   - Estimand: $\tau_{RDD}$ for log(price/estimate) at threshold
   - Direction: Negative
   - Significance: $\alpha = 0.05$, two-tailed

3. **H3 (Mechanism):** Effects operate through increased information disclosure
   - Estimand: $\tau_{indirect}$ via Mechanism Index
   - Direction: Positive
   - Significance: $\alpha = 0.05$, one-tailed

### 7.2 Secondary/Exploratory Hypotheses

- Heterogeneity by sector
- Heterogeneity by contracting authority
- Long-term dynamic effects (DiD event study)
- Cross-country variation in effect sizes

---

## References

- Callaway, B., & Sant'Anna, P. H. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica*, 82(6), 2295-2326.
- Goodman-Bacon, A. (2021). Difference-in-differences with variation in treatment timing. *Journal of Econometrics*, 225(2), 254-277.
- Imbens, G., & Angrist, J. (1994). Identification and estimation of local average treatment effects. *Econometrica*, 62(2), 467-475.
- Imbens, G., & Kalyanaraman, K. (2012). Optimal bandwidth choice for the regression discontinuity estimator. *Review of Economic Studies*, 79(3), 933-959.
- McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. *Journal of Econometrics*, 142(2), 698-714.
- Roth, J. (2022). Pretest with caution: Event-study estimates after testing for parallel trends. *AER: Insights*, 4(3), 305-322.
- Rubin, D. B. (1974). Estimating causal effects of treatments in randomized and nonrandomized studies. *Journal of Educational Psychology*, 66(5), 688-701.
