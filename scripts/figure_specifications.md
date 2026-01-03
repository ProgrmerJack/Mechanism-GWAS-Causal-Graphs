# Figure Specifications for Nature Sustainability

## Overview

This document specifies the "killer figure" set for submission to **Nature Sustainability** (Research Article format). Each figure is designed to:

1. **Headline sustainability outcomes** - carbon footprint, resilience, green procurement
2. Meet the ≤8 display items constraint (5 main + 3 extended data)
3. Tell the story: Rules → Market Structure → **Sustainability Outcomes**

**Core Narrative**: Procurement transparency rules causally reshape market structure and innovation diffusion in ways that change the carbon/resilience profile of public spending—at scale—across countries.

---

## Figure Format Requirements

### Journal Standards

| Specification | Nature | Science |
|--------------|--------|---------|
| Width (single column) | 88 mm | 85 mm |
| Width (double column) | 180 mm | 178 mm |
| Resolution | 300 dpi minimum | 300 dpi minimum |
| Format | PDF, EPS, TIFF | PDF, EPS, TIFF |
| Font | Arial/Helvetica | Arial |
| Font size (labels) | 7-9 pt | 7-9 pt |
| Line width | 0.5-1.5 pt | 0.5-1.5 pt |

### Color Palette

Use colorblind-safe palette (Okabe-Ito):

| Color | Hex | Use |
|-------|-----|-----|
| Blue | #0072B2 | Treatment/Positive |
| Orange | #E69F00 | Control/Negative |
| Green | #009E73 | Secondary positive |
| Red | #D55E00 | Warnings/Fails |
| Purple | #CC79A7 | Neutral |
| Gray | #999999 | Background |

---

## Figure 1: Control Surface Schematic

### Purpose
Illustrate the multi-layered causal identification strategy combining RDD and DiD.

### Layout
Three-panel layout showing:
1. **Panel A**: RDD threshold visualization
2. **Panel B**: DiD event study
3. **Panel C**: Mechanism pathway diagram

### Panel A Specifications

**Title:** "Transparency Threshold Effects"

**Elements:**
- X-axis: Contract Value (normalized to threshold)
- Y-axis: Outcome (e.g., Number of Bidders)
- Scatter plot with local polynomial fit
- Vertical dashed line at cutoff (c=0)
- Shaded confidence bands

**Annotations:**
- τ_RDD estimate with SE in box
- McCrary p-value in corner
- Optimal bandwidth range shaded

**Code Template:**
```python
fig, ax = plt.subplots(figsize=(88/25.4, 70/25.4))  # Single column

# Scatter
ax.scatter(X[below], Y[below], alpha=0.3, c='#0072B2', s=10)
ax.scatter(X[above], Y[above], alpha=0.3, c='#E69F00', s=10)

# Local polynomial fit
ax.plot(x_fit_below, y_fit_below, c='#0072B2', lw=1.5)
ax.plot(x_fit_above, y_fit_above, c='#E69F00', lw=1.5)

# Confidence bands
ax.fill_between(x_fit_below, ci_lower, ci_upper, alpha=0.2)

# Cutoff line
ax.axvline(0, color='black', linestyle='--', lw=0.8)

# Annotation box
ax.text(0.05, 0.95, f'τ = {tau:.3f}***\n(SE = {se:.3f})',
        transform=ax.transAxes, fontsize=8,
        bbox=dict(boxstyle='round', facecolor='white'))
```

### Panel B Specifications

**Title:** "Dynamic Treatment Effects"

**Elements:**
- X-axis: Event Time (periods relative to reform)
- Y-axis: Treatment Effect Coefficient
- Point estimates with 95% CI bars
- Horizontal line at y=0
- Vertical line at event time = -0.5 (treatment onset)
- Pre-treatment period shaded gray

**Annotations:**
- Pre-trend test p-value
- ATT estimate with SE

**Code Template:**
```python
ax.errorbar(event_times, coefficients,
            yerr=[coef - ci_lower, ci_upper - coef],
            fmt='o', capsize=3, color='#0072B2', markersize=5)

ax.axhline(0, color='black', linestyle='--', lw=0.8)
ax.axvline(-0.5, color='#D55E00', linestyle='--', lw=0.8)
ax.axvspan(min(event_times)-0.5, -0.5, alpha=0.1, color='gray')
```

### Panel C Specifications

**Title:** "Causal Pathway"

**Elements:**
- DAG (Directed Acyclic Graph) showing:
  - Treatment → Mechanism Index → Outcome
  - Direct effect arrow
  - Confounders blocked

**Style:**
- Nodes: Rounded rectangles
- Arrows: Curved with heads
- Labels: Clear variable names

---

## Figure 2: Cross-Country Comparison

### Purpose
Show consistency of effects across Ukraine, Colombia, and UK.

### Layout
3x3 grid: Countries (rows) × Outcomes (columns)

### Specifications

**Grid Layout:**
```
         Competition    Price       Quality
Ukraine  [RDD plot]    [RDD plot]  [RDD plot]
Colombia [RDD plot]    [RDD plot]  [RDD plot]
UK       [RDD plot]    [RDD plot]  [RDD plot]
```

**Each subplot:**
- Minimal axes (spine off)
- Shared y-axis scale within outcome
- Effect size (Cohen's d) displayed
- Significance stars: * p<0.05, ** p<0.01, *** p<0.001

**Summary panel:**
- Forest plot of meta-analytic summary
- Heterogeneity I² statistic

---

## Figure 3: Mechanism Index Decomposition

### Purpose
Show which transparency dimensions drive effects.

### Layout
**Panel A:** Stacked bar chart of Mechanism Index components by country
**Panel B:** Mediation analysis results

### Panel A Specifications

**Elements:**
- X-axis: Countries
- Y-axis: Mechanism Index score (0-100)
- Stacked bars by dimension (5 colors)
- Error bars on total score

**Dimensions:**
1. Disclosure (blue)
2. Participation (green)
3. Oversight (orange)
4. Integrity (purple)
5. Efficiency (gray)

### Panel B Specifications

**Elements:**
- Path diagram showing:
  - Treatment → Mechanism Index (a path)
  - Mechanism Index → Outcome (b path)
  - Treatment → Outcome (c' path)
- Path coefficients with SEs
- Indirect effect: a×b with bootstrap CI

---

## Figure 4: Sustainability Outcomes (HEADLINE FIGURE)

### Purpose
**Core Nature Sustainability contribution.** Show procurement rule changes causally affect sustainability—not just competition/prices.

### Layout
Four panels (2×2):
- **Panel A**: Carbon Footprint Intensity (kg CO₂e per $)
- **Panel B**: Supply Chain Resilience Index  
- **Panel C**: Green Procurement Share
- **Panel D**: Composite Sustainability Score

### Panel A: Carbon Footprint RDD

**Title:** "Procurement Carbon Intensity at Transparency Threshold"

**Method:**
1. Map CPV codes → EXIOBASE sectors
2. Apply sector-level CO₂e intensity (kg CO₂e/$)
3. RDD on weighted-average intensity

**Elements:**
- X-axis: Contract Value (normalized to threshold)
- Y-axis: Carbon Intensity (kg CO₂e/$ spent)
- Local polynomial fit both sides
- Effect annotation: "Δ = -X.XX kg CO₂e/$ (95% CI)"

### Panel B: Resilience Metrics

**Title:** "Supply Chain Resilience Effects"

**Metrics (coefficient plot):**
- Supplier HHI change (lower = more diversified)
- Single-bid rate change (lower = more competition)
- Completion rate change (higher = more reliable)
- Standardized effect sizes with 95% CIs

### Panel C: Green Procurement Share

**Title:** "Green Procurement Adoption"

**Method:** Classify via CPV green codes + NLP keywords

**Elements:**
- Bar chart: Green share below vs. above threshold
- RDD estimate with CI

### Panel D: Composite Score

**Title:** "Sustainability Index Change"

**Components:**
- Carbon Score (40%): intensity → 0-100
- Resilience Score (35%): HHI, competition
- Green Score (25%): green share

**Elements:**
- Stacked bar comparison
- Net change annotation

---

## Figure 5: Mechanisms + Cross-Country Replication

### Purpose
Show (a) why sustainability improves and (b) effects replicate.

### Layout
- **Row 1**: Mechanism pathway (Panels A-B)
- **Row 2**: Country replication (Panels C-E)

### Panel A: Mediation Analysis

**Elements:**
- Path: Rules → Mechanism Index → Carbon Intensity
- Indirect effect with bootstrap CI

### Panel B: Supplier Composition Channel

**Elements:**
- Winner composition at threshold
- New entrants vs. incumbents
- Innovators (patent-holders) share

### Panels C-E: Country Replication

**One mini-RDD per country:**
- Ukraine, Colombia, UK
- Carbon intensity effect ± CI
- Meta-analytic heterogeneity I²

---

## Extended Data Figures

### ED Figure 1: Data Coverage Map

**Elements:**
- World map with countries colored by data availability
- Inset: Timeline of data collection periods
- Table: N observations by country

### ED Figure 2: OCDS Field Completeness

**Elements:**
- Heatmap: Fields (rows) × Countries (columns)
- Color scale: % complete (0-100%)
- Hierarchical clustering by similarity

### ED Figure 3: McCrary Density Tests

**Elements:**
- One subplot per threshold
- Density histogram with fitted local polynomial
- Discontinuity estimate and p-value

### ED Figure 4: Placebo Cutoffs

**Elements:**
- Multiple RDD estimates at false cutoffs
- True cutoff highlighted
- Distribution of placebo estimates

### ED Figure 5: Heterogeneity Analysis

**Elements:**
- Forest plot of subgroup effects
- Categories: Sector, Authority type, Time period
- Overall pooled estimate

### ED Figure 6: Entity Resolution Diagnostics

**Elements:**
- Precision-recall curve
- Match score distribution
- Examples of matches by confidence level

---

## Production Notes

### File Naming Convention

```
Figure_1_Control_Surface_v{version}.{format}
ED_Figure_1_Data_Coverage_v{version}.{format}
```

### Version Control

- All figures tracked in `outputs/figures/`
- Version notes in `outputs/figures/CHANGELOG.md`
- Final versions tagged with submission ID

### Review Checklist

- [ ] All text legible at print size
- [ ] Color scheme colorblind-safe
- [ ] Axes labels complete with units
- [ ] Statistical annotations include sample size
- [ ] Legend placement doesn't obscure data
- [ ] Resolution ≥ 300 dpi
- [ ] Vector format for line art

---

## Implementation

See `scripts/generate_figures.py` for automated figure generation:

```bash
python scripts/generate_figures.py --output outputs/figures/ --format pdf
```

Figure generation integrated into Snakemake workflow:

```bash
snakemake outputs/figures/all_figures.done
```
