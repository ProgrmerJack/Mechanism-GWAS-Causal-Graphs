# Nature Portfolio Reporting Summary

## For: Nature Sustainability Research Article

**Title:** Procurement Transparency Rules Causally Reshape Market Structure and Sustainability Outcomes: Evidence from Three Countries

**Author:** Abduxoliq Ashuraliyev

---

## Life Sciences Reporting Summary

> This checklist is used when submissions contain life sciences data. Leave blank or mark N/A for social sciences submissions.

☐ N/A - This is a social sciences/economics study

---

## Social Sciences Reporting Summary

### Policy Information

**1. For studies involving human research participants:**

☑ Our study uses **administrative data only** (publicly released government procurement records)
☐ No human subjects research approval required
☑ Data is from official government transparency portals (ProZorro, SECOP II, Contracts Finder)
☑ All data is aggregated/anonymized at contracting authority level

**2. Describe any IRB/ethics approval:**

☐ Not required - study uses only publicly available administrative records
☐ No individual-level data; no personally identifiable information
☐ All procurement data is released under open data licenses by respective governments

---

### Data Availability

**3. Field-specific reporting:**

| Data Source | Country | Format | Access | License |
|------------|---------|--------|--------|---------|
| ProZorro | Ukraine | OCDS JSON | api.prozorro.gov.ua | CC BY 4.0 |
| SECOP II | Colombia | OCDS JSON | colombiacompra.gov.co | Open Government |
| Contracts Finder | UK | OCDS JSON | contracts-finder.service.gov.uk | OGL v3.0 |

**4. Data processing:**

- Raw data downloaded: [dates TBD]
- Processing scripts: `workflow/rules/` (Snakemake)
- Intermediate files: `data/processed/`
- Final analysis datasets: `data/final/`

---

### Code Availability

**5. Software and algorithms:**

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| Python | 3.11+ | Primary analysis | python.org |
| R | 4.3+ | RDD/DiD robustness | r-project.org |
| rdrobust | 2.1+ | RDD estimation | Calonico et al. (2014) |
| did | 2.1+ | Callaway-Sant'Anna | Callaway & Sant'Anna (2021) |
| Snakemake | 8.0+ | Workflow management | Mölder et al. (2021) |

**6. Custom code:**

☑ All analysis code available at: [GitHub repository TBD]
☑ DOI for code: 10.5281/zenodo.17802423
☑ Code is documented with README and docstrings
☑ Environment specification: `environment.yml`

---

### Statistics

**7. Statistical analyses:**

| Analysis | Method | Software | Citation |
|----------|--------|----------|----------|
| Primary RDD | Local polynomial regression | rdrobust | Calonico et al. (2014) |
| Bandwidth selection | MSE-optimal (Imbens-Kalyanaraman) | rdrobust | Imbens & Kalyanaraman (2012) |
| DiD | Callaway-Sant'Anna | did (R) | Callaway & Sant'Anna (2021) |
| Standard errors | Robust heteroskedasticity-consistent | HC2 | MacKinnon & White (1985) |
| Clustering | At contracting authority level | cluster.vcov | Abadie et al. (2023) |

**8. Statistical tests:**

- Two-sided tests throughout unless otherwise specified
- Significance levels: * p<0.05, ** p<0.01, *** p<0.001
- Multiple comparison corrections: Benjamini-Hochberg FDR where applicable

**9. Sample sizes:**

| Country | N contracts | N authorities | Period |
|---------|-------------|---------------|--------|
| Ukraine | TBD | TBD | 2016-2023 |
| Colombia | TBD | TBD | 2016-2023 |
| UK | TBD | TBD | 2015-2023 |

**10. Effect sizes:**

- Primary effects reported with:
  - Point estimate ± SE
  - 95% confidence intervals
  - Cohen's d (standardized)
- Heterogeneity: I² statistic for cross-country meta-analysis

---

### Reproducibility

**11. Robustness checks:**

☑ Bandwidth sensitivity (0.25h* to 2h*)
☑ McCrary density test (manipulation)
☑ Donut-hole estimates (excluding near-cutoff)
☑ Placebo cutoffs
☑ Callaway-Sant'Anna (staggered DiD)
☑ Bacon decomposition (TWFE diagnostics)
☑ Pre-trend tests

**12. Pre-registration:**

☐ This study was not pre-registered
☐ Pre-registration not standard in applied economics
☑ Analysis plan documented in `docs/analysis_plan.md`

---

### Data Collection

**13. Data collection procedures:**

1. **API access**: Automated download via OCDS endpoints
2. **Date range**: 2015-2023 (varies by country)
3. **Fields extracted**: 42 OCDS fields (see `docs/OCDS_field_audit.md`)
4. **Quality filters**: Minimum 80% field completeness per record
5. **Entity resolution**: Fuzzy matching + manual validation (κ > 0.85)

**14. Timing:**

- Data downloaded: [TBD]
- Analysis completed: [TBD]
- No time-sensitive data (historical administrative records)

---

### Research Involving Human Participants

**15. Recruitment:**

☐ N/A - No human participants; administrative data only

**16. Consent:**

☐ N/A - Public administrative data; no individual consent required

---

### Sustainability Metrics (Study-Specific)

**17. Carbon footprint methodology:**

| Component | Source | Method |
|-----------|--------|--------|
| Sector classification | CPV → EXIOBASE | Author crosswalk |
| Emissions intensity | EXIOBASE 3 | Stadler et al. (2018) |
| Uncertainty | Monte Carlo | 1000 draws from intensity distributions |

**18. Resilience metrics:**

| Metric | Definition | Source |
|--------|------------|--------|
| Supplier HHI | Herfindahl-Hirschman Index | Computed from contract data |
| Single-bid rate | % contracts with 1 bidder | OCDS field |
| Completion rate | % contracts completed | OCDS status field |

**19. Green procurement classification:**

- CPV code matching (EU GPP criteria)
- NLP keyword detection (renewable, solar, etc.)
- Inter-annotator agreement: κ > 0.80 (target)

---

### AI/Machine Learning Disclosure

**20. Use of AI tools:**

☐ LLMs used for code documentation assistance
☐ No LLM-generated content in manuscript text
☐ All AI assistance disclosed per Nature policy
☐ Human verification of all AI-assisted outputs

---

## Checklist Summary

### Required for Submission

- [ ] Completed Reporting Summary
- [ ] Data availability statement
- [ ] Code availability statement
- [ ] Statistics section complete
- [ ] Reproducibility documentation
- [ ] Ethics statement (if applicable)

### Files to Include

1. `reporting_summary.pdf` - This document
2. `data_availability.md` - Data sources and access
3. `code_availability.md` - Code repository and DOI
4. `supplementary_information.pdf` - Extended methods

---

*Document Version: 1.0*
*Last Updated: [Date]*
*Corresponding Author: Jack00040008@outlook.com*
