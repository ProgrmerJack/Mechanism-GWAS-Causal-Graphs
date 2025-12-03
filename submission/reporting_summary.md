# Nature Genetics Reporting Summary

This document addresses all items required by the Nature Research Reporting Summary
for computational studies.

## Statistical Methods

### Sample sizes
- **GWAS loci analyzed:** 2,134 genome-wide significant loci across 8 cardiometabolic traits
- **Benchmark genes (Tier 1):** 47 Mendelian genes with verified absence from training sets
- **Benchmark genes (Tier 2):** 89 approved drug targets
- **CRISPRi validation pairs:** 847 enhancer-gene pairs

### Statistical tests used
- **Calibration:** Expected Calibration Error (ECE) with 10 equal-width bins
- **Correlation:** Pearson correlation coefficient (r) for effect size comparisons
- **Confidence intervals:** 1,000 bootstrap replicates, 2.5th-97.5th percentiles
- **Multiple testing:** Benjamini-Hochberg FDR correction where applicable

### Reproducibility
- All analyses are reproducible via Snakemake workflow
- Random seeds fixed for reproducibility: `numpy.random.seed(42)`
- Full audit log provided in REPRODUCE.md

## Data Availability

### Data sources and accessibility
| Dataset | Source | Accession/URL | Access Type |
|---------|--------|---------------|-------------|
| GLGC GWAS | Public | http://csg.sph.umich.edu/willer/public/lipids2013/ | Open |
| CARDIoGRAMplusC4D | Public | http://www.cardiogramplusc4d.org/ | Open |
| DIAGRAM T2D | Public | https://diagram-consortium.org/ | Open |
| GTEx v8 eQTLs | dbGaP | phs000424.v8.p2 | Open summary statistics |
| eQTL Catalogue | EBI | https://www.ebi.ac.uk/eqtl/ | Open |
| ABC Model | Engreitz Lab | https://www.engreitzlab.org/resources/ | Open |
| PCHi-C (Jung) | GEO | GSE118752 | Open |
| PCHi-C (Javierre) | OSF | https://osf.io/u8tzp/ | Open |
| 1000 Genomes | IGSR | https://www.internationalgenome.org/ | Open |
| GENCODE v40 | GENCODE | https://www.gencodegenes.org/ | Open |

### No controlled-access data
All analyses use publicly available summary-level data.
No individual-level genotypes or controlled-access datasets were used.

### Data manifests
Versioned data manifests with SHA256 checksums for all input files
are provided in `data/manifests/`.

## Code Availability

- **GitHub:** https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs
- **Zenodo:** DOI: 10.5281/zenodo.17798899
- **License:** MIT (permissive, allows commercial use)

### Software versions
| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.11 | Core runtime |
| R | 4.3 | Statistical computing |
| coloc | 5.2.3 | Colocalization |
| susieR | 0.12.27 | Fine-mapping |
| snakemake | 7.32 | Workflow |
| pandas | 2.0 | Data manipulation |
| numpy | 1.24 | Numerical computing |
| scipy | 1.10 | Statistical functions |
| scikit-learn | 1.2 | ML utilities |

## Methods

### Fine-mapping
- Tool: SuSiE-RSS (susieR v0.12.27)
- Parameters: L=10 signals, 95% credible set coverage
- LD reference: 1000 Genomes Phase 3 European (n=503)

### Colocalization
- Tool: coloc.susie (coloc v5.2.3)
- Priors: p1=10^-4, p2=10^-4, p12=5×10^-6
- Threshold: PP.H4 ≥ 0.8 for colocalization

### Enhancer-gene linking
- ABC: Score threshold ≥ 0.015, 131 biosamples
- PCHi-C: CHiCAGO score ≥ 5, 44 cell types
- Ensemble: Weighted logistic regression, cross-validated

### Calibration assessment
- Metric: Expected Calibration Error (ECE)
- Bins: 10 equal-width probability bins
- Confidence: 1,000 bootstrap replicates

## Benchmark Design

### Anti-leak provisions (Tier 1)
- Genes verified absent from Open Targets L2G v22.09 training set
- Genes verified absent from PoPS training set
- 500 kb exclusion zone around any training gene
- Two independent curators verified inclusion/exclusion
- Full provenance in `data/manifests/benchmark_genes.yaml`

### Drug target benchmark (Tier 2)
- Source: ChEMBL v32, max_phase = 4
- Indication: EFO cardiovascular/metabolic
- Filtered to exclude L2G training overlap

### CRISPR validation (Tier 3)
- Fulco et al. (2019): 664 pairs in K562
- Gasperini et al. (2019): 183 pairs in K562/WTC11
- Independent of genetic association

## Ethics and Competing Interests

- **Ethics approval:** Not applicable (publicly available summary statistics only)
- **Competing interests:** None declared
- **Funding:** None

## Author Contributions

A.A. conceived the study, developed the methodology, performed analyses,
and wrote the manuscript.

---

*Completed by: Abduxoliq Ashuraliyev*
*Date: December 3, 2025*
