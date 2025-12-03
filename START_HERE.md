# START HERE: Quick Guide for Editors and Reviewers

**Paper:** Path-Probability Models Outperform Point-Estimate Scores for Noncoding GWAS Gene Prioritization

**Author:** Abduxoliq Ashuraliyev | Jack00040008@outlook.com

**Resources:**
- **GitHub:** https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs
- **Zenodo:** https://doi.org/10.5281/zenodo.17798899
- **License:** MIT (open, permissive)

---

## The One-Sentence Summary

We replace single-number gene scores with **mechanism graphs**—probabilistic paths from variants → enhancers → genes → tissues → traits—and show this improves prioritization accuracy (+18% recall) while providing calibrated probabilities at every step.

---

## Key Claims and Where to Verify Them

| Claim | Evidence Location | How to Verify |
|-------|-------------------|---------------|
| 76% recall at rank 20 (vs 58% L2G) | Fig 3a, `results/benchmarks/` | `snakemake results/figures/fig3_benchmark.pdf` |
| ECE < 0.05 per module | Table 1, Fig 4 | `snakemake results/figures/fig4_calibration.pdf` |
| r = 0.89 eQTL Catalogue correlation | Fig 5b | `snakemake results/figures/fig5_replication.pdf` |
| ABC+PCHi-C AUPRC 0.71 vs 0.54 distance | Fig 2a | `snakemake results/figures/fig2_bridge.pdf` |
| 78% cross-study replication | Fig 5a | `results/replication/summary.tsv` |

---

## Reproduce Figures 1-2 in 15 Minutes

```bash
# 1. Clone
git clone https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs.git
cd Mechanism-GWAS-Causal-Graphs

# 2. Setup (requires conda)
conda env create -f environment.yml
conda activate mechanism-gwas

# 3. Download test data (5 min)
snakemake --cores 1 download_test_data

# 4. Generate figures (10 min)
snakemake --cores 8 results/figures/fig1_overview.pdf results/figures/fig2_bridge.pdf
```

---

## What Makes This Different from Existing Methods

| Aspect | Open Targets L2G / PoPS | This Work |
|--------|-------------------------|-----------|
| Output | Single score per gene | Full probability path |
| Calibration | Not calibrated | ECE < 0.05 per module |
| Mechanism | Collapsed into score | Explicit: variant → cCRE → gene → tissue |
| Multi-signal | Single-causal assumption | coloc.susie handles multiple signals |
| Enhancer-gene | Distance-weighted | ABC + PCHi-C ensemble |
| Validation | Overlapping benchmarks | Anti-leak tiers with provenance |

---

## Benchmark Anti-Leak Verification

Tier 1 benchmark genes were verified absent from training data of:
- Open Targets L2G v22.09
- PoPS
- MAGMA

Provenance documented in `data/manifests/benchmark_genes.yaml` with:
- OMIM accession numbers
- PMIDs for each gene
- Curation date
- Anti-overlap verification script: `scripts/verify_antileak.py`

---

## Directory Structure

```
Mechanism-GWAS-Causal-Graphs/
├── src/                      # Python source code
│   ├── finemapping/         # SuSiE wrapper
│   ├── colocalization/      # coloc.susie interface
│   ├── enhancer_gene_linking/  # ABC + PCHi-C
│   ├── mechanism_graph/     # Noisy-OR inference
│   └── calibration/         # ECE computation
├── workflow/
│   └── Snakefile            # Reproducible pipeline
├── data/manifests/          # Data sources with checksums
├── manuscript/
│   ├── main.tex             # Full manuscript
│   └── references.bib       # Bibliography
├── results/                 # Generated outputs
├── REPRODUCE.md             # Full reproduction guide
└── environment.yml          # Conda dependencies
```

---

## For Any Questions

- **Issues:** https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs/issues
- **Email:** Jack00040008@outlook.com

---

*This guide was generated to facilitate editorial review. All materials are publicly available without login or access restrictions.*
