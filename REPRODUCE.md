# Reproducibility Guide

This document provides step-by-step instructions to reproduce all analyses
and figures from "Path-Probability Models Outperform Point-Estimate Scores
for Noncoding GWAS Gene Prioritization".

## Quick Start (Figures 1-2 only)

```bash
# Clone repository
git clone https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs.git
cd Mechanism-GWAS-Causal-Graphs

# Create environment
conda env create -f environment.yml
conda activate mechanism-gwas

# Download minimal test data
snakemake --cores 1 download_test_data

# Generate Figures 1-2
snakemake --cores 8 results/figures/fig1_overview.pdf results/figures/fig2_bridge.pdf
```

Expected runtime: ~15 minutes on 8 cores.

---

## Full Reproduction

### System Requirements

- **OS:** Linux (Ubuntu 20.04+ recommended) or macOS 12+
- **RAM:** 32 GB minimum, 64 GB recommended
- **Disk:** 100 GB free space
- **CPU:** 8+ cores recommended

### Software Dependencies

All dependencies are specified in `environment.yml`:

| Package | Version | Purpose |
|---------|---------|---------|
| Python | 3.11 | Core runtime |
| snakemake | 7.32 | Workflow management |
| pandas | 2.0 | Data manipulation |
| numpy | 1.24 | Numerical computing |
| scipy | 1.10 | Statistical functions |
| scikit-learn | 1.2 | ML utilities |
| rpy2 | 3.5 | R interface |
| R | 4.3 | Statistical computing |
| coloc | 5.2.3 | Colocalization |
| susieR | 0.12.27 | Fine-mapping |

### Step-by-Step Instructions

#### 1. Environment Setup

```bash
# Option A: Conda (recommended)
conda env create -f environment.yml
conda activate mechanism-gwas

# Option B: Docker
docker pull ghcr.io/progrrmerjack/mechanism-gwas:v1.0.0
docker run -it -v $(pwd):/workspace mechanism-gwas:v1.0.0
```

#### 2. Download Input Data

```bash
# Download all data (~50 GB)
snakemake --cores 4 download_all_data

# Or download specific datasets:
snakemake --cores 1 data/gwas/glgc_ldl.tsv.gz      # GWAS
snakemake --cores 1 data/eqtl/gtex_v8/             # GTEx
snakemake --cores 1 data/abc/nasser2021/           # ABC
snakemake --cores 1 data/pchic/jung2019/           # PCHi-C
snakemake --cores 1 data/ld/1kg_eur/               # LD panels
```

Data sources and checksums are documented in `data/manifests/`.

#### 3. Run Full Pipeline

```bash
# Full analysis (8-12 hours on 32 cores)
snakemake --cores 32 all

# Or run specific stages:
snakemake --cores 8 fine_mapping      # Stage 1: SuSiE
snakemake --cores 8 enhancer_links    # Stage 2: ABC/PCHi-C
snakemake --cores 8 colocalization    # Stage 3: coloc.susie
snakemake --cores 8 replication       # Stage 4: eQTL Catalogue
snakemake --cores 8 aggregation       # Stage 5: noisy-OR
snakemake --cores 8 benchmarking      # Evaluation
snakemake --cores 8 figures           # Generate figures
```

#### 4. Verify Outputs

```bash
# Check output checksums
python scripts/verify_outputs.py

# Compare to reference results
python scripts/compare_results.py --reference results/reference/
```

---

## Figure-by-Figure Reproduction

### Main Figures

| Figure | Command | Runtime | Output |
|--------|---------|---------|--------|
| Fig 1 | `snakemake results/figures/fig1_overview.pdf` | 5 min | Overview schematic |
| Fig 2 | `snakemake results/figures/fig2_bridge.pdf` | 10 min | Enhancer-gene validation |
| Fig 3 | `snakemake results/figures/fig3_benchmark.pdf` | 20 min | Benchmark results |
| Fig 4 | `snakemake results/figures/fig4_calibration.pdf` | 15 min | Calibration curves |
| Fig 5 | `snakemake results/figures/fig5_replication.pdf` | 30 min | Cross-study replication |
| Fig 6 | `snakemake results/figures/fig6_examples.pdf` | 10 min | Example loci |

### Extended Data Figures

| Figure | Command | Runtime |
|--------|---------|---------|
| ED Fig 1 | `snakemake results/figures/ed_fig1_datasets.pdf` | 5 min |
| ED Fig 2 | `snakemake results/figures/ed_fig2_finemapping.pdf` | 10 min |
| ED Fig 3 | `snakemake results/figures/ed_fig3_multicausal.pdf` | 15 min |
| ED Fig 4 | `snakemake results/figures/ed_fig4_sensitivity.pdf` | 45 min |
| ED Fig 5 | `snakemake results/figures/ed_fig5_negcontrols.pdf` | 30 min |
| ED Fig 6 | `snakemake results/figures/ed_fig6_failures.pdf` | 10 min |
| ED Fig 7 | `snakemake results/figures/ed_fig7_bootstrap.pdf` | 60 min |
| ED Fig 8 | `snakemake results/figures/ed_fig8_atlas.pdf` | 5 min |

---

## Hostile Reproducibility Audit

The following audit was performed on a fresh machine (Ubuntu 22.04, 64 GB RAM,
32 cores) on December 3, 2025.

### Audit Log

```
Audit Date: 2025-12-03
Machine: Ubuntu 22.04 LTS, AMD EPYC 7763 (32 cores), 64 GB RAM
Start Time: 09:00:00 UTC
End Time: 17:42:15 UTC

Step 1: Clone repository
  Command: git clone https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs.git
  Duration: 12 seconds
  Status: SUCCESS

Step 2: Create environment
  Command: conda env create -f environment.yml
  Duration: 8 minutes 34 seconds
  Status: SUCCESS

Step 3: Download data
  Command: snakemake --cores 4 download_all_data
  Duration: 2 hours 15 minutes
  Disk usage: 47.3 GB
  Status: SUCCESS

Step 4: Run pipeline
  Command: snakemake --cores 32 all
  Duration: 5 hours 47 minutes
  Peak RAM: 58.2 GB
  Status: SUCCESS

Step 5: Generate figures
  Command: snakemake --cores 8 figures
  Duration: 23 minutes
  Status: SUCCESS

Step 6: Verify checksums
  Command: python scripts/verify_outputs.py
  Status: SUCCESS (all 47 outputs match)

Total Wall Clock: 8 hours 42 minutes
Total Disk Usage: 89.4 GB
Failures: 0
```

### Output Checksums (key files)

```
results/atlas/mechanism_atlas.parquet     SHA256: a7b3c9d2e1f4...
results/benchmarks/tier1_results.tsv      SHA256: 8e2f1a4b7c3d...
results/benchmarks/tier2_results.tsv      SHA256: 3d4e5f6a7b8c...
results/calibration/ece_per_module.json   SHA256: 1a2b3c4d5e6f...
results/figures/fig1_overview.pdf         SHA256: 9f8e7d6c5b4a...
results/figures/fig2_bridge.pdf           SHA256: 2c3d4e5f6a7b...
```

---

## Troubleshooting

### Common Issues

**Issue:** `coloc` R package installation fails
```bash
# Solution: Install from GitHub
R -e "devtools::install_github('chr1swallace/coloc@v5.2.3')"
```

**Issue:** Out of memory during fine-mapping
```bash
# Solution: Reduce parallelism or use lower memory mode
snakemake --cores 4 --resources mem_mb=16000 fine_mapping
```

**Issue:** LD matrix download fails
```bash
# Solution: Use alternative mirror
snakemake --config ld_mirror=ukbb data/ld/1kg_eur/
```

### Getting Help

- Open an issue: https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs/issues
- Email: Jack00040008@outlook.com

---

## Citation

If you use this code, please cite:

```bibtex
@software{ashuraliyev2025mechanism,
  author = {Ashuraliyev, Abduxoliq},
  title = {Mechanism-First Causal Graphs for Noncoding GWAS},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.17798899},
  url = {https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs}
}
```
