# RegulatoryBench: Task-Stratified Benchmark for Regulatory Element-to-Gene Linking

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Overview

RegulatoryBench is a task-stratified evaluation framework for regulatory genomics methods. It separates three distinct prediction tasks that have been conflated in prior benchmarking efforts:

- **Task A**: GWAS Credible Set → Causal Gene
- **Task B**: Regulatory Element → Target Gene  
- **Task C**: Regulatory Variant → Affected Gene

## Key Findings

1. **Distance dominates GWAS-to-gene mapping (Task A)**: Simple distance ranking achieves AUROC=0.93, vastly outperforming all enhancer-gene methods including ABC (0.60) and PoPS (0.79).

2. **ABC performs comparably to distance on E2G (Task B)**: AUROC=0.88 vs 0.87, but evaluation is confounded by 78% training data leakage.

3. **Systematic leakage plagues published benchmarks**: We identify three types of contamination affecting the majority of benchmark pairs.

## Repository Structure

```
regulatorybench/
├── task_taxonomy.py          # Task definitions and evaluation contract
├── split_benchmarks.py       # Benchmark splitter by task type
├── evaluate_task_a.py        # Task A evaluation (23 methods)
├── evaluate_task_b.py        # Task B evaluation
├── leakage_audit.py          # Leakage detection framework
├── generate_figures.py       # Publication figure generator
├── MANUSCRIPT_DRAFT.md       # Nature Genetics Analysis draft
├── RESULTS_SUMMARY.md        # Key findings summary
├── benchmarks/
│   ├── task_a_gwas_to_gene.parquet       # 14,016 pairs
│   ├── task_b_enhancer_to_gene.parquet   # 19,825 pairs
│   ├── task_a_evaluation_results.json
│   ├── task_b_baseline_results.json
│   ├── task_a_leakage_audit.json
│   └── task_b_leakage_audit.json
└── figures/
    ├── figure1_task_taxonomy.pdf
    ├── figure2_task_a_results.pdf
    ├── figure3_task_b_results.pdf
    ├── figure4_leakage_audit.pdf
    ├── table1_method_comparison.csv
    └── table2_leakage_controls.csv
```

## Quick Start

```bash
# Clone repository
git clone https://github.com/[user]/regulatorybench.git
cd regulatorybench

# Install dependencies
pip install pandas numpy scikit-learn matplotlib seaborn pyarrow

# Run full evaluation
python evaluate_task_a.py
python evaluate_task_b.py
python leakage_audit.py
python generate_figures.py
```

## Benchmark Data

### Task A: GWAS → Gene (14,016 pairs)

| Statistic | Value |
|-----------|-------|
| Total pairs | 14,016 |
| Unique loci | 560 |
| Positives | 569 (4.1%) |
| Source | UK Biobank via E2G framework |

### Task B: Enhancer → Gene (19,825 pairs)

| Statistic | Value |
|-----------|-------|
| Total pairs | 19,825 |
| Unique loci | 6,927 |
| Positives | 863 (4.4%) |
| Sources | ENCODE CRISPRi, Fulco 2019 |

## Citation

If you use RegulatoryBench in your research, please cite:

```bibtex
@article{regulatorybench2024,
  title={RegulatoryBench: A Task-Stratified Benchmark Reveals Distance Dominance in GWAS-to-Gene Mapping},
  author={[Authors]},
  journal={Nature Genetics},
  year={2024},
  doi={10.5281/zenodo.XXXXXXX}
}
```

## License

MIT License

## Contact

[Contact information]
