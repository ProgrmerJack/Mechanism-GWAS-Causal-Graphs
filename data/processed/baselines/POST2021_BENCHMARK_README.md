# Post-2021 Independent GWAS Benchmark

**Version**: 1.0
**Created**: 2025-12-23
**Purpose**: Fair comparison of gene prioritization methods

## Why This Benchmark?

The standard Open Targets gold standards (gwas_gold_standards.191108.tsv)
were used to **train the L2G model** (Mountjoy et al. 2021). Using them
for benchmarking L2G creates 100% data leakage.

This benchmark contains **independent loci** from publications after 2021,
ensuring fair comparison.

## Curation Criteria

**Inclusion Criteria**:
1. Publication date: 2022 or later
2. Genome-wide significant association (P < 5e-8)
3. High-confidence causal gene identification via:
   - Mendelian disease (ClinGen Definitive/Strong)
   - FDA-approved drug target
   - CRISPR/functional validation
   - Coding variant in gene
   - Multiple lines of converging evidence

**Exclusion Criteria**:
1. Publication before 2022
2. Ambiguous gene assignment
3. Only eQTL evidence (no functional validation)

## Evidence Tiers

- **Tier1_CRISPR**: CRISPR screen or base editing validation
- **Tier1_Mendelian**: Monogenic disease with ClinGen evidence
- **Tier1_Drug**: FDA-approved drug target
- **Tier1_Coding**: Coding variant with functional impact
- **Tier2_MultiEvidence**: Multiple converging lines of evidence

## Benchmark Statistics

- **Total loci**: 17
- **Unique genes**: 17
- **Trait categories**: 5
- **L2G training overlap**: 0 loci (verified independent)

## Data Sources

1. **GWAS Catalog** (NHGRI-EBI)
2. **ClinGen** Gene-Disease Validity Curations
3. **Recent Publications**:
   - Graham et al. 2022 Nature (CAD)
   - Mahajan et al. 2022 Nat Genet (T2D)
   - Koyama et al. 2024 Nat Genet (CAD)
   - Zoodsma et al. 2025 Nat Genet (Metabolism)

## Files

- `post2021_independent_benchmark.tsv`: Main benchmark file
- `POST2021_BENCHMARK_README.md`: This file

## Column Descriptions

- `locus_id`: Unique locus identifier
- `chr`: Chromosome
- `pos_hg38`: Position (GRCh38)
- `lead_snp`: Lead GWAS variant (rsID)
- `gene_symbol`: Causal gene symbol
- `gene_id`: Ensembl gene ID
- `trait`: GWAS trait
- `trait_category`: Trait category (Cardiovascular, Metabolic, Immune)
- `gwas_pmid`: PubMed ID of GWAS publication
- `gwas_pvalue`: GWAS P-value
- `evidence_tier`: Evidence strength (Tier1/Tier2)
- `validation_type`: Type of functional validation
- `validation_pmid`: PubMed ID of validation study
- `curated_date`: Date of manual curation
- `curator`: Curator name
- `notes`: Additional notes

## Usage

```python
import pandas as pd

# Load benchmark
benchmark = pd.read_csv('post2021_independent_benchmark.tsv', sep='\t')

# Filter by evidence tier
tier1_only = benchmark[benchmark['evidence_tier'].str.startswith('Tier1')]
```

---
**Generated**: 2025-12-23 16:48:30