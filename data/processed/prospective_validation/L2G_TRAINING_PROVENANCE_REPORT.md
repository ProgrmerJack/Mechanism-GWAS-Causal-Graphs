
# L2G Training Data Provenance Verification Report

## Summary
- **STING-seq benchmark genes**: 128
- **L2G gold standard genes**: 526
- **Overlap**: 29 genes

## Conclusion

⚠️ **POTENTIAL TRAINING LEAKAGE DETECTED**

29 genes overlap between STING-seq benchmark and L2G training data.

These genes should be EXCLUDED from AUROC calculation or the validation
should be explicitly labeled as "retrospective" for these loci.

Overlapping genes:
- ENSG00000000971 (CFH)
- ENSG00000010704 (HFE)
- ENSG00000021826 (CPS1)
- ENSG00000090534 (THPO)
- ENSG00000091513 (TF)
- ENSG00000096968 (JAK2)
- ENSG00000113302 (IL12B)
- ENSG00000113721 (PDGFRB)
- ENSG00000119866 (BCL11A)
- ENSG00000130203 (APOE)
- ENSG00000134242 (PTPN22)
- ENSG00000134460 (IL2RA)
- ENSG00000135111 (TBX3)
- ENSG00000136997 (MYC)
- ENSG00000142208 (AKT1)
- ENSG00000146648 (EGFR)
- ENSG00000153162 (BMP6)
- ENSG00000159216 (RUNX1)
- ENSG00000160791 (CCR5)
- ENSG00000169442 (CD52)
- ENSG00000171522 (PTGER4)
- ENSG00000179348 (GATA2)
- ENSG00000180210 (F2)
- ENSG00000182578 (CSF1R)
- ENSG00000183873 (SCN5A)
- ENSG00000187045 (TMPRSS6)
- ENSG00000187266 (EPOR)
- ENSG00000198734 (F5)
- ENSG00000198851 (CD3E)

## Data Sources

### L2G Training Data (v22.09)
- OTG Gold Standards: https://github.com/opentargets/genetics-gold-standards
- ChEMBL drug-target pairs (Phase II-IV)
- ClinVar/UniProt genetic associations  
- Gene2Phenotype, PanelApp, ClinGen

### STING-seq Benchmark
- Morris et al. (2024) Cell Genomics
- 132 CRE-gene pairs from CRISPR perturbation
- Blood cell traits (K562, lymphocytes)

## Verification Method
1. Fetched gold standards from OTG GitHub repository
2. Extracted all Ensembl gene IDs from training data
3. Mapped STING-seq gene symbols to Ensembl IDs
4. Computed intersection to find overlap
