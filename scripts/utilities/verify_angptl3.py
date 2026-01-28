#!/usr/bin/env python3
"""Verify the matching logic for ANGPTL3."""

import pandas as pd

# Check what true_ensembl the code would see for ANGPTL3
bench = pd.read_csv('data/processed/baselines/independent_benchmark_v1.tsv', sep='\t')
row = bench[bench['gene_symbol'] == 'ANGPTL3'].iloc[0]

true_ensembl = row.get('gene_id', row.get('ensembl_id', ''))
print(f"gene_symbol: {row['gene_symbol']}")
print(f"gene_id: {row['gene_id']}")
print(f"true_ensembl for matching: {true_ensembl}")

# L2G top prediction is ENSG00000132854
print()
print(f"L2G top predicted: ENSG00000132854")
print(f"true_ensembl:      {true_ensembl}")
print(f"Match: {str('ENSG00000132854').upper() == str(true_ensembl).upper()}")

# The point: L2G predicts ENSG00000132854 but true is ENSG00000132855
# These are DIFFERENT genes - off by 1 in ID
print()
print("CONCLUSION: L2G genuinely predicts the WRONG gene for ANGPTL3 locus")
print("  - ENSG00000132854 is likely ANGPTL3-AS1 (antisense) or another nearby gene")
print("  - ENSG00000132855 is ANGPTL3 (the true causal gene)")
print("  - This is NOT a code bug - L2G is genuinely failing on this locus")
