#!/usr/bin/env python3
"""
Analyze the ENCODE CRISPRi benchmark data for RegulatoryBench v3.
"""

import gzip
import pandas as pd
from pathlib import Path
import os

os.chdir(Path(__file__).parent.parent)

# Load training data
data_path = Path("data/external/crispr_benchmark/resources/crispr_data")
training_file = data_path / "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
heldout_file = data_path / "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz"

print("="*60)
print("ENCODE CRISPRi BENCHMARK ANALYSIS")
print("="*60)

# Load training data
with gzip.open(training_file, 'rt') as f:
    df_train = pd.read_csv(f, sep='\t')

print(f"\n=== K562 Training Dataset ===")
print(f"Total element-gene pairs: {len(df_train)}")
print(f"Regulated (positive): {df_train['Regulated'].sum()}")
print(f"Not regulated (negative): {(~df_train['Regulated']).sum()}")

# Load heldout data
with gzip.open(heldout_file, 'rt') as f:
    df_heldout = pd.read_csv(f, sep='\t')

print(f"\n=== Heldout 5 Cell Types Dataset ===")
print(f"Total element-gene pairs: {len(df_heldout)}")
print(f"Regulated (positive): {df_heldout['Regulated'].sum()}")
print(f"Not regulated (negative): {(~df_heldout['Regulated']).sum()}")

# Combined stats
df_all = pd.concat([df_train, df_heldout], ignore_index=True)
print(f"\n=== Combined Dataset ===")
print(f"Total element-gene pairs: {len(df_all)}")
print(f"Regulated (positive): {df_all['Regulated'].sum()}")
print(f"Unique genes: {df_all['measuredGeneSymbol'].nunique()}")
print(f"Unique Ensembl IDs: {df_all['measuredGeneEnsemblId'].nunique()}")

print(f"\n=== Cell Type Distribution ===")
print(df_all['CellType'].value_counts())

print(f"\n=== Dataset Sources ===")
print(df_all['Dataset'].value_counts())

# Effect size distribution for positives
positives = df_all[df_all['Regulated'] == True]
print(f"\n=== Effect Size for Positives ===")
print(f"Min: {positives['EffectSize'].min():.3f}")
print(f"Mean: {positives['EffectSize'].mean():.3f}")
print(f"Median: {positives['EffectSize'].median():.3f}")
print(f"Max: {positives['EffectSize'].max():.3f}")

# Distance distribution
print(f"\n=== Distance to TSS ===")
print(f"Min: {positives['distanceToTSS'].min()}")
print(f"Mean: {positives['distanceToTSS'].mean():.0f}")
print(f"Median: {positives['distanceToTSS'].median():.0f}")
print(f"Max: {positives['distanceToTSS'].max()}")

# Show sample of positive pairs
print(f"\n=== Sample of Positive Regulated Pairs ===")
sample_cols = ['measuredGeneSymbol', 'measuredGeneEnsemblId', 'CellType', 
               'EffectSize', 'distanceToTSS', 'Dataset']
print(positives[sample_cols].head(20).to_string())

# Save the positive pairs for benchmark construction
output_path = Path("data/processed/crispr_positives.tsv")
positives.to_csv(output_path, sep='\t', index=False)
print(f"\n✓ Saved {len(positives)} positive pairs to {output_path}")

# Count unique validated target genes - this is what we need for RegulatoryBench
print(f"\n=== SUMMARY FOR REGULATORYBENCH V3 ===")
unique_genes = positives['measuredGeneEnsemblId'].nunique()
unique_symbols = positives['measuredGeneSymbol'].nunique()
print(f"Unique CRISPRi-validated target genes: {unique_genes} Ensembl IDs / {unique_symbols} symbols")
print(f"These are EXPERIMENTAL GROUND TRUTH with direct perturbation evidence")
