#!/usr/bin/env python3
"""Analyze how positive labels are defined in Task A benchmark."""

import pandas as pd
from pathlib import Path

# Load benchmark
df = pd.read_parquet('benchmarks/task_a_gwas_to_gene.parquet')
pos = df[df.truth == True]

print("=" * 80)
print("TASK A POSITIVE LABEL ANALYSIS")
print("=" * 80)
print(f"\nTotal pairs: {len(df)}")
print(f"Positive pairs: {len(pos)} ({len(pos)/len(df)*100:.2f}%)")
print(f"Unique loci: {df.locus_id.nunique()}")

print("\n" + "=" * 80)
print("EVIDENCE BREAKDOWN FOR 569 POSITIVES")
print("=" * 80)

print("\n1. CODING EVIDENCE:")
print(f"   AnyCoding: {pos.AnyCoding.sum()} ({pos.AnyCoding.sum()/len(pos)*100:.1f}%)")
print(f"   AnySpliceSite: {pos.AnySpliceSite.sum()} ({pos.AnySpliceSite.sum()/len(pos)*100:.1f}%)")
print(f"   AnyPromoter: {pos.AnyPromoter.sum()} ({pos.AnyPromoter.sum()/len(pos)*100:.1f}%)")

print("\n2. ABC PREDICTION:")
print(f"   MaxABC > 0: {(pos.MaxABC > 0).sum()} ({(pos.MaxABC > 0).sum()/len(pos)*100:.1f}%)")
print(f"   MaxABC > 0.015: {(pos.MaxABC > 0.015).sum()} ({(pos.MaxABC > 0.015).sum()/len(pos)*100:.1f}%)")
print(f"   ABCPrediction (any): {(pos.ABCPrediction > 0).sum()} ({(pos.ABCPrediction > 0).sum()/len(pos)*100:.1f}%)")

print("\n3. DISTANCE:")
print(f"   DistanceRank == 1 (nearest gene): {(pos.DistanceRank == 1).sum()} ({(pos.DistanceRank == 1).sum()/len(pos)*100:.1f}%)")
print(f"   IsClosestGeneToBestSNP: {pos.IsClosestGeneToBestSNP.sum()} ({pos.IsClosestGeneToBestSNP.sum()/len(pos)*100:.1f}%)")

print("\n4. CORRECT LABEL COLUMNS (these indicate whether each method got it right):")
abc_cols = ['Correct.ABCAny', 'Correct.FMEnriched', 'Correct.LDSCEnriched', 'Correct.Distance']
for col in abc_cols:
    count = pos[col].sum()
    print(f"   {col}: {count} ({count/len(pos)*100:.1f}%)")

print("\n5. OVERLAP ANALYSIS:")
abc_pos = pos[pos.MaxABC > 0.015]
dist_pos = pos[pos.DistanceRank == 1]
both = pos[(pos.MaxABC > 0.015) & (pos.DistanceRank == 1)]
abc_only = pos[(pos.MaxABC > 0.015) & (pos.DistanceRank != 1)]
dist_only = pos[(pos.MaxABC <= 0.015) & (pos.DistanceRank == 1)]
neither = pos[(pos.MaxABC <= 0.015) & (pos.DistanceRank != 1)]

print(f"   Both ABC>0.015 AND DistanceRank=1: {len(both)} ({len(both)/len(pos)*100:.1f}%)")
print(f"   ABC>0.015 only: {len(abc_only)} ({len(abc_only)/len(pos)*100:.1f}%)")
print(f"   DistanceRank=1 only: {len(dist_only)} ({len(dist_only)/len(pos)*100:.1f}%)")
print(f"   Neither ABC nor Distance: {len(neither)} ({len(neither)/len(pos)*100:.1f}%)")

print("\n6. WHAT DEFINES THESE POSITIVES?")
print("   Checking if truth==True requires ANY of these conditions...")

# Check correlation
print(f"\n   Correlation analysis:")
print(f"   - truth==True & MaxABC>0: {((pos.MaxABC > 0).sum())} / {len(pos)}")
print(f"   - truth==True & DistanceRank==1: {(pos.DistanceRank == 1).sum()} / {len(pos)}")
print(f"   - truth==True & (coding OR ABC OR distance=1): {((pos.AnyCoding | (pos.MaxABC > 0.015) | (pos.DistanceRank == 1)).sum())} / {len(pos)}")

# Look for external validation
if 'ABCTrainingExample' in pos.columns:
    print(f"\n7. ABC TRAINING CONTAMINATION:")
    print(f"   Positives that are ABC training examples: {pos.ABCTrainingExample.sum()} ({pos.ABCTrainingExample.sum()/len(pos)*100:.1f}%)")

print("\n" + "=" * 80)
print("CRITICAL FINDING:")
print("=" * 80)
print(f"""
The positive labels contain {(pos.MaxABC > 0.015).sum()} pairs ({(pos.MaxABC > 0.015).sum()/len(pos)*100:.1f}%) 
where ABC score > 0.015.

This means the benchmark uses ABC predictions as part of the gold standard,
creating circularity when evaluating ABC performance.

For a truly independent benchmark, we need labels based ONLY on:
1. Coding/LoF variants mapping unambiguously to genes
2. Rare variant burden tests (ExWAS gene hits)
3. Experimental perturbation data (CRISPRi at GWAS loci)
4. Clinical genetic evidence (OMIM, ClinVar, Gene2Phenotype)
""")
