#!/usr/bin/env python3
"""Debug ANGPTL3 locus L2G predictions."""

import pandas as pd
import pyarrow.parquet as pq
from pathlib import Path

bench = pd.read_csv('data/processed/baselines/independent_benchmark_v1.tsv', sep='\t')
angptl3 = bench[bench['gene_symbol'] == 'ANGPTL3'].iloc[0]
chrom = '1'
pos = int(angptl3['pos_hg38'])  # 62597513
true_gene_id = angptl3['gene_id']  # ENSG00000132855

print(f'Looking for L2G predictions near chr1:{pos} (ANGPTL3 = {true_gene_id})')
print()

# Search all parquet files
l2g_dir = Path('data/external/open_targets/l2g_bulk')
parquet_files = sorted(l2g_dir.glob('*.parquet'))

all_nearby = []
for f in parquet_files:
    try:
        df = pq.read_table(f).to_pandas()
        df['chrom_clean'] = df['chrom'].astype(str).str.replace('chr', '')
        nearby = df[(df['chrom_clean'] == chrom) & (df['pos'].between(pos - 500000, pos + 500000))]
        if len(nearby) > 0:
            all_nearby.append(nearby)
    except Exception as e:
        pass

if all_nearby:
    combined = pd.concat(all_nearby).sort_values('y_proba_full_model', ascending=False)
    print(f'Total nearby: {len(combined)}')
    print()
    print('Top 10 predictions:')
    for i, row in combined.head(10).iterrows():
        is_true = 'TRUE' if row['gene_id'] == true_gene_id else '    '
        print(f"  {is_true} {row['gene_id']} score={row['y_proba_full_model']:.3f}")
    
    # Check if true gene appears anywhere
    has_true = combined[combined['gene_id'] == true_gene_id]
    if len(has_true) > 0:
        print(f"\nTrue gene {true_gene_id} appears at rank:")
        for rank, (i, row) in enumerate(combined.iterrows(), 1):
            if row['gene_id'] == true_gene_id:
                print(f"  Rank {rank}, score={row['y_proba_full_model']:.3f}")
                break
    else:
        print(f"\nTrue gene {true_gene_id} NOT FOUND in L2G predictions!")
else:
    print('No L2G predictions found near this locus!')
