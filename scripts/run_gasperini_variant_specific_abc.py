#!/usr/bin/env python3
"""
Variant-Specific ABC Validation on Gasperini 2019 CRISPRi Screen

This script validates our path-probability framework using the large-scale
Gasperini 2019 CRISPRi screen (514+ significant enhancer-gene pairs).

Key insight: Gene-level ABC fails because it includes noise from unrelated loci.
Variant-specific ABC (enhancers overlapping the tested CRISPR target) works.

Author: Mechanism-GWAS-Causal-Graphs
Date: 2024
"""

import pandas as pd
import numpy as np
import gzip
import json
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
ABC_FILE = PROJECT_ROOT / "data/external/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
GASPERINI_FILE = PROJECT_ROOT / "data/external/crispr_validation/gasperini_2019_results.txt.gz"
L2G_FILE = PROJECT_ROOT / "data/external/opentargets_l2g/l2g_scores_v22.09.tsv.gz"
OUTPUT_DIR = PROJECT_ROOT / "data/processed/prospective_validation"

# Blood-relevant cell types (K562 is the Gasperini cell line)
BLOOD_CELL_TYPES = [
    'K562-Roadmap',
    'erythroblast-Corces2016', 
    'GM12878-Roadmap',
    'CD14-positive_monocyte-ENCODE',
    'CD14-positive_monocyte_treated_with_LPS_for_1_hour-Novakovic2016',
    'CD14-positive_monocyte_treated_with_BG_for_1_hour-Novakovic2016',
    'CD4-positive_helper_T_cell-ENCODE',
    'CD4-positive_helper_T_cell-Corces2016',
    'CD8-positive_alpha-beta_T_cell-ENCODE',
    'CD8-positive_alpha-beta_T_cell-Corces2016',
    'B_cell-ENCODE',
    'natural_killer_cell-ENCODE',
    'natural_killer_cell-Corces2016',
    'monocyte-ENCODE',
    'neutrophil-Corces2016',
    'megakaryocyte-Corces2016',
    'common_myeloid_progenitor-Corces2016',
    'hematopoietic_multipotent_progenitor_cell-Corces2016',
    'lymphocyte-Roadmap',
    'spleen-ENCODE',
    'spleen-Roadmap',
    'thymus-ENCODE'
]


def load_gasperini_data():
    """Load and filter Gasperini 2019 significant enhancer-gene pairs."""
    print("Loading Gasperini 2019 data...")
    df = pd.read_csv(gzip.open(GASPERINI_FILE, 'rt'), sep='\t', low_memory=False)
    
    # Filter for enhancer sites (DHS = DNase hypersensitive sites = enhancers)
    enhancers = df[df['site_type'] == 'DHS'].copy()
    print(f"  Total enhancer pairs: {len(enhancers)}")
    
    # Get significant hits
    def is_significant(x, threshold=0.05):
        try:
            return float(x) < threshold
        except:
            return False
    
    significant = enhancers[enhancers['pvalue.empirical.adjusted'].apply(lambda x: is_significant(x, 0.05))].copy()
    print(f"  Significant at padj < 0.05: {len(significant)}")
    
    # Also get looser threshold for more power
    significant_10 = enhancers[enhancers['pvalue.empirical.adjusted'].apply(lambda x: is_significant(x, 0.10))].copy()
    print(f"  Significant at padj < 0.10: {len(significant_10)}")
    
    return enhancers, significant, significant_10


def load_abc_blood_matched():
    """Load ABC predictions filtered to blood cell types."""
    print("\nLoading blood-matched ABC predictions...")
    
    abc_chunks = []
    for chunk in pd.read_csv(gzip.open(ABC_FILE, 'rt'), sep='\t', chunksize=500000):
        blood_chunk = chunk[chunk['CellType'].isin(BLOOD_CELL_TYPES)]
        if len(blood_chunk) > 0:
            abc_chunks.append(blood_chunk[['chr', 'start', 'end', 'TargetGene', 'ABC.Score', 'CellType']])
    
    abc_df = pd.concat(abc_chunks, ignore_index=True)
    print(f"  Blood-matched ABC predictions: {len(abc_df)}")
    print(f"  Cell types represented: {abc_df['CellType'].nunique()}")
    
    return abc_df


def find_overlapping_abc(enhancer_chr, enhancer_start, enhancer_end, abc_df, window=0):
    """
    Find ABC enhancer predictions overlapping a Gasperini enhancer.
    
    This is the key insight: we only use ABC scores for enhancers that
    directly overlap the tested CRISPR target, not genome-wide.
    """
    # Handle NaN or invalid coordinates
    try:
        enhancer_start = int(enhancer_start)
        enhancer_end = int(enhancer_end)
    except (ValueError, TypeError):
        return pd.DataFrame()
    
    # Normalize chromosome format
    if pd.isna(enhancer_chr):
        return pd.DataFrame()
    chrom = str(enhancer_chr) if str(enhancer_chr).startswith('chr') else f'chr{enhancer_chr}'
    abc_chrom = abc_df[abc_df['chr'] == chrom]
    
    if len(abc_chrom) == 0:
        return pd.DataFrame()
    
    # Find overlapping enhancers (with optional window extension)
    overlapping = abc_chrom[
        (abc_chrom['start'] <= enhancer_end + window) & 
        (abc_chrom['end'] >= enhancer_start - window)
    ]
    
    return overlapping


def compute_noisy_or(probs, epsilon=0.01):
    """
    Compute noisy-OR aggregation of probabilities.
    P = 1 - (1-epsilon) * product((1-p_i) for p_i in probs)
    """
    if len(probs) == 0:
        return epsilon
    
    product = (1 - epsilon)
    for p in probs:
        product *= (1 - p)
    
    return 1 - product


def run_validation(gasperini_df, abc_df, threshold_type='strict'):
    """
    Run variant-specific ABC validation on Gasperini data.
    
    For each enhancer-gene pair:
    1. Find ABC predictions overlapping the enhancer
    2. Get ABC scores for the tested gene
    3. Apply noisy-OR with background probability
    4. Compare validated vs non-validated pairs
    """
    print(f"\n{'='*60}")
    print(f"Running validation ({threshold_type} threshold)")
    print(f"{'='*60}")
    
    results = []
    enhancers_with_abc = 0
    
    for idx, row in gasperini_df.iterrows():
        chrom = row['target_site.chr']
        start = row['target_site.start']
        end = row['target_site.stop']
        target_gene = row['gene_short_name']
        ensg = row['ENSG']
        
        # Is this a validated (significant) pair?
        try:
            validated = float(row['pvalue.empirical.adjusted']) < (0.05 if threshold_type == 'strict' else 0.10)
        except:
            validated = False
        
        # Get effect size for validated pairs
        try:
            fold_change = float(row['fold_change.transcript_remaining'])
            effect_validated = validated and fold_change < 0.9  # >10% knockdown
        except:
            effect_validated = validated
        
        # Find overlapping ABC predictions
        overlapping = find_overlapping_abc(chrom, start, end, abc_df, window=0)
        
        if len(overlapping) == 0:
            continue
        
        enhancers_with_abc += 1
        
        # Get ABC scores for the target gene
        gene_abc = overlapping[overlapping['TargetGene'] == target_gene]
        
        if len(gene_abc) > 0:
            abc_score = gene_abc['ABC.Score'].max()
        else:
            abc_score = 0.0
        
        # Compute path-probability via noisy-OR
        # Using ABC as the main signal (since we don't have L2G for all genes)
        path_prob = compute_noisy_or([abc_score], epsilon=0.01)
        
        results.append({
            'enhancer_chr': chrom,
            'enhancer_start': start,
            'enhancer_end': end,
            'gene': target_gene,
            'ensg': ensg,
            'validated': effect_validated,
            'abc_score': abc_score,
            'path_prob': path_prob,
            'n_overlapping_abc': len(overlapping),
            'n_gene_abc': len(gene_abc)
        })
    
    results_df = pd.DataFrame(results)
    print(f"\nPairs with ABC overlap: {enhancers_with_abc}")
    print(f"Total results: {len(results_df)}")
    
    if len(results_df) == 0:
        return None
    
    # Compute metrics
    validated_df = results_df[results_df['validated']]
    non_validated_df = results_df[~results_df['validated']]
    
    print(f"Validated pairs: {len(validated_df)}")
    print(f"Non-validated pairs: {len(non_validated_df)}")
    
    if len(validated_df) < 5 or len(non_validated_df) < 5:
        print("Insufficient samples for AUROC calculation")
        return results_df
    
    # AUROC for ABC score
    y_true = results_df['validated'].astype(int)
    y_score_abc = results_df['abc_score']
    y_score_pp = results_df['path_prob']
    
    auroc_abc = roc_auc_score(y_true, y_score_abc)
    auroc_pp = roc_auc_score(y_true, y_score_pp)
    
    print(f"\n=== RESULTS ===")
    print(f"ABC Score AUROC: {auroc_abc:.4f}")
    print(f"Path-prob AUROC: {auroc_pp:.4f}")
    
    # Mean scores
    print(f"\nMean ABC (validated): {validated_df['abc_score'].mean():.4f}")
    print(f"Mean ABC (non-validated): {non_validated_df['abc_score'].mean():.4f}")
    
    # Mann-Whitney U test
    stat, pval = stats.mannwhitneyu(
        validated_df['abc_score'],
        non_validated_df['abc_score'],
        alternative='greater'
    )
    print(f"Mann-Whitney U p-value: {pval:.2e}")
    
    # Bootstrap confidence intervals
    print("\nBootstrapping confidence intervals...")
    n_bootstrap = 1000
    aurocs_abc = []
    aurocs_pp = []
    
    np.random.seed(42)
    for _ in range(n_bootstrap):
        indices = np.random.choice(len(results_df), size=len(results_df), replace=True)
        boot_df = results_df.iloc[indices]
        
        if boot_df['validated'].sum() < 2 or (~boot_df['validated']).sum() < 2:
            continue
        
        try:
            aurocs_abc.append(roc_auc_score(boot_df['validated'].astype(int), boot_df['abc_score']))
            aurocs_pp.append(roc_auc_score(boot_df['validated'].astype(int), boot_df['path_prob']))
        except:
            continue
    
    aurocs_abc = np.array(aurocs_abc)
    aurocs_pp = np.array(aurocs_pp)
    
    print(f"\nABC Score AUROC: {auroc_abc:.4f} [95% CI: {np.percentile(aurocs_abc, 2.5):.4f} - {np.percentile(aurocs_abc, 97.5):.4f}]")
    print(f"Path-prob AUROC: {auroc_pp:.4f} [95% CI: {np.percentile(aurocs_pp, 2.5):.4f} - {np.percentile(aurocs_pp, 97.5):.4f}]")
    
    return results_df


def run_validation_against_all_pairs(significant_df, all_enhancers_df, abc_df):
    """
    More rigorous validation: compare validated pairs against ALL tested pairs.
    This is the proper AUROC calculation.
    """
    print(f"\n{'='*60}")
    print("RIGOROUS VALIDATION: Significant vs All Tested Pairs")
    print(f"{'='*60}")
    
    # Create set of validated enhancer-gene pairs
    validated_pairs = set()
    for _, row in significant_df.iterrows():
        key = (row['target_site.chr'], row['target_site.start'], row['target_site.stop'], row['gene_short_name'])
        validated_pairs.add(key)
    
    print(f"Validated pairs: {len(validated_pairs)}")
    
    results = []
    
    for idx, row in all_enhancers_df.iterrows():
        if idx % 10000 == 0:
            print(f"  Processing {idx}/{len(all_enhancers_df)}...")
        
        chrom = row['target_site.chr']
        start = row['target_site.start']
        end = row['target_site.stop']
        target_gene = row['gene_short_name']
        
        # Check if validated
        key = (chrom, start, end, target_gene)
        validated = key in validated_pairs
        
        # Find overlapping ABC predictions
        overlapping = find_overlapping_abc(chrom, start, end, abc_df, window=0)
        
        if len(overlapping) == 0:
            continue
        
        # Get ABC scores for the target gene
        gene_abc = overlapping[overlapping['TargetGene'] == target_gene]
        
        if len(gene_abc) > 0:
            abc_score = gene_abc['ABC.Score'].max()
        else:
            abc_score = 0.0
        
        path_prob = compute_noisy_or([abc_score], epsilon=0.01)
        
        results.append({
            'gene': target_gene,
            'validated': validated,
            'abc_score': abc_score,
            'path_prob': path_prob
        })
    
    results_df = pd.DataFrame(results)
    print(f"\nTotal pairs with ABC: {len(results_df)}")
    print(f"Validated: {results_df['validated'].sum()}")
    print(f"Non-validated: {(~results_df['validated']).sum()}")
    
    # AUROC
    y_true = results_df['validated'].astype(int)
    auroc_abc = roc_auc_score(y_true, results_df['abc_score'])
    auroc_pp = roc_auc_score(y_true, results_df['path_prob'])
    
    print(f"\n=== RIGOROUS AUROC ===")
    print(f"ABC Score AUROC: {auroc_abc:.4f}")
    print(f"Path-prob AUROC: {auroc_pp:.4f}")
    
    # Bootstrap
    print("\nBootstrapping...")
    n_bootstrap = 1000
    aurocs = []
    np.random.seed(42)
    
    for _ in range(n_bootstrap):
        indices = np.random.choice(len(results_df), size=len(results_df), replace=True)
        boot_df = results_df.iloc[indices]
        if boot_df['validated'].sum() >= 2:
            try:
                aurocs.append(roc_auc_score(boot_df['validated'].astype(int), boot_df['abc_score']))
            except:
                pass
    
    aurocs = np.array(aurocs)
    ci_low = np.percentile(aurocs, 2.5)
    ci_high = np.percentile(aurocs, 97.5)
    
    print(f"ABC AUROC: {auroc_abc:.4f} [95% CI: {ci_low:.4f} - {ci_high:.4f}]")
    
    return results_df, auroc_abc, (ci_low, ci_high)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load data
    all_enhancers, significant_05, significant_10 = load_gasperini_data()
    abc_df = load_abc_blood_matched()
    
    # Run quick validation on significant pairs only
    print("\n" + "="*60)
    print("QUICK VALIDATION: Among significant pairs")
    print("="*60)
    results_strict = run_validation(significant_05, abc_df, 'strict')
    
    # Run rigorous validation (all pairs)
    print("\n" + "="*60)
    print("FULL VALIDATION: All enhancer-gene pairs")
    print("="*60)
    
    # Sample for speed (full dataset is 84k pairs)
    sample_size = min(20000, len(all_enhancers))
    print(f"Sampling {sample_size} pairs for efficiency...")
    
    np.random.seed(42)
    sample_idx = np.random.choice(len(all_enhancers), size=sample_size, replace=False)
    sampled_enhancers = all_enhancers.iloc[sample_idx]
    
    full_results, auroc, ci = run_validation_against_all_pairs(
        significant_05, sampled_enhancers, abc_df
    )
    
    # Save results
    output_file = OUTPUT_DIR / "gasperini_variant_specific_abc_results.json"
    results_summary = {
        'dataset': 'Gasperini_2019',
        'cell_type': 'K562',
        'n_total_enhancer_pairs': len(all_enhancers),
        'n_significant_05': len(significant_05),
        'n_significant_10': len(significant_10),
        'n_sampled': sample_size,
        'n_with_abc_overlap': len(full_results),
        'n_validated_with_abc': int(full_results['validated'].sum()),
        'auroc_abc': auroc,
        'auroc_ci_low': ci[0],
        'auroc_ci_high': ci[1],
        'method': 'variant_specific_abc',
        'abc_cell_types': BLOOD_CELL_TYPES[:5]  # Store sample
    }
    
    with open(output_file, 'w') as f:
        json.dump(results_summary, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    
    return full_results, results_summary


if __name__ == "__main__":
    results, summary = main()
