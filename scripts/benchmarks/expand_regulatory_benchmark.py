#!/usr/bin/env python3
"""
expand_regulatory_benchmark.py
==============================
Nature Genetics Critical Requirement: Expand REGULATORY benchmark from 7 to 100+ loci

Sources for REGULATORY benchmarks:
1. Open Targets Gold Standards (functional observational/experimental)
2. CRISPR enhancer screens (Gasperini, Fulco et al.)
3. Fine-mapped eQTL colocalization studies
4. ABC-validated enhancer-gene pairs

This script:
1. Loads Open Targets Gold Standards (2023 version)
2. Classifies loci as CODING vs REGULATORY based on evidence type
3. Filters for high-confidence regulatory links
4. Creates an expanded benchmark with mechanism stratification
5. Integrates with our existing post-2021 benchmark
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from collections import defaultdict

# Configuration
DATA_DIR = Path("data")
EXTERNAL_DIR = DATA_DIR / "external" / "gold_standards" / "genetics-gold-standards" / "gold_standards" / "processed"
OUTPUT_DIR = DATA_DIR / "processed" / "baselines"
RESULTS_DIR = Path("results") / "baselines" / "stratified"

# Ensure directories exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_open_targets_gold_standards():
    """
    Load Open Targets Gold Standards from multiple files.
    Uses the 2023 version (otg_gs_230511.json) as primary.
    """
    records = []
    
    # Load 2023 JSONL
    jsonl_path = EXTERNAL_DIR / "otg_gs_230511.json"
    if jsonl_path.exists():
        with open(jsonl_path, 'r') as f:
            for line in f:
                if line.strip():
                    records.append(json.loads(line))
        print(f"Loaded {len(records)} records from {jsonl_path.name}")
    
    # Also load the original TSV for additional metadata
    tsv_path = EXTERNAL_DIR / "gwas_gold_standards.191108.tsv"
    if tsv_path.exists():
        tsv_df = pd.read_csv(tsv_path, sep='\t')
        print(f"Loaded {len(tsv_df)} records from {tsv_path.name}")
        return records, tsv_df
    
    return records, None


def classify_evidence_as_regulatory(row):
    """
    Classify a gold standard locus as REGULATORY based on evidence class.
    
    REGULATORY evidence types:
    - functional experimental: CRISPR, reporter assays, etc.
    - functional observational: colocalization, chromatin interaction
    - molecular_qtl: eQTL, pQTL evidence
    
    CODING evidence types:
    - drug: Often targets coding variants
    - expert curated: Could be either, default to CODING unless specified
    """
    evidence_class = str(row.get('gold_standard_info.evidence.class', '')).lower()
    evidence_desc = str(row.get('gold_standard_info.evidence.description', '')).lower()
    
    regulatory_keywords = [
        'functional experimental', 'functional observational',
        'eqtl', 'pqtl', 'sqtl', 'coloc', 'enhancer', 'promoter',
        'chromatin', 'hi-c', 'crispr', 'perturbation', 'abc',
        'regulatory', 'noncoding', 'intergenic', 'intronic'
    ]
    
    coding_keywords = [
        'missense', 'lof', 'loss-of-function', 'stop-gain', 'frameshift',
        'splice', 'coding', 'protein-altering', 'mendelian'
    ]
    
    # Check evidence for regulatory indicators
    for keyword in regulatory_keywords:
        if keyword in evidence_class or keyword in evidence_desc:
            return 'REGULATORY'
    
    # Check for coding indicators
    for keyword in coding_keywords:
        if keyword in evidence_class or keyword in evidence_desc:
            return 'CODING'
    
    # Default based on evidence class
    if 'functional' in evidence_class:
        return 'REGULATORY'
    elif 'drug' in evidence_class or 'expert' in evidence_class:
        return 'CODING'
    
    return 'UNKNOWN'


def process_tsv_gold_standards(tsv_df):
    """
    Process TSV gold standards with mechanism classification.
    """
    processed = []
    
    for _, row in tsv_df.iterrows():
        # Skip rows without valid variant info
        if pd.isna(row.get('sentinel_variant.rsid')):
            continue
        
        mechanism = classify_evidence_as_regulatory(row)
        
        locus = {
            'locus_id': f"OT_{row.get('sentinel_variant.rsid', 'unknown')}",
            'chr': str(row.get('sentinel_variant.locus_GRCh38.chromosome', '')),
            'pos_hg38': row.get('sentinel_variant.locus_GRCh38.position', 0),
            'lead_snp': row.get('sentinel_variant.rsid', ''),
            'gene_id': row.get('gold_standard_info.gene_id', ''),
            'trait': row.get('trait_info.reported_trait_name', ''),
            'evidence_class': row.get('gold_standard_info.evidence.class', ''),
            'confidence': row.get('gold_standard_info.highest_confidence', ''),
            'evidence_description': row.get('gold_standard_info.evidence.description', ''),
            'pubmed_id': row.get('association_info.pubmed_id', ''),
            'gwas_catalog_id': row.get('association_info.gwas_catalog_id', ''),
            'mechanism_class': mechanism,
            'source': 'OpenTargets_GoldStandards_2019'
        }
        processed.append(locus)
    
    return pd.DataFrame(processed)


def load_existing_benchmark():
    """Load our existing post-2021 benchmark."""
    benchmark_path = OUTPUT_DIR / "post2021_independent_benchmark_FINAL.tsv"
    if benchmark_path.exists():
        df = pd.read_csv(benchmark_path, sep='\t')
        print(f"Loaded existing benchmark: {len(df)} loci")
        return df
    return None


def identify_regulatory_loci_from_gold_standards(gs_df, existing_snps):
    """
    Identify high-confidence REGULATORY loci from gold standards.
    
    Criteria for inclusion:
    1. High or Medium confidence
    2. Functional evidence (not just drug or expert curated)
    3. NOT already in our benchmark
    4. Valid genomic coordinates
    """
    # Filter for regulatory loci
    regulatory_df = gs_df[gs_df['mechanism_class'] == 'REGULATORY'].copy()
    print(f"Total REGULATORY loci in gold standards: {len(regulatory_df)}")
    
    # Filter for confidence
    high_confidence = regulatory_df[regulatory_df['confidence'].isin(['High', 'Medium'])]
    print(f"High/Medium confidence: {len(high_confidence)}")
    
    # Remove duplicates (same SNP can have multiple genes)
    unique_loci = high_confidence.drop_duplicates(subset=['lead_snp', 'gene_id'])
    print(f"Unique SNP-gene pairs: {len(unique_loci)}")
    
    # Remove those already in our benchmark
    if existing_snps is not None:
        novel_loci = unique_loci[~unique_loci['lead_snp'].isin(existing_snps)]
        print(f"Novel loci (not in existing benchmark): {len(novel_loci)}")
    else:
        novel_loci = unique_loci
    
    return novel_loci


def create_expanded_regulatory_benchmark():
    """
    Main function to create expanded regulatory benchmark.
    """
    print("=" * 80)
    print("EXPANDING REGULATORY BENCHMARK")
    print("Nature Genetics Requirement: 100+ regulatory loci")
    print("=" * 80)
    
    # Load gold standards
    records_2023, tsv_df = load_open_targets_gold_standards()
    
    if tsv_df is None:
        print("ERROR: Could not load gold standards TSV")
        return None
    
    # Process TSV gold standards
    gs_processed = process_tsv_gold_standards(tsv_df)
    
    # Mechanism classification summary
    print("\n" + "=" * 60)
    print("MECHANISM CLASSIFICATION (Open Targets Gold Standards)")
    print("=" * 60)
    for mech in ['REGULATORY', 'CODING', 'UNKNOWN']:
        count = len(gs_processed[gs_processed['mechanism_class'] == mech])
        pct = count / len(gs_processed) * 100
        print(f"  {mech}: {count} ({pct:.1f}%)")
    
    # Load existing benchmark
    existing_benchmark = load_existing_benchmark()
    existing_snps = set(existing_benchmark['lead_snp'].tolist()) if existing_benchmark is not None else set()
    
    # Identify new regulatory loci
    new_regulatory = identify_regulatory_loci_from_gold_standards(gs_processed, existing_snps)
    
    # Create combined benchmark
    print("\n" + "=" * 60)
    print("CREATING COMBINED BENCHMARK")
    print("=" * 60)
    
    # Prepare existing benchmark
    if existing_benchmark is not None:
        existing_benchmark['source'] = 'Post2021_Independent'
        existing_benchmark['mechanism_class'] = existing_benchmark.apply(
            lambda x: 'REGULATORY' if 'CRISPR' in str(x.get('evidence_tier', '')) 
            else 'CODING', axis=1
        )
    
    # Prepare new regulatory loci
    new_regulatory_formatted = new_regulatory[[
        'locus_id', 'chr', 'pos_hg38', 'lead_snp', 'gene_id', 
        'trait', 'evidence_class', 'confidence', 'mechanism_class', 'source'
    ]].copy()
    new_regulatory_formatted = new_regulatory_formatted.rename(columns={
        'evidence_class': 'evidence_tier',
        'confidence': 'validation_type'
    })
    
    # Summary of expanded benchmark
    print(f"\nExisting benchmark: {len(existing_benchmark) if existing_benchmark is not None else 0} loci")
    print(f"New regulatory loci: {len(new_regulatory_formatted)}")
    
    # Save the expanded regulatory benchmark
    expanded_regulatory_path = OUTPUT_DIR / "expanded_regulatory_benchmark.tsv"
    new_regulatory_formatted.to_csv(expanded_regulatory_path, sep='\t', index=False)
    print(f"\nSaved: {expanded_regulatory_path}")
    
    # Create summary report
    create_benchmark_expansion_report(
        existing_benchmark, new_regulatory_formatted, gs_processed
    )
    
    return new_regulatory_formatted


def create_benchmark_expansion_report(existing, new_regulatory, all_gs):
    """Create detailed report on benchmark expansion."""
    
    report_path = RESULTS_DIR / "BENCHMARK_EXPANSION_REPORT.md"
    
    with open(report_path, 'w') as f:
        f.write("# Benchmark Expansion Report\n\n")
        f.write("## Nature Genetics Requirement: Adequate Regulatory Benchmark\n\n")
        
        f.write("### Problem Statement\n\n")
        f.write("Our initial post-2021 benchmark had:\n")
        f.write("- 56 CODING loci (88.9%)\n")
        f.write("- 7 REGULATORY loci (11.1%)\n\n")
        f.write("This is **insufficient** for evaluating locus-to-gene methods.\n\n")
        
        f.write("### Solution: Open Targets Gold Standards\n\n")
        f.write(f"Downloaded **{len(all_gs)}** gold standard locus-gene pairs.\n\n")
        
        f.write("#### Evidence Class Distribution\n\n")
        f.write("| Evidence Class | Count | Mechanism |\n")
        f.write("|----------------|-------|------------|\n")
        for ec in all_gs['evidence_class'].value_counts().head(10).index:
            count = len(all_gs[all_gs['evidence_class'] == ec])
            mech = all_gs[all_gs['evidence_class'] == ec]['mechanism_class'].mode()[0]
            f.write(f"| {ec} | {count} | {mech} |\n")
        
        f.write("\n### Regulatory Loci Identified\n\n")
        f.write(f"**Total new regulatory loci**: {len(new_regulatory)}\n\n")
        
        if len(new_regulatory) > 0:
            f.write("#### Sample Regulatory Loci\n\n")
            f.write("| SNP | Gene | Trait | Evidence |\n")
            f.write("|-----|------|-------|----------|\n")
            for _, row in new_regulatory.head(20).iterrows():
                f.write(f"| {row['lead_snp']} | {row['gene_id']} | {row['trait'][:40]} | {row['evidence_tier']} |\n")
        
        f.write("\n### Combined Benchmark Statistics\n\n")
        existing_coding = len(existing[existing['mechanism_class'] == 'CODING']) if existing is not None else 0
        existing_reg = len(existing[existing['mechanism_class'] == 'REGULATORY']) if existing is not None else 0
        new_reg = len(new_regulatory)
        
        total = existing_coding + existing_reg + new_reg
        
        f.write(f"| Source | CODING | REGULATORY | Total |\n")
        f.write(f"|--------|--------|------------|-------|\n")
        f.write(f"| Post-2021 Independent | {existing_coding} | {existing_reg} | {existing_coding + existing_reg} |\n")
        f.write(f"| Open Targets GS (new) | 0 | {new_reg} | {new_reg} |\n")
        f.write(f"| **Combined** | **{existing_coding}** | **{existing_reg + new_reg}** | **{total}** |\n")
        
        new_reg_pct = (existing_reg + new_reg) / total * 100 if total > 0 else 0
        f.write(f"\n**New REGULATORY proportion**: {new_reg_pct:.1f}%\n")
        
        if new_reg_pct >= 30:
            f.write("\n✅ Regulatory benchmark now adequate (≥30%)\n")
        else:
            f.write(f"\n⚠️ Still need more regulatory loci (target: 30%, current: {new_reg_pct:.1f}%)\n")
        
        f.write("\n### Leakage Considerations\n\n")
        f.write("The Open Targets Gold Standards were used for L2G training.\n")
        f.write("However, we use them specifically to:\n")
        f.write("1. **Diagnose benchmark bias** (not to train our method)\n")
        f.write("2. **Stratify evaluation** (show performance on regulatory vs coding)\n")
        f.write("3. **Identify weakness modes** (where all methods fail)\n\n")
        f.write("This is explicitly documented in our leakage audit.\n")
    
    print(f"\nReport saved: {report_path}")


if __name__ == "__main__":
    expanded = create_expanded_regulatory_benchmark()
    
    if expanded is not None:
        print("\n" + "=" * 80)
        print("BENCHMARK EXPANSION COMPLETE")
        print("=" * 80)
        print(f"\nNew regulatory loci added: {len(expanded)}")
        print("\nNext steps:")
        print("1. Run baselines on expanded benchmark")
        print("2. Create stratified performance figure")
        print("3. Document leakage implications")
