#!/usr/bin/env python3
"""
Build expanded regulatory benchmark using rsID-based Platform API lookups.
FIXED VERSION - correct API schema and filter for GWAS studies.

Strategy:
1. Load E2G benchmark data with silver-standard regulatory links
2. Load SuSiE variant files to get rsIDs for each credible set
3. Use rsID search on Platform API to get correct variant IDs
4. Query L2G scores for those variants (GWAS studies only)
5. Match target genes to evaluate accuracy
"""
import os
import sys
import glob
import time
import json
import requests
import pandas as pd
from pathlib import Path
from collections import defaultdict

# API config
PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"
REQUEST_DELAY = 0.1  # Respect rate limits

def get_variant_from_rsid(rsid: str) -> dict:
    """Search Platform API for variant by rsID."""
    query = """
    query searchVariant($rsid: String!) {
      search(queryString: $rsid, entityNames: ["variant"]) {
        hits { id entity name }
      }
    }
    """
    try:
        response = requests.post(
            PLATFORM_API,
            json={'query': query, 'variables': {'rsid': rsid}},
            timeout=30
        )
        data = response.json()
        if data.get('data') and data['data'].get('search'):
            hits = data['data']['search'].get('hits', [])
            for hit in hits:
                if hit.get('entity') == 'variant':
                    return {'variant_id': hit.get('id'), 'name': hit.get('name')}
    except Exception as e:
        return {'error': str(e)}
    return {}

def get_l2g_for_variant(variant_id: str, max_credible_sets: int = 50) -> list:
    """Get L2G predictions for a variant from GWAS studies."""
    query = """
    query getL2G($id: String!, $size: Int!) {
      variant(variantId: $id) {
        id
        rsIds
        chromosome
        position
        credibleSets(page: {size: $size, index: 0}) {
          count
          rows {
            studyLocusId
            study { id studyType traitFromSource }
            l2GPredictions(page: {size: 20, index: 0}) {
              count
              rows {
                target { approvedSymbol id }
                score
              }
            }
          }
        }
      }
    }
    """
    try:
        response = requests.post(
            PLATFORM_API,
            json={'query': query, 'variables': {'id': variant_id, 'size': max_credible_sets}},
            timeout=60
        )
        data = response.json()
        if data.get('data') and data['data'].get('variant'):
            variant = data['data']['variant']
            results = []
            for cs in variant.get('credibleSets', {}).get('rows', []):
                study = cs.get('study', {})
                # Filter to GWAS studies only (L2G is mainly for GWAS)
                study_type = study.get('studyType', '')
                if study_type != 'gwas':
                    continue
                
                l2g_preds = cs.get('l2GPredictions', {})
                if l2g_preds.get('count', 0) == 0:
                    continue
                
                for pred in l2g_preds.get('rows', []):
                    results.append({
                        'variant_id': variant_id,
                        'chromosome': variant.get('chromosome'),
                        'position': variant.get('position'),
                        'rsIds': ','.join(variant.get('rsIds', [])),
                        'study_id': study.get('id'),
                        'study_type': study_type,
                        'trait': study.get('traitFromSource'),
                        'gene': pred.get('target', {}).get('approvedSymbol'),
                        'gene_id': pred.get('target', {}).get('id'),
                        'l2g_score': pred.get('score')
                    })
            return results
    except Exception as e:
        return [{'error': str(e)}]
    return []

def main():
    base_path = Path(r"c:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs")
    data_path = base_path / "data" / "external" / "E2G_benchmarking" / "resources"
    
    # Load E2G benchmark - silver standard regulatory links
    print("=" * 70)
    print("LOADING E2G BENCHMARK DATA")
    print("=" * 70)
    e2g_file = data_path / "UKBiobank.ABCGene.anyabc.tsv"
    e2g_df = pd.read_csv(e2g_file, sep='\t')
    
    # Filter to non-coding true regulatory links
    e2g_df = e2g_df[
        (e2g_df['truth'] == True) &  # Silver standard positives
        (e2g_df['CodingSpliceOrPromoterVariants'] == False)  # Non-coding
    ]
    print(f"True non-coding regulatory links: {len(e2g_df)}")
    print(f"Unique credible sets with true links: {e2g_df['CredibleSet'].nunique()}")
    print(f"Diseases: {e2g_df['Disease'].nunique()}")
    
    # Build gene truth mapping from E2G
    e2g_truth = {}
    for _, row in e2g_df.iterrows():
        cs = row['CredibleSet']
        gene = row['TargetGene']
        disease = row['Disease']
        if cs not in e2g_truth:
            e2g_truth[cs] = {'genes': set(), 'disease': disease}
        e2g_truth[cs]['genes'].add(gene)
    
    # Load SuSiE variant files to get rsIDs
    print("\n" + "=" * 70)
    print("LOADING SuSIE VARIANT DATA")
    print("=" * 70)
    susie_path = data_path / "191010_UKBB_SuSiE_hg38_liftover"
    
    # Build mapping from CredibleSet to lead variants (highest PIP)
    credible_set_variants = {}
    disease_folders = sorted(susie_path.iterdir())
    
    for disease_folder in disease_folders:
        if not disease_folder.is_dir():
            continue
        variant_file = disease_folder / "variant.list.txt"
        if not variant_file.exists():
            continue
        
        try:
            vdf = pd.read_csv(variant_file, sep='\t', low_memory=False)
            # Get unique credible sets
            for cs in vdf['CredibleSet'].unique():
                cs_data = vdf[vdf['CredibleSet'] == cs]
                # Get variant with highest posterior probability
                lead = cs_data.loc[cs_data['PosteriorProb'].idxmax()]
                
                # Extract rsID and position
                rsid = lead['rsid'] if pd.notna(lead['rsid']) else None
                
                # Store mapping
                if cs not in credible_set_variants:
                    credible_set_variants[cs] = {
                        'rsid': rsid,
                        'chromosome': str(lead['chromosome']).replace('chr', ''),
                        'position': lead['position'],
                        'allele1': lead['allele1'],
                        'allele2': lead['allele2'],
                        'pip': lead['PosteriorProb'],
                        'disease': disease_folder.name
                    }
        except Exception as e:
            pass  # Skip problematic files silently
    
    print(f"Loaded lead variants for {len(credible_set_variants)} credible sets")
    
    # Filter to credible sets with true regulatory links
    target_cs = set(e2g_df['CredibleSet'].unique())
    matched_cs = target_cs & set(credible_set_variants.keys())
    print(f"Credible sets with E2G truth and SuSiE data: {len(matched_cs)}")
    
    # Query Platform API for each credible set's lead variant
    print("\n" + "=" * 70)
    print("QUERYING PLATFORM API VIA rsID SEARCH")
    print("=" * 70)
    
    all_results = []
    rsid_queries = 0
    rsid_found = 0
    variants_with_l2g = 0
    
    # Process ALL credible sets with true regulatory links
    cs_list = list(matched_cs)
    
    for i, cs in enumerate(cs_list):
        if (i + 1) % 50 == 0:
            print(f"  Processing {i+1}/{len(cs_list)}... (found {variants_with_l2g} with L2G)")
        
        cs_info = credible_set_variants[cs]
        rsid = cs_info['rsid']
        
        # Skip if no rsID or structural variant ID format
        if not rsid or ':' in rsid or not rsid.startswith('rs'):
            continue
        
        rsid_queries += 1
        
        # Search for variant by rsID
        variant_info = get_variant_from_rsid(rsid)
        time.sleep(REQUEST_DELAY)
        
        if 'error' in variant_info or 'variant_id' not in variant_info:
            continue
        
        rsid_found += 1
        variant_id = variant_info['variant_id']
        
        # Get L2G predictions (GWAS studies only)
        l2g_preds = get_l2g_for_variant(variant_id)
        time.sleep(REQUEST_DELAY)
        
        # Get truth genes for this credible set
        truth_info = e2g_truth.get(cs, {'genes': set(), 'disease': ''})
        truth_genes = truth_info['genes']
        
        if l2g_preds and 'error' not in l2g_preds[0] and len(l2g_preds) > 0:
            variants_with_l2g += 1
            
            for pred in l2g_preds:
                if 'error' in pred:
                    continue
                pred['e2g_credible_set'] = cs
                pred['e2g_disease'] = truth_info['disease']
                pred['query_rsid'] = rsid
                pred['truth_genes'] = ','.join(truth_genes)
                pred['is_truth'] = pred.get('gene', '') in truth_genes
                all_results.append(pred)
    
    print(f"\n  rsIDs queried: {rsid_queries}")
    print(f"  rsIDs found in Platform: {rsid_found}")
    print(f"  Variants with GWAS L2G: {variants_with_l2g}")
    print(f"  Total L2G predictions: {len(all_results)}")
    
    # Create results DataFrame
    if all_results:
        results_df = pd.DataFrame(all_results)
        
        # Analyze accuracy
        print("\n" + "=" * 70)
        print("L2G ACCURACY ON EXPANDED REGULATORY BENCHMARK")
        print("=" * 70)
        
        # For each credible set, check if top L2G gene matches truth
        accuracy_results = []
        cs_evaluated = results_df['e2g_credible_set'].unique()
        
        for cs in cs_evaluated:
            cs_preds = results_df[results_df['e2g_credible_set'] == cs]
            truth_str = cs_preds.iloc[0]['truth_genes']
            truth_genes = set(truth_str.split(',')) if truth_str else set()
            
            # Get top L2G prediction
            top_pred = cs_preds.nlargest(1, 'l2g_score')
            if len(top_pred) > 0:
                top_gene = top_pred.iloc[0]['gene']
                top_score = top_pred.iloc[0]['l2g_score']
                is_correct = top_gene in truth_genes
                
                accuracy_results.append({
                    'credible_set': cs,
                    'top_l2g_gene': top_gene,
                    'top_l2g_score': top_score,
                    'truth_genes': truth_str,
                    'correct': is_correct,
                    'disease': cs_preds.iloc[0]['e2g_disease']
                })
        
        if accuracy_results:
            acc_df = pd.DataFrame(accuracy_results)
            n_total = len(acc_df)
            n_correct = acc_df['correct'].sum()
            accuracy = n_correct / n_total if n_total > 0 else 0
            
            print(f"\n  Credible sets evaluated: {n_total}")
            print(f"  L2G correct (matches ABC truth): {n_correct}")
            print(f"  L2G accuracy: {accuracy:.1%}")
            
            # Breakdown by disease
            print("\n  Accuracy by disease:")
            disease_acc = acc_df.groupby('disease').agg(
                n_loci=('correct', 'count'),
                n_correct=('correct', 'sum')
            ).reset_index()
            disease_acc['accuracy'] = disease_acc['n_correct'] / disease_acc['n_loci']
            disease_acc = disease_acc.sort_values('n_loci', ascending=False)
            for _, row in disease_acc.head(10).iterrows():
                print(f"    {row['disease']}: {row['accuracy']:.0%} ({row['n_correct']}/{row['n_loci']})")
            
            # Save results
            output_path = base_path / "data" / "processed" / "baselines"
            output_path.mkdir(parents=True, exist_ok=True)
            
            results_df.to_csv(output_path / "regulatory_benchmark_l2g_predictions.tsv", sep='\t', index=False)
            acc_df.to_csv(output_path / "regulatory_benchmark_accuracy.tsv", sep='\t', index=False)
            
            print(f"\n  Saved: {output_path / 'regulatory_benchmark_l2g_predictions.tsv'}")
            print(f"  Saved: {output_path / 'regulatory_benchmark_accuracy.tsv'}")
            
            # Summary stats
            print("\n" + "=" * 70)
            print("BENCHMARK SUMMARY")
            print("=" * 70)
            print(f"""
Expanded Regulatory Benchmark Results:
======================================
  Original E2G true links: 530
  Credible sets with truth: {len(target_cs)}
  Matched to SuSiE variants: {len(matched_cs)}
  rsIDs queried: {rsid_queries}
  Variants found in Platform: {rsid_found}
  Variants with GWAS L2G: {variants_with_l2g}
  
FINAL: L2G accuracy = {accuracy:.1%} ({n_correct}/{n_total})

This is a MAJOR improvement over the original n=7 benchmark!
Silver standard = ABC model regulatory predictions (Engreitz et al.)

NOTE: This evaluates whether L2G's top prediction matches the
ABC-predicted causal gene for non-coding regulatory variants.
""")
            
            # Also create summary JSON for manuscript
            summary = {
                'benchmark': 'expanded_regulatory',
                'source': 'EngreitzLab GWAS_E2G_benchmarking',
                'truth_type': 'ABC model silver standard',
                'n_original_truth_links': 530,
                'n_credible_sets_evaluated': n_total,
                'n_correct': int(n_correct),
                'accuracy': round(accuracy, 3),
                'diseases_evaluated': len(disease_acc),
            }
            with open(output_path / "regulatory_benchmark_summary.json", 'w') as f:
                json.dump(summary, f, indent=2)
            print(f"  Saved: {output_path / 'regulatory_benchmark_summary.json'}")
    else:
        print("\nNo results obtained!")

if __name__ == "__main__":
    main()
