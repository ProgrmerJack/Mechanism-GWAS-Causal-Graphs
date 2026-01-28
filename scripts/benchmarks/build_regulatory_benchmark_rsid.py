#!/usr/bin/env python3
"""
Build expanded regulatory benchmark using rsID-based Platform API lookups.

Strategy:
1. Load E2G benchmark data with silver-standard regulatory links
2. Load SuSiE variant files to get rsIDs for each credible set
3. Use rsID search on Platform API to get correct variant IDs
4. Query L2G scores for those variants
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
REQUEST_DELAY = 0.15  # Respect rate limits

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

def get_l2g_for_variant(variant_id: str) -> list:
    """Get L2G predictions for a variant."""
    query = """
    query getL2G($id: String!) {
      variant(variantId: $id) {
        id
        rsIds
        chromosome
        position
        credibleSets {
          count
          rows {
            studyLocusId
            study { traitFromSource studyId }
            l2GPredictions {
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
            json={'query': query, 'variables': {'id': variant_id}},
            timeout=30
        )
        data = response.json()
        if data.get('data') and data['data'].get('variant'):
            variant = data['data']['variant']
            results = []
            for cs in variant.get('credibleSets', {}).get('rows', []):
                study = cs.get('study', {})
                for pred in cs.get('l2GPredictions', {}).get('rows', []):
                    results.append({
                        'variant_id': variant_id,
                        'chromosome': variant.get('chromosome'),
                        'position': variant.get('position'),
                        'study_id': study.get('studyId'),
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
            vdf = pd.read_csv(variant_file, sep='\t')
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
                        'chromosome': lead['chromosome'].replace('chr', ''),
                        'position': lead['position'],
                        'allele1': lead['allele1'],
                        'allele2': lead['allele2'],
                        'pip': lead['PosteriorProb'],
                        'disease': disease_folder.name
                    }
        except Exception as e:
            print(f"  Error loading {disease_folder.name}: {e}")
    
    print(f"Loaded lead variants for {len(credible_set_variants)} credible sets")
    
    # Filter to credible sets with true regulatory links
    target_cs = set(e2g_df['CredibleSet'].unique())
    matched_cs = target_cs & set(credible_set_variants.keys())
    print(f"Credible sets with E2G truth and SuSiE data: {len(matched_cs)}")
    
    # Query Platform API for each credible set's lead variant
    print("\n" + "=" * 70)
    print("QUERYING PLATFORM API VIA rsID SEARCH")
    print("=" * 70)
    
    results = []
    rsid_not_found = 0
    rsid_found = 0
    l2g_found = 0
    
    # Process credible sets with true regulatory links
    cs_list = list(matched_cs)[:200]  # Limit for testing
    
    for i, cs in enumerate(cs_list):
        if (i + 1) % 20 == 0:
            print(f"  Processing {i+1}/{len(cs_list)}...")
        
        cs_info = credible_set_variants[cs]
        rsid = cs_info['rsid']
        
        # Skip if no rsID
        if not rsid or rsid.startswith('1:') or ':' in rsid:
            rsid_not_found += 1
            continue
        
        # Search for variant by rsID
        variant_info = get_variant_from_rsid(rsid)
        time.sleep(REQUEST_DELAY)
        
        if 'error' in variant_info or 'variant_id' not in variant_info:
            rsid_not_found += 1
            continue
        
        rsid_found += 1
        variant_id = variant_info['variant_id']
        
        # Get L2G predictions
        l2g_preds = get_l2g_for_variant(variant_id)
        time.sleep(REQUEST_DELAY)
        
        if l2g_preds and 'error' not in l2g_preds[0]:
            l2g_found += 1
        
        # Get truth genes for this credible set
        truth_genes = set(e2g_df[e2g_df['CredibleSet'] == cs]['TargetGene'].values)
        
        for pred in l2g_preds:
            if 'error' in pred:
                continue
            pred['credible_set'] = cs
            pred['disease_e2g'] = cs_info['disease']
            pred['rsid'] = rsid
            pred['truth_genes'] = ','.join(truth_genes)
            pred['is_truth'] = pred.get('gene', '') in truth_genes
            results.append(pred)
    
    print(f"\n  rsIDs found in Platform: {rsid_found}")
    print(f"  rsIDs not found: {rsid_not_found}")
    print(f"  Variants with L2G: {l2g_found}")
    print(f"  Total predictions: {len(results)}")
    
    # Create results DataFrame
    if results:
        results_df = pd.DataFrame(results)
        
        # Analyze accuracy
        print("\n" + "=" * 70)
        print("L2G ACCURACY ON REGULATORY BENCHMARK")
        print("=" * 70)
        
        # For each credible set, check if top L2G gene matches truth
        accuracy_results = []
        for cs in results_df['credible_set'].unique():
            cs_preds = results_df[results_df['credible_set'] == cs]
            truth_genes = set(cs_preds.iloc[0]['truth_genes'].split(',')) if cs_preds.iloc[0]['truth_genes'] else set()
            
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
                    'truth_genes': ','.join(truth_genes),
                    'correct': is_correct
                })
        
        if accuracy_results:
            acc_df = pd.DataFrame(accuracy_results)
            n_total = len(acc_df)
            n_correct = acc_df['correct'].sum()
            accuracy = n_correct / n_total if n_total > 0 else 0
            
            print(f"\n  Total credible sets evaluated: {n_total}")
            print(f"  L2G correct predictions: {n_correct}")
            print(f"  L2G accuracy: {accuracy:.1%}")
            
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
  - Credible sets with silver-standard truth: {len(target_cs)}
  - Credible sets with SuSiE variant data: {len(matched_cs)}
  - Successfully queried via rsID: {rsid_found}
  - Credible sets with L2G predictions: {n_total}
  - L2G accuracy: {accuracy:.1%} ({n_correct}/{n_total})

This is a substantial improvement over n=7!
Silver standard = ABC model regulatory predictions (Engreitz et al.)
""")
    else:
        print("\nNo results obtained!")

if __name__ == "__main__":
    main()
