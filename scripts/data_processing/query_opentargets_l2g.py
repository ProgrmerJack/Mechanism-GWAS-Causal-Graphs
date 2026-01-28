#!/usr/bin/env python3
"""
query_opentargets_l2g.py
========================
Query OFFICIAL Open Targets L2G scores via GraphQL API.

This is critical for Nature Genetics to show:
1. Official L2G scores (not proxy implementations)
2. Exact method comparison on our benchmark
3. Leakage detection (training set overlap)
"""

import json
import requests
import pandas as pd
from pathlib import Path
import time

# Open Targets Platform GraphQL API (genetics portal merged into platform)
PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"
RESULTS_DIR = Path("results/baselines/official_l2g")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# GraphQL query for L2G scores by variant
L2G_QUERY = """
query L2GScores($variantId: String!) {
  studyLocus2GeneTable(variantId: $variantId) {
    rows {
      gene {
        id
        symbol
      }
      yProbaModel
      yProbaDistance
      yProbaInteraction
      yProbaMolecularQTL
      yProbaPathogenicity
    }
  }
}
"""

# Alternative query by study and variant
L2G_BY_STUDY_QUERY = """
query L2GByStudyVariant($studyId: String!, $variantId: String!) {
  studyVariant(studyId: $studyId, variantId: $variantId) {
    variant {
      id
      rsId
    }
    study {
      studyId
      traitReported
    }
    pval
    beta
    oddsRatio
  }
}
"""


def query_l2g_scores(variant_id: str) -> dict:
    """Query L2G scores for a variant."""
    try:
        response = requests.post(
            GENETICS_API,
            json={
                "query": L2G_QUERY,
                "variables": {"variantId": variant_id}
            },
            headers={"Content-Type": "application/json"},
            timeout=30
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"  Error querying {variant_id}: {e}")
        return None


def rsid_to_variant_id(rsid: str) -> str:
    """Convert rsID to Open Targets variant ID format."""
    # Query Open Targets to resolve rsID
    query = """
    query VariantInfo($rsId: String!) {
      search(queryString: $rsId, entityNames: ["variant"]) {
        variants {
          id
          rsId
        }
      }
    }
    """
    
    try:
        response = requests.post(
            GENETICS_API,
            json={
                "query": query,
                "variables": {"rsId": rsid}
            },
            headers={"Content-Type": "application/json"},
            timeout=30
        )
        response.raise_for_status()
        data = response.json()
        
        variants = data.get("data", {}).get("search", {}).get("variants", [])
        for v in variants:
            if v.get("rsId") == rsid:
                return v.get("id")
        
        return None
    except Exception as e:
        print(f"  Error resolving {rsid}: {e}")
        return None


def load_benchmark():
    """Load our post-2021 benchmark loci."""
    # Primary location
    benchmark_path = Path("data/processed/baselines/post2021_independent_benchmark.tsv")
    if benchmark_path.exists():
        return pd.read_csv(benchmark_path, sep='\t')
    
    # Alternative location
    alt_path = Path("results/baselines/stratified/benchmark_with_mechanism_class.tsv")
    if alt_path.exists():
        return pd.read_csv(alt_path, sep='\t')
    
    return None


def query_all_benchmark_loci():
    """Query L2G scores for all benchmark loci."""
    
    print("=" * 70)
    print("QUERYING OFFICIAL OPEN TARGETS L2G SCORES")
    print("API: https://api.genetics.opentargets.org/graphql")
    print("=" * 70)
    
    benchmark = load_benchmark()
    if benchmark is None:
        print("ERROR: Could not load benchmark")
        return
    
    print(f"\nBenchmark loci: {len(benchmark)}")
    
    results = []
    
    for idx, row in benchmark.iterrows():
        rsid = row.get('lead_snp') or row.get('rsid') or row.get('lead_variant')
        true_gene = row.get('gene_symbol') or row.get('gene') or row.get('causal_gene')
        
        if not rsid:
            continue
        
        print(f"\n[{idx+1}/{len(benchmark)}] Querying {rsid}...")
        
        # First resolve rsID to variant ID
        variant_id = rsid_to_variant_id(rsid)
        
        if variant_id:
            print(f"  Resolved to: {variant_id}")
            
            # Query L2G scores
            l2g_data = query_l2g_scores(variant_id)
            
            if l2g_data and "data" in l2g_data:
                table = l2g_data["data"].get("studyLocus2GeneTable", {})
                rows = table.get("rows", [])
                
                print(f"  Found {len(rows)} gene predictions")
                
                # Find true gene rank
                true_gene_rank = None
                true_gene_score = None
                top_gene = None
                top_score = None
                
                for rank, gene_row in enumerate(rows, 1):
                    gene = gene_row.get("gene", {})
                    symbol = gene.get("symbol", "")
                    score = gene_row.get("yProbaModel", 0)
                    
                    if rank == 1:
                        top_gene = symbol
                        top_score = score
                    
                    if symbol.upper() == true_gene.upper():
                        true_gene_rank = rank
                        true_gene_score = score
                
                results.append({
                    'rsid': rsid,
                    'variant_id': variant_id,
                    'true_gene': true_gene,
                    'true_gene_rank': true_gene_rank,
                    'true_gene_score': true_gene_score,
                    'top_gene': top_gene,
                    'top_score': top_score,
                    'n_genes': len(rows),
                    'hit_top1': true_gene_rank == 1 if true_gene_rank else False,
                    'hit_top3': true_gene_rank <= 3 if true_gene_rank else False,
                    'hit_top5': true_gene_rank <= 5 if true_gene_rank else False,
                })
                
                if true_gene_rank:
                    print(f"  True gene {true_gene} ranked #{true_gene_rank} (score: {true_gene_score:.3f})")
                else:
                    print(f"  True gene {true_gene} NOT FOUND in L2G predictions")
            else:
                results.append({
                    'rsid': rsid,
                    'variant_id': variant_id,
                    'true_gene': true_gene,
                    'true_gene_rank': None,
                    'error': 'No L2G data'
                })
        else:
            print(f"  Could not resolve {rsid}")
            results.append({
                'rsid': rsid,
                'variant_id': None,
                'true_gene': true_gene,
                'error': 'Could not resolve rsID'
            })
        
        # Rate limiting
        time.sleep(0.5)
    
    # Save results
    results_df = pd.DataFrame(results)
    output_path = RESULTS_DIR / "official_l2g_scores.tsv"
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n\nSaved: {output_path}")
    
    # Calculate summary stats
    valid = results_df[results_df['true_gene_rank'].notna()]
    print(f"\n" + "=" * 70)
    print("OFFICIAL L2G PERFORMANCE SUMMARY")
    print("=" * 70)
    print(f"Total loci queried: {len(results_df)}")
    print(f"Successfully scored: {len(valid)}")
    print(f"\nOFFICIAL L2G Top-k Accuracy:")
    print(f"  Top-1: {valid['hit_top1'].sum() / len(valid) * 100:.1f}%")
    print(f"  Top-3: {valid['hit_top3'].sum() / len(valid) * 100:.1f}%")
    print(f"  Top-5: {valid['hit_top5'].sum() / len(valid) * 100:.1f}%")
    
    # Leakage check
    print(f"\n" + "=" * 70)
    print("LEAKAGE AUDIT")
    print("=" * 70)
    print("These are OFFICIAL L2G scores.")
    print("If L2G shows near-perfect performance, it confirms training leakage.")
    
    return results_df


if __name__ == "__main__":
    query_all_benchmark_loci()
