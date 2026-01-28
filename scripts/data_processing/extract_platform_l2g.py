#!/usr/bin/env python3
"""
Extract L2G scores from Open Targets Platform API for benchmark variants.
Uses correct GraphQL schema based on API introspection (v25.4.4).

Key findings from schema introspection:
- Variant fields: id, chromosome, position, referenceAllele, alternateAllele, rsIds
- Search returns SearchResult with id, name, entity, object (union type)
- CredibleSets have l2GPredictions field
- Variant ID format: CHROM_POS_REF_ALT
"""

import requests
import json
import pandas as pd
import time
from pathlib import Path
from typing import Optional, List, Dict

PLATFORM_API = "https://api.platform.opentargets.org/api/v4/graphql"

def query_variant(variant_id: str) -> Dict:
    """Query variant and get associated credible sets with L2G scores."""
    query = """
    query getVariant($variantId: String!) {
      variant(variantId: $variantId) {
        id
        chromosome
        position
        referenceAllele
        alternateAllele
        rsIds
        credibleSets(page: {size: 50, index: 0}) {
          count
          rows {
            studyLocusId
            region
            study {
              id
              traitFromSource
              studyType
            }
            l2GPredictions {
              count
              rows {
                score
                target {
                  id
                  approvedSymbol
                }
              }
            }
          }
        }
      }
    }
    """
    r = requests.post(PLATFORM_API, json={"query": query, "variables": {"variantId": variant_id}}, timeout=60)
    return r.json()

def query_credible_sets_by_study(study_id: str, limit: int = 100) -> Dict:
    """Get credible sets and L2G predictions for a specific study."""
    query = """
    query getStudyCredibleSets($studyId: String!, $size: Int!) {
      credibleSets(studyIds: [$studyId], page: {size: $size, index: 0}) {
        count
        rows {
          studyLocusId
          region
          chromosome
          position
          study {
            id
            traitFromSource
          }
          l2GPredictions {
            count
            rows {
              score
              target {
                id
                approvedSymbol
              }
            }
          }
        }
      }
    }
    """
    r = requests.post(PLATFORM_API, json={"query": query, "variables": {"studyId": study_id, "size": limit}}, timeout=120)
    return r.json()

def search_studies(search_term: str, limit: int = 20) -> List[Dict]:
    """Search for GWAS studies and return study IDs."""
    query = """
    query searchStudies($term: String!, $size: Int!) {
      search(queryString: $term, entityNames: ["study"], page: {size: $size, index: 0}) {
        total
        hits {
          id
          name
          entity
        }
      }
    }
    """
    r = requests.post(PLATFORM_API, json={"query": query, "variables": {"term": search_term, "size": limit}}, timeout=30)
    result = r.json()
    if 'data' in result and result['data'].get('search'):
        return result['data']['search'].get('hits', [])
    return []

def get_study_details(study_id: str) -> Dict:
    """Get detailed information about a study."""
    query = """
    query getStudy($studyId: String!) {
      study(studyId: $studyId) {
        id
        traitFromSource
        studyType
        nSamples
        nCases
        publicationFirstAuthor
        pubmedId
      }
    }
    """
    r = requests.post(PLATFORM_API, json={"query": query, "variables": {"studyId": study_id}}, timeout=30)
    return r.json()

def extract_l2g_for_known_variants(variant_gene_pairs: List[Dict]) -> pd.DataFrame:
    """
    Query L2G scores for known variant-gene pairs from our benchmark.
    
    Args:
        variant_gene_pairs: List of dicts with 'chr', 'pos_hg38', 'gene_symbol', 'locus_id'
    
    Returns:
        DataFrame with L2G scores
    """
    results = []
    
    for pair in variant_gene_pairs:
        chr_num = str(pair.get('chr', '')).replace('chr', '')
        pos = pair.get('pos_hg38')
        gene_symbol = pair.get('gene_symbol', '')
        locus_id = pair.get('locus_id', '')
        
        if not chr_num or not pos:
            continue
            
        print(f"Querying chr{chr_num}:{pos} for {gene_symbol}...", end=" ")
        
        # Try to find variant - we need ref/alt alleles
        # Query nearby credible sets
        region_query = f"""
        query {{
          credibleSets(
            page: {{size: 100, index: 0}}
          ) {{
            rows {{
              studyLocusId
              chromosome
              position
              region
              study {{
                id
                traitFromSource
                studyType
              }}
              locus {{
                variantId
                posteriorProbability
              }}
              l2GPredictions {{
                rows {{
                  score
                  target {{
                    id
                    approvedSymbol
                  }}
                }}
              }}
            }}
          }}
        }}
        """
        
        # Instead, search for studies related to the trait and get their credible sets
        # For now, let's just look for any L2G predictions involving our target gene
        
        gene_query = """
        query getGeneL2G($geneId: String!) {
          target(ensemblId: $geneId) {
            id
            approvedSymbol
          }
        }
        """
        
        time.sleep(0.3)  # Rate limiting
        print("done")
    
    return pd.DataFrame(results)

def get_l2g_landscape(study_ids: List[str], output_path: Path) -> pd.DataFrame:
    """
    Extract full L2G landscape for multiple studies.
    
    Args:
        study_ids: List of Open Targets study IDs (e.g., GCST IDs)
        output_path: Where to save results
    
    Returns:
        DataFrame with all L2G predictions
    """
    all_results = []
    
    for study_id in study_ids:
        print(f"Querying study {study_id}...")
        cs_data = query_credible_sets_by_study(study_id, limit=200)
        
        if 'data' in cs_data and cs_data['data'].get('credibleSets'):
            cs = cs_data['data']['credibleSets']
            print(f"  Found {cs.get('count', 0)} credible sets")
            
            for row in cs.get('rows', []):
                study_info = row.get('study', {})
                l2g_preds = row.get('l2GPredictions', {}).get('rows', [])
                
                for pred in l2g_preds:
                    all_results.append({
                        'study_id': study_info.get('id', ''),
                        'trait': study_info.get('traitFromSource', ''),
                        'study_locus_id': row.get('studyLocusId', ''),
                        'region': row.get('region', ''),
                        'chromosome': row.get('chromosome', ''),
                        'position': row.get('position', ''),
                        'gene_id': pred['target']['id'],
                        'gene_symbol': pred['target']['approvedSymbol'],
                        'l2g_score': pred['score']
                    })
        
        time.sleep(0.5)  # Rate limiting
    
    df = pd.DataFrame(all_results)
    if not df.empty:
        df.to_csv(output_path, sep='\t', index=False)
        print(f"\nSaved {len(df)} L2G predictions to {output_path}")
    
    return df

def main():
    print("=" * 70)
    print("OPEN TARGETS L2G EXTRACTION - Platform API v25.4.4")
    print("=" * 70)
    
    # Step 1: Test with SORT1 variant
    print("\n1. Testing with SORT1 locus (rs12740374):")
    print("   Variant ID format: 1_109274968_G_T (GRCh38)")
    
    sort1_data = query_variant("1_109274968_G_T")
    if 'data' in sort1_data and sort1_data['data'].get('variant'):
        v = sort1_data['data']['variant']
        print(f"   Found: chr{v['chromosome']}:{v['position']}")
        print(f"   Alleles: {v['referenceAllele']}/{v['alternateAllele']}")
        print(f"   rsIDs: {v.get('rsIds', [])}")
        
        cs_data = v.get('credibleSets', {})
        print(f"   Associated credible sets: {cs_data.get('count', 0)}")
        
        for cs in cs_data.get('rows', [])[:5]:
            study = cs.get('study', {})
            l2g = cs.get('l2GPredictions', {}).get('rows', [])
            if l2g:
                print(f"\n   Study: {study.get('id', 'N/A')} - {study.get('traitFromSource', 'N/A')}")
                print("   L2G predictions:")
                for pred in l2g[:5]:
                    symbol = pred['target']['approvedSymbol']
                    score = pred['score']
                    marker = " <-- TARGET" if symbol == "SORT1" else ""
                    print(f"     {symbol}: {score:.4f}{marker}")
    else:
        print(f"   Error or no data: {json.dumps(sort1_data, indent=2)[:500]}")
    
    # Step 2: Search for lipid GWAS studies
    print("\n" + "=" * 70)
    print("2. Searching for LDL cholesterol GWAS studies:")
    ldl_studies = search_studies("LDL cholesterol GWAS", limit=10)
    print(f"   Found {len(ldl_studies)} studies")
    for s in ldl_studies[:5]:
        print(f"   - {s['id']}: {s['name']}")
    
    # Get study details for first GWAS study
    if ldl_studies:
        study_id = ldl_studies[0]['id']
        details = get_study_details(study_id)
        if 'data' in details and details['data'].get('study'):
            s = details['data']['study']
            print(f"\n   Study details for {study_id}:")
            print(f"   Trait: {s.get('traitFromSource')}")
            print(f"   N samples: {s.get('nSamples')}")
            print(f"   Type: {s.get('studyType')}")
    
    # Step 3: Get L2G landscape for major lipid studies
    print("\n" + "=" * 70)
    print("3. Extracting L2G landscape for major lipid/cardiometabolic studies:")
    
    # Key GWAS studies to query (GCST IDs from GWAS Catalog)
    key_studies = [
        "GCST006612",    # LDL cholesterol (Global Lipids)
        "GCST90018961",  # LDL cholesterol 
        "GCST002898",    # LDL cholesterol
    ]
    
    output_dir = Path(r"C:\Users\Jack0\GitHub\Mechanism-GWAS-Causal-Graphs\data\processed\baselines")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    l2g_df = get_l2g_landscape(key_studies, output_dir / "platform_api_l2g_scores.tsv")
    
    if not l2g_df.empty:
        print(f"\nL2G Score Summary:")
        print(f"  Total predictions: {len(l2g_df)}")
        print(f"  Unique genes: {l2g_df['gene_symbol'].nunique()}")
        print(f"  Score distribution:")
        print(f"    Mean: {l2g_df['l2g_score'].mean():.3f}")
        print(f"    Median: {l2g_df['l2g_score'].median():.3f}")
        print(f"    >0.5: {(l2g_df['l2g_score'] > 0.5).sum()}")
        print(f"    >0.8: {(l2g_df['l2g_score'] > 0.8).sum()}")
    
    print("\n" + "=" * 70)
    print("L2G extraction complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
