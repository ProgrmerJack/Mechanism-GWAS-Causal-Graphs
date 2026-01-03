"""
Analyze baseline comparison case studies
Identify representative loci where different methods succeed/fail
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent
RESULTS_DIR = BASE_DIR / "results" / "baselines"
BENCHMARK_FILE = BASE_DIR / "data" / "processed" / "baselines" / "post2021_independent_benchmark_FINAL.tsv"
PREDICTIONS_FILE = RESULTS_DIR / "post2021_predictions_all_methods.tsv"
OUTPUT_FILE = RESULTS_DIR / "CASE_STUDY_ANALYSIS.md"

def load_data():
    """Load benchmark and predictions"""
    benchmark = pd.read_csv(BENCHMARK_FILE, sep='\t')
    predictions = pd.read_csv(PREDICTIONS_FILE, sep='\t')
    return benchmark, predictions

def find_distance_success_functional_fail(predictions):
    """Find loci where Distance succeeds but functional methods fail"""
    cases = []
    
    for locus_id in predictions['locus_id'].unique():
        locus_preds = predictions[predictions['locus_id'] == locus_id]
        
        # Get ranks
        distance_rank = locus_preds[locus_preds['method'] == 'Distance']['true_gene_rank'].values[0]
        abc_rank = locus_preds[locus_preds['method'] == 'ABC_Only']['true_gene_rank'].values[0]
        eqtl_rank = locus_preds[locus_preds['method'] == 'eQTL_Only']['true_gene_rank'].values[0]
        
        # Distance Top-1, functional methods Top-10 or worse
        if distance_rank == 1 and abc_rank > 5 and eqtl_rank > 5:
            evidence_tier = locus_preds['evidence_tier'].values[0]
            trait = locus_preds['trait_category'].values[0]
            true_gene = locus_preds['true_gene'].values[0]
            
            cases.append({
                'locus_id': locus_id,
                'true_gene': true_gene,
                'evidence_tier': evidence_tier,
                'trait_category': trait,
                'distance_rank': distance_rank,
                'abc_rank': abc_rank,
                'eqtl_rank': eqtl_rank,
                'case_type': 'Distance Success / Functional Fail'
            })
    
    return pd.DataFrame(cases)

def find_functional_success_distance_fail(predictions):
    """Find loci where functional methods succeed but Distance fails"""
    cases = []
    
    for locus_id in predictions['locus_id'].unique():
        locus_preds = predictions[predictions['locus_id'] == locus_id]
        
        # Get ranks
        distance_rank = locus_preds[locus_preds['method'] == 'Distance']['true_gene_rank'].values[0]
        abc_rank = locus_preds[locus_preds['method'] == 'ABC_Only']['true_gene_rank'].values[0]
        eqtl_rank = locus_preds[locus_preds['method'] == 'eQTL_Only']['true_gene_rank'].values[0]
        
        # Distance fails Top-3, at least one functional method Top-3
        if distance_rank > 3 and (abc_rank <= 3 or eqtl_rank <= 3):
            evidence_tier = locus_preds['evidence_tier'].values[0]
            trait = locus_preds['trait_category'].values[0]
            true_gene = locus_preds['true_gene'].values[0]
            
            cases.append({
                'locus_id': locus_id,
                'true_gene': true_gene,
                'evidence_tier': evidence_tier,
                'trait_category': trait,
                'distance_rank': distance_rank,
                'abc_rank': abc_rank,
                'eqtl_rank': eqtl_rank,
                'case_type': 'Functional Success / Distance Fail'
            })
    
    return pd.DataFrame(cases)

def find_all_methods_fail(predictions):
    """Find loci where all methods fail Top-5"""
    cases = []
    
    for locus_id in predictions['locus_id'].unique():
        locus_preds = predictions[predictions['locus_id'] == locus_id]
        
        # Get ranks for all methods
        distance_rank = locus_preds[locus_preds['method'] == 'Distance']['true_gene_rank'].values[0]
        abc_rank = locus_preds[locus_preds['method'] == 'ABC_Only']['true_gene_rank'].values[0]
        eqtl_rank = locus_preds[locus_preds['method'] == 'eQTL_Only']['true_gene_rank'].values[0]
        pops_rank = locus_preds[locus_preds['method'] == 'PoPS']['true_gene_rank'].values[0]
        
        # All methods fail Top-5
        if min(distance_rank, abc_rank, eqtl_rank, pops_rank) > 5:
            evidence_tier = locus_preds['evidence_tier'].values[0]
            trait = locus_preds['trait_category'].values[0]
            true_gene = locus_preds['true_gene'].values[0]
            total_candidates = locus_preds['total_candidates'].values[0]
            
            cases.append({
                'locus_id': locus_id,
                'true_gene': true_gene,
                'evidence_tier': evidence_tier,
                'trait_category': trait,
                'distance_rank': distance_rank,
                'abc_rank': abc_rank,
                'eqtl_rank': eqtl_rank,
                'pops_rank': pops_rank,
                'total_candidates': total_candidates,
                'case_type': 'All Methods Fail'
            })
    
    return pd.DataFrame(cases)

def find_all_methods_succeed(predictions):
    """Find loci where all methods succeed Top-3"""
    cases = []
    
    for locus_id in predictions['locus_id'].unique():
        locus_preds = predictions[predictions['locus_id'] == locus_id]
        
        # Get ranks for all methods
        distance_rank = locus_preds[locus_preds['method'] == 'Distance']['true_gene_rank'].values[0]
        abc_rank = locus_preds[locus_preds['method'] == 'ABC_Only']['true_gene_rank'].values[0]
        eqtl_rank = locus_preds[locus_preds['method'] == 'eQTL_Only']['true_gene_rank'].values[0]
        pops_rank = locus_preds[locus_preds['method'] == 'PoPS']['true_gene_rank'].values[0]
        
        # All methods succeed Top-3
        if max(distance_rank, abc_rank, eqtl_rank, pops_rank) <= 3:
            evidence_tier = locus_preds['evidence_tier'].values[0]
            trait = locus_preds['trait_category'].values[0]
            true_gene = locus_preds['true_gene'].values[0]
            
            cases.append({
                'locus_id': locus_id,
                'true_gene': true_gene,
                'evidence_tier': evidence_tier,
                'trait_category': trait,
                'distance_rank': distance_rank,
                'abc_rank': abc_rank,
                'eqtl_rank': eqtl_rank,
                'pops_rank': pops_rank,
                'case_type': 'All Methods Succeed'
            })
    
    return pd.DataFrame(cases)

def get_locus_details(locus_id, benchmark):
    """Get detailed information about a locus from benchmark"""
    locus = benchmark[benchmark['locus_id'] == locus_id].iloc[0]
    
    details = {
        'locus_id': locus['locus_id'],
        'true_gene': locus['gene_symbol'],  # Column is gene_symbol
        'trait': locus['trait'],
        'chr': locus['chr'],
        'pos': locus['pos_hg38'],
        'rsid': locus.get('lead_snp', 'N/A'),
        'evidence_tier': locus['evidence_tier'],
        'validation_evidence': locus['validation_type'],
        'trait_category': locus['trait_category'],
        'pubmed_id': locus['validation_pmid'],
        'publication_year': locus['gwas_pmid']  # Using gwas_pmid instead
    }
    
    return details

def create_case_study_report(benchmark, predictions):
    """Generate comprehensive case study report"""
    
    # Find different case types
    print("Finding case studies...")
    distance_success = find_distance_success_functional_fail(predictions)
    functional_success = find_functional_success_distance_fail(predictions)
    all_fail = find_all_methods_fail(predictions)
    all_succeed = find_all_methods_succeed(predictions)
    
    print(f"Distance success / Functional fail: {len(distance_success)} loci")
    print(f"Functional success / Distance fail: {len(functional_success)} loci")
    print(f"All methods fail: {len(all_fail)} loci")
    print(f"All methods succeed: {len(all_succeed)} loci")
    
    # Create markdown report
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write("# Baseline Comparison Case Study Analysis\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now()}\n\n")
        
        f.write("## Summary\n\n")
        f.write(f"- **Distance success / Functional fail**: {len(distance_success)} loci\n")
        f.write(f"- **Functional success / Distance fail**: {len(functional_success)} loci\n")
        f.write(f"- **All methods fail (Top-5)**: {len(all_fail)} loci\n")
        f.write(f"- **All methods succeed (Top-3)**: {len(all_succeed)} loci\n\n")
        
        # Case Type 1: Distance Success / Functional Fail
        f.write("## Case Type 1: Distance Success / Functional Fail\n\n")
        f.write("These loci represent **Mendelian/coding variants** where the causal gene is typically ")
        f.write("the nearest gene. Functional methods (ABC, eQTL) fail because they prioritize regulatory ")
        f.write("mechanisms over coding/loss-of-function variants.\n\n")
        
        if len(distance_success) > 0:
            # Select top 5 examples
            examples = distance_success.nlargest(5, 'abc_rank')
            
            f.write("### Representative Examples\n\n")
            for idx, row in examples.iterrows():
                details = get_locus_details(row['locus_id'], benchmark)
                
                f.write(f"#### {details['locus_id']}\n\n")
                f.write(f"- **True Gene**: {details['true_gene']}\n")
                f.write(f"- **Trait**: {details['trait']}\n")
                f.write(f"- **Evidence Tier**: {details['evidence_tier']}\n")
                f.write(f"- **Validation**: {details['validation_evidence']}\n")
                f.write(f"- **PubMed**: {details['pubmed_id']}\n\n")
                
                f.write("**Method Rankings:**\n")
                f.write(f"- Distance: **Rank {row['distance_rank']}** ✓\n")
                f.write(f"- ABC: Rank {row['abc_rank']}\n")
                f.write(f"- eQTL: Rank {row['eqtl_rank']}\n\n")
                
                f.write("**Interpretation:**\n")
                if 'Mendelian' in details['evidence_tier']:
                    f.write(f"- Mendelian disease gene ({details['true_gene']}) with definitive causative evidence\n")
                    f.write("- Likely harbors coding/LOF variant disrupting protein function\n")
                    f.write("- Distance prioritizes correctly (nearest gene = causal gene)\n")
                    f.write("- ABC/eQTL prioritize distal regulatory targets incorrectly\n\n")
                elif 'Coding' in details['evidence_tier']:
                    f.write(f"- Damaging coding variant in {details['true_gene']}\n")
                    f.write("- Direct protein-level mechanism, not regulatory\n")
                    f.write("- Distance succeeds, functional methods fail\n\n")
                
                f.write("---\n\n")
        
        # Case Type 2: Functional Success / Distance Fail
        f.write("## Case Type 2: Functional Success / Distance Fail\n\n")
        f.write("These loci represent **regulatory variants** acting through distal enhancers ")
        f.write("to modulate gene expression. The causal gene is NOT the nearest gene, requiring ")
        f.write("functional annotations (ABC, eQTL) to identify the correct target.\n\n")
        
        if len(functional_success) > 0:
            f.write("### Representative Examples\n\n")
            for idx, row in functional_success.iterrows():
                details = get_locus_details(row['locus_id'], benchmark)
                
                f.write(f"#### {details['locus_id']}\n\n")
                f.write(f"- **True Gene**: {details['true_gene']}\n")
                f.write(f"- **Trait**: {details['trait']}\n")
                f.write(f"- **Evidence Tier**: {details['evidence_tier']}\n")
                f.write(f"- **Validation**: {details['validation_evidence']}\n")
                f.write(f"- **PubMed**: {details['pubmed_id']}\n\n")
                
                f.write("**Method Rankings:**\n")
                f.write(f"- Distance: Rank {row['distance_rank']}\n")
                if row['abc_rank'] <= 3:
                    f.write(f"- ABC: **Rank {row['abc_rank']}** ✓\n")
                else:
                    f.write(f"- ABC: Rank {row['abc_rank']}\n")
                
                if row['eqtl_rank'] <= 3:
                    f.write(f"- eQTL: **Rank {row['eqtl_rank']}** ✓\n")
                else:
                    f.write(f"- eQTL: Rank {row['eqtl_rank']}\n")
                f.write("\n")
                
                f.write("**Interpretation:**\n")
                if 'CRISPR' in details['evidence_tier']:
                    f.write(f"- CRISPR screen validated {details['true_gene']} as causal\n")
                    f.write("- Likely regulatory variant affecting enhancer activity\n")
                    f.write("- Functional methods prioritize correctly\n")
                    f.write("- Distance fails (true gene not nearest)\n\n")
                else:
                    f.write(f"- Regulatory mechanism affecting {details['true_gene']}\n")
                    f.write("- Functional annotations needed to identify target\n\n")
                
                f.write("---\n\n")
        
        # Case Type 3: All Methods Fail
        f.write("## Case Type 3: All Methods Fail\n\n")
        f.write("These are **challenging loci** requiring mechanistic causal chain integration ")
        f.write("beyond simple functional scores. Possible explanations:\n")
        f.write("- Complex regulatory architecture (multiple enhancers)\n")
        f.write("- Long-range interactions not captured by ABC/eQTL\n")
        f.write("- Context-dependent mechanisms (cell type, developmental stage)\n")
        f.write("- Pathway-level effects requiring systems biology\n\n")
        
        if len(all_fail) > 0:
            f.write("### Representative Examples\n\n")
            for idx, row in all_fail.iterrows():
                details = get_locus_details(row['locus_id'], benchmark)
                
                f.write(f"#### {details['locus_id']}\n\n")
                f.write(f"- **True Gene**: {details['true_gene']}\n")
                f.write(f"- **Trait**: {details['trait']}\n")
                f.write(f"- **Evidence Tier**: {details['evidence_tier']}\n")
                f.write(f"- **Validation**: {details['validation_evidence']}\n")
                f.write(f"- **Candidate Genes**: {row['total_candidates']}\n")
                f.write(f"- **PubMed**: {details['pubmed_id']}\n\n")
                
                f.write("**Method Rankings (All Fail):**\n")
                f.write(f"- Distance: Rank {row['distance_rank']} / {row['total_candidates']}\n")
                f.write(f"- ABC: Rank {row['abc_rank']} / {row['total_candidates']}\n")
                f.write(f"- eQTL: Rank {row['eqtl_rank']} / {row['total_candidates']}\n")
                f.write(f"- PoPS: Rank {row['pops_rank']} / {row['total_candidates']}\n\n")
                
                f.write("**Why This is Challenging:**\n")
                f.write("- True gene not ranked highly by any baseline method\n")
                f.write("- Requires integrated mechanistic evidence (variant → cCRE → gene → pathway)\n")
                f.write("- Our causal chain approach needed to prioritize correctly\n\n")
                
                f.write("---\n\n")
        
        # Case Type 4: All Methods Succeed
        f.write("## Case Type 4: All Methods Succeed\n\n")
        f.write("These loci have **strong convergent evidence** where multiple independent ")
        f.write("signals agree on the causal gene. High confidence predictions.\n\n")
        
        if len(all_succeed) > 0:
            f.write("### Representative Examples\n\n")
            for idx, row in all_succeed.iterrows():
                details = get_locus_details(row['locus_id'], benchmark)
                
                f.write(f"#### {details['locus_id']}\n\n")
                f.write(f"- **True Gene**: {details['true_gene']}\n")
                f.write(f"- **Trait**: {details['trait']}\n")
                f.write(f"- **Evidence Tier**: {details['evidence_tier']}\n")
                f.write(f"- **Validation**: {details['validation_evidence']}\n\n")
                
                f.write("**Method Rankings (All Succeed):**\n")
                f.write(f"- Distance: Rank {row['distance_rank']}\n")
                f.write(f"- ABC: Rank {row['abc_rank']}\n")
                f.write(f"- eQTL: Rank {row['eqtl_rank']}\n")
                f.write(f"- PoPS: Rank {row['pops_rank']}\n\n")
                
                f.write("**Interpretation:**\n")
                f.write(f"- Strong convergent evidence for {details['true_gene']}\n")
                f.write("- High confidence across multiple prioritization strategies\n\n")
                
                f.write("---\n\n")
        
        # Recommendations
        f.write("## Recommendations for Manuscript\n\n")
        f.write("### Case Study Selection\n\n")
        f.write("We recommend including 3-5 detailed case studies in the manuscript:\n\n")
        
        f.write("1. **Mendelian Example** (Distance success): Show coding/LOF variant in nearest gene\n")
        f.write("2. **Regulatory Example** (Functional success): Show distal enhancer targeting non-nearest gene\n")
        f.write("3. **Complex Example** (All fail): Demonstrate need for mechanistic causal chains\n")
        f.write("4. **Convergent Example** (All succeed): Show high-confidence multi-method agreement\n\n")
        
        f.write("### Mechanistic Interpretation\n\n")
        f.write("For each case study, we should provide:\n")
        f.write("- Variant → cCRE mapping (candidate regulatory elements)\n")
        f.write("- cCRE → Gene links (ABC enhancer-gene predictions, eQTL colocalization)\n")
        f.write("- Gene → Pathway connections (biological mechanism)\n")
        f.write("- Network diagram visualizing the causal chain\n\n")
        
        f.write("This demonstrates how our method integrates multiple evidence types ")
        f.write("that baseline methods use in isolation.\n\n")
    
    print(f"\n✓ Case study analysis saved to: {OUTPUT_FILE}")
    
    return {
        'distance_success': distance_success,
        'functional_success': functional_success,
        'all_fail': all_fail,
        'all_succeed': all_succeed
    }

def main():
    """Main execution"""
    print("Loading data...")
    benchmark, predictions = load_data()
    
    print(f"Benchmark: {len(benchmark)} loci")
    print(f"Predictions: {len(predictions)} rows ({len(predictions['locus_id'].unique())} loci × {len(predictions['method'].unique())} methods)")
    
    # Create case study report
    case_studies = create_case_study_report(benchmark, predictions)
    
    print("\nCase study analysis complete!")
    print(f"\nSummary:")
    print(f"  Distance success: {len(case_studies['distance_success'])} loci")
    print(f"  Functional success: {len(case_studies['functional_success'])} loci")
    print(f"  All methods fail: {len(case_studies['all_fail'])} loci")
    print(f"  All methods succeed: {len(case_studies['all_succeed'])} loci")

if __name__ == "__main__":
    main()
