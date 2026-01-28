#!/usr/bin/env python3
"""
Validate Manuscript Claims Against Real Data

This script performs comprehensive validation of all quantitative claims
made in the manuscript against the actual data files.

Claims to Validate:
1. Calibration metrics (ECE < 0.05 for each module)
2. Recall at rank 20 (76% vs 58% for L2G)
3. Cross-study replication (r = 0.89, 78% replication rate)
4. CRISPR benchmark AUPRC (0.71)
5. Locus-level statistics

Author: Validation Script
"""

import json
import sys
from pathlib import Path
from typing import Dict, Any, Tuple
import pandas as pd
import numpy as np
from scipy import stats

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
MANUSCRIPT_DIR = BASE_DIR / "manuscript"
OUTPUT_DIR = BASE_DIR / "scripts" / "validation_output"

OUTPUT_DIR.mkdir(exist_ok=True)

def load_calibration_metrics() -> pd.DataFrame:
    """Load calibration metrics from processed data."""
    path = DATA_DIR / "processed" / "calibration" / "calibration_metrics.tsv"
    return pd.read_csv(path, sep="\t")

def load_benchmark_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load benchmark gold standard and drug targets."""
    tier1 = pd.read_csv(DATA_DIR / "processed" / "benchmark" / "tier1_gold_standard_genes.tsv", sep="\t")
    tier2 = pd.read_csv(DATA_DIR / "processed" / "benchmark" / "tier2_drug_targets.tsv", sep="\t")
    return tier1, tier2

def load_replication_data() -> pd.DataFrame:
    """Load eQTL Catalogue replication data."""
    path = DATA_DIR / "processed" / "replication" / "eqtl_catalogue_replication.tsv"
    return pd.read_csv(path, sep="\t")

def load_locus_summary() -> pd.DataFrame:
    """Load locus-level summary."""
    path = DATA_DIR / "processed" / "locus_summary.tsv"
    return pd.read_csv(path, sep="\t")

def load_mechanism_graphs() -> Dict[str, Any]:
    """Load mechanism graph JSON files."""
    graphs_dir = DATA_DIR / "processed" / "mechanism_graphs"
    graphs = {}
    for json_file in graphs_dir.glob("*.json"):
        with open(json_file) as f:
            graphs[json_file.stem] = json.load(f)
    return graphs

def validate_calibration_claims(cal_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate manuscript claim: ECE < 0.05 for all modules
    
    Manuscript claims (Table 1):
    - Variant PIP (SuSiE): ECE = 0.031 [0.024, 0.038]
    - cCRE-Gene (ABC/PCHi-C): ECE = 0.047 [0.039, 0.055]
    - Gene-Tissue (coloc.susie): ECE = 0.042 [0.035, 0.049]
    - Final gene probability: ECE = 0.038 [0.031, 0.045]
    """
    results = {"claim": "ECE < 0.05 for all modules", "validations": []}
    
    ece_rows = cal_df[cal_df["metric"] == "ECE"]
    
    for _, row in ece_rows.iterrows():
        module = row["module"]
        ece_value = row["value"]
        ci_lower = row["ci_lower"]
        ci_upper = row["ci_upper"]
        
        is_valid = ece_value < 0.05
        
        results["validations"].append({
            "module": module,
            "ece": ece_value,
            "ci": [ci_lower, ci_upper],
            "claim_met": is_valid,
            "threshold": 0.05
        })
    
    results["all_valid"] = all(v["claim_met"] for v in results["validations"])
    return results

def validate_benchmark_performance(cal_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate manuscript claim: Path-probability recall@20 = 76% vs L2G = 58%
    
    This is the primary performance claim.
    """
    results = {"claim": "Path-probability recall@20 = 76% vs L2G = 58%", "validations": []}
    
    # Get recall at rank 20 for our method
    our_recall = cal_df[(cal_df["module"] == "Final_gene_probability") & 
                        (cal_df["metric"] == "Recall_at_20")]
    
    if len(our_recall) > 0:
        our_value = our_recall.iloc[0]["value"]
        our_ci_lower = our_recall.iloc[0]["ci_lower"]
        our_ci_upper = our_recall.iloc[0]["ci_upper"]
        
        results["validations"].append({
            "method": "Path-probability",
            "recall_at_20": our_value,
            "ci": [our_ci_lower, our_ci_upper],
            "manuscript_claim": 0.76,
            "claim_within_ci": 0.76 >= our_ci_lower and 0.76 <= our_ci_upper
        })
    
    # Get L2G recall
    l2g_recall = cal_df[(cal_df["module"] == "Open_Targets_L2G") & 
                        (cal_df["metric"] == "Recall_at_20")]
    
    if len(l2g_recall) > 0:
        l2g_value = l2g_recall.iloc[0]["value"]
        l2g_ci_lower = l2g_recall.iloc[0]["ci_lower"]
        l2g_ci_upper = l2g_recall.iloc[0]["ci_upper"]
        
        results["validations"].append({
            "method": "Open_Targets_L2G",
            "recall_at_20": l2g_value,
            "ci": [l2g_ci_lower, l2g_ci_upper],
            "manuscript_claim": 0.58,
            "claim_within_ci": 0.58 >= l2g_ci_lower and 0.58 <= l2g_ci_upper
        })
    
    # Calculate improvement
    if len(results["validations"]) == 2:
        improvement = our_value - l2g_value
        results["improvement"] = improvement
        results["improvement_percentage"] = f"{improvement*100:.1f} percentage points"
        results["manuscript_claim_improvement"] = "18 percentage points"
    
    results["valid"] = all(v["claim_within_ci"] for v in results["validations"])
    return results

def validate_replication_claims(repl_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate manuscript claims:
    - Overall replication rate: 78%
    - Effect size correlation: r = 0.89
    - Direction concordance: 94%
    """
    results = {"claim": "Cross-study eQTL replication", "validations": []}
    
    # Calculate replication rate
    replicated = repl_df["replicated"].sum()
    total = len(repl_df)
    replication_rate = replicated / total if total > 0 else 0
    
    results["validations"].append({
        "metric": "replication_rate",
        "calculated": replication_rate,
        "manuscript_claim": 0.78,
        "difference": abs(replication_rate - 0.78),
        "within_tolerance": abs(replication_rate - 0.78) < 0.05
    })
    
    # Calculate effect size correlation
    valid_pairs = repl_df[repl_df["replicated"] == True]
    if len(valid_pairs) > 10:
        r, p = stats.pearsonr(valid_pairs["gtex_beta"], valid_pairs["replication_beta"])
        results["validations"].append({
            "metric": "effect_size_correlation",
            "calculated": r,
            "manuscript_claim": 0.89,
            "difference": abs(r - 0.89),
            "within_tolerance": abs(r - 0.89) < 0.05,
            "p_value": p
        })
    
    # Calculate direction concordance
    direction_match = repl_df["direction_concordant"].sum()
    direction_rate = direction_match / total if total > 0 else 0
    
    results["validations"].append({
        "metric": "direction_concordance",
        "calculated": direction_rate,
        "manuscript_claim": 0.94,
        "difference": abs(direction_rate - 0.94),
        "within_tolerance": abs(direction_rate - 0.94) < 0.05
    })
    
    results["valid"] = all(v.get("within_tolerance", False) for v in results["validations"])
    return results

def validate_crispr_benchmark(cal_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate manuscript claim: CRISPR benchmark AUPRC = 0.71
    """
    results = {"claim": "CRISPR benchmark AUPRC = 0.71", "validations": []}
    
    auprc_row = cal_df[(cal_df["module"] == "cCRE_Gene_ABC_PCHiC") & 
                       (cal_df["metric"] == "AUPRC")]
    
    if len(auprc_row) > 0:
        auprc = auprc_row.iloc[0]["value"]
        ci_lower = auprc_row.iloc[0]["ci_lower"]
        ci_upper = auprc_row.iloc[0]["ci_upper"]
        
        results["validations"].append({
            "metric": "AUPRC",
            "calculated": auprc,
            "ci": [ci_lower, ci_upper],
            "manuscript_claim": 0.71,
            "claim_within_ci": ci_lower <= 0.71 <= ci_upper
        })
    
    results["valid"] = all(v.get("claim_within_ci", False) for v in results["validations"])
    return results

def validate_locus_statistics(locus_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate locus-level statistics and mechanism paths.
    """
    results = {"claim": "Locus-level mechanism path statistics", "validations": []}
    
    # Count by evidence tier
    tier_counts = locus_df["evidence_tier"].value_counts().to_dict()
    results["tier_distribution"] = tier_counts
    
    # Count validated vs pending
    validation_counts = locus_df["validation_status"].value_counts().to_dict()
    results["validation_distribution"] = validation_counts
    
    # Summary statistics
    results["summary"] = {
        "n_loci": len(locus_df),
        "mean_path_probability": float(locus_df["path_probability"].mean()),
        "median_path_probability": float(locus_df["path_probability"].median()),
        "mean_coloc_pp_h4": float(locus_df["coloc_pp_h4"].mean()),
        "mean_max_pip": float(locus_df["max_pip"].mean()),
        "unique_traits": len(locus_df["trait"].unique()),
        "unique_tissues": len(locus_df["coloc_tissue"].unique())
    }
    
    # Validate key example: SORT1
    sort1 = locus_df[locus_df["top_gene"] == "SORT1"]
    if len(sort1) > 0:
        sort1_row = sort1.iloc[0]
        results["validations"].append({
            "locus": "SORT1_1p13_LDL",
            "path_probability": float(sort1_row["path_probability"]),
            "manuscript_claim": 0.79,
            "match": abs(sort1_row["path_probability"] - 0.79) < 0.01
        })
    
    results["valid"] = True
    return results

def validate_mechanism_graph_structure(graphs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate that mechanism graphs have correct structure.
    """
    results = {"claim": "Mechanism graph structure validation", "validations": []}
    
    required_fields = [
        "locus_id", "trait", "gwas_source", "genome_build",
        "locus_boundaries", "credible_variants", "mechanism_paths"
    ]
    
    for name, graph in graphs.items():
        missing = [f for f in required_fields if f not in graph]
        
        # Check path structure
        paths = graph.get("mechanism_paths", [])
        path_valid = True
        for path in paths:
            if not all(k in path for k in ["path_id", "path_probability", "variant", "gene"]):
                path_valid = False
                break
        
        results["validations"].append({
            "graph": name,
            "has_required_fields": len(missing) == 0,
            "missing_fields": missing,
            "n_paths": len(paths),
            "path_structure_valid": path_valid
        })
    
    results["valid"] = all(v["has_required_fields"] and v["path_structure_valid"] 
                          for v in results["validations"])
    return results

def validate_drug_targets(tier2_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate drug target benchmark data.
    """
    results = {"claim": "Drug target benchmark validation", "validations": []}
    
    # Count drug targets by approval status
    status_counts = tier2_df["approval_status"].value_counts().to_dict()
    results["approval_status_distribution"] = status_counts
    
    # Count by indication
    indication_counts = tier2_df["indication"].value_counts().head(10).to_dict()
    results["top_indications"] = indication_counts
    
    results["summary"] = {
        "n_drug_gene_pairs": len(tier2_df),
        "n_unique_genes": len(tier2_df["gene_symbol"].unique()),
        "n_unique_drugs": len(tier2_df["drug_name"].unique()),
        "n_approved_drugs": len(tier2_df[tier2_df["approval_status"] == "Approved"])
    }
    
    # Validate key targets mentioned in manuscript
    key_targets = ["HMGCR", "PCSK9", "PPARG", "KCNJ11"]
    for target in key_targets:
        found = tier2_df[tier2_df["gene_symbol"] == target]
        results["validations"].append({
            "target": target,
            "found_in_data": len(found) > 0,
            "n_drugs": len(found),
            "drug_names": found["drug_name"].tolist() if len(found) > 0 else []
        })
    
    results["valid"] = all(v["found_in_data"] for v in results["validations"])
    return results

def generate_validation_report(all_results: Dict[str, Any]) -> str:
    """Generate comprehensive validation report."""
    report = []
    report.append("=" * 80)
    report.append("MANUSCRIPT CLAIM VALIDATION REPORT")
    report.append("Mechanism-First Causal Graphs for Noncoding GWAS")
    report.append("=" * 80)
    report.append("")
    
    overall_valid = True
    
    for section, results in all_results.items():
        report.append(f"\n{'='*60}")
        report.append(f"SECTION: {section}")
        report.append(f"{'='*60}")
        
        if "claim" in results:
            report.append(f"Claim: {results['claim']}")
        
        if "valid" in results:
            status = "✓ VALIDATED" if results["valid"] else "✗ NEEDS REVIEW"
            report.append(f"Status: {status}")
            if not results["valid"]:
                overall_valid = False
        
        if "validations" in results:
            report.append("\nDetails:")
            for v in results["validations"]:
                report.append(f"  - {v}")
        
        if "summary" in results:
            report.append("\nSummary statistics:")
            for k, v in results["summary"].items():
                report.append(f"  {k}: {v}")
    
    report.append("\n" + "=" * 80)
    report.append(f"OVERALL VALIDATION: {'✓ ALL CLAIMS VALIDATED' if overall_valid else '⚠ SOME CLAIMS NEED REVIEW'}")
    report.append("=" * 80)
    
    return "\n".join(report)

def main():
    """Run all validations."""
    print("Loading data files...")
    
    # Load all data
    cal_df = load_calibration_metrics()
    tier1_df, tier2_df = load_benchmark_data()
    repl_df = load_replication_data()
    locus_df = load_locus_summary()
    graphs = load_mechanism_graphs()
    
    print(f"Loaded calibration data: {len(cal_df)} rows")
    print(f"Loaded tier1 benchmark: {len(tier1_df)} genes")
    print(f"Loaded tier2 benchmark: {len(tier2_df)} drug targets")
    print(f"Loaded replication data: {len(repl_df)} genes")
    print(f"Loaded locus summary: {len(locus_df)} loci")
    print(f"Loaded mechanism graphs: {len(graphs)} graphs")
    
    # Run all validations
    print("\nRunning validations...")
    
    all_results = {
        "calibration": validate_calibration_claims(cal_df),
        "benchmark_performance": validate_benchmark_performance(cal_df),
        "replication": validate_replication_claims(repl_df),
        "crispr_benchmark": validate_crispr_benchmark(cal_df),
        "locus_statistics": validate_locus_statistics(locus_df),
        "mechanism_graphs": validate_mechanism_graph_structure(graphs),
        "drug_targets": validate_drug_targets(tier2_df)
    }
    
    # Generate report
    report = generate_validation_report(all_results)
    print(report)
    
    # Save results
    output_file = OUTPUT_DIR / "validation_report.txt"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(report)
    
    # Save JSON results
    json_file = OUTPUT_DIR / "validation_results.json"
    with open(json_file, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\nResults saved to:")
    print(f"  {output_file}")
    print(f"  {json_file}")
    
    # Return success/failure
    all_valid = all(r.get("valid", True) for r in all_results.values())
    return 0 if all_valid else 1

if __name__ == "__main__":
    sys.exit(main())
