"""
Atlas Assembly Script

Assembles all mechanism graphs into a unified atlas for publication.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd


def load_trait_results(
    results_dir: str,
    trait: str,
) -> Dict[str, Any]:
    """Load all results for a trait."""
    trait_dir = Path(results_dir) / trait
    
    results = {"trait": trait}
    
    # Load mechanism graphs
    graphs_file = trait_dir / "mechanism_graphs.json"
    if graphs_file.exists():
        with open(graphs_file) as f:
            results["graphs"] = json.load(f)
    
    # Load gene rankings
    rankings_file = trait_dir / "gene_rankings.tsv"
    if rankings_file.exists():
        rankings = pd.read_csv(rankings_file, sep="\t")
        results["rankings"] = rankings.to_dict("records")
    
    # Load calibration report
    calibration_file = trait_dir / "calibration_report.json"
    if calibration_file.exists():
        with open(calibration_file) as f:
            results["calibration"] = json.load(f)
    
    return results


def assemble_atlas(
    results_dir: str,
    traits: List[str],
) -> Dict[str, Any]:
    """Assemble atlas from all trait results."""
    atlas = {
        "version": "1.0.0",
        "genome_build": "GRCh38",
        "traits": {},
        "metadata": {
            "n_traits": len(traits),
            "pipeline": "Mechanism-First Causal Graphs",
        },
    }
    
    all_genes = {}
    
    for trait in traits:
        print(f"Loading {trait}...")
        results = load_trait_results(results_dir, trait)
        atlas["traits"][trait] = results
        
        # Collect genes across traits
        for gene_info in results.get("rankings", []):
            gene = gene_info.get("gene", "")
            if gene:
                if gene not in all_genes:
                    all_genes[gene] = {"traits": [], "scores": []}
                all_genes[gene]["traits"].append(trait)
                all_genes[gene]["scores"].append(
                    gene_info.get("combined_probability", 0)
                )
    
    # Summary statistics
    atlas["summary"] = {
        "n_genes_prioritized": len(all_genes),
        "genes_by_n_traits": {},
    }
    
    for gene, info in all_genes.items():
        n = len(info["traits"])
        atlas["summary"]["genes_by_n_traits"][n] = (
            atlas["summary"]["genes_by_n_traits"].get(n, 0) + 1
        )
    
    return atlas


def generate_summary_table(
    atlas: Dict[str, Any],
) -> pd.DataFrame:
    """Generate summary table for all genes."""
    rows = []
    
    for trait, trait_data in atlas["traits"].items():
        for gene_info in trait_data.get("rankings", []):
            rows.append({
                "trait": trait,
                "gene": gene_info.get("gene", ""),
                "combined_probability": gene_info.get("combined_probability", 0),
                "ci_lower": gene_info.get("ci_lower", 0),
                "ci_upper": gene_info.get("ci_upper", 0),
                "n_loci": gene_info.get("n_loci", 0),
                "best_tissue": gene_info.get("best_tissue", ""),
            })
    
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Assemble mechanism graph atlas"
    )
    parser.add_argument("--results-dir", required=True, help="Results directory")
    parser.add_argument("--traits", nargs="+", required=True, help="Traits to include")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--summary-output", required=True, help="Summary TSV file")
    
    args = parser.parse_args()
    
    # Assemble atlas
    atlas = assemble_atlas(args.results_dir, args.traits)
    
    # Save atlas
    with open(args.output, "w") as f:
        json.dump(atlas, f, indent=2)
    
    print(f"Saved atlas to {args.output}")
    
    # Generate summary
    summary = generate_summary_table(atlas)
    summary.to_csv(args.summary_output, sep="\t", index=False)
    
    print(f"Saved summary to {args.summary_output}")
    print(f"Total genes: {atlas['summary']['n_genes_prioritized']}")


if __name__ == "__main__":
    main()
