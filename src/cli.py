"""
CLI Module for Mechanism-First Causal Graphs
"""

import argparse
import sys


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Mechanism-First Causal Graphs for Noncoding GWAS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run full pipeline
    mechanism-gwas run --config config/config.yaml
    
    # Fine-map a single trait
    mechanism-gwas finemap --trait ldl_cholesterol --output results/
    
    # Build mechanism graphs
    mechanism-gwas graph --pips results/pips.tsv --coloc results/coloc.tsv
    
    # Run calibration analysis
    mechanism-gwas calibrate --rankings results/rankings.tsv --trait ldl_cholesterol
        """
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Run command
    run_parser = subparsers.add_parser("run", help="Run full pipeline")
    run_parser.add_argument("--config", required=True, help="Configuration file")
    run_parser.add_argument("--cores", type=int, default=4, help="Number of cores")
    run_parser.add_argument("--dry-run", action="store_true", help="Dry run mode")
    
    # Finemap command
    finemap_parser = subparsers.add_parser("finemap", help="Run fine-mapping")
    finemap_parser.add_argument("--trait", required=True, help="Trait name")
    finemap_parser.add_argument("--sumstats", required=True, help="Summary statistics")
    finemap_parser.add_argument("--ld-dir", required=True, help="LD reference directory")
    finemap_parser.add_argument("--output", required=True, help="Output directory")
    
    # Graph command
    graph_parser = subparsers.add_parser("graph", help="Build mechanism graphs")
    graph_parser.add_argument("--pips", required=True, help="PIP file from fine-mapping")
    graph_parser.add_argument("--coloc", required=True, help="Colocalization results")
    graph_parser.add_argument("--ccre", required=True, help="cCRE annotations")
    graph_parser.add_argument("--output", required=True, help="Output JSON file")
    
    # Calibrate command
    calibrate_parser = subparsers.add_parser("calibrate", help="Run calibration")
    calibrate_parser.add_argument("--rankings", required=True, help="Gene rankings")
    calibrate_parser.add_argument("--trait", required=True, help="Trait name")
    calibrate_parser.add_argument("--output", required=True, help="Output file")
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    
    if args.command == "run":
        from .pipeline import run_pipeline
        run_pipeline(args.config, args.cores, args.dry_run)
    
    elif args.command == "finemap":
        from .finemapping import run_finemapping
        run_finemapping(args.trait, args.sumstats, args.ld_dir, args.output)
    
    elif args.command == "graph":
        from .mechanism_graph import build_graphs
        build_graphs(args.pips, args.coloc, args.ccre, args.output)
    
    elif args.command == "calibrate":
        from .calibration import run_calibration
        run_calibration(args.rankings, args.trait, args.output)


if __name__ == "__main__":
    main()
