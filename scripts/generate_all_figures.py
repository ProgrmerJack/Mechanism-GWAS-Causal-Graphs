#!/usr/bin/env python3
"""
FLAMES Figure Generation Pipeline
==================================

Main entry point for generating all publication figures.

This script orchestrates the generation of:
- Main Figures 1-4
- Extended Data Figures 1-10
- Composite preview sheet

All figures follow Nature Genetics specifications:
- Okabe-Ito colorblind-safe palette
- Constrained layout (no hardcoded positions)
- 600 DPI for print, PDF with embedded TrueType fonts
- Proper panel labeling (8pt bold lowercase)

Usage:
    python generate_all_figures.py [--figures] [--extended] [--preview] [--all]

Author: FLAMES Project
"""

import sys

# Increase recursion limit for Python 3.14 dateutil import chain
sys.setrecursionlimit(3000)

import argparse
import time
from pathlib import Path
from datetime import datetime

# Add parent to path
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent  # Go up one level to project root
sys.path.insert(0, str(SCRIPT_DIR))

import matplotlib
matplotlib.use('Agg')

# Import figure modules
from figures.style import setup_nature_style, validate_figure_dimensions
from figures.utils import load_all_data, generate_preview_sheet, check_overlaps
from figures.figure1 import create_figure_1
from figures.figure2 import create_figure_2
from figures.figure3 import create_figure_3
from figures.figure4 import create_figure_4
from figures.extended_data import create_all_ed_figures


def print_header():
    """Print script header."""
    print("\n" + "=" * 70)
    print("  FLAMES Figure Generation Pipeline")
    print("  Nature Genetics Submission")
    print("=" * 70)
    print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)


def print_summary(start_time: float, figures_generated: list, output_dir: Path):
    """Print generation summary."""
    elapsed = time.time() - start_time
    
    print("\n" + "=" * 70)
    print("  GENERATION COMPLETE")
    print("=" * 70)
    print(f"  Time elapsed: {elapsed:.1f} seconds")
    print(f"  Figures generated: {len(figures_generated)}")
    print(f"  Output directory: {output_dir}")
    print()
    print("  Generated figures:")
    for fig in figures_generated:
        print(f"    ✓ {fig}")
    print("=" * 70 + "\n")


def validate_outputs(output_dir: Path) -> dict:
    """
    Validate all generated figures.
    
    Returns
    -------
    dict
        Validation results with any issues found
    """
    print("\n  Validating outputs...")
    
    results = {
        'total': 0,
        'valid': 0,
        'issues': []
    }
    
    # Check for expected files
    expected_pdfs = [
        'fig1_calibration_overview.pdf',
        'fig2_stress_tests.pdf',
        'fig3_case_studies.pdf',
        'fig4_benchmark_performance.pdf',
    ] + [f'ed_fig{i}*.pdf' for i in range(1, 11)]
    
    for pattern in expected_pdfs:
        matches = list(output_dir.glob(pattern))
        if matches:
            results['total'] += len(matches)
            for match in matches:
                # Check file size
                size_kb = match.stat().st_size / 1024
                if size_kb < 10:
                    results['issues'].append(f"{match.name}: File too small ({size_kb:.1f} KB)")
                else:
                    results['valid'] += 1
        else:
            results['issues'].append(f"Missing: {pattern}")
    
    # Report results
    if results['issues']:
        print(f"  ⚠ Found {len(results['issues'])} issue(s):")
        for issue in results['issues'][:5]:  # Show first 5
            print(f"    - {issue}")
        if len(results['issues']) > 5:
            print(f"    ... and {len(results['issues']) - 5} more")
    else:
        print(f"  ✓ All {results['valid']} figures validated successfully")
    
    return results


def generate_main_figures(data: dict, output_dir: Path, base_path: Path) -> list:
    """
    Generate main manuscript figures 1-4.
    
    Returns
    -------
    list
        Names of generated figures
    """
    generated = []
    
    print("\n" + "-" * 50)
    print("  MAIN FIGURES (1-4)")
    print("-" * 50)
    
    try:
        create_figure_1(data, output_dir, base_path)
        generated.append('Figure 1 - Calibration Overview')
    except Exception as e:
        print(f"  ✗ Figure 1 failed: {e}")
    
    try:
        create_figure_2(data, output_dir, base_path)
        generated.append('Figure 2 - Stress Tests')
    except Exception as e:
        print(f"  ✗ Figure 2 failed: {e}")
    
    try:
        create_figure_3(data, output_dir, base_path)
        generated.append('Figure 3 - Case Studies')
    except Exception as e:
        print(f"  ✗ Figure 3 failed: {e}")
    
    try:
        create_figure_4(data, output_dir, base_path)
        generated.append('Figure 4 - Benchmark Performance')
    except Exception as e:
        print(f"  ✗ Figure 4 failed: {e}")
    
    return generated


def generate_extended_data(data: dict, output_dir: Path, base_path: Path) -> list:
    """
    Generate Extended Data figures 1-10.
    
    Returns
    -------
    list
        Names of generated figures
    """
    generated = []
    
    print("\n" + "-" * 50)
    print("  EXTENDED DATA FIGURES (1-10)")
    print("-" * 50)
    
    try:
        create_all_ed_figures(data, output_dir, base_path)
        generated.extend([f'ED Figure {i}' for i in range(1, 11)])
    except Exception as e:
        print(f"  ✗ Extended Data figures failed: {e}")
    
    return generated


def generate_composite_preview(output_dir: Path) -> str:
    """
    Generate composite preview sheet.
    
    Returns
    -------
    str
        Path to preview file
    """
    print("\n" + "-" * 50)
    print("  COMPOSITE PREVIEW SHEET")
    print("-" * 50)
    
    # Collect all PNG files
    png_files = sorted(output_dir.glob('*.png'))
    
    if not png_files:
        print("  ⚠ No PNG files found for preview")
        return None
    
    preview_path = output_dir / 'composite_preview.png'
    
    try:
        generate_preview_sheet(
            png_files,  # Pass Path objects, not strings
            preview_path,
            title='FLAMES Manuscript Figures - Preview Sheet'
        )
        print(f"  ✓ Preview generated: {preview_path.name}")
        return str(preview_path)
    except Exception as e:
        print(f"  ✗ Preview generation failed: {e}")
        return None


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Generate FLAMES publication figures'
    )
    parser.add_argument(
        '--figures', '-f', action='store_true',
        help='Generate main figures only (1-4)'
    )
    parser.add_argument(
        '--extended', '-e', action='store_true',
        help='Generate Extended Data figures only (1-10)'
    )
    parser.add_argument(
        '--preview', '-p', action='store_true',
        help='Generate composite preview only'
    )
    parser.add_argument(
        '--all', '-a', action='store_true',
        help='Generate all figures (default)'
    )
    parser.add_argument(
        '--validate', '-v', action='store_true',
        help='Validate outputs after generation'
    )
    parser.add_argument(
        '--output', '-o', type=str, default=None,
        help='Output directory (default: figures_professional)'
    )
    
    args = parser.parse_args()
    
    # Default to all if nothing specified
    if not any([args.figures, args.extended, args.preview]):
        args.all = True
    
    # Setup paths - PROJECT_ROOT contains results/
    base_path = PROJECT_ROOT
    
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = base_path / 'figures_professional'
    
    output_dir.mkdir(exist_ok=True)
    
    # Print header
    print_header()
    print(f"\n  Output directory: {output_dir}")
    
    # Track time and generated figures
    start_time = time.time()
    all_generated = []
    
    # Load data
    print("\n  Loading data...")
    data = load_all_data(base_path)
    print(f"  ✓ Loaded {len(data)} data sources")
    
    # Generate figures
    if args.all or args.figures:
        generated = generate_main_figures(data, output_dir, base_path)
        all_generated.extend(generated)
    
    if args.all or args.extended:
        generated = generate_extended_data(data, output_dir, base_path)
        all_generated.extend(generated)
    
    if args.all or args.preview:
        preview = generate_composite_preview(output_dir)
        if preview:
            all_generated.append('Composite Preview')
    
    # Validate if requested
    if args.validate:
        validate_outputs(output_dir)
    
    # Print summary
    print_summary(start_time, all_generated, output_dir)
    
    # Copy manifest
    _create_manifest(output_dir, all_generated)


def _create_manifest(output_dir: Path, figures: list):
    """Create manifest file listing all generated figures."""
    manifest_path = output_dir / 'MANIFEST.txt'
    
    with open(manifest_path, 'w') as f:
        f.write("FLAMES Manuscript Figures\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("Main Figures:\n")
        for fig in figures:
            if not fig.startswith('ED'):
                f.write(f"  - {fig}\n")
        
        f.write("\nExtended Data Figures:\n")
        for fig in figures:
            if fig.startswith('ED'):
                f.write(f"  - {fig}\n")
        
        f.write("\nFiles:\n")
        for pdf in sorted(output_dir.glob('*.pdf')):
            size_kb = pdf.stat().st_size / 1024
            f.write(f"  {pdf.name} ({size_kb:.1f} KB)\n")
    
    print(f"\n  Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
