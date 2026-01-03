#!/usr/bin/env python3
"""
FLAMES Figure Utilities Module
==============================

Helper functions for figure generation:
- Data loading
- Figure saving with metadata
- Overlap detection
- Text wrapping and annotation

Author: FLAMES Project
"""

import json
import os
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from datetime import datetime

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from .style import PUB_DPI, COLORS


# =============================================================================
# DATA LOADING
# =============================================================================

def load_all_data(base_path: Path) -> Dict[str, Any]:
    """
    Load all data files required for FLAMES figures.
    
    Parameters
    ----------
    base_path : Path
        Root directory of the project (contains 'results/' folder)
    
    Returns
    -------
    Dict[str, Any]
        Dictionary with loaded data keyed by type
    """
    data = {}
    results_path = base_path / 'results'
    
    print("\n" + "=" * 60)
    print("Loading Data Files")
    print("=" * 60)
    
    # 1. Leave-family-out stress test results
    lfo_path = results_path / 'stress_test' / 'leave_family_out_results.json'
    if lfo_path.exists():
        with open(lfo_path, 'r') as f:
            data['lfo'] = json.load(f)
        n_families = len(data['lfo'].get('results', []))
        print(f"✓ Leave-family-out stress test: {n_families} disease families")
    else:
        print(f"✗ Missing: {lfo_path.relative_to(base_path)}")
    
    # 2. Cross-validation ECE results
    cv_path = results_path / 'calibration_validation' / 'cv_ece_results.json'
    if cv_path.exists():
        with open(cv_path, 'r') as f:
            data['cv_ece'] = json.load(f)
        ece = data['cv_ece'].get('cv_mean_ece', 0)
        print(f"✓ Cross-validation ECE: mean={ece:.4f}")
    else:
        print(f"✗ Missing: {cv_path.relative_to(base_path)}")
    
    # 3. Expected discoveries (calibration validation)
    exp_path = results_path / 'decision_curve' / 'expected_discoveries.json'
    if exp_path.exists():
        with open(exp_path, 'r') as f:
            data['expected'] = json.load(f)
        n_k = len(data['expected'])
        print(f"✓ Expected discoveries: {n_k} k-values")
    else:
        print(f"✗ Missing: {exp_path.relative_to(base_path)}")
    
    # 4. Benchmark comparison metrics
    bench_path = results_path / 'baselines' / 'post2021_comparison_metrics.tsv'
    if bench_path.exists():
        data['benchmark'] = pd.read_csv(bench_path, sep='\t')
        n_methods = len(data['benchmark'])
        print(f"✓ Benchmark comparison: {n_methods} methods")
    else:
        print(f"✗ Missing: {bench_path.relative_to(base_path)}")
    
    # 5. Case studies
    case_path = results_path / 'case_studies' / 'case_studies_detailed.json'
    if case_path.exists():
        with open(case_path, 'r') as f:
            data['cases'] = json.load(f)
        n_cases = len(data['cases'])
        print(f"✓ Case studies: {n_cases} examples")
    else:
        print(f"✗ Missing: {case_path.relative_to(base_path)}")
    
    # 6. Ablation study results
    ablation_path = results_path / 'ablation'
    if ablation_path.exists():
        ablation_files = list(ablation_path.glob('*.json'))
        if ablation_files:
            data['ablation'] = {}
            for f in ablation_files:
                with open(f, 'r') as fp:
                    data['ablation'][f.stem] = json.load(fp)
            print(f"✓ Ablation study: {len(ablation_files)} files")
    
    # 7. Disease-specific calibration
    disease_cal_path = results_path / 'calibration_validation' / 'disease_calibration.tsv'
    if disease_cal_path.exists():
        data['disease_calibration'] = pd.read_csv(disease_cal_path, sep='\t')
        print(f"✓ Disease calibration: {len(data['disease_calibration'])} diseases")
    
    print("=" * 60 + "\n")
    return data


# =============================================================================
# FIGURE SAVING
# =============================================================================

def save_figure(
    fig,
    output_dir: Path,
    name: str,
    title: Optional[str] = None,
    author: str = "FLAMES Project",
    subject: Optional[str] = None,
    formats: Tuple[str, ...] = ('pdf', 'png'),
    copy_to_manuscript: bool = True,
    base_path: Optional[Path] = None,
) -> List[Path]:
    """
    Save figure in multiple formats with proper metadata.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save
    output_dir : Path
        Output directory
    name : str
        Base filename (without extension)
    title : str, optional
        PDF metadata title
    author : str
        PDF metadata author
    subject : str, optional
        PDF metadata subject
    formats : Tuple[str, ...]
        Output formats ('pdf', 'png', 'svg', 'tiff')
    copy_to_manuscript : bool
        Also save to manuscript/figures/
    base_path : Path, optional
        Project root for manuscript directory
    
    Returns
    -------
    List[Path]
        Paths to saved files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_files = []
    
    for fmt in formats:
        output_path = output_dir / f"{name}.{fmt}"
        
        if fmt == 'pdf':
            # Save with metadata using PdfPages
            with PdfPages(output_path) as pdf:
                pdf.savefig(fig, bbox_inches='tight', pad_inches=0.02)
                
                # Add metadata
                d = pdf.infodict()
                d['Title'] = title or name.replace('_', ' ').title()
                d['Author'] = author
                d['Subject'] = subject or 'FLAMES Nature Genetics Figure'
                d['Creator'] = 'FLAMES Figure Generator v2.0'
                d['CreationDate'] = datetime.now()
                d['ModDate'] = datetime.now()
        
        elif fmt == 'png':
            fig.savefig(
                output_path,
                format='png',
                dpi=PUB_DPI,
                bbox_inches='tight',
                pad_inches=0.02,
                facecolor='white',
                edgecolor='none',
            )
        
        elif fmt == 'tiff':
            fig.savefig(
                output_path,
                format='tiff',
                dpi=PUB_DPI,
                bbox_inches='tight',
                pad_inches=0.02,
                facecolor='white',
                pil_kwargs={'compression': 'tiff_lzw'},
            )
        
        elif fmt == 'svg':
            fig.savefig(
                output_path,
                format='svg',
                bbox_inches='tight',
                pad_inches=0.02,
            )
        
        saved_files.append(output_path)
        print(f"  → Saved {output_path.name}")
    
    # Copy to manuscript directory if requested
    if copy_to_manuscript and base_path is not None:
        manuscript_dir = base_path / 'manuscript' / 'figures'
        if manuscript_dir.exists():
            import shutil
            for src in saved_files:
                if src.suffix in ('.pdf', '.png'):
                    dst = manuscript_dir / src.name
                    shutil.copy(src, dst)
            print(f"  → Copied to manuscript/figures/")
    
    return saved_files


# =============================================================================
# OVERLAP DETECTION
# =============================================================================

def check_overlaps(
    fig,
    renderer=None,
    verbose: bool = True,
) -> List[Tuple[str, str]]:
    """
    Check for overlapping text elements in a figure.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to check
    renderer : matplotlib.backend_bases.RendererBase, optional
        Renderer for getting extents
    verbose : bool
        Print warnings for detected overlaps
    
    Returns
    -------
    List[Tuple[str, str]]
        List of (text1, text2) pairs that overlap
    """
    if renderer is None:
        # Create a dummy renderer
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    
    overlaps = []
    
    # Collect all text objects
    texts = []
    for ax in fig.axes:
        texts.extend(ax.texts)
        if ax.title.get_text():
            texts.append(ax.title)
        if ax.xaxis.label.get_text():
            texts.append(ax.xaxis.label)
        if ax.yaxis.label.get_text():
            texts.append(ax.yaxis.label)
        # Tick labels
        texts.extend(ax.get_xticklabels())
        texts.extend(ax.get_yticklabels())
    
    # Add figure-level texts
    texts.extend(fig.texts)
    
    # Check pairwise overlaps
    for i, t1 in enumerate(texts):
        if not t1.get_text():
            continue
        try:
            bbox1 = t1.get_window_extent(renderer)
        except Exception:
            continue
        
        for t2 in texts[i+1:]:
            if not t2.get_text():
                continue
            try:
                bbox2 = t2.get_window_extent(renderer)
            except Exception:
                continue
            
            # Check intersection
            if bbox1.overlaps(bbox2):
                overlaps.append((t1.get_text()[:20], t2.get_text()[:20]))
    
    if verbose and overlaps:
        print(f"  ⚠ Detected {len(overlaps)} text overlaps:")
        for t1, t2 in overlaps[:5]:  # Show first 5
            print(f"    '{t1}' ↔ '{t2}'")
    
    return overlaps


# =============================================================================
# TEXT UTILITIES
# =============================================================================

def wrap_text(
    text: str,
    max_width: int = 20,
    max_lines: int = 3,
    ellipsis: bool = True,
) -> str:
    """
    Wrap long text for figure labels.
    
    Parameters
    ----------
    text : str
        Text to wrap
    max_width : int
        Maximum characters per line
    max_lines : int
        Maximum number of lines
    ellipsis : bool
        Add '...' if truncated
    
    Returns
    -------
    str
        Wrapped text with newlines
    """
    lines = textwrap.wrap(text, width=max_width)
    
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        if ellipsis:
            lines[-1] = lines[-1][:max_width-3] + '...'
    
    return '\n'.join(lines)


def adaptive_annotation(
    ax,
    x: float,
    y: float,
    text: str,
    existing_positions: List[Tuple[float, float]] = None,
    offset_step: float = 0.05,
    max_attempts: int = 8,
    fontsize: int = 6,
    **kwargs,
) -> Tuple[float, float]:
    """
    Add annotation with adaptive positioning to avoid overlaps.
    
    Tries different offsets to find a non-overlapping position.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to annotate
    x, y : float
        Data coordinates of point to annotate
    text : str
        Annotation text
    existing_positions : List[Tuple[float, float]], optional
        List of already-used annotation positions
    offset_step : float
        Step size for offset attempts (in axes fraction)
    max_attempts : int
        Maximum positioning attempts
    fontsize : int
        Text font size
    **kwargs
        Additional arguments for ax.annotate
    
    Returns
    -------
    Tuple[float, float]
        Final (x_offset, y_offset) used
    """
    if existing_positions is None:
        existing_positions = []
    
    # Direction offsets to try (in axes coordinates)
    offsets = [
        (1, 1), (1, -1), (-1, 1), (-1, -1),
        (0, 1), (0, -1), (1, 0), (-1, 0),
    ]
    
    # Transform data to axes coordinates
    transform = ax.transData + ax.transAxes.inverted()
    ax_x, ax_y = transform.transform((x, y))
    
    # Try different offsets
    for attempt in range(max_attempts):
        # Calculate offset
        dir_idx = attempt % len(offsets)
        scale = 1 + (attempt // len(offsets))
        
        dx = offsets[dir_idx][0] * offset_step * scale
        dy = offsets[dir_idx][1] * offset_step * scale
        
        new_pos = (ax_x + dx, ax_y + dy)
        
        # Check for conflicts with existing positions
        conflict = False
        for ex, ey in existing_positions:
            if abs(new_pos[0] - ex) < 0.08 and abs(new_pos[1] - ey) < 0.05:
                conflict = True
                break
        
        if not conflict:
            # Found good position
            ax.annotate(
                text,
                xy=(x, y),
                xytext=(dx * 100, dy * 100),  # Points offset
                textcoords='offset points',
                fontsize=fontsize,
                ha='center',
                va='center',
                arrowprops=dict(
                    arrowstyle='-',
                    color=COLORS['gray'],
                    linewidth=0.5,
                ) if abs(dx) + abs(dy) > offset_step else None,
                **kwargs,
            )
            existing_positions.append(new_pos)
            return (dx, dy)
    
    # Fallback: just place it
    ax.annotate(
        text,
        xy=(x, y),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=fontsize,
        **kwargs,
    )
    return (0.05, 0.05)


# =============================================================================
# CALIBRATION UTILITIES
# =============================================================================

def compute_calibration_curve(
    predictions: np.ndarray,
    outcomes: np.ndarray,
    n_bins: int = 10,
    strategy: str = 'uniform',
) -> Dict[str, np.ndarray]:
    """
    Compute calibration curve data with confidence intervals.
    
    Parameters
    ----------
    predictions : np.ndarray
        Predicted probabilities
    outcomes : np.ndarray
        Binary outcomes (0 or 1)
    n_bins : int
        Number of calibration bins
    strategy : str
        'uniform' for equal-width bins, 'quantile' for equal-count bins
    
    Returns
    -------
    Dict with keys:
        'bin_centers': center of each bin
        'bin_means': observed frequency in each bin
        'bin_counts': number of samples per bin
        'bin_ci_lower': lower 95% CI
        'bin_ci_upper': upper 95% CI
        'ece': expected calibration error
    """
    predictions = np.asarray(predictions)
    outcomes = np.asarray(outcomes)
    
    # Determine bin edges
    if strategy == 'quantile':
        # Equal-count bins
        quantiles = np.linspace(0, 100, n_bins + 1)
        bin_edges = np.percentile(predictions, quantiles)
        bin_edges = np.unique(bin_edges)  # Remove duplicates
    else:
        # Equal-width bins
        bin_edges = np.linspace(0, 1, n_bins + 1)
    
    n_bins_actual = len(bin_edges) - 1
    
    bin_centers = []
    bin_means = []
    bin_counts = []
    bin_ci_lower = []
    bin_ci_upper = []
    
    for i in range(n_bins_actual):
        if i < n_bins_actual - 1:
            mask = (predictions >= bin_edges[i]) & (predictions < bin_edges[i+1])
        else:
            mask = (predictions >= bin_edges[i]) & (predictions <= bin_edges[i+1])
        
        count = mask.sum()
        bin_counts.append(count)
        
        if count >= 5:
            center = predictions[mask].mean()
            mean = outcomes[mask].mean()
            
            # Wilson score interval for confidence
            z = 1.96
            n = count
            p = mean
            
            denom = 1 + z**2/n
            center_adj = (p + z**2/(2*n)) / denom
            margin = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
            
            bin_centers.append(center)
            bin_means.append(mean)
            bin_ci_lower.append(max(0, center_adj - margin))
            bin_ci_upper.append(min(1, center_adj + margin))
        else:
            bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
            bin_means.append(np.nan)
            bin_ci_lower.append(np.nan)
            bin_ci_upper.append(np.nan)
    
    # Compute ECE
    valid = ~np.isnan(bin_means)
    bin_centers = np.array(bin_centers)
    bin_means = np.array(bin_means)
    bin_counts = np.array(bin_counts)
    
    if valid.any():
        ece = np.sum(
            bin_counts[valid] * np.abs(bin_means[valid] - bin_centers[valid])
        ) / bin_counts[valid].sum()
    else:
        ece = np.nan
    
    return {
        'bin_centers': bin_centers,
        'bin_means': bin_means,
        'bin_counts': bin_counts,
        'bin_ci_lower': np.array(bin_ci_lower),
        'bin_ci_upper': np.array(bin_ci_upper),
        'ece': ece,
    }


def generate_preview_sheet(
    figure_paths: List[Path],
    output_path: Path,
    cols: int = 3,
    title: str = "FLAMES Figures Preview",
) -> None:
    """
    Generate a composite preview sheet of all figures.
    
    Parameters
    ----------
    figure_paths : List[Path]
        Paths to figure PNG files
    output_path : Path
        Output path for preview sheet
    cols : int
        Number of columns in grid
    title : str
        Title for preview sheet
    """
    from PIL import Image
    
    # Load images
    images = []
    names = []
    for p in figure_paths:
        if p.suffix == '.png' and p.exists():
            images.append(Image.open(p))
            names.append(p.stem)
    
    if not images:
        print("  No images found for preview sheet")
        return
    
    # Calculate grid
    n = len(images)
    rows = (n + cols - 1) // cols
    
    # Get max dimensions
    max_width = max(img.width for img in images)
    max_height = max(img.height for img in images)
    
    # Scale down if too large
    scale = min(1.0, 2400 / max_width, 1600 / max_height)
    thumb_width = int(max_width * scale)
    thumb_height = int(max_height * scale)
    
    # Create canvas
    margin = 40
    label_height = 30
    canvas_width = cols * thumb_width + (cols + 1) * margin
    canvas_height = rows * (thumb_height + label_height) + (rows + 1) * margin + 60
    
    canvas = Image.new('RGB', (canvas_width, canvas_height), 'white')
    
    # Paste images
    for idx, (img, name) in enumerate(zip(images, names)):
        row = idx // cols
        col = idx % cols
        
        x = margin + col * (thumb_width + margin)
        y = 60 + margin + row * (thumb_height + label_height + margin)
        
        # Resize and paste
        img_resized = img.resize(
            (int(img.width * scale), int(img.height * scale)),
            Image.Resampling.LANCZOS
        )
        canvas.paste(img_resized, (x, y))
    
    # Save
    canvas.save(output_path, 'PNG', dpi=(150, 150))
    print(f"  → Preview sheet saved: {output_path.name}")


if __name__ == '__main__':
    # Test data loading
    base = Path(__file__).parent.parent.parent
    data = load_all_data(base)
    print(f"Loaded {len(data)} data types")
