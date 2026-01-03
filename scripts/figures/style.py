#!/usr/bin/env python3
"""
FLAMES Figure Style Module
==========================

Nature Genetics publication standards for all figures.
Implements Okabe-Ito colorblind-safe palette and precise journal dimensions.

Standards Reference:
- Nature Genetics Author Guidelines (2024)
- Figure dimensions: 89mm (single), 183mm (double), 170mm max height
- Font: Arial/Helvetica, 5-8pt for labels, 7pt max body text
- Line widths: 0.25-1.0pt
- DPI: 600 for print

Author: FLAMES Project
"""

from typing import Tuple, Optional
import textwrap

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np

# =============================================================================
# DIMENSION CONSTANTS (Nature Genetics Specifications)
# =============================================================================

MM_TO_INCH = 1 / 25.4

# Column widths
SINGLE_COL_MM = 89   # Single column width
DOUBLE_COL_MM = 183  # Double column width (full page width)
MAX_HEIGHT_MM = 170  # Maximum figure height

# Convert to inches for matplotlib
SINGLE_COL = SINGLE_COL_MM * MM_TO_INCH   # ~3.5 inches
DOUBLE_COL = DOUBLE_COL_MM * MM_TO_INCH   # ~7.2 inches
MAX_HEIGHT = MAX_HEIGHT_MM * MM_TO_INCH   # ~6.7 inches

# Resolution
PUB_DPI = 600  # Required for Nature print quality

# =============================================================================
# OKABE-ITO COLORBLIND-SAFE PALETTE
# =============================================================================
# This palette provides 8 high-contrast colors distinguishable by
# individuals with all types of color vision deficiency.
# Reference: Okabe & Ito (2008) Color Universal Design

OKABE_ITO = {
    'orange':       '#E69F00',  # Primary warm color
    'sky_blue':     '#56B4E9',  # Primary cool color
    'bluish_green': '#009E73',  # Teal/green
    'yellow':       '#F0E442',  # Bright accent
    'blue':         '#0072B2',  # Primary blue
    'vermillion':   '#D55E00',  # Red-orange
    'reddish_purple': '#CC79A7',  # Pink/purple
    'black':        '#000000',  # Neutral
}

# Extended palette for consistent usage
COLORS = {
    # Primary Okabe-Ito colors
    'blue': OKABE_ITO['blue'],
    'orange': OKABE_ITO['orange'],
    'green': OKABE_ITO['bluish_green'],
    'vermillion': OKABE_ITO['vermillion'],
    'purple': OKABE_ITO['reddish_purple'],
    'skyblue': OKABE_ITO['sky_blue'],
    'yellow': OKABE_ITO['yellow'],
    'black': OKABE_ITO['black'],
    # Neutrals for backgrounds/annotations
    'gray': '#666666',
    'lightgray': '#E5E5E5',
    'white': '#FFFFFF',
    'darkgray': '#333333',
    'dark': '#222222',  # Near-black for text/lines
}

# Method-specific colors for benchmark plots
METHOD_COLORS = {
    'FLAMES':               COLORS['blue'],
    'Distance':             COLORS['orange'],
    'L2G':                  COLORS['green'],
    'cS2G':                 COLORS['purple'],
    'cS2G_LocusAware_max':  COLORS['purple'],
    'ABC_Only':             COLORS['skyblue'],
    'eQTL_Only':            COLORS['vermillion'],
    'PoPS':                 COLORS['yellow'],
}

# Validation tier colors
TIER_COLORS = {
    # Simple tier keys (used in figure3.py)
    'tier1':           COLORS['green'],
    'tier2':           COLORS['orange'],
    'tier3':           COLORS['yellow'],
    # Detailed tier keys (used elsewhere)
    'Tier1_CRISPR':    COLORS['green'],
    'Tier1_Drug':      COLORS['blue'],
    'Tier1_Mendelian': COLORS['blue'],
    'Tier2_MultiEvidence': COLORS['orange'],
    'Unknown':         COLORS['gray'],
}

# =============================================================================
# STYLE CONFIGURATION
# =============================================================================

def setup_nature_style() -> None:
    """
    Configure matplotlib for Nature Genetics publication standards.
    
    This function sets all rcParams to meet journal requirements:
    - TrueType font embedding (pdf.fonttype = 42)
    - Arial/Helvetica fonts
    - Appropriate font sizes (5-8pt)
    - Line widths meeting minimum requirements (0.25pt)
    - Clean, minimal aesthetics
    
    Call this at the start of any figure generation script.
    """
    # Reset to defaults first
    plt.style.use('default')
    
    plt.rcParams.update({
        # =================================================================
        # CRITICAL: Font embedding for journal submission
        # =================================================================
        'pdf.fonttype': 42,    # TrueType fonts (REQUIRED by Nature)
        'ps.fonttype': 42,     # TrueType for PostScript
        
        # =================================================================
        # Figure settings
        # =================================================================
        'figure.dpi': 150,               # Screen display DPI
        'figure.facecolor': 'white',     # White background
        'figure.edgecolor': 'white',
        'figure.constrained_layout.use': True,  # Auto-prevent overlaps
        'figure.constrained_layout.h_pad': 0.04167,  # 3pt padding
        'figure.constrained_layout.w_pad': 0.04167,
        
        # =================================================================
        # Save settings
        # =================================================================
        'savefig.dpi': PUB_DPI,
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'none',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'savefig.transparent': False,
        
        # =================================================================
        # Font settings (Arial/Helvetica, Nature standard)
        # =================================================================
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica Neue', 'Helvetica', 
                           'DejaVu Sans', 'Liberation Sans'],
        'font.size': 7,               # Default 7pt (Nature max)
        'font.weight': 'normal',
        
        # =================================================================
        # Axes settings
        # =================================================================
        'axes.labelsize': 7,          # Axis label size
        'axes.titlesize': 8,          # Subplot title size
        'axes.titleweight': 'bold',   # Bold titles
        'axes.titlepad': 4,           # Title padding
        'axes.labelpad': 3,           # Label padding from axis
        'axes.labelweight': 'normal',
        'axes.linewidth': 0.5,        # Spine width (min 0.25pt)
        'axes.spines.top': False,     # Remove top spine
        'axes.spines.right': False,   # Remove right spine
        'axes.unicode_minus': False,  # Use proper minus sign
        'axes.axisbelow': True,       # Grid below data
        'axes.facecolor': 'white',    # White plot background
        'axes.edgecolor': COLORS['darkgray'],
        'axes.prop_cycle': plt.cycler('color', list(OKABE_ITO.values())),
        
        # =================================================================
        # Tick settings
        # =================================================================
        'xtick.labelsize': 6,         # X tick labels
        'ytick.labelsize': 6,         # Y tick labels
        'xtick.major.width': 0.5,     # Major tick width
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,     # Minor tick width
        'ytick.minor.width': 0.3,
        'xtick.major.size': 3,        # Major tick length
        'ytick.major.size': 3,
        'xtick.minor.size': 1.5,      # Minor tick length
        'ytick.minor.size': 1.5,
        'xtick.major.pad': 2,         # Tick-to-label padding
        'ytick.major.pad': 2,
        'xtick.direction': 'out',     # Ticks point outward
        'ytick.direction': 'out',
        'xtick.color': COLORS['darkgray'],
        'ytick.color': COLORS['darkgray'],
        
        # =================================================================
        # Legend settings (clean, no box)
        # =================================================================
        'legend.fontsize': 6,
        'legend.frameon': False,      # No frame (Nature style)
        'legend.framealpha': 0,
        'legend.edgecolor': 'none',
        'legend.fancybox': False,
        'legend.borderpad': 0.4,
        'legend.labelspacing': 0.3,
        'legend.handlelength': 1.2,
        'legend.handletextpad': 0.4,
        'legend.columnspacing': 1.0,
        'legend.loc': 'best',
        
        # =================================================================
        # Line and marker settings
        # =================================================================
        'lines.linewidth': 0.75,      # Default line width
        'lines.markersize': 4,        # Default marker size
        'lines.markeredgewidth': 0.5,
        'lines.solid_capstyle': 'round',
        'lines.solid_joinstyle': 'round',
        
        # =================================================================
        # Patch settings (bars, etc.)
        # =================================================================
        'patch.linewidth': 0.5,
        'patch.edgecolor': 'white',
        'patch.facecolor': COLORS['blue'],
        
        # =================================================================
        # Grid settings (subtle if used)
        # =================================================================
        'grid.linewidth': 0.25,       # Minimum visible
        'grid.alpha': 0.3,
        'grid.linestyle': '-',
        'grid.color': COLORS['lightgray'],
        
        # =================================================================
        # Error bar settings
        # =================================================================
        'errorbar.capsize': 2,
        
        # =================================================================
        # Hatch settings
        # =================================================================
        'hatch.linewidth': 0.5,
        'hatch.color': COLORS['darkgray'],
    })


def get_figure_size(
    width: str = 'double',
    aspect_ratio: float = 0.618,  # Golden ratio
    height: Optional[float] = None,
    n_rows: int = 1,
    n_cols: int = 1,
    row_height: Optional[float] = None,
) -> Tuple[float, float]:
    """
    Calculate figure dimensions meeting Nature Genetics requirements.
    
    Parameters
    ----------
    width : str
        'single' (~3.5") or 'double' (~7.2")
    aspect_ratio : float
        Height/width ratio (default golden ratio)
    height : float, optional
        Explicit height in inches (overrides aspect_ratio)
    n_rows : int
        Number of subplot rows
    n_cols : int
        Number of subplot columns  
    row_height : float, optional
        Height per row in inches (for dynamic sizing)
    
    Returns
    -------
    Tuple[float, float]
        (width, height) in inches
    
    Examples
    --------
    >>> get_figure_size('double')  # Standard double-column
    (7.2, 4.45)
    >>> get_figure_size('single', height=4)  # Single col, fixed height
    (3.5, 4.0)
    """
    # Determine width
    fig_width = DOUBLE_COL if width == 'double' else SINGLE_COL
    
    # Determine height
    if height is not None:
        fig_height = height
    elif row_height is not None:
        fig_height = row_height * n_rows
    else:
        fig_height = fig_width * aspect_ratio
    
    # Enforce maximum height
    fig_height = min(fig_height, MAX_HEIGHT)
    
    return (fig_width, fig_height)


def add_panel_letter(
    ax,
    letter: str,
    loc: str = 'upper left',
    fontsize: int = 8,
    fontweight: str = 'bold',
    offset: Tuple[float, float] = (-0.12, 1.08),  # Moved further outside to avoid tick labels
) -> None:
    """
    Add a panel letter (a, b, c, etc.) to an axes.
    
    Nature Genetics style uses lowercase bold 8pt letters.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to label
    letter : str
        The panel letter (e.g., 'a', 'b', 'c')
    loc : str
        Location ('upper left', 'upper right', etc.)
    fontsize : int
        Font size in points (Nature standard: 8pt)
    fontweight : str
        Font weight ('bold' for Nature)
    offset : Tuple[float, float]
        (x, y) offset from corner in axes coordinates
    """
    # Position letter using axes transform
    ax.text(
        offset[0], offset[1],
        letter.lower(),  # Nature uses lowercase
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight=fontweight,
        verticalalignment='top',
        horizontalalignment='left',
        family='sans-serif',
    )


def create_panel_letter_anchored(
    letter: str,
    fontsize: int = 8,
    loc: str = 'upper left',
) -> AnchoredText:
    """
    Create an AnchoredText object for panel letters.
    
    Alternative to add_panel_letter that uses matplotlib's
    AnchoredText for more reliable positioning.
    
    Parameters
    ----------
    letter : str
        Panel letter to display
    fontsize : int
        Font size (8pt for Nature)
    loc : str
        Location string for AnchoredText
    
    Returns
    -------
    AnchoredText
        Object to add to axes via ax.add_artist()
    """
    return AnchoredText(
        letter.lower(),
        loc=loc,
        frameon=False,
        pad=0.3,
        prop=dict(
            size=fontsize,
            weight='bold',
            family='sans-serif',
        ),
    )


# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

def validate_figure_dimensions(fig) -> dict:
    """
    Check that a figure meets Nature Genetics dimension requirements.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to validate
    
    Returns
    -------
    dict
        Validation results with 'passed', 'width', 'height', 'issues'
    """
    width, height = fig.get_size_inches()
    issues = []
    
    # Check width
    if width > DOUBLE_COL + 0.01:
        issues.append(f"Width {width:.2f}in exceeds max {DOUBLE_COL:.2f}in")
    
    # Check height
    if height > MAX_HEIGHT + 0.01:
        issues.append(f"Height {height:.2f}in exceeds max {MAX_HEIGHT:.2f}in")
    
    # Check if near standard widths
    near_single = abs(width - SINGLE_COL) < 0.1
    near_double = abs(width - DOUBLE_COL) < 0.1
    if not (near_single or near_double):
        issues.append(f"Width {width:.2f}in not near standard ({SINGLE_COL:.2f} or {DOUBLE_COL:.2f})")
    
    return {
        'passed': len(issues) == 0,
        'width_inches': width,
        'height_inches': height,
        'width_mm': width / MM_TO_INCH,
        'height_mm': height / MM_TO_INCH,
        'issues': issues,
    }


def print_style_summary() -> None:
    """Print current style configuration summary."""
    print("\n" + "=" * 60)
    print("FLAMES Figure Style Configuration")
    print("=" * 60)
    print(f"\nDimensions:")
    print(f"  Single column: {SINGLE_COL_MM}mm ({SINGLE_COL:.2f}in)")
    print(f"  Double column: {DOUBLE_COL_MM}mm ({DOUBLE_COL:.2f}in)")
    print(f"  Max height: {MAX_HEIGHT_MM}mm ({MAX_HEIGHT:.2f}in)")
    print(f"\nResolution: {PUB_DPI} DPI")
    print(f"\nFont: {plt.rcParams['font.sans-serif'][0]}")
    print(f"  Default size: {plt.rcParams['font.size']}pt")
    print(f"  Label size: {plt.rcParams['axes.labelsize']}pt")
    print(f"  Tick size: {plt.rcParams['xtick.labelsize']}pt")
    print(f"\nPDF fonttype: {plt.rcParams['pdf.fonttype']} (TrueType)")
    print(f"\nOkabe-Ito palette: {len(OKABE_ITO)} colors")
    print("=" * 60 + "\n")


if __name__ == '__main__':
    setup_nature_style()
    print_style_summary()
