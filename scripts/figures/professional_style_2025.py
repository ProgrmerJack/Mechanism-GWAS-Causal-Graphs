#!/usr/bin/env python3
"""
FLAMES Professional Visualization Style Module (2025 Edition)
=============================================================

State-of-the-art scientific visualization for Nature Genetics submission.

Implements cutting-edge 2025 visualization best practices:
- Edward Tufte's data-ink ratio principles (minimize chartjunk)
- Gradient fills with elegant transparency for confidence intervals
- Sophisticated color harmonics with depth perception
- Clean, minimal aesthetics that emphasize the data
- Raincloud-inspired distribution visualizations
- Professional shadows and subtle glow effects
- Modern typography hierarchy

Research sources:
- Tufte, E.R. "The Visual Display of Quantitative Information"
- Nature Genetics Author Guidelines (2024-2025)
- Allen et al. (2021) "Raincloud plots" - cited 1676 times
- SciencePlots matplotlib styling
- "Error Bars Considered Harmful" - gradient uncertainty encoding

Author: FLAMES Project - Professional Edition
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Shadow, PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
import numpy as np
from typing import Tuple, List, Optional, Union, Dict
from dataclasses import dataclass

# =============================================================================
# PROFESSIONAL 2025 COLOR SYSTEM
# =============================================================================
# Enhanced Okabe-Ito with gradient-ready variants and depth perception

@dataclass
class ColorGradient:
    """Professional gradient with light/main/dark variants for depth."""
    light: str    # Highlight/glow
    main: str     # Primary color
    dark: str     # Shadow/depth
    alpha: float = 0.9  # Default alpha

# Professional gradient palettes
PROFESSIONAL_GRADIENTS = {
    'blue': ColorGradient(
        light='#7CB9E8',   # Light azure
        main='#0072B2',    # Okabe-Ito blue
        dark='#004C78'     # Deep navy
    ),
    'orange': ColorGradient(
        light='#FFD699',   # Light amber
        main='#E69F00',    # Okabe-Ito orange
        dark='#B37700'     # Deep gold
    ),
    'green': ColorGradient(
        light='#66D4B3',   # Light mint
        main='#009E73',    # Okabe-Ito green
        dark='#006B4F'     # Deep forest
    ),
    'vermillion': ColorGradient(
        light='#FF9B80',   # Light coral
        main='#D55E00',    # Okabe-Ito vermillion
        dark='#8B3D00'     # Deep rust
    ),
    'purple': ColorGradient(
        light='#E6B8D9',   # Light lavender
        main='#CC79A7',    # Okabe-Ito purple
        dark='#8A4B72'     # Deep plum
    ),
    'skyblue': ColorGradient(
        light='#B8E0F7',   # Light sky
        main='#56B4E9',    # Okabe-Ito sky blue
        dark='#2E8BC0'     # Deep cerulean
    ),
    'neutral': ColorGradient(
        light='#E8E8E8',   # Light gray
        main='#6B7280',    # Medium gray
        dark='#374151'     # Dark charcoal
    ),
}

# Semantic color assignments (using gradient system)
SEMANTIC_COLORS = {
    'primary': PROFESSIONAL_GRADIENTS['blue'],
    'secondary': PROFESSIONAL_GRADIENTS['orange'],
    'success': PROFESSIONAL_GRADIENTS['green'],
    'danger': PROFESSIONAL_GRADIENTS['vermillion'],
    'accent': PROFESSIONAL_GRADIENTS['purple'],
    'info': PROFESSIONAL_GRADIENTS['skyblue'],
    'neutral': PROFESSIONAL_GRADIENTS['neutral'],
}

# =============================================================================
# NATURE GENETICS PROFESSIONAL DIMENSIONS
# =============================================================================

MM_TO_INCH = 1 / 25.4

# Journal specifications
SINGLE_COL_MM = 89   # Single column width
DOUBLE_COL_MM = 183  # Double column width
MAX_HEIGHT_MM = 170  # Maximum figure height

# Convert to inches
SINGLE_COL = SINGLE_COL_MM * MM_TO_INCH
DOUBLE_COL = DOUBLE_COL_MM * MM_TO_INCH
MAX_HEIGHT = MAX_HEIGHT_MM * MM_TO_INCH

# Professional DPI settings
SCREEN_DPI = 150   # Preview quality
PRINT_DPI = 600    # Publication quality

# =============================================================================
# PROFESSIONAL TYPOGRAPHY HIERARCHY
# =============================================================================

TYPOGRAPHY = {
    'title': {
        'fontsize': 10,
        'fontweight': 'bold',
        'color': '#1a1a1a',
        'family': 'Arial',
    },
    'subtitle': {
        'fontsize': 8,
        'fontweight': 'semibold',
        'color': '#333333',
        'family': 'Arial',
    },
    'body': {
        'fontsize': 7,
        'fontweight': 'normal',
        'color': '#4a4a4a',
        'family': 'Arial',
    },
    'label': {
        'fontsize': 7,
        'fontweight': 'normal',
        'color': '#666666',
        'family': 'Arial',
    },
    'caption': {
        'fontsize': 6,
        'fontweight': 'normal',
        'color': '#888888',
        'family': 'Arial',
    },
    'annotation': {
        'fontsize': 6,
        'fontweight': 'normal',
        'color': '#666666',
        'style': 'italic',
        'family': 'Arial',
    },
    'mono': {
        'fontsize': 6.5,
        'fontweight': 'normal',
        'color': '#333333',
        'family': 'DejaVu Sans Mono',
    },
}

# =============================================================================
# PROFESSIONAL RCPARAMS CONFIGURATION
# =============================================================================

def setup_professional_style() -> None:
    """
    Configure matplotlib for professional 2025 publication standards.
    
    Implements Edward Tufte's principles:
    - Maximize data-ink ratio
    - Minimize chartjunk
    - Clean, elegant aesthetics
    - Proportional and honest representation
    """
    plt.style.use('default')
    
    plt.rcParams.update({
        # =====================================================================
        # FONT EMBEDDING (CRITICAL FOR JOURNAL SUBMISSION)
        # =====================================================================
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        
        # =====================================================================
        # PROFESSIONAL FONT SETTINGS
        # =====================================================================
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
        
        # =====================================================================
        # FIGURE SETTINGS
        # =====================================================================
        'figure.dpi': SCREEN_DPI,
        'figure.facecolor': 'white',
        'figure.edgecolor': 'none',
        'figure.constrained_layout.use': True,
        'figure.constrained_layout.h_pad': 0.04,
        'figure.constrained_layout.w_pad': 0.04,
        'figure.autolayout': False,
        
        # =====================================================================
        # SAVE SETTINGS
        # =====================================================================
        'savefig.dpi': PRINT_DPI,
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'none',
        'savefig.transparent': False,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        
        # =====================================================================
        # AXES (TUFTE-INSPIRED MINIMAL)
        # =====================================================================
        'axes.facecolor': 'white',
        'axes.edgecolor': '#333333',
        'axes.linewidth': 0.6,
        'axes.grid': False,
        'axes.axisbelow': True,
        'axes.labelsize': 7,
        'axes.titlesize': 8,
        'axes.titleweight': 'bold',
        'axes.labelweight': 'normal',
        'axes.labelcolor': '#333333',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.prop_cycle': plt.cycler('color', [
            '#0072B2', '#E69F00', '#009E73', '#D55E00', 
            '#CC79A7', '#56B4E9', '#F0E442', '#666666'
        ]),
        
        # =====================================================================
        # TICKS (MINIMAL AND PRECISE)
        # =====================================================================
        'xtick.major.size': 3,
        'xtick.minor.size': 1.5,
        'xtick.major.width': 0.6,
        'xtick.minor.width': 0.4,
        'xtick.labelsize': 6.5,
        'xtick.color': '#333333',
        'xtick.direction': 'out',
        'xtick.major.pad': 2,
        
        'ytick.major.size': 3,
        'ytick.minor.size': 1.5,
        'ytick.major.width': 0.6,
        'ytick.minor.width': 0.4,
        'ytick.labelsize': 6.5,
        'ytick.color': '#333333',
        'ytick.direction': 'out',
        'ytick.major.pad': 2,
        
        # =====================================================================
        # LINES (PROFESSIONAL THICKNESS)
        # =====================================================================
        'lines.linewidth': 1.2,
        'lines.markersize': 4,
        'lines.markeredgewidth': 0.8,
        'lines.solid_capstyle': 'round',
        'lines.solid_joinstyle': 'round',
        
        # =====================================================================
        # PATCHES (FOR BARS, BOXES)
        # =====================================================================
        'patch.linewidth': 0.6,
        'patch.edgecolor': '#333333',
        'patch.facecolor': '#0072B2',
        
        # =====================================================================
        # LEGEND (CLEAN AND MINIMAL)
        # =====================================================================
        'legend.fontsize': 6.5,
        'legend.frameon': False,
        'legend.borderpad': 0.3,
        'legend.labelspacing': 0.3,
        'legend.handlelength': 1.5,
        'legend.handletextpad': 0.4,
        'legend.columnspacing': 1.0,
        
        # =====================================================================
        # GRID (SUBTLE WHEN USED)
        # =====================================================================
        'grid.color': '#E5E5E5',
        'grid.linewidth': 0.4,
        'grid.alpha': 0.7,
        
        # =====================================================================
        # ERROR BARS
        # =====================================================================
        'errorbar.capsize': 2,
    })


# =============================================================================
# PROFESSIONAL GRADIENT AND SHADING FUNCTIONS
# =============================================================================

def create_gradient_colormap(color_gradient: ColorGradient, n_colors: int = 256) -> mcolors.LinearSegmentedColormap:
    """
    Create a professional gradient colormap from light to dark.
    
    Args:
        color_gradient: ColorGradient object with light/main/dark
        n_colors: Number of color steps
        
    Returns:
        LinearSegmentedColormap for smooth gradients
    """
    colors = [color_gradient.light, color_gradient.main, color_gradient.dark]
    return mcolors.LinearSegmentedColormap.from_list('gradient', colors, N=n_colors)


def add_gradient_fill(ax, x: np.ndarray, y_lower: np.ndarray, y_upper: np.ndarray,
                      color: str = '#0072B2', n_levels: int = 50, 
                      alpha_max: float = 0.4) -> None:
    """
    Add a sophisticated gradient fill for confidence intervals.
    
    Creates a smooth gradient that fades from center (most opaque) to edges
    (most transparent), following 2025 best practices for uncertainty visualization.
    
    Args:
        ax: Matplotlib axes
        x: X-axis values
        y_lower: Lower bound of interval
        y_upper: Upper bound of interval
        color: Base color (will create gradient)
        n_levels: Number of gradient layers (more = smoother)
        alpha_max: Maximum opacity at center
    """
    y_center = (y_lower + y_upper) / 2
    half_width = (y_upper - y_lower) / 2
    
    for i in range(n_levels):
        # Create gradient from center outward
        frac = (i + 1) / n_levels
        alpha = alpha_max * (1 - frac * 0.8)  # Fade toward edges
        
        y_lo = y_center - half_width * frac
        y_hi = y_center + half_width * frac
        
        ax.fill_between(x, y_lo, y_hi, color=color, alpha=alpha, 
                       linewidth=0, zorder=1)


def add_confidence_band(ax, x: np.ndarray, y: np.ndarray, 
                        ci_lower: np.ndarray, ci_upper: np.ndarray,
                        color: str = '#0072B2', line_alpha: float = 1.0,
                        fill_alpha: float = 0.25, label: Optional[str] = None,
                        linewidth: float = 1.5, linestyle: str = '-') -> None:
    """
    Add a line with confidence band using modern gradient fill.
    
    Args:
        ax: Matplotlib axes
        x: X values
        y: Y values (main line)
        ci_lower: Lower confidence interval
        ci_upper: Upper confidence interval
        color: Line and fill color
        line_alpha: Line opacity
        fill_alpha: Fill opacity
        label: Legend label
        linewidth: Line thickness
        linestyle: Line style
    """
    # Gradient confidence band
    add_gradient_fill(ax, x, ci_lower, ci_upper, color=color, 
                      alpha_max=fill_alpha)
    
    # Main line on top
    ax.plot(x, y, color=color, alpha=line_alpha, linewidth=linewidth,
            linestyle=linestyle, label=label, zorder=10)


def create_professional_bar(ax, x: float, height: float, width: float = 0.6,
                           color: Union[str, ColorGradient] = '#0072B2',
                           edge_color: Optional[str] = None,
                           alpha: float = 0.9, add_shadow: bool = True,
                           shadow_offset: Tuple[float, float] = (0.02, -0.02),
                           error: Optional[float] = None,
                           error_color: str = '#333333') -> mpatches.Rectangle:
    """
    Create a professional bar with optional gradient and shadow.
    
    Args:
        ax: Matplotlib axes
        x: X position (center)
        height: Bar height
        width: Bar width
        color: Fill color or ColorGradient
        edge_color: Edge color (None for auto-darken)
        alpha: Fill opacity
        add_shadow: Whether to add subtle shadow
        shadow_offset: Shadow offset (x, y)
        error: Error bar value (optional)
        error_color: Error bar color
        
    Returns:
        Rectangle patch
    """
    # Handle gradient colors
    if isinstance(color, ColorGradient):
        fill_color = color.main
        edge_color = edge_color or color.dark
    else:
        fill_color = color
        if edge_color is None:
            # Auto-darken edge
            edge_color = mcolors.to_hex([c * 0.7 for c in mcolors.to_rgb(color)])
    
    # Shadow first (behind bar)
    if add_shadow:
        shadow = mpatches.Rectangle(
            (x - width/2 + shadow_offset[0], shadow_offset[1]),
            width, height,
            facecolor='#000000', alpha=0.08, linewidth=0, zorder=1
        )
        ax.add_patch(shadow)
    
    # Main bar
    bar = mpatches.Rectangle(
        (x - width/2, 0), width, height,
        facecolor=fill_color, edgecolor=edge_color,
        alpha=alpha, linewidth=0.8, zorder=2
    )
    ax.add_patch(bar)
    
    # Error bar
    if error is not None:
        ax.errorbar(x, height, yerr=error, color=error_color, 
                   capsize=3, capthick=0.8, linewidth=0.8, zorder=3)
    
    return bar


def draw_professional_box(ax, x: float, y: float, width: float, height: float,
                         text: str, color: Union[str, ColorGradient] = '#0072B2',
                         text_color: str = 'white', fontsize: float = 7,
                         fontweight: str = 'bold', alpha: float = 0.9,
                         rounded: bool = True, shadow: bool = True,
                         shadow_offset: Tuple[float, float] = (0.03, -0.03),
                         sublabel: Optional[str] = None) -> FancyBboxPatch:
    """
    Draw a professional labeled box with optional shadow and sublabel.
    
    Args:
        ax: Matplotlib axes
        x, y: Center position
        width, height: Box dimensions
        text: Main label
        color: Fill color or ColorGradient
        text_color: Label color
        fontsize: Label font size
        fontweight: Label font weight
        alpha: Fill opacity
        rounded: Use rounded corners
        shadow: Add drop shadow
        shadow_offset: Shadow offset
        sublabel: Secondary label below main
        
    Returns:
        FancyBboxPatch
    """
    # Handle gradient colors
    if isinstance(color, ColorGradient):
        fill_color = color.main
        edge_color = color.dark
    else:
        fill_color = color
        edge_color = mcolors.to_hex([c * 0.7 for c in mcolors.to_rgb(color)])
    
    boxstyle = 'round,pad=0.03,rounding_size=0.15' if rounded else 'square,pad=0.02'
    
    # Shadow
    if shadow:
        shadow_box = FancyBboxPatch(
            (x - width/2 + shadow_offset[0], y - height/2 + shadow_offset[1]),
            width, height, boxstyle=boxstyle,
            facecolor='#000000', alpha=0.12, linewidth=0, zorder=1
        )
        ax.add_patch(shadow_box)
    
    # Main box
    box = FancyBboxPatch(
        (x - width/2, y - height/2), width, height,
        boxstyle=boxstyle, facecolor=fill_color, edgecolor=edge_color,
        alpha=alpha, linewidth=0.8, zorder=2
    )
    ax.add_patch(box)
    
    # Main label
    if sublabel:
        ax.text(x, y + height * 0.12, text, ha='center', va='center',
               fontsize=fontsize, fontweight=fontweight, color=text_color, zorder=3)
        ax.text(x, y - height * 0.18, sublabel, ha='center', va='center',
               fontsize=fontsize - 1, color=text_color, alpha=0.85, zorder=3)
    else:
        ax.text(x, y, text, ha='center', va='center',
               fontsize=fontsize, fontweight=fontweight, color=text_color, zorder=3)
    
    return box


def draw_connection_line(ax, x1: float, y1: float, x2: float, y2: float,
                        color: str = '#333333', linewidth: float = 1.2,
                        alpha: float = 0.8, style: str = 'solid',
                        show_endpoint: bool = True, 
                        endpoint_size: float = 20,
                        endpoint_marker: str = 'o') -> None:
    """
    Draw a professional connection line with optional endpoint marker.
    
    Args:
        ax: Matplotlib axes
        x1, y1: Start position
        x2, y2: End position
        color: Line color
        linewidth: Line width
        alpha: Line opacity
        style: 'solid', 'dashed', 'dotted'
        show_endpoint: Show marker at endpoint
        endpoint_size: Marker size
        endpoint_marker: Marker style
    """
    linestyles = {
        'solid': '-',
        'dashed': '--',
        'dotted': ':',
    }
    ls = linestyles.get(style, '-')
    
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth,
           linestyle=ls, alpha=alpha, solid_capstyle='round', zorder=5)
    
    if show_endpoint:
        ax.scatter([x2], [y2], c=color, s=endpoint_size, marker=endpoint_marker,
                  alpha=alpha, zorder=6, linewidths=0)


def draw_direction_indicator(ax, x: float, y: float, direction: str = 'right',
                            color: str = '#333333', size: float = 50,
                            alpha: float = 0.9) -> None:
    """
    Draw a direction indicator (triangle/chevron).
    
    Args:
        ax: Matplotlib axes
        x, y: Position
        direction: 'right', 'left', 'up', 'down'
        color: Marker color
        size: Marker size
        alpha: Opacity
    """
    markers = {
        'right': '>',
        'left': '<',
        'up': '^',
        'down': 'v',
    }
    marker = markers.get(direction, '>')
    ax.scatter([x], [y], c=color, s=size, marker=marker, alpha=alpha, 
              zorder=10, linewidths=0)


# =============================================================================
# PROFESSIONAL CALIBRATION DIAGRAM FUNCTIONS
# =============================================================================

def draw_calibration_diagram(ax, predicted: np.ndarray, observed: np.ndarray,
                            bins: int = 10, color: str = '#0072B2',
                            show_histogram: bool = True,
                            show_confidence: bool = True,
                            ci_level: float = 0.95) -> Dict:
    """
    Draw a professional calibration/reliability diagram.
    
    State-of-the-art visualization following latest ML calibration literature.
    
    Args:
        ax: Matplotlib axes
        predicted: Predicted probabilities
        observed: Observed outcomes (0/1)
        bins: Number of calibration bins
        color: Main color
        show_histogram: Show prediction distribution
        show_confidence: Show confidence intervals
        ci_level: Confidence level for intervals
        
    Returns:
        Dict with calibration metrics
    """
    # Calculate calibration curve
    bin_edges = np.linspace(0, 1, bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    bin_means = []
    bin_trues = []
    bin_counts = []
    bin_ci_lower = []
    bin_ci_upper = []
    
    for i in range(bins):
        mask = (predicted >= bin_edges[i]) & (predicted < bin_edges[i+1])
        if mask.sum() > 0:
            mean_pred = predicted[mask].mean()
            true_frac = observed[mask].mean()
            count = mask.sum()
            
            # Wilson score interval for binomial proportion
            from scipy import stats
            z = stats.norm.ppf(1 - (1 - ci_level) / 2)
            n = count
            p = true_frac
            
            denominator = 1 + z**2 / n
            center = (p + z**2 / (2*n)) / denominator
            spread = z * np.sqrt(p * (1 - p) / n + z**2 / (4*n**2)) / denominator
            
            ci_low = max(0, center - spread)
            ci_high = min(1, center + spread)
        else:
            mean_pred = bin_centers[i]
            true_frac = np.nan
            count = 0
            ci_low = np.nan
            ci_high = np.nan
        
        bin_means.append(mean_pred)
        bin_trues.append(true_frac)
        bin_counts.append(count)
        bin_ci_lower.append(ci_low)
        bin_ci_upper.append(ci_high)
    
    bin_means = np.array(bin_means)
    bin_trues = np.array(bin_trues)
    bin_counts = np.array(bin_counts)
    bin_ci_lower = np.array(bin_ci_lower)
    bin_ci_upper = np.array(bin_ci_upper)
    
    # Perfect calibration line
    ax.plot([0, 1], [0, 1], color='#888888', linestyle='--', linewidth=1,
           alpha=0.7, label='Perfect calibration', zorder=1)
    
    # Confidence intervals with gradient
    if show_confidence:
        valid = ~np.isnan(bin_trues)
        if valid.any():
            add_gradient_fill(ax, bin_means[valid], bin_ci_lower[valid], 
                            bin_ci_upper[valid], color=color, alpha_max=0.3)
    
    # Calibration line
    valid = ~np.isnan(bin_trues)
    ax.plot(bin_means[valid], bin_trues[valid], color=color, linewidth=2,
           marker='o', markersize=5, markerfacecolor='white', 
           markeredgecolor=color, markeredgewidth=1.5, label='FLAMES', zorder=10)
    
    # Histogram inset for prediction distribution
    if show_histogram:
        ax_hist = ax.inset_axes([0.55, 0.05, 0.4, 0.15])
        ax_hist.hist(predicted, bins=20, color=color, alpha=0.6, edgecolor='white',
                    linewidth=0.5)
        ax_hist.set_xlim(0, 1)
        ax_hist.set_xticks([0, 0.5, 1])
        ax_hist.set_xticklabels(['0', '0.5', '1'], fontsize=5)
        ax_hist.set_yticks([])
        ax_hist.set_xlabel('Predicted prob.', fontsize=5)
        ax_hist.spines['top'].set_visible(False)
        ax_hist.spines['right'].set_visible(False)
        ax_hist.spines['left'].set_visible(False)
    
    # Labels and formatting
    ax.set_xlabel('Predicted probability', fontsize=7)
    ax.set_ylabel('Observed frequency', fontsize=7)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect('equal')
    ax.legend(loc='upper left', fontsize=6)
    
    # Calculate ECE
    ece = np.nansum(bin_counts * np.abs(bin_trues - bin_means)) / np.sum(bin_counts)
    
    return {
        'ece': ece,
        'bin_means': bin_means,
        'bin_trues': bin_trues,
        'bin_counts': bin_counts,
    }


# =============================================================================
# PROFESSIONAL BAR CHART FUNCTIONS
# =============================================================================

def draw_professional_grouped_bars(ax, categories: List[str], 
                                   groups: Dict[str, List[float]],
                                   errors: Optional[Dict[str, List[float]]] = None,
                                   colors: Optional[Dict[str, str]] = None,
                                   bar_width: float = 0.25,
                                   show_values: bool = True,
                                   value_format: str = '{:.2f}',
                                   ylabel: str = 'Value') -> None:
    """
    Draw professional grouped bar chart.
    
    Args:
        ax: Matplotlib axes
        categories: Category labels (x-axis)
        groups: Dict of group_name -> values
        errors: Optional error bars for each group
        colors: Custom colors per group
        bar_width: Width of individual bars
        show_values: Show value labels on bars
        value_format: Format string for values
        ylabel: Y-axis label
    """
    n_groups = len(groups)
    n_cats = len(categories)
    
    # Default colors
    if colors is None:
        default_colors = ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7']
        colors = {name: default_colors[i % len(default_colors)] 
                 for i, name in enumerate(groups.keys())}
    
    x = np.arange(n_cats)
    offset = (n_groups - 1) * bar_width / 2
    
    for i, (group_name, values) in enumerate(groups.items()):
        pos = x - offset + i * bar_width
        color = colors.get(group_name, '#0072B2')
        err = errors.get(group_name) if errors else None
        
        # Darken edge color
        edge_color = mcolors.to_hex([c * 0.7 for c in mcolors.to_rgb(color)])
        
        bars = ax.bar(pos, values, bar_width * 0.85, label=group_name,
                     color=color, edgecolor=edge_color, linewidth=0.6,
                     alpha=0.9, zorder=2)
        
        if err is not None:
            ax.errorbar(pos, values, yerr=err, fmt='none', color='#333333',
                       capsize=2, capthick=0.6, linewidth=0.6, zorder=3)
        
        if show_values:
            for bar, val in zip(bars, values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, height + 0.01,
                       value_format.format(val), ha='center', va='bottom',
                       fontsize=5, color='#333333')
    
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=6.5)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.legend(fontsize=6, frameon=False, ncol=min(3, n_groups))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# =============================================================================
# PROFESSIONAL ANNOTATIONS
# =============================================================================

def add_comparison_badge(ax, x: float, y: float, value: str, 
                        direction: str = 'up', color: str = '#009E73',
                        fontsize: float = 7) -> None:
    """
    Add a professional comparison badge (e.g., "+18%").
    
    Args:
        ax: Matplotlib axes
        x, y: Position
        value: Text to display (e.g., "+18%")
        direction: 'up' (positive) or 'down' (negative)
        color: Badge color
        fontsize: Text size
    """
    # Background badge
    bbox_props = dict(
        boxstyle='round,pad=0.3,rounding_size=0.2',
        facecolor=color, edgecolor='none', alpha=0.15
    )
    
    # Direction indicator
    prefix = '↑ ' if direction == 'up' else '↓ '
    
    ax.text(x, y, prefix + value, ha='center', va='center',
           fontsize=fontsize, fontweight='bold', color=color,
           bbox=bbox_props, zorder=10)


def add_panel_label(ax, label: str, x: float = -0.12, y: float = 1.08,
                   fontsize: float = 11, fontweight: str = 'bold') -> None:
    """
    Add professional panel label (a, b, c, etc.).
    
    Args:
        ax: Matplotlib axes
        label: Panel label
        x, y: Position in axes coordinates
        fontsize: Label size
        fontweight: Label weight
    """
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
           fontweight=fontweight, va='top', ha='left', color='#1a1a1a')


def add_significance_indicator(ax, x1: float, x2: float, y: float, 
                              height: float = 0.02, 
                              significance: str = '***',
                              color: str = '#333333') -> None:
    """
    Add significance bracket between two groups.
    
    Args:
        ax: Matplotlib axes
        x1, x2: X positions of groups
        y: Y position of bracket base
        height: Bracket height
        significance: Significance text (*, **, ***, ns)
        color: Bracket color
    """
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y],
           color=color, linewidth=0.8, solid_capstyle='round')
    ax.text((x1 + x2) / 2, y + height + 0.005, significance,
           ha='center', va='bottom', fontsize=7, color=color)


# =============================================================================
# HELPER UTILITIES
# =============================================================================

def compute_grid_positions(n: int, start: float, end: float) -> np.ndarray:
    """Compute evenly-spaced positions for n items."""
    return np.linspace(start, end, n)


def lighten_color(color: str, amount: float = 0.3) -> str:
    """Lighten a color by the given amount."""
    rgb = mcolors.to_rgb(color)
    white = np.array([1, 1, 1])
    return mcolors.to_hex(rgb + (white - rgb) * amount)


def darken_color(color: str, amount: float = 0.3) -> str:
    """Darken a color by the given amount."""
    rgb = np.array(mcolors.to_rgb(color))
    return mcolors.to_hex(rgb * (1 - amount))


# =============================================================================
# QUICK STYLE APPLICATION
# =============================================================================

# Apply style on import for convenience
setup_professional_style()

if __name__ == '__main__':
    # Demo figure
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, 3))
    
    # Demo 1: Confidence band
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    ci_lower = y - 0.3
    ci_upper = y + 0.3
    
    add_confidence_band(axes[0], x, y, ci_lower, ci_upper, color='#0072B2')
    axes[0].set_title('Confidence Band', fontsize=8)
    axes[0].set_xlabel('X', fontsize=7)
    axes[0].set_ylabel('Y', fontsize=7)
    
    # Demo 2: Calibration diagram
    np.random.seed(42)
    pred = np.random.beta(2, 5, 500)
    obs = (np.random.random(500) < pred).astype(float)
    draw_calibration_diagram(axes[1], pred, obs, color='#009E73')
    axes[1].set_title('Calibration Diagram', fontsize=8)
    
    # Demo 3: Bar chart
    draw_professional_grouped_bars(
        axes[2],
        categories=['A', 'B', 'C'],
        groups={
            'Method 1': [0.8, 0.7, 0.9],
            'Method 2': [0.6, 0.8, 0.7],
        },
        errors={
            'Method 1': [0.05, 0.06, 0.04],
            'Method 2': [0.07, 0.05, 0.06],
        },
        ylabel='Performance'
    )
    axes[2].set_title('Grouped Bars', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('professional_style_demo.pdf', dpi=600)
    print("Demo saved as professional_style_demo.pdf")
