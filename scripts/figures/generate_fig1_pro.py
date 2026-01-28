#!/usr/bin/env python3
"""
FIGURE 1: Mechanism Graph Framework Overview
Nature Genetics Professional Quality

Caption requirements:
a) Mechanism graph for SORT1 locus showing causal path from rs12740374 through liver enhancer
b) Five-stage inference pipeline diagram
c) Formal probabilistic model showing path multiplication and noisy-OR
d) Comparison with L2G showing opaque score vs mechanism decomposition
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, ConnectionPatch
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path
import json

# ============================================================================
# NATURE GENETICS SPECIFICATIONS
# ============================================================================
SINGLE_COL_MM = 89
DOUBLE_COL_MM = 183
MM_TO_INCH = 1 / 25.4
DPI = 600

# Okabe-Ito colorblind-safe palette
COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00', 
    'green': '#009E73',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'skyblue': '#56B4E9',
    'yellow': '#F0E442',
    'black': '#000000',
    'gray': '#999999',
    'light_gray': '#E5E5E5',
    'white': '#FFFFFF'
}

def setup_nature_style():
    """Configure matplotlib for Nature Genetics specifications."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'axes.titlesize': 8,
        'axes.labelsize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,
        'figure.dpi': 150,
        'savefig.dpi': DPI,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def draw_mechanism_graph(ax):
    """Panel a: SORT1 locus mechanism graph."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Title
    ax.text(0.02, 0.98, 'a', transform=ax.transAxes, fontsize=10, 
            fontweight='bold', va='top', ha='left')
    
    # Node positions
    nodes = {
        'variant': (2, 6.5),
        'enhancer': (5, 6.5),
        'gene': (8, 6.5),
        'tissue': (5, 4.5),
    }
    
    # Draw variant node (diamond shape)
    variant_box = FancyBboxPatch((nodes['variant'][0]-0.8, nodes['variant'][1]-0.4), 
                                  1.6, 0.8, boxstyle="round,pad=0.05",
                                  facecolor=COLORS['orange'], edgecolor='black',
                                  linewidth=1.2, alpha=0.9)
    ax.add_patch(variant_box)
    ax.text(nodes['variant'][0], nodes['variant'][1], 'rs12740374\nPIP = 0.94',
            ha='center', va='center', fontsize=6, fontweight='bold')
    
    # Draw enhancer node (hexagon-like)
    enhancer_box = FancyBboxPatch((nodes['enhancer'][0]-1.0, nodes['enhancer'][1]-0.4),
                                   2.0, 0.8, boxstyle="round,pad=0.05",
                                   facecolor=COLORS['green'], edgecolor='black',
                                   linewidth=1.2, alpha=0.9)
    ax.add_patch(enhancer_box)
    ax.text(nodes['enhancer'][0], nodes['enhancer'][1], 'Liver Enhancer\nABC = 0.84',
            ha='center', va='center', fontsize=6, fontweight='bold', color='white')
    
    # Draw gene node
    gene_box = FancyBboxPatch((nodes['gene'][0]-0.7, nodes['gene'][1]-0.4),
                               1.4, 0.8, boxstyle="round,pad=0.05",
                               facecolor=COLORS['blue'], edgecolor='black',
                               linewidth=1.2, alpha=0.9)
    ax.add_patch(gene_box)
    ax.text(nodes['gene'][0], nodes['gene'][1], 'SORT1\nP = 0.79',
            ha='center', va='center', fontsize=6, fontweight='bold', color='white')
    
    # Draw tissue context
    tissue_box = FancyBboxPatch((nodes['tissue'][0]-1.2, nodes['tissue'][1]-0.35),
                                 2.4, 0.7, boxstyle="round,pad=0.05",
                                 facecolor=COLORS['light_gray'], edgecolor=COLORS['gray'],
                                 linewidth=1, linestyle='--', alpha=0.7)
    ax.add_patch(tissue_box)
    ax.text(nodes['tissue'][0], nodes['tissue'][1], 'Hepatocyte context',
            ha='center', va='center', fontsize=6, style='italic', color=COLORS['gray'])
    
    # Draw connecting lines with probability labels
    # Variant -> Enhancer (clean connector line with chevron)
    v_end = nodes['variant'][0]+0.8
    e_start = nodes['enhancer'][0]-1.0
    ax.plot([v_end, e_start], [nodes['variant'][1], nodes['enhancer'][1]], 
            color=COLORS['black'], lw=2, solid_capstyle='round')
    ax.text((v_end + e_start)/2, nodes['variant'][1], '▸', ha='center', va='center', 
            fontsize=12, color=COLORS['black'], fontweight='bold')
    ax.text(3.5, 7.0, 'cCRE\noverlap', ha='center', va='bottom', fontsize=5)
    
    # Enhancer -> Gene (clean connector line with chevron)
    e_end = nodes['enhancer'][0]+1.0
    g_start = nodes['gene'][0]-0.7
    ax.plot([e_end, g_start], [nodes['enhancer'][1], nodes['gene'][1]], 
            color=COLORS['black'], lw=2, solid_capstyle='round')
    ax.text((e_end + g_start)/2, nodes['enhancer'][1], '▸', ha='center', va='center', 
            fontsize=12, color=COLORS['black'], fontweight='bold')
    ax.text(6.5, 7.0, 'coloc\nPP.H4=0.96', ha='center', va='bottom', fontsize=5)
    
    # Tissue connection (dashed)
    ax.plot([nodes['tissue'][0], nodes['enhancer'][0]], 
            [nodes['tissue'][1]+0.35, nodes['enhancer'][1]-0.4],
            color=COLORS['gray'], linestyle='--', linewidth=1)
    
    # Path probability annotation
    ax.text(5, 2.5, 'Path Probability: P = 0.94 × 0.84 × 0.96 = 0.79',
            ha='center', va='center', fontsize=7, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=COLORS['skyblue'], 
                     edgecolor=COLORS['blue'], alpha=0.3))
    
    # Locus label
    ax.text(5, 1.2, '1p13.3 (LDL-C associated locus)',
            ha='center', va='center', fontsize=6, style='italic')

def draw_pipeline(ax):
    """Panel b: Five-stage inference pipeline."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    ax.text(0.02, 0.98, 'b', transform=ax.transAxes, fontsize=10,
            fontweight='bold', va='top', ha='left')
    
    stages = [
        ('SuSiE\nFine-mapping', COLORS['orange'], 'Variant\nPIP'),
        ('ENCODE\ncCRE', COLORS['yellow'], 'Regulatory\nelement'),
        ('ABC/PCHi-C\nensemble', COLORS['green'], 'E-G\nlinking'),
        ('coloc.susie', COLORS['blue'], 'Tissue\ncoloc'),
        ('Noisy-OR', COLORS['purple'], 'Gene\nprob'),
    ]
    
    y_center = 5.5
    box_width = 1.4
    box_height = 1.2
    gap = 0.35
    
    total_width = len(stages) * box_width + (len(stages) - 1) * gap
    start_x = (10 - total_width) / 2
    
    for i, (name, color, output) in enumerate(stages):
        x = start_x + i * (box_width + gap)
        
        # Main box
        box = FancyBboxPatch((x, y_center - box_height/2), box_width, box_height,
                              boxstyle="round,pad=0.08", facecolor=color,
                              edgecolor='black', linewidth=1, alpha=0.85)
        ax.add_patch(box)
        
        # Determine text color based on background
        text_color = 'white' if color in [COLORS['green'], COLORS['blue'], COLORS['purple']] else 'black'
        ax.text(x + box_width/2, y_center, name, ha='center', va='center',
                fontsize=5.5, fontweight='bold', color=text_color)
        
        # Output label below
        ax.text(x + box_width/2, y_center - box_height/2 - 0.4, output,
                ha='center', va='top', fontsize=5, color=COLORS['gray'])
        
        # Connector to next stage (line with chevron)
        if i < len(stages) - 1:
            ax.plot([x + box_width + 0.05, x + box_width + gap - 0.05], [y_center, y_center],
                   color='black', lw=1, solid_capstyle='round')
            ax.text(x + box_width + gap/2, y_center, '▸', ha='center', va='center',
                   fontsize=8, color='black', fontweight='bold')
    
    # Data source labels at top
    data_sources = ['GWAS\nsummary', 'ENCODE\nV3', 'ABC/\nPCHi-C', 'eQTL\nCatalogue', '']
    for i, source in enumerate(data_sources[:-1]):
        x = start_x + i * (box_width + gap) + box_width/2
        ax.text(x, y_center + box_height/2 + 0.6, source, ha='center', va='bottom',
                fontsize=4.5, color=COLORS['gray'], style='italic')
        # Clean vertical connector line with downward chevron
        ax.plot([x, x], [y_center + box_height/2 + 0.5, y_center + box_height/2 + 0.15],
               color=COLORS['gray'], lw=0.5, solid_capstyle='round')
        ax.text(x, y_center + box_height/2 + 0.08, '▾', ha='center', va='center',
               fontsize=6, color=COLORS['gray'])
    
    # Final output
    ax.text(5, 2.8, 'Calibrated gene-level probabilities\nwith mechanism decomposition',
            ha='center', va='center', fontsize=6,
            bbox=dict(boxstyle='round,pad=0.3', facecolor=COLORS['light_gray'],
                     edgecolor=COLORS['gray'], alpha=0.5))

def draw_probabilistic_model(ax):
    """Panel c: Formal probabilistic model."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    ax.text(0.02, 0.98, 'c', transform=ax.transAxes, fontsize=10,
            fontweight='bold', va='top', ha='left')
    
    # Path multiplication formula
    ax.text(5, 6.5, 'Path Probability (single path):', ha='center', va='center',
            fontsize=7, fontweight='bold')
    ax.text(5, 5.5, r'$P_{path} = P_{variant} \times P_{E \rightarrow G} \times P_{coloc}$',
            ha='center', va='center', fontsize=9)
    
    # Noisy-OR formula
    ax.text(5, 4.0, 'Gene Probability (multiple paths via Noisy-OR):', ha='center', va='center',
            fontsize=7, fontweight='bold')
    ax.text(5, 3.0, r'$P_{gene} = 1 - \prod_{i}(1 - P_{path_i})$',
            ha='center', va='center', fontsize=9)
    
    # Example calculation
    ax.text(5, 1.5, 'Example: SORT1 = 1 - (1 - 0.79) = 0.79 (single dominant path)',
            ha='center', va='center', fontsize=6, style='italic',
            bbox=dict(boxstyle='round,pad=0.2', facecolor=COLORS['light_gray'],
                     edgecolor=COLORS['gray'], alpha=0.5))

def draw_l2g_comparison(ax):
    """Panel d: Comparison with L2G."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    ax.text(0.02, 0.98, 'd', transform=ax.transAxes, fontsize=10,
            fontweight='bold', va='top', ha='left')
    
    # L2G side (left)
    ax.text(2.5, 7.2, 'L2G (baseline)', ha='center', va='center',
            fontsize=7, fontweight='bold', color=COLORS['vermillion'])
    
    # L2G black box
    l2g_box = FancyBboxPatch((1.2, 4.5), 2.6, 2.0, boxstyle="round,pad=0.1",
                              facecolor=COLORS['gray'], edgecolor='black',
                              linewidth=1.5, alpha=0.8)
    ax.add_patch(l2g_box)
    ax.text(2.5, 5.5, 'Opaque\nML Model', ha='center', va='center',
            fontsize=6, fontweight='bold', color='white')
    
    # L2G output
    ax.text(2.5, 3.5, 'SORT1 = 0.82', ha='center', va='center', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.2', facecolor=COLORS['vermillion'],
                     edgecolor='black', alpha=0.3))
    ax.text(2.5, 2.5, '- No mechanism\n- No calibration\n- No uncertainty',
            ha='center', va='center', fontsize=5, color=COLORS['vermillion'])
    
    # Divider
    ax.axvline(x=5, ymin=0.15, ymax=0.85, color=COLORS['gray'], linestyle='--', linewidth=1)
    ax.text(5, 7.5, 'vs', ha='center', va='center', fontsize=8, fontweight='bold')
    
    # Mechanism Graph side (right)
    ax.text(7.5, 7.2, 'Mechanism Graphs', ha='center', va='center',
            fontsize=7, fontweight='bold', color=COLORS['blue'])
    
    # Transparent pipeline
    mg_box = FancyBboxPatch((6.2, 4.5), 2.6, 2.0, boxstyle="round,pad=0.1",
                             facecolor=COLORS['skyblue'], edgecolor=COLORS['blue'],
                             linewidth=1.5, alpha=0.4)
    ax.add_patch(mg_box)
    ax.text(7.5, 5.5, 'Interpretable\nCausal Path', ha='center', va='center',
            fontsize=6, fontweight='bold', color=COLORS['blue'])
    
    # Mechanism graph output
    ax.text(7.5, 3.5, 'SORT1 = 0.79 [0.71-0.86]', ha='center', va='center', fontsize=6,
            bbox=dict(boxstyle='round,pad=0.2', facecolor=COLORS['blue'],
                     edgecolor='black', alpha=0.3))
    ax.text(7.5, 2.5, '+ Full mechanism\n+ Calibrated probability\n+ 95% CI quantified',
            ha='center', va='center', fontsize=5, color=COLORS['green'])

def create_figure_1():
    """Generate Figure 1: Mechanism Graph Framework Overview."""
    setup_nature_style()
    
    # Double column width, appropriate height
    fig_width = DOUBLE_COL_MM * MM_TO_INCH
    fig_height = fig_width * 0.65  # Golden ratio inspired
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Create 2x2 grid
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.2,
                          left=0.05, right=0.98, top=0.95, bottom=0.05)
    
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    
    # Draw each panel
    draw_mechanism_graph(ax_a)
    draw_pipeline(ax_b)
    draw_probabilistic_model(ax_c)
    draw_l2g_comparison(ax_d)
    
    return fig

def main():
    """Generate and save Figure 1."""
    base_path = Path(__file__).parent.parent
    output_dir = base_path / 'figures_professional'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("GENERATING FIGURE 1: Mechanism Graph Framework Overview")
    print("=" * 70)
    
    fig = create_figure_1()
    
    # Save in multiple formats
    for fmt, dpi_val in [('pdf', DPI), ('png', DPI), ('tiff', DPI)]:
        output_path = output_dir / f'fig1_calibration_overview.{fmt}'
        fig.savefig(output_path, format=fmt, dpi=dpi_val, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        print(f"  ✓ Saved {output_path.name}")
    
    # Also copy to manuscript/figures
    manuscript_dir = base_path / 'manuscript' / 'figures'
    manuscript_dir.mkdir(exist_ok=True)
    fig.savefig(manuscript_dir / 'fig1_calibration_overview.pdf',
               format='pdf', dpi=DPI, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print(f"  ✓ Copied to manuscript/figures/")
    
    plt.close(fig)
    
    print("\n" + "=" * 70)
    print("FIGURE 1 COMPLETE")
    print("=" * 70)

if __name__ == '__main__':
    main()
