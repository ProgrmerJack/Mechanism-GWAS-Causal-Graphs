"""
FLAMES Figure Generation Package
================================

Publication-quality figures for Nature Genetics submission.
Implements all journal guidelines with modern matplotlib best practices.

Main Figures:
    - Figure 1: Calibration Overview (3 panels)
    - Figure 2: Stress Tests & Generalization (2 panels)
    - Figure 3: Case Studies (3 panels)
    - Figure 4: Benchmark Performance (2 panels)

Extended Data Figures:
    - ED Figure 1-10: Various supplementary analyses

Author: FLAMES Project
Version: 2.0.0
"""

from .style import (
    setup_nature_style,
    get_figure_size,
    add_panel_letter,
    validate_figure_dimensions,
    COLORS,
    OKABE_ITO,
    METHOD_COLORS,
    TIER_COLORS,
    SINGLE_COL,
    DOUBLE_COL,
    MAX_HEIGHT,
    PUB_DPI,
)
from .utils import (
    load_all_data,
    save_figure,
    check_overlaps,
    wrap_text,
    adaptive_annotation,
    compute_calibration_curve,
    generate_preview_sheet,
)
from .figure1 import create_figure_1
from .figure2 import create_figure_2
from .figure3 import create_figure_3
from .figure4 import create_figure_4
from .extended_data import (
    create_all_ed_figures,
    create_ed_figure_1,
    create_ed_figure_2,
    create_ed_figure_3,
    create_ed_figure_4,
    create_ed_figure_5,
    create_ed_figure_6,
    create_ed_figure_7,
    create_ed_figure_8,
    create_ed_figure_9,
    create_ed_figure_10,
)

__all__ = [
    # Style
    'setup_nature_style',
    'get_figure_size',
    'add_panel_letter',
    'validate_figure_dimensions',
    'COLORS',
    'OKABE_ITO', 
    'METHOD_COLORS',
    'TIER_COLORS',
    'SINGLE_COL',
    'DOUBLE_COL',
    'MAX_HEIGHT',
    'PUB_DPI',
    # Utils
    'load_all_data',
    'save_figure',
    'check_overlaps',
    'wrap_text',
    'adaptive_annotation',
    'compute_calibration_curve',
    'generate_preview_sheet',
    # Main figures
    'create_figure_1',
    'create_figure_2',
    'create_figure_3',
    'create_figure_4',
    # Extended Data
    'create_all_ed_figures',
    'create_ed_figure_1',
    'create_ed_figure_2',
    'create_ed_figure_3',
    'create_ed_figure_4',
    'create_ed_figure_5',
    'create_ed_figure_6',
    'create_ed_figure_7',
    'create_ed_figure_8',
    'create_ed_figure_9',
    'create_ed_figure_10',
]
