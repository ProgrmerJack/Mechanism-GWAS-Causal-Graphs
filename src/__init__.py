"""
Mechanism-GWAS-Causal-Graphs

A framework for constructing calibrated mechanism graphs linking 
variants → regulatory elements → genes → tissues → traits.
"""

__version__ = "1.0.0"
__author__ = "Research Team"
__email__ = "contact@example.com"

from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# Submodule imports
from . import harmonization
from . import finemapping
from . import colocalization
from . import mechanism_graphs
from . import calibration
from . import utils

__all__ = [
    "harmonization",
    "finemapping", 
    "colocalization",
    "mechanism_graphs",
    "calibration",
    "utils",
    "PROJECT_ROOT",
    "CONFIG_DIR",
    "DATA_DIR",
    "RESULTS_DIR",
]
