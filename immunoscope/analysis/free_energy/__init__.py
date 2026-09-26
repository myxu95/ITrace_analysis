"""
Free Energy Landscape (FEL) Analysis Module

This module provides tools for calculating and visualizing free energy landscapes
from molecular dynamics collective variables.

Core Components:
    - FELCalculator: Core FEL computation
    - FELVisualizer: 2D/3D visualization
    - FELComparator: Difference FEL analysis
    - FELConvergenceAnalyzer: Convergence assessment
"""

from .fel_calculator import FELCalculator
from .fel_visualizer import FELVisualizer

__all__ = [
    'FELCalculator',
    'FELVisualizer',
]

__version__ = '0.1.0'
