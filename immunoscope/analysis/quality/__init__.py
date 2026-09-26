"""
ImmunoScope Quality Analysis Module

This module provides comprehensive quality control and validation tools for MD simulation data.
Includes completeness checking, structure validation, batch tracking, quality reporting,
energy-based quality analysis, and trajectory convergence assessment.
"""

# Import modules without scipy dependency
from .md_completeness import MDCompletenessChecker
from .structure_validator import StructureValidator
from .batch_tracker import BatchTracker
from .quality_reporter import QualityReporter
from .post_pbc_validator import PostPBCValidator

# Optional imports for modules with scipy dependency
try:
    from .energy_quality import EnergyQualityChecker
    from .convergence_checker import TrajectoryConvergenceChecker
    HAS_SCIPY_MODULES = True
except ImportError:
    EnergyQualityChecker = None
    TrajectoryConvergenceChecker = None
    HAS_SCIPY_MODULES = False

__all__ = [
    'MDCompletenessChecker',
    'StructureValidator',
    'BatchTracker',
    'QualityReporter',
    'PostPBCValidator',
    'EnergyQualityChecker',
    'TrajectoryConvergenceChecker'
]