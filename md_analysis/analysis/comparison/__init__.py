"""Two-system result comparison analysis module."""

from .comparison_schema import ComparisonArtifacts, SingleCaseArtifacts
from .single_case_loader import SingleCaseLoader
from .comparison_builder import (
    ComparisonBuildResult,
    SystemComparisonBuilder,
)
from .comparison_preflight import ComparisonPreflightBuilder, ComparisonPreflightResult
from .comparison_plotter import (
    write_flexibility_comparison_plot,
    write_interaction_family_comparison_plot,
    write_quality_interface_comparison_plot,
    write_rrcs_comparison_plot,
)
from .statistical_tests import StatisticalComparator, StatisticalComparisonResult
from .fel_comparison import FELComparator
from .fel_visualizer import FELVisualizer

__all__ = [
    "ComparisonArtifacts",
    "SingleCaseArtifacts",
    "SingleCaseLoader",
    "ComparisonBuildResult",
    "ComparisonPreflightBuilder",
    "ComparisonPreflightResult",
    "SystemComparisonBuilder",
    "write_flexibility_comparison_plot",
    "write_interaction_family_comparison_plot",
    "write_quality_interface_comparison_plot",
    "write_rrcs_comparison_plot",
    "StatisticalComparator",
    "StatisticalComparisonResult",
    "FELComparator",
    "FELComparisonResult",
    "FELVisualizer",
]
