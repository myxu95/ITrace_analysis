"""
Standard Pipelines - Predefined pipeline templates.

This module provides commonly used pipeline configurations.
"""

from .base_pipeline import Pipeline
from .nodes import PreprocessNode, RMSDNode, RMSDPlotNode, RMSDQualityNode

class PreprocessOnlyPipeline(Pipeline):
    """
    Preprocessing-only pipeline.

    Use this to only perform PBC correction without analysis.

    Example:
        >>> pipeline = PreprocessOnlyPipeline(method="3step")
        >>> result = pipeline.execute(context)
        >>> print(result.trajectory_processed)
    """

    def __init__(self,
                 method: str = "3step",
                 dt: float = None,
                 gmx_executable: str = "gmx",
                 fit_group: str = None,
                 output_group: str = None,
                 center_group: str = None,
                 keep_temp_files: bool = False):
        """
        Initialize preprocess-only pipeline.

        Args:
            method: PBC processing method ('2step' or '3step')
            dt: Time interval for frame sampling in ps
            gmx_executable: GROMACS executable command
            fit_group: Optional fit group override
            output_group: Optional output group override
            center_group: Optional center group override
            keep_temp_files: Whether to keep temporary files
        """
        nodes = [
            PreprocessNode(
                method=method,
                dt=dt,
                gmx_executable=gmx_executable,
                fit_group=fit_group,
                output_group=output_group,
                center_group=center_group,
                keep_temp_files=keep_temp_files,
            ),
        ]
        super().__init__(nodes=nodes)


class PreprocessQualityPipeline(Pipeline):
    """
    Preprocessing plus RMSD quality-assessment pipeline.

    This pipeline follows the mature pbc_rmsd design while staying in the new architecture:
    1. PBC preprocessing
    2. RMSD calculation
    3. RMSD convergence analysis and quality report
    """

    def __init__(self,
                 method: str = "2step",
                 dt: float = None,
                 gmx_executable: str = "gmx",
                 fit_group: str = None,
                 output_group: str = None,
                 center_group: str = None,
                 keep_temp_files: bool = False,
                 rmsd_selection: str = "backbone"):
        nodes = [
            PreprocessNode(
                method=method,
                dt=dt,
                gmx_executable=gmx_executable,
                fit_group=fit_group,
                output_group=output_group,
                center_group=center_group,
                keep_temp_files=keep_temp_files,
            ),
            RMSDNode(
                selection=rmsd_selection,
                output_filename="rmsd.xvg",
                output_subdir="analysis/rmsd",
            ),
            RMSDPlotNode(),
            RMSDQualityNode(),
        ]
        super().__init__(nodes=nodes)
