"""
Docking Angle Analysis Pipeline Node

Pipeline node for calculating TCR-pMHC docking angles, integrating
the standardized DockingAngleAnalyzer into the ITrace pipeline framework.

Author: ITrace Development Team
Date: 2026-03-18
"""

from typing import Optional
from pathlib import Path

from md_analysis.core.base_node import PipelineNode
from md_analysis.core.context import PipelineContext
from md_analysis.core.exceptions import PipelineError
from md_analysis.analysis.angles import (
    DockingAngleAnalyzer,
    DockingAngleInput,
    DockingAngleResult
)


class DockingAngleNode(PipelineNode):
    """
    Pipeline node for docking angle analysis.

    This node calculates TCR-pMHC docking angles (Crossing and Incident)
    using the standardized DockingAngleAnalyzer. It integrates seamlessly
    with the ITrace pipeline framework.

    Responsibilities:
    1. Validate inputs from context
    2. Build DockingAngleInput from context
    3. Call DockingAngleAnalyzer
    4. Record results to context
    5. Handle errors gracefully

    Parameters
    ----------
    stride : int, default=1
        Frame stride for trajectory analysis
    output_subdir : str, default="angles"
        Subdirectory name for output files
    auto_identify_chains : bool, default=True
        Whether to use automatic chain identification
        If False, requires selections in context.selections

    Examples
    --------
    Basic usage in a pipeline:
    >>> from md_analysis.pipeline import Pipeline
    >>> from md_analysis.core import PipelineContext
    >>> from md_analysis.pipeline.nodes import DockingAngleNode
    >>>
    >>> pipeline = Pipeline([
    ...     PBCNode(),
    ...     DockingAngleNode(stride=10),
    ...     ReportNode()
    ... ])
    >>>
    >>> context = PipelineContext(
    ...     system_id="1ao7",
    ...     topology="md.tpr",
    ...     trajectory_raw="md.xtc"
    ... )
    >>> result = pipeline.execute(context)
    >>> print(result.results['docking_angles']['statistics'])

    With manual chain selections:
    >>> context.selections['mhc'] = 'chainID A'
    >>> context.selections['tcr_alpha'] = 'chainID D'
    >>> context.selections['tcr_beta'] = 'chainID E'
    >>> node = DockingAngleNode(auto_identify_chains=False)
    """

    def __init__(
        self,
        stride: int = 1,
        output_subdir: str = "angles",
        auto_identify_chains: bool = True,
        print_each_frame: bool = False,
        name: Optional[str] = None
    ):
        """
        Initialize docking angle node.

        Parameters
        ----------
        stride : int, default=1
            Frame stride for trajectory analysis
        output_subdir : str, default="angles"
            Subdirectory for output files
        auto_identify_chains : bool, default=True
            Whether to use automatic chain identification
        print_each_frame : bool, default=False
            Whether to log per-frame angle results
        name : str, optional
            Node name (defaults to class name)
        """
        super().__init__(name=name)
        self.stride = stride
        self.output_subdir = output_subdir
        self.auto_identify_chains = auto_identify_chains
        self.print_each_frame = print_each_frame
        self.analyzer = DockingAngleAnalyzer()

    def validate_inputs(self, context: PipelineContext):
        """
        Validate required inputs from context.

        Parameters
        ----------
        context : PipelineContext
            Pipeline context

        Raises
        ------
        PipelineError
            If required inputs are missing
        """
        # Prefer processed trajectories for angle analysis, but allow raw trajectories.
        has_processed = bool(context.metadata.get('trajectory_processed') or context.trajectory_processed)
        has_raw = bool(context.trajectory_raw)
        if not has_processed and not has_raw:
            raise PipelineError(
                node_name=self.name,
                reason="Missing trajectory input. Provide processed or raw trajectory first."
            )

        # If manual mode, check for chain selections
        if not self.auto_identify_chains:
            required_selections = ['mhc', 'tcr_alpha', 'tcr_beta']
            missing = [s for s in required_selections if s not in context.selections]

            if missing:
                raise PipelineError(
                    node_name=self.name,
                    reason=f"Manual chain mode requires selections: {missing}"
                )

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute docking angle analysis.

        Parameters
        ----------
        context : PipelineContext
            Pipeline context

        Returns
        -------
        PipelineContext
            Updated context with angle results
        """
        try:
            self.logger.info(f"Starting docking angle analysis for {context.system_id}")

            # 1. Validate inputs
            self.validate_inputs(context)

            # 2. Determine trajectory path
            trajectory_path = (
                context.metadata.get('trajectory_processed') or
                context.trajectory_processed
            )

            if not trajectory_path:
                # Use raw trajectory if no processed version available
                self.logger.warning("No processed trajectory found, using raw trajectory")
                trajectory_path = context.trajectory_raw

            # 3. Build input parameters
            input_params = self._build_input_params(context, trajectory_path)

            # 4. Execute analysis
            self.logger.info(f"  Stride: {self.stride}")
            self.logger.info(f"  Auto-identify chains: {self.auto_identify_chains}")
            self.logger.info(f"  Print each frame: {self.print_each_frame}")

            result = self.analyzer.analyze(input_params)

            # 5. Handle results
            if result.success:
                self._record_results(context, result)
                self.logger.info("Docking angle analysis completed successfully")
            else:
                error_msg = f"Angle analysis failed: {result.error_message}"
                context.add_error(error_msg)
                context.should_stop = True
                self.logger.error(error_msg)

        except PipelineError:
            # Re-raise pipeline errors
            raise

        except Exception as e:
            error_msg = f"Unexpected error in {self.name}: {str(e)}"
            self.logger.exception(error_msg)
            context.add_error(error_msg)
            context.should_stop = True

        return context

    def _build_input_params(
        self,
        context: PipelineContext,
        trajectory_path: str
    ) -> DockingAngleInput:
        """
        Build DockingAngleInput from context.

        Parameters
        ----------
        context : PipelineContext
            Pipeline context
        trajectory_path : str
            Path to trajectory file

        Returns
        -------
        DockingAngleInput
            Standardized input parameters
        """
        # Trajectory analysis must prefer topology matched to XTC; single-frame analysis may fall back to structure files.
        topology_path = context.topology if trajectory_path else (context.structure_pdb or context.topology)

        chain_mapping = context.metadata.get('chain_mapping') or {}
        mhc_chain = chain_mapping.get('mhc_alpha')
        tcr_alpha_chain = chain_mapping.get('tcr_alpha')
        tcr_beta_chain = chain_mapping.get('tcr_beta')
        use_manual_selections = bool(mhc_chain and tcr_alpha_chain and tcr_beta_chain)

        # Output directory
        output_dir = context.get_subdir_path(self.output_subdir, "")

        # Build input
        input_params = DockingAngleInput(
            topology=topology_path,
            trajectory=trajectory_path,
            stride=self.stride,
            print_each_frame=self.print_each_frame,
            output_dir=output_dir,
            auto_identify_chains=(False if use_manual_selections else self.auto_identify_chains)
        )

        # Prefer upstream chain identification results to avoid duplicate detection inside the angle analyzer.
        if use_manual_selections:
            input_params.mhc_selection = f'chainID {mhc_chain}'
            input_params.tcr_alpha_selection = f'chainID {tcr_alpha_chain}'
            input_params.tcr_beta_selection = f'chainID {tcr_beta_chain}'
        elif not self.auto_identify_chains:
            input_params.mhc_selection = context.selections.get('mhc')
            input_params.tcr_alpha_selection = context.selections.get('tcr_alpha')
            input_params.tcr_beta_selection = context.selections.get('tcr_beta')

        return input_params

    def _record_results(
        self,
        context: PipelineContext,
        result: DockingAngleResult
    ):
        """
        Record analysis results to context.

        Parameters
        ----------
        context : PipelineContext
            Pipeline context
        result : DockingAngleResult
            Analysis result
        """
        # Add to context.results
        context.results['docking_angles'] = {
            'crossing_angle': result.crossing_angle,
            'incident_angle': result.incident_angle,
            'statistics': result.statistics,
            'output_files': result.output_files,
            'metadata': result.metadata,
        }

        # Add metadata
        if 'angle_analysis' not in context.metadata:
            context.metadata['angle_analysis'] = {}

        context.metadata['angle_analysis'].update(result.metadata)

        # Log summary
        if result.statistics:
            stats = result.statistics
            self.logger.info(
                f"  Crossing: {stats['crossing_mean']:.2f} ± {stats['crossing_std']:.2f}° "
                f"[{stats['crossing_min']:.2f}, {stats['crossing_max']:.2f}]"
            )
            self.logger.info(
                f"  Incident: {stats['incident_mean']:.2f} ± {stats['incident_std']:.2f}° "
                f"[{stats['incident_min']:.2f}, {stats['incident_max']:.2f}]"
            )
        elif result.crossing_angle is not None:
            self.logger.info(
                f"  Crossing: {result.crossing_angle:.2f}°, "
                f"Incident: {result.incident_angle:.2f}°"
            )
