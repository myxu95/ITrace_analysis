"""
Preprocess Pipeline Node - PBC processing wrapper.
"""

from typing import Optional
import logging

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError
from ....analysis.trajectory.pbc import PBCProcessor

logger = logging.getLogger(__name__)


class PreprocessNode(PipelineNode):
    """
    PBC preprocessing pipeline node.

    This node wraps the PBCProcessor and performs periodic boundary
    condition (PBC) artifact removal. It:
    - Centers trajectory on shortest chain
    - Makes molecules whole
    - Fits rotational and translational motion

    The processed trajectory is saved to context.trajectory_processed.

    Example:
        >>> node = PreprocessNode(method="3step")
        >>> context = node.execute(context)
        >>> print(context.trajectory_processed)
        './results/1ao7/md_processed.xtc'
    """

    def __init__(self,
                 method: str = "3step",
                 dt: Optional[float] = None,
                 gmx_executable: str = "gmx",
                 fit_group: Optional[str] = None,
                 output_group: Optional[str] = None,
                 center_group: Optional[str] = None,
                 keep_temp_files: bool = False,
                 name: Optional[str] = None):
        """
        Initialize preprocess node.

        Args:
            method: PBC processing method ('2step' or '3step')
            dt: Time interval for frame sampling in ps (optional)
            gmx_executable: GROMACS executable command
            fit_group: Optional fit group override
            output_group: Optional output group override
            center_group: Optional center group override
            keep_temp_files: Whether to keep temporary files
            name: Node name (defaults to 'PreprocessNode')
        """
        super().__init__(name=name)
        self.method = method
        self.dt = dt
        self.gmx_executable = gmx_executable
        self.fit_group = fit_group
        self.output_group = output_group
        self.center_group = center_group
        self.keep_temp_files = keep_temp_files

    def validate_inputs(self, context: PipelineContext):
        """
        Validate that required inputs exist in context.

        Args:
            context: Pipeline context

        Raises:
            PipelineError: If required files are missing
        """
        if not context.trajectory_raw:
            raise PipelineError(
                node_name=self.name,
                reason="trajectory_raw not found in context",
                context_state={'system_id': context.system_id}
            )

        if not context.topology:
            raise PipelineError(
                node_name=self.name,
                reason="topology not found in context",
                context_state={'system_id': context.system_id}
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute PBC processing.

        Args:
            context: Pipeline context

        Returns:
            Updated context with trajectory_processed

        Raises:
            PipelineError: If execution fails
        """
        self.logger.info(f"Executing {self.name} for system {context.system_id}")

        # Validate inputs
        try:
            self.validate_inputs(context)
        except PipelineError as e:
            context.add_error(str(e))
            context.should_stop = True
            return context

        # Create output directory
        output_dir = context.get_output_path("")  # Get directory
        output_file = context.get_output_path("md_processed.xtc")

        try:
            # Initialize PBC processor
            pbc_processor = PBCProcessor(
                gmx_executable=self.gmx_executable,
                keep_temp_files=self.keep_temp_files,
            )

            # Run PBC processing
            result = pbc_processor.comprehensive_pbc_process(
                trajectory=context.trajectory_raw,
                topology=context.topology,
                output_dir=output_dir,
                method=self.method,
                dt=self.dt,
                fit_group=self.fit_group,
                output_group=self.output_group,
                center_group=self.center_group,
            )

            # Update context
            if result.get('success', False):
                context.trajectory_processed = result.get('processed', output_file)

                self.logger.info(
                    f"PBC processing completed. "
                    f"Output: {context.trajectory_processed}"
                )

                # Add metadata
                context.metadata['pbc_processing'] = {
                    'method': self.method,
                    'dt': self.dt,
                    'fit_group': self.fit_group,
                    'output_group': self.output_group,
                    'center_group': self.center_group,
                    'shortest_chain': result.get('shortest_chain_id')
                }
            else:
                error_msg = f"PBC processing failed: {result.get('error', 'unknown error')}"
                context.add_error(error_msg)
                context.should_stop = True
                self.logger.error(error_msg)

        except Exception as e:
            error_msg = f"PBC processing exception: {str(e)}"
            context.add_error(error_msg)
            context.should_stop = True
            self.logger.exception(error_msg)

        return context

    def __repr__(self) -> str:
        return f"PreprocessNode(method='{self.method}')"
