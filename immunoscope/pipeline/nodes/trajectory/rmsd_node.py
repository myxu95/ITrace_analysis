"""
RMSD Pipeline Node - Wrapper around RMSDCalculator.
"""

from typing import Optional
import logging

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError
from ....analysis.trajectory.rmsd_refactored import RMSDCalculator, RMSDInput

logger = logging.getLogger(__name__)


class RMSDNode(PipelineNode):
    """
    RMSD calculation pipeline node.

    This node wraps the RMSDCalculator and integrates it into the pipeline.
    It is responsible for:
    - Validating that processed trajectory exists in context
    - Calling RMSDCalculator with proper inputs
    - Recording results to context
    - Handling errors gracefully

    Example:
        >>> pipeline = Pipeline(nodes=[
        ...     PreprocessNode(),
        ...     RMSDNode(selection="protein and name CA"),
        ...     ReportNode()
        ... ])
    """

    def __init__(self,
                 selection: str = "protein and name CA",
                 reference_frame: int = 0,
                 method: str = "mdanalysis",
                 output_filename: str = "rmsd.xvg",
                 output_subdir: Optional[str] = None,
                 name: Optional[str] = None):
        """
        Initialize RMSD node.

        Args:
            selection: Atom selection string for RMSD calculation
            reference_frame: Reference frame index
            method: Calculation method ('mdanalysis' or 'gromacs')
            output_filename: RMSD output file name
            output_subdir: Optional output subdirectory
            name: Node name (defaults to 'RMSDNode')
        """
        super().__init__(name=name)
        self.selection = selection
        self.reference_frame = reference_frame
        self.method = method
        self.output_filename = output_filename
        self.output_subdir = output_subdir

    def validate_inputs(self, context: PipelineContext):
        """
        Validate that required inputs exist in context.

        Args:
            context: Pipeline context

        Raises:
            PipelineError: If trajectory_processed is missing
        """
        if not context.trajectory_processed:
            raise PipelineError(
                node_name=self.name,
                reason="trajectory_processed not found in context",
                context_state={
                    'system_id': context.system_id,
                    'has_trajectory_raw': bool(context.trajectory_raw),
                    'has_trajectory_processed': False
                }
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute RMSD calculation.

        Args:
            context: Pipeline context

        Returns:
            Updated context with RMSD results

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

        # Prepare output file path
        if self.output_subdir:
            output_file = context.get_subdir_path(self.output_subdir, self.output_filename)
        else:
            output_file = context.get_output_path(self.output_filename)

        # Create input parameters
        input_params = RMSDInput(
            topology=context.topology,
            trajectory=context.trajectory_processed,
            selection=self.selection,
            reference_frame=self.reference_frame,
            output_file=output_file,
            method=self.method
        )

        # Call RMSDCalculator
        calculator = RMSDCalculator()
        result = calculator.calculate(input_params)

        # Handle result
        if result.success:
            # Record results to context
            context.add_result('rmsd', result.to_dict())

            self.logger.info(
                f"RMSD calculation completed. "
                f"Mean: {result.mean_rmsd:.3f} nm, Std: {result.std_rmsd:.3f} nm"
            )
        else:
            # Record error
            error_msg = f"RMSD calculation failed: {result.error_message}"
            context.add_error(error_msg)
            self.logger.error(error_msg)

            # Optionally stop pipeline
            # context.should_stop = True

        return context

    def __repr__(self) -> str:
        return (f"RMSDNode(selection='{self.selection}', "
                f"method='{self.method}')")
