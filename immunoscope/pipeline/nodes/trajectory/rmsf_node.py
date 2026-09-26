"""
RMSF Pipeline Node - Wrapper around GROMACS rmsf calculation.

This node calculates Root Mean Square Fluctuation for specified atom groups.
"""

from typing import Optional, List
import logging
import subprocess
from pathlib import Path

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError

logger = logging.getLogger(__name__)


class RMSFNode(PipelineNode):
    """
    RMSF calculation pipeline node.

    This node wraps GROMACS rmsf and integrates it into the pipeline.
    It can work with:
    - Predefined atom selections (e.g., 'Protein', 'C-alpha')
    - CDR regions detected by CDRDetectionNode
    - Custom index groups

    Example:
        >>> # Standard protein RMSF
        >>> pipeline = Pipeline(nodes=[
        ...     RMSFNode(selection='Protein')
        ... ])

        >>> # CDR-specific RMSF
        >>> pipeline = Pipeline(nodes=[
        ...     CDRDetectionNode(),
        ...     RMSFNode(use_cdr_regions=True)
        ... ])
    """

    def __init__(self,
                 selection: Optional[str] = None,
                 use_cdr_regions: bool = False,
                 index_file: Optional[str] = None,
                 output_basename: str = "rmsf",
                 residue_averaging: bool = True,
                 name: Optional[str] = None):
        """
        Initialize RMSF node.

        Args:
            selection: Atom selection string (e.g., 'Protein', 'C-alpha')
                      Ignored if use_cdr_regions=True
            use_cdr_regions: Use CDR regions from context (requires CDRDetectionNode)
            index_file: Path to custom index file (optional)
            output_basename: Base name for output files
            residue_averaging: Average by residue (-res flag)
            name: Node name (defaults to 'RMSFNode')
        """
        super().__init__(name=name)
        self.selection = selection
        self.use_cdr_regions = use_cdr_regions
        self.index_file = index_file
        self.output_basename = output_basename
        self.residue_averaging = residue_averaging

    def validate_inputs(self, context: PipelineContext):
        """
        Validate that required inputs exist in context.

        Args:
            context: Pipeline context

        Raises:
            PipelineError: If required inputs are missing
        """
        if not context.topology:
            raise PipelineError(
                node_name=self.name,
                reason="topology not found in context",
                context_state={'system_id': context.system_id}
            )

        if not context.trajectory_processed and not context.trajectory_raw:
            raise PipelineError(
                node_name=self.name,
                reason="No trajectory found in context",
                context_state={'system_id': context.system_id}
            )

        if self.use_cdr_regions and 'cdr_detection' not in context.metadata:
            raise PipelineError(
                node_name=self.name,
                reason="use_cdr_regions=True but CDR detection not found in context. "
                       "Ensure CDRDetectionNode runs before RMSFNode.",
                context_state={'system_id': context.system_id}
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute RMSF calculation.

        Args:
            context: Pipeline context

        Returns:
            Updated context with RMSF results
        """
        self.validate_inputs(context)

        # Use processed trajectory if available, otherwise raw
        trajectory = context.trajectory_processed or context.trajectory_raw

        logger.info(f"[{self.name}] Calculating RMSF for {context.system_id}")

        try:
            if self.use_cdr_regions:
                # Calculate RMSF for each detected CDR region
                self._calculate_cdr_rmsf(context, trajectory)
            else:
                # Calculate RMSF for single selection
                self._calculate_single_rmsf(context, trajectory)

        except Exception as e:
            raise PipelineError(
                node_name=self.name,
                reason=f"RMSF calculation failed: {str(e)}",
                context_state={
                    'system_id': context.system_id,
                    'trajectory': trajectory
                }
            ) from e

        return context

    def _calculate_single_rmsf(self, context: PipelineContext, trajectory: str):
        """Calculate RMSF for a single selection."""
        output_file = context.get_output_path(f"{self.output_basename}.xvg")

        cmd = ['gmx', 'rmsf',
               '-s', context.topology,
               '-f', trajectory,
               '-o', output_file]

        if self.residue_averaging:
            cmd.append('-res')

        if self.index_file:
            cmd.extend(['-n', self.index_file])

        # Prepare input (selection)
        gmx_input = (self.selection or 'Protein') + '\n'

        result = subprocess.run(
            cmd,
            input=gmx_input,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"gmx rmsf failed: {result.stderr}")

        context.results['rmsf'] = {
            'output_file': output_file,
            'selection': self.selection or 'Protein'
        }

        logger.info(f"  ✓ RMSF calculated: {output_file}")

    def _calculate_cdr_rmsf(self, context: PipelineContext, trajectory: str):
        """Calculate RMSF for each detected CDR region."""
        cdr_info = context.metadata['cdr_detection']
        index_file = cdr_info['index_file']

        if not Path(index_file).exists():
            raise FileNotFoundError(f"CDR index file not found: {index_file}")

        rmsf_results = []

        # Iterate through detected CDRs
        for chain_name, chain_data in cdr_info['chains'].items():
            for cdr_num, cdr_data in chain_data.get('cdrs', {}).items():
                group_name = cdr_data['index_group']
                output_file = context.get_output_path(
                    f"rmsf_cdr{cdr_num}_{chain_name}.xvg"
                )

                cmd = ['gmx', 'rmsf',
                       '-s', context.topology,
                       '-f', trajectory,
                       '-n', index_file,
                       '-o', output_file]

                if self.residue_averaging:
                    cmd.append('-res')

                result = subprocess.run(
                    cmd,
                    input=group_name + '\n',
                    capture_output=True,
                    text=True
                )

                if result.returncode == 0:
                    rmsf_results.append({
                        'chain': chain_name,
                        'cdr': cdr_num,
                        'sequence': cdr_data['sequence'],
                        'output_file': output_file
                    })
                    logger.info(f"  ✓ CDR{cdr_num}_{chain_name}: {output_file}")
                else:
                    logger.warning(f"  ✗ CDR{cdr_num}_{chain_name}: gmx rmsf failed")

        context.results['cdr_rmsf'] = {
            'n_cdrs': len(rmsf_results),
            'cdrs': rmsf_results
        }

        logger.info(f"[{self.name}] RMSF calculated for {len(rmsf_results)} CDR regions")
