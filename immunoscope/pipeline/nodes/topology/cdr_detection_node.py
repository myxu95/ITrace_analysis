"""
CDR Detection Pipeline Node - Wrapper around CDRManager.

This node uses ANARCI to automatically detect CDR regions in TCR chains
and stores the results in the pipeline context.
"""

from typing import Optional, Dict
import logging
from pathlib import Path

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError
from ....analysis import CDRManager

logger = logging.getLogger(__name__)


class CDRDetectionNode(PipelineNode):
    """
    CDR region detection pipeline node.

    This node wraps the CDRManager and integrates ANARCI-based CDR detection
    into the pipeline. It is responsible for:
    - Using ANARCI to detect CDR1, CDR2, CDR3 regions
    - Generating GROMACS-compatible index files
    - Storing CDR information in context for downstream nodes

    Example:
        >>> pipeline = Pipeline(nodes=[
        ...     CDRDetectionNode(chains={'TCR_alpha': 'chainID D'}),
        ...     RMSFNode(use_cdr_index=True),
        ... ])
    """

    def __init__(self,
                 chains: Optional[Dict[str, str]] = None,
                 include_cdrs: list = [1, 2, 3],
                 allow_fallback: bool = True,
                 numbering_scheme: str = "imgt",
                 name: Optional[str] = None):
        """
        Initialize CDR detection node.

        Args:
            chains: Dictionary mapping chain names to selection strings
                   e.g., {'TCR_alpha': 'chainID D', 'TCR_beta': 'chainID E'}
                   If None, uses default TCR chains (D, E)
            include_cdrs: List of CDR numbers to include (1, 2, 3)
            allow_fallback: Allow regex fallback if ANARCI unavailable
            numbering_scheme: ANARCI numbering scheme (imgt/kabat/chothia)
            name: Node name (defaults to 'CDRDetectionNode')
        """
        super().__init__(name=name)
        self.chains = chains
        self.include_cdrs = include_cdrs
        self.allow_fallback = allow_fallback
        self.numbering_scheme = numbering_scheme

    def validate_inputs(self, context: PipelineContext):
        """
        Validate that required inputs exist in context.

        Args:
            context: Pipeline context

        Raises:
            PipelineError: If topology is missing
        """
        if not context.topology and not context.structure_pdb:
            raise PipelineError(
                node_name=self.name,
                reason="topology/structure_pdb not found in context",
                context_state={
                    'system_id': context.system_id,
                    'has_topology': bool(context.topology),
                    'has_structure_pdb': bool(context.structure_pdb),
                }
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute CDR detection.

        Args:
            context: Pipeline context

        Returns:
            Updated context with CDR information
        """
        self.validate_inputs(context)

        logger.info(f"[{self.name}] Detecting CDR regions for {context.system_id}")

        try:
            # Initialize CDR manager
            cdr_topology = context.structure_pdb or context.topology
            cdr_manager = CDRManager(
                topology_file=cdr_topology,
                allow_fallback=self.allow_fallback,
                numbering_scheme=self.numbering_scheme
            )

            chains = self._resolve_chains(context)

            # Detect CDR regions
            cdr_data = cdr_manager.detect_all_cdr_regions(chains)

            # Generate index file
            output_dir = Path(context.get_output_path("cdr_regions.ndx")).parent
            output_dir.mkdir(parents=True, exist_ok=True)

            index_file = context.get_output_path("cdr_regions.ndx")
            cdr_manager.generate_index_file(
                output_file=index_file,
                include_cdrs=self.include_cdrs,
                chains=list(chains.keys())
            )

            # Save metadata
            metadata_file = context.get_output_path("cdr_metadata.json")
            cdr_manager.save_metadata(metadata_file, task_name=context.system_id)

            # Store results in context
            context.metadata['cdr_detection'] = {
                'method': cdr_data['detection_method'],
                'numbering_scheme': cdr_data['numbering_scheme'],
                'index_file': index_file,
                'metadata_file': metadata_file,
                'chains': {}
            }

            # Extract CDR information
            n_cdrs_detected = 0
            for chain_name, chain_data in cdr_data['chains'].items():
                if 'error' in chain_data:
                    logger.warning(f"  {chain_name}: Failed - {chain_data['error']}")
                    continue

                chain_info = {'cdrs': {}}
                for cdr_num in self.include_cdrs:
                    cdr_key = f'cdr{cdr_num}'
                    if cdr_key in chain_data:
                        cdr_info = chain_data[cdr_key]
                        chain_info['cdrs'][cdr_num] = {
                            'sequence': cdr_info['sequence'],
                            'residue_range': cdr_info['residue_range'],
                            'index_group': f"CDR{cdr_num}_{chain_name}_CA"
                        }
                        n_cdrs_detected += 1
                        logger.info(f"  ✓ {chain_name} CDR{cdr_num}: {cdr_info['sequence']}")
                chain_info['chain_id'] = self._selection_to_chain_id(chains.get(chain_name))
                context.metadata['cdr_detection']['chains'][chain_name] = chain_info

            # Summary
            context.results['cdr_detection'] = {
                'n_chains': len(chains),
                'n_cdrs_detected': n_cdrs_detected,
                'detection_method': cdr_data['detection_method']
            }

            logger.info(f"[{self.name}] Detected {n_cdrs_detected} CDR regions")

        except Exception as e:
            raise PipelineError(
                node_name=self.name,
                reason=f"CDR detection failed: {str(e)}",
                context_state={
                    'system_id': context.system_id,
                    'topology': context.topology,
                    'structure_pdb': context.structure_pdb
                }
            ) from e

        return context

    def _resolve_chains(self, context: PipelineContext) -> Dict[str, str]:
        if self.chains is not None:
            return self.chains

        chain_mapping = context.metadata.get('chain_mapping', {})
        tcr_alpha = chain_mapping.get('tcr_alpha')
        tcr_beta = chain_mapping.get('tcr_beta')
        if tcr_alpha and tcr_beta:
            return {
                'TCR_alpha': f'chainID {tcr_alpha}',
                'TCR_beta': f'chainID {tcr_beta}'
            }

        return {
            'TCR_alpha': 'chainID D',
            'TCR_beta': 'chainID E'
        }

    @staticmethod
    def _selection_to_chain_id(selection: Optional[str]) -> Optional[str]:
        if not selection:
            return None
        if selection.startswith("chainID "):
            return selection.split(" ", 1)[1].strip()
        return selection
