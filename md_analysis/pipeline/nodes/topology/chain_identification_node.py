"""
Chain Identification Node - Intelligent chain identification for pHLA-TCR complexes.

This node wraps IntelligentChainIdentifier and adds it to the pipeline context.
"""

import logging
from pathlib import Path
from typing import Optional

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....analysis.topology import IntelligentChainIdentifier

logger = logging.getLogger(__name__)


class ChainIdentificationNode(PipelineNode):
    """
    Chain identification node using ANARCI-based intelligent identification.

    This node:
    1. Uses IntelligentChainIdentifier to identify chains
    2. Stores chain mapping in context.metadata['chain_mapping']
    3. Handles identification failures gracefully

    Example:
        >>> node = ChainIdentificationNode(method="anarci")
        >>> context = node.execute(context)
        >>> chain_mapping = context.metadata['chain_mapping']
        >>> # {'mhc_alpha': 'A', 'b2m': 'B', 'peptide': 'C',
        >>> #  'tcr_alpha': 'D', 'tcr_beta': 'E'}
    """

    def __init__(self, method: str = "anarci", fallback_to_heuristic: bool = True):
        """
        Initialize chain identification node.

        Args:
            method: Identification method ('anarci' or 'heuristic')
            fallback_to_heuristic: If True, fall back to heuristic method on failure
        """
        super().__init__(name="ChainIdentification")
        self.method = method
        self.fallback_to_heuristic = fallback_to_heuristic
        self.use_anarci = (method == "anarci")

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute chain identification.

        Args:
            context: Pipeline context

        Returns:
            Updated context with chain_mapping in metadata

        Raises:
            ValueError: If structure_pdb is not provided
        """
        logger.info(f"[{context.system_id}] Starting chain identification")

        # Validate input
        if not context.structure_pdb:
            error_msg = "structure_pdb is required for chain identification"
            logger.error(f"[{context.system_id}] {error_msg}")
            context.add_error(error_msg)
            context.should_stop = True
            return context

        if not Path(context.structure_pdb).exists():
            error_msg = f"Structure file not found: {context.structure_pdb}"
            logger.error(f"[{context.system_id}] {error_msg}")
            context.add_error(error_msg)
            context.should_stop = True
            return context

        try:
            # Initialize identifier
            identifier = IntelligentChainIdentifier(use_anarci=self.use_anarci)

            # Identify chains
            chain_results = identifier.identify_chains(context.structure_pdb)

            # Build chain mapping from results
            chain_mapping = {}
            for chain_id, chain_info in chain_results.items():
                chain_type = chain_info.chain_type

                # Map chain_type to standard names
                if chain_type == 'HLA_alpha':
                    chain_mapping['mhc_alpha'] = chain_id
                elif chain_type == 'beta2m':
                    chain_mapping['b2m'] = chain_id
                elif chain_type == 'peptide':
                    chain_mapping['peptide'] = chain_id
                elif chain_type == 'TCR_alpha':
                    chain_mapping['tcr_alpha'] = chain_id
                elif chain_type == 'TCR_beta':
                    chain_mapping['tcr_beta'] = chain_id

            # Validate we have all required chains
            required_chains = ['mhc_alpha', 'b2m', 'peptide', 'tcr_alpha', 'tcr_beta']
            missing_chains = [c for c in required_chains if c not in chain_mapping]

            if missing_chains:
                warning_msg = f"Missing chains: {missing_chains}"
                logger.warning(f"[{context.system_id}] {warning_msg}")
                context.add_warning(warning_msg)

                # Apply fallback if enabled
                if self.fallback_to_heuristic:
                    logger.info(f"[{context.system_id}] Applying fallback chain assignment")
                    chain_mapping = self._fallback_chain_assignment()
                    context.add_warning("Using fallback chain assignment (A-E)")

            # Store in context
            context.metadata['chain_mapping'] = chain_mapping
            context.metadata['chain_identification_method'] = self.method

            logger.info(f"[{context.system_id}] Chain identification successful")
            logger.debug(f"[{context.system_id}] Chain mapping: {chain_mapping}")

        except Exception as e:
            error_msg = f"Chain identification failed: {str(e)}"
            logger.exception(f"[{context.system_id}] {error_msg}")
            context.add_error(error_msg)

            if self.fallback_to_heuristic:
                logger.info(f"[{context.system_id}] Applying fallback chain assignment")
                context.metadata['chain_mapping'] = self._fallback_chain_assignment()
                context.add_warning("Using fallback chain assignment due to error")
            else:
                context.should_stop = True

        return context

    @staticmethod
    def _fallback_chain_assignment() -> dict:
        """
        Fallback chain assignment (standard A-E mapping).

        Returns:
            Standard chain mapping
        """
        return {
            'mhc_alpha': 'A',
            'b2m': 'B',
            'peptide': 'C',
            'tcr_alpha': 'D',
            'tcr_beta': 'E'
        }
