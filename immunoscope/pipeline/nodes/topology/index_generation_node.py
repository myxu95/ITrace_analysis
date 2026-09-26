"""
Index Generation Node - Generate GROMACS index files for group selections.

This node creates index files based on chain mapping and group definitions.
"""

import logging
from pathlib import Path
from typing import Dict, Optional
import MDAnalysis as mda

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext

logger = logging.getLogger(__name__)


class IndexGenerationNode(PipelineNode):
    """
    Generate GROMACS index file based on chain mapping.

    This node:
    1. Reads chain_mapping from context.metadata
    2. Generates index file with user-defined groups
    3. Stores index file path in context

    Example:
        >>> node = IndexGenerationNode(group_definitions={
        ...     'pHLA': 'chainID {mhc_alpha} {b2m} {peptide}',
        ...     'TCR': 'chainID {tcr_alpha} {tcr_beta}'
        ... })
        >>> context = node.execute(context)
        >>> print(context.index_file)
    """

    def __init__(self, group_definitions: Dict[str, str], output_name: str = "groups.ndx"):
        """
        Initialize index generation node.

        Args:
            group_definitions: Dictionary mapping group name to MDAnalysis selection template
                              Template can use {mhc_alpha}, {b2m}, {peptide}, {tcr_alpha}, {tcr_beta}
            output_name: Name of output index file

        Example:
            >>> group_definitions = {
            ...     'pHLA': 'chainID {mhc_alpha} {b2m} {peptide}',
            ...     'TCR': 'chainID {tcr_alpha} {tcr_beta}',
            ...     'TCR_CA': 'chainID {tcr_alpha} {tcr_beta} and name CA'
            ... }
        """
        super().__init__(name="IndexGeneration")
        self.group_definitions = group_definitions
        self.output_name = output_name

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute index file generation.

        Args:
            context: Pipeline context (requires structure_pdb and chain_mapping)

        Returns:
            Updated context with index_file path
        """
        logger.info(f"[{context.system_id}] Starting index file generation")

        # Validate input
        if not context.structure_pdb:
            error_msg = "structure_pdb is required for index generation"
            logger.error(f"[{context.system_id}] {error_msg}")
            context.add_error(error_msg)
            context.should_stop = True
            return context

        chain_mapping = context.metadata.get('chain_mapping')
        if not chain_mapping:
            error_msg = "chain_mapping not found in metadata (run ChainIdentificationNode first)"
            logger.error(f"[{context.system_id}] {error_msg}")
            context.add_error(error_msg)
            context.should_stop = True
            return context

        try:
            # Determine output path
            if context.output_dir:
                output_path = Path(context.output_dir) / self.output_name
            else:
                output_path = Path(context.structure_pdb).parent / self.output_name

            # Load structure
            u = mda.Universe(context.structure_pdb)

            # Generate index file
            with open(output_path, 'w') as f:
                for group_name, selection_template in self.group_definitions.items():
                    try:
                        # Fill template with chain mapping
                        selection = selection_template.format(**chain_mapping)

                        # Select atoms
                        atoms = u.select_atoms(selection)

                        if atoms.n_atoms == 0:
                            warning_msg = f"Group '{group_name}' selected 0 atoms (selection: {selection})"
                            logger.warning(f"[{context.system_id}] {warning_msg}")
                            context.add_warning(warning_msg)
                            continue

                        # Write group
                        f.write(f"[ {group_name} ]\n")
                        indices = [a.index + 1 for a in atoms]  # GROMACS uses 1-based indexing

                        # Write indices (10 per line)
                        for i in range(0, len(indices), 10):
                            f.write(" ".join(map(str, indices[i:i+10])) + "\n")
                        f.write("\n")

                        logger.info(f"[{context.system_id}] Group '{group_name}': {len(indices)} atoms")

                    except KeyError as e:
                        error_msg = f"Group '{group_name}': Missing chain in mapping: {e}"
                        logger.error(f"[{context.system_id}] {error_msg}")
                        context.add_error(error_msg)
                        context.should_stop = True
                        return context

                    except Exception as e:
                        error_msg = f"Group '{group_name}': Selection failed: {str(e)}"
                        logger.error(f"[{context.system_id}] {error_msg}")
                        context.add_error(error_msg)
                        context.should_stop = True
                        return context

            # Store index file path
            context.index_file = str(output_path)
            logger.info(f"[{context.system_id}] Index file generated: {output_path}")

        except Exception as e:
            error_msg = f"Index generation failed: {str(e)}"
            logger.exception(f"[{context.system_id}] {error_msg}")
            context.add_error(error_msg)
            context.should_stop = True

        return context
