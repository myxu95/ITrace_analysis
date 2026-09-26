"""
Chain Identification Adapter for pHLA-TCR Complexes

Unified interface for automatic chain identification from both PDB and TPR/GRO files.
Routes to appropriate identifier (IntelligentChainIdentifier or TPRChainSequenceExtractor).

Author: ITrace Development Team
Date: 2026-03-10
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from .intelligent_chain_identifier import IntelligentChainIdentifier, ChainIdentification
from .tpr_chain_extractor import TPRChainSequenceExtractor, UnifiedChainInfo

logger = logging.getLogger(__name__)


class ChainIdentificationAdapter:
    """
    Unified adapter for automatic chain identification.

    Supports both PDB and TPR/GRO file formats, automatically routing
    to the appropriate identifier based on file extension.
    """

    def __init__(self, use_anarci: bool = True):
        """
        Initialize chain identification adapter.

        Parameters
        ----------
        use_anarci : bool
            Whether to use ANARCI for TCR identification
        """
        self.use_anarci = use_anarci
        self.pdb_identifier = IntelligentChainIdentifier(use_anarci=use_anarci)
        self.tpr_extractor = TPRChainSequenceExtractor(use_anarci=use_anarci)

    def identify_chains(
        self,
        topology_file: str,
        structure_file: Optional[str] = None
    ) -> Dict[str, UnifiedChainInfo]:
        """
        Identify all protein chains from topology file.

        Automatically detects file format and routes to appropriate identifier.

        Parameters
        ----------
        topology_file : str
            Path to topology/structure file (PDB, TPR, GRO)
        structure_file : str, optional
            Path to structure file (for TPR + GRO combination)

        Returns
        -------
        identifications : Dict[str, UnifiedChainInfo]
            Dictionary mapping chain_id to UnifiedChainInfo

        Example
        -------
        >>> adapter = ChainIdentificationAdapter()
        >>> # PDB file
        >>> results_pdb = adapter.identify_chains('1ao7.pdb')
        >>> # TPR file
        >>> results_tpr = adapter.identify_chains('md.tpr', 'md.gro')
        """
        topology_path = Path(topology_file)
        suffix = topology_path.suffix.lower()

        logger.info(f"Identifying chains from: {topology_file}")

        # Route based on file type
        if suffix == '.pdb':
            return self._identify_from_pdb(topology_file)
        elif suffix in ['.tpr', '.gro']:
            return self._identify_from_tpr_gro(topology_file, structure_file)
        else:
            raise ValueError(
                f"Unsupported file format: {suffix}. "
                f"Expected: .pdb, .tpr, or .gro"
            )

    def _identify_from_pdb(self, pdb_file: str) -> Dict[str, UnifiedChainInfo]:
        """
        Identify chains from PDB file using IntelligentChainIdentifier.

        Parameters
        ----------
        pdb_file : str
            Path to PDB file

        Returns
        -------
        identifications : Dict[str, UnifiedChainInfo]
            Converted to UnifiedChainInfo format
        """
        logger.info("Using IntelligentChainIdentifier for PDB file")

        # Run identification
        pdb_results = self.pdb_identifier.identify_chains(pdb_file)

        # Convert ChainIdentification to UnifiedChainInfo
        unified_results = {}
        for chain_id, pdb_info in pdb_results.items():
            unified_results[chain_id] = UnifiedChainInfo(
                chain_id=chain_id,
                chain_type=pdb_info.chain_type,
                length=pdb_info.length,
                sequence=pdb_info.sequence,
                confidence=pdb_info.confidence,
                residue_range=(1, pdb_info.length),  # PDB doesn't store this
                naming_convention='chainID',
                anarci_result=pdb_info.anarci_result
            )

        return unified_results

    def _identify_from_tpr_gro(
        self,
        topology_file: str,
        structure_file: Optional[str]
    ) -> Dict[str, UnifiedChainInfo]:
        """
        Identify chains from TPR/GRO using TPRChainSequenceExtractor.

        Parameters
        ----------
        topology_file : str
            Path to TPR or GRO file
        structure_file : str, optional
            Path to GRO file (if topology is TPR)

        Returns
        -------
        identifications : Dict[str, UnifiedChainInfo]
            Already in UnifiedChainInfo format
        """
        logger.info("Using TPRChainSequenceExtractor for TPR/GRO file")

        # Run identification
        results = self.tpr_extractor.identify_chains(topology_file, structure_file)

        return results

    def validate_identification(
        self,
        identifications: Dict[str, UnifiedChainInfo],
        strict: bool = False
    ) -> Tuple[bool, List[str]]:
        """
        Validate chain identification results.

        Validation modes:
        - Strict mode: Requires all 5 components (HLA-α, β2m, peptide, TCR-α, TCR-β)
        - Graceful mode (default): Requires HLA-α + at least 1 TCR chain

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Results from identify_chains()
        strict : bool, default=False
            If True, require all 5 components. If False, allow missing peptide/β2m

        Returns
        -------
        is_valid : bool
            Whether identification meets requirements
        messages : List[str]
            Warning/error messages

        Example
        -------
        >>> adapter = ChainIdentificationAdapter()
        >>> results = adapter.identify_chains('md.tpr')
        >>> is_valid, messages = adapter.validate_identification(results, strict=False)
        >>> if not is_valid:
        ...     for msg in messages:
        ...         print(f"Error: {msg}")
        """
        messages = []

        # Count identified chain types
        chain_types = [info.chain_type for info in identifications.values()]
        has_hla_alpha = 'HLA_alpha' in chain_types
        has_tcr_alpha = 'TCR_alpha' in chain_types
        has_tcr_beta = 'TCR_beta' in chain_types
        has_peptide = 'peptide' in chain_types
        has_beta2m = 'beta2m' in chain_types
        has_unknown = 'unknown' in chain_types

        # Strict mode: require all 5 components
        if strict:
            if not has_hla_alpha:
                messages.append("Missing HLA-alpha chain")
            if not has_tcr_alpha:
                messages.append("Missing TCR-alpha chain")
            if not has_tcr_beta:
                messages.append("Missing TCR-beta chain")
            if not has_peptide:
                messages.append("Missing peptide chain")
            if not has_beta2m:
                messages.append("Missing beta2-microglobulin chain")

            is_valid = len(messages) == 0

            if has_unknown:
                messages.append("Warning: Some chains could not be identified")

        # Graceful mode: require HLA-α + at least 1 TCR chain
        else:
            # Critical requirements
            if not has_hla_alpha:
                messages.append("ERROR: Missing HLA-alpha chain (required)")
            if not has_tcr_alpha and not has_tcr_beta:
                messages.append("ERROR: Missing both TCR chains (at least 1 required)")

            # Optional components (warnings only)
            if not has_peptide:
                messages.append("WARNING: Missing peptide chain (optional in graceful mode)")
            if not has_beta2m:
                messages.append("WARNING: Missing beta2-microglobulin chain (optional in graceful mode)")
            if has_unknown:
                messages.append("WARNING: Some chains could not be identified")

            # Determine validity (only critical requirements)
            critical_errors = [msg for msg in messages if msg.startswith("ERROR")]
            is_valid = len(critical_errors) == 0

        # Log validation results
        if is_valid:
            logger.info("✓ Chain identification validation passed")
            for msg in messages:
                if msg.startswith("WARNING"):
                    logger.warning(msg)
        else:
            logger.error("✗ Chain identification validation failed")
            for msg in messages:
                if msg.startswith("ERROR"):
                    logger.error(msg)
                else:
                    logger.warning(msg)

        return is_valid, messages

    def get_chain_summary(self, identifications: Dict[str, UnifiedChainInfo]) -> str:
        """
        Generate a human-readable summary of identified chains.

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Results from identify_chains()

        Returns
        -------
        summary : str
            Formatted summary text

        Example
        -------
        >>> adapter = ChainIdentificationAdapter()
        >>> results = adapter.identify_chains('md.tpr')
        >>> print(adapter.get_chain_summary(results))
        """
        lines = ["Chain Identification Summary:", "=" * 50]

        # Sort by chain type for consistent ordering
        type_order = {
            'HLA_alpha': 1,
            'beta2m': 2,
            'peptide': 3,
            'TCR_alpha': 4,
            'TCR_beta': 5,
            'unknown': 6
        }

        sorted_items = sorted(
            identifications.items(),
            key=lambda x: type_order.get(x[1].chain_type, 99)
        )

        for chain_id, info in sorted_items:
            confidence_str = f"{info.confidence:.2f}"
            anarci_str = " (ANARCI)" if info.anarci_result else ""

            lines.append(
                f"  {chain_id:6s}: {info.chain_type:12s} | "
                f"{info.length:3d} AA | conf={confidence_str}{anarci_str}"
            )

        lines.append("=" * 50)

        return "\n".join(lines)
