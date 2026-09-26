"""
Selection String Builder for MDAnalysis

Generates MDAnalysis selection strings from chain identification results.
Handles both chainID (PDB) and segname (TPR/GRO) naming conventions.

Author: ImmunoScope Development Team
Date: 2026-03-10
"""

import logging
from typing import Dict, Optional

from ..analysis.topology import UnifiedChainInfo

logger = logging.getLogger(__name__)


class SelectionStringBuilder:
    """
    Build MDAnalysis selection strings from chain identification results.

    Automatically handles:
    - chainID format (PDB files): 'chainID A'
    - segname format (TPR/GRO files): 'segname PROA'
    """

    # Standard residue ranges for MHC groove
    MHC_GROOVE_ALPHA1 = (50, 86)
    MHC_GROOVE_ALPHA2 = (140, 176)

    # Standard residue range for TCR V domain
    TCR_V_DOMAIN = (1, 115)

    def __init__(self):
        """Initialize selection string builder."""
        pass

    def build_all_selections(
        self,
        identifications: Dict[str, UnifiedChainInfo]
    ) -> Dict[str, str]:
        """
        Build all selection strings for docking angle analysis.

        Generates selections for:
        - mhc_selection: MHC alpha chain (full or groove region)
        - tcr_alpha_selection: TCR alpha chain (if present)
        - tcr_beta_selection: TCR beta chain (if present)
        - peptide_selection: Peptide chain (if present)
        - beta2m_selection: Beta2-microglobulin chain (if present)

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Results from ChainIdentificationAdapter.identify_chains()

        Returns
        -------
        selections : Dict[str, str]
            Dictionary of selection strings

        Example
        -------
        >>> builder = SelectionStringBuilder()
        >>> selections = builder.build_all_selections(identifications)
        >>> print(selections['mhc_selection'])
        'segname PROA'
        >>> print(selections['tcr_alpha_selection'])
        'segname PROD'
        """
        selections = {}

        # Get naming convention (chainID or segname)
        naming_convention = self._detect_naming_convention(identifications)

        # Build selections for each component
        for chain_id, info in identifications.items():
            if info.chain_type == 'HLA_alpha':
                selections['mhc_selection'] = self._build_chain_selection(
                    chain_id, naming_convention
                )
                # Also create groove-specific selection
                selections['mhc_groove_selection'] = self._build_mhc_groove_selection(
                    chain_id, naming_convention, info.residue_range
                )

            elif info.chain_type == 'TCR_alpha':
                selections['tcr_alpha_selection'] = self._build_chain_selection(
                    chain_id, naming_convention
                )

            elif info.chain_type == 'TCR_beta':
                selections['tcr_beta_selection'] = self._build_chain_selection(
                    chain_id, naming_convention
                )

            elif info.chain_type == 'peptide':
                selections['peptide_selection'] = self._build_chain_selection(
                    chain_id, naming_convention
                )

            elif info.chain_type == 'beta2m':
                selections['beta2m_selection'] = self._build_chain_selection(
                    chain_id, naming_convention
                )

        # Log generated selections
        logger.info("Generated MDAnalysis selection strings:")
        for key, value in selections.items():
            logger.info(f"  {key}: {value}")

        return selections

    def _detect_naming_convention(
        self,
        identifications: Dict[str, UnifiedChainInfo]
    ) -> str:
        """
        Detect naming convention from identification results.

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Chain identification results

        Returns
        -------
        convention : str
            Either 'chainID' or 'segname'
        """
        # Use the naming convention from the first entry
        if identifications:
            first_info = next(iter(identifications.values()))
            return first_info.naming_convention
        else:
            # Default to chainID
            return 'chainID'

    def _build_chain_selection(self, chain_id: str, convention: str) -> str:
        """
        Build basic chain selection string.

        Parameters
        ----------
        chain_id : str
            Chain identifier (chainID or segname)
        convention : str
            'chainID' or 'segname'

        Returns
        -------
        selection : str
            MDAnalysis selection string

        Example
        -------
        >>> builder._build_chain_selection('A', 'chainID')
        'chainID A'
        >>> builder._build_chain_selection('PROA', 'segname')
        'segname PROA'
        """
        if convention == 'chainID':
            return f"chainID {chain_id}"
        elif convention == 'segname':
            return f"segname {chain_id}"
        else:
            raise ValueError(f"Unknown naming convention: {convention}")

    def _build_mhc_groove_selection(
        self,
        chain_id: str,
        convention: str,
        residue_range: tuple
    ) -> str:
        """
        Build MHC groove region selection (α1 + α2 domains).

        Uses standard ranges:
        - α1 domain: residues 50-86
        - α2 domain: residues 140-176

        Parameters
        ----------
        chain_id : str
            MHC chain identifier
        convention : str
            'chainID' or 'segname'
        residue_range : tuple
            (first_resid, last_resid) for the chain

        Returns
        -------
        selection : str
            MDAnalysis selection string for groove CA atoms

        Example
        -------
        >>> builder._build_mhc_groove_selection('A', 'chainID', (1, 300))
        'chainID A and (resid 50:86 or resid 140:176) and name CA'
        """
        first_resid, last_resid = residue_range
        alpha1_start, alpha1_end = self.MHC_GROOVE_ALPHA1
        alpha2_start, alpha2_end = self.MHC_GROOVE_ALPHA2

        # Check if residue range covers groove regions
        if last_resid < alpha2_end:
            logger.warning(
                f"MHC chain {chain_id} may not fully cover α2 domain "
                f"(ends at {last_resid}, expected up to {alpha2_end})"
            )

        # Build chain part
        if convention == 'chainID':
            chain_part = f"chainID {chain_id}"
        else:
            chain_part = f"segname {chain_id}"

        # Build selection
        selection = (
            f"{chain_part} and "
            f"(resid {alpha1_start}:{alpha1_end} or resid {alpha2_start}:{alpha2_end}) and "
            f"name CA"
        )

        return selection

    def build_mhc_selection(
        self,
        identifications: Dict[str, UnifiedChainInfo],
        use_groove_only: bool = False
    ) -> Optional[str]:
        """
        Build MHC selection string.

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Chain identification results
        use_groove_only : bool, default=False
            If True, select only groove region (α1 + α2)
            If False, select entire MHC chain

        Returns
        -------
        selection : str or None
            MDAnalysis selection string, or None if no MHC found
        """
        for chain_id, info in identifications.items():
            if info.chain_type == 'HLA_alpha':
                convention = info.naming_convention

                if use_groove_only:
                    return self._build_mhc_groove_selection(
                        chain_id, convention, info.residue_range
                    )
                else:
                    return self._build_chain_selection(chain_id, convention)

        logger.warning("No HLA-alpha chain found in identifications")
        return None

    def build_tcr_selection(
        self,
        identifications: Dict[str, UnifiedChainInfo],
        chain_type: str = 'alpha'
    ) -> Optional[str]:
        """
        Build TCR selection string.

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Chain identification results
        chain_type : str
            'alpha' or 'beta'

        Returns
        -------
        selection : str or None
            MDAnalysis selection string, or None if chain not found
        """
        target_type = f'TCR_{chain_type}'

        for chain_id, info in identifications.items():
            if info.chain_type == target_type:
                return self._build_chain_selection(chain_id, info.naming_convention)

        logger.warning(f"No TCR-{chain_type} chain found in identifications")
        return None

    def build_peptide_selection(
        self,
        identifications: Dict[str, UnifiedChainInfo]
    ) -> Optional[str]:
        """
        Build peptide selection string.

        Parameters
        ----------
        identifications : Dict[str, UnifiedChainInfo]
            Chain identification results

        Returns
        -------
        selection : str or None
            MDAnalysis selection string, or None if peptide not found
        """
        for chain_id, info in identifications.items():
            if info.chain_type == 'peptide':
                return self._build_chain_selection(chain_id, info.naming_convention)

        logger.warning("No peptide chain found in identifications")
        return None
