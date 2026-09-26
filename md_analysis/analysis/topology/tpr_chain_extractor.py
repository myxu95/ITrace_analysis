"""
TPR/GRO Chain Sequence Extractor for pHLA-TCR Complexes

Extracts protein chain sequences from TPR/GRO files using MDAnalysis,
then identifies chains using ANARCI + length heuristics (same strategy as IntelligentChainIdentifier).

Author: ITrace Development Team
Date: 2026-03-10
"""

import logging
from typing import Dict, Optional
from dataclasses import dataclass
import MDAnalysis as mda

from .cdr_manager import ANARCIWrapper

logger = logging.getLogger(__name__)


THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


@dataclass
class UnifiedChainInfo:
    """Unified chain information structure for both PDB and TPR/GRO files."""
    chain_id: str
    chain_type: str
    length: int
    sequence: str
    confidence: float
    residue_range: tuple
    naming_convention: str
    anarci_result: Optional[Dict] = None


class TPRChainSequenceExtractor:
    """
    Extract and identify protein chains from TPR/GRO files.

    Uses the same identification strategy as IntelligentChainIdentifier:
    1. Peptide: <=20 AA (definitive)
    2. Beta2m: 90-110 AA (definitive)
    3. ANARCI for TCR alpha/beta identification
    4. HLA-alpha by elimination
    """

    # Length thresholds (same as IntelligentChainIdentifier)
    PEPTIDE_MAX_LENGTH = 20
    BETA2M_MIN_LENGTH = 90
    BETA2M_MAX_LENGTH = 110
    HLA_ALPHA_MIN_LENGTH = 230
    HLA_ALPHA_MAX_LENGTH = 310

    def __init__(self, use_anarci: bool = True):
        """
        Initialize TPR chain extractor.

        Parameters
        ----------
        use_anarci : bool
            Whether to use ANARCI for TCR identification
        """
        self.use_anarci = use_anarci

        if use_anarci:
            try:
                self.anarci_wrapper = ANARCIWrapper(allow_fallback=False)
                logger.info("ANARCI available for TCR chain identification")
            except RuntimeError as e:
                logger.warning(f"ANARCI not available: {e}")
                logger.warning("Will use fallback length-based identification")
                self.use_anarci = False
                self.anarci_wrapper = None
        else:
            self.anarci_wrapper = None

    def identify_chains(
        self,
        topology_file: str,
        structure_file: Optional[str] = None
    ) -> Dict[str, UnifiedChainInfo]:
        """
        Identify all protein chains from TPR/GRO file.

        Parameters
        ----------
        topology_file : str
            Path to topology file (TPR)
        structure_file : str, optional
            Path to structure file (GRO). If None, use topology as structure

        Returns
        -------
        identifications : Dict[str, UnifiedChainInfo]
            Dictionary mapping segname to chain information

        Example
        -------
        >>> extractor = TPRChainSequenceExtractor()
        >>> results = extractor.identify_chains('md.tpr', 'md.gro')
        >>> for segname, info in results.items():
        ...     print(f"{segname}: {info.chain_type} ({info.length} AA)")
        """
        # Load universe
        if structure_file:
            universe = mda.Universe(topology_file, structure_file)
        else:
            universe = mda.Universe(topology_file)

        logger.info(f"Loaded topology: {topology_file}")

        # Extract protein chains by segname
        chains_data = self._extract_protein_chains(universe)

        if len(chains_data) != 5:
            logger.warning(
                f"Expected 5 chains in pHLA-TCR complex, found {len(chains_data)}"
            )

        # Sort by length
        sorted_chains = sorted(
            chains_data.items(),
            key=lambda x: x[1]['length']
        )

        identifications = {}

        # Step 1: Identify peptide (<=20 AA)
        segname, chain_data = sorted_chains[0]
        if chain_data['length'] <= self.PEPTIDE_MAX_LENGTH:
            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='peptide',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=1.0,
                residue_range=chain_data['residue_range'],
                naming_convention='segname'
            )
            logger.info(f"[DEFINITIVE] {segname} = peptide ({chain_data['length']} AA, <=20)")
        else:
            logger.warning(
                f"Shortest chain {segname} has {chain_data['length']} AA, "
                f"exceeds peptide threshold (<=20 AA)"
            )
            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='unknown',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=0.5,
                residue_range=chain_data['residue_range'],
                naming_convention='segname'
            )

        # Step 2: Identify beta2m (90-110 AA)
        beta2m_found = False
        for segname, chain_data in sorted_chains[1:]:
            if self.BETA2M_MIN_LENGTH <= chain_data['length'] <= self.BETA2M_MAX_LENGTH:
                identifications[segname] = UnifiedChainInfo(
                    chain_id=segname,
                    chain_type='beta2m',
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    confidence=1.0,
                    residue_range=chain_data['residue_range'],
                    naming_convention='segname'
                )
                logger.info(f"[DEFINITIVE] {segname} = beta2m ({chain_data['length']} AA, 90-110 range)")
                beta2m_found = True
                break

        if not beta2m_found:
            logger.warning("Could not identify beta2m chain (expected 90-110 AA range)")

        # Step 3: Identify remaining chains using ANARCI
        remaining_chains = {
            sid: cdata for sid, cdata in chains_data.items()
            if sid not in identifications
        }

        if len(remaining_chains) == 3:
            tcr_hla_identifications = self._identify_tcr_and_hla(remaining_chains)
            identifications.update(tcr_hla_identifications)
        else:
            logger.warning(
                f"Expected 3 remaining chains (2 TCR + 1 HLA-alpha), "
                f"found {len(remaining_chains)}"
            )
            # Fallback: mark as unknown
            for segname, chain_data in remaining_chains.items():
                identifications[segname] = UnifiedChainInfo(
                    chain_id=segname,
                    chain_type='unknown',
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    confidence=0.0,
                    residue_range=chain_data['residue_range'],
                    naming_convention='segname'
                )

        return identifications

    def _extract_protein_chains(self, universe: mda.Universe) -> Dict[str, dict]:
        """
        Extract protein chain sequences from MDAnalysis Universe.

        Uses segname (segment IDs) to group chains.

        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis Universe object

        Returns
        -------
        chains_data : Dict[str, dict]
            Dictionary mapping segname to chain data:
            - 'sequence': single-letter amino acid sequence
            - 'length': number of residues
            - 'residue_range': (first_resid, last_resid)
        """
        protein = universe.select_atoms("protein")

        if len(protein) == 0:
            raise ValueError("No protein atoms found in structure")

        chains_data = {}

        # Group by segments (segids)
        for segment in protein.segments:
            segname = segment.segid

            # Select protein residues in this segment
            seg_residues = protein.select_atoms(f"segid {segname}").residues

            if len(seg_residues) == 0:
                continue

            # Extract sequence (convert three-letter to one-letter)
            sequence = ""
            for residue in seg_residues:
                resname = residue.resname
                if resname in THREE_TO_ONE:
                    sequence += THREE_TO_ONE[resname]
                else:
                    logger.warning(f"Unknown residue {resname} in segment {segname}, using 'X'")
                    sequence += 'X'

            # Get residue range
            first_resid = seg_residues[0].resid
            last_resid = seg_residues[-1].resid

            chains_data[segname] = {
                'sequence': sequence,
                'length': len(sequence),
                'residue_range': (first_resid, last_resid)
            }

            logger.debug(f"Extracted {segname}: {len(sequence)} AA (resid {first_resid}-{last_resid})")

        logger.info(f"Extracted {len(chains_data)} protein chains from TPR/GRO")

        return chains_data

    def _identify_tcr_and_hla(
        self,
        chains: Dict[str, Dict]
    ) -> Dict[str, UnifiedChainInfo]:
        """
        Identify TCR alpha, TCR beta, and HLA-alpha from 3 remaining chains.

        Strategy:
        1. Send all 3 chains to ANARCI for TCR identification
        2. ANARCI identifies 2 TCR chains (alpha and beta)
        3. Remaining chain is HLA-alpha

        Parameters
        ----------
        chains : Dict[str, Dict]
            Dictionary of 3 remaining chains

        Returns
        -------
        identifications : Dict[str, UnifiedChainInfo]
            Dictionary of chain identifications
        """
        identifications = {}

        if not self.use_anarci or self.anarci_wrapper is None:
            return self._identify_by_length_fallback(chains)

        # Run ANARCI on all 3 chains
        logger.info("Running ANARCI on 3 remaining chains to identify TCR...")
        anarci_results = {}
        for segname, chain_data in chains.items():
            try:
                result = self.anarci_wrapper.run_anarci(
                    sequence=chain_data['sequence'],
                    chain_type='TCR'
                )
                anarci_results[segname] = result
                chain_type_detected = result.get('chain_type', 'unknown')
                logger.info(
                    f"  {segname} ({chain_data['length']} AA): "
                    f"ANARCI result = {chain_type_detected}"
                )
            except Exception as e:
                logger.warning(f"  {segname}: ANARCI failed ({e})")
                anarci_results[segname] = None

        # Classify chains based on ANARCI results
        tcr_alpha_candidates = []
        tcr_beta_candidates = []
        hla_alpha_candidates = []

        for segname, result in anarci_results.items():
            chain_data = chains[segname]

            if result and 'chain_type' in result:
                chain_type_anarci = result['chain_type'].lower()

                # TCR alpha-like: TCR_alpha or TCR_delta
                if ('alpha' in chain_type_anarci or 'delta' in chain_type_anarci) and 'tcr' in chain_type_anarci:
                    tcr_alpha_candidates.append((segname, chain_data, result))
                # TCR beta-like: TCR_beta or TCR_gamma
                elif ('beta' in chain_type_anarci or 'gamma' in chain_type_anarci) and 'tcr' in chain_type_anarci:
                    tcr_beta_candidates.append((segname, chain_data, result))
                else:
                    # Likely HLA-alpha
                    hla_alpha_candidates.append((segname, chain_data, result))
            else:
                # No TCR recognition - likely HLA-alpha
                hla_alpha_candidates.append((segname, chain_data, None))

        # Assign TCR-alpha
        if len(tcr_alpha_candidates) == 1:
            segname, chain_data, result = tcr_alpha_candidates[0]
            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='TCR_alpha',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=0.95,
                residue_range=chain_data['residue_range'],
                naming_convention='segname',
                anarci_result=result
            )
            logger.info(f"[ANARCI] {segname} = TCR-alpha")
        elif len(tcr_alpha_candidates) > 1:
            logger.warning(f"Multiple TCR-alpha candidates: {[c[0] for c in tcr_alpha_candidates]}")
            # Use shortest as alpha
            tcr_alpha_candidates.sort(key=lambda x: x[1]['length'])
            segname, chain_data, result = tcr_alpha_candidates[0]
            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='TCR_alpha',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=0.7,
                residue_range=chain_data['residue_range'],
                naming_convention='segname',
                anarci_result=result
            )
            logger.info(f"[ANARCI] {segname} = TCR-alpha (multiple candidates, using shortest)")

        # Assign TCR-beta
        if len(tcr_beta_candidates) == 1:
            segname, chain_data, result = tcr_beta_candidates[0]
            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='TCR_beta',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=0.95,
                residue_range=chain_data['residue_range'],
                naming_convention='segname',
                anarci_result=result
            )
            logger.info(f"[ANARCI] {segname} = TCR-beta")
        elif len(tcr_beta_candidates) > 1:
            logger.warning(f"Multiple TCR-beta candidates: {[c[0] for c in tcr_beta_candidates]}")
            # Use longest as beta
            tcr_beta_candidates.sort(key=lambda x: x[1]['length'], reverse=True)
            segname, chain_data, result = tcr_beta_candidates[0]
            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='TCR_beta',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=0.7,
                residue_range=chain_data['residue_range'],
                naming_convention='segname',
                anarci_result=result
            )
            logger.info(f"[ANARCI] {segname} = TCR-beta (multiple candidates, using longest)")

        # Assign HLA-alpha (remaining chain)
        remaining = [
            (sid, cdata) for sid, cdata in chains.items()
            if sid not in identifications
        ]

        if len(remaining) == 1:
            segname, chain_data = remaining[0]

            # Validate HLA-alpha length
            confidence = self._validate_hla_alpha_length(chain_data['length'], segname)

            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='HLA_alpha',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=confidence,
                residue_range=chain_data['residue_range'],
                naming_convention='segname'
            )
            logger.info(
                f"[ELIMINATION] {segname} = HLA-alpha "
                f"(remaining after TCR identification, confidence={confidence:.2f})"
            )
        elif len(remaining) > 1:
            # Multiple remaining - assign longest as HLA-alpha
            remaining.sort(key=lambda x: x[1]['length'], reverse=True)
            segname, chain_data = remaining[0]

            confidence = self._validate_hla_alpha_length(chain_data['length'], segname)
            confidence = min(confidence, 0.70)

            identifications[segname] = UnifiedChainInfo(
                chain_id=segname,
                chain_type='HLA_alpha',
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                confidence=confidence,
                residue_range=chain_data['residue_range'],
                naming_convention='segname'
            )
            logger.warning(
                f"Multiple unidentified chains, assigned longest {segname} as HLA-alpha "
                f"(length={chain_data['length']} AA, confidence={confidence:.2f})"
            )

            # Mark others as unknown
            for segname, chain_data in remaining[1:]:
                identifications[segname] = UnifiedChainInfo(
                    chain_id=segname,
                    chain_type='unknown',
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    confidence=0.3,
                    residue_range=chain_data['residue_range'],
                    naming_convention='segname'
                )

        return identifications

    def _validate_hla_alpha_length(self, length: int, chain_id: str) -> float:
        """
        Validate HLA-alpha length and return confidence score.

        Typical: 250-290 AA (high confidence)
        Extended: 230-310 AA (medium confidence)
        Outside: low confidence

        Parameters
        ----------
        length : int
            Chain length in amino acids
        chain_id : str
            Chain ID for logging

        Returns
        -------
        confidence : float
            Confidence score (0.0-1.0)
        """
        if 250 <= length <= 290:
            logger.info(
                f"  {chain_id} length {length} AA is within typical HLA-alpha range (250-290 AA)"
            )
            return 0.90
        elif 230 <= length < 250:
            logger.warning(
                f"  {chain_id} length {length} AA is slightly shorter than typical HLA-alpha (250-290 AA)"
            )
            return 0.70
        elif 290 < length <= 310:
            logger.warning(
                f"  {chain_id} length {length} AA is slightly longer than typical HLA-alpha (250-290 AA)"
            )
            return 0.70
        else:
            logger.error(
                f"  {chain_id} length {length} AA is VERY UNUSUAL for HLA-alpha! "
                f"Expected range: 230-310 AA"
            )
            return 0.30

    def _identify_by_length_fallback(
        self,
        chains: Dict[str, Dict]
    ) -> Dict[str, UnifiedChainInfo]:
        """
        Fallback identification based on length alone (without ANARCI).

        Assumes:
        - Longest chain: HLA-alpha
        - Middle 2 chains: TCR (shorter=alpha, longer=beta)

        Parameters
        ----------
        chains : Dict[str, Dict]
            Dictionary of 3 chains

        Returns
        -------
        identifications : Dict[str, UnifiedChainInfo]
            Dictionary of chain identifications
        """
        identifications = {}
        sorted_chains = sorted(chains.items(), key=lambda x: x[1]['length'])

        if len(sorted_chains) != 3:
            logger.error(f"Expected 3 chains, got {len(sorted_chains)}")
            return identifications

        # Shortest: TCR-alpha (usually)
        segname, chain_data = sorted_chains[0]
        identifications[segname] = UnifiedChainInfo(
            chain_id=segname,
            chain_type='TCR_alpha',
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            confidence=0.6,
            residue_range=chain_data['residue_range'],
            naming_convention='segname'
        )
        logger.info(f"Identified {segname} as TCR-alpha (length fallback)")

        # Middle: TCR-beta (usually)
        segname, chain_data = sorted_chains[1]
        identifications[segname] = UnifiedChainInfo(
            chain_id=segname,
            chain_type='TCR_beta',
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            confidence=0.6,
            residue_range=chain_data['residue_range'],
            naming_convention='segname'
        )
        logger.info(f"Identified {segname} as TCR-beta (length fallback)")

        # Longest: HLA-alpha
        segname, chain_data = sorted_chains[2]
        confidence = self._validate_hla_alpha_length(chain_data['length'], segname)
        confidence = min(confidence, 0.70)

        identifications[segname] = UnifiedChainInfo(
            chain_id=segname,
            chain_type='HLA_alpha',
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            confidence=confidence,
            residue_range=chain_data['residue_range'],
            naming_convention='segname'
        )
        logger.info(
            f"Identified {segname} as HLA-alpha "
            f"(length fallback, {chain_data['length']} AA, confidence={confidence:.2f})"
        )

        return identifications
