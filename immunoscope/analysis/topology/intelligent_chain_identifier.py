"""
Intelligent Chain Identifier for pHLA-TCR Complexes

Uses ANARCI to identify TCR alpha/beta chains, combined with length-based
heuristics for peptide and beta2-microglobulin identification.

Author: ImmunoScope Development Team
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from ..structure import PDBSequenceExtractor
from .cdr_manager import ANARCIWrapper

logger = logging.getLogger(__name__)


@dataclass
class ChainIdentification:
    """Result of chain identification."""
    chain_id: str
    length: int
    sequence: str
    chain_type: str  # 'peptide', 'beta2m', 'TCR_alpha', 'TCR_beta', 'HLA_alpha'
    confidence: float  # 0.0 to 1.0
    anarci_result: Optional[Dict] = None


class IntelligentChainIdentifier:
    """
    Intelligent chain identifier using ANARCI and length heuristics.

    NEW IDENTIFICATION STRATEGY (2026-01-22):
    1. Identify peptide by length (<=20 AA) - DEFINITIVE
    2. Identify beta2-microglobulin by length (90-110 AA) - DEFINITIVE
    3. Send remaining 3 chains to ANARCI:
       - ANARCI identifies TCR-alpha and TCR-beta
       - Remaining chain is HLA-alpha

    Standard chain assignment:
    - Chain A: HLA-alpha (identified by elimination)
    - Chain B: beta2-microglobulin (90-110 AA)
    - Chain C: Peptide (<=20 AA)
    - Chain D: TCR-alpha (ANARCI identification)
    - Chain E: TCR-beta (ANARCI identification)
    """

    # Length thresholds for identification (NEW STRATEGY)
    PEPTIDE_MAX_LENGTH = 20      # Peptide: <=20 AA (definitive)
    BETA2M_MIN_LENGTH = 90       # Beta2m: ~100 AA (90-110 range)
    BETA2M_MAX_LENGTH = 110

    # HLA-alpha length validation ranges
    HLA_ALPHA_TYPICAL_MIN = 250  # Typical range: 250-290 AA
    HLA_ALPHA_TYPICAL_MAX = 290
    HLA_ALPHA_EXTENDED_MIN = 230  # Extended range: 230-310 AA (with warning)
    HLA_ALPHA_EXTENDED_MAX = 310

    def __init__(self, use_anarci: bool = True):
        """
        Initialize intelligent chain identifier.

        Args:
            use_anarci: Whether to use ANARCI for TCR identification
        """
        self.use_anarci = use_anarci
        self.sequence_extractor = PDBSequenceExtractor()

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
        pdb_file: str
    ) -> Dict[str, ChainIdentification]:
        """
        Identify all chains in PDB file using intelligent strategy.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary mapping original chain_id to ChainIdentification

        Example:
            >>> identifier = IntelligentChainIdentifier()
            >>> results = identifier.identify_chains("1ao7.pdb")
            >>> for chain_id, info in results.items():
            ...     print(f"{chain_id}: {info.chain_type} ({info.length} aa)")
        """
        # Step 1: Extract sequences
        chains_data = self.sequence_extractor.extract_sequences_from_pdb(pdb_file)

        if len(chains_data) != 5:
            logger.warning(
                f"Expected 5 chains in pHLA-TCR complex, found {len(chains_data)}"
            )

        # Step 2: Sort by length
        sorted_chains = sorted(
            chains_data.items(),
            key=lambda x: x[1]['length']
        )

        identifications = {}

        # Step 3: Identify peptide (<=20 AA, definitive)
        peptide_chain_id, peptide_data = sorted_chains[0]
        if peptide_data['length'] <= self.PEPTIDE_MAX_LENGTH:
            identifications[peptide_chain_id] = ChainIdentification(
                chain_id=peptide_chain_id,
                length=peptide_data['length'],
                sequence=peptide_data['sequence'],
                chain_type='peptide',
                confidence=1.0
            )
            logger.info(f"[DEFINITIVE] {peptide_chain_id} = peptide ({peptide_data['length']} AA, <=20)")
        else:
            logger.warning(
                f"Shortest chain {peptide_chain_id} has {peptide_data['length']} AA, "
                f"exceeds peptide threshold (<=20 AA)"
            )
            identifications[peptide_chain_id] = ChainIdentification(
                chain_id=peptide_chain_id,
                length=peptide_data['length'],
                sequence=peptide_data['sequence'],
                chain_type='unknown',
                confidence=0.5
            )

        # Step 4: Identify beta2-microglobulin (90-110 AA, definitive)
        beta2m_found = False
        for chain_id, chain_data in sorted_chains[1:]:
            if self.BETA2M_MIN_LENGTH <= chain_data['length'] <= self.BETA2M_MAX_LENGTH:
                identifications[chain_id] = ChainIdentification(
                    chain_id=chain_id,
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    chain_type='beta2m',
                    confidence=1.0
                )
                logger.info(f"[DEFINITIVE] {chain_id} = beta2m ({chain_data['length']} AA, 90-110 range)")
                beta2m_found = True
                break

        if not beta2m_found:
            logger.warning("Could not identify beta2m chain (expected 90-110 AA range)")

        # Step 5: Identify remaining chains using ANARCI
        remaining_chains = {
            cid: cdata for cid, cdata in chains_data.items()
            if cid not in identifications
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
            for chain_id, chain_data in remaining_chains.items():
                identifications[chain_id] = ChainIdentification(
                    chain_id=chain_id,
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    chain_type='unknown',
                    confidence=0.0
                )

        return identifications

    def _validate_hla_alpha_length(self, length: int, chain_id: str) -> float:
        """
        Validate HLA-alpha length and return appropriate confidence score.

        Typical HLA-alpha length: 270-280 AA
        Acceptable range: 250-290 AA (high confidence)
        Extended range: 230-310 AA (medium confidence with warning)
        Outside range: low confidence with error

        Args:
            length: Chain length in amino acids
            chain_id: Chain ID for logging

        Returns:
            Confidence score (0.0-1.0)
        """
        if self.HLA_ALPHA_TYPICAL_MIN <= length <= self.HLA_ALPHA_TYPICAL_MAX:
            logger.info(
                f"  {chain_id} length {length} AA is within typical HLA-alpha range "
                f"({self.HLA_ALPHA_TYPICAL_MIN}-{self.HLA_ALPHA_TYPICAL_MAX} AA)"
            )
            return 0.90
        elif self.HLA_ALPHA_EXTENDED_MIN <= length < self.HLA_ALPHA_TYPICAL_MIN:
            logger.warning(
                f"  {chain_id} length {length} AA is slightly shorter than typical HLA-alpha "
                f"(typical: {self.HLA_ALPHA_TYPICAL_MIN}-{self.HLA_ALPHA_TYPICAL_MAX} AA)"
            )
            return 0.70
        elif self.HLA_ALPHA_TYPICAL_MAX < length <= self.HLA_ALPHA_EXTENDED_MAX:
            logger.warning(
                f"  {chain_id} length {length} AA is slightly longer than typical HLA-alpha "
                f"(typical: {self.HLA_ALPHA_TYPICAL_MIN}-{self.HLA_ALPHA_TYPICAL_MAX} AA)"
            )
            return 0.70
        else:
            logger.error(
                f"  {chain_id} length {length} AA is VERY UNUSUAL for HLA-alpha! "
                f"Expected range: {self.HLA_ALPHA_EXTENDED_MIN}-{self.HLA_ALPHA_EXTENDED_MAX} AA. "
                f"This may indicate: (1) truncated HLA-alpha, (2) misidentified chain, "
                f"or (3) non-standard pHLA-TCR structure."
            )
            return 0.30

    def _identify_tcr_and_hla(
        self,
        chains: Dict[str, Dict]
    ) -> Dict[str, ChainIdentification]:
        """
        Identify TCR alpha, TCR beta, and HLA-alpha from 3 remaining chains.

        NEW STRATEGY (2026-01-22):
        1. Send all 3 chains to ANARCI for TCR identification
        2. ANARCI should identify 2 TCR chains (alpha and beta)
        3. The remaining chain (not identified as TCR) is HLA-alpha

        Args:
            chains: Dictionary of 3 remaining chains (after peptide and beta2m removed)

        Returns:
            Dictionary of ChainIdentification results
        """
        identifications = {}

        if not self.use_anarci or self.anarci_wrapper is None:
            # Fallback: length-based identification
            return self._identify_by_length_fallback(chains)

        # Run ANARCI on all 3 chains
        logger.info("Running ANARCI on 3 remaining chains to identify TCR...")
        anarci_results = {}
        for chain_id, chain_data in chains.items():
            try:
                result = self.anarci_wrapper.run_anarci(
                    sequence=chain_data['sequence'],
                    chain_type='TCR'
                )
                anarci_results[chain_id] = result
                chain_type_detected = result.get('chain_type', 'unknown')
                logger.info(
                    f"  {chain_id} ({chain_data['length']} AA): "
                    f"ANARCI result = {chain_type_detected}"
                )
            except Exception as e:
                logger.warning(f"  {chain_id}: ANARCI failed ({e})")
                anarci_results[chain_id] = None

        # Identify chains based on ANARCI results
        tcr_alpha_candidates = []
        tcr_beta_candidates = []
        hla_alpha_candidates = []

        for chain_id, result in anarci_results.items():
            chain_data = chains[chain_id]

            if result and 'chain_type' in result:
                chain_type_anarci = result['chain_type'].lower()

                # TCR alpha-like: TCR_alpha or TCR_delta
                if ('alpha' in chain_type_anarci or 'delta' in chain_type_anarci) and 'tcr' in chain_type_anarci:
                    tcr_alpha_candidates.append((chain_id, chain_data, result))
                # TCR beta-like: TCR_beta or TCR_gamma
                elif ('beta' in chain_type_anarci or 'gamma' in chain_type_anarci) and 'tcr' in chain_type_anarci:
                    tcr_beta_candidates.append((chain_id, chain_data, result))
                else:
                    # Likely HLA-alpha (ANARCI doesn't recognize it as TCR)
                    hla_alpha_candidates.append((chain_id, chain_data, result))
            else:
                # No TCR recognition - likely HLA-alpha
                hla_alpha_candidates.append((chain_id, chain_data, None))

        # Assign TCR-alpha
        if len(tcr_alpha_candidates) == 1:
            chain_id, chain_data, result = tcr_alpha_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_alpha',
                confidence=0.95,
                anarci_result=result
            )
            logger.info(f"[ANARCI] {chain_id} = TCR-alpha")
        elif len(tcr_alpha_candidates) > 1:
            logger.warning(f"Multiple TCR-alpha candidates: {[c[0] for c in tcr_alpha_candidates]}")
            # Use shortest as alpha
            tcr_alpha_candidates.sort(key=lambda x: x[1]['length'])
            chain_id, chain_data, result = tcr_alpha_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_alpha',
                confidence=0.7,
                anarci_result=result
            )
            logger.info(f"[ANARCI] {chain_id} = TCR-alpha (multiple candidates, using shortest)")

        # Assign TCR-beta
        if len(tcr_beta_candidates) == 1:
            chain_id, chain_data, result = tcr_beta_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_beta',
                confidence=0.95,
                anarci_result=result
            )
            logger.info(f"[ANARCI] {chain_id} = TCR-beta")
        elif len(tcr_beta_candidates) > 1:
            logger.warning(f"Multiple TCR-beta candidates: {[c[0] for c in tcr_beta_candidates]}")
            # Use longest as beta
            tcr_beta_candidates.sort(key=lambda x: x[1]['length'], reverse=True)
            chain_id, chain_data, result = tcr_beta_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_beta',
                confidence=0.7,
                anarci_result=result
            )
            logger.info(f"[ANARCI] {chain_id} = TCR-beta (multiple candidates, using longest)")

        # Assign HLA-alpha (remaining chain after TCR identification)
        remaining = [
            (cid, cdata) for cid, cdata in chains.items()
            if cid not in identifications
        ]

        if len(remaining) == 1:
            chain_id, chain_data = remaining[0]

            # Validate HLA-alpha length and get appropriate confidence
            confidence = self._validate_hla_alpha_length(chain_data['length'], chain_id)

            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='HLA_alpha',
                confidence=confidence
            )
            logger.info(
                f"[ELIMINATION] {chain_id} = HLA-alpha "
                f"(remaining after TCR identification, confidence={confidence:.2f})"
            )
        elif len(remaining) > 1:
            # Multiple remaining - assign longest as HLA-alpha
            remaining.sort(key=lambda x: x[1]['length'], reverse=True)
            chain_id, chain_data = remaining[0]

            # Validate HLA-alpha length
            confidence = self._validate_hla_alpha_length(chain_data['length'], chain_id)
            # Lower confidence due to ambiguity
            confidence = min(confidence, 0.70)

            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='HLA_alpha',
                confidence=confidence
            )
            logger.warning(
                f"Multiple unidentified chains, assigned longest {chain_id} as HLA-alpha "
                f"(length={chain_data['length']} AA, confidence={confidence:.2f})"
            )

            # Mark others as unknown
            for chain_id, chain_data in remaining[1:]:
                identifications[chain_id] = ChainIdentification(
                    chain_id=chain_id,
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    chain_type='unknown',
                    confidence=0.3
                )

        return identifications

    def _identify_by_length_fallback(
        self,
        chains: Dict[str, Dict]
    ) -> Dict[str, ChainIdentification]:
        """
        Fallback identification based on length alone (without ANARCI).

        Assumes:
        - Longest chain: HLA-alpha
        - Middle 2 chains: TCR (shorter=alpha, longer=beta)

        Args:
            chains: Dictionary of 3 chains

        Returns:
            Dictionary of ChainIdentification results
        """
        identifications = {}
        sorted_chains = sorted(chains.items(), key=lambda x: x[1]['length'])

        if len(sorted_chains) != 3:
            logger.error(f"Expected 3 chains, got {len(sorted_chains)}")
            return identifications

        # Shortest of the 3: TCR-alpha (usually)
        chain_id, chain_data = sorted_chains[0]
        identifications[chain_id] = ChainIdentification(
            chain_id=chain_id,
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            chain_type='TCR_alpha',
            confidence=0.6  # Lower confidence without ANARCI
        )
        logger.info(f"Identified {chain_id} as TCR-alpha (length fallback)")

        # Middle: TCR-beta (usually)
        chain_id, chain_data = sorted_chains[1]
        identifications[chain_id] = ChainIdentification(
            chain_id=chain_id,
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            chain_type='TCR_beta',
            confidence=0.6
        )
        logger.info(f"Identified {chain_id} as TCR-beta (length fallback)")

        # Longest: HLA-alpha (assumption - may not always be true)
        chain_id, chain_data = sorted_chains[2]

        # Validate HLA-alpha length
        confidence = self._validate_hla_alpha_length(chain_data['length'], chain_id)
        # Lower confidence due to length-only fallback
        confidence = min(confidence, 0.70)

        identifications[chain_id] = ChainIdentification(
            chain_id=chain_id,
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            chain_type='HLA_alpha',
            confidence=confidence
        )
        logger.info(
            f"Identified {chain_id} as HLA-alpha "
            f"(length fallback, {chain_data['length']} AA, confidence={confidence:.2f})"
        )

        return identifications

    def create_standardization_mapping(
        self,
        identifications: Dict[str, ChainIdentification]
    ) -> Dict[str, str]:
        """
        Create chain ID mapping for standardization.

        Standard ordering:
        - Chain A: HLA-alpha
        - Chain B: beta2m
        - Chain C: Peptide
        - Chain D: TCR-alpha
        - Chain E: TCR-beta

        Args:
            identifications: Dictionary from identify_chains()

        Returns:
            Dictionary mapping old_chain_id -> new_chain_id

        Example:
            >>> mapping = identifier.create_standardization_mapping(identifications)
            >>> print(mapping)  # {'X': 'C', 'Y': 'B', 'Z': 'D', 'W': 'E', 'V': 'A'}
        """
        mapping = {}

        for chain_id, info in identifications.items():
            if info.chain_type == 'peptide':
                mapping[chain_id] = 'C'
            elif info.chain_type == 'beta2m':
                mapping[chain_id] = 'B'
            elif info.chain_type == 'TCR_alpha':
                mapping[chain_id] = 'D'
            elif info.chain_type == 'TCR_beta':
                mapping[chain_id] = 'E'
            elif info.chain_type == 'HLA_alpha':
                mapping[chain_id] = 'A'
            else:
                logger.warning(f"Unknown chain type for {chain_id}: {info.chain_type}")
                mapping[chain_id] = chain_id  # Keep original

        return mapping
