"""
GROMACS Index Manager - Unified index file management

This module provides centralized, enumerable, and traceable index file management
for all analysis modules. It eliminates fragmented index generation and provides
automatic group name to group ID mapping.

Key Features:
- Enum-based component definitions (StandardComponent)
- Automatic group registry (name -> ID mapping via .ndx parsing)
- Lazy loading with caching for static and dynamic components
- Thread-safe through instance isolation
- Seamless integration with PipelineContext

Author: ImmunoScope Development Team
Date: 2026-03-17
"""

import logging
import re
import threading
from pathlib import Path
from typing import Dict, Optional, List
from enum import Enum
from dataclasses import dataclass

from immunoscope.analysis.topology import ComponentIndexGenerator

logger = logging.getLogger(__name__)


class StandardComponent(Enum):
    """
    Standard component definitions for pHLA-TCR complexes.

    Components are divided into two categories:

    1. Static Components (chain-based, always available):
       - HLA, HLA_ALPHA, HLA_BETA
       - PHLA, PEPTIDE
       - TCR, TCR_ALPHA, TCR_BETA

    2. Dynamic Components (require additional information):
       - CDR loops: CDR1_ALPHA, CDR2_ALPHA, CDR3_ALPHA,
                    CDR1_BETA, CDR2_BETA, CDR3_BETA
       - Interface: INTERFACE (TCR-pMHC contact region)

    Chain standardization:
    - Chain A: HLA-alpha (~270-280 residues)
    - Chain B: HLA-beta/beta2m (~99-100 residues)
    - Chain C: Peptide (5-22 residues)
    - Chain D: TCR-alpha (115-215 residues)
    - Chain E: TCR-beta (209-252 residues)

    Note: Dynamic components use tuple values to avoid Enum aliasing
    """

    # Static components (chain-based)
    HLA = ['A', 'B']
    HLA_ALPHA = ['A']
    HLA_BETA = ['B']
    PHLA = ['A', 'B', 'C']
    PEPTIDE = ['C']
    TCR = ['D', 'E']
    TCR_ALPHA = ['D']
    TCR_BETA = ['E']

    # Dynamic components - CDR loops (use tuples to avoid aliasing)
    CDR1_ALPHA = (['D'], 'CDR1')  # TCR alpha CDR1 loop
    CDR2_ALPHA = (['D'], 'CDR2')  # TCR alpha CDR2 loop
    CDR3_ALPHA = (['D'], 'CDR3')  # TCR alpha CDR3 loop
    CDR1_BETA = (['E'], 'CDR1')   # TCR beta CDR1 loop
    CDR2_BETA = (['E'], 'CDR2')   # TCR beta CDR2 loop
    CDR3_BETA = (['E'], 'CDR3')   # TCR beta CDR3 loop

    # Dynamic components - Interface (require distance calculation)
    INTERFACE = ['A', 'B', 'C', 'D', 'E']  # TCR-pMHC interface residues

    @property
    def chains(self) -> List[str]:
        """Get chain IDs for this component."""
        value = self.value
        # Handle both list and tuple formats
        if isinstance(value, tuple):
            return value[0]
        return value

    @property
    def component_name(self) -> str:
        """Get standardized component name (e.g., 'pHLA', 'TCR', 'CDR3_alpha')."""
        name = self.name
        if name == 'PHLA':
            return 'pHLA'
        elif '_' in name:
            # HLA_ALPHA -> HLA_alpha
            # CDR3_ALPHA -> CDR3_alpha
            parts = name.split('_')
            if len(parts) == 2:
                return parts[0] + '_' + parts[1].lower()
            else:
                # Handle cases like CDR_LOOPS_ALPHA if needed
                return '_'.join([parts[0]] + [p.lower() for p in parts[1:]])
        else:
            return name

    @property
    def is_dynamic(self) -> bool:
        """Check if this component requires dynamic generation."""
        name = self.name
        return (name.startswith('CDR') or
                name == 'INTERFACE')

    @property
    def requires_sequence(self) -> bool:
        """Check if this component requires sequence information."""
        return self.name.startswith('CDR')

    @property
    def requires_distance_calculation(self) -> bool:
        """Check if this component requires distance-based calculation."""
        return self.name == 'INTERFACE'


@dataclass
class GroupInfo:
    """Information about a group in the index file."""
    group_id: int          # Group number in .ndx file (0-based in file, displayed as is)
    group_name: str        # Group name (e.g., 'pHLA', 'TCR')
    atom_count: int        # Number of atoms in this group
    component: Optional[StandardComponent] = None  # Associated component enum


class IndexManager:
    """
    GROMACS Index Manager - Unified index file management.

    This class provides centralized index generation and group ID tracking
    for all analysis modules. It ensures that:
    1. Index files are generated once and cached
    2. Group names are automatically mapped to group IDs
    3. Analysis modules can request group IDs by component name
    4. GROMACS commands receive group IDs directly (e.g., "0\n4\n")

    Usage:
        # Initialize with context
        index_mgr = IndexManager(context)

        # Ensure base index exists (static components)
        base_index = index_mgr.ensure_base_index()

        # Get group ID for a component
        phla_group_id = index_mgr.get_group_id('pHLA')  # Returns: "20"
        tcr_group_id = index_mgr.get_group_id('TCR')    # Returns: "21"

        # Use in GROMACS command
        stdin_input = f"{phla_group_id}\n{tcr_group_id}\n"

    Thread Safety:
        Each task gets its own IndexManager instance through PipelineContext,
        ensuring thread isolation. Internal locks protect shared operations.
    """

    def __init__(self, context: 'PipelineContext'):
        """
        Initialize IndexManager.

        Args:
            context: PipelineContext instance for this task
        """
        self.context = context
        self.base_index_file: Optional[Path] = None
        self.dynamic_index_files: Dict[str, Path] = {}  # component_key -> index_file
        self.group_registry: Dict[str, GroupInfo] = {}  # group_name -> GroupInfo
        self._lock = threading.Lock()

        # Initialize IndexGenerator (will be created when needed)
        self._index_generator: Optional[ComponentIndexGenerator] = None

        # Chain standardization status
        self._chain_standardized: Optional[bool] = None
        self._chain_validation_result: Optional[Dict] = None

        logger.debug(f"IndexManager initialized for {context.system_id}")

    @property
    def index_generator(self) -> ComponentIndexGenerator:
        """Lazy initialization of the component index generator."""
        if self._index_generator is None:
            # Use structure_pdb if available, otherwise topology
            topology_file = self.context.structure_pdb or self.context.topology
            self._index_generator = ComponentIndexGenerator(
                topology=topology_file,
                auto_standardize=True
            )

            # Verify chain standardization
            self._verify_chain_standardization()

        return self._index_generator

    def _verify_chain_standardization(self) -> None:
        """
        Verify that chains are properly standardized.

        This method checks:
        1. Whether chain standardization was performed
        2. Whether the standardized topology exists
        3. Chain count and ordering correctness

        Sets:
            self._chain_standardized: True if standardization successful
            self._chain_validation_result: Detailed validation results
        """
        try:
            import MDAnalysis as mda

            # Check if standardized PDB exists
            if self._index_generator.standardized_pdb is None:
                logger.warning("⚠️  Chain standardization was not performed!")
                logger.warning("    Using original topology - chain IDs may not be standardized")
                self._chain_standardized = False
                self._chain_validation_result = {
                    'standardized': False,
                    'reason': 'No standardized PDB generated',
                    'recommendation': 'Provide a PDB file with correct chain ordering'
                }
                return

            # Load standardized structure
            std_pdb = str(self._index_generator.standardized_pdb)
            u = mda.Universe(std_pdb)

            # Get protein chains
            protein_chains = {}
            for chain_id in set(u.select_atoms('protein').chainIDs):
                chain_atoms = u.select_atoms(f'protein and chainID {chain_id}')
                residue_count = len(chain_atoms.residues)
                protein_chains[chain_id] = residue_count

            # Validate chain count
            if len(protein_chains) != 5:
                logger.warning(f"⚠️  Expected 5 chains, found {len(protein_chains)}")
                logger.warning(f"    Chains: {protein_chains}")
                self._chain_standardized = False
                self._chain_validation_result = {
                    'standardized': True,
                    'valid': False,
                    'reason': f'Incorrect chain count: {len(protein_chains)} (expected 5)',
                    'chains': protein_chains
                }
                return

            # Expected chain ordering (by residue count)
            expected_order = {
                'C': (5, 22),      # Peptide: shortest
                'B': (90, 110),    # β2m: ~100 residues
                'D': (100, 220),   # TCR-α: variable
                'E': (200, 260),   # TCR-β: longer
                'A': (260, 290)    # HLA-α: longest
            }

            # Validate each chain
            validation_issues = []
            for chain_id, (min_res, max_res) in expected_order.items():
                if chain_id not in protein_chains:
                    validation_issues.append(f"Chain {chain_id} missing")
                    continue

                res_count = protein_chains[chain_id]
                if not (min_res <= res_count <= max_res):
                    validation_issues.append(
                        f"Chain {chain_id}: {res_count} residues (expected {min_res}-{max_res})"
                    )

            if validation_issues:
                logger.warning("⚠️  Chain validation warnings:")
                for issue in validation_issues:
                    logger.warning(f"    - {issue}")
                logger.warning("    Chain ordering may not match pHLA-TCR standard")

                self._chain_standardized = True
                self._chain_validation_result = {
                    'standardized': True,
                    'valid': False,
                    'reason': 'Chain residue counts outside expected ranges',
                    'chains': protein_chains,
                    'issues': validation_issues
                }
            else:
                logger.info("✅ Chain standardization verified:")
                logger.info(f"    Chain A (HLA-α):  {protein_chains['A']} residues")
                logger.info(f"    Chain B (β2m):    {protein_chains['B']} residues")
                logger.info(f"    Chain C (peptide): {protein_chains['C']} residues")
                logger.info(f"    Chain D (TCR-α):  {protein_chains['D']} residues")
                logger.info(f"    Chain E (TCR-β):  {protein_chains['E']} residues")

                self._chain_standardized = True
                self._chain_validation_result = {
                    'standardized': True,
                    'valid': True,
                    'chains': protein_chains
                }

        except ImportError:
            logger.warning("MDAnalysis not available - cannot verify chain standardization")
            self._chain_standardized = None
            self._chain_validation_result = {
                'standardized': None,
                'reason': 'MDAnalysis not available'
            }
        except Exception as e:
            logger.error(f"Chain validation failed: {e}")
            self._chain_standardized = None
            self._chain_validation_result = {
                'standardized': None,
                'reason': f'Validation error: {str(e)}'
            }

    def get_chain_validation_report(self) -> Dict:
        """
        Get chain standardization validation report.

        Returns:
            Dictionary with validation results

        Example:
            >>> report = index_mgr.get_chain_validation_report()
            >>> print(report['standardized'])  # True/False/None
            >>> print(report['chains'])        # {'A': 275, 'B': 99, ...}
        """
        if self._chain_validation_result is None:
            # Trigger validation by accessing index_generator
            _ = self.index_generator

        return self._chain_validation_result

    def is_chain_standardized(self) -> bool:
        """
        Check if chains are properly standardized and validated.

        Returns:
            True if chains are standardized and validated, False otherwise

        Example:
            >>> if not index_mgr.is_chain_standardized():
            ...     print("Warning: Chains may not be standardized!")
        """
        if self._chain_standardized is None:
            # Trigger validation
            _ = self.index_generator

        return self._chain_standardized is True and \
               (self._chain_validation_result or {}).get('valid', False)

    def ensure_base_index(self, force_regenerate: bool = False) -> Path:
        """
        Ensure base index file exists (contains all static components).

        The base index includes all StandardComponent enums:
        - HLA, HLA_alpha, HLA_beta
        - pHLA, peptide
        - TCR, TCR_alpha, TCR_beta

        Args:
            force_regenerate: If True, regenerate even if cached

        Returns:
            Path to base index file

        Raises:
            RuntimeError: If chain standardization failed and validation is critical

        Example:
            >>> index_mgr = IndexManager(context)
            >>> base_index = index_mgr.ensure_base_index()
            >>> print(base_index)
            ./results/1ao7/indices/base_components.ndx
        """
        with self._lock:
            # Check chain standardization status before generating index
            if not self.is_chain_standardized():
                report = self.get_chain_validation_report()

                if report.get('standardized') is False:
                    logger.error("❌ Chain standardization failed!")
                    logger.error(f"   Reason: {report.get('reason', 'Unknown')}")
                    logger.error("   ⚠️  Generated indices may use incorrect chain IDs!")
                    logger.error("   ⚠️  Analysis results will be INCORRECT!")

                    # Add error to context
                    self.context.add_error(
                        f"Chain standardization failed: {report.get('reason')}"
                    )

                    # Add warning to metadata
                    self.context.metadata['chain_standardization'] = 'FAILED'

                elif report.get('valid') is False:
                    logger.warning("⚠️  Chain validation warnings detected:")
                    if 'issues' in report:
                        for issue in report['issues']:
                            logger.warning(f"    - {issue}")
                    logger.warning("   Proceeding with caution - verify results manually!")

                    self.context.metadata['chain_standardization'] = 'WARNING'

            if self.base_index_file is not None and not force_regenerate:
                if self.base_index_file.exists():
                    logger.debug(f"Using cached base index: {self.base_index_file}")
                    return self.base_index_file

            # Generate base index with all standard components
            output_file = self.context.get_index_path("base_components.ndx")
            output_path = Path(output_file)

            logger.info(f"Generating base index with {len(StandardComponent)} components")

            # Get all static component names (exclude dynamic components)
            static_components = [
                comp.component_name for comp in StandardComponent
                if not comp.is_dynamic
            ]

            logger.info(f"Static components: {', '.join(static_components)}")

            # Use IndexGenerator to create multi-component index
            self.index_generator.generate_multi_component_index(
                components=static_components,
                output_file=str(output_path)
            )

            self.base_index_file = output_path

            # Parse the index file to build group registry
            self._parse_index_file(output_path)

            logger.info(f"Base index created with {len(self.group_registry)} groups")

            return self.base_index_file

    def get_group_id(self, component_name: str) -> Optional[str]:
        """
        Get group ID for a component by name.

        Args:
            component_name: Component name (e.g., 'pHLA', 'TCR', 'HLA_alpha')

        Returns:
            Group ID as string (e.g., "20"), or None if not found

        Example:
            >>> phla_id = index_mgr.get_group_id('pHLA')
            >>> print(phla_id)
            20
            >>> tcr_id = index_mgr.get_group_id('TCR')
            >>> print(tcr_id)
            21
        """
        # Ensure base index exists
        if self.base_index_file is None:
            self.ensure_base_index()

        group_info = self.group_registry.get(component_name)

        if group_info is None:
            logger.warning(f"Component '{component_name}' not found in group registry")
            return None

        return str(group_info.group_id)

    def get_group_info(self, component_name: str) -> Optional[GroupInfo]:
        """
        Get complete group information.

        Args:
            component_name: Component name

        Returns:
            GroupInfo instance, or None if not found
        """
        if self.base_index_file is None:
            self.ensure_base_index()

        return self.group_registry.get(component_name)

    def list_available_components(self) -> List[str]:
        """
        List all available component names in the registry.

        Returns:
            List of component names

        Example:
            >>> components = index_mgr.list_available_components()
            >>> print(components)
            ['HLA', 'HLA_alpha', 'HLA_beta', 'pHLA', 'peptide', 'TCR', 'TCR_alpha', 'TCR_beta']
        """
        if self.base_index_file is None:
            self.ensure_base_index()

        return sorted(self.group_registry.keys())

    def ensure_cdr3_index(self,
                         chain: str,
                         sequence: str,
                         ca_only: bool = True) -> Optional[str]:
        """
        Ensure CDR3 index exists for a specific chain and sequence.

        Args:
            chain: TCR chain ('alpha' or 'beta')
            sequence: CDR3 amino acid sequence
            ca_only: If True, select only CA atoms (default: True)

        Returns:
            Group ID as string, or None if sequence not found

        Example:
            >>> cdr3_id = index_mgr.ensure_cdr3_index('beta', 'CASSLGQAYEQYF')
            >>> print(cdr3_id)
            25
        """
        if chain.lower() not in ['alpha', 'beta']:
            raise ValueError("Chain must be 'alpha' or 'beta'")

        component_key = f"CDR3_{chain.lower()}"

        with self._lock:
            # Check if already generated
            if component_key in self.dynamic_index_files:
                existing_file = self.dynamic_index_files[component_key]
                if existing_file.exists():
                    logger.debug(f"Using cached CDR3 index: {existing_file}")
                    # Return group ID from registry
                    group_info = self.group_registry.get(component_key)
                    return str(group_info.group_id) if group_info else None

            # Generate CDR3 index
            output_file = self.context.get_index_path(f"cdr3_{chain.lower()}.ndx")

            logger.info(f"Generating CDR3 index for {chain} chain")
            logger.info(f"  Sequence: {sequence}")

            result_file = self.index_generator.generate_cdr3_index(
                chain=chain,
                cdr3_sequence=sequence,
                output_file=output_file,
                ca_only=ca_only
            )

            if result_file is None:
                logger.warning(f"CDR3 sequence not found in {chain} chain")
                return None

            result_path = Path(result_file)
            self.dynamic_index_files[component_key] = result_path

            # Parse the CDR3 index to get group info
            self._parse_cdr3_index(result_path, component_key)

            group_info = self.group_registry.get(component_key)
            if group_info:
                logger.info(f"CDR3 index created: group {group_info.group_id} with {group_info.atom_count} atoms")
                return str(group_info.group_id)
            else:
                return None

    def ensure_cdr_loops_index(self,
                              chain: str,
                              cdr1_seq: Optional[str] = None,
                              cdr2_seq: Optional[str] = None,
                              cdr3_seq: Optional[str] = None,
                              ca_only: bool = True) -> Dict[str, Optional[str]]:
        """
        Ensure CDR loop indices exist for a specific chain.

        Generates indices for CDR1, CDR2, and CDR3 loops based on sequences.

        Args:
            chain: TCR chain ('alpha' or 'beta')
            cdr1_seq: CDR1 amino acid sequence (optional)
            cdr2_seq: CDR2 amino acid sequence (optional)
            cdr3_seq: CDR3 amino acid sequence (optional)
            ca_only: If True, select only CA atoms (default: True)

        Returns:
            Dictionary mapping CDR loop names to group IDs
            e.g., {'CDR1_alpha': '25', 'CDR2_alpha': '26', 'CDR3_alpha': '27'}

        Example:
            >>> group_ids = index_mgr.ensure_cdr_loops_index(
            ...     chain='beta',
            ...     cdr1_seq='SGHNS',
            ...     cdr2_seq='SYNSPL',
            ...     cdr3_seq='CASSLGQAYEQYF'
            ... )
            >>> print(group_ids)
            {'CDR1_beta': '25', 'CDR2_beta': '26', 'CDR3_beta': '27'}
        """
        if chain.lower() not in ['alpha', 'beta']:
            raise ValueError("Chain must be 'alpha' or 'beta'")

        chain_id = 'D' if chain.lower() == 'alpha' else 'E'
        result_group_ids = {}

        # Process each CDR loop
        cdr_sequences = {
            'CDR1': cdr1_seq,
            'CDR2': cdr2_seq,
            'CDR3': cdr3_seq
        }

        for cdr_name, cdr_seq in cdr_sequences.items():
            if cdr_seq is None:
                result_group_ids[f"{cdr_name}_{chain.lower()}"] = None
                continue

            component_key = f"{cdr_name}_{chain.lower()}"

            with self._lock:
                # Check if already generated
                if component_key in self.dynamic_index_files:
                    existing_file = self.dynamic_index_files[component_key]
                    if existing_file.exists():
                        logger.debug(f"Using cached {cdr_name} index: {existing_file}")
                        group_info = self.group_registry.get(component_key)
                        result_group_ids[component_key] = str(group_info.group_id) if group_info else None
                        continue

                # Generate CDR index using sequence matching
                output_file = self.context.get_index_path(f"{cdr_name.lower()}_{chain.lower()}.ndx")

                logger.info(f"Generating {cdr_name} index for {chain} chain")
                logger.info(f"  Sequence: {cdr_seq}")

                # Use the same method as CDR3 (generate_cdr3_index works for any CDR)
                result_file = self._generate_cdr_index_by_sequence(
                    chain_id=chain_id,
                    cdr_sequence=cdr_seq,
                    component_name=component_key,
                    output_file=output_file,
                    ca_only=ca_only
                )

                if result_file is None:
                    logger.warning(f"{cdr_name} sequence not found in {chain} chain")
                    result_group_ids[component_key] = None
                    continue

                result_path = Path(result_file)
                self.dynamic_index_files[component_key] = result_path

                # Parse the CDR index
                self._parse_cdr3_index(result_path, component_key)

                group_info = self.group_registry.get(component_key)
                if group_info:
                    logger.info(f"{cdr_name} index created: group {group_info.group_id} with {group_info.atom_count} atoms")
                    result_group_ids[component_key] = str(group_info.group_id)
                else:
                    result_group_ids[component_key] = None

        return result_group_ids

    def ensure_interface_index(self,
                              cutoff: float = 4.5,
                              ca_only: bool = True) -> Optional[str]:
        """
        Ensure Interface index exists (TCR-pMHC contact residues).

        Generates an index containing residues at the TCR-pMHC interface
        based on distance cutoff.

        Args:
            cutoff: Distance cutoff in Angstroms (default: 4.5)
            ca_only: If True, select only CA atoms (default: True)

        Returns:
            Group ID as string, or None if generation failed

        Example:
            >>> interface_id = index_mgr.ensure_interface_index(cutoff=5.0)
            >>> print(interface_id)
            28
        """
        component_key = "INTERFACE"

        with self._lock:
            # Check if already generated
            if component_key in self.dynamic_index_files:
                existing_file = self.dynamic_index_files[component_key]
                if existing_file.exists():
                    logger.debug(f"Using cached Interface index: {existing_file}")
                    group_info = self.group_registry.get(component_key)
                    return str(group_info.group_id) if group_info else None

            # Generate Interface index
            output_file = self.context.get_index_path("interface.ndx")

            logger.info(f"Generating Interface index (cutoff={cutoff} Å)")

            try:
                result_file = self._generate_interface_index(
                    cutoff=cutoff,
                    output_file=output_file,
                    ca_only=ca_only
                )

                if result_file is None:
                    logger.warning("Interface index generation failed")
                    return None

                result_path = Path(result_file)
                self.dynamic_index_files[component_key] = result_path

                # Parse the interface index
                self._parse_cdr3_index(result_path, component_key)

                group_info = self.group_registry.get(component_key)
                if group_info:
                    logger.info(f"Interface index created: group {group_info.group_id} with {group_info.atom_count} atoms")
                    return str(group_info.group_id)
                else:
                    return None

            except Exception as e:
                logger.error(f"Failed to generate interface index: {e}")
                return None

    def _generate_cdr_index_by_sequence(self,
                                       chain_id: str,
                                       cdr_sequence: str,
                                       component_name: str,
                                       output_file: str,
                                       ca_only: bool = True) -> Optional[str]:
        """
        Generate CDR index by sequence matching (internal method).

        Args:
            chain_id: Chain ID ('D' or 'E')
            cdr_sequence: CDR amino acid sequence
            component_name: Component name for index group
            output_file: Output index file path
            ca_only: If True, select only CA atoms

        Returns:
            Path to generated index file, or None if failed
        """
        try:
            import MDAnalysis as mda
            from MDAnalysis.lib.util import convert_aa_code

            topology_file = self.context.structure_pdb or self.context.topology
            u = mda.Universe(str(topology_file))

            # Select chain
            chain_atoms = u.select_atoms(f"chainID {chain_id} and protein and name CA")
            chain_sequence = ''.join([
                convert_aa_code(r.resname) for r in chain_atoms.residues
            ])

            # Find sequence position
            pos = chain_sequence.find(cdr_sequence)
            if pos == -1:
                logger.warning(f"Sequence '{cdr_sequence}' not found in chain {chain_id}")
                return None

            start_idx = pos
            end_idx = pos + len(cdr_sequence)

            cdr_residues = chain_atoms.residues[start_idx:end_idx]

            # Build selection
            if ca_only:
                selection = f"chainID {chain_id} and protein and resid {cdr_residues[0].resid}-{cdr_residues[-1].resid} and name CA"
            else:
                selection = f"chainID {chain_id} and protein and resid {cdr_residues[0].resid}-{cdr_residues[-1].resid}"

            cdr_atoms = u.select_atoms(selection)
            atom_indices = cdr_atoms.indices + 1  # GROMACS uses 1-based indexing

            # Write index file
            with open(output_file, 'w') as f:
                f.write(f"[ {component_name} ]\n")
                for i in range(0, len(atom_indices), 10):
                    line = ' '.join(str(idx) for idx in atom_indices[i:i+10])
                    f.write(f"{line}\n")

            logger.info(f"{component_name} index created: {len(atom_indices)} atoms")
            return output_file

        except ImportError:
            logger.error("MDAnalysis is required for CDR index generation")
            return None
        except Exception as e:
            logger.error(f"Failed to generate CDR index: {e}")
            return None

    def _generate_interface_index(self,
                                  cutoff: float,
                                  output_file: str,
                                  ca_only: bool = True) -> Optional[str]:
        """
        Generate Interface index by distance calculation (internal method).

        Args:
            cutoff: Distance cutoff in Angstroms
            output_file: Output index file path
            ca_only: If True, select only CA atoms

        Returns:
            Path to generated index file, or None if failed
        """
        try:
            import MDAnalysis as mda
            from MDAnalysis.lib.distances import distance_array
            import numpy as np

            topology_file = self.context.structure_pdb or self.context.topology
            u = mda.Universe(str(topology_file))

            # Define TCR and pMHC selections
            tcr_sel = "chainID D E and protein"
            pmhc_sel = "chainID A B C and protein"

            if ca_only:
                tcr_sel += " and name CA"
                pmhc_sel += " and name CA"

            tcr_atoms = u.select_atoms(tcr_sel)
            pmhc_atoms = u.select_atoms(pmhc_sel)

            # Calculate distance matrix
            dist_matrix = distance_array(
                tcr_atoms.positions,
                pmhc_atoms.positions,
                box=u.dimensions
            )

            # Find interface residues (atoms within cutoff)
            tcr_interface_mask = np.any(dist_matrix <= cutoff, axis=1)
            pmhc_interface_mask = np.any(dist_matrix <= cutoff, axis=0)

            tcr_interface_atoms = tcr_atoms[tcr_interface_mask]
            pmhc_interface_atoms = pmhc_atoms[pmhc_interface_mask]

            # Combine interface atoms
            interface_atom_indices = np.concatenate([
                pmhc_interface_atoms.indices + 1,  # 1-based
                tcr_interface_atoms.indices + 1    # 1-based
            ])
            interface_atom_indices = np.sort(interface_atom_indices)

            # Write index file
            with open(output_file, 'w') as f:
                f.write("[ INTERFACE ]\n")
                for i in range(0, len(interface_atom_indices), 10):
                    line = ' '.join(str(idx) for idx in interface_atom_indices[i:i+10])
                    f.write(f"{line}\n")

            logger.info(f"Interface index created: {len(interface_atom_indices)} atoms "
                       f"({len(pmhc_interface_atoms)} from pMHC, {len(tcr_interface_atoms)} from TCR)")
            return output_file

        except ImportError:
            logger.error("MDAnalysis is required for interface index generation")
            return None
        except Exception as e:
            logger.error(f"Failed to generate interface index: {e}")
            return None

    def get_combined_index(self, *component_names: str) -> Path:
        """
        Get a combined index file containing multiple components.

        If components are all in base_index, returns base_index.
        If dynamic components (CDR loops, Interface) are included, merges base + dynamic indices.

        Args:
            *component_names: Variable number of component names

        Returns:
            Path to combined index file

        Example:
            >>> combined_index = index_mgr.get_combined_index('pHLA', 'TCR')
            >>> print(combined_index)
            ./results/1ao7/indices/base_components.ndx

            >>> combined_index = index_mgr.get_combined_index('pHLA', 'TCR', 'CDR3_beta')
            >>> print(combined_index)
            ./results/1ao7/indices/combined_with_dynamic.ndx
        """
        # Check if any components are dynamic (CDR loops, Interface)
        dynamic_components = [
            name for name in component_names
            if any(keyword in name for keyword in ['CDR', 'INTERFACE', 'Interface'])
        ]

        if not dynamic_components:
            # All components are static (in base index)
            self.ensure_base_index()
            return self.base_index_file
        else:
            # Need to merge base + dynamic indices
            with self._lock:
                combined_file = self.context.get_index_path("combined_with_dynamic.ndx")
                combined_path = Path(combined_file)

                # Check if already exists and is up-to-date
                if combined_path.exists():
                    logger.debug(f"Using cached combined index: {combined_path}")
                    return combined_path

                # Collect all index files to merge
                index_files = [str(self.base_index_file)]

                for comp_name in component_names:
                    if comp_name in self.dynamic_index_files:
                        index_files.append(str(self.dynamic_index_files[comp_name]))

                # Merge index files
                logger.info(f"Merging {len(index_files)} index files")
                self.index_generator.merge_index_files(
                    index_files=index_files,
                    output_file=str(combined_path)
                )

                logger.info(f"Combined index created: {combined_path}")
                return combined_path

    def _parse_index_file(self, index_file: Path) -> None:
        """
        Parse .ndx file to build group_name -> GroupInfo mapping.

        GROMACS .ndx file format:
        [ group_name ]
        atom_id1 atom_id2 atom_id3 ...

        Group IDs are determined by order in the file (0-indexed internally,
        but displayed/used as-is in GROMACS commands).

        Args:
            index_file: Path to .ndx file
        """
        if not index_file.exists():
            logger.error(f"Index file not found: {index_file}")
            return

        logger.debug(f"Parsing index file: {index_file}")

        current_group_id = 0
        current_group_name = None
        current_atom_ids = []

        with open(index_file, 'r') as f:
            for line in f:
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Check for group header: [ group_name ]
                match = re.match(r'\[\s*(.+?)\s*\]', line)
                if match:
                    # Save previous group if exists
                    if current_group_name is not None:
                        self._register_group(
                            current_group_id,
                            current_group_name,
                            len(current_atom_ids)
                        )
                        current_group_id += 1

                    # Start new group
                    current_group_name = match.group(1)
                    current_atom_ids = []
                else:
                    # Atom ID line
                    atom_ids = line.split()
                    current_atom_ids.extend(atom_ids)

            # Save last group
            if current_group_name is not None:
                self._register_group(
                    current_group_id,
                    current_group_name,
                    len(current_atom_ids)
                )

        logger.debug(f"Parsed {len(self.group_registry)} groups from {index_file.name}")

    def _parse_cdr3_index(self, index_file: Path, component_key: str) -> None:
        """
        Parse CDR3 index file (single group).

        Args:
            index_file: Path to CDR3 .ndx file
            component_key: Component key (e.g., 'CDR3_TCR_beta')
        """
        if not index_file.exists():
            logger.error(f"CDR3 index file not found: {index_file}")
            return

        logger.debug(f"Parsing CDR3 index: {index_file}")

        group_name = None
        atom_ids = []

        with open(index_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Check for group header
                match = re.match(r'\[\s*(.+?)\s*\]', line)
                if match:
                    group_name = match.group(1)
                else:
                    # Atom ID line
                    atom_ids.extend(line.split())

        if group_name:
            # CDR3 groups start at a high ID (e.g., 25+)
            # For simplicity, assign sequential IDs starting from 25
            cdr3_start_id = 25
            existing_cdr3_count = sum(1 for k in self.dynamic_index_files.keys() if k != component_key)
            group_id = cdr3_start_id + existing_cdr3_count

            self._register_group(group_id, component_key, len(atom_ids))
            logger.debug(f"CDR3 group registered: {component_key} (ID {group_id}, {len(atom_ids)} atoms)")

    def _register_group(self, group_id: int, group_name: str, atom_count: int) -> None:
        """
        Register a group in the registry.

        Args:
            group_id: Group ID (0-indexed in file)
            group_name: Group name
            atom_count: Number of atoms
        """
        # Try to match with StandardComponent
        component = None
        try:
            # Normalize name for matching
            normalized_name = group_name.upper().replace('_', '_')
            if normalized_name == 'PHLA':
                component = StandardComponent.PHLA
            elif normalized_name in [comp.name for comp in StandardComponent]:
                component = StandardComponent[normalized_name]
        except (KeyError, ValueError):
            pass

        group_info = GroupInfo(
            group_id=group_id,
            group_name=group_name,
            atom_count=atom_count,
            component=component
        )

        self.group_registry[group_name] = group_info
        logger.debug(f"Registered group: {group_name} (ID {group_id}, {atom_count} atoms)")

    def identify_chains_from_topology(
        self,
        topology_file: Optional[str] = None,
        output_index_file: Optional[str] = None
    ) -> 'TopologyChainMapping':
        """
        Identify chains directly from topology file (TPR/GRO) using length-based analysis.

        This method solves the PDB-TPR chain ID mismatch problem by working directly
        with the topology file. It uses chain splitting and residue counting to
        identify components independently of chain ID labels.

        Strategy:
        1. Split topology into chains using 'gmx make_ndx -splitch'
        2. Count atoms/residues for each chain
        3. Sort chains by length
        4. Assign components:
           - Shortest (5-22 res)     → peptide
           - 2nd shortest (90-110)   → beta2-microglobulin
           - Middle (100-220)        → TCR-alpha
           - 2nd longest (200-260)   → TCR-beta
           - Longest (260-290)       → HLA-alpha
        5. Generate renamed index file with standard component names

        Args:
            topology_file: Path to topology file (.tpr, .gro). If None, uses context.topology
            output_index_file: Optional output index file with renamed groups.
                             If None, uses default: results/{system_id}/indices/topology_components.ndx

        Returns:
            TopologyChainMapping with:
            - success: bool
            - chains: Dict[int, TopologyChainInfo]  # group_id -> chain_info
            - component_map: Dict[str, int]         # component_name -> group_id
            - error_message: Optional[str]

        Example:
            >>> # Identify chains from topology
            >>> mapping = index_mgr.identify_chains_from_topology()
            >>> if mapping.success:
            ...     peptide_group = mapping.component_map['peptide']
            ...     tcr_group_alpha = mapping.component_map['TCR_alpha']
            ...     print(f"Peptide is group {peptide_group}")
            ...     print(f"TCR-alpha is group {tcr_group_alpha}")

        Note:
            This method is preferred over PDB-based standardization when:
            - PDB file may have been standardized but TPR was not
            - Original MD simulation used non-standard chain ordering
            - Chain IDs in PDB and TPR/trajectory do not match
        """
        from immunoscope.analysis.topology import (
            TopologyChainIdentifier,
            TopologyChainMapping
        )

        # Use context topology if not provided
        if topology_file is None:
            topology_file = self.context.topology

        # Default output file
        if output_index_file is None:
            output_index_file = self.context.get_index_path("topology_components.ndx")

        logger.info(f"Identifying chains from topology: {topology_file}")

        # Create identifier and run analysis
        identifier = TopologyChainIdentifier()
        result = identifier.identify_chains_from_topology(
            topology_file=str(topology_file),
            output_index_file=str(output_index_file)
        )

        if result.success:
            logger.info("Topology-based chain identification successful")
            logger.info(f"Renamed index file: {output_index_file}")

            # Store the mapping in context metadata for future reference
            self.context.metadata['topology_chain_mapping'] = {
                comp_name: group_id
                for comp_name, group_id in result.component_map.items()
            }

            # Update chain standardization status
            self._chain_standardized = True
            self._chain_validation_result = {
                'standardized': True,
                'valid': True,
                'method': 'topology_based',
                'component_map': result.component_map,
                'chains': {
                    chain.assigned_component: chain.residue_count
                    for chain in result.chains.values()
                }
            }
        else:
            logger.error(f"Topology-based chain identification failed: {result.error_message}")
            self._chain_standardized = False
            self._chain_validation_result = {
                'standardized': False,
                'method': 'topology_based',
                'reason': result.error_message
            }

        return result

    def ensure_base_index_from_topology(
        self,
        force_regenerate: bool = False
    ) -> Path:
        """
        Alternative to ensure_base_index() using topology-based chain identification.

        This method works directly with TPR/GRO files to avoid PDB-TPR chain ID mismatch.
        It generates a base index file with standardized component names, but the
        underlying group IDs correspond to the original topology chain ordering.

        Args:
            force_regenerate: If True, regenerate even if cached

        Returns:
            Path to topology-based base index file

        Example:
            >>> # Use topology-based identification instead of PDB standardization
            >>> base_index = index_mgr.ensure_base_index_from_topology()
            >>> phla_id = index_mgr.get_group_id('pHLA')
            >>> # Now phla_id corresponds to correct groups in TPR, not PDB

        Note:
            This method is recommended when:
            - You suspect PDB and TPR have different chain orderings
            - You want to work directly with simulation files without PDB dependencies
            - Chain standardization warnings appear in logs
        """
        with self._lock:
            # Check if already generated
            topology_index = self.context.get_index_path("topology_components.ndx")
            topology_index_path = Path(topology_index)

            if self.base_index_file is not None and not force_regenerate:
                if self.base_index_file.exists():
                    # Check if it's the topology-based index
                    if self.base_index_file == topology_index_path:
                        logger.debug(f"Using cached topology-based index: {self.base_index_file}")
                        return self.base_index_file

            # Identify chains from topology
            result = self.identify_chains_from_topology(
                output_index_file=str(topology_index_path)
            )

            if not result.success:
                raise RuntimeError(
                    f"Failed to identify chains from topology: {result.error_message}"
                )

            # Set as base index
            self.base_index_file = topology_index_path

            # Parse the index file to build group registry
            self._parse_index_file(topology_index_path)

            # Register topology-identified components
            for comp_name, group_id in result.component_map.items():
                if comp_name in self.group_registry:
                    logger.debug(f"Component {comp_name} mapped to group {group_id}")
                else:
                    # If not found in registry, register manually
                    chain_info = result.chains.get(group_id)
                    if chain_info:
                        self._register_group(
                            group_id=group_id,
                            group_name=comp_name,
                            atom_count=chain_info.atom_count
                        )

            logger.info(f"Topology-based index created with {len(self.group_registry)} groups")
            logger.info("Using topology chain ordering - avoids PDB-TPR mismatch")

            return self.base_index_file

    def __repr__(self) -> str:
        """String representation for debugging."""
        n_groups = len(self.group_registry)
        has_base = self.base_index_file is not None
        n_dynamic = len(self.dynamic_index_files)
        return (f"IndexManager(system={self.context.system_id}, "
                f"base_index={has_base}, groups={n_groups}, dynamic={n_dynamic})")
