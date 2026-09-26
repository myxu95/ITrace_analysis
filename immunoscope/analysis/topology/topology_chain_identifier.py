"""
Topology-based Chain Identifier for GROMACS Files

Identifies chain components directly from TPR/GRO topology files
based on chain splitting and atom/residue counting, avoiding
chain ID label dependencies.

This approach solves the PDB-TPR chain ID mismatch problem:
- PDB may be standardized (A/B/C/D/E)
- TPR/trajectory may have original chain IDs
- This module works directly with TPR, using physical properties (chain length)

Author: ImmunoScope Development Team
Date: 2026-03-17
"""

import logging
import subprocess
import tempfile
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class TopologyChainInfo:
    """Information about a chain from topology file."""
    original_group_id: int      # Original group ID from splitch
    original_group_name: str    # Original group name (e.g., "ch0_Protein")
    atom_count: int             # Number of atoms
    residue_count: int          # Number of residues
    assigned_component: str     # Assigned component name (e.g., "peptide", "HLA_alpha")
    confidence: float           # Assignment confidence (0.0-1.0)


@dataclass
class TopologyChainMapping:
    """Complete chain mapping result."""
    success: bool
    chains: Dict[int, TopologyChainInfo]  # group_id -> chain_info
    component_map: Dict[str, int]         # component_name -> group_id
    error_message: Optional[str] = None


class TopologyChainIdentifier:
    """
    Identify chain components from topology files based on atom/residue counts.

    Strategy:
    1. Use 'gmx make_ndx -splitch' to split topology by chains
    2. Extract atom count and residue count for each chain
    3. Sort chains by length and assign components:
       - Shortest (5-22 res)     → peptide
       - 2nd shortest (90-110)   → beta2-microglobulin (β2m)
       - Middle (100-220)        → TCR-alpha
       - 2nd longest (200-260)   → TCR-beta
       - Longest (260-290)       → HLA-alpha

    This method:
    - Works directly with TPR/GRO files
    - Independent of chain ID labels
    - Avoids PDB-TPR chain ID mismatch issues
    """

    # Length ranges for component identification
    LENGTH_RANGES = {
        'peptide': {'min_res': 5, 'max_res': 22, 'typical': 9},
        'beta2m': {'min_res': 90, 'max_res': 110, 'typical': 99},
        'TCR_alpha': {'min_res': 100, 'max_res': 220, 'typical': 180},
        'TCR_beta': {'min_res': 200, 'max_res': 260, 'typical': 240},
        'HLA_alpha': {'min_res': 260, 'max_res': 290, 'typical': 275}
    }

    def __init__(self, gmx_executable: str = "gmx"):
        """
        Initialize topology chain identifier.

        Args:
            gmx_executable: GROMACS executable command
        """
        self.gmx = gmx_executable

    def identify_chains_from_topology(
        self,
        topology_file: str,
        output_index_file: Optional[str] = None
    ) -> TopologyChainMapping:
        """
        Identify chains from topology file and generate renamed index.

        Args:
            topology_file: Path to topology file (.tpr, .gro, .pdb)
            output_index_file: Optional output index file with renamed groups

        Returns:
            TopologyChainMapping with chain identification results

        Example:
            >>> identifier = TopologyChainIdentifier()
            >>> result = identifier.identify_chains_from_topology("md.tpr")
            >>> if result.success:
            ...     peptide_group = result.component_map['peptide']
            ...     print(f"Peptide is group {peptide_group}")
        """
        logger.info(f"Identifying chains from topology: {topology_file}")

        try:
            # Step 1: Split chains using gmx make_ndx
            chain_groups = self._split_chains_gmx(topology_file)

            if not chain_groups:
                return TopologyChainMapping(
                    success=False,
                    chains={},
                    component_map={},
                    error_message="No protein chains found via splitch"
                )

            # Step 2: Extract chain information
            chains_info = self._extract_chain_info(chain_groups, topology_file)

            if len(chains_info) != 5:
                logger.warning(f"Expected 5 chains, found {len(chains_info)}")
                logger.warning("Chain identification may be incorrect")

            # Step 3: Assign components based on length
            assigned_chains = self._assign_components_by_length(chains_info)

            # Step 4: Build component map
            component_map = {
                chain.assigned_component: group_id
                for group_id, chain in assigned_chains.items()
            }

            # Step 5: Generate renamed index file (optional)
            if output_index_file:
                self._generate_renamed_index(
                    assigned_chains,
                    output_index_file,
                    topology_file
                )

            logger.info("✅ Chain identification successful:")
            for comp_name, group_id in sorted(component_map.items()):
                chain = assigned_chains[group_id]
                logger.info(f"  {comp_name:12s} → Group {group_id:2d} "
                           f"({chain.residue_count} residues, {chain.atom_count} atoms)")

            return TopologyChainMapping(
                success=True,
                chains=assigned_chains,
                component_map=component_map
            )

        except Exception as e:
            logger.error(f"Chain identification failed: {e}")
            return TopologyChainMapping(
                success=False,
                chains={},
                component_map={},
                error_message=str(e)
            )

    def _split_chains_gmx(self, topology_file: str) -> List[Dict]:
        """
        Split chains using 'gmx make_ndx -splitch'.

        Args:
            topology_file: Path to topology file

        Returns:
            List of chain group dictionaries with group_id, name, etc.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_index = Path(tmpdir) / "chains.ndx"

            # Run gmx make_ndx with splitch
            cmd = [
                self.gmx, 'make_ndx',
                '-f', topology_file,
                '-o', str(temp_index)
            ]

            logger.debug(f"Running: {' '.join(cmd)}")

            try:
                process = subprocess.run(
                    cmd,
                    input="splitch\nq\n",  # splitch command, then quit
                    text=True,
                    capture_output=True,
                    timeout=30
                )

                # Parse stderr for chain group information
                chain_groups = self._parse_splitch_output(process.stderr)

                return chain_groups

            except subprocess.TimeoutExpired:
                logger.error("gmx make_ndx timed out")
                return []
            except Exception as e:
                logger.error(f"Failed to split chains: {e}")
                return []

    def _parse_splitch_output(self, stderr_output: str) -> List[Dict]:
        """
        Parse gmx make_ndx stderr to extract chain groups.

        Expected format in stderr:
        "Splitting Protein into chains (using column 21 of the pdb file).
         There are 5 chains and 0 residues with unknown chain ID."

        Then group list shows:
        "18 ch0_Protein   : 1234 atoms"
        "19 ch1_Protein   : 567 atoms"
        ...

        Args:
            stderr_output: stderr from gmx make_ndx

        Returns:
            List of dicts with group_id, name, atom_count
        """
        chain_groups = []

        # Pattern: group_id group_name : atom_count atoms
        # Example: "18 ch0_Protein   : 1234 atoms"
        pattern = r'^\s*(\d+)\s+(ch\d+_Protein)\s*:\s*(\d+)\s+atoms'

        for line in stderr_output.split('\n'):
            match = re.match(pattern, line.strip())
            if match:
                group_id = int(match.group(1))
                group_name = match.group(2)
                atom_count = int(match.group(3))

                chain_groups.append({
                    'group_id': group_id,
                    'name': group_name,
                    'atom_count': atom_count
                })

                logger.debug(f"Found chain group: {group_name} (ID {group_id}, {atom_count} atoms)")

        logger.info(f"Detected {len(chain_groups)} protein chains")
        return chain_groups

    def _extract_chain_info(
        self,
        chain_groups: List[Dict],
        topology_file: str
    ) -> List[TopologyChainInfo]:
        """
        Extract detailed chain information including residue count.

        Args:
            chain_groups: List of chain groups from splitch
            topology_file: Topology file path

        Returns:
            List of TopologyChainInfo objects
        """
        chains_info = []

        for group in chain_groups:
            # Estimate residue count from atom count
            # Typical: ~8-10 atoms per residue for backbone+sidechain
            # Use conservative estimate: atom_count / 8
            estimated_residues = group['atom_count'] // 8

            chain_info = TopologyChainInfo(
                original_group_id=group['group_id'],
                original_group_name=group['name'],
                atom_count=group['atom_count'],
                residue_count=estimated_residues,
                assigned_component='',  # Will be assigned later
                confidence=0.0
            )

            chains_info.append(chain_info)

        return chains_info

    def _assign_components_by_length(
        self,
        chains_info: List[TopologyChainInfo]
    ) -> Dict[int, TopologyChainInfo]:
        """
        Assign components to chains based on residue count.

        Strategy:
        1. Sort chains by residue count (ascending)
        2. Assign in order: peptide, β2m, TCR-α, TCR-β, HLA-α
        3. Validate against expected ranges
        4. Set confidence based on range match

        Args:
            chains_info: List of chain information

        Returns:
            Dictionary mapping group_id to assigned TopologyChainInfo
        """
        # Sort chains by residue count
        sorted_chains = sorted(chains_info, key=lambda c: c.residue_count)

        # Expected order: peptide < β2m < TCR-α < TCR-β < HLA-α
        component_order = ['peptide', 'beta2m', 'TCR_alpha', 'TCR_beta', 'HLA_alpha']

        assigned_chains = {}

        for idx, chain in enumerate(sorted_chains):
            if idx < len(component_order):
                component = component_order[idx]
                chain.assigned_component = component

                # Calculate confidence based on length range match
                length_range = self.LENGTH_RANGES[component]
                min_res = length_range['min_res']
                max_res = length_range['max_res']

                if min_res <= chain.residue_count <= max_res:
                    chain.confidence = 1.0
                    logger.info(f"✓ {component:12s}: {chain.residue_count} residues (expected {min_res}-{max_res})")
                else:
                    # Out of range - lower confidence
                    chain.confidence = 0.5
                    logger.warning(f"⚠ {component:12s}: {chain.residue_count} residues "
                                 f"(expected {min_res}-{max_res})")

                assigned_chains[chain.original_group_id] = chain

        return assigned_chains

    def _generate_renamed_index(
        self,
        assigned_chains: Dict[int, TopologyChainInfo],
        output_file: str,
        topology_file: str
    ) -> None:
        """
        Generate index file with renamed groups based on assignments.

        Creates groups like: pHLA, peptide, TCR, HLA_alpha, etc.

        Args:
            assigned_chains: Dictionary of assigned chains
            output_file: Output index file path
            topology_file: Topology file path
        """
        # Build gmx make_ndx commands to create renamed groups
        commands = []

        # Get group IDs by component
        component_groups = {}
        for group_id, chain in assigned_chains.items():
            component_groups[chain.assigned_component] = group_id

        # Create individual component groups
        for component, group_id in component_groups.items():
            commands.append(f"{group_id}")  # Select group
            commands.append(f"name {20 + len(commands)//2} {component}")

        # Create combined groups
        # pHLA = HLA_alpha + beta2m + peptide
        if all(c in component_groups for c in ['HLA_alpha', 'beta2m', 'peptide']):
            hla_id = component_groups['HLA_alpha']
            b2m_id = component_groups['beta2m']
            pep_id = component_groups['peptide']
            commands.append(f"{hla_id} | {b2m_id} | {pep_id}")
            commands.append(f"name {20 + len(commands)//2} pHLA")

        # TCR = TCR_alpha + TCR_beta
        if all(c in component_groups for c in ['TCR_alpha', 'TCR_beta']):
            tcr_a = component_groups['TCR_alpha']
            tcr_b = component_groups['TCR_beta']
            commands.append(f"{tcr_a} | {tcr_b}")
            commands.append(f"name {20 + len(commands)//2} TCR")

        commands.append("q")  # Quit

        # Run gmx make_ndx with commands
        cmd = [self.gmx, 'make_ndx', '-f', topology_file, '-o', output_file]

        stdin_input = '\n'.join(commands) + '\n'

        try:
            process = subprocess.run(
                cmd,
                input=stdin_input,
                text=True,
                capture_output=True,
                timeout=30
            )

            if process.returncode == 0:
                logger.info(f"✅ Renamed index file created: {output_file}")
            else:
                logger.error(f"Failed to create renamed index: {process.stderr}")

        except Exception as e:
            logger.error(f"Error generating renamed index: {e}")

    def get_component_selection_string(
        self,
        component_name: str,
        mapping: TopologyChainMapping
    ) -> Optional[str]:
        """
        Get GROMACS selection string for a component.

        Args:
            component_name: Component name (e.g., 'peptide', 'TCR')
            mapping: TopologyChainMapping result

        Returns:
            Group ID as string, or None if not found

        Example:
            >>> group_id = identifier.get_component_selection_string('peptide', mapping)
            >>> print(group_id)  # "18"
        """
        if not mapping.success:
            return None

        # Handle combined components
        if component_name == 'pHLA':
            # Combine HLA_alpha, beta2m, peptide
            group_ids = []
            for comp in ['HLA_alpha', 'beta2m', 'peptide']:
                if comp in mapping.component_map:
                    group_ids.append(str(mapping.component_map[comp]))
            return ' | '.join(group_ids) if group_ids else None

        elif component_name == 'TCR':
            # Combine TCR_alpha, TCR_beta
            group_ids = []
            for comp in ['TCR_alpha', 'TCR_beta']:
                if comp in mapping.component_map:
                    group_ids.append(str(mapping.component_map[comp]))
            return ' | '.join(group_ids) if group_ids else None

        elif component_name == 'HLA':
            # HLA_alpha + beta2m
            group_ids = []
            for comp in ['HLA_alpha', 'beta2m']:
                if comp in mapping.component_map:
                    group_ids.append(str(mapping.component_map[comp]))
            return ' | '.join(group_ids) if group_ids else None

        else:
            # Individual component
            if component_name in mapping.component_map:
                return str(mapping.component_map[component_name])

        return None
