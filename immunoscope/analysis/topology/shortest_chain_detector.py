"""
Shortest Chain Detector for GROMACS MD Systems

This module provides functionality to detect the shortest protein chain
from GROMACS structure files (primarily .gro format) and generate index
files for PBC centering operations.

Designed for pHLA-TCR complexes where the peptide is typically the
shortest chain and should be used as the centering reference.
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class ShortestChainDetector:
    """
    Detect shortest protein chain from GRO files and generate GROMACS index files.

    This detector parses GRO format structure files to identify protein chains
    based on residue continuity, then generates a GROMACS index file containing
    the shortest chain for use in PBC correction and trajectory centering.

    Attributes:
        gro_file: Path to GROMACS .gro structure file
        topology_file: Path to topology file (.tpr, .gro, or .pdb)
        gmx: GROMACS executable command (default: 'gmx')
        protein_residues: Set of recognized protein residue names
    """

    def __init__(self, gro_file: str, topology_file: str, gmx_executable: str = "gmx"):
        """
        Initialize ShortestChainDetector.

        Args:
            gro_file: Path to .gro structure file
            topology_file: Path to topology file
            gmx_executable: GROMACS executable command (default: 'gmx')
        """
        self.gro_file = gro_file
        self.topology_file = topology_file
        self.gmx = gmx_executable

        # Standard and modified protein residue types
        self.protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'ACE', 'NME', 'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN'
        }

    def generate_shortest_chain_index(self, output_dir: str) -> Optional[str]:
        """
        Generate GROMACS index file containing the shortest protein chain.

        This is the main workflow:
        1. Detect all protein chains from GRO file
        2. Identify the shortest chain (typically peptide)
        3. Generate index file with standard groups + shortest chain

        Args:
            output_dir: Output directory for index file

        Returns:
            Path to generated index file, or None if failed

        Example:
            >>> detector = ShortestChainDetector("md.gro", "md.tpr")
            >>> index_file = detector.generate_shortest_chain_index("./output")
        """
        try:
            # 1. Detect chains from GRO file
            chains = self._detect_chains_from_gro()
            if not chains:
                logger.warning("No valid protein chains detected")
                return None

            # 2. Find shortest chain
            shortest_chain_atoms = self._find_shortest_chain_atoms(chains)
            if not shortest_chain_atoms:
                logger.warning("No shortest chain found")
                return None

            # 3. Generate index file
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            index_file = output_path / "shortest_chain.ndx"
            self._create_index_file_with_shortest_chain(str(index_file), shortest_chain_atoms)

            logger.info(f"Generated shortest chain index file: {index_file}")
            logger.info(f"Shortest chain contains {len(shortest_chain_atoms)} atoms")

            return str(index_file)

        except Exception as e:
            logger.error(f"Failed to generate index file: {e}")
            return None

    def _detect_chains_from_gro(self) -> Dict[int, List[int]]:
        """
        Detect protein chains from GRO file based on residue continuity.

        Parses GRO format to extract protein atoms, then splits into chains
        based on residue ID gaps (discontinuity > 1 indicates new chain).

        Returns:
            Dictionary mapping chain_id to list of atom IDs
            Example: {0: [1, 2, 3, ...], 1: [150, 151, ...], ...}

        Note:
            Only chains with >= 20 atoms are considered valid.
        """
        try:
            with open(self.gro_file, 'r') as f:
                lines = f.readlines()

            # Skip header (line 0), atom count (line 1), and box vector (last line)
            atom_lines = lines[2:-1]

            # Parse protein atoms
            protein_atoms = []
            for line in atom_lines:
                if len(line.strip()) < 20:
                    continue

                try:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    match = re.match(r'(\d+)([A-Za-z]+)', parts[0])
                    if match:
                        resid = int(match.group(1))
                        resname = match.group(2).upper()
                        atom_id = int(parts[2])
                    else:
                        resid = int(line[0:5].strip())
                        resname = line[5:10].strip().upper()
                        atom_id = int(line[15:20].strip())
                    coord_key = tuple(parts[3:6]) if len(parts) >= 6 else None

                    if resname in self.protein_residues:
                        protein_atoms.append({
                            'resid': resid,
                            'atom_id': atom_id,
                            'coord_key': coord_key,
                        })

                except (ValueError, IndexError):
                    continue

            # Sort by residue ID
            protein_atoms.sort(key=lambda x: x['resid'])

            # Split into chains based on residue ID continuity
            chains = {}
            current_chain = 0
            current_atoms = []
            last_resid = None
            last_coord_key = None

            for atom in protein_atoms:
                coord_changed = (
                    last_coord_key is not None
                    and atom.get('coord_key') is not None
                    and atom.get('coord_key') != last_coord_key
                )
                large_gap = last_resid is not None and atom['resid'] - last_resid > 10
                if last_resid is not None and (large_gap or coord_changed):
                    # New chain detected (gap in residue numbering)
                    if current_atoms:
                        chains[current_chain] = current_atoms
                        current_chain += 1
                    current_atoms = []

                current_atoms.append(atom['atom_id'])
                last_resid = atom['resid']
                last_coord_key = atom.get('coord_key')

            # Add last chain
            if current_atoms:
                chains[current_chain] = current_atoms

            logger.info(f"Detected {len(chains)} protein chains")
            for chain_id, atoms in chains.items():
                logger.info(f"  Chain {chain_id}: {len(atoms)} atoms")

            return chains

        except Exception as e:
            logger.error(f"Failed to detect chains from GRO file: {e}")
            return {}

    def _find_shortest_chain_atoms(self, chains: Dict[int, List[int]]) -> Optional[List[int]]:
        """
        Find the shortest chain from detected chains.

        Args:
            chains: Dictionary mapping chain_id to atom lists

        Returns:
            List of atom IDs for the shortest chain, or None if no chains
        """
        if not chains:
            return None

        # Find chain with minimum atom count
        shortest_chain_id = min(chains.keys(), key=lambda x: len(chains[x]))
        shortest_atoms = chains[shortest_chain_id]

        logger.info(f"Shortest chain: Chain {shortest_chain_id} ({len(shortest_atoms)} atoms)")
        return shortest_atoms

    def _create_index_file_with_shortest_chain(self, output_file: str, shortest_atoms: List[int]):
        """
        Create GROMACS index file with standard groups plus shortest chain.

        Workflow:
        1. Try to generate base index from topology file
        2. If fails, try from GRO file
        3. If fails, use minimal default groups
        4. Append shortest chain group

        Args:
            output_file: Path to output index file
            shortest_atoms: List of atom IDs for shortest chain

        Raises:
            Exception if file creation fails
        """

        # 1. Try to generate base index file from topology
        temp_base = output_file + ".temp"
        base_content = ""

        try:
            cmd = [self.gmx, "make_ndx", "-f", self.topology_file, "-o", temp_base]
            result = subprocess.run(
                cmd,
                input="q\n",
                text=True,
                capture_output=True,
                check=True,
                timeout=30
            )

            # Read generated base content
            if Path(temp_base).exists():
                with open(temp_base, 'r') as f:
                    base_content = f.read()
                Path(temp_base).unlink(missing_ok=True)

        except subprocess.CalledProcessError as e:
            logger.warning(f"Failed to generate base index from topology: {e}")
            logger.info("Trying to generate from GRO file instead")

            # Fallback: Generate from GRO file
            try:
                cmd = [self.gmx, "make_ndx", "-f", self.gro_file, "-o", temp_base]
                result = subprocess.run(
                    cmd,
                    input="q\n",
                    text=True,
                    capture_output=True,
                    check=True,
                    timeout=30
                )

                if Path(temp_base).exists():
                    with open(temp_base, 'r') as f:
                        base_content = f.read()
                    Path(temp_base).unlink(missing_ok=True)

            except subprocess.CalledProcessError as e2:
                logger.warning(f"Failed to generate from GRO file too: {e2}")
                # Last resort: Use minimal default index content
                base_content = self._create_minimal_index_content()

        # 2. Create final index file
        try:
            # Ensure output directory exists
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Format shortest chain group
            shortest_content = self._format_shortest_chain_group(shortest_atoms)

            # Write final file
            with open(output_file, 'w') as f:
                if base_content:
                    f.write(base_content)
                    if not base_content.endswith('\n'):
                        f.write("\n")
                f.write(shortest_content)

            logger.info(f"Successfully created index file with shortest chain: {output_file}")

        except Exception as e:
            logger.error(f"Failed to create final index file: {e}")
            # Cleanup temporary file if exists
            Path(temp_base).unlink(missing_ok=True)
            raise

    def _create_minimal_index_content(self) -> str:
        """
        Create minimal default index content as fallback.

        Returns:
            String containing basic GROMACS index groups
        """
        return """[ System ]
   1

[ Protein ]
   1

[ Protein-H ]
   1

[ C-alpha ]
   1

[ Backbone ]
   1

[ MainChain ]
   1

[ MainChain+Cb ]
   1

[ MainChain+H ]
   1

[ SideChain ]
   1

[ SideChain-H ]
   1

"""

    def _format_shortest_chain_group(self, atom_ids: List[int]) -> str:
        """
        Format shortest chain group in GROMACS index format.

        Args:
            atom_ids: List of atom IDs to include

        Returns:
            Formatted index group string with 12 atoms per line

        Format:
            [ Shortest_Chain ]
            1    2    3 ...  (12 atoms per line, right-aligned 5 chars each)
        """
        lines = ["[ Shortest_Chain ]"]

        # Format 12 atoms per line, right-aligned in 5-character fields
        for i in range(0, len(atom_ids), 12):
            line_atoms = atom_ids[i:i + 12]
            formatted_line = "".join(f"{atom_id:>5}" for atom_id in line_atoms)
            lines.append(formatted_line)

        return "\n".join(lines) + "\n"


def create_shortest_chain_index(md_gro_path: str, topology_path: str,
                               output_dir: str, gmx_executable: str = "gmx") -> Optional[str]:
    """
    Convenience function to create shortest chain index from md.gro file.

    This is the main entry point for external code to quickly generate
    a shortest chain index file without instantiating the detector class.

    Args:
        md_gro_path: Path to md.gro structure file
        topology_path: Path to topology file (.tpr, .gro, or .pdb)
        output_dir: Output directory for index file
        gmx_executable: GROMACS executable command (default: 'gmx')

    Returns:
        Path to generated index file, or None if failed

    Example:
        >>> from immunoscope.analysis import create_shortest_chain_index
        >>> index_file = create_shortest_chain_index(
        ...     "md.gro", "md.tpr", "./output"
        ... )
    """
    if not Path(md_gro_path).exists():
        logger.error(f"GRO file does not exist: {md_gro_path}")
        return None

    detector = ShortestChainDetector(md_gro_path, topology_path, gmx_executable)
    return detector.generate_shortest_chain_index(output_dir)
