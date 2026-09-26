#!/usr/bin/env python3
"""
Strict Shortest Chain Detector

For complex systems, ONLY the shortest chain (peptide) can be used as center group.
No fallback strategies allowed - if shortest chain detection fails, processing must stop.
"""

import re
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import tempfile

logger = logging.getLogger(__name__)


class StrictShortestChainDetector:
    """
    Strict detector that ONLY identifies shortest chain for centering.
    No fallback strategies - either find shortest chain or fail safely.
    """

    def __init__(self, topology_file: str, gro_file: Optional[str] = None, gmx_executable: str = "gmx"):
        self.topology_file = topology_file
        self.gro_file = gro_file
        self.gmx = gmx_executable

        # Extended protein residue types
        self.protein_residues = {
            # Standard amino acids
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # Modified residues
            'ACE', 'NME', 'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN',
            # Alternative names and caps
            'HSD', 'HSE', 'HSP', 'CSS', 'CYM', 'NTER', 'CTER'
        }

    def detect_shortest_chain_or_fail(self) -> Tuple[str, str]:
        """
        Simplified shortest chain detection using gmx make_ndx + splitch.

        Process:
        1. Use 'gmx make_ndx' to load standard groups
        2. Execute 'splitch 1' to split Protein group by chains
        3. Find the shortest newly created chain group (peptide)
        4. Return the group number for use in trjconv -center

        Returns:
            Tuple of (group_id, index_file_path)

        Raises:
            RuntimeError: If shortest chain cannot be detected
        """
        logger.info("Starting simplified strict shortest-chain detection with gmx make_ndx + splitch")

        try:
            # Execute simplified detection using gmx make_ndx + splitch
            result = self._detect_shortest_chain_splitch()

            if result:
                shortest_group_id, index_file = result
                logger.info("Detected shortest chain (peptide)")
                logger.info(f"   Group index: {shortest_group_id}")
                logger.info(f"   Index file: {index_file}")
                return shortest_group_id, index_file
            else:
                raise RuntimeError("pHLA-TCR splitch detection failed; peptide chain could not be identified")

        except Exception as e:
            error_msg = (
                f"pHLA-TCR shortest-chain detection failed: {e}\n"
                "   A pHLA-TCR complex is expected to contain HLA heavy chain, beta2m, peptide, and TCR alpha/beta chains.\n"
                "   Check:\n"
                "   1. Whether the TPR file contains all chain information.\n"
                "   2. Whether the GROMACS version supports the splitch command.\n"
                "   3. Whether chain boundaries are clearly represented.\n"
                "   4. Whether manual validation of gmx make_ndx output is needed."
            )
            logger.error(error_msg)
            raise RuntimeError("pHLA-TCR peptide detection failed; stopping to preserve scientific correctness")

    def _detect_shortest_chain_splitch(self) -> Optional[Tuple[str, str]]:
        """
        Detect the shortest chain with gmx make_ndx + splitch 1.

        Returns:
            Tuple of (group_id, index_file_path) or None if failed
        """
        try:
            # Create temporary index file
            temp_dir = tempfile.mkdtemp(prefix="immunoscope_splitch_")
            index_file = Path(temp_dir) / "chains.ndx"

            logger.info(f"Created temporary index file: {index_file}")

            # Use confirmed splitch syntax for the server's GROMACS version
            commands = "splitch 1\nq\n"

            # Execute gmx make_ndx
            logger.debug(f"Running command: {self.gmx} make_ndx -f {self.topology_file} -o {index_file}")
            logger.debug(f"Input commands: {repr(commands)}")

            result = subprocess.run(
                [self.gmx, "make_ndx", "-f", self.topology_file, "-o", str(index_file)],
                input=commands,
                text=True,
                capture_output=True,
                timeout=60
            )

            logger.info(f"Command return code: {result.returncode}")

            # Log the complete output for debugging pHLA-TCR systems
            logger.info("=== COMPLETE GMX MAKE_NDX OUTPUT ===")
            logger.info("STDERR:")
            logger.info(result.stderr)
            logger.info("STDOUT:")
            logger.info(result.stdout)
            logger.info("=== END OUTPUT ===")

            if result.returncode != 0:
                logger.warning(f"gmx make_ndx returned non-zero code {result.returncode}, but continuing to parse output")

            # Try to parse the output - check both stderr and stdout
            # In GROMACS 2023.2, chain info appears in stdout after splitch command

            # Parse stderr first (traditional location)
            chain_groups = self._parse_splitch_output(result.stderr)

            # If not found in stderr, try stdout (GROMACS 2023.2 behavior)
            if not chain_groups:
                logger.info("No chain groups found in stderr; trying stdout")
                chain_groups = self._parse_splitch_output(result.stdout)

            if not chain_groups:
                logger.error("No chain groups generated by splitch were found")
                return None

            # Find shortest chain group
            shortest_group = self._find_shortest_chain_group(chain_groups)

            if not shortest_group:
                logger.error("Unable to determine the shortest chain group")
                return None

            group_id, atom_count = shortest_group
            logger.info(f"Detected shortest chain: group {group_id}, {atom_count} atoms")

            # Validate the index file exists and is readable
            if not index_file.exists():
                logger.error("Index file was not created")
                return None

            return str(group_id), str(index_file)

        except Exception as e:
            logger.error(f"splitch detection failed: {e}")
            return None

    def _detect_from_gro_multiformat(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """Parse GRO files with multiple supported formats."""
        if not self.gro_file or not Path(self.gro_file).exists():
            logger.warning("GRO file does not exist; skipping GRO detection")
            return None

        try:
            # Method 1: Standard format parsing
            result = self._parse_gro_standard_format()
            if result and self._is_valid_chain_detection(result):
                return result

            # Method 2: Flexible format parsing
            result = self._parse_gro_flexible_format()
            if result and self._is_valid_chain_detection(result):
                return result

            # Method 3: Residue-gap based parsing
            result = self._parse_gro_residue_gaps()
            if result and self._is_valid_chain_detection(result):
                return result

            return None

        except Exception as e:
            logger.error(f"Multi-format GRO parsing failed: {e}")
            return None

    def _detect_from_combined_analysis(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """Analyze combined TPR and GRO information."""
        try:
            # Combine information from both sources
            # This would implement a more sophisticated approach
            # that cross-validates findings from TPR and GRO
            logger.info("Combined analysis method is not implemented yet")
            return None

        except Exception as e:
            logger.error(f"Combined analysis failed: {e}")
            return None

    def _parse_gro_standard_format(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """Parse standard GRO format with stricter format checks."""
        try:
            with open(self.gro_file, 'r') as f:
                lines = f.readlines()

            if len(lines) < 3:
                logger.warning("Invalid GRO file format")
                return None

            # Parse header
            title = lines[0].strip()
            try:
                atom_count = int(lines[1].strip())
            except ValueError:
                logger.warning("Unable to parse atom count")
                return None

            # Ensure we have the right number of lines
            expected_lines = atom_count + 3  # title + count + atoms + box
            if len(lines) < expected_lines:
                logger.warning(f"Insufficient GRO lines: expected {expected_lines}, actual {len(lines)}")
                return None

            atom_lines = lines[2:2+atom_count]
            chains = {}
            current_chain_id = 0
            current_atoms = []
            last_resid = None

            for line_num, line in enumerate(atom_lines, 1):
                try:
                    atom_info = self._parse_gro_line_strict(line, line_num)
                    if not atom_info:
                        continue

                    resid, resname, atom_id = atom_info

                    if resname in self.protein_residues:
                        # Chain boundary detection
                        if last_resid is not None and resid != last_resid + 1:
                            # New chain detected
                            if len(current_atoms) >= 5:  # Minimum peptide size
                                chains[current_chain_id] = current_atoms
                                current_chain_id += 1
                            current_atoms = []

                        current_atoms.append(atom_id)
                        last_resid = resid

                except Exception as e:
                    logger.debug(f"Failed to parse line {line_num}: {e}")
                    continue

            # Add final chain
            if len(current_atoms) >= 5:
                chains[current_chain_id] = current_atoms

            if len(chains) < 2:
                logger.warning(f"Too few chains detected: {len(chains)}")
                return None

            # Find shortest chain
            shortest_chain_atoms = min(chains.values(), key=len)

            logger.info(f"Standard format parsing detected {len(chains)} chains")
            for i, atoms in chains.items():
                logger.info(f"  Chain {i}: {len(atoms)} atoms")

            return chains, shortest_chain_atoms

        except Exception as e:
            logger.error(f"Standard GRO parsing failed: {e}")
            return None

    def _parse_gro_line_strict(self, line: str, line_num: int) -> Optional[Tuple[int, str, int]]:
        """Strict GRO line parser supporting multiple format variants."""

        if len(line.strip()) < 30:
            return None

        try:
            # GRO format: positions are fixed
            # %5d%-5s%5s%5d%8.3f%8.3f%8.3f
            # resid resname atomname atomid x y z

            # Method 1: Fixed position parsing (most common)
            if len(line) >= 44:
                resid_resname = line[0:5].strip()
                atom_name = line[5:10].strip()
                atom_id_str = line[10:15].strip()

                # Parse residue ID and name
                patterns = [
                    r'^(\d+)([A-Z]{3,4})$',      # 1ALA, 123PRO
                    r'^(\d+)\s+([A-Z]{3,4})$',   # "1 ALA"
                    r'^(\d+)([A-Z]+)$'           # 1ACE, 1NTER
                ]

                for pattern in patterns:
                    match = re.match(pattern, resid_resname)
                    if match:
                        resid = int(match.group(1))
                        resname = match.group(2)
                        atom_id = int(atom_id_str)
                        return resid, resname, atom_id

            # Method 2: Space-separated parsing (fallback)
            parts = line.split()
            if len(parts) >= 6:
                # Try to parse first field as "123ALA" format
                first_field = parts[0]
                match = re.match(r'^(\d+)([A-Z]+)', first_field)
                if match:
                    resid = int(match.group(1))
                    resname = match.group(2)
                    # Atom ID usually in 3rd position for space-separated
                    atom_id = int(parts[2])
                    return resid, resname, atom_id

            return None

        except (ValueError, IndexError) as e:
            logger.debug(f"Line {line_num} parsing failed: {e}")
            return None

    def _parse_gro_flexible_format(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """Flexible GRO parser for format variants."""
        # Implementation similar to standard but with more flexible parsing
        logger.info("Flexible format parsing is not fully implemented yet")
        return None

    def _parse_gro_residue_gaps(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """Detect chains based on residue-number gaps."""
        # Implementation that looks for larger gaps in residue numbering
        logger.info("Residue-gap detection is not fully implemented yet")
        return None

    def _parse_splitch_output(self, stderr_output: str) -> List[Tuple[int, int]]:
        """
        Parse gmx make_ndx splitch output and extract chain-group information.

        Args:
            stderr_output: stderr or stdout output from gmx make_ndx

        Returns:
            List of (group_id, atom_count) tuples for newly created chain groups
        """
        chain_groups = []

        try:
            lines = stderr_output.split('\n')

            # Focus on parsing - detailed output already logged above

            # Multiple patterns to try, ordered by specificity
            patterns = [
                # Pattern 1: Actual GROMACS 2023.2 splitch output format
                r'^(\d+):\s*(\d+)\s+atoms\s*\(\d+\s+to\s+\d+\)',
                # Pattern 2: Standard splitch output with "Protein_chain_"
                r'^\s*(\d+)\s+Protein_chain_\w+\s*:\s*(\d+)\s+atoms',
                # Pattern 3: General chain pattern (more flexible)
                r'^\s*(\d+)\s+\w*[Cc]hain\w*\s*:\s*(\d+)\s+atoms',
                # Pattern 4: Very flexible pattern for any group >= 18
                r'^\s*(\d+)\s+(\S+)\s*:\s*(\d+)\s+atoms'
            ]

            # Try each pattern
            for pattern_num, pattern in enumerate(patterns, 1):
                pattern_chains = []

                for line in lines:
                    match = re.search(pattern, line)
                    if match:
                        if pattern_num == 1:  # Pattern 1: GROMACS 2023.2 format "1: 4343 atoms"
                            chain_id = int(match.group(1))
                            atom_count = int(match.group(2))
                            # Convert chain ID to group ID (add 17 to get group number)
                            group_id = chain_id + 17
                            pattern_chains.append((group_id, atom_count))
                            logger.debug(f"Pattern {pattern_num} matched: chain {chain_id}->group {group_id}, {atom_count} atoms")

                        elif pattern_num <= 3:  # Patterns 2-3: direct group_id, atom_count
                            group_id = int(match.group(1))
                            atom_count = int(match.group(2))
                            # Only accept newly created groups (>= 18)
                            if group_id >= 18:
                                pattern_chains.append((group_id, atom_count))
                                logger.debug(f"Pattern {pattern_num} matched: group {group_id}, {atom_count} atoms")

                        else:  # Pattern 4: group_id, group_name, atom_count
                            group_id = int(match.group(1))
                            group_name = match.group(2)
                            atom_count = int(match.group(3))

                            # For pattern 4, only accept groups >= 18 with chain-like names
                            if group_id < 18:
                                continue
                            if not any(keyword in group_name.lower() for keyword in ['chain', 'prot']):
                                continue

                            pattern_chains.append((group_id, atom_count))
                            logger.debug(f"Pattern {pattern_num} matched: group {group_id}, {atom_count} atoms")

                if pattern_chains:
                    chain_groups = pattern_chains
                    logger.info(f"Pattern {pattern_num} parsed {len(chain_groups)} chain groups")
                    break

            # If still no chains found, show diagnostic info
            if not chain_groups:
                logger.warning("No chain groups found; showing diagnostics:")

                # Look for "Splitting" message
                splitting_found = False
                for line in lines:
                    if 'splitting' in line.lower() or 'splitch' in line.lower():
                        logger.warning(f"Found splitting message: {line.strip()}")
                        splitting_found = True

                if not splitting_found:
                    logger.warning("No splitting message found; splitch may not have executed")

                # Show any lines with numbers that might be groups
                logger.warning("Potential group lines:")
                for line in lines:
                    if re.search(r'^\s*\d+\s+\S+.*atoms', line):
                        logger.warning(f"  {line.strip()}")

            logger.info(f"Parsed {len(chain_groups)} chain groups")
            return chain_groups

        except Exception as e:
            logger.error(f"Failed to parse splitch output: {e}")
            return []

    def _find_shortest_chain_group(self, chain_groups: List[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        """
        Find the shortest chain (peptide) from chain groups.

        Args:
            chain_groups: List of (group_id, atom_count) tuples

        Returns:
            (group_id, atom_count) of shortest chain, or None if invalid
        """
        if not chain_groups:
            logger.error("No chain groups available for selection")
            return None

        if len(chain_groups) < 2:
            logger.error(f"Too few chains detected: {len(chain_groups)} < 2; not a complex system")
            return None

        # Sort by atom count to find shortest
        sorted_chains = sorted(chain_groups, key=lambda x: x[1])
        shortest_group_id, shortest_atom_count = sorted_chains[0]

        # Validate shortest chain is reasonable for peptide (≤20 AAs = ~300 atoms)
        if shortest_atom_count < 10:
            logger.error(f"Shortest chain is too short: {shortest_atom_count} < 10 atoms; likely not a valid peptide")
            return None

        if shortest_atom_count > 300:
            logger.error(f"Shortest chain is too long: {shortest_atom_count} > 300 atoms; outside peptide range (<=20 AAs)")
            return None

        # Check size distribution - ensure meaningful difference
        longest_atom_count = sorted_chains[-1][1]
        if shortest_atom_count * 2 > longest_atom_count:
            logger.warning(f"Chain size separation is weak: shortest {shortest_atom_count}, longest {longest_atom_count}")

        logger.info(f"Selected shortest chain as peptide: group {shortest_group_id}, {shortest_atom_count} atoms")
        for group_id, atom_count in sorted_chains:
            logger.info(f"  Chain group {group_id}: {atom_count} atoms")

        return shortest_group_id, shortest_atom_count

    def _is_valid_chain_detection(self, result: Tuple[Dict[int, List[int]], List[int]]) -> bool:
        """Validate chain detection result."""
        chains, shortest_chain_atoms = result

        # Must have at least 2 chains (complex system)
        if len(chains) < 2:
            logger.warning(f"Too few chains detected: {len(chains)} < 2")
            return False

        # Shortest chain must be reasonable peptide size
        if len(shortest_chain_atoms) < 5:
            logger.warning(f"Shortest chain is too short: {len(shortest_chain_atoms)} < 5")
            return False

        # Shortest chain shouldn't be too long (likely not a peptide)
        if len(shortest_chain_atoms) > 1000:
            logger.warning(f"Shortest chain is too long: {len(shortest_chain_atoms)} > 1000")
            return False

        # Check chain size distribution
        chain_sizes = [len(atoms) for atoms in chains.values()]
        min_size = min(chain_sizes)
        max_size = max(chain_sizes)

        # Ensure there's meaningful size difference
        if min_size * 3 > max_size:
            logger.warning(f"Chain size separation is weak: min {min_size}, max {max_size}")
            return False

        logger.info("Chain detection result validation passed")
        logger.info(f"   Chain count: {len(chains)}")
        logger.info(f"   Shortest chain: {min_size} atoms")
        logger.info(f"   Longest chain: {max_size} atoms")

        return True

    def _validate_shortest_chain_result(self, result: Tuple[Dict[int, List[int]], List[int]]) -> bool:
        """Validate shortest-chain detection result."""
        return self._is_valid_chain_detection(result)


    def _validate_index_file(self, index_file: str, group_id: str) -> bool:
        """Validate that the index file contains the specified group."""
        try:
            with open(index_file, 'r') as f:
                content = f.read()

            # Check if the group ID exists in the file
            # Look for group headers like "[ Protein_chain_A ]" or similar
            group_pattern = rf'\[\s*\w*\s*\]'
            groups_found = re.findall(group_pattern, content)

            if len(groups_found) >= int(group_id):
                logger.info(f"Index file validation passed: contains group {group_id}")
                return True
            else:
                logger.error(f"Index file validation failed: group {group_id} not found")
                return False

        except Exception as e:
            logger.error(f"Index file validation failed: {e}")
            return False


def detect_shortest_chain_strict(topology_file: str,
                               gro_file: Optional[str] = None,
                               gmx_executable: str = "gmx") -> Tuple[str, str]:
    """
    Simplified strict shortest-chain detection using gmx make_ndx + splitch.

    Process:
    1. Run 'gmx make_ndx -f topology.tpr' to load standard groups.
    2. Run 'splitch 1' to split the Protein group by chain.
    3. Parse the newly generated chain groups and select the shortest one.
    4. Return the group index for later 'gmx trjconv -center' use.

    Args:
        topology_file: TPR file path.
        gro_file: GRO file path; currently unused and kept for API compatibility.
        gmx_executable: GROMACS executable path.

    Returns:
        Tuple of (group_id, index_file_path)
        - group_id: shortest-chain group index, such as "18" or "19".
        - index_file_path: generated index file path.

    Raises:
        RuntimeError: If the shortest chain cannot be detected.

    Example:
        group_id, index_file = detect_shortest_chain_strict("system.tpr")
        # Use group_id for centering in trjconv:
        # gmx trjconv -center -pbc mol -n index_file
    """
    detector = StrictShortestChainDetector(topology_file, gro_file, gmx_executable)
    return detector.detect_shortest_chain_or_fail()
