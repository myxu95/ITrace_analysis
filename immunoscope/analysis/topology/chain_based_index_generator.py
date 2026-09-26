"""
Chain-Based Index Generator for GROMACS MD Systems

This module provides functionality to generate GROMACS index files based on
intelligent chain identification using ANARCI for TCR recognition.

Chain identification strategy for Class I pHLA-TCR complexes:
1. Shortest chain (5-22 AA) -> Peptide
2. ~100 AA chain (90-110 AA) -> Beta-2-microglobulin
3. Remaining 3 chains -> Use ANARCI to identify TCR chains
4. Longest non-TCR chain -> MHC alpha

Author: ImmunoScope Development Team
Date: 2026-01-21
"""

import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import tempfile

logger = logging.getLogger(__name__)


class ChainBasedIndexGenerator:
    """
    Generate GROMACS index files with intelligent chain identification.

    Supports flexible input:
    - PDB file (protein structure)
    - TPR file (extract PDB temporarily using gmx trjconv)

    Chain identification strategy for Class I pHLA-TCR:
    1. Shortest chain (5-22 AA) -> Peptide
    2. ~100 AA chain (90-110 AA) -> Beta-2-microglobulin
    3. Remaining chains -> Use ANARCI to identify TCR-alpha and TCR-beta
    4. Longest non-TCR chain -> MHC-alpha

    Attributes:
        structure_source: Path to input structure file (PDB or TPR)
        input_type: Type of input file ('pdb' or 'tpr')
        gmx: GROMACS executable command
        pdb_file: Path to PDB file (direct or extracted)
        anarci_wrapper: Optional ANARCIWrapper for TCR identification
    """

    def __init__(
        self,
        structure_source: str,
        gmx_executable: str = "gmx",
        use_anarci: bool = True
    ):
        """
        Initialize ChainBasedIndexGenerator.

        Args:
            structure_source: Path to structure file (.pdb or .tpr)
            gmx_executable: GROMACS executable command (default: 'gmx')
            use_anarci: Use ANARCI for automatic chain identification (default: True)

        Raises:
            ValueError: If file type is not .pdb or .tpr
            FileNotFoundError: If structure file does not exist
        """
        self.structure_source = Path(structure_source)
        self.gmx = gmx_executable
        self._temp_files = []
        self.use_anarci = use_anarci
        self.anarci_wrapper = None

        # Initialize ANARCI wrapper if requested
        if self.use_anarci:
            try:
                from .cdr_manager import ANARCIWrapper
                self.anarci_wrapper = ANARCIWrapper(allow_fallback=False)
                logger.info("ANARCI available for chain identification")
            except Exception as e:
                logger.warning(f"ANARCI not available: {e}. Chain identification will require manual specification.")
                self.use_anarci = False

        if not self.structure_source.exists():
            raise FileNotFoundError(f"Structure file not found: {structure_source}")

        # Auto-detect input type
        suffix = self.structure_source.suffix.lower()
        if suffix == '.pdb':
            self.input_type = 'pdb'
            self.pdb_file = str(self.structure_source)
            logger.info(f"Input type: PDB file")
        elif suffix == '.tpr':
            self.input_type = 'tpr'
            self.pdb_file = None  # Will be extracted
            logger.info(f"Input type: TPR file (will extract PDB)")
        else:
            raise ValueError(
                f"Unsupported file type: {suffix}. "
                f"Expected .pdb or .tpr"
            )

    def generate_peptide_index(
        self,
        output_dir: str,
        peptide_chain_id: str = 'C',
        index_name: str = "peptide_chain.ndx"
    ) -> str:
        """
        Generate GROMACS index file for peptide chain.

        Args:
            output_dir: Output directory for index file
            peptide_chain_id: Chain ID of peptide (default: 'C' for standardized PDB)
            index_name: Name of output index file

        Returns:
            Path to generated index file

        Raises:
            ValueError: If peptide chain not found in structure
            RuntimeError: If index generation fails

        Example:
            >>> generator = ChainBasedIndexGenerator("complex_std_protein.pdb")
            >>> index_file = generator.generate_peptide_index("./")
            >>> print(f"Generated: {index_file}")
        """
        try:
            # 1. Get PDB file (existing or extract from TPR)
            pdb_file = self._get_pdb_file()

            # 2. Read peptide chain information from PDB
            peptide_info = self._read_chain_info(pdb_file, peptide_chain_id)

            logger.info(
                f"Peptide chain '{peptide_chain_id}': "
                f"{peptide_info['num_residues']} residues, "
                f"{peptide_info['num_atoms']} atoms"
            )

            # 3. Generate index file
            index_file = self._create_index_file(
                pdb_file,
                peptide_info,
                output_dir,
                index_name
            )

            logger.info(f"Successfully generated index file: {index_file}")

            return index_file

        finally:
            # Always cleanup temporary files
            self._cleanup_temp_files()

    def _get_pdb_file(self) -> str:
        """
        Get PDB file path (use existing or extract from TPR).

        Returns:
            Path to PDB file

        Raises:
            RuntimeError: If extraction from TPR fails
        """
        if self.input_type == 'pdb':
            # Use existing PDB file
            logger.debug(f"Using existing PDB: {self.pdb_file}")
            return self.pdb_file

        elif self.input_type == 'tpr':
            # Extract PDB from TPR temporarily
            logger.info("Extracting PDB from TPR file...")

            # Create temporary PDB file
            temp_pdb = tempfile.NamedTemporaryFile(
                mode='w',
                suffix='.pdb',
                delete=False,
                dir=self.structure_source.parent
            )
            temp_pdb_path = temp_pdb.name
            temp_pdb.close()

            self._temp_files.append(Path(temp_pdb_path))

            # Extract structure using gmx trjconv
            cmd = [
                self.gmx, "trjconv",
                "-s", str(self.structure_source),
                "-dump", "0",
                "-o", temp_pdb_path
            ]

            try:
                result = subprocess.run(
                    cmd,
                    input="0\n",  # Select "System" group
                    text=True,
                    capture_output=True,
                    timeout=60
                )

                if result.returncode != 0:
                    raise RuntimeError(
                        f"Failed to extract PDB from TPR: {result.stderr}"
                    )

                self.pdb_file = temp_pdb_path
                logger.info(f"Extracted temporary PDB: {temp_pdb_path}")

                return temp_pdb_path

            except subprocess.TimeoutExpired:
                raise RuntimeError("Extraction from TPR timed out (>60s)")
            except Exception as e:
                raise RuntimeError(f"Failed to extract PDB: {str(e)}")

    def _extract_chain_sequences(self, pdb_file: str) -> Dict[str, Tuple[str, int]]:
        """
        Extract sequences and lengths for all chains in PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary: {chain_id: (sequence, num_residues)}
        """
        chains = {}
        current_chain = None
        residue_seq = {}  # {chain_id: {res_num: res_name}}

        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()

                    if chain_id not in residue_seq:
                        residue_seq[chain_id] = {}

                    # Store residue (using dict to avoid duplicates)
                    residue_seq[chain_id][res_num] = res_name

        # Convert residue names to single-letter sequence
        aa_mapping = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }

        for chain_id, residues in residue_seq.items():
            # Sort by residue number
            sorted_residues = sorted(residues.items())
            sequence = ''.join([
                aa_mapping.get(res_name, 'X')
                for _, res_name in sorted_residues
            ])
            chains[chain_id] = (sequence, len(residues))

        logger.info(f"Extracted {len(chains)} chains from PDB")
        for chain_id, (seq, length) in chains.items():
            logger.debug(f"  Chain {chain_id}: {length} residues")

        return chains

    def identify_chains(self, pdb_file: Optional[str] = None) -> Dict[str, str]:
        """
        Automatically identify chain types using ANARCI for Class I pHLA-TCR complexes.

        Strategy:
        1. Extract all chain sequences and lengths
        2. Identify peptide (shortest, 5-22 AA)
        3. Identify beta2m (~100 AA, 90-110 range)
        4. Use ANARCI to identify TCR chains from remaining 3 chains
        5. Longest non-TCR chain is MHC-alpha

        Args:
            pdb_file: Optional PDB file path (uses self.pdb_file if not provided)

        Returns:
            Dictionary: {
                'peptide': chain_id,
                'beta2m': chain_id,
                'mhc_alpha': chain_id,
                'tcr_alpha': chain_id,
                'tcr_beta': chain_id
            }

        Raises:
            ValueError: If chain identification fails
        """
        if pdb_file is None:
            pdb_file = self._get_pdb_file()

        logger.info("Starting automatic chain identification...")

        # Step 1: Extract all chains
        chains = self._extract_chain_sequences(pdb_file)

        if len(chains) != 5:
            raise ValueError(
                f"Expected 5 chains for pHLA-TCR complex, found {len(chains)}. "
                f"Chains: {list(chains.keys())}"
            )

        # Sort chains by length
        sorted_chains = sorted(chains.items(), key=lambda x: x[1][1])
        chain_ids = [c[0] for c in sorted_chains]
        chain_lengths = [c[1][1] for c in sorted_chains]

        logger.info(f"Chains sorted by length: {dict(zip(chain_ids, chain_lengths))}")

        chain_mapping = {}

        # Step 2: Identify peptide (shortest chain, 5-22 AA)
        peptide_chain = chain_ids[0]
        peptide_length = chain_lengths[0]

        if not (5 <= peptide_length <= 22):
            logger.warning(
                f"Shortest chain ({peptide_chain}) has {peptide_length} residues, "
                f"expected 5-22 for peptide"
            )

        chain_mapping['peptide'] = peptide_chain
        logger.info(f"Identified peptide: Chain {peptide_chain} ({peptide_length} AA)")

        # Step 3: Identify beta2m (~100 AA)
        beta2m_chain = None
        for i in range(1, len(sorted_chains)):
            chain_id = chain_ids[i]
            length = chain_lengths[i]
            if 90 <= length <= 110:
                beta2m_chain = chain_id
                chain_mapping['beta2m'] = beta2m_chain
                logger.info(f"Identified beta2m: Chain {beta2m_chain} ({length} AA)")
                break

        if beta2m_chain is None:
            raise ValueError(
                f"Could not identify beta2m chain. "
                f"Expected ~100 AA chain in range 90-110. "
                f"Found lengths: {chain_lengths}"
            )

        # Step 4: Identify TCR chains using ANARCI (remaining 3 chains)
        remaining_chains = [
            (cid, chains[cid])
            for cid in chain_ids
            if cid not in [peptide_chain, beta2m_chain]
        ]

        if len(remaining_chains) != 3:
            raise ValueError(
                f"Expected 3 remaining chains (MHC-α, TCR-α, TCR-β), "
                f"found {len(remaining_chains)}"
            )

        # Use ANARCI to identify TCR chains
        if self.use_anarci and self.anarci_wrapper:
            tcr_chains = []
            non_tcr_chain = None

            for chain_id, (sequence, length) in remaining_chains:
                try:
                    logger.info(f"Running ANARCI on chain {chain_id} ({length} AA)...")
                    result = self.anarci_wrapper.run_anarci(sequence, chain_type="TCR")

                    # Check if ANARCI successfully numbered the sequence
                    if result.get('cdr3_range') is not None:
                        tcr_chains.append((chain_id, length))
                        logger.info(f"  Chain {chain_id}: Identified as TCR (CDR3 detected)")
                    else:
                        non_tcr_chain = (chain_id, length)
                        logger.info(f"  Chain {chain_id}: Not TCR (no CDR3 detected)")

                except Exception as e:
                    logger.warning(f"  Chain {chain_id}: ANARCI failed ({e})")
                    # Assume it's not TCR if ANARCI fails
                    non_tcr_chain = (chain_id, length)

            # Validate results
            if len(tcr_chains) != 2:
                logger.warning(
                    f"ANARCI identified {len(tcr_chains)} TCR chains, expected 2. "
                    f"Falling back to length-based heuristic."
                )
                # Fallback: two shorter chains are TCR, longest is MHC-alpha
                sorted_remaining = sorted(remaining_chains, key=lambda x: x[1][1])
                tcr_chains = sorted_remaining[:2]
                non_tcr_chain = sorted_remaining[2]

            # Assign TCR-alpha and TCR-beta (shorter is usually alpha)
            tcr_sorted = sorted(tcr_chains, key=lambda x: x[1][1])
            chain_mapping['tcr_alpha'] = tcr_sorted[0][0]
            chain_mapping['tcr_beta'] = tcr_sorted[1][0]
            chain_mapping['mhc_alpha'] = non_tcr_chain[0]

            logger.info(f"Identified TCR-α: Chain {chain_mapping['tcr_alpha']} ({tcr_sorted[0][1]} AA)")
            logger.info(f"Identified TCR-β: Chain {chain_mapping['tcr_beta']} ({tcr_sorted[1][1]} AA)")
            logger.info(f"Identified MHC-α: Chain {chain_mapping['mhc_alpha']} ({non_tcr_chain[1]} AA)")

        else:
            # Fallback without ANARCI: use length heuristic
            logger.warning("ANARCI not available, using length-based heuristic")
            sorted_remaining = sorted(remaining_chains, key=lambda x: x[1][1])

            # Assume: shortest two are TCR, longest is MHC-alpha
            chain_mapping['tcr_alpha'] = sorted_remaining[0][0]
            chain_mapping['tcr_beta'] = sorted_remaining[1][0]
            chain_mapping['mhc_alpha'] = sorted_remaining[2][0]

            logger.info(f"Heuristic TCR-α: Chain {chain_mapping['tcr_alpha']} ({sorted_remaining[0][1][1]} AA)")
            logger.info(f"Heuristic TCR-β: Chain {chain_mapping['tcr_beta']} ({sorted_remaining[1][1][1]} AA)")
            logger.info(f"Heuristic MHC-α: Chain {chain_mapping['mhc_alpha']} ({sorted_remaining[2][1][1]} AA)")

        logger.info("Chain identification completed successfully")
        logger.info(f"Final mapping: {chain_mapping}")

        return chain_mapping

    def _read_chain_info(self, pdb_file: str, chain_id: str) -> Dict:
        """
        Read chain information from PDB file.

        Args:
            pdb_file: Path to PDB file
            chain_id: Chain identifier to read

        Returns:
            Dictionary containing:
            {
                'chain_id': str,
                'residues': List[int],
                'atoms': List[int],
                'num_residues': int,
                'num_atoms': int
            }

        Raises:
            ValueError: If chain not found in PDB
        """
        residues = set()
        atoms = []

        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_chain = line[21]

                    if atom_chain == chain_id:
                        try:
                            atom_num = int(line[6:11].strip())
                            res_num = int(line[22:26].strip())

                            atoms.append(atom_num)
                            residues.add(res_num)
                        except ValueError:
                            # Skip malformed lines
                            continue

        if not atoms:
            raise ValueError(
                f"Chain '{chain_id}' not found in PDB file: {pdb_file}\n"
                f"Please ensure the PDB has been standardized with "
                f"PDBChainStandardizer."
            )

        return {
            'chain_id': chain_id,
            'residues': sorted(residues),
            'atoms': atoms,
            'num_residues': len(residues),
            'num_atoms': len(atoms)
        }

    def _create_index_file(
        self,
        pdb_file: str,
        peptide_info: Dict,
        output_dir: str,
        index_name: str
    ) -> str:
        """
        Create GROMACS index file with peptide group.

        Args:
            pdb_file: PDB file path
            peptide_info: Peptide chain information dictionary
            output_dir: Output directory
            index_name: Index file name

        Returns:
            Path to created index file
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        index_file = output_path / index_name

        # Try to generate base index groups
        base_index_content = self._generate_base_index_content(pdb_file)

        # Write index file
        with open(index_file, 'w') as f:
            # Write base groups if available
            if base_index_content:
                f.write(base_index_content)
                f.write('\n')

            # Write peptide group
            f.write(f"[ Peptide_Chain_{peptide_info['chain_id']} ]\n")

            # Format atom numbers (15 per line)
            atoms = peptide_info['atoms']
            for i in range(0, len(atoms), 15):
                line_atoms = atoms[i:i+15]
                f.write(' '.join(f"{a:6d}" for a in line_atoms) + '\n')

        logger.info(
            f"Created index file with peptide group:\n"
            f"  File: {index_file}\n"
            f"  Chain: {peptide_info['chain_id']}\n"
            f"  Residues: {peptide_info['num_residues']}\n"
            f"  Atoms: {peptide_info['num_atoms']}"
        )

        return str(index_file)

    def _generate_base_index_content(self, structure_file: str) -> Optional[str]:
        """
        Generate base index groups (System, Protein, etc.) using gmx make_ndx.

        Args:
            structure_file: PDB file path

        Returns:
            Index file content as string, or None if generation fails
        """
        try:
            # Create temporary index file
            temp_index = tempfile.NamedTemporaryFile(
                mode='w',
                suffix='.ndx',
                delete=False
            )
            temp_index_path = temp_index.name
            temp_index.close()

            self._temp_files.append(Path(temp_index_path))

            # Generate index using gmx make_ndx
            cmd = [self.gmx, "make_ndx", "-f", structure_file, "-o", temp_index_path]

            result = subprocess.run(
                cmd,
                input="q\n",  # Quit immediately (just generate default groups)
                text=True,
                capture_output=True,
                timeout=30
            )

            if result.returncode == 0 and Path(temp_index_path).exists():
                with open(temp_index_path) as f:
                    content = f.read()
                logger.debug("Generated base index groups using gmx make_ndx")
                return content
            else:
                logger.warning(
                    "Failed to generate base index groups, "
                    "index will contain only peptide group"
                )
                return None

        except Exception as e:
            logger.warning(f"Error generating base index: {e}")
            return None

    def _cleanup_temp_files(self):
        """Clean up temporary files created during processing."""
        if not self._temp_files:
            return

        for temp_file in self._temp_files:
            try:
                if temp_file.exists():
                    temp_file.unlink()
                    logger.debug(f"Cleaned up temporary file: {temp_file}")
            except Exception as e:
                logger.warning(f"Failed to clean up {temp_file}: {e}")

        self._temp_files.clear()


# Convenience function
def generate_peptide_index_from_pdb(
    pdb_file: str,
    output_dir: str = "./",
    peptide_chain_id: str = 'C'
) -> str:
    """
    Convenience function to generate peptide index from PDB file.

    Args:
        pdb_file: Path to standardized PDB file
        output_dir: Output directory for index file
        peptide_chain_id: Chain ID of peptide (default: 'C')

    Returns:
        Path to generated index file

    Example:
        >>> index_file = generate_peptide_index_from_pdb(
        ...     "complex_std_protein.pdb",
        ...     output_dir="./"
        ... )
    """
    generator = ChainBasedIndexGenerator(pdb_file)
    return generator.generate_peptide_index(output_dir, peptide_chain_id)


def generate_peptide_index_from_tpr(
    tpr_file: str,
    output_dir: str = "./",
    peptide_chain_id: str = 'C'
) -> str:
    """
    Convenience function to generate peptide index from TPR file.

    Args:
        tpr_file: Path to GROMACS TPR file
        output_dir: Output directory for index file
        peptide_chain_id: Chain ID of peptide (default: 'C')

    Returns:
        Path to generated index file

    Example:
        >>> index_file = generate_peptide_index_from_tpr(
        ...     "md.tpr",
        ...     output_dir="./"
        ... )
    """
    generator = ChainBasedIndexGenerator(tpr_file)
    return generator.generate_peptide_index(output_dir, peptide_chain_id)
