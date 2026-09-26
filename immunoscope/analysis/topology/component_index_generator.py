"""
Index Generator for GROMACS Analysis

This module provides a unified interface for generating GROMACS index files
for different protein components in pHLA-TCR complex simulations.

Supported components:
- HLA: MHC alpha and beta chains
- pHLA: MHC + peptide complex
- peptide: Antigenic peptide
- TCR: T-cell receptor alpha and beta chains
- TCR_alpha: TCR alpha chain only
- TCR_beta: TCR beta chain only
- CDR3_alpha: CDR3 region of TCR alpha (requires sequence)
- CDR3_beta: CDR3 region of TCR beta (requires sequence)
- custom: User-defined selections

Standard chain ordering (after standardization):
- Chain A: HLA-alpha (longest, ~270-280 residues)
- Chain B: HLA-beta/beta2m (~99-100 residues)
- Chain C: Peptide (shortest, 5-22 residues)
- Chain D: TCR-alpha (115-215 residues)
- Chain E: TCR-beta (209-252 residues)

Author: ImmunoScope Development Team
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict, List, Union
import MDAnalysis as mda

from ..structure import PDBChainStandardizer

logger = logging.getLogger(__name__)


class ComponentIndexGenerator:
    """
    Generate GROMACS index files for different protein components.

    This class provides a unified interface to create index files for
    various structural components of pHLA-TCR complexes.
    """

    COMPONENT_DEFINITIONS = {
        'HLA': {
            'description': 'MHC alpha and beta chains',
            'chains': ['A', 'B'],
            'command_template': 'chain {chains}\nname {group_id} HLA\n'
        },
        'HLA_alpha': {
            'description': 'MHC alpha chain',
            'chains': ['A'],
            'command_template': 'chain {chains}\nname {group_id} HLA_alpha\n'
        },
        'HLA_beta': {
            'description': 'MHC beta chain (beta2-microglobulin)',
            'chains': ['B'],
            'command_template': 'chain {chains}\nname {group_id} HLA_beta\n'
        },
        'pHLA': {
            'description': 'pMHC complex (MHC + peptide)',
            'chains': ['A', 'B', 'C'],
            'command_template': 'chain {chains}\nname {group_id} pHLA\n'
        },
        'peptide': {
            'description': 'Antigenic peptide',
            'chains': ['C'],
            'command_template': 'chain {chains}\nname {group_id} peptide\n'
        },
        'TCR': {
            'description': 'T-cell receptor alpha and beta chains',
            'chains': ['D', 'E'],
            'command_template': 'chain {chains}\nname {group_id} TCR\n'
        },
        'TCR_alpha': {
            'description': 'TCR alpha chain',
            'chains': ['D'],
            'command_template': 'chain {chains}\nname {group_id} TCR_alpha\n'
        },
        'TCR_beta': {
            'description': 'TCR beta chain',
            'chains': ['E'],
            'command_template': 'chain {chains}\nname {group_id} TCR_beta\n'
        }
    }

    def __init__(self,
                 topology: str,
                 gmx_executable: str = "gmx",
                 auto_standardize: bool = True,
                 standardized_pdb: Optional[str] = None):
        """
        Initialize ComponentIndexGenerator.

        Args:
            topology: Path to topology file (.tpr, .gro, .pdb)
            gmx_executable: GROMACS executable command
            auto_standardize: If True, automatically standardize chain ordering
            standardized_pdb: Path to pre-standardized PDB file (if available)
        """
        self.topology = Path(topology)
        self.gmx = gmx_executable
        self.auto_standardize = auto_standardize
        self.standardized_pdb = None
        self.chain_standardizer = PDBChainStandardizer()

        if not self.topology.exists():
            raise FileNotFoundError(f"Topology file not found: {self.topology}")

        logger.info(f"ComponentIndexGenerator initialized with topology: {self.topology.name}")

        if auto_standardize:
            if standardized_pdb and Path(standardized_pdb).exists():
                self.standardized_pdb = Path(standardized_pdb)
                logger.info(f"Using pre-standardized PDB: {self.standardized_pdb.name}")
            else:
                self.standardized_pdb = self._prepare_standardized_pdb()
                if self.standardized_pdb:
                    logger.info(f"Chain standardization complete: {self.standardized_pdb.name}")
                else:
                    logger.warning("Chain standardization failed, using original topology")

    def _prepare_standardized_pdb(self) -> Optional[Path]:
        """
        Prepare a standardized PDB file for index generation.

        Returns:
            Path to standardized PDB file, or None if failed
        """
        try:
            if self.topology.suffix.lower() == '.pdb':
                pdb_file = self.topology
            else:
                pdb_file = self._convert_to_pdb()
                if not pdb_file:
                    return None

            standardized_file = pdb_file.parent / f"{pdb_file.stem}_standardized.pdb"

            result = self.chain_standardizer.process_single(
                input_pdb=pdb_file,
                output_pdb=standardized_file,
                task_name=self.topology.stem,
                skip_if_standard=True
            )

            if result.status in ['OK', 'ALREADY_STANDARD']:
                if result.status == 'OK':
                    return standardized_file
                else:
                    return pdb_file
            else:
                logger.warning(f"Chain standardization failed: {result.error_message}")
                return None

        except Exception as e:
            logger.error(f"Failed to prepare standardized PDB: {e}")
            return None

    def _convert_to_pdb(self) -> Optional[Path]:
        """
        Convert TPR/GRO file to PDB format using gmx editconf.

        Returns:
            Path to converted PDB file, or None if failed
        """
        try:
            output_pdb = self.topology.parent / f"{self.topology.stem}_converted.pdb"

            cmd = [
                self.gmx, 'editconf',
                '-f', str(self.topology),
                '-o', str(output_pdb)
            ]

            logger.info(f"Converting {self.topology.suffix} to PDB format...")

            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60
            )

            if process.returncode != 0:
                logger.error(f"PDB conversion failed: {process.stderr}")
                return None

            if output_pdb.exists():
                logger.info(f"PDB conversion successful: {output_pdb.name}")
                return output_pdb
            else:
                logger.error("PDB file was not created")
                return None

        except subprocess.TimeoutExpired:
            logger.error("PDB conversion timed out")
            return None
        except Exception as e:
            logger.error(f"PDB conversion error: {e}")
            return None

    def _get_topology_for_indexing(self) -> Path:
        """
        Get the appropriate topology file for index generation.

        Returns:
            Path to topology file (standardized PDB if available, else original)
        """
        if self.standardized_pdb and self.standardized_pdb.exists():
            return self.standardized_pdb
        else:
            return self.topology

    def generate_component_index(self,
                                 component: str,
                                 output_file: Optional[str] = None,
                                 custom_chains: Optional[List[str]] = None) -> str:
        """
        Generate index file for a specific component.

        Args:
            component: Component name (HLA/pHLA/peptide/TCR/TCR_alpha/TCR_beta)
            output_file: Output index file path (default: {component}.ndx)
            custom_chains: Custom chain list (overrides default)

        Returns:
            Path to generated index file

        Raises:
            ValueError: If component is not recognized
            RuntimeError: If index generation fails
        """
        if component not in self.COMPONENT_DEFINITIONS:
            raise ValueError(
                f"Unknown component: {component}. "
                f"Available: {', '.join(self.COMPONENT_DEFINITIONS.keys())}"
            )

        component_def = self.COMPONENT_DEFINITIONS[component]
        chains = custom_chains if custom_chains else component_def['chains']

        if output_file is None:
            output_file = f"{component.lower()}.ndx"

        output_path = Path(output_file)

        logger.info(f"Generating {component} index (chains: {', '.join(chains)})")

        chain_str = ' '.join(chains)
        commands = f"chain {chain_str}\nname 20 {component}\nq\n"

        topology_file = self._get_topology_for_indexing()
        cmd = [self.gmx, 'make_ndx', '-f', str(topology_file), '-o', str(output_path)]

        try:
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate(input=commands, timeout=30)

            if process.returncode != 0:
                raise RuntimeError(f"gmx make_ndx failed: {stderr}")

            if not output_path.exists():
                raise RuntimeError("Index file was not created")

            with open(output_path, 'r') as f:
                content = f.read()
                if component not in content:
                    raise RuntimeError(f"{component} group not found in index file")

            logger.info(f"Index file created: {output_path}")
            return str(output_path)

        except subprocess.TimeoutExpired:
            raise RuntimeError("Index generation timed out")
        except Exception as e:
            logger.error(f"Failed to generate {component} index: {e}")
            raise

    def generate_multi_component_index(self,
                                       components: List[str],
                                       output_file: str = "combined.ndx") -> str:
        """
        Generate a single index file with multiple components.

        Args:
            components: List of component names
            output_file: Output index file path

        Returns:
            Path to generated index file

        Example:
            generator.generate_multi_component_index(
                ['pHLA', 'TCR', 'peptide'],
                'phla_tcr_peptide.ndx'
            )
        """
        output_path = Path(output_file)

        logger.info(f"Generating multi-component index with {len(components)} components")

        commands = []

        # Get the starting group number (after default groups + chain groups)
        # Default groups: 0-17, chain groups from PDB: 18-25 (chA, chB, chAB, chC, chDE, chD, chE)
        # So our custom groups start at 18 (we'll overwrite the chain groups with our named versions)
        group_id = 18

        for component in components:
            if component not in self.COMPONENT_DEFINITIONS:
                raise ValueError(f"Unknown component: {component}")

            component_def = self.COMPONENT_DEFINITIONS[component]
            chains = component_def['chains']

            # Use 'chain X' command for gmx make_ndx
            # Multiple chains: "chain A B" (space-separated)
            chain_str = ' '.join(chains)
            commands.append(f"chain {chain_str}")
            # IMPORTANT: Need blank line after "chain" to save the selection
            commands.append("")
            # Then rename the saved group
            commands.append(f"name {group_id} {component}")
            group_id += 1

        commands.append("q")
        command_string = '\n'.join(commands) + '\n'

        # Use standardized PDB for index generation (not TPR)
        topology_file = self._get_topology_for_indexing()
        logger.info(f"Using topology file for indexing: {topology_file.name}")

        cmd = [self.gmx, 'make_ndx', '-f', str(topology_file), '-o', str(output_path)]

        logger.debug(f"gmx make_ndx commands:\n{command_string}")

        try:
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate(input=command_string, timeout=30)

            logger.debug(f"gmx make_ndx stdout:\n{stdout}")
            if stderr:
                logger.debug(f"gmx make_ndx stderr:\n{stderr}")

            if process.returncode != 0:
                raise RuntimeError(f"gmx make_ndx failed: {stderr}")

            logger.info(f"Multi-component index created: {output_path}")
            logger.info(f"  Components: {', '.join(components)}")

            return str(output_path)

        except Exception as e:
            logger.error(f"Failed to generate multi-component index: {e}")
            raise

    def generate_cdr3_index(self,
                           chain: str,
                           cdr3_sequence: str,
                           output_file: Optional[str] = None,
                           ca_only: bool = True) -> Optional[str]:
        """
        Generate index file for CDR3 region based on sequence.

        Args:
            chain: TCR chain ('alpha' or 'beta')
            cdr3_sequence: CDR3 amino acid sequence
            output_file: Output index file path
            ca_only: If True, select only CA atoms (default: True)

        Returns:
            Path to generated index file, or None if sequence not found

        Note:
            This method uses MDAnalysis to locate CDR3 residues by sequence matching.
        """
        if chain.lower() not in ['alpha', 'beta']:
            raise ValueError("Chain must be 'alpha' or 'beta'")

        chain_id = 'D' if chain.lower() == 'alpha' else 'E'
        component_name = f"CDR3_TCR_{chain.lower()}"

        if output_file is None:
            output_file = f"cdr3_{chain.lower()}.ndx"

        output_path = Path(output_file)

        logger.info(f"Generating {component_name} index")
        logger.info(f"  CDR3 sequence: {cdr3_sequence}")

        try:
            topology_file = self._get_topology_for_indexing()
            u = mda.Universe(str(topology_file))

            chain_selection = f"chainID {chain_id} and protein"
            chain_atoms = u.select_atoms(f"{chain_selection} and name CA")

            chain_sequence = ''.join([
                mda.lib.util.convert_aa_code(r.resname)
                for r in chain_atoms.residues
            ])

            logger.info(f"  Full chain length: {len(chain_sequence)} residues")

            pos = chain_sequence.find(cdr3_sequence)

            if pos == -1:
                logger.warning(f"  CDR3 sequence not found in chain {chain_id}")
                return None

            start_idx = pos
            end_idx = pos + len(cdr3_sequence)

            logger.info(f"  CDR3 position: sequence index {start_idx}-{end_idx}")

            cdr3_residues = chain_atoms.residues[start_idx:end_idx]

            if ca_only:
                selection = f"chainID {chain_id} and protein and resid {cdr3_residues[0].resid}-{cdr3_residues[-1].resid} and name CA"
                component_name_full = f"{component_name}_CA"
            else:
                selection = f"chainID {chain_id} and protein and resid {cdr3_residues[0].resid}-{cdr3_residues[-1].resid}"
                component_name_full = component_name

            cdr3_atoms = u.select_atoms(selection)
            atom_indices = cdr3_atoms.indices + 1

            with open(output_path, 'w') as f:
                f.write(f"[ {component_name_full} ]\n")

                for i in range(0, len(atom_indices), 10):
                    line = ' '.join(str(idx) for idx in atom_indices[i:i+10])
                    f.write(f"{line}\n")

            logger.info(f"  CDR3 index created: {output_path}")
            logger.info(f"  Atom count: {len(atom_indices)}")

            return str(output_path)

        except Exception as e:
            logger.error(f"Failed to generate CDR3 index: {e}")
            return None

    def generate_combined_cdr3_index(self,
                                    cdr3_alpha_seq: Optional[str] = None,
                                    cdr3_beta_seq: Optional[str] = None,
                                    output_file: str = "cdr3_combined.ndx",
                                    ca_only: bool = True) -> str:
        """
        Generate combined CDR3 index for both alpha and beta chains.

        Args:
            cdr3_alpha_seq: CDR3 alpha sequence
            cdr3_beta_seq: CDR3 beta sequence
            output_file: Output index file path
            ca_only: If True, select only CA atoms

        Returns:
            Path to generated index file
        """
        output_path = Path(output_file)

        logger.info("Generating combined CDR3 index")

        temp_files = []

        try:
            if cdr3_alpha_seq:
                alpha_file = self.generate_cdr3_index(
                    'alpha', cdr3_alpha_seq, 'temp_cdr3_alpha.ndx', ca_only
                )
                if alpha_file:
                    temp_files.append(alpha_file)

            if cdr3_beta_seq:
                beta_file = self.generate_cdr3_index(
                    'beta', cdr3_beta_seq, 'temp_cdr3_beta.ndx', ca_only
                )
                if beta_file:
                    temp_files.append(beta_file)

            if not temp_files:
                raise RuntimeError("No CDR3 sequences were found")

            with open(output_path, 'w') as outf:
                for temp_file in temp_files:
                    with open(temp_file, 'r') as inf:
                        outf.write(inf.read())
                        outf.write('\n')

            for temp_file in temp_files:
                Path(temp_file).unlink()

            logger.info(f"Combined CDR3 index created: {output_path}")

            return str(output_path)

        except Exception as e:
            for temp_file in temp_files:
                if Path(temp_file).exists():
                    Path(temp_file).unlink()
            raise

    def merge_index_files(self,
                         index_files: List[str],
                         output_file: str = "merged.ndx") -> str:
        """
        Merge multiple index files into one.

        Args:
            index_files: List of index file paths
            output_file: Output merged index file

        Returns:
            Path to merged index file
        """
        output_path = Path(output_file)

        logger.info(f"Merging {len(index_files)} index files")

        with open(output_path, 'w') as outf:
            for idx_file in index_files:
                idx_path = Path(idx_file)
                if not idx_path.exists():
                    logger.warning(f"Index file not found: {idx_file}")
                    continue

                with open(idx_path, 'r') as inf:
                    outf.write(inf.read())
                    outf.write('\n')

        logger.info(f"Merged index file created: {output_path}")

        return str(output_path)

    @staticmethod
    def list_available_components() -> Dict[str, str]:
        """
        List all available component definitions.

        Returns:
            Dictionary of component names and descriptions
        """
        return {
            name: info['description']
            for name, info in IndexGenerator.COMPONENT_DEFINITIONS.items()
        }


# Temporary internal alias while downstream modules are migrated.
IndexGenerator = ComponentIndexGenerator
