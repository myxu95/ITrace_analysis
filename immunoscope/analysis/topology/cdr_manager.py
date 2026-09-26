#!/usr/bin/env python3
"""
CDR Recognition and Management Module for pHLA-TCR Complexes

This module provides comprehensive CDR (Complementarity Determining Region)
recognition capabilities with ANARCI integration, GROMACS index file generation,
and metadata management.

Author: ImmunoScope Development Team
Date: 2025-12-26
"""

import subprocess
import re
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
import MDAnalysis as mda

logger = logging.getLogger(__name__)


class ANARCIWrapper:
    """
    Wrapper for ANARCI (Antibody Numbering and Receptor ClassIfication) tool.

    Provides three-tier CDR detection strategy:
    1. ANARCI command-line execution (primary)
    2. ANARCI Python package import (secondary)
    3. Regex-based CDR3 fallback (tertiary)
    """

    def __init__(self, allow_fallback: bool = True, numbering_scheme: str = "imgt"):
        """
        Initialize ANARCI wrapper.

        Args:
            allow_fallback: Allow regex-based CDR3 detection if ANARCI unavailable
            numbering_scheme: ANARCI numbering scheme (imgt/kabat/chothia)
        """
        self.allow_fallback = allow_fallback
        self.numbering_scheme = numbering_scheme

        # Fix PATH to include conda environment bin directory
        self._setup_environment()

        self.anarci_available = self.validate_installation()

        if not self.anarci_available:
            if allow_fallback:
                logger.warning("ANARCI not available, will use regex fallback for CDR3")
            else:
                raise RuntimeError(
                    "ANARCI not installed and fallback disabled. "
                    "Install ANARCI: pip install anarci"
                )

    def _setup_environment(self):
        """Setup environment variables for ANARCI execution."""
        import os
        import sys

        # Add Python executable's directory to PATH (for hmmscan, ANARCI, etc.)
        python_bin_dir = Path(sys.executable).parent
        current_path = os.environ.get('PATH', '')

        if str(python_bin_dir) not in current_path:
            os.environ['PATH'] = str(python_bin_dir) + os.pathsep + current_path
            logger.debug(f"Added {python_bin_dir} to PATH for ANARCI dependencies")

    def validate_installation(self) -> bool:
        """
        Validate ANARCI installation via multiple methods.

        Returns:
            True if ANARCI is available
        """
        # Method 1: Try command-line ANARCI
        try:
            result = subprocess.run(
                ["ANARCI", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.info("ANARCI command-line tool detected")
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        # Method 2: Try Python package import
        try:
            import anarci
            logger.info("ANARCI Python package detected")
            return True
        except ImportError:
            pass

        logger.warning("ANARCI not found via command-line or Python package")
        return False

    def run_anarci(self, sequence: str, chain_type: str = "TCR") -> Dict[str, Any]:
        """
        Run ANARCI on amino acid sequence.

        Args:
            sequence: Amino acid sequence (single-letter codes)
            chain_type: Chain type ('TCR' or 'Ig')

        Returns:
            Dictionary with CDR regions and numbering information
        """
        if not self.anarci_available:
            # Use fallback for CDR3 only
            logger.info("Using regex fallback for CDR3 detection")
            return self._fallback_cdr3_only(sequence, chain_type)

        try:
            # Try Python package first
            return self._run_anarci_python(sequence, chain_type)
        except Exception as e:
            logger.warning(f"ANARCI Python execution failed: {e}")

            # Try command-line as fallback
            try:
                return self._run_anarci_command(sequence, chain_type)
            except Exception as e2:
                logger.error(f"ANARCI command-line execution failed: {e2}")

                # Use regex fallback
                if self.allow_fallback:
                    logger.info("Falling back to regex CDR3 detection")
                    return self._fallback_cdr3_only(sequence, chain_type)
                else:
                    raise RuntimeError(f"ANARCI execution failed: {e}, {e2}")

    def _run_anarci_python(self, sequence: str, chain_type: str) -> Dict[str, Any]:
        """Run ANARCI via Python package."""
        try:
            from anarci import anarci, number

            # Run ANARCI
            # Use 'human' as species for TCR sequences (not 'TCR')
            results = anarci(
                [("sequence", sequence)],
                scheme=self.numbering_scheme,
                allowed_species=['human']
            )

            if not results or not results[0] or not results[0][0] or not results[0][0][0]:
                raise ValueError("ANARCI returned no results")

            # ANARCI returns: ([[(numbering, start, end)]], [hit_tables], [details])
            numbering, start_pos, end_pos = results[0][0][0]
            hit_table = results[2][0] if len(results) > 2 and results[2] else None

            # Parse CDR regions from numbering
            cdr_regions = self._parse_anarci_numbering(numbering, None)

            # Determine chain type from hit table
            chain_type_detected = self._determine_chain_type_from_hits(hit_table)

            return {
                'method': 'anarci_python',
                'numbering_scheme': self.numbering_scheme,
                'cdr1_range': cdr_regions.get('cdr1'),
                'cdr2_range': cdr_regions.get('cdr2'),
                'cdr3_range': cdr_regions.get('cdr3'),
                'numbering': numbering,
                'start_pos': start_pos,
                'end_pos': end_pos,
                'chain_type': chain_type_detected,
                'hit_table': hit_table,
                'raw_output': str(results)
            }

        except ImportError:
            raise RuntimeError("anarci package not available")

    def _run_anarci_command(self, sequence: str, chain_type: str) -> Dict[str, Any]:
        """Run ANARCI via command-line."""
        # Create temporary FASTA file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">sequence\n{sequence}\n")
            fasta_file = f.name

        try:
            # Run ANARCI
            cmd = [
                "ANARCI",
                "-i", fasta_file,
                "--scheme", self.numbering_scheme,
                "--assign_germline"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode != 0:
                raise RuntimeError(f"ANARCI command failed: {result.stderr}")

            # Parse output
            cdr_regions = self._parse_anarci_output(result.stdout)

            return {
                'method': 'anarci_command',
                'numbering_scheme': self.numbering_scheme,
                'cdr1_range': cdr_regions.get('cdr1'),
                'cdr2_range': cdr_regions.get('cdr2'),
                'cdr3_range': cdr_regions.get('cdr3'),
                'raw_output': result.stdout
            }

        finally:
            # Cleanup
            Path(fasta_file).unlink(missing_ok=True)

    def _determine_chain_type_from_hits(self, hit_table: List) -> str:
        """
        Determine TCR chain type from ANARCI hit table.

        Args:
            hit_table: ANARCI hit table with format:
                [['id', 'description', 'evalue', ...],  # Header
                 ['human_A', '', 1e-38, ...],           # Best hit
                 ['mouse_B', '', 2e-37, ...]]           # Second hit

        Returns:
            Chain type: 'TCR_alpha', 'TCR_beta', 'TCR_delta', 'TCR_gamma', or 'unknown'
        """
        if not hit_table or len(hit_table) < 2:
            return 'unknown'

        # Get best hit (second row, first column)
        best_hit_id = hit_table[1][0] if len(hit_table[1]) > 0 else ''

        # Parse chain type from hit ID (e.g., 'human_A' -> A, 'mouse_B' -> B)
        chain_letter = best_hit_id.split('_')[-1] if '_' in best_hit_id else ''

        chain_type_mapping = {
            'A': 'TCR_alpha',
            'B': 'TCR_beta',
            'D': 'TCR_delta',
            'G': 'TCR_gamma'
        }

        detected_type = chain_type_mapping.get(chain_letter, 'unknown')

        logger.debug(f"Hit ID '{best_hit_id}' -> chain type '{detected_type}'")

        return detected_type

    def _parse_anarci_numbering(self, numbering: List, alignment_details: Any) -> Dict:
        """
        Parse ANARCI numbering to extract CDR regions.

        IMGT CDR definitions:
        - CDR1: 27-38
        - CDR2: 56-65
        - CDR3: 105-117

        Args:
            numbering: List of ((position, insertion), amino_acid) tuples
            alignment_details: Not used (kept for compatibility)
        """
        cdr_regions = {}

        # Extract positions from numbering.
        # numbering format: [((pos, insertion), amino_acid), ...]
        #
        # IMPORTANT:
        # The numbering list can include alignment gaps. The enumerate index is
        # therefore NOT guaranteed to match the ungapped sequence index. We must
        # track the original sequence index by incrementing only when ANARCI
        # emits a real amino acid.
        positions = []
        seq_idx = 0
        for item in numbering:
            if isinstance(item, tuple) and len(item) == 2:
                pos_info, aa = item
                if isinstance(pos_info, tuple) and len(pos_info) == 2:
                    pos, insertion = pos_info
                    if pos is not None and aa != '-':  # Skip gaps
                        positions.append((seq_idx, pos, insertion))
                        seq_idx += 1

        if not positions:
            return cdr_regions

        # IMGT CDR ranges
        imgt_cdr_ranges = {
            'cdr1': (27, 38),
            'cdr2': (56, 65),
            'cdr3': (105, 117)
        }

        for cdr_name, (start_pos, end_pos) in imgt_cdr_ranges.items():
            # Find sequence indices corresponding to IMGT positions
            indices = [
                seq_i for seq_i, pos, ins in positions
                if start_pos <= pos <= end_pos
            ]
            if indices:
                # Return as (start_index, end_index) in original sequence
                cdr_regions[cdr_name] = (min(indices), max(indices) + 1)

        return cdr_regions

    def _parse_anarci_output(self, output: str) -> Dict:
        """Parse ANARCI command-line output."""
        # This is a simplified parser - actual implementation would be more robust
        cdr_regions = {}

        # Look for CDR regions in output
        # ANARCI output format varies, this is a placeholder
        # In practice, you'd parse the numbered output more carefully

        return cdr_regions

    def _fallback_cdr3_only(self, sequence: str, chain_type: str) -> Dict[str, Any]:
        """
        Fallback regex-based CDR3 detection.

        TCR CDR3 pattern: C-X(8-17)-[FW]
        - Starts with Cysteine (C)
        - Followed by 8-17 variable amino acids
        - Ends with Phenylalanine (F) or Tryptophan (W)
        """
        logger.info("Using regex fallback for CDR3 detection")

        # TCR CDR3 pattern
        pattern = r'C.{8,17}[FW]'

        match = re.search(pattern, sequence)

        if match:
            start = match.start()
            end = match.end()
            cdr3_seq = match.group()

            logger.info(f"CDR3 detected via regex: {cdr3_seq} (pos {start}-{end})")

            return {
                'method': 'regex_fallback',
                'numbering_scheme': 'none',
                'cdr1_range': None,
                'cdr2_range': None,
                'cdr3_range': (start, end),
                'cdr3_sequence': cdr3_seq,
                'raw_output': f"Regex match: {cdr3_seq}"
            }
        else:
            logger.warning("CDR3 not detected via regex fallback")
            return {
                'method': 'regex_fallback',
                'numbering_scheme': 'none',
                'cdr1_range': None,
                'cdr2_range': None,
                'cdr3_range': None,
                'raw_output': "No CDR3 match"
            }


class CDRIndexGenerator:
    """
    Generator for GROMACS-compatible index files (.ndx) with CDR regions.

    Creates index groups for CDR regions compatible with GROMACS tools.
    """

    def __init__(self):
        """Initialize CDR index generator."""
        self.groups = {}

    def add_cdr_group(self,
                     group_name: str,
                     atom_indices: List[int],
                     ca_indices: Optional[List[int]] = None):
        """
        Add CDR group to index.

        Args:
            group_name: Group name (e.g., 'CDR3_alpha')
            atom_indices: All atom indices for the CDR
            ca_indices: CA atom indices (optional, creates separate group)
        """
        # Add all-atom group
        self.groups[group_name] = atom_indices

        # Add CA-only group if provided
        if ca_indices:
            ca_group_name = f"{group_name}_CA"
            self.groups[ca_group_name] = ca_indices

    def generate_index_file(self, output_file: str):
        """
        Generate GROMACS .ndx index file.

        Args:
            output_file: Output file path
        """
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w') as f:
            for group_name, indices in self.groups.items():
                f.write(f"[ {group_name} ]\n")

                # Write indices (10 per line, GROMACS convention)
                for i in range(0, len(indices), 10):
                    line_indices = indices[i:i+10]
                    f.write(" ".join(map(str, line_indices)) + "\n")

                f.write("\n")

        logger.info(f"Index file generated: {output_file}")
        logger.info(f"  Groups: {list(self.groups.keys())}")

    def clear(self):
        """Clear all groups."""
        self.groups = {}


class CDRMetadataManager:
    """
    Manager for CDR metadata in JSON format.

    Handles saving and loading of CDR information with validation.
    """

    def __init__(self):
        """Initialize metadata manager."""
        pass

    def save_metadata(self,
                     cdr_data: Dict[str, Any],
                     output_file: str,
                     task_name: str,
                     topology_file: str):
        """
        Save CDR metadata to JSON file.

        Args:
            cdr_data: CDR detection results
            output_file: Output JSON file path
            task_name: Task identifier
            topology_file: Topology file path
        """
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Build metadata structure
        metadata = {
            'task_name': task_name,
            'topology_file': str(topology_file),
            'timestamp': datetime.now().isoformat(),
            'anarci_method': cdr_data.get('detection_method', 'unknown'),
            'numbering_scheme': cdr_data.get('numbering_scheme', 'none'),
            'chains': cdr_data.get('chains', {})
        }

        # Add index file references
        if 'index_files' in cdr_data:
            metadata['index_files'] = cdr_data['index_files']

        # Save to file
        with open(output_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"Metadata saved: {output_file}")

    def load_metadata(self, metadata_file: str) -> Dict[str, Any]:
        """
        Load CDR metadata from JSON file.

        Args:
            metadata_file: Metadata JSON file path

        Returns:
            Metadata dictionary
        """
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        logger.info(f"Metadata loaded: {metadata_file}")
        return metadata

    def validate_schema(self, metadata: Dict[str, Any]) -> bool:
        """
        Validate metadata schema.

        Args:
            metadata: Metadata dictionary

        Returns:
            True if valid
        """
        required_fields = ['task_name', 'topology_file', 'timestamp', 'chains']

        for field in required_fields:
            if field not in metadata:
                logger.error(f"Missing required field: {field}")
                return False

        return True


class CDRManager:
    """
    Main coordinator for CDR recognition and management.

    Integrates ANARCI detection, index file generation, and metadata management.
    """

    def __init__(self,
                 topology_file: str,
                 allow_fallback: bool = True,
                 numbering_scheme: str = "imgt"):
        """
        Initialize CDR manager.

        Args:
            topology_file: Path to topology file (.tpr/.pdb/.gro)
            allow_fallback: Allow regex fallback if ANARCI unavailable
            numbering_scheme: ANARCI numbering scheme
        """
        self.topology_file = Path(topology_file)

        if not self.topology_file.exists():
            raise FileNotFoundError(f"Topology file not found: {topology_file}")

        # Initialize components
        self.anarci = ANARCIWrapper(
            allow_fallback=allow_fallback,
            numbering_scheme=numbering_scheme
        )
        self.index_generator = CDRIndexGenerator()
        self.metadata_manager = CDRMetadataManager()

        # Load topology
        try:
            self.universe = mda.Universe(str(self.topology_file))
            logger.info(f"Topology loaded: {self.topology_file}")
            logger.info(f"Total atoms: {self.universe.atoms.n_atoms}")
        except Exception as e:
            raise RuntimeError(f"Failed to load topology: {e}")

        # Storage for CDR data
        self.cdr_data = {
            'topology_file': str(self.topology_file),
            'detection_method': None,
            'numbering_scheme': numbering_scheme,
            'chains': {},
            'index_files': {}
        }

    def detect_all_cdr_regions(self, chains: Dict[str, str]) -> Dict[str, Any]:
        """
        Detect CDR regions for all specified chains.

        Args:
            chains: Dictionary mapping chain names to selection strings
                   e.g., {'TCR_alpha': 'chainID D', 'TCR_beta': 'chainID E'}

        Returns:
            Dictionary with CDR data for each chain
        """
        logger.info(f"Detecting CDR regions for {len(chains)} chains")

        for chain_name, selection in chains.items():
            logger.info(f"\nProcessing chain: {chain_name}")
            logger.info(f"  Selection: {selection}")

            try:
                chain_data = self._detect_chain_cdrs(chain_name, selection)
                self.cdr_data['chains'][chain_name] = chain_data
            except Exception as e:
                logger.error(f"Failed to process {chain_name}: {e}")
                self.cdr_data['chains'][chain_name] = {
                    'error': str(e),
                    'cdr3_detected': False
                }

        return self.cdr_data

    def _detect_chain_cdrs(self, chain_name: str, selection: str) -> Dict[str, Any]:
        """
        Detect CDR regions for a single chain.

        Args:
            chain_name: Chain name
            selection: MDAnalysis selection string

        Returns:
            Chain CDR data
        """
        # Select chain atoms
        atoms = self.universe.select_atoms(selection)

        if len(atoms) == 0:
            raise ValueError(f"No atoms found for selection: {selection}")

        logger.info(f"  Selected {len(atoms)} atoms")

        # Extract sequence
        sequence = self._extract_sequence(atoms)
        logger.info(f"  Sequence length: {len(sequence)} residues")
        logger.info(f"  Sequence: {sequence[:50]}..." if len(sequence) > 50 else f"  Sequence: {sequence}")

        # Detect CDR regions with ANARCI
        anarci_result = self.anarci.run_anarci(sequence, chain_type="TCR")

        # Store detection method
        if self.cdr_data['detection_method'] is None:
            self.cdr_data['detection_method'] = anarci_result['method']

        # Build chain data
        chain_data = {
            'chain_name': chain_name,
            'selection': selection,
            'sequence': sequence,
            'numbering_scheme': anarci_result['numbering_scheme']
        }

        # Process each CDR
        ca_atoms = atoms.select_atoms("name CA")

        for cdr_num in [1, 2, 3]:
            cdr_key = f'cdr{cdr_num}'
            cdr_range = anarci_result.get(f'{cdr_key}_range')

            if cdr_range:
                start_idx, end_idx = cdr_range

                # Get CDR sequence
                cdr_sequence = sequence[start_idx:end_idx]

                # Get CA atoms for CDR
                cdr_ca_atoms = ca_atoms[start_idx:end_idx]

                # Get residue range (PDB numbering)
                if len(cdr_ca_atoms) > 0:
                    residue_range = (
                        int(cdr_ca_atoms.resids[0]),
                        int(cdr_ca_atoms.resids[-1])
                    )

                    # Get all atoms for these residues
                    resid_selection = f"resid {residue_range[0]}:{residue_range[1]}"
                    cdr_all_atoms = atoms.select_atoms(resid_selection)

                    # Get atom indices (1-based for GROMACS)
                    all_atom_indices = [int(a.index) + 1 for a in cdr_all_atoms]
                    ca_atom_indices = [int(a.index) + 1 for a in cdr_ca_atoms]

                    chain_data[cdr_key] = {
                        'sequence': cdr_sequence,
                        'residue_range': residue_range,
                        'residue_count': len(cdr_ca_atoms),
                        'atom_count': len(all_atom_indices),
                        'ca_count': len(ca_atom_indices),
                        'atom_indices': all_atom_indices,
                        'ca_indices': ca_atom_indices,
                        'index_group_name': f"{cdr_key}_{chain_name}"
                    }

                    logger.info(f"  {cdr_key.upper()}: {cdr_sequence}")
                    logger.info(f"    Residues: {residue_range[0]}-{residue_range[1]}")
                    logger.info(f"    Atoms: {len(all_atom_indices)} total, {len(ca_atom_indices)} CA")
                else:
                    logger.warning(f"  {cdr_key.upper()}: No CA atoms found")
            else:
                logger.info(f"  {cdr_key.upper()}: Not detected")

        # Mark if CDR3 was detected
        chain_data['cdr3_detected'] = 'cdr3' in chain_data

        return chain_data

    def _extract_sequence(self, atoms: mda.AtomGroup) -> str:
        """
        Extract amino acid sequence from atoms.

        Args:
            atoms: MDAnalysis atom group

        Returns:
            Single-letter amino acid sequence
        """
        # Three-letter to one-letter mapping
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        # Select CA atoms (one per residue)
        ca_atoms = atoms.select_atoms("name CA")

        # Convert to single-letter sequence
        sequence = ''.join([aa_map.get(res, 'X') for res in ca_atoms.resnames])

        return sequence

    def generate_index_file(self,
                           output_file: str,
                           include_cdrs: List[int] = [3],
                           chains: Optional[List[str]] = None):
        """
        Generate GROMACS index file for CDR regions.

        Args:
            output_file: Output .ndx file path
            include_cdrs: List of CDR numbers to include (e.g., [1, 2, 3])
            chains: List of chain names to include (None = all)
        """
        logger.info(f"Generating index file: {output_file}")
        logger.info(f"  Including CDR{include_cdrs}")

        self.index_generator.clear()

        # Add groups for each chain and CDR
        for chain_name, chain_data in self.cdr_data['chains'].items():
            if chains and chain_name not in chains:
                continue

            if 'error' in chain_data:
                logger.warning(f"Skipping {chain_name}: {chain_data['error']}")
                continue

            for cdr_num in include_cdrs:
                cdr_key = f'cdr{cdr_num}'
                if cdr_key in chain_data:
                    cdr_info = chain_data[cdr_key]
                    group_name = f"CDR{cdr_num}_{chain_name}"

                    self.index_generator.add_cdr_group(
                        group_name=group_name,
                        atom_indices=cdr_info['atom_indices'],
                        ca_indices=cdr_info['ca_indices']
                    )

        self.index_generator.generate_index_file(output_file)

        # Store in metadata
        self.cdr_data['index_files'][f"cdr{'_'.join(map(str, include_cdrs))}"] = str(output_file)

    def save_metadata(self, output_file: str, task_name: str):
        """
        Save CDR metadata to JSON file.

        Args:
            output_file: Output JSON file path
            task_name: Task identifier
        """
        self.metadata_manager.save_metadata(
            cdr_data=self.cdr_data,
            output_file=output_file,
            task_name=task_name,
            topology_file=str(self.topology_file)
        )

    def get_cdr_summary(self) -> Dict[str, Any]:
        """
        Get summary of detected CDR regions.

        Returns:
            Summary dictionary
        """
        summary = {
            'total_chains': len(self.cdr_data['chains']),
            'detection_method': self.cdr_data['detection_method'],
            'chains': {}
        }

        for chain_name, chain_data in self.cdr_data['chains'].items():
            if 'error' in chain_data:
                summary['chains'][chain_name] = {'status': 'failed'}
            else:
                cdrs_detected = [f'CDR{i}' for i in [1, 2, 3] if f'cdr{i}' in chain_data]
                summary['chains'][chain_name] = {
                    'status': 'success',
                    'cdrs_detected': cdrs_detected,
                    'cdr3_sequence': chain_data.get('cdr3', {}).get('sequence', 'N/A')
                }

        return summary
