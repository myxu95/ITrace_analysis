"""
PDB Sequence Extractor Module

Extracts amino acid sequences from PDB files and exports to CSV format.
Essential for chain identification and ANARCI analysis in pHLA-TCR complexes.

Author: ImmunoScope Development Team
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import csv
import MDAnalysis as mda

logger = logging.getLogger(__name__)


class PDBSequenceExtractor:
    """
    Extract amino acid sequences from PDB structures.

    Provides sequence extraction for protein chains with support for:
    - Standard and modified amino acids
    - Multiple chain systems
    - CSV export format
    - Integration with ANARCI for chain identification

    Attributes:
        standard_residues: Standard 20 amino acids (3-letter codes)
        three_to_one: Mapping from 3-letter to 1-letter codes
    """

    # Standard amino acid mapping
    THREE_TO_ONE = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # Modified residues
        'MSE': 'M',  # Selenomethionine
        'HSD': 'H', 'HSE': 'H', 'HSP': 'H',  # Histidine variants
        'HIE': 'H', 'HID': 'H', 'HIP': 'H',
        'CYX': 'C', 'CYM': 'C',  # Cysteine variants
        'ASH': 'D', 'GLH': 'E',  # Protonated variants
        'LYN': 'K',  # Deprotonated lysine
    }

    def __init__(self):
        """Initialize PDB sequence extractor."""
        self.three_to_one = self.THREE_TO_ONE

    def extract_sequences_from_pdb(
        self,
        pdb_file: str
    ) -> Dict[str, Dict[str, any]]:
        """
        Extract sequences from all protein chains in PDB file.

        Uses MDAnalysis to parse PDB structure and extract sequences
        for each chain, along with metadata like residue count.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary mapping chain_id to chain information:
            {
                'A': {
                    'sequence': 'MKTAYIAK...',
                    'length': 275,
                    'residue_ids': [1, 2, 3, ...],
                    'residue_names': ['MET', 'LYS', ...]
                },
                ...
            }

        Example:
            >>> extractor = PDBSequenceExtractor()
            >>> chains = extractor.extract_sequences_from_pdb("1ao7.pdb")
            >>> print(f"Chain A: {chains['A']['length']} residues")
        """
        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        try:
            # Load structure with MDAnalysis
            u = mda.Universe(str(pdb_path))

            # Extract chains
            chains_data = {}

            # Get all protein atoms
            protein = u.select_atoms("protein")

            if len(protein) == 0:
                logger.warning(f"No protein atoms found in {pdb_file}")
                return {}

            # Group by chain (segid)
            chain_ids = set(protein.segments.segids)

            for chain_id in sorted(chain_ids):
                # Select this chain
                chain_atoms = u.select_atoms(f"protein and segid {chain_id}")

                if len(chain_atoms) == 0:
                    continue

                # Extract residues (unique by resid)
                residues = chain_atoms.residues

                # Build sequence
                sequence_parts = []
                residue_ids = []
                residue_names = []

                for residue in residues:
                    resname = residue.resname
                    resid = residue.resid

                    # Convert to single-letter code
                    if resname in self.three_to_one:
                        sequence_parts.append(self.three_to_one[resname])
                        residue_ids.append(resid)
                        residue_names.append(resname)
                    else:
                        # Unknown residue - use X
                        logger.debug(f"Unknown residue {resname} at position {resid} in chain {chain_id}")
                        sequence_parts.append('X')
                        residue_ids.append(resid)
                        residue_names.append(resname)

                sequence = ''.join(sequence_parts)

                chains_data[chain_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'residue_ids': residue_ids,
                    'residue_names': residue_names
                }

                logger.info(f"Chain {chain_id}: {len(sequence)} residues")

            return chains_data

        except Exception as e:
            logger.error(f"Failed to extract sequences from {pdb_file}: {e}")
            raise

    def save_sequences_to_csv(
        self,
        chains_data: Dict[str, Dict],
        output_csv: str,
        task_name: Optional[str] = None
    ) -> None:
        """
        Save extracted sequences to CSV file.

        CSV format:
        TaskName,ChainID,Length,Sequence,ResidueIDs,ResidueNames

        Args:
            chains_data: Dictionary from extract_sequences_from_pdb()
            output_csv: Path to output CSV file
            task_name: Optional task/sample name for the CSV

        Example:
            >>> extractor.save_sequences_to_csv(
            ...     chains_data,
            ...     "sequences.csv",
            ...     task_name="1ao7_run1"
            ... )
        """
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Prepare data for CSV
        rows = []
        for chain_id, chain_info in sorted(chains_data.items()):
            # Convert lists to comma-separated strings
            residue_ids_str = ','.join(map(str, chain_info['residue_ids']))
            residue_names_str = ','.join(chain_info['residue_names'])

            rows.append({
                'TaskName': task_name or '',
                'ChainID': chain_id,
                'Length': chain_info['length'],
                'Sequence': chain_info['sequence'],
                'ResidueIDs': residue_ids_str,
                'ResidueNames': residue_names_str
            })

        # Write CSV
        fieldnames = ['TaskName', 'ChainID', 'Length', 'Sequence', 'ResidueIDs', 'ResidueNames']

        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        logger.info(f"Saved {len(rows)} chain sequences to {output_csv}")

    def extract_and_save(
        self,
        pdb_file: str,
        output_csv: str,
        task_name: Optional[str] = None
    ) -> Dict[str, Dict]:
        """
        Convenience method: extract sequences and save to CSV in one step.

        Args:
            pdb_file: Path to PDB file
            output_csv: Path to output CSV file
            task_name: Optional task name

        Returns:
            Dictionary of extracted chain data

        Example:
            >>> extractor = PDBSequenceExtractor()
            >>> chains = extractor.extract_and_save(
            ...     "1ao7.pdb",
            ...     "1ao7_sequences.csv",
            ...     "1ao7"
            ... )
        """
        chains_data = self.extract_sequences_from_pdb(pdb_file)

        if chains_data:
            if task_name is None:
                task_name = Path(pdb_file).stem
            self.save_sequences_to_csv(chains_data, output_csv, task_name)
        else:
            logger.warning(f"No sequences extracted from {pdb_file}")

        return chains_data

    def batch_extract(
        self,
        pdb_files: List[str],
        output_csv: str
    ) -> None:
        """
        Extract sequences from multiple PDB files to single CSV.

        Args:
            pdb_files: List of PDB file paths
            output_csv: Path to output CSV file

        Example:
            >>> extractor = PDBSequenceExtractor()
            >>> extractor.batch_extract(
            ...     ["1ao7.pdb", "2vlk.pdb"],
            ...     "all_sequences.csv"
            ... )
        """
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        all_rows = []

        for pdb_file in pdb_files:
            pdb_path = Path(pdb_file)
            task_name = pdb_path.stem

            try:
                chains_data = self.extract_sequences_from_pdb(pdb_file)

                for chain_id, chain_info in sorted(chains_data.items()):
                    residue_ids_str = ','.join(map(str, chain_info['residue_ids']))
                    residue_names_str = ','.join(chain_info['residue_names'])

                    all_rows.append({
                        'TaskName': task_name,
                        'ChainID': chain_id,
                        'Length': chain_info['length'],
                        'Sequence': chain_info['sequence'],
                        'ResidueIDs': residue_ids_str,
                        'ResidueNames': residue_names_str
                    })

            except Exception as e:
                logger.error(f"Failed to process {pdb_file}: {e}")
                continue

        # Write all to CSV
        fieldnames = ['TaskName', 'ChainID', 'Length', 'Sequence', 'ResidueIDs', 'ResidueNames']

        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)

        logger.info(f"Batch extracted {len(all_rows)} chains from {len(pdb_files)} PDB files")


def extract_sequences_from_pdb(pdb_file: str, output_csv: str) -> Dict[str, Dict]:
    """
    Convenience function to extract sequences from PDB to CSV.

    Args:
        pdb_file: Path to PDB file
        output_csv: Path to output CSV file

    Returns:
        Dictionary of extracted chain sequences

    Example:
        >>> from immunoscope.utils import extract_sequences_from_pdb
        >>> chains = extract_sequences_from_pdb("1ao7.pdb", "sequences.csv")
    """
    extractor = PDBSequenceExtractor()
    return extractor.extract_and_save(pdb_file, output_csv)
