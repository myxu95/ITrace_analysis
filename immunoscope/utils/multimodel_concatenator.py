"""
Multi-Model PDB Concatenator

Merge multi-model PDB files into single PDB with continuous atom numbering.
No alignment or coordinate changes - just concatenate and renumber atoms.

Author: ImmunoScope Development Team
Date: 2026-01
"""

from pathlib import Path
from typing import List, Dict, Tuple
import logging


class MultiModelConcatenator:
    """Merge multi-model PDB into single PDB with continuous atom numbering.

    This class handles PDB files containing multiple MODEL records (common in
    biological assemblies) and concatenates them into a single PDB file with:
    - Continuous atom numbering (1, 2, 3, ...)
    - All coordinates preserved unchanged
    - MODEL/ENDMDL records removed
    - Original metadata (HEADER, CRYST1) preserved

    Usage:
        concatenator = MultiModelConcatenator()
        result = concatenator.concatenate_models("input.pdb", "output.pdb")
    """

    def __init__(self):
        """Initialize concatenator."""
        self.logger = logging.getLogger(__name__)

    def count_models(self, pdb_file: str) -> int:
        """Count number of MODEL records in PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Number of MODEL records found (0 if no MODEL records)
        """
        count = 0
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('MODEL'):
                        count += 1
        except Exception as e:
            self.logger.error(f"Failed to read {pdb_file}: {e}")
            raise

        return count

    def concatenate_models(
        self,
        input_pdb: str,
        output_pdb: str,
        preserve_remarks: bool = True
    ) -> Dict:
        """Concatenate multi-model PDB into single PDB with continuous atom numbering.

        Args:
            input_pdb: Input multi-model PDB file path
            output_pdb: Output single PDB file path
            preserve_remarks: Keep HEADER/TITLE/CRYST1/REMARK lines

        Returns:
            Dictionary with processing results:
            {
                'input_file': str,
                'output_file': str,
                'n_models': int,
                'total_atoms': int,
                'atoms_per_model': List[int],
                'status': str
            }
        """
        # Check input file
        if not Path(input_pdb).exists():
            raise FileNotFoundError(f"Input file not found: {input_pdb}")

        # Count models
        n_models = self.count_models(input_pdb)
        self.logger.info(f"Processing {input_pdb}: {n_models} model(s) found")

        # Handle single-model case
        if n_models <= 1:
            self.logger.info("Single model or no MODEL records - copying file")
            import shutil
            shutil.copy(input_pdb, output_pdb)
            return {
                'input_file': input_pdb,
                'output_file': output_pdb,
                'n_models': n_models,
                'total_atoms': self._count_atoms(input_pdb),
                'atoms_per_model': [self._count_atoms(input_pdb)],
                'status': 'single_model_copied'
            }

        # Parse multi-model PDB
        header_lines, atom_lines, atoms_per_model = self._parse_multimodel_pdb(
            input_pdb, preserve_remarks
        )

        # Renumber atoms continuously
        renumbered_lines = self._renumber_atoms(atom_lines)

        # Write output
        self._write_pdb(output_pdb, header_lines, renumbered_lines)

        total_atoms = len(atom_lines)
        self.logger.info(f"Merged {n_models} models -> {total_atoms} atoms")
        self.logger.info(f"  Atoms per model: {atoms_per_model}")

        return {
            'input_file': input_pdb,
            'output_file': output_pdb,
            'n_models': n_models,
            'total_atoms': total_atoms,
            'atoms_per_model': atoms_per_model,
            'status': 'success'
        }

    def _parse_multimodel_pdb(
        self,
        pdb_file: str,
        preserve_remarks: bool
    ) -> Tuple[List[str], List[str], List[int]]:
        """Parse multi-model PDB file.

        Returns:
            (header_lines, atom_lines, atoms_per_model)
        """
        header_lines = []
        atom_lines = []
        atoms_per_model = []
        current_model_atoms = 0
        in_model = False
        before_first_model = True

        with open(pdb_file, 'r') as f:
            for line in f:
                # Skip MODEL record
                if line.startswith('MODEL'):
                    if current_model_atoms > 0:
                        atoms_per_model.append(current_model_atoms)
                    current_model_atoms = 0
                    in_model = True
                    before_first_model = False
                    continue

                # Skip ENDMDL record
                if line.startswith('ENDMDL'):
                    in_model = False
                    continue

                # Collect header lines (before first MODEL)
                if before_first_model and preserve_remarks:
                    if line.startswith(('HEADER', 'TITLE', 'CRYST1', 'REMARK')):
                        header_lines.append(line)
                    continue

                # Skip lines after first model but before entering any model
                if not in_model and not before_first_model:
                    continue

                # Collect ATOM/HETATM lines
                if line.startswith(('ATOM', 'HETATM')):
                    atom_lines.append(line)
                    current_model_atoms += 1

        # Record last model's atom count
        if current_model_atoms > 0:
            atoms_per_model.append(current_model_atoms)

        return header_lines, atom_lines, atoms_per_model

    def _renumber_atoms(self, atom_lines: List[str]) -> List[str]:
        """Renumber atoms continuously starting from 1.

        Args:
            atom_lines: List of ATOM/HETATM record strings

        Returns:
            List of renumbered ATOM/HETATM records
        """
        renumbered = []
        atom_number = 1

        for line in atom_lines:
            # PDB format specifications (columns are 1-indexed in spec, 0-indexed in Python)
            # Columns 1-6: Record type (ATOM or HETATM)
            # Columns 7-11: Atom serial number (right-justified)
            # Column 12: Space
            # Columns 13-16: Atom name
            # Column 17: Alternate location indicator
            # Columns 18-20: Residue name
            # Column 21: Space
            # Column 22: Chain identifier
            # Columns 23-26: Residue sequence number
            # Column 27: Insertion code
            # Columns 28-30: Space
            # Columns 31-38: X coordinate
            # Columns 39-46: Y coordinate
            # Columns 47-54: Z coordinate
            # Columns 55-60: Occupancy
            # Columns 61-66: Temperature factor
            # Columns 77-78: Element symbol
            # Columns 79-80: Charge

            # Extract fields (using 0-indexed slicing)
            record_type = line[0:6]       # "ATOM  " or "HETATM"
            # Skip old atom number [6:11]
            atom_name = line[12:16]       # Atom name
            alt_loc = line[16]            # Alternate location
            res_name = line[17:20]        # Residue name
            chain_id = line[21]           # Chain ID
            res_seq = line[22:26]         # Residue number
            icode = line[26]              # Insertion code
            coords = line[27:54]          # X, Y, Z coordinates (27 chars)
            occupancy = line[54:60]       # Occupancy (6 chars)
            temp_factor = line[60:66]     # Temperature factor (6 chars)

            # Handle optional element/charge (may be missing)
            remainder = line[66:].rstrip('\n')

            # Reconstruct line with new atom number
            new_line = (
                f"{record_type}"                    # ATOM or HETATM (6 chars)
                f"{atom_number:5d}"                 # Atom serial (5 chars, right-aligned)
                f" {atom_name}"                     # Space + atom name (5 chars total)
                f"{alt_loc}{res_name} {chain_id}"  # Alt loc + res name + space + chain (6 chars)
                f"{res_seq}{icode}"                 # Res seq + icode (5 chars)
                f"{coords}"                         # Coordinates (27 chars)
                f"{occupancy}{temp_factor}"         # Occupancy + B-factor (12 chars)
                f"{remainder}\n"                    # Element, charge, etc.
            )

            renumbered.append(new_line)
            atom_number += 1

            # Check for overflow (PDB format limit: 99999)
            if atom_number > 99999:
                self.logger.warning(
                    f"Atom number exceeds 99999 - PDB format may be invalid"
                )

        return renumbered

    def _write_pdb(
        self,
        output_file: str,
        header_lines: List[str],
        atom_lines: List[str]
    ):
        """Write PDB file with header and atom records.

        Args:
            output_file: Output PDB file path
            header_lines: Header records (HEADER, CRYST1, etc.)
            atom_lines: ATOM/HETATM records
        """
        # Create output directory if needed
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w') as f:
            # Write header
            f.writelines(header_lines)

            # Write atoms
            f.writelines(atom_lines)

            # Write END record
            f.write("END\n")

    def _count_atoms(self, pdb_file: str) -> int:
        """Count ATOM/HETATM records in PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Number of ATOM/HETATM records
        """
        count = 0
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        count += 1
        except Exception as e:
            self.logger.error(f"Failed to count atoms in {pdb_file}: {e}")
            raise

        return count

    def batch_concatenate(
        self,
        input_dir: str,
        output_dir: str,
        pattern: str = "*.pdb",
        preserve_remarks: bool = True
    ) -> List[Dict]:
        """Batch process multi-model PDB files in directory.

        Args:
            input_dir: Input directory containing PDB files
            output_dir: Output directory for concatenated PDBs
            pattern: Glob pattern for matching files (default: "*.pdb")
            preserve_remarks: Keep HEADER/TITLE/CRYST1/REMARK lines

        Returns:
            List of processing results (one dict per file)
        """
        input_path = Path(input_dir)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Find PDB files
        pdb_files = sorted(input_path.glob(pattern))
        self.logger.info(f"Found {len(pdb_files)} PDB files in {input_dir}")

        results = []
        for pdb_file in pdb_files:
            output_file = output_path / f"{pdb_file.stem}_concatenated.pdb"

            self.logger.info(f"Processing {pdb_file.name}...")
            try:
                result = self.concatenate_models(
                    str(pdb_file),
                    str(output_file),
                    preserve_remarks=preserve_remarks
                )
                results.append(result)
                print(f"  ✓ {result['n_models']} models -> {result['total_atoms']} atoms")
            except Exception as e:
                self.logger.error(f"  ✗ Failed: {e}")
                results.append({
                    'input_file': str(pdb_file),
                    'status': 'failed',
                    'error': str(e)
                })

        return results
