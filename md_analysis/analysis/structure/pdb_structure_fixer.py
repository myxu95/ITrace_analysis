#!/usr/bin/env python3
"""
PDB Structure Fixer using PDBFixer

This module provides functionality to repair PDB structures:
- Add missing residues
- Add missing atoms
- Add missing hydrogens
- Remove heterogens (optional)
- Fix non-standard residues

Author: ITrace Development Team
Date: 2026-01-22
"""

import logging
from pathlib import Path
from typing import Optional, List, Dict
import tempfile
import shutil

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False

logger = logging.getLogger(__name__)


class PDBStructureFixer:
    """
    PDB structure fixer using PDBFixer.

    Repairs common issues in PDB structures:
    - Missing residues
    - Missing heavy atoms
    - Missing hydrogens
    - Non-standard residues
    """

    def __init__(self,
                 add_missing_residues: bool = True,
                 add_missing_atoms: bool = True,
                 add_hydrogens: bool = True,
                 remove_heterogens: bool = False,
                 keep_water: bool = False,
                 skip_terminal_residues: bool = True,
                 pH: float = 7.0):
        """
        Initialize PDB structure fixer.

        Args:
            add_missing_residues: Add missing residues based on SEQRES records
            add_missing_atoms: Add missing heavy atoms
            add_hydrogens: Add missing hydrogen atoms
            remove_heterogens: Remove all heterogens (ligands, ions, etc.)
            keep_water: Keep water molecules (only if remove_heterogens=True)
            skip_terminal_residues: Skip missing residues at chain termini (default: True)
                                   Terminal residues are often disordered and adding them
                                   may introduce unreliable structure. Only internal
                                   missing residues (loops) will be added.
            pH: pH for protonation state (default: 7.0)
        """
        if not PDBFIXER_AVAILABLE:
            raise ImportError(
                "PDBFixer is not available. Install with: conda install -c conda-forge pdbfixer"
            )

        self.add_missing_residues = add_missing_residues
        self.add_missing_atoms = add_missing_atoms
        self.add_hydrogens = add_hydrogens
        self.remove_heterogens = remove_heterogens
        self.keep_water = keep_water
        self.skip_terminal_residues = skip_terminal_residues
        self.pH = pH

    def _filter_terminal_missing_residues(self, fixer) -> int:
        """
        Filter out terminal missing residues from fixer.missingResidues.

        Only keeps internal missing residues (gaps/loops), removes N- and C-terminal
        missing residues which are often disordered.

        Args:
            fixer: PDBFixer object with missingResidues already identified

        Returns:
            Number of terminal residues filtered out
        """
        if not fixer.missingResidues:
            return 0

        # Get existing residue ranges for each chain
        # Note: missingResidues keys use (chain_index, res_id), not (chain_id, res_id)
        chain_ranges = {}  # {chain_index: (min_resid, max_resid, chain_id)}

        for chain in fixer.topology.chains():
            chain_idx = chain.index
            chain_id = chain.id
            residue_ids = [res.id for res in chain.residues()]

            if residue_ids:
                # Extract numeric part of residue ID
                numeric_ids = []
                for res_id in residue_ids:
                    try:
                        # Handle both numeric and string residue IDs
                        numeric_ids.append(int(res_id.split()[0]))
                    except (ValueError, AttributeError):
                        pass

                if numeric_ids:
                    chain_ranges[chain_idx] = (min(numeric_ids), max(numeric_ids), chain_id)

        # Filter missing residues
        original_missing = dict(fixer.missingResidues)
        filtered_missing = {}
        n_terminal_filtered = 0

        for (chain_idx, res_id), res_names in original_missing.items():
            if chain_idx not in chain_ranges:
                # Unknown chain, keep it
                filtered_missing[(chain_idx, res_id)] = res_names
                continue

            min_res, max_res, chain_id = chain_ranges[chain_idx]

            # Check if this is a terminal residue
            is_n_terminal = (res_id < min_res)
            is_c_terminal = (res_id > max_res)

            if is_n_terminal or is_c_terminal:
                # Skip terminal residues
                n_terminal_filtered += len(res_names)
                term_type = "N-terminal" if is_n_terminal else "C-terminal"
                logger.info(f"  Skipping {term_type}: Chain {chain_id} Res {res_id} ({len(res_names)} residues)")
            else:
                # Keep internal missing residues
                filtered_missing[(chain_idx, res_id)] = res_names

        # Update fixer's missing residues
        fixer.missingResidues = filtered_missing

        return n_terminal_filtered

    def fix_structure(self,
                     input_pdb: str,
                     output_pdb: str,
                     verbose: bool = True) -> Dict:
        """
        Fix a single PDB structure.

        Args:
            input_pdb: Input PDB file path
            output_pdb: Output PDB file path
            verbose: Print detailed information

        Returns:
            Dictionary with fixing statistics
        """
        input_path = Path(input_pdb)
        output_path = Path(output_pdb)

        if not input_path.exists():
            raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

        if verbose:
            logger.info(f"Fixing structure: {input_path.name}")

        # Initialize fixer
        fixer = PDBFixer(filename=str(input_path))

        stats = {
            'input_file': str(input_path),
            'output_file': str(output_path),
            'missing_residues': 0,
            'missing_heavy_atoms': 0,
            'missing_hydrogens': 0,
            'removed_chains': [],
            'fixed': True
        }

        try:
            # Find missing residues
            if self.add_missing_residues:
                fixer.findMissingResidues()
                n_missing_residues_total = sum(len(residues) for residues in fixer.missingResidues.values())

                # Filter out terminal residues if requested
                n_terminal_filtered = 0
                if self.skip_terminal_residues and n_missing_residues_total > 0:
                    n_terminal_filtered = self._filter_terminal_missing_residues(fixer)

                n_missing_residues = sum(len(residues) for residues in fixer.missingResidues.values())
                stats['missing_residues'] = n_missing_residues
                stats['terminal_residues_skipped'] = n_terminal_filtered

                if verbose and n_missing_residues_total > 0:
                    logger.info(f"  Found {n_missing_residues_total} missing residues")
                    if n_terminal_filtered > 0:
                        logger.info(f"  Skipped {n_terminal_filtered} terminal residues")
                        logger.info(f"  Will add {n_missing_residues} internal missing residues")
                    else:
                        logger.info(f"  Will add {n_missing_residues} missing residues")

                    for chain_id, residues in fixer.missingResidues.items():
                        if residues:
                            logger.info(f"    Chain {chain_id}: {len(residues)} residues")

            # Find non-standard residues
            fixer.findNonstandardResidues()
            if verbose and fixer.nonstandardResidues:
                logger.info(f"  Found {len(fixer.nonstandardResidues)} non-standard residues")

            # Replace non-standard residues with standard ones
            fixer.replaceNonstandardResidues()

            # Remove heterogens if requested
            if self.remove_heterogens:
                fixer.removeHeterogens(keepWater=self.keep_water)
                if verbose:
                    logger.info(f"  Removed heterogens (keep_water={self.keep_water})")

            # Find missing atoms
            if self.add_missing_atoms:
                fixer.findMissingAtoms()
                n_missing_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
                stats['missing_heavy_atoms'] = n_missing_atoms

                if verbose and n_missing_atoms > 0:
                    logger.info(f"  Found {n_missing_atoms} missing heavy atoms")

                # Add missing atoms
                fixer.addMissingAtoms()

            # Add hydrogens
            if self.add_hydrogens:
                fixer.addMissingHydrogens(pH=self.pH)
                if verbose:
                    logger.info(f"  Added missing hydrogens (pH={self.pH})")
                stats['added_hydrogens'] = True

            # Write output
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

            if verbose:
                logger.info(f"  ✓ Fixed structure written to: {output_path}")

        except Exception as e:
            logger.error(f"  ✗ Error fixing structure: {e}")
            stats['fixed'] = False
            stats['error'] = str(e)
            raise

        return stats

    def fix_structure_safe(self,
                          input_pdb: str,
                          output_pdb: str,
                          backup: bool = True,
                          verbose: bool = True) -> Dict:
        """
        Fix structure with backup and error recovery.

        If input and output are the same file, creates a backup first.

        Args:
            input_pdb: Input PDB file path
            output_pdb: Output PDB file path
            backup: Create backup of original file
            verbose: Print detailed information

        Returns:
            Dictionary with fixing statistics
        """
        input_path = Path(input_pdb).resolve()
        output_path = Path(output_pdb).resolve()

        is_same_file = (input_path == output_path)

        if is_same_file and backup:
            backup_path = input_path.with_suffix('.pdb.backup')
            shutil.copy2(input_path, backup_path)
            if verbose:
                logger.info(f"Created backup: {backup_path}")

        try:
            # Use temporary file if writing to same location
            if is_same_file:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb',
                                                 dir=output_path.parent,
                                                 delete=False) as tmp:
                    temp_output = Path(tmp.name)

                stats = self.fix_structure(str(input_path), str(temp_output), verbose=verbose)

                # Verify output is not empty
                if temp_output.stat().st_size == 0:
                    raise ValueError("Fixed PDB file is empty")

                # Move temp file to final location
                shutil.move(str(temp_output), str(output_path))
            else:
                stats = self.fix_structure(str(input_path), str(output_path), verbose=verbose)

            return stats

        except Exception as e:
            logger.error(f"Error during structure fixing: {e}")
            # Restore from backup if same file
            if is_same_file and backup:
                backup_path = input_path.with_suffix('.pdb.backup')
                if backup_path.exists():
                    shutil.copy2(backup_path, input_path)
                    logger.info(f"Restored from backup: {backup_path}")
            raise

    def batch_fix(self,
                  input_files: List[str],
                  output_dir: str,
                  suffix: str = "_fixed",
                  verbose: bool = True) -> List[Dict]:
        """
        Fix multiple PDB structures.

        Args:
            input_files: List of input PDB file paths
            output_dir: Output directory
            suffix: Suffix to add to output filenames (default: "_fixed")
            verbose: Print detailed information

        Returns:
            List of dictionaries with fixing statistics for each file
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        results = []

        for i, input_file in enumerate(input_files, 1):
            input_path = Path(input_file)
            output_file = output_path / f"{input_path.stem}{suffix}.pdb"

            if verbose:
                logger.info(f"\n[{i}/{len(input_files)}] Processing: {input_path.name}")

            try:
                stats = self.fix_structure(str(input_file), str(output_file), verbose=verbose)
                stats['status'] = 'success'
                results.append(stats)
            except Exception as e:
                logger.error(f"Failed to fix {input_path.name}: {e}")
                results.append({
                    'input_file': str(input_file),
                    'status': 'failed',
                    'error': str(e)
                })

        return results


def main():
    """Command-line interface for PDB structure fixer."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Fix PDB structures using PDBFixer"
    )
    parser.add_argument(
        "input",
        help="Input PDB file or directory"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output PDB file or directory (default: input_fixed.pdb or ./fixed/)"
    )
    parser.add_argument(
        "--no-missing-residues",
        action="store_true",
        help="Do not add missing residues"
    )
    parser.add_argument(
        "--no-missing-atoms",
        action="store_true",
        help="Do not add missing heavy atoms"
    )
    parser.add_argument(
        "--no-hydrogens",
        action="store_true",
        help="Do not add hydrogen atoms"
    )
    parser.add_argument(
        "--add-terminal-residues",
        action="store_true",
        help="Add terminal missing residues (default: skip terminals, only add internal gaps)"
    )
    parser.add_argument(
        "--remove-heterogens",
        action="store_true",
        help="Remove all heterogens (ligands, ions, etc.)"
    )
    parser.add_argument(
        "--keep-water",
        action="store_true",
        help="Keep water molecules (only with --remove-heterogens)"
    )
    parser.add_argument(
        "--pH",
        type=float,
        default=7.0,
        help="pH for protonation state (default: 7.0)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format='%(levelname)s: %(message)s'
    )

    # Initialize fixer
    fixer = PDBStructureFixer(
        add_missing_residues=not args.no_missing_residues,
        add_missing_atoms=not args.no_missing_atoms,
        add_hydrogens=not args.no_hydrogens,
        remove_heterogens=args.remove_heterogens,
        keep_water=args.keep_water,
        skip_terminal_residues=not args.add_terminal_residues,
        pH=args.pH
    )

    input_path = Path(args.input)

    # Single file or batch?
    if input_path.is_file():
        # Single file
        if args.output:
            output_file = args.output
        else:
            output_file = input_path.with_stem(f"{input_path.stem}_fixed")

        print("=" * 80)
        print("PDB Structure Fixer")
        print("=" * 80)
        print(f"Input:  {input_path}")
        print(f"Output: {output_file}")
        print("=" * 80)

        stats = fixer.fix_structure(str(input_path), str(output_file), verbose=True)

        print("\n" + "=" * 80)
        print("Summary:")
        print("=" * 80)
        print(f"Missing residues added:     {stats['missing_residues']}")
        print(f"Terminal residues skipped:  {stats.get('terminal_residues_skipped', 0)}")
        print(f"Missing heavy atoms added:  {stats['missing_heavy_atoms']}")
        print(f"Hydrogens added:            {stats.get('added_hydrogens', False)}")
        print("=" * 80)

    elif input_path.is_dir():
        # Batch processing
        pdb_files = sorted(input_path.glob("*.pdb"))

        if not pdb_files:
            print(f"No PDB files found in {input_path}")
            return 1

        output_dir = args.output if args.output else input_path / "fixed"

        print("=" * 80)
        print("Batch PDB Structure Fixer")
        print("=" * 80)
        print(f"Input directory:  {input_path}")
        print(f"Output directory: {output_dir}")
        print(f"Files to process: {len(pdb_files)}")
        print("=" * 80)

        results = fixer.batch_fix(
            [str(f) for f in pdb_files],
            str(output_dir),
            verbose=True
        )

        # Summary
        success = sum(1 for r in results if r.get('status') == 'success')
        failed = len(results) - success

        print("\n" + "=" * 80)
        print("Batch Summary:")
        print("=" * 80)
        print(f"Total files:  {len(results)}")
        print(f"Successful:   {success}")
        print(f"Failed:       {failed}")
        print("=" * 80)
    else:
        print(f"Error: {input_path} is not a file or directory")
        return 1

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
