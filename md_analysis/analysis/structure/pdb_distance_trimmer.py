#!/usr/bin/env python3
"""
PDB Distance Trimmer

Trim PDB structures by removing residues beyond a distance threshold from peptide center.
This serves as a safety measure to remove potentially incorrect terminal residues.

Author: ITrace Development Team
Date: 2026-01-22
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array
    MDA_AVAILABLE = True
except ImportError:
    MDA_AVAILABLE = False

logger = logging.getLogger(__name__)


class PDBDistanceTrimmer:
    """
    Trim PDB structures by distance from peptide center.

    Removes residues that are beyond a specified distance threshold from
    the center of the peptide chain. This is useful as a safety measure
    to remove potentially incorrect terminal residues that may have been
    incorrectly added during structure repair.
    """

    def __init__(self, distance_cutoff_nm: float = 8.0):
        """
        Initialize PDB distance trimmer.

        Args:
            distance_cutoff_nm: Distance cutoff in nanometers (default: 8.0 nm)
        """
        if not MDA_AVAILABLE:
            raise ImportError(
                "MDAnalysis is required for PDB trimming. "
                "Install with: conda install -c conda-forge mdanalysis"
            )

        self.distance_cutoff_nm = distance_cutoff_nm
        self.distance_cutoff_angstrom = distance_cutoff_nm * 10.0

    def find_peptide_chain(self, universe) -> str:
        """
        Find peptide chain (shortest protein chain).

        Args:
            universe: MDAnalysis Universe object

        Returns:
            Segment ID of the peptide chain

        Raises:
            ValueError: If no protein chains found
        """
        chains = {}
        for seg in universe.segments:
            # Count protein residues only
            protein_residues = universe.select_atoms(
                f"segid {seg.segid} and protein"
            ).residues
            if len(protein_residues) > 0:
                chains[seg.segid] = len(protein_residues)

        # Find shortest chain
        if not chains:
            raise ValueError("No protein chains found!")

        peptide_segid = min(chains.items(), key=lambda x: x[1])[0]
        peptide_length = chains[peptide_segid]

        logger.info(f"  Peptide chain: {peptide_segid} ({peptide_length} residues)")
        return peptide_segid

    def get_peptide_center_residue(self, universe, peptide_segid):
        """
        Get the central residue of peptide chain.

        Args:
            universe: MDAnalysis Universe object
            peptide_segid: Segment ID of peptide chain

        Returns:
            MDAnalysis Residue object representing the center residue
        """
        peptide_residues = universe.select_atoms(
            f"segid {peptide_segid} and protein"
        ).residues

        n_residues = len(peptide_residues)
        center_idx = n_residues // 2
        center_residue = peptide_residues[center_idx]

        logger.info(
            f"  Center residue: {center_residue.resname}{center_residue.resid} "
            f"(position {center_idx+1}/{n_residues})"
        )
        return center_residue

    def trim_structure(self,
                      input_pdb: str,
                      output_pdb: str,
                      verbose: bool = True) -> Dict:
        """
        Trim PDB structure by distance from peptide center.

        Args:
            input_pdb: Input PDB file path
            output_pdb: Output PDB file path
            verbose: Print detailed information

        Returns:
            Dictionary with trimming statistics
        """
        input_path = Path(input_pdb)
        output_path = Path(output_pdb)

        if not input_path.exists():
            raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

        if verbose:
            logger.info(f"Trimming structure: {input_path.name}")
            logger.info(f"  Distance cutoff: {self.distance_cutoff_nm} nm")

        # Load structure
        u = mda.Universe(str(input_path))

        # Find peptide chain
        peptide_segid = self.find_peptide_chain(u)

        # Get peptide center residue
        center_residue = self.get_peptide_center_residue(u, peptide_segid)
        center_com = center_residue.atoms.center_of_mass()

        # Calculate distances for all residues
        all_residues = u.select_atoms("protein").residues

        residues_to_keep = []
        residues_to_remove = []

        for residue in all_residues:
            residue_com = residue.atoms.center_of_mass()
            distance = np.linalg.norm(residue_com - center_com)

            if distance <= self.distance_cutoff_angstrom:
                residues_to_keep.append(residue)
            else:
                residues_to_remove.append({
                    'segid': residue.segment.segid,
                    'resid': residue.resid,
                    'resname': residue.resname,
                    'distance_nm': distance / 10.0
                })

        # Print removed residues if verbose
        if verbose and residues_to_remove:
            logger.info(f"  Removing {len(residues_to_remove)} residues beyond {self.distance_cutoff_nm} nm:")
            # Group by chain
            by_chain = {}
            for res in residues_to_remove:
                segid = res['segid']
                if segid not in by_chain:
                    by_chain[segid] = []
                by_chain[segid].append(res)

            for segid in sorted(by_chain.keys()):
                chain_residues = by_chain[segid]
                logger.info(f"    Chain {segid}: {len(chain_residues)} residues")

        # Create atom selection for kept residues
        atom_indices = []
        for res in residues_to_keep:
            atom_indices.extend(res.atoms.indices)

        if not atom_indices:
            raise ValueError("No residues left after trimming!")

        # Select and write trimmed structure
        trimmed_atoms = u.atoms[atom_indices]

        output_path.parent.mkdir(parents=True, exist_ok=True)
        trimmed_atoms.write(str(output_path))

        if verbose:
            retention_rate = len(residues_to_keep) / len(all_residues) * 100
            logger.info(f"  Kept {len(residues_to_keep)}/{len(all_residues)} residues ({retention_rate:.1f}%)")
            logger.info(f"  ✓ Trimmed structure written to: {output_path}")

        # Return statistics
        return {
            'input_file': str(input_path),
            'output_file': str(output_path),
            'cutoff_nm': self.distance_cutoff_nm,
            'total_residues': len(all_residues),
            'kept_residues': len(residues_to_keep),
            'removed_residues': len(residues_to_remove),
            'retention_rate': len(residues_to_keep) / len(all_residues),
            'removed_details': residues_to_remove
        }

    def batch_trim(self,
                   input_files: List[str],
                   output_dir: str,
                   suffix: str = "_trimmed",
                   verbose: bool = True) -> List[Dict]:
        """
        Trim multiple PDB structures.

        Args:
            input_files: List of input PDB file paths
            output_dir: Output directory
            suffix: Suffix to add to output filenames (default: "_trimmed")
            verbose: Print detailed information

        Returns:
            List of dictionaries with trimming statistics for each file
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
                stats = self.trim_structure(str(input_file), str(output_file), verbose=verbose)
                stats['status'] = 'success'
                results.append(stats)
            except Exception as e:
                logger.error(f"Failed to trim {input_path.name}: {e}")
                results.append({
                    'input_file': str(input_file),
                    'status': 'failed',
                    'error': str(e)
                })

        return results


def main():
    """Command-line interface for PDB distance trimmer."""
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Trim PDB by removing residues beyond distance from peptide center"
    )
    parser.add_argument(
        "input",
        help="Input PDB file or directory"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output PDB file or directory (default: input_trimmed.pdb or ./trimmed/)"
    )
    parser.add_argument(
        "-d", "--distance",
        type=float,
        default=8.0,
        help="Distance cutoff in nanometers (default: 8.0 nm)"
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

    # Initialize trimmer
    trimmer = PDBDistanceTrimmer(distance_cutoff_nm=args.distance)

    input_path = Path(args.input)

    # Single file or batch?
    if input_path.is_file():
        # Single file
        if args.output:
            output_file = args.output
        else:
            output_file = input_path.with_stem(f"{input_path.stem}_trimmed")

        print("=" * 80)
        print("PDB Distance Trimmer")
        print("=" * 80)
        print(f"Input:  {input_path}")
        print(f"Output: {output_file}")
        print(f"Cutoff: {args.distance} nm")
        print("=" * 80)

        stats = trimmer.trim_structure(str(input_path), str(output_file), verbose=True)

        print("\n" + "=" * 80)
        print("Summary:")
        print("=" * 80)
        print(f"Total residues:      {stats['total_residues']}")
        print(f"Kept residues:       {stats['kept_residues']}")
        print(f"Removed residues:    {stats['removed_residues']}")
        print(f"Retention rate:      {stats['retention_rate']*100:.1f}%")
        print("=" * 80)

    elif input_path.is_dir():
        # Batch processing
        pdb_files = sorted(input_path.glob("*.pdb"))

        if not pdb_files:
            print(f"No PDB files found in {input_path}")
            return 1

        output_dir = args.output if args.output else input_path / "trimmed"

        print("=" * 80)
        print("Batch PDB Distance Trimmer")
        print("=" * 80)
        print(f"Input directory:  {input_path}")
        print(f"Output directory: {output_dir}")
        print(f"Files to process: {len(pdb_files)}")
        print(f"Cutoff: {args.distance} nm")
        print("=" * 80)

        results = trimmer.batch_trim(
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
