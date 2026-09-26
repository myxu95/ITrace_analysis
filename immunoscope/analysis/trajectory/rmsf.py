#!/usr/bin/env python3
"""
RMSF (Root Mean Square Fluctuation) Analysis Module

Provides comprehensive RMSF analysis capabilities for MD trajectories,
with specialized support for pHLA-TCR complexes including CDR3 regions.

Author: ImmunoScope Development Team
Date: 2025-12-09
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
from typing import Dict, List, Tuple, Optional
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class RMSFAnalyzer:
    """
    RMSF analysis for MD trajectories.

    Provides three levels of analysis:
    - Level 1: Global RMSF (entire system)
    - Level 2: Per-chain RMSF (individual chains)
    - Level 3: CDR3 RMSF (specific regions like CDR3 loops)
    """

    def __init__(self, universe: mda.Universe, output_dir: Optional[Path] = None):
        """
        Initialize RMSF analyzer.

        Args:
            universe: MDAnalysis Universe object
            output_dir: Output directory for results
        """
        self.u = universe
        self.output_dir = Path(output_dir) if output_dir else Path("./rmsf_results")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Results storage
        self.results = {}

    def calculate_global_rmsf(self,
                             selection: str = "protein and name CA",
                             save: bool = True) -> Dict:
        """
        Calculate global RMSF for entire system.

        Args:
            selection: Atom selection string
            save: Whether to save results to CSV

        Returns:
            Dictionary with RMSF statistics and data
        """
        logger.info("Calculating global RMSF...")

        # Select atoms
        atoms = self.u.select_atoms(selection)
        n_atoms = len(atoms)

        if n_atoms == 0:
            logger.error(f"No atoms found for selection: {selection}")
            return {}

        logger.info(f"  Selected {n_atoms} atoms for RMSF calculation")

        # Calculate RMSF
        rmsf_analyzer = RMSF(atoms).run()
        rmsf_values = rmsf_analyzer.results.rmsf

        # Create DataFrame
        df = pd.DataFrame({
            'residue_id': atoms.resids,
            'residue_name': atoms.resnames,
            'chain_id': atoms.chainIDs,
            'residue_seq_index': range(len(atoms)),
            'rmsf_angstrom': rmsf_values
        })

        # Save if requested
        if save:
            output_file = self.output_dir / 'global_rmsf.csv'
            df.to_csv(output_file, index=False)
            logger.info(f"  Global RMSF saved to {output_file}")

        # Calculate statistics
        stats = {
            'mean': float(np.mean(rmsf_values)),
            'std': float(np.std(rmsf_values)),
            'min': float(np.min(rmsf_values)),
            'max': float(np.max(rmsf_values)),
            'median': float(np.median(rmsf_values)),
            'n_residues': n_atoms,
            'dataframe': df
        }

        logger.info(f"  Mean RMSF: {stats['mean']:.2f} Å")
        logger.info(f"  Max RMSF: {stats['max']:.2f} Å")

        self.results['global_rmsf'] = stats

        return stats

    def calculate_chain_rmsf(self,
                            chain_selections: Dict[str, str],
                            save: bool = True) -> Dict:
        """
        Calculate RMSF for individual chains.

        Args:
            chain_selections: Dictionary mapping chain names to selection strings
            save: Whether to save results to CSV

        Returns:
            Dictionary with per-chain RMSF statistics
        """
        logger.info("Calculating per-chain RMSF...")

        chain_results = {}

        for chain_name, selection in chain_selections.items():
            logger.info(f"  Processing {chain_name}...")

            # Select CA atoms
            ca_selection = f"{selection} and name CA"
            atoms = self.u.select_atoms(ca_selection)

            if len(atoms) == 0:
                logger.warning(f"    No atoms found for {chain_name}")
                continue

            # Calculate RMSF
            rmsf_analyzer = RMSF(atoms).run()
            rmsf_values = rmsf_analyzer.results.rmsf

            logger.info(f"    {len(atoms)} residues, Mean RMSF: {np.mean(rmsf_values):.2f} Å")

            # Create DataFrame
            df = pd.DataFrame({
                'residue_id': atoms.resids,
                'residue_name': atoms.resnames,
                'residue_seq_index': range(len(atoms)),
                'rmsf_angstrom': rmsf_values
            })

            # Save individual chain
            if save:
                output_file = self.output_dir / f'{chain_name}_rmsf.csv'
                df.to_csv(output_file, index=False)

            # Store statistics
            chain_results[chain_name] = {
                'mean': float(np.mean(rmsf_values)),
                'std': float(np.std(rmsf_values)),
                'min': float(np.min(rmsf_values)),
                'max': float(np.max(rmsf_values)),
                'median': float(np.median(rmsf_values)),
                'n_residues': len(atoms),
                'dataframe': df
            }

        # Save combined results
        if save and chain_results:
            combined_df = pd.DataFrame()
            for chain_name, data in chain_results.items():
                temp_df = data['dataframe'].copy()
                temp_df['chain'] = chain_name
                combined_df = pd.concat([combined_df, temp_df], ignore_index=True)

            output_file = self.output_dir / 'chain_rmsf.csv'
            combined_df.to_csv(output_file, index=False)
            logger.info(f"  Chain RMSF saved to {output_file}")

        self.results['chain_rmsf'] = chain_results

        return chain_results

    def calculate_region_rmsf(self,
                             region_name: str,
                             selection: str,
                             save: bool = True) -> Dict:
        """
        Calculate RMSF for a specific region.

        Args:
            region_name: Name of the region (e.g., "CDR3_alpha")
            selection: Atom selection string
            save: Whether to save results to CSV

        Returns:
            Dictionary with RMSF statistics
        """
        logger.info(f"Calculating RMSF for {region_name}...")

        # Select atoms
        atoms = self.u.select_atoms(selection)

        if len(atoms) == 0:
            logger.warning(f"  No atoms found for {region_name}")
            return {}

        # Calculate RMSF
        rmsf_analyzer = RMSF(atoms).run()
        rmsf_values = rmsf_analyzer.results.rmsf

        # Create DataFrame
        df = pd.DataFrame({
            'residue_id': atoms.resids,
            'residue_name': atoms.resnames,
            'residue_seq_index': range(len(atoms)),
            'rmsf_angstrom': rmsf_values
        })

        # Save if requested
        if save:
            output_file = self.output_dir / f'{region_name}_rmsf.csv'
            df.to_csv(output_file, index=False)
            logger.info(f"  {region_name} RMSF saved to {output_file}")

        # Calculate statistics
        stats = {
            'mean': float(np.mean(rmsf_values)),
            'std': float(np.std(rmsf_values)),
            'min': float(np.min(rmsf_values)),
            'max': float(np.max(rmsf_values)),
            'median': float(np.median(rmsf_values)),
            'n_residues': len(atoms),
            'residue_range': (int(atoms.resids[0]), int(atoms.resids[-1])),
            'dataframe': df
        }

        logger.info(f"  Mean RMSF: {stats['mean']:.2f} ± {stats['std']:.2f} Å")

        return stats


def extract_sequence_from_topology(universe: mda.Universe,
                                   selection: str) -> str:
    """
    Extract amino acid sequence from MD topology.

    Args:
        universe: MDAnalysis Universe object
        selection: Atom selection string

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
    atoms = universe.select_atoms(f"{selection} and name CA")

    # Convert to single-letter sequence
    sequence = ''.join([aa_map.get(res, 'X') for res in atoms.resnames])

    return sequence


def find_subsequence_position(full_seq: str, subseq: str) -> Optional[Tuple[int, int]]:
    """
    Find subsequence position in full sequence.

    Args:
        full_seq: Full amino acid sequence
        subseq: Subsequence to find

    Returns:
        Tuple of (start_index, end_index) or None if not found
    """
    pos = full_seq.find(subseq)

    if pos == -1:
        return None

    return (pos, pos + len(subseq))


def calculate_cdr3_rmsf(universe: mda.Universe,
                       chain_selection: str,
                       cdr3_sequence: str,
                       chain_name: str = "TCR",
                       output_dir: Optional[Path] = None) -> Optional[Dict]:
    """
    Calculate RMSF for a CDR3 region.

    Args:
        universe: MDAnalysis Universe object
        chain_selection: Selection string for the TCR chain
        cdr3_sequence: CDR3 amino acid sequence
        chain_name: Name of the chain (for labeling)
        output_dir: Output directory

    Returns:
        Dictionary with CDR3 RMSF statistics or None if sequence not found
    """
    logger.info(f"Calculating CDR3 RMSF for {chain_name}...")
    logger.info(f"  CDR3 sequence: {cdr3_sequence}")

    # Extract full chain sequence
    full_seq = extract_sequence_from_topology(universe, chain_selection)
    logger.info(f"  Full chain length: {len(full_seq)} residues")

    # Find CDR3 position
    pos = find_subsequence_position(full_seq, cdr3_sequence)

    if pos is None:
        logger.warning(f"  CDR3 sequence not found in {chain_name}")
        return None

    start_idx, end_idx = pos
    logger.info(f"  CDR3 position: sequence index {start_idx}-{end_idx}")

    # Select CDR3 CA atoms
    ca_atoms = universe.select_atoms(f"{chain_selection} and name CA")
    cdr3_atoms = ca_atoms[start_idx:end_idx]

    logger.info(f"  Residue range: {cdr3_atoms.resids[0]}-{cdr3_atoms.resids[-1]}")

    # Calculate RMSF
    rmsf_analyzer = RMSF(cdr3_atoms).run()
    rmsf_values = rmsf_analyzer.results.rmsf

    # Create DataFrame
    df = pd.DataFrame({
        'residue_id': cdr3_atoms.resids,
        'residue_name': cdr3_atoms.resnames,
        'residue_seq_index': range(len(cdr3_atoms)),
        'rmsf_angstrom': rmsf_values
    })

    # Save if output directory provided
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / f'{chain_name}_cdr3_rmsf.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"  CDR3 RMSF saved to {output_file}")

    # Calculate statistics
    stats = {
        'mean': float(np.mean(rmsf_values)),
        'std': float(np.std(rmsf_values)),
        'min': float(np.min(rmsf_values)),
        'max': float(np.max(rmsf_values)),
        'median': float(np.median(rmsf_values)),
        'sequence': cdr3_sequence,
        'residue_range': (int(cdr3_atoms.resids[0]), int(cdr3_atoms.resids[-1])),
        'n_residues': len(cdr3_atoms),
        'dataframe': df
    }

    logger.info(f"  Mean RMSF: {stats['mean']:.2f} ± {stats['std']:.2f} Å")

    return stats


# Convenience function for pHLA-TCR complexes
def analyze_phla_tcr_rmsf(universe: mda.Universe,
                         output_dir: Path,
                         cdr3_sequences: Optional[Dict[str, str]] = None) -> Dict:
    """
    Complete RMSF analysis for pHLA-TCR complex.

    Args:
        universe: MDAnalysis Universe object
        output_dir: Output directory
        cdr3_sequences: Optional dict with 'cdr3_alpha' and 'cdr3_beta' keys

    Returns:
        Dictionary with all RMSF results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # Initialize analyzer
    analyzer = RMSFAnalyzer(universe, output_dir)

    # Level 1: Global RMSF
    logger.info("\n=== Level 1: Global RMSF ===")
    global_stats = analyzer.calculate_global_rmsf()
    results['global_rmsf'] = {k: v for k, v in global_stats.items() if k != 'dataframe'}

    # Level 2: Per-chain RMSF
    logger.info("\n=== Level 2: Per-Chain RMSF ===")
    chain_selections = {
        'HLA_alpha': 'chainID A and protein',
        'HLA_beta': 'chainID B and protein',
        'Peptide': 'chainID C and protein',
        'TCR_alpha': 'chainID D and protein',
        'TCR_beta': 'chainID E and protein'
    }

    chain_stats = analyzer.calculate_chain_rmsf(chain_selections)
    results['chain_rmsf'] = {
        chain: {k: v for k, v in stats.items() if k != 'dataframe'}
        for chain, stats in chain_stats.items()
    }

    # Level 3: CDR3 RMSF (optional)
    if cdr3_sequences:
        logger.info("\n=== Level 3: CDR3 RMSF ===")
        results['cdr3_rmsf'] = {}

        if 'cdr3_alpha' in cdr3_sequences:
            alpha_stats = calculate_cdr3_rmsf(
                universe,
                'chainID D and protein',
                cdr3_sequences['cdr3_alpha'],
                'CDR3_alpha',
                output_dir
            )
            if alpha_stats:
                results['cdr3_rmsf']['CDR3_alpha'] = {
                    k: v for k, v in alpha_stats.items() if k != 'dataframe'
                }

        if 'cdr3_beta' in cdr3_sequences:
            beta_stats = calculate_cdr3_rmsf(
                universe,
                'chainID E and protein',
                cdr3_sequences['cdr3_beta'],
                'CDR3_beta',
                output_dir
            )
            if beta_stats:
                results['cdr3_rmsf']['CDR3_beta'] = {
                    k: v for k, v in beta_stats.items() if k != 'dataframe'
                }

    return results
