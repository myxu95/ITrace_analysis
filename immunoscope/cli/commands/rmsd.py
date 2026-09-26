#!/usr/bin/env python3
"""
IMS RMSD - RMSD Calculation and Analysis

Usage:
    ims rmsd -f trajectory.xtc -s topology.tpr --selection SELECTION [-o output.xvg]

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import argparse
import sys
import logging
from pathlib import Path

from immunoscope.analysis.trajectory import RMSDCalculator


def create_parser():
    """Create argument parser for rmsd subcommand"""
    parser = argparse.ArgumentParser(
        prog='ims rmsd',
        description='RMSD calculation and analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Backbone RMSD:
    ims rmsd -f processed.xtc -s md.tpr --selection backbone

  C-alpha RMSD:
    ims rmsd -f processed.xtc -s md.tpr --selection "protein and name CA"

  Save to file:
    ims rmsd -f processed.xtc -s md.tpr --selection backbone -o rmsd.xvg

  Use GROMACS backend:
    ims rmsd -f processed.xtc -s md.tpr --selection backbone --backend gromacs

Common selections:
  backbone          Protein backbone (N, CA, C)
  protein and name CA    C-alpha atoms only
  resid 1-50        Specific residue range
  chainID A         Specific chain
        """
    )

    # Required arguments
    parser.add_argument(
        '-f', '--trajectory',
        required=True,
        metavar='FILE',
        help='Input trajectory file'
    )

    parser.add_argument(
        '-s', '--topology',
        required=True,
        metavar='FILE',
        help='Structure/topology file (e.g., md.tpr) - GROMACS standard'
    )

    parser.add_argument(
        '--selection',
        required=True,
        metavar='STR',
        help='Atom selection (MDAnalysis syntax)'
    )

    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        metavar='FILE',
        help='Output file (default: rmsd.xvg)'
    )

    parser.add_argument(
        '--backend',
        choices=['mdanalysis', 'gromacs'],
        default='mdanalysis',
        help='Calculation backend (default: mdanalysis)'
    )

    parser.add_argument(
        '--ref-frame',
        type=int,
        default=0,
        metavar='N',
        help='Reference frame index (default: 0)'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    return parser


def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(levelname)s: %(message)s'
    )


def main(argv=None):
    """Main function for rmsd subcommand"""
    parser = create_parser()
    args = parser.parse_args(argv)

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Validate input files
    if not Path(args.trajectory).exists():
        logger.error(f"Trajectory file not found: {args.trajectory}")
        return 1

    if not Path(args.topology).exists():
        logger.error(f"Topology file not found: {args.topology}")
        return 1

    # Set default output file
    output_file = args.output if args.output else 'rmsd.xvg'

    # Print processing information
    logger.info("=" * 60)
    logger.info("IMS RMSD - RMSD Calculation")
    logger.info("=" * 60)
    logger.info(f"Trajectory: {args.trajectory}")
    logger.info(f"Topology: {args.topology}")
    logger.info(f"Selection: {args.selection}")
    logger.info(f"Backend: {args.backend}")
    logger.info(f"Output: {output_file}")
    logger.info("")

    # Calculate RMSD
    try:
        logger.info("Calculating RMSD...")

        calc = RMSDCalculator(args.topology, args.trajectory)

        if args.backend == 'mdanalysis':
            times, rmsd_values = calc.calculate_mdanalysis(
                selection=args.selection,
                ref_frame=args.ref_frame
            )
        else:  # gromacs
            times, rmsd_values = calc.calculate_gromacs(
                selection=args.selection,
                ref_frame=args.ref_frame
            )

        # Save results
        calc.save_xvg(times, rmsd_values, output_file,
                     xlabel='Time (ps)', ylabel='RMSD (nm)',
                     title=f'RMSD: {args.selection}')

        # Calculate statistics
        import numpy as np
        mean_rmsd = np.mean(rmsd_values)
        std_rmsd = np.std(rmsd_values)
        min_rmsd = np.min(rmsd_values)
        max_rmsd = np.max(rmsd_values)

        logger.info("")
        logger.info("=" * 60)
        logger.info("RMSD Statistics")
        logger.info("=" * 60)
        logger.info(f"Mean RMSD: {mean_rmsd:.3f} nm")
        logger.info(f"Std RMSD:  {std_rmsd:.3f} nm")
        logger.info(f"Min RMSD:  {min_rmsd:.3f} nm")
        logger.info(f"Max RMSD:  {max_rmsd:.3f} nm")
        logger.info(f"N frames:  {len(rmsd_values)}")
        logger.info("")
        logger.info(f"Results saved to: {output_file}")

        return 0

    except Exception as e:
        logger.error("")
        logger.error(f"RMSD calculation failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
