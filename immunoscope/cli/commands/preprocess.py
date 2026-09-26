#!/usr/bin/env python3
"""
IMS Preprocess - PBC Correction and Trajectory Processing

Usage:
    ims preprocess -f trajectory.xtc -s topology.tpr -o processed.xtc
    ims preprocess -f trajectory.xtc -s topology.tpr -o processed.xtc -m 2step

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import argparse
import sys
import logging
import shutil
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import PreprocessQualityPipeline


def create_parser():
    """Create argument parser for preprocess subcommand"""
    parser = argparse.ArgumentParser(
        prog='ims preprocess',
        description='PBC correction and trajectory preprocessing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic PBC processing (interactive method selection):
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc

  Specify 2-step method:
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc -m 2step

  Specify 3-step method:
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc -m 3step

  With custom fit group:
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc -m 2step --fit "protein and name CA"

  With time downsampling:
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc -m 2step --dt 100

  Verbose output:
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc -v

PBC Methods:
  2step     nojump + fit rot+trans (recommended, faster)
  3step     center + whole + fit (standard, more thorough)

GROMACS-style Parameters:
  -f        Input trajectory file (GROMACS -f)
  -s        Structure/topology file (GROMACS -s)
  -o        Output file (GROMACS -o)
  -m        Method selection (similar to GROMACS -method)
        """
    )

    # Required arguments (GROMACS style)
    parser.add_argument(
        '-f', '--trajectory',
        required=True,
        metavar='FILE',
        help='Input trajectory file (e.g., md.xtc, md.trr)'
    )

    parser.add_argument(
        '-s', '--topology',
        required=True,
        metavar='FILE',
        help='Structure/topology file (e.g., md.tpr) - GROMACS standard'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='FILE',
        help='Output processed trajectory file'
    )

    # Method selection (GROMACS style -m)
    parser.add_argument(
        '-m', '--method',
        choices=['2step', '3step'],
        metavar='METHOD',
        help='PBC correction method: 2step (recommended) or 3step (standard). If not specified, interactive menu will be shown.'
    )

    # Optional processing arguments
    parser.add_argument(
        '--fit',
        metavar='GROUP',
        help='Fit group for alignment (default: Backbone for 2step, auto-detect for 3step)'
    )

    parser.add_argument(
        '--output-group',
        metavar='GROUP',
        help='Output group selection (default: System)'
    )

    parser.add_argument(
        '--dt',
        type=float,
        metavar='PS',
        help='Time step for output trajectory in ps (for downsampling)'
    )

    parser.add_argument(
        '--rmsd-selection',
        default='backbone',
        metavar='SELECTION',
        help='RMSD atom selection for preprocessing quality assessment (default: backbone)'
    )

    parser.add_argument(
        '--gmx',
        default='gmx',
        metavar='CMD',
        help='GROMACS executable (default: gmx)'
    )

    parser.add_argument(
        '--keep-temp',
        action='store_true',
        help='Keep temporary intermediate files'
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


def interactive_method_selection():
    """
    Interactive menu for PBC method selection

    Returns:
        Selected method ('2step' or '3step')
    """
    print()
    print("=" * 70)
    print("PBC Correction Method Selection")
    print("=" * 70)
    print()
    print("Available methods:")
    print()
    print("  1. 2-step method (RECOMMENDED)")
    print("     - nojump: Remove atom jumps across periodic boundaries")
    print("     - fit: Align trajectory (remove rotation and translation)")
    print("     - Faster and more stable for large complexes")
    print("     - Best for: TCR-pMHC, antibody-antigen complexes")
    print()
    print("  2. 3-step method (STANDARD)")
    print("     - center: Center system on shortest chain")
    print("     - whole: Make molecules whole")
    print("     - fit: Align trajectory")
    print("     - More thorough, traditional approach")
    print("     - Best for: Standard protein systems")
    print()
    print("-" * 70)

    while True:
        try:
            choice = input("Select method [1-2] (default: 1): ").strip()

            # Default to 2step if user just presses Enter
            if choice == '':
                choice = '1'

            if choice == '1':
                print()
                print("Selected: 2-step method (nojump + fit)")
                print()
                return '2step'
            elif choice == '2':
                print()
                print("Selected: 3-step method (center + whole + fit)")
                print()
                return '3step'
            else:
                print("Invalid choice. Please enter 1 or 2.")
        except (KeyboardInterrupt, EOFError):
            print()
            print("Operation cancelled by user")
            sys.exit(1)


def main(argv=None):
    """Main function for preprocess subcommand"""
    parser = create_parser()
    args = parser.parse_args(argv)

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Validate input files
    trajectory_path = Path(args.trajectory)
    topology_path = Path(args.topology)

    if not trajectory_path.exists():
        logger.error(f"Trajectory file not found: {args.trajectory}")
        return 1

    if not topology_path.exists():
        logger.error(f"Topology file not found: {args.topology}")
        return 1

    # Interactive method selection if not specified
    if args.method is None:
        args.method = interactive_method_selection()

    # Print processing information
    logger.info("=" * 60)
    logger.info("IMS Preprocess - PBC Correction + RMSD Quality Assessment")
    logger.info("=" * 60)
    logger.info(f"Input trajectory: {args.trajectory}")
    logger.info(f"Topology file: {args.topology}")
    logger.info(f"PBC method: {args.method}")
    logger.info(f"Output file: {args.output}")
    if args.fit:
        logger.info(f"Fit group: {args.fit}")
    if args.output_group:
        logger.info(f"Output group: {args.output_group}")
        if args.dt:
            logger.info(f"Time step: {args.dt} ps")
        logger.info(f"RMSD selection: {args.rmsd_selection}")
        logger.info("")

    try:
        logger.info("Starting PBC correction...")
        logger.info("")

        output_path = Path(args.output).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

        context = PipelineContext(
            system_id=trajectory_path.stem,
            topology=str(topology_path.resolve()),
            trajectory_raw=str(trajectory_path.resolve()),
            output_dir=str(output_path.parent),
        )

        pipeline = PreprocessQualityPipeline(
            method=args.method,
            dt=args.dt,
            gmx_executable=args.gmx,
            fit_group=args.fit,
            output_group=args.output_group,
            center_group=args.fit if args.method == '3step' else None,
            keep_temp_files=args.keep_temp,
            rmsd_selection=args.rmsd_selection,
        )

        result = pipeline.execute(context)

        if result.has_errors():
            raise RuntimeError('; '.join(result.errors))

        generated_processed = Path(result.trajectory_processed).resolve()
        if generated_processed != output_path:
            if output_path.exists():
                output_path.unlink()
            shutil.move(str(generated_processed), str(output_path))
            result.trajectory_processed = str(output_path)

        generated_pdb = generated_processed.with_name(f"{generated_processed.stem}_converted.pdb")
        renamed_pdb = output_path.with_name(f"{output_path.stem}_converted.pdb")
        if generated_pdb.exists() and generated_pdb.resolve() != renamed_pdb.resolve():
            if renamed_pdb.exists():
                renamed_pdb.unlink()
            shutil.move(str(generated_pdb), str(renamed_pdb))

        logger.info("")
        logger.info("=" * 60)
        logger.info("Success!")
        logger.info("=" * 60)
        logger.info(f"Processed trajectory: {output_path}")
        if renamed_pdb.exists():
            logger.info(f"Converted PDB: {renamed_pdb}")
        rmsd_result = result.get_result('rmsd')
        if rmsd_result and rmsd_result.get('output_file'):
            logger.info(f"RMSD file: {rmsd_result['output_file']}")
        quality_result = result.get_result('preprocess_quality')
        if quality_result:
            metrics = quality_result['metrics']
            logger.info(f"RMSD full variation: {metrics['full_variation_nm']:.3f} nm")
            logger.info(f"RMSD tail90 variation: {metrics['tail90_variation_nm']:.3f} nm")
            logger.info(f"Quality report: {quality_result['reports']['markdown']}")

        # Check output file size
        if output_path.exists():
            size_mb = output_path.stat().st_size / (1024 * 1024)
            logger.info(f"Output file size: {size_mb:.2f} MB")

            # Calculate compression ratio if possible
            input_size_mb = trajectory_path.stat().st_size / (1024 * 1024)
            if args.dt:
                logger.info(f"Input file size: {input_size_mb:.2f} MB")
                logger.info(f"Note: File size may differ due to downsampling (dt={args.dt} ps)")

        logger.info("")
        logger.info("Next steps:")
        logger.info(f"  1. Check RMSD report in {output_path.parent / 'quality'}")
        logger.info(f"  2. Run deeper quality checks if needed: ims quality {output_path.parent}")
        logger.info(f"  3. Visualization: Open {args.output} in VMD or PyMOL")

        return 0

    except Exception as e:
        logger.error("")
        logger.error("=" * 60)
        logger.error("Processing failed!")
        logger.error("=" * 60)
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
