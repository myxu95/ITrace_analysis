#!/usr/bin/env python3
"""
IMS Quality - Quality Assessment and Validation

Usage:
    ims quality -f trajectory.xtc -s topology.tpr [-o output_dir]

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import argparse
import sys
import logging
from pathlib import Path

from immunoscope.pipeline import QualityAssessmentPipeline


def create_parser():
    """Create argument parser for quality subcommand"""
    parser = argparse.ArgumentParser(
        prog='ims quality',
        description='Quality assessment and validation for MD trajectories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic quality assessment:
    ims quality -f processed.xtc -s md.tpr

  Save report to specific directory:
    ims quality -f processed.xtc -s md.tpr -o ./quality_reports

  Custom RMSD selection:
    ims quality -f processed.xtc -s md.tpr --selection "protein and name CA"

  Fast mode (skip validation):
    ims quality -f processed.xtc -s md.tpr --skip-validation

Quality Grades:
  Grade A (Excellent):   Mean RMSD < 0.3 nm, Std < 0.05 nm
  Grade B (Good):        Mean RMSD < 0.5 nm, Std < 0.1 nm
  Grade C (Acceptable):  Mean RMSD < 0.8 nm, Std < 0.15 nm
  Grade D (Poor):        Mean RMSD > 0.8 nm or Std > 0.15 nm
        """
    )

    # Required arguments
    parser.add_argument(
        '-f', '--trajectory',
        required=True,
        metavar='FILE',
        help='Input trajectory file (e.g., processed.xtc)'
    )

    parser.add_argument(
        '-s', '--topology',
        required=True,
        metavar='FILE',
        help='Structure/topology file (e.g., md.tpr) - GROMACS standard'
    )

    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        metavar='DIR',
        help='Output directory for reports (default: current directory)'
    )

    parser.add_argument(
        '--selection',
        default='backbone',
        metavar='STR',
        help='MDAnalysis selection for RMSD (default: backbone)'
    )

    parser.add_argument(
        '--skip-validation',
        action='store_true',
        help='Skip post-PBC validation (faster)'
    )

    parser.add_argument(
        '--stride',
        type=int,
        default=1,
        metavar='N',
        help='Validation stride (use >1 for faster checks, default: 1)'
    )

    parser.add_argument(
        '--format',
        choices=['markdown', 'json', 'both'],
        default='both',
        help='Report format (default: both)'
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
    """Main function for quality subcommand"""
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

    # Set output directory
    output_dir = args.output if args.output else './quality_assessment'
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Print processing information
    logger.info("=" * 60)
    logger.info("IMS Quality - Quality Assessment")
    logger.info("=" * 60)
    logger.info(f"Trajectory: {args.trajectory}")
    logger.info(f"Topology: {args.topology}")
    logger.info(f"RMSD selection: {args.selection}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("")

    # Initialize quality assessment pipeline
    try:
        pipeline = QualityAssessmentPipeline()
    except Exception as e:
        logger.error(f"Failed to initialize quality pipeline: {e}")
        return 1

    # Run quality assessment
    try:
        logger.info("Running quality assessment...")

        results = pipeline.run_comprehensive_assessment(
            trajectory=str(trajectory_path),
            topology=str(topology_path),
            rmsd_selection=args.selection,
            output_dir=output_dir,
            validation_stride=args.stride,
            enable_post_pbc_validation=not args.skip_validation
        )

        # Generate reports
        logger.info("Generating reports...")

        if args.format in ['markdown', 'both']:
            report_file = Path(output_dir) / 'quality_report.md'
            pipeline.generate_quality_report(
                results=results,
                output_file=str(report_file),
                format='markdown'
            )
            logger.info(f"Markdown report: {report_file}")

        if args.format in ['json', 'both']:
            import json
            json_file = Path(output_dir) / 'quality_report.json'
            with open(json_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            logger.info(f"JSON report: {json_file}")

        # Print summary
        logger.info("")
        logger.info("=" * 60)
        logger.info("Quality Assessment Results")
        logger.info("=" * 60)
        logger.info(f"Overall Grade: {results['overall_grade']}")
        logger.info(f"Qualification: {'PASS' if results['is_qualified'] else 'FAIL'}")

        if 'rmsd_metrics' in results:
            rmsd = results['rmsd_metrics']
            logger.info(f"Mean RMSD: {rmsd['mean_rmsd']:.3f} nm")
            logger.info(f"Std RMSD: {rmsd['std_rmsd']:.3f} nm")
            if 'is_converged' in rmsd:
                logger.info(f"Converged: {'Yes' if rmsd['is_converged'] else 'No'}")

        logger.info("")
        logger.info(f"Full report saved to: {output_dir}")

        return 0

    except Exception as e:
        logger.error("")
        logger.error("=" * 60)
        logger.error("Quality assessment failed!")
        logger.error("=" * 60)
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
