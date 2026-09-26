#!/usr/bin/env python3
"""IMS dssp - Per-residue secondary structure (DSSP) analysis."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import DSSPPipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="ims dssp",
        description=(
            "Per-residue secondary structure aggregated over the trajectory. "
            "Reports the dominant DSSP code (H/E/L) plus per-frame propensities "
            "for mutation-design risk scoring."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims dssp -f md_processed.xtc -s md_processed_converted.pdb -o ./output/dssp_demo
  ims dssp -f md_processed.xtc -s md.tpr --stride 10 -o ./output/dssp_demo
        """,
    )
    parser.add_argument("-f", "--trajectory", required=True, metavar="FILE", help="Input trajectory file")
    parser.add_argument("-s", "--topology", required=True, metavar="FILE", help="Topology file")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
    parser.add_argument("--selection", default="protein", metavar="SEL",
                        help="Atom selection passed to DSSP (default: protein)")
    parser.add_argument("--stride", type=int, default=1, metavar="N",
                        help="Frame stride (default: 1; use 10+ for long trajectories)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser


def setup_logging(verbose: bool = False):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    trajectory_path = Path(args.trajectory)
    topology_path = Path(args.topology)
    if not trajectory_path.exists():
        logger.error(f"Trajectory file not found: {trajectory_path}")
        return 1
    if not topology_path.exists():
        logger.error(f"Topology file not found: {topology_path}")
        return 1

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    system_id = trajectory_path.stem

    logger.info("=" * 60)
    logger.info("IMS dssp - Per-Residue Secondary Structure")
    logger.info("=" * 60)
    logger.info(f"Trajectory: {trajectory_path}")
    logger.info(f"Topology: {topology_path}")
    logger.info(f"Output dir: {output_dir}")
    logger.info(f"Selection: {args.selection}")
    logger.info(f"Stride: {args.stride}")
    logger.info("")

    context = PipelineContext(
        system_id=system_id,
        topology=str(topology_path.resolve()),
        trajectory_raw=str(trajectory_path.resolve()),
        trajectory_processed=str(trajectory_path.resolve()),
        output_dir=str(output_dir),
    )

    pipeline = DSSPPipeline(selection=args.selection, stride=args.stride)
    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("DSSP analysis failed")
        for error in result.errors:
            logger.error(error)
        return 1

    dssp_result = result.results.get("dssp", {})
    logger.info("DSSP analysis completed")
    logger.info(f"Per-residue CSV: {dssp_result.get('residue_csv')}")
    logger.info(f"Summary: {dssp_result.get('summary_json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
