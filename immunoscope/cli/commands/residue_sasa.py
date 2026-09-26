#!/usr/bin/env python3
"""IMS residue-sasa - Per-residue bound/unbound SASA analysis."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import ResidueSasaPipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="ims residue-sasa",
        description=(
            "Per-residue solvent accessible surface area in bound (complex) and "
            "unbound (own group only) states. Exposes per-residue exposure metrics "
            "for mutation design."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Auto-identify pHLA and TCR:
    ims residue-sasa -f md_processed.xtc -s md_processed_converted.pdb \\
      --structure md_processed_converted.pdb -o ./output/residue_sasa

  Manual selections:
    ims residue-sasa -f md_processed.xtc -s md.tpr \\
      --sel-a "chainID A or chainID B or chainID C" --sel-b "chainID D or chainID E" \\
      -o ./output/residue_sasa
        """,
    )
    parser.add_argument("-f", "--trajectory", required=True, metavar="FILE", help="Input trajectory file")
    parser.add_argument("-s", "--topology", required=True, metavar="FILE", help="Topology file")
    parser.add_argument("--structure", metavar="FILE", help="Structure PDB for automatic chain identification")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
    parser.add_argument("--sel-a", metavar="SEL", help="Manual selection for molecule A (e.g. pHLA)")
    parser.add_argument("--sel-b", metavar="SEL", help="Manual selection for molecule B (e.g. TCR)")
    parser.add_argument("--probe-radius", type=float, default=1.4, metavar="A",
                        help="Probe radius in Angstrom (default: 1.4)")
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
    structure_path = Path(args.structure) if args.structure else None

    if not trajectory_path.exists():
        logger.error(f"Trajectory file not found: {trajectory_path}")
        return 1
    if not topology_path.exists():
        logger.error(f"Topology file not found: {topology_path}")
        return 1
    if structure_path and not structure_path.exists():
        logger.error(f"Structure file not found: {structure_path}")
        return 1
    if (args.sel_a is None) ^ (args.sel_b is None):
        logger.error("Manual mode requires both --sel-a and --sel-b")
        return 1
    if args.sel_a is None and structure_path is None:
        logger.error("Auto mode requires --structure so chain identification can run")
        return 1

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    system_id = trajectory_path.stem

    logger.info("=" * 60)
    logger.info("IMS residue-sasa - Per-Residue SASA Analysis")
    logger.info("=" * 60)
    logger.info(f"Trajectory: {trajectory_path}")
    logger.info(f"Topology: {topology_path}")
    logger.info(f"Output dir: {output_dir}")
    logger.info(f"Probe radius: {args.probe_radius} A")
    logger.info(f"Stride: {args.stride}")
    logger.info("")

    context = PipelineContext(
        system_id=system_id,
        topology=str(topology_path.resolve()),
        trajectory_raw=str(trajectory_path.resolve()),
        trajectory_processed=str(trajectory_path.resolve()),
        structure_pdb=str(structure_path.resolve()) if structure_path else None,
        output_dir=str(output_dir),
    )

    pipeline = ResidueSasaPipeline(
        probe_radius=args.probe_radius,
        stride=args.stride,
        selection_a=args.sel_a,
        selection_b=args.sel_b,
        auto_identify_chains=(args.sel_a is None and args.sel_b is None),
    )
    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("Residue SASA analysis failed")
        for error in result.errors:
            logger.error(error)
        return 1

    sasa_result = result.results.get("residue_sasa", {})
    logger.info("Residue SASA analysis completed")
    logger.info(f"Per-residue CSV: {sasa_result.get('residue_csv')}")
    logger.info(f"Summary: {sasa_result.get('summary_json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
