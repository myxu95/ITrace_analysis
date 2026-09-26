#!/usr/bin/env python3
"""IMS RRCS - key interface RRCS analysis."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import RRCSPipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="ims rrcs",
        description="Residue-residue contact score analysis for pHLA-TCR interfaces",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims rrcs --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/rrcs_demo

  ims rrcs --structure md_processed_converted.pdb --trajectory md_processed.xtc -o ./output/rrcs_demo --stride 5
        """,
    )
    parser.add_argument("--structure", required=True, metavar="FILE", help="Representative structure PDB")
    parser.add_argument("--topology", metavar="FILE", help="Topology file for trajectory loading; if omitted, structure is used")
    parser.add_argument("--trajectory", required=True, metavar="FILE", help="Trajectory file")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
    parser.add_argument("--stride", type=int, default=1, metavar="N", help="Frame stride for RRCS sampling (default: 1)")
    parser.add_argument("--radius-min", type=float, default=3.23, metavar="A", help="RRCS minimum full-score distance in Angstrom (default: 3.23)")
    parser.add_argument("--radius-max", type=float, default=4.63, metavar="A", help="RRCS zero-score distance in Angstrom (default: 4.63)")
    parser.add_argument(
        "--pair-scope",
        choices=["interface", "cdr3_peptide", "cdr3_groove", "tcr_peptide", "tcr_groove", "tcr_interface"],
        default="interface",
        help="Default residue-pair scope when no pair file is provided (default: interface)",
    )
    parser.add_argument(
        "--pair-file",
        metavar="FILE",
        help="Optional custom residue-pair file. Format per line: 'D:95-100 $ C:4-6'",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser


def setup_logging(verbose: bool = False):
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(levelname)s: %(message)s")


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    structure_path = Path(args.structure)
    if not structure_path.exists():
        logger.error(f"Structure file not found: {structure_path}")
        return 1

    topology_path = Path(args.topology) if args.topology else structure_path
    if not topology_path.exists():
        logger.error(f"Topology file not found: {topology_path}")
        return 1

    trajectory_path = Path(args.trajectory)
    if not trajectory_path.exists():
        logger.error(f"Trajectory file not found: {trajectory_path}")
        return 1
    if args.pair_file and not Path(args.pair_file).exists():
        logger.error(f"Pair file not found: {args.pair_file}")
        return 1

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    system_id = structure_path.stem

    context = PipelineContext(
        system_id=system_id,
        topology=str(topology_path.resolve()),
        trajectory_raw=str(trajectory_path.resolve()),
        structure_pdb=str(structure_path.resolve()),
        output_dir=str(output_dir),
    )

    pipeline = RRCSPipeline(
        radius_min=args.radius_min,
        radius_max=args.radius_max,
        stride=args.stride,
        pair_scope=args.pair_scope,
        pair_file=str(Path(args.pair_file).resolve()) if args.pair_file else None,
        auto_identify_chains=True,
        auto_detect_cdr=True,
    )
    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("RRCS analysis failed")
        for error in result.errors:
            logger.error(error)
        return 1

    rrcs_result = result.results.get("rrcs", {})
    logger.info("RRCS analysis completed")
    logger.info(f"RRCS time series: {rrcs_result.get('timeseries_file')}")
    logger.info(f"RRCS pair summary: {rrcs_result.get('pair_summary_file')}")
    logger.info(f"RRCS region summary: {rrcs_result.get('region_summary_file')}")
    logger.info(f"Pair scope: {rrcs_result.get('pair_scope')}")
    if rrcs_result.get("pair_file"):
        logger.info(f"Pair file: {rrcs_result.get('pair_file')}")

    summary_path = rrcs_result.get("summary_file")
    if summary_path:
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        logger.info(
            json.dumps(
                {
                    "n_pairs": summary.get("n_pairs"),
                    "n_nonzero_pairs": summary.get("n_nonzero_pairs"),
                    "top_pairs": summary.get("top_pairs", [])[:5],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
