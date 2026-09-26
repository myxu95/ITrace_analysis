#!/usr/bin/env python3
"""IMS RMSF - legacy selection mode and region-aware flexibility analysis."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis.rms import RMSF

from ...core.context import PipelineContext
from ...pipeline import AnnotatedRMSFPipeline


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ims rmsf",
        description="RMSF calculation and region-aware flexibility analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Legacy low-level selection mode
  ims rmsf -f processed.xtc -s md.tpr --selection "protein and name CA" -o rmsf.xvg

  # Region-aware flexibility analysis
  ims rmsf -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb \
    --annotated -d output/rmsf_demo_1OGA --stride 10
        """,
    )
    parser.add_argument("-f", "--trajectory", required=True, help="Trajectory file")
    parser.add_argument("-s", "--topology", required=True, help="Topology file")
    parser.add_argument("--structure", default=None, help="Structure PDB for semantic annotation")
    parser.add_argument("--selection", default=None, help="Atom selection for legacy mode")
    parser.add_argument("-o", "--output", default="rmsf.xvg", help="Legacy mode output file")
    parser.add_argument("-d", "--output-dir", default=None, help="Annotated mode output directory")
    parser.add_argument("--annotated", action="store_true", help="Run region-aware RMSF pipeline")
    parser.add_argument(
        "--split-atom-sets",
        action="store_true",
        help="Annotated mode only: also emit per-residue backbone/sidechain RMSF",
    )
    parser.add_argument("--stride", type=int, default=1, help="Frame stride")
    parser.add_argument("--time-unit", choices=["ps", "ns"], default="ps", help="Time unit for summaries")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def _write_legacy_xvg(trajectory: str, topology: str, selection: str, output: Path) -> dict:
    universe = mda.Universe(topology, trajectory)
    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        raise ValueError(f"No atoms found for selection: {selection}")
    values = RMSF(atoms).run().results.rmsf
    frame = pd.DataFrame({
        "resid": atoms.resids,
        "resname": atoms.resnames,
        "rmsf_angstrom": values,
    })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("# residue_index rmsf_angstrom\n")
        for idx, (_, row) in enumerate(frame.iterrows(), start=1):
            handle.write(f"{idx:6d} {float(row['rmsf_angstrom']):.6f}\n")
    return {
        "n_atoms": int(len(atoms)),
        "mean_rmsf_angstrom": float(frame["rmsf_angstrom"].mean()),
        "max_rmsf_angstrom": float(frame["rmsf_angstrom"].max()),
        "output_file": str(output),
    }


def _run_annotated_mode(args, logger: logging.Logger) -> int:
    if not args.structure:
        logger.error("--annotated mode requires --structure")
        return 1

    output_dir = args.output_dir or f"./output/{Path(args.trajectory).stem}_rmsf"
    context = PipelineContext(
        system_id=Path(output_dir).name,
        topology=args.topology,
        trajectory_raw=args.trajectory,
        structure_pdb=args.structure,
        trajectory_processed=args.trajectory,
        output_dir=output_dir,
    )
    pipeline = AnnotatedRMSFPipeline(
        stride=args.stride,
        time_unit=args.time_unit,
        selection=args.selection,
        split_atom_sets=args.split_atom_sets,
    )
    result = pipeline.execute(context)
    if result.should_stop:
        logger.error("Annotated RMSF pipeline stopped unexpectedly")
        for err in result.errors:
            logger.error(err)
        return 1

    rmsf_result = result.results.get("residue_rmsf", {})
    logger.info("Annotated RMSF analysis completed")
    logger.info(f"Residue CSV: {rmsf_result.get('residue_csv')}")
    logger.info(f"Region summary: {rmsf_result.get('region_summary_csv')}")
    logger.info(f"Summary JSON: {rmsf_result.get('summary_json')}")
    return 0


def main(argv=None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info("=" * 60)
    logger.info("IMS RMSF")
    logger.info("=" * 60)
    logger.info(f"Trajectory: {args.trajectory}")
    logger.info(f"Topology: {args.topology}")
    if args.structure:
        logger.info(f"Structure: {args.structure}")
    logger.info("")

    try:
        if args.annotated:
            return _run_annotated_mode(args, logger)

        if not args.selection:
            logger.error("Legacy mode requires --selection, or use --annotated")
            return 1

        result = _write_legacy_xvg(
            trajectory=args.trajectory,
            topology=args.topology,
            selection=args.selection,
            output=Path(args.output),
        )
        logger.info("Legacy RMSF calculation completed")
        logger.info(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except Exception as exc:
        logger.error(f"RMSF calculation failed: {exc}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
