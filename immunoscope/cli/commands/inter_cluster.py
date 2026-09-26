#!/usr/bin/env python3
"""IMS InterCluster - key interface-state clustering analysis."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import InterfaceClusteringPipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="ims inter_cluster",
        description="Interface-aware clustering for pHLA-TCR key binding states",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims inter_cluster -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb -o ./output/interface_cluster

  ims inter_cluster -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb -o ./output/interface_cluster --stride 5 --distance-cutoff 0.30
        """,
    )
    parser.add_argument("-f", "--trajectory", required=True, metavar="FILE", help="Processed trajectory file")
    parser.add_argument("-s", "--topology", required=True, metavar="FILE", help="Topology file")
    parser.add_argument("--structure", required=True, metavar="FILE", help="Reference structure PDB")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
    parser.add_argument(
        "--stride",
        type=int,
        default=10,
        metavar="N",
        help="Geometry/sidechain frame stride; contact fingerprint uses ceil(stride/2) (default: 10)",
    )
    parser.add_argument("--contact-cutoff", type=float, default=4.5, metavar="A", help="Coarse contact cutoff in Angstrom (default: 4.5)")
    parser.add_argument("--distance-cutoff", type=float, default=0.35, metavar="X", help="Hierarchical clustering cutoff on joint distance (default: 0.35)")
    parser.add_argument("--linkage-method", default="average", metavar="NAME", help="Linkage method (default: average)")
    parser.add_argument("--geometry-weight", type=float, default=0.45, metavar="W", help="Backbone geometry weight (default: 0.45)")
    parser.add_argument("--sidechain-weight", type=float, default=0.25, metavar="W", help="Sidechain orientation weight (default: 0.25)")
    parser.add_argument("--interaction-weight", type=float, default=0.30, metavar="W", help="Interface interaction fingerprint weight (default: 0.30)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser


def setup_logging(verbose: bool = False):
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(levelname)s: %(message)s")


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    topology_path = Path(args.topology)
    trajectory_path = Path(args.trajectory)
    structure_path = Path(args.structure)
    for label, path in [("topology", topology_path), ("trajectory", trajectory_path), ("structure", structure_path)]:
        if not path.exists():
            logger.error(f"{label} file not found: {path}")
            return 1

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    system_id = structure_path.stem

    context = PipelineContext(
        system_id=system_id,
        topology=str(topology_path.resolve()),
        trajectory_raw=str(trajectory_path.resolve()),
        trajectory_processed=str(trajectory_path.resolve()),
        structure_pdb=str(structure_path.resolve()),
        output_dir=str(output_dir),
    )

    pipeline = InterfaceClusteringPipeline(
        stride=args.stride,
        contact_cutoff_angstrom=args.contact_cutoff,
        distance_cutoff=args.distance_cutoff,
        linkage_method=args.linkage_method,
        geometry_weight=args.geometry_weight,
        sidechain_weight=args.sidechain_weight,
        interaction_weight=args.interaction_weight,
        auto_identify_chains=True,
        auto_detect_cdr=True,
    )
    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("Interface clustering failed")
        for error in result.errors:
            logger.error(error)
        return 1

    cluster_result = result.results.get("interface_clustering", {})
    logger.info("Interface clustering completed")
    logger.info(f"Frame assignments: {cluster_result.get('frame_assignments_csv')}")
    logger.info(f"Summary table: {cluster_result.get('summary_table_csv')}")
    logger.info(f"Summary JSON: {cluster_result.get('summary_json')}")

    summary_path = cluster_result.get("summary_json")
    if summary_path:
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        logger.info(
            json.dumps(
                {
                    "n_frames": summary.get("n_frames"),
                    "n_clusters": summary.get("n_clusters"),
                    "largest_cluster": summary.get("largest_cluster"),
                },
                ensure_ascii=False,
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
