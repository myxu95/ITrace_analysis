#!/usr/bin/env python3
"""IMS Landscape - energy landscape analysis."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import LandscapeAnalysisPipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="ims landscape",
        description="Energy landscape analysis for pHLA-TCR complexes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic landscape analysis:
    ims landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/landscape_demo

  Disable interactive HTML:
    ims landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/landscape_demo --no-interactive

  Use UMAP projection:
    ims landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc \
      -o ./output/landscape_umap --reducer umap --umap-n-neighbors 20 --umap-min-dist 0.05

  Use internal TICA slow-coordinate projection:
    ims landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc \
      -o ./output/landscape_tica --reducer tica --tica-lag 10

  Inject precomputed RMSD curves:
    ims landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc \
      --rmsd-file peptide_rmsd=analysis/rmsd/peptide_rmsd.xvg \
      --rmsd-file cdr3a_rmsd=analysis/rmsd/cdr3a_rmsd.csv \
      -o ./output/landscape_demo
        """,
    )
    parser.add_argument("--structure", required=True, metavar="FILE", help="Representative structure PDB")
    parser.add_argument("--topology", metavar="FILE", help="Topology file for trajectory loading; if omitted, structure is used")
    parser.add_argument("--trajectory", required=True, metavar="FILE", help="Trajectory file")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
    parser.add_argument("--stride", type=int, default=1, metavar="N", help="Trajectory stride for upstream analyses (default: 1)")
    parser.add_argument("--hotspot-threshold", type=float, default=0.5, metavar="F", help="Hotspot contact frequency threshold (default: 0.5)")
    parser.add_argument("--reducer", choices=["pca", "umap", "tica"], default="pca", help="Dimensionality reducer for landscape projection (default: pca)")
    parser.add_argument("--n-components", type=int, default=3, metavar="N", help="Number of retained reduced dimensions (default: 3)")
    parser.add_argument("--umap-n-neighbors", type=int, default=15, metavar="N", help="UMAP neighbor count when --reducer umap (default: 15)")
    parser.add_argument("--umap-min-dist", type=float, default=0.1, metavar="F", help="UMAP minimum distance when --reducer umap (default: 0.1)")
    parser.add_argument("--umap-metric", default="euclidean", metavar="NAME", help="UMAP metric when --reducer umap (default: euclidean)")
    parser.add_argument("--random-state", type=int, default=42, metavar="N", help="Random seed for stochastic reducers (default: 42)")
    parser.add_argument("--tica-lag", type=int, default=10, metavar="N", help="TICA lag in frames when --reducer tica (default: 10)")
    parser.add_argument("--tica-regularization", type=float, default=1e-6, metavar="F", help="TICA covariance regularization when --reducer tica (default: 1e-6)")
    parser.add_argument("--no-tica-kinetic-map", action="store_true", help="Disable kinetic-map scaling for internal TICA")
    parser.add_argument("--no-tica-stabilize-sign", action="store_true", help="Disable deterministic TIC sign orientation")
    parser.add_argument("--temperature", type=float, default=300.0, metavar="K", help="FEL temperature in Kelvin (default: 300)")
    parser.add_argument("--bins", default="auto", metavar="SPEC", help="FEL bins setting, e.g. auto / 20 / 20,20 (default: auto)")
    parser.add_argument("--no-contact", action="store_true", help="Skip contact/contact-annotation features")
    parser.add_argument("--no-angles", action="store_true", help="Skip docking-angle features")
    parser.add_argument("--no-com-distance", action="store_true", help="Skip TCR-MHC COM distance features")
    parser.add_argument("--no-bsa", action="store_true", help="Skip buried-surface-area features")
    parser.add_argument("--no-interactive", action="store_true", help="Do not generate plotly interactive HTML")
    parser.add_argument(
        "--rmsd-file",
        action="append",
        default=[],
        metavar="NAME=FILE",
        help="Optional precomputed RMSD file to inject as extra conformation feature; repeatable",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser


def setup_logging(verbose: bool = False):
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(levelname)s: %(message)s")


def _parse_bins(raw_value: str):
    text = str(raw_value).strip()
    if text.lower() == "auto":
        return "auto"
    if "," in text:
        parts = [part.strip() for part in text.split(",") if part.strip()]
        if len(parts) != 2:
            raise ValueError("Comma-formatted bins must contain two integers, for example 20,20")
        return (int(parts[0]), int(parts[1]))
    return int(text)


def _parse_rmsd_file_args(raw_items: list[str]) -> dict[str, str]:
    result = {}
    for item in raw_items:
        if "=" not in item:
            raise ValueError(f"Invalid --rmsd-file argument, expected NAME=FILE: {item}")
        feature_name, file_path = item.split("=", 1)
        feature_name = feature_name.strip()
        file_path = file_path.strip()
        if not feature_name:
            raise ValueError(f"Invalid RMSD feature name: {item}")
        if not file_path:
            raise ValueError(f"Invalid RMSD file path: {item}")
        result[feature_name] = file_path
    return result


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

    try:
        bins = _parse_bins(args.bins)
        rmsd_files = _parse_rmsd_file_args(args.rmsd_file)
    except ValueError as exc:
        logger.error(str(exc))
        return 1

    for feature_name, file_path in rmsd_files.items():
        if not Path(file_path).exists():
            logger.error(f"RMSD file not found for {feature_name}: {file_path}")
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

    rmsd_result_keys = {}
    for index, (feature_name, file_path) in enumerate(rmsd_files.items()):
        result_key = f"landscape_rmsd_{index}"
        rmsd_result_keys[feature_name] = result_key
        context.results[result_key] = {
            "output_file": str(Path(file_path).resolve()),
        }

    pipeline = LandscapeAnalysisPipeline(
        stride=args.stride,
        auto_identify_chains=True,
        auto_detect_cdr=not args.no_contact,
        include_contact=not args.no_contact,
        include_angles=not args.no_angles,
        include_com_distance=not args.no_com_distance,
        include_bsa=not args.no_bsa,
        rmsd_result_keys=rmsd_result_keys or None,
        hotspot_threshold=args.hotspot_threshold,
        reducer=args.reducer,
        n_components=args.n_components,
        umap_n_neighbors=args.umap_n_neighbors,
        umap_min_dist=args.umap_min_dist,
        umap_metric=args.umap_metric,
        random_state=args.random_state,
        tica_lag=args.tica_lag,
        tica_regularization=args.tica_regularization,
        tica_kinetic_map=not args.no_tica_kinetic_map,
        tica_stabilize_sign=not args.no_tica_stabilize_sign,
        temperature=args.temperature,
        bins=bins,
        include_interactive=not args.no_interactive,
    )

    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("Landscape analysis failed")
        for error in result.errors:
            logger.error(error)
        return 1

    landscape_result = result.results.get("landscape", {})
    logger.info("Landscape analysis completed")
    logger.info(f"Reducer: {landscape_result.get('reducer')}")
    logger.info(f"Feature matrix: {landscape_result.get('feature_matrix')}")
    logger.info(f"Reduced coordinates: {landscape_result.get('pca_coordinates')}")
    logger.info(f"Loadings: {landscape_result.get('loadings')}")
    logger.info(f"Summary: {landscape_result.get('summary')}")
    logger.info(f"2D landscape: {landscape_result.get('landscape_2d')}")
    logger.info(f"Trajectory overlay: {landscape_result.get('trajectory_overlay')}")
    if landscape_result.get("landscape_interactive"):
        logger.info(f"Interactive HTML: {landscape_result.get('landscape_interactive')}")

    summary_path = landscape_result.get("summary")
    if summary_path and Path(summary_path).exists():
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        logger.info(
            json.dumps(
                {
                    "n_frames": summary.get("n_frames"),
                    "n_features": summary.get("n_features"),
                    "reducer": summary.get("reducer"),
                    "coordinate_labels": summary.get("coordinate_labels"),
                    "explained_variance": summary.get("explained_variance", []),
                    "reducer_metadata": summary.get("reducer_metadata", {}),
                    "top_contributors_pc1": summary.get("top_contributors_pc1", [])[:5],
                    "top_contributors_pc2": summary.get("top_contributors_pc2", [])[:5],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
