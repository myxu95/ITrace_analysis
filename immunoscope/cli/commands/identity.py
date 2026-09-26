#!/usr/bin/env python3
"""IMS Identity - basic biological identity annotation."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.pipeline import BiologicalIdentityPipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="ims identity",
        description="Biological identity annotation for pHLA-TCR complexes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims identity --structure md_processed_converted.pdb -o ./output/identity_demo

  ims identity --structure md_processed_converted.pdb --cdr-metadata cdr_metadata.json -o ./output/identity_demo
        """,
    )
    parser.add_argument("--structure", required=True, metavar="FILE", help="Structure PDB used for annotation")
    parser.add_argument("--cdr-metadata", default=None, metavar="FILE", help="Optional existing cdr_metadata.json")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
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

    cdr_metadata_path = Path(args.cdr_metadata) if args.cdr_metadata else None
    if cdr_metadata_path and not cdr_metadata_path.exists():
        logger.error(f"CDR metadata not found: {cdr_metadata_path}")
        return 1

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    system_id = structure_path.stem

    context = PipelineContext(
        system_id=system_id,
        topology=str(structure_path.resolve()),
        trajectory_raw=str(structure_path.resolve()),
        structure_pdb=str(structure_path.resolve()),
        output_dir=str(output_dir),
    )

    pipeline = BiologicalIdentityPipeline(
        cdr_metadata_path=str(cdr_metadata_path.resolve()) if cdr_metadata_path else None,
        auto_identify_chains=True,
        auto_detect_cdr=cdr_metadata_path is None,
    )
    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("Biological identity annotation failed")
        for error in result.errors:
            logger.error(error)
        return 1

    identity_result = result.results.get("biological_identity", {})
    logger.info("Biological identity annotation completed")
    logger.info(f"JSON: {identity_result.get('json_file')}")
    logger.info(
        json.dumps(
            {
                "hla_identity": identity_result.get("hla_identity", {}),
                "peptide_identity": identity_result.get("peptide_identity", {}),
                "tcr_identity": identity_result.get("tcr_identity", {}),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
