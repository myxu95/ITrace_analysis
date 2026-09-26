#!/usr/bin/env python3
"""
IMS PDB - PDB File Processing

Usage:
    ims pdb download <pdb_id> [pdb_id...]
    ims pdb fix -f input.pdb -o output.pdb
    ims pdb standardize -f input.pdb -o output.pdb

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import argparse
import sys
import logging
from pathlib import Path


def create_parser():
    """Create argument parser for pdb subcommand"""
    parser = argparse.ArgumentParser(
        prog='ims pdb',
        description='PDB file processing and management',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Download PDB files:
    ims pdb download 1ao7 2ckb 6vma

  Fix PDB structure:
    ims pdb fix -f input.pdb -o fixed.pdb

  Standardize chain IDs:
    ims pdb standardize -f input.pdb -o standardized.pdb

  Batch processing:
    ims pdb batch --pdb-list pdbs.txt --workers 4
        """
    )

    subparsers = parser.add_subparsers(dest='action', help='PDB action')

    # Download subcommand
    download_parser = subparsers.add_parser(
        'download',
        help='Download PDB files from RCSB'
    )
    download_parser.add_argument(
        'pdb_ids',
        nargs='+',
        help='PDB IDs to download'
    )
    download_parser.add_argument(
        '-o', '--output',
        default='./pdbs',
        help='Output directory (default: ./pdbs)'
    )

    # Fix subcommand
    fix_parser = subparsers.add_parser(
        'fix',
        help='Fix PDB structure defects'
    )
    fix_parser.add_argument(
        '-f', '--input',
        required=True,
        help='Input PDB file'
    )
    fix_parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output PDB file'
    )

    # Standardize subcommand
    std_parser = subparsers.add_parser(
        'standardize',
        help='Standardize PDB chain IDs'
    )
    std_parser.add_argument(
        '-f', '--input',
        required=True,
        help='Input PDB file'
    )
    std_parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output PDB file'
    )

    # Batch subcommand
    batch_parser = subparsers.add_parser(
        'batch',
        help='Batch process multiple PDBs'
    )
    batch_parser.add_argument(
        '--pdb-list',
        required=True,
        help='File containing PDB IDs (one per line)'
    )
    batch_parser.add_argument(
        '--workers',
        type=int,
        default=4,
        help='Number of parallel workers (default: 4)'
    )
    batch_parser.add_argument(
        '-o', '--output',
        default='./processed_pdbs',
        help='Output directory (default: ./processed_pdbs)'
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


def handle_download(args, logger):
    """Handle PDB download action"""
    from immunoscope.utils import PDBDownloader

    downloader = PDBDownloader()
    Path(args.output).mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading {len(args.pdb_ids)} PDB(s) to {args.output}")

    success = 0
    failed = 0

    for pdb_id in args.pdb_ids:
        try:
            result = downloader.download_pdb(pdb_id, output_dir=args.output)
            if result['status'] == 'success':
                logger.info(f"✓ {pdb_id.upper()}: {result['file_path']}")
                success += 1
            else:
                logger.error(f"✗ {pdb_id.upper()}: {result.get('error', 'Unknown error')}")
                failed += 1
        except Exception as e:
            logger.error(f"✗ {pdb_id.upper()}: {e}")
            failed += 1

    logger.info(f"\nCompleted: {success} succeeded, {failed} failed")
    return 0 if failed == 0 else 1


def handle_fix(args, logger):
    """Handle PDB fix action"""
    from immunoscope.analysis import PDBStructureFixer

    if not Path(args.input).exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    logger.info(f"Fixing PDB structure: {args.input}")

    fixer = PDBStructureFixer()
    stats = fixer.fix_structure(args.input, args.output)

    logger.info(f"✓ Fixed PDB saved to: {args.output}")
    logger.info(f"  Internal residues added: {stats['internal_added']}")
    logger.info(f"  Terminal residues skipped: {stats['terminal_skipped']}")

    return 0


def handle_standardize(args, logger):
    """Handle PDB standardize action"""
    from immunoscope.analysis import PDBChainStandardizer

    if not Path(args.input).exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    logger.info(f"Standardizing PDB chains: {args.input}")

    standardizer = PDBChainStandardizer()
    result = standardizer.process_single(args.input, args.output)

    if result.status == 'OK':
        logger.info(f"✓ Standardized PDB saved to: {args.output}")
        logger.info(f"  Chain mapping: {result.chain_mapping}")
    else:
        logger.error(f"✗ Standardization failed: {result.status}")
        return 1

    return 0


def handle_batch(args, logger):
    """Handle PDB batch processing"""
    from immunoscope.cli.batch_pdb import PDBProcessingPipeline

    # Read PDB IDs from file
    with open(args.pdb_list, 'r') as f:
        pdb_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    logger.info(f"Processing {len(pdb_ids)} PDBs from {args.pdb_list}")

    pipeline = PDBProcessingPipeline(base_dir=args.output)
    results = pipeline.process_batch(pdb_ids, workers=args.workers)

    logger.info(f"\nCompleted: {results['succeeded']} succeeded, {results['failed']} failed")

    return 0 if results['failed'] == 0 else 1


def main(argv=None):
    """Main function for pdb subcommand"""
    parser = create_parser()
    args = parser.parse_args(argv)

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Check if action is specified
    if not args.action:
        parser.print_help()
        return 1

    # Dispatch to appropriate handler
    if args.action == 'download':
        return handle_download(args, logger)
    elif args.action == 'fix':
        return handle_fix(args, logger)
    elif args.action == 'standardize':
        return handle_standardize(args, logger)
    elif args.action == 'batch':
        return handle_batch(args, logger)
    else:
        logger.error(f"Unknown action: {args.action}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
