#!/usr/bin/env python3
"""
PDB Batch Processing Pipeline - CLI Module

Complete workflow from PDB ID to MD-ready structure:
1. Download: Fetch PDB files from RCSB
2. Concatenate: Merge multi-model PDBs (auto-detect)
3. Fix: Repair structure defects (add missing residues)
4. Standardize: Standardize chain IDs to A/B/C/D/E
5. Trim: Distance-based trimming (8nm default)

CLI Usage:
    immunoscope-batch-pdb 1ao7 2ckb 6vma
    immunoscope-batch-pdb --pdb-list my_pdbs.txt
    immunoscope-batch-pdb --pdb-list pdbs.txt --workers 4
"""

import argparse
import json
import shutil
import sys
import time
import multiprocessing
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import MDAnalysis as mda

from immunoscope.utils import (
    PDBDownloader,
    MultiModelConcatenator,
    PDBStructureFixer,
    PDBChainStandardizer,
    PDBDistanceTrimmer
)


class PDBProcessingPipeline:
    """
    PDB batch processing pipeline

    Complete workflow:
    1. Download: Download PDB from RCSB
    2. Concatenate: Merge multi-model PDBs (auto-detect)
    3. Fix: Repair structure defects
    4. Standardize: Standardize chain IDs
    5. Trim: Distance-based trimming
    """

    def __init__(
        self,
        base_dir: str = "input/standardizedpdbs",
        skip_existing: bool = True,
        verbose: bool = False
    ):
        """
        Initialize pipeline

        Args:
            base_dir: Base directory for outputs
            skip_existing: Skip existing files
            verbose: Verbose output
        """
        self.base_dir = Path(base_dir)
        self.skip_existing = skip_existing
        self.verbose = verbose

        # Initialize processing tools
        self.downloader = PDBDownloader()
        self.concatenator = MultiModelConcatenator()
        self.fixer = PDBStructureFixer()
        self.standardizer = PDBChainStandardizer()
        self.trimmer = PDBDistanceTrimmer()

        # Directory structure
        self.dirs = {
            'raw': self.base_dir / 'pdbs_raw',
            'concat': self.base_dir / 'pdbs_concatenated',
            'fixed': self.base_dir / 'fix_pdbs',
            'standardized': self.base_dir / 'protein_only',
            'trimmed': self.base_dir / 'trimmed_sd_pdbs'
        }

        # Create directories
        for dir_path in self.dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

    def process_batch(
        self,
        pdb_ids: List[str],
        steps: List[str] = ['download', 'fix', 'standardize', 'trim'],
        workers: int = 1
    ) -> Dict:
        """
        Batch process PDB list

        Args:
            pdb_ids: List of PDB IDs
            steps: Steps to execute
            workers: Number of parallel workers (1=serial, >1=parallel, 'auto'=CPU count)

        Returns:
            Processing results dictionary
        """
        if workers == 'auto':
            workers = multiprocessing.cpu_count()

        if workers == 1:
            return self._process_batch_serial(pdb_ids, steps)
        else:
            return self._process_batch_parallel(pdb_ids, steps, workers)

    def _process_batch_serial(self, pdb_ids: List[str], steps: List[str]) -> Dict:
        """Serial processing (clear output)"""
        results = {
            'total': len(pdb_ids),
            'succeeded': 0,
            'failed': 0,
            'details': []
        }

        for i, pdb_id in enumerate(pdb_ids, 1):
            if not self.verbose:
                print(f"[{i}/{len(pdb_ids)}] {pdb_id.upper()}...", end=' ', flush=True)

            try:
                result = self._process_single(pdb_id, steps)
                results['details'].append(result)

                if result['status'] == 'success':
                    results['succeeded'] += 1
                    status_msg = self._format_status_message(result)
                    if not self.verbose:
                        print(f"✓{status_msg}")
                else:
                    results['failed'] += 1
                    if not self.verbose:
                        print("✗")

            except Exception as e:
                results['failed'] += 1
                results['details'].append({
                    'pdb_id': pdb_id.upper(),
                    'status': 'error',
                    'error': str(e)
                })
                if not self.verbose:
                    print(f"✗ ({str(e)})")

        return results

    def _process_batch_parallel(self, pdb_ids: List[str], steps: List[str], workers: int) -> Dict:
        """Parallel processing (high speed)"""
        results = {
            'total': len(pdb_ids),
            'succeeded': 0,
            'failed': 0,
            'details': []
        }

        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_pdb = {
                executor.submit(self._process_single, pdb_id, steps): pdb_id
                for pdb_id in pdb_ids
            }

            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result = future.result()
                    results['details'].append(result)

                    if result['status'] == 'success':
                        results['succeeded'] += 1
                        status_msg = self._format_status_message(result)
                        print(f"✓ {pdb_id.upper()}{status_msg}")
                    else:
                        results['failed'] += 1
                        print(f"✗ {pdb_id.upper()}")

                except Exception as e:
                    results['failed'] += 1
                    results['details'].append({
                        'pdb_id': pdb_id.upper(),
                        'status': 'error',
                        'error': str(e)
                    })
                    print(f"✗ {pdb_id.upper()}: {e}")

        return results

    def _format_status_message(self, result: Dict) -> str:
        """Format special status message"""
        messages = []

        # Check multi-model concatenation
        if 'concatenate' in result['steps']:
            concat_info = result['steps']['concatenate']
            if concat_info.get('status') == 'success' and concat_info.get('models', 1) > 1:
                n_models = concat_info['models']
                messages.append(f"multi-model merge: {n_models} → 1")

        # Check chain remapping
        if 'standardize' in result['steps']:
            std_info = result['steps']['standardize']
            if std_info.get('status') == 'success' and std_info.get('note') != 'chains_already_standard':
                chain_mapping = std_info.get('chain_mapping', {})
                if chain_mapping:
                    old_chains = '/'.join(sorted(chain_mapping.keys()))
                    new_chains = '/'.join([chain_mapping[k] for k in sorted(chain_mapping.keys())])
                    messages.append(f"chain remap: {old_chains} → {new_chains}")

        if messages:
            return f" ({'; '.join(messages)})"
        return ""

    def _process_single(self, pdb_id: str, steps: List[str]) -> Dict:
        """Process single PDB"""
        pdb_id = pdb_id.upper()
        result = {'pdb_id': pdb_id, 'steps': {}}

        if self.verbose:
            print(f"\n{'='*80}")
            print(f"Processing: {pdb_id}")
            print('='*80)

        # Step 1: Download
        if 'download' in steps:
            if self.verbose:
                print("\n[Step 1] Downloading...")
            result['steps']['download'] = self._step_download(pdb_id)
            if self.verbose and result['steps']['download']['status'] == 'success':
                print(f"  ✓ Downloaded: {result['steps']['download'].get('file')}")

        # Step 2: Concatenate (auto-detect)
        if self.verbose:
            print("\n[Step 2] Checking multi-model...")
        result['steps']['concatenate'] = self._step_concatenate(pdb_id)
        if self.verbose:
            concat_info = result['steps']['concatenate']
            if concat_info['status'] == 'success':
                print(f"  ✓ Merged {concat_info.get('models')} models")
            elif concat_info['status'] == 'single_model':
                print("  ✓ Single model (no merge needed)")

        # Step 3: Fix
        if 'fix' in steps:
            if self.verbose:
                print("\n[Step 3] Fixing structure...")
            result['steps']['fix'] = self._step_fix(pdb_id)
            if self.verbose and result['steps']['fix']['status'] == 'success':
                stats = result['steps']['fix'].get('stats', {})
                print(f"  ✓ Fixed: {stats.get('internal_added', 0)} internal residues added, "
                      f"{stats.get('terminal_skipped', 0)} terminal skipped")

        # Step 4: Standardize
        if 'standardize' in steps:
            if self.verbose:
                print("\n[Step 4] Standardizing chains...")
            result['steps']['standardize'] = self._step_standardize(pdb_id)
            if self.verbose and result['steps']['standardize']['status'] == 'success':
                chain_mapping = result['steps']['standardize'].get('chain_mapping', {})
                if chain_mapping:
                    print(f"  ✓ Chain mapping: {chain_mapping}")
                else:
                    print("  ✓ Chains already standard")

        # Step 5: Trim
        if 'trim' in steps:
            if self.verbose:
                print("\n[Step 5] Trimming by distance...")
            result['steps']['trim'] = self._step_trim(pdb_id)
            if self.verbose and result['steps']['trim']['status'] == 'success':
                kept = result['steps']['trim'].get('kept_residues')
                total = result['steps']['trim'].get('total_residues')
                rate = result['steps']['trim'].get('retention_rate', 0)
                print(f"  ✓ Trimmed: {kept}/{total} residues ({rate:.1f}%)")

        # Overall status
        result['status'] = 'success' if all(
            s.get('status') in ['success', 'single_model', 'skipped']
            for s in result['steps'].values()
        ) else 'partial'

        return result

    def _step_download(self, pdb_id: str) -> Dict:
        """Download step"""
        output_file = self.dirs['raw'] / f"{pdb_id}.pdb"

        if self.skip_existing and output_file.exists():
            return {'status': 'skipped', 'file': str(output_file)}

        result = self.downloader.download_pdb(pdb_id, output_dir=str(self.dirs['raw']))
        return {
            'status': 'success' if result['status'] == 'success' else 'failed',
            'file': result.get('file_path'),
            'assembly': result.get('assembly_id')
        }

    def _step_concatenate(self, pdb_id: str) -> Dict:
        """Multi-model concatenation step (auto-detect)"""
        raw_file = self.dirs['raw'] / f"{pdb_id}.pdb"

        if not raw_file.exists():
            return {'status': 'skipped', 'reason': 'no_raw_file'}

        # Detect multi-model
        n_models = self.concatenator.count_models(str(raw_file))

        if n_models <= 1:
            # Single model, copy to concat directory
            concat_file = self.dirs['concat'] / f"{pdb_id}.pdb"
            if not concat_file.exists() or not self.skip_existing:
                shutil.copy(raw_file, concat_file)
            return {'status': 'single_model', 'file': str(concat_file), 'models': 1}

        # Multi-model, concatenate
        concat_file = self.dirs['concat'] / f"{pdb_id}.pdb"

        if self.skip_existing and concat_file.exists():
            return {'status': 'skipped', 'file': str(concat_file)}

        result = self.concatenator.concatenate_models(
            str(raw_file),
            str(concat_file)
        )
        return {
            'status': 'success',
            'file': str(concat_file),
            'models': n_models,
            'total_atoms': result['total_atoms']
        }

    def _step_fix(self, pdb_id: str) -> Dict:
        """Structure fixing step"""
        input_file = self.dirs['concat'] / f"{pdb_id}.pdb"
        output_file = self.dirs['fixed'] / f"{pdb_id}.pdb"

        if not input_file.exists():
            return {'status': 'failed', 'error': 'input_file_not_found'}

        if self.skip_existing and output_file.exists():
            return {'status': 'skipped', 'file': str(output_file)}

        stats = self.fixer.fix_structure(
            str(input_file),
            str(output_file),
            verbose=False
        )
        return {
            'status': 'success',
            'file': str(output_file),
            'stats': stats
        }

    def _step_standardize(self, pdb_id: str) -> Dict:
        """Chain standardization step"""
        input_file = self.dirs['fixed'] / f"{pdb_id}.pdb"
        output_file = self.dirs['standardized'] / f"{pdb_id}_sd.pdb"

        if not input_file.exists():
            return {'status': 'failed', 'error': 'input_file_not_found'}

        if self.skip_existing and output_file.exists():
            return {'status': 'skipped', 'file': str(output_file)}

        # Use intelligent chain identification
        result = self.standardizer.process_single(
            input_pdb=str(input_file),
            output_pdb=str(output_file),
            task_name=pdb_id
        )

        # Handle "already standard" case (need manual copy)
        if result.status == 'ALREADY_STANDARD':
            if not output_file.exists():
                shutil.copy(input_file, output_file)
            return {
                'status': 'success',
                'file': str(output_file),
                'chain_mapping': result.chain_mapping,
                'note': 'chains_already_standard'
            }

        return {
            'status': 'success' if result.status == 'OK' else 'failed',
            'file': str(output_file),
            'chain_mapping': result.chain_mapping
        }

    def _step_trim(self, pdb_id: str) -> Dict:
        """Distance trimming step"""
        input_file = self.dirs['standardized'] / f"{pdb_id}_sd.pdb"
        output_file = self.dirs['trimmed'] / f"{pdb_id}_sd.pdb"

        if not input_file.exists():
            return {'status': 'failed', 'error': 'input_file_not_found'}

        if self.skip_existing and output_file.exists():
            return {'status': 'skipped', 'file': str(output_file)}

        stats = self.trimmer.trim_structure(
            str(input_file),
            str(output_file),
            verbose=False
        )
        return {
            'status': 'success',
            'file': str(output_file),
            'kept_residues': stats['kept_residues'],
            'total_residues': stats['total_residues'],
            'retention_rate': stats['retention_rate']
        }


def main():
    """CLI entry point for batch PDB processing"""
    parser = argparse.ArgumentParser(
        prog='immunoscope-batch-pdb',
        description='PDB Batch Processing Pipeline: From PDB ID to MD-ready structure',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:

1. Complete workflow for multiple PDBs:
   immunoscope-batch-pdb 1ao7 2ckb 6vma 6vmc

2. Read PDB IDs from file:
   immunoscope-batch-pdb --pdb-list my_pdbs.txt

3. Run partial steps only:
   immunoscope-batch-pdb 1ao7 2ckb --steps fix,standardize,trim

4. Custom output directory:
   immunoscope-batch-pdb 1ao7 --output-dir /path/to/output

5. Generate detailed report:
   immunoscope-batch-pdb --pdb-list pdbs.txt --report processing_report.json

6. Parallel processing (4 workers):
   immunoscope-batch-pdb --pdb-list pdbs.txt --workers 4

7. Auto parallel (use all CPU cores):
   immunoscope-batch-pdb --pdb-list pdbs.txt --workers auto

File format (my_pdbs.txt):
  1ao7
  2ckb
  3dxa
  6vma
  6vmc
        """
    )

    # PDB input
    parser.add_argument(
        'pdb_ids',
        nargs='*',
        help='PDB IDs (e.g., 1ao7 2ckb 6vma)'
    )
    parser.add_argument(
        '--pdb-list',
        type=str,
        help='Text file containing PDB IDs (one per line)'
    )

    # Processing options
    parser.add_argument(
        '--steps',
        type=str,
        default='download,fix,standardize,trim',
        help='Steps to execute (default: download,fix,standardize,trim)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='input/standardizedpdbs',
        help='Output directory (default: input/standardizedpdbs)'
    )
    parser.add_argument(
        '--workers',
        type=str,
        default='1',
        help='Number of parallel workers (1=serial, 2-8=parallel, auto=CPU count, default: 1)'
    )
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        default=True,
        help='Skip existing files (default: True)'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force reprocess all files'
    )

    # Output options
    parser.add_argument(
        '--report',
        type=str,
        help='JSON report output path'
    )
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Collect PDB IDs
    pdb_ids = []
    if args.pdb_ids:
        pdb_ids.extend(args.pdb_ids)
    if args.pdb_list:
        with open(args.pdb_list, 'r') as f:
            pdb_ids.extend([line.strip() for line in f if line.strip() and not line.startswith('#')])

    if not pdb_ids:
        print("Error: No PDB IDs provided")
        parser.print_help()
        sys.exit(1)

    # Parse steps
    steps = [s.strip() for s in args.steps.split(',')]

    # Parse workers
    if args.workers == 'auto':
        workers = 'auto'
    else:
        workers = int(args.workers)

    # Print header
    print("\nPDB Batch Processing Pipeline")
    print("=" * 80)
    print(f"Input: {len(pdb_ids)} PDB(s) ({', '.join([p.upper() for p in pdb_ids[:5]])}{', ...' if len(pdb_ids) > 5 else ''})")
    print(f"Steps: {' → '.join(steps)}")
    print(f"Output: {args.output_dir}/trimmed_sd_pdbs/")
    if workers != 1:
        worker_str = f"{multiprocessing.cpu_count()} (auto)" if workers == 'auto' else workers
        print(f"Workers: {worker_str}")
    print()

    # Initialize pipeline
    pipeline = PDBProcessingPipeline(
        base_dir=args.output_dir,
        skip_existing=not args.force,
        verbose=args.verbose
    )

    # Execute processing
    print("Processing:")
    start_time = time.time()

    results = pipeline.process_batch(
        pdb_ids=pdb_ids,
        steps=steps,
        workers=workers
    )

    elapsed_time = time.time() - start_time

    # Output results
    print("\n" + "=" * 80)
    print("Processing Complete")
    print("=" * 80)
    print(f"Total: {results['total']} | Succeeded: {results['succeeded']} | "
          f"Failed: {results['failed']} | Time: {elapsed_time:.1f}s")

    # Save report
    if args.report:
        results['elapsed_time'] = elapsed_time
        with open(args.report, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nDetailed report saved: {args.report}")

    # Display final files
    print("\nFinal files:")
    for detail in results['details']:
        if detail['status'] == 'success' and 'trim' in detail.get('steps', {}):
            trim_info = detail['steps']['trim']
            if trim_info.get('status') == 'success':
                pdb_id = detail['pdb_id']
                final_file = Path(trim_info['file'])
                if final_file.exists():
                    size_kb = final_file.stat().st_size / 1024
                    kept = trim_info.get('kept_residues', 0)
                    print(f"  ✓ {pdb_id}_sd.pdb ({size_kb:.1f} KB, {kept} residues)")

    print("=" * 80 + "\n")

    # Exit with appropriate code
    sys.exit(0 if results['failed'] == 0 else 1)


if __name__ == '__main__':
    main()
