#!/usr/bin/env python3
"""
Batch workflow wrappers over the manifest-based discovery and batch executor contracts.
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.task_discovery import discover_tasks

logger = logging.getLogger(__name__)


def process_md_tasks(simulations_path: str,
                    output_dir: Optional[str] = None,
                    dt: Optional[float] = None,
                    max_workers: Optional[int] = None,
                    gmx_executable: str = 'gmx',
                    method: str = '2step') -> Dict[str, Any]:
    """
    Process multiple MD simulation tasks with automatic file discovery.
    """
    sims_path = Path(simulations_path)
    if not sims_path.exists():
        raise ValueError(f'Simulations path does not exist: {simulations_path}')
    if not sims_path.is_dir():
        raise ValueError(f'Simulations path is not a directory: {simulations_path}')

    resolved_output_dir = output_dir or str(sims_path.parent / f'{sims_path.name}_processed')
    worker_count = max_workers or 4

    logger.info('Starting batch MD processing')
    logger.info(f'Input directory: {simulations_path}')
    logger.info(f'Output directory: {resolved_output_dir}')
    if dt is not None:
        logger.info(f'Frame sampling: {dt} ps')
    logger.info(f'Max workers: {worker_count}')
    logger.info(f'PBC method: {method}')

    report = discover_tasks(
        sims_path,
        required_files=['topology', 'trajectory'],
    )
    if report.num_valid == 0:
        logger.error('No valid MD tasks found for batch processing')
        return {
            'total_tasks': 0,
            'successful': 0,
            'failed': 0,
            'results': [],
            'output_directory': resolved_output_dir,
        }

    from . import BatchExecutor, PreprocessOnlyPipeline

    pipeline = PreprocessOnlyPipeline(
        method=method,
        dt=dt,
        gmx_executable=gmx_executable,
    )
    executor = BatchExecutor(max_workers=worker_count)
    contexts = executor.execute_pipeline(
        report,
        pipeline,
        show_progress=False,
        output_base_dir=resolved_output_dir,
    )

    results = []
    for context in contexts:
        status = 'failed' if context.has_errors() else 'success'
        entry = {
            'task_name': context.system_id,
            'status': status,
            'processed': context.trajectory_processed,
            'output_dir': context.output_dir,
        }
        if context.has_errors():
            entry['error'] = '; '.join(context.errors)
        results.append(entry)

    successful = sum(1 for item in results if item['status'] == 'success')
    failed = len(results) - successful
    summary = {
        'total_tasks': len(results),
        'successful': successful,
        'failed': failed,
        'results': results,
        'output_directory': resolved_output_dir,
    }

    logger.info(f'Batch processing completed: {successful}/{len(results)} tasks processed successfully')
    return summary


def discover_md_tasks(simulations_path: str) -> Dict[str, tuple]:
    """Discover valid MD tasks using the task discovery contract."""
    sims_path = Path(simulations_path)
    if not sims_path.exists():
        raise ValueError(f'Simulations path does not exist: {simulations_path}')

    report = discover_tasks(
        sims_path,
        required_files=['topology', 'trajectory'],
    )
    discovered_tasks = {
        task.task_id: (task.input_files.trajectory_path, task.input_files.topology_path)
        for task in report.valid_tasks
    }

    logger.info(f'Discovered {len(discovered_tasks)} valid MD tasks in {simulations_path}')
    return discovered_tasks


def check_task_status(simulations_path: str) -> Dict[str, Dict[str, Any]]:
    """Check task validity using the task discovery contract."""
    sims_path = Path(simulations_path)
    if not sims_path.exists():
        raise ValueError(f'Simulations path does not exist: {simulations_path}')

    report = discover_tasks(
        sims_path,
        required_files=['topology', 'trajectory'],
    )
    task_status = {}

    for task in report.all_tasks:
        task_status[task.task_id] = {
            'status': task.validation_status,
            'reason': '; '.join(task.validation_messages) if task.validation_messages else 'Task is valid',
            'trajectory': task.input_files.trajectory_path,
            'topology': task.input_files.topology_path,
            'location': str(Path(task.input_files.trajectory_path).parent) if task.input_files.trajectory_path else None,
        }

    logger.info(f'Checked {len(task_status)} tasks in {simulations_path}')
    return task_status


def main():
    """Command line interface for batch processing."""
    import argparse

    parser = argparse.ArgumentParser(
        description='ImmunoScope Batch Processing - Process multiple MD simulation tasks',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m immunoscope.pipeline.batch_workflow /data/simulations
    python -m immunoscope.pipeline.batch_workflow /data/simulations --output /data/processed
    python -m immunoscope.pipeline.batch_workflow /data/simulations --dt 10.0 --workers 4
    python -m immunoscope.pipeline.batch_workflow /data/simulations --check-only
        """
    )

    parser.add_argument('simulations_path', help='Path containing multiple task folders with MD files')
    parser.add_argument('--output', '-o', help="Output directory (default: simulations_path + '_processed')")
    parser.add_argument('--dt', type=float, help='Time interval for frame sampling in ps')
    parser.add_argument('--workers', '-w', type=int, help='Maximum number of parallel workers')
    parser.add_argument('--gmx', default='gmx', help='GROMACS executable command (default: gmx)')
    parser.add_argument('--check-only', action='store_true', help="Only check task status, don't process")
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose logging')

    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        if args.check_only:
            print(f'Checking MD tasks in: {args.simulations_path}')
            print()

            status = check_task_status(args.simulations_path)
            valid_count = sum(1 for info in status.values() if info['status'] == 'valid')
            total_count = len(status)

            print(f'Task Status Summary: {valid_count}/{total_count} valid tasks')
            print('=' * 60)

            for task_name, info in status.items():
                status_icon = '✅' if info['status'] == 'valid' else '❌'
                print(f"{status_icon} {task_name}: {info['reason']}")
                if info['status'] == 'valid':
                    print(f"   📁 {info['trajectory']}")
                    print(f"   📁 {info['topology']}")
                print()
        else:
            results = process_md_tasks(
                simulations_path=args.simulations_path,
                output_dir=args.output,
                dt=args.dt,
                max_workers=args.workers,
                gmx_executable=args.gmx,
                method='2step'
            )

            print()
            print('Batch Processing Results:')
            print('=' * 40)
            print(f"Total tasks: {results['total_tasks']}")
            print(f"Successful: {results['successful']}")
            print(f"Failed: {results['failed']}")
            print(f"Output directory: {results['output_directory']}")
    except Exception as e:
        logger.error(f'Error: {e}')
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
