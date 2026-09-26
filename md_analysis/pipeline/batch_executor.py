"""
Batch Executor - Universal batch processing for any pipeline.

This module provides a generic batch executor that can run any pipeline
on multiple tasks in parallel.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import logging
from tqdm import tqdm

from ..core.context import PipelineContext
from ..core.models import DiscoveryReport
from ..core.task_adapters import discovery_report_to_contexts
from .base_pipeline import Pipeline

logger = logging.getLogger(__name__)


class BatchExecutor:
    """
    Universal batch executor for pipelines.

    This executor can run any pipeline on multiple tasks in parallel.
    It handles:
    - Parallel execution with configurable workers
    - Error isolation (one task failure doesn't affect others)
    - Progress tracking
    - Result aggregation

    Example:
        >>> from pathlib import Path
        >>> from md_analysis.pipeline import BatchExecutor, PreprocessQualityPipeline
        >>> from md_analysis.core import discover_tasks
        >>>
        >>> report = discover_tasks(Path("/data/simulations"))
        >>> executor = BatchExecutor(max_workers=4)
        >>> results = executor.execute_pipeline(
        ...     report,
        ...     PreprocessQualityPipeline()
        ... )
        >>> summary = executor.summarize_results(results)
        >>> print(f"Success: {summary['successful']}/{summary['total_tasks']}")
    """

    def __init__(self, max_workers: int = 4):
        """
        Initialize batch executor.

        Args:
            max_workers: Maximum number of parallel workers
        """
        self.max_workers = max_workers

    def prepare_tasks(self,
                      tasks: List[PipelineContext] | DiscoveryReport,
                      output_base_dir: Optional[str] = None) -> List[PipelineContext]:
        """
        Normalize execution input into PipelineContext objects.

        Args:
            tasks: Either a list of PipelineContext objects or a DiscoveryReport.
            output_base_dir: Optional output directory root used when adapting a
                DiscoveryReport.

        Returns:
            List of PipelineContext instances ready for execution.
        """
        if isinstance(tasks, DiscoveryReport):
            return discovery_report_to_contexts(
                tasks,
                output_base_dir=output_base_dir,
            )
        return list(tasks)

    def execute_pipeline(self,
                        tasks: List[PipelineContext] | DiscoveryReport,
                        pipeline: Pipeline,
                        show_progress: bool = True,
                        output_base_dir: Optional[str] = None) -> List[PipelineContext]:
        """
        Execute pipeline on multiple tasks.

        Args:
            tasks: Either a list of PipelineContext instances or a DiscoveryReport.
            pipeline: Pipeline to execute
            show_progress: Whether to show progress bar
            output_base_dir: Optional output root used when ``tasks`` is a
                DiscoveryReport.

        Returns:
            List of PipelineContext instances with results

        Example:
            >>> results = executor.execute_pipeline(report, pipeline)
        """
        prepared_tasks = self.prepare_tasks(tasks, output_base_dir=output_base_dir)
        if not prepared_tasks:
            logger.warning("No tasks to process")
            return []

        logger.info(f"Starting batch execution: {len(prepared_tasks)} tasks, {self.max_workers} workers")

        if self.max_workers == 1:
            results = []
            progress_bar = tqdm(prepared_tasks, desc="Processing", disable=not show_progress)
            for task in progress_bar:
                result = self._execute_single_task(pipeline, task)
                results.append(result)
                if show_progress:
                    status = "✓" if not result.has_errors() else "✗"
                    progress_bar.set_postfix_str(f"{status} {result.system_id}")
            logger.info(f"Batch execution completed: {len(results)} results")
            return results

        results = []

        # Use multiprocessing for parallel execution
        with ProcessPoolExecutor(max_workers=self.max_workers) as pool:
            # Submit tasks
            futures = {
                pool.submit(self._execute_single_task, pipeline, task): task
                for task in prepared_tasks
            }

            # Process results as they complete
            progress_bar = tqdm(
                as_completed(futures),
                total=len(prepared_tasks),
                desc="Processing",
                disable=not show_progress
            )

            for future in progress_bar:
                task = futures[future]
                try:
                    result = future.result()
                    results.append(result)

                    # Update progress bar description
                    if show_progress:
                        status = "✓" if not result.has_errors() else "✗"
                        progress_bar.set_postfix_str(f"{status} {result.system_id}")

                except Exception as e:
                    logger.error(f"Task {task.system_id} failed with exception: {e}")
                    task.add_error(f"Executor exception: {str(e)}")
                    results.append(task)

        logger.info(f"Batch execution completed: {len(results)} results")
        return results


    @staticmethod
    def _execute_single_task(pipeline: Pipeline, task: PipelineContext) -> PipelineContext:
        """
        Execute pipeline on a single task.

        This is a static method to enable pickling for multiprocessing.

        Args:
            pipeline: Pipeline to execute
            task: Task context

        Returns:
            Result context
        """
        try:
            return pipeline.execute(task)
        except Exception as e:
            logger.exception(f"Pipeline execution failed for {task.system_id}: {e}")
            task.add_error(f"Pipeline exception: {str(e)}")
            task.should_stop = True
            return task

    def summarize_results(self, results: List[PipelineContext]) -> Dict[str, Any]:
        """
        Summarize batch execution results.

        Args:
            results: List of result contexts

        Returns:
            Summary dictionary

        Example:
            >>> summary = executor.summarize_results(results)
            >>> print(summary)
            {
                'total_tasks': 10,
                'successful': 8,
                'failed': 2,
                'success_rate': 80.0,
                'results': [...]
            }
        """
        total = len(results)
        successful = sum(1 for r in results if not r.has_errors())
        failed = total - successful
        success_rate = (successful / total * 100) if total > 0 else 0

        summary = {
            'total_tasks': total,
            'successful': successful,
            'failed': failed,
            'success_rate': success_rate,
            'results': results
        }

        logger.info(
            f"Summary: {successful}/{total} successful ({success_rate:.1f}%), "
            f"{failed} failed"
        )

        return summary

    def get_failed_tasks(self, results: List[PipelineContext]) -> List[PipelineContext]:
        """
        Get list of failed tasks.

        Args:
            results: List of result contexts

        Returns:
            List of failed task contexts
        """
        failed = [r for r in results if r.has_errors()]
        logger.info(f"Found {len(failed)} failed tasks")
        return failed

    def get_successful_tasks(self, results: List[PipelineContext]) -> List[PipelineContext]:
        """
        Get list of successful tasks.

        Args:
            results: List of result contexts

        Returns:
            List of successful task contexts
        """
        successful = [r for r in results if not r.has_errors()]
        logger.info(f"Found {len(successful)} successful tasks")
        return successful

    def save_summary(self, results: List[PipelineContext], output_file: str):
        """
        Save batch execution summary to file.

        Args:
            results: List of result contexts
            output_file: Output file path (CSV or JSON)

        Example:
            >>> executor.save_summary(results, "batch_summary.csv")
        """
        import pandas as pd
        import json
        from pathlib import Path

        summary_data = []

        for result in results:
            data = {
                'system_id': result.system_id,
                'success': not result.has_errors(),
                'n_errors': len(result.errors),
                'n_warnings': len(result.warnings),
                'n_results': len(result.results)
            }

            # Add result summaries
            for key, value in result.results.items():
                if isinstance(value, dict):
                    for k, v in value.items():
                        if isinstance(v, (int, float, str, bool)):
                            data[f"{key}_{k}"] = v

            summary_data.append(data)

        # Save based on extension
        if output_file.endswith('.csv'):
            df = pd.DataFrame(summary_data)
            df.to_csv(output_file, index=False)
            logger.info(f"Summary saved to CSV: {output_file}")

        elif output_file.endswith('.json'):
            with open(output_file, 'w') as f:
                json.dump(summary_data, f, indent=2)
            logger.info(f"Summary saved to JSON: {output_file}")

        else:
            raise ValueError(f"Unsupported file format: {output_file}")

    def __repr__(self) -> str:
        return f"BatchExecutor(max_workers={self.max_workers})"
