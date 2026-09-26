"""
Task adapters between the manifest-based discovery contract and legacy pipeline contexts.
"""

from pathlib import Path
from typing import List

from .context import PipelineContext
from .models import DiscoveredTask, DiscoveryReport


def discovered_task_to_context(
    task: DiscoveredTask,
    output_base_dir: str | None = None,
) -> PipelineContext:
    """
    Convert a discovered task into a PipelineContext.

    This adapter exists so legacy PipelineContext-based pipelines can consume
    the modern manifest-based discovery contract without importing the
    legacy discovery module directly.
    """
    if output_base_dir:
        output_dir = str(Path(output_base_dir) / task.task_id)
    else:
        output_dir = f"./output/{task.task_id}"

    metadata = dict(task.metadata)
    metadata.update({
        'task_root': task.task_root,
        'source_root': task.source_root,
        'task_file_path': task.task_file_path,
        'tags': list(task.tags),
        'validation_status': task.validation_status,
        'validation_messages': list(task.validation_messages),
        'discovered_at': task.discovered_at,
        'discovery_contract': 'manifest',
    })

    return PipelineContext(
        system_id=task.task_id,
        topology=task.input_files.topology_path,
        trajectory_raw=task.input_files.trajectory_path,
        structure_pdb=task.input_files.structure_path,
        output_dir=output_dir,
        metadata=metadata,
    )


def discovery_report_to_contexts(
    report: DiscoveryReport,
    output_base_dir: str | None = None,
    include_invalid: bool = False,
    include_ambiguous: bool = False,
) -> List[PipelineContext]:
    """
    Convert tasks from a DiscoveryReport into PipelineContext objects.

    By default, only valid tasks are converted.
    """
    tasks = list(report.valid_tasks)
    if include_invalid:
        tasks.extend(report.invalid_tasks)
    if include_ambiguous:
        tasks.extend(report.ambiguous_tasks)

    return [
        discovered_task_to_context(task, output_base_dir=output_base_dir)
        for task in tasks
    ]
