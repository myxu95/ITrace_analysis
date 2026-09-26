"""
ITrace Core Module - Core processing functionality

Contains contracts, discovery, and execution primitives.
"""

# Core infrastructure (new architecture)
from .context import PipelineContext, ProcessingResult
from .exceptions import (
    ITraceError,
    InputValidationError,
    ProcessingError,
    PipelineError,
    ConfigurationError,
    TaskDiscoveryError,
    TaskValidationError,
    ManifestWriteError,
    TaskFileParseError
)
from .base_node import PipelineNode

# Task discovery (new manifest-based)
from .models import TaskInputFiles, DiscoveredTask, DiscoveryReport
from .task_discovery import (
    TaskDiscoverer,
    ManifestWriter,
    discover_tasks,
    discover_tasks_from_list,
    write_manifest
)
from .task_adapters import (
    discovered_task_to_context,
    discovery_report_to_contexts,
)

# Public task-discovery contract: manifest-based discovery returning DiscoveryReport.
discover_tasks_new = discover_tasks

__all__ = [
    # New architecture
    'PipelineContext',
    'ProcessingResult',
    'PipelineNode',
    # Task discovery (new manifest-based)
    'TaskInputFiles',
    'DiscoveredTask',
    'DiscoveryReport',
    'TaskDiscoverer',
    'ManifestWriter',
    'discover_tasks',
    'discover_tasks_new',
    'discover_tasks_from_list',
    'write_manifest',
    'discovered_task_to_context',
    'discovery_report_to_contexts',
    # Exceptions
    'ITraceError',
    'InputValidationError',
    'ProcessingError',
    'PipelineError',
    'ConfigurationError',
    'TaskDiscoveryError',
    'TaskValidationError',
    'ManifestWriteError',
    'TaskFileParseError',
    'discover_tasks',
]
