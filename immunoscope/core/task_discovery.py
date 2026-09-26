"""
ImmunoScope Task Discovery and Manifest Generation

This module provides tools for discovering single-task directories in a batch root,
validating their completeness, and generating structured manifests.

Design Principles:
- Single is a special case of batch
- Standard directory structure enables discoverability
- Manifest enables executability
- Scanning discovers candidate tasks, batch execution relies on manifest

Usage:
    >>> from immunoscope.core import discover_tasks, write_manifest
    >>> report = discover_tasks(Path("/data/batch"))
    >>> write_manifest(report, Path("manifest.jsonl"), format="jsonl")
    >>> print(f"Found {report.num_valid} valid tasks")
"""

from pathlib import Path
from typing import Optional, List, Dict, Any, Iterable
from datetime import datetime
import logging
import yaml
import json
import csv

from .models import TaskInputFiles, DiscoveredTask, DiscoveryReport
from .exceptions import (
    TaskDiscoveryError,
    TaskFileParseError,
    ManifestWriteError
)

logger = logging.getLogger(__name__)


# ============================================================================
# File Extensions
# ============================================================================

STRUCTURE_EXTENSIONS = ['.pdb', '.gro']
TOPOLOGY_EXTENSIONS = ['.tpr', '.top']
TRAJECTORY_EXTENSIONS = ['.xtc', '.trr']

TASK_FILE_NAMES = ['task.yaml', 'immunoscope_task.yaml']
REQUIRED_FILE_TYPES = ('structure', 'topology', 'trajectory')


# ============================================================================
# Task Discoverer
# ============================================================================

class TaskDiscoverer:
    """
    Discovers tasks in a batch root directory.

    Scans first-level subdirectories for candidate single-task directories,
    validates their completeness, and generates discovery reports.
    """

    def __init__(self):
        """Initialize task discoverer."""
        pass

    def discover_tasks(self, root: Path, task_depth: int = 1, required_files: Optional[List[str]] = None) -> DiscoveryReport:
        """
        Discover all tasks in a batch root directory.

        Args:
            root: Path to batch root directory

        Returns:
            DiscoveryReport with all discovered tasks

        Raises:
            TaskDiscoveryError: If root directory does not exist or is not accessible
        """
        # Validate root directory
        if not root.exists():
            raise TaskDiscoveryError(str(root), "directory does not exist")

        if not root.is_dir():
            raise TaskDiscoveryError(str(root), "path is not a directory")

        logger.info(f"Discovering tasks in: {root}")

        # Scan candidate task directories at the requested depth
        discovered_tasks = []
        for task_dir in self._iter_candidate_task_dirs(root, task_depth=task_depth):
            logger.debug(f"Scanning candidate task: {task_dir}")
            task = self._discover_single_task(task_dir, root, required_files=required_files)
            discovered_tasks.append(task)

        # Classify tasks
        valid_tasks = [t for t in discovered_tasks if t.is_valid()]
        invalid_tasks = [t for t in discovered_tasks if t.is_invalid()]
        ambiguous_tasks = [t for t in discovered_tasks if t.is_ambiguous()]

        # Create report
        report = DiscoveryReport(
            source_root=str(root.absolute()),
            discovered_at=datetime.now().isoformat(),
            total_tasks=len(discovered_tasks),
            valid_tasks=valid_tasks,
            invalid_tasks=invalid_tasks,
            ambiguous_tasks=ambiguous_tasks,
            all_tasks=discovered_tasks
        )

        logger.info(f"Discovery complete: {report.num_valid} valid, "
                   f"{report.num_invalid} invalid, {report.num_ambiguous} ambiguous")

        return report

    def _iter_candidate_task_dirs(self, root: Path, task_depth: int = 1) -> Iterable[Path]:
        """Yield candidate task directories at the requested depth."""
        if task_depth < 1:
            raise TaskDiscoveryError(str(root), "task_depth must be >= 1")

        for candidate in sorted(root.rglob('*')):
            if not candidate.is_dir():
                continue
            if len(candidate.relative_to(root).parts) != task_depth:
                continue
            yield candidate

    def discover_tasks_from_list(
        self,
        task_list: List[Dict[str, Any]],
        source_root: str = '<task_list>',
        required_files: Optional[List[str]] = None
    ) -> DiscoveryReport:
        """Build a DiscoveryReport from a manual list of task dictionaries."""
        discovered_tasks = []
        source_root_path = Path(source_root).absolute() if source_root != '<task_list>' else Path.cwd()

        for index, entry in enumerate(task_list, start=1):
            task = self._discover_list_entry(entry, source_root_path, index=index, required_files=required_files)
            discovered_tasks.append(task)

        valid_tasks = [t for t in discovered_tasks if t.is_valid()]
        invalid_tasks = [t for t in discovered_tasks if t.is_invalid()]
        ambiguous_tasks = [t for t in discovered_tasks if t.is_ambiguous()]

        return DiscoveryReport(
            source_root=str(source_root_path),
            discovered_at=datetime.now().isoformat(),
            total_tasks=len(discovered_tasks),
            valid_tasks=valid_tasks,
            invalid_tasks=invalid_tasks,
            ambiguous_tasks=ambiguous_tasks,
            all_tasks=discovered_tasks,
        )

    def _discover_list_entry(self, entry: Dict[str, Any], source_root: Path, index: int, required_files: Optional[List[str]] = None) -> DiscoveredTask:
        """Convert a manual task-list entry into a discovered task."""
        task_id = entry.get('task_id') or entry.get('system_id') or f'task_{index:03d}'
        structure_path = entry.get('structure') or entry.get('structure_path') or entry.get('structure_pdb')
        topology_path = entry.get('topology') or entry.get('topology_path')
        trajectory_path = entry.get('trajectory') or entry.get('trajectory_path') or entry.get('trajectory_raw')

        def _resolve_path(value: Optional[str]) -> Optional[str]:
            if not value:
                return None
            path = Path(value)
            if not path.is_absolute():
                path = source_root / path
            return str(path.absolute())

        task_root = entry.get('task_root')
        if task_root:
            task_root_path = Path(task_root)
            if not task_root_path.is_absolute():
                task_root_path = source_root / task_root_path
        else:
            task_root_path = source_root / task_id

        task = DiscoveredTask(
            task_id=task_id,
            task_root=str(task_root_path.absolute()),
            source_root=str(source_root),
            input_files=TaskInputFiles(
                structure_path=_resolve_path(structure_path),
                topology_path=_resolve_path(topology_path),
                trajectory_path=_resolve_path(trajectory_path),
            ),
            task_file_path=None,
            validation_status='pending',
            validation_messages=[],
            tags=list(entry.get('tags', [])),
            metadata=dict(entry.get('metadata', {})),
        )

        self._validate_list_task(task, required_files=required_files)
        return task

    def _validate_list_task(self, task: DiscoveredTask, required_files: Optional[List[str]] = None) -> None:
        """Validate a task that was created from a manual list entry."""
        messages = []
        required_types = self._normalize_required_files(required_files)
        missing = [file_type for file_type in task.input_files.missing_files() if file_type in required_types]
        for file_type in missing:
            messages.append(f'Missing {file_type} file')

        for label, path_value in [
            ('structure', task.input_files.structure_path),
            ('topology', task.input_files.topology_path),
            ('trajectory', task.input_files.trajectory_path),
        ]:
            if label not in required_types:
                continue
            if path_value and not Path(path_value).exists():
                messages.append(f'{label.capitalize()} path does not exist: {path_value}')

        task.validation_status = 'invalid' if messages else 'valid'
        task.validation_messages = messages

    def _discover_single_task(self, task_dir: Path, source_root: Path, required_files: Optional[List[str]] = None) -> DiscoveredTask:
        """
        Discover and validate a single task directory.

        Args:
            task_dir: Path to task directory
            source_root: Path to batch root directory

        Returns:
            DiscoveredTask object
        """
        task_id = task_dir.name
        task_root = str(task_dir.absolute())

        # Try to find task.yaml
        task_file = self._find_task_file(task_dir)

        # Parse task file if exists
        task_config = {}
        metadata = {}
        tags = []

        if task_file:
            logger.debug(f"Found task file: {task_file.name}")
            try:
                task_config = self._parse_task_file(task_file)
                # Extract task_id override
                if 'task_id' in task_config:
                    task_id = task_config['task_id']
                # Extract metadata and tags
                metadata = task_config.get('metadata', {})
                tags = task_config.get('tags', [])
            except TaskFileParseError as e:
                logger.warning(f"Failed to parse task file for {task_id}: {e}")

        # Discover input files
        if task_file and 'input_files' in task_config:
            # Use explicit paths from task.yaml
            input_files = self._resolve_explicit_files(
                task_dir,
                task_config['input_files']
            )
        else:
            # Auto-discover files
            input_files = self._auto_discover_files(task_dir)

        # Create task object
        task = DiscoveredTask(
            task_id=task_id,
            task_root=task_root,
            source_root=str(source_root.absolute()),
            input_files=input_files,
            task_file_path=str(task_file.absolute()) if task_file else None,
            validation_status="pending",  # Will be updated by validation
            validation_messages=[],
            tags=tags,
            metadata=metadata
        )

        # Validate task
        self._validate_task(task, required_files=required_files)

        return task

    def _find_task_file(self, task_dir: Path) -> Optional[Path]:
        """
        Find task.yaml or immunoscope_task.yaml in task directory.

        Args:
            task_dir: Path to task directory

        Returns:
            Path to task file if found, None otherwise
        """
        for task_file_name in TASK_FILE_NAMES:
            task_file = task_dir / task_file_name
            if task_file.exists() and task_file.is_file():
                return task_file
        return None

    def _parse_task_file(self, task_file: Path) -> Dict[str, Any]:
        """
        Parse task.yaml file.

        Args:
            task_file: Path to task file

        Returns:
            Dictionary with task configuration

        Raises:
            TaskFileParseError: If YAML cannot be parsed
        """
        try:
            with open(task_file, 'r') as f:
                config = yaml.safe_load(f)

            if not isinstance(config, dict):
                raise TaskFileParseError(
                    str(task_file),
                    "task file must contain a dictionary"
                )

            return config

        except yaml.YAMLError as e:
            raise TaskFileParseError(str(task_file), f"invalid YAML: {e}")
        except Exception as e:
            raise TaskFileParseError(str(task_file), str(e))

    def _resolve_explicit_files(
        self,
        task_dir: Path,
        input_files_config: Dict[str, str]
    ) -> TaskInputFiles:
        """
        Resolve file paths from explicit task.yaml configuration.

        Args:
            task_dir: Path to task directory
            input_files_config: Dictionary with file paths from task.yaml

        Returns:
            TaskInputFiles object
        """
        input_files = TaskInputFiles()

        # Resolve structure
        if 'structure' in input_files_config:
            structure_path = task_dir / input_files_config['structure']
            if structure_path.exists():
                input_files.structure_path = str(structure_path.absolute())

        # Resolve topology
        if 'topology' in input_files_config:
            topology_path = task_dir / input_files_config['topology']
            if topology_path.exists():
                input_files.topology_path = str(topology_path.absolute())

        # Resolve trajectory
        if 'trajectory' in input_files_config:
            trajectory_path = task_dir / input_files_config['trajectory']
            if trajectory_path.exists():
                input_files.trajectory_path = str(trajectory_path.absolute())

        return input_files

    def _auto_discover_files(self, task_dir: Path) -> TaskInputFiles:
        """
        Auto-discover input files in task directory.

        Search order:
        1. input/ subdirectory
        2. Task root directory

        Args:
            task_dir: Path to task directory

        Returns:
            TaskInputFiles object
        """
        input_files = TaskInputFiles()

        # Define search locations
        search_locations = []
        input_dir = task_dir / 'input'
        prod_dir = task_dir / 'prod'
        if input_dir.exists() and input_dir.is_dir():
            search_locations.append(input_dir)
        search_locations.append(task_dir)
        if prod_dir.exists() and prod_dir.is_dir():
            search_locations.append(prod_dir)

        # Auto-discover structure file
        input_files.structure_path = self._find_file_by_extension(
            search_locations,
            STRUCTURE_EXTENSIONS,
            role='structure',
            task_dir=task_dir,
        )

        # Auto-discover topology file
        input_files.topology_path = self._find_file_by_extension(
            search_locations,
            TOPOLOGY_EXTENSIONS,
            role='topology',
            task_dir=task_dir,
        )

        # Auto-discover trajectory file
        input_files.trajectory_path = self._find_file_by_extension(
            search_locations,
            TRAJECTORY_EXTENSIONS,
            role='trajectory',
            task_dir=task_dir,
        )

        return input_files

    def _find_file_by_extension(
        self,
        search_locations: List[Path],
        extensions: List[str],
        role: str | None = None,
        task_dir: Path | None = None,
    ) -> Optional[str]:
        """
        Find file with given extensions in search locations.

        Returns the first unique match. If multiple candidates exist,
        returns None (ambiguous).

        Args:
            search_locations: List of directories to search
            extensions: List of file extensions (e.g., ['.pdb', '.gro'])

        Returns:
            Absolute path to file if unique match found, None otherwise
        """
        candidates = []

        for location in search_locations:
            for ext in extensions:
                matches = list(location.glob(f"*{ext}"))
                candidates.extend(matches)

        # Remove duplicates
        candidates = list(set(candidates))

        if len(candidates) == 0:
            return None
        elif len(candidates) == 1:
            return str(candidates[0].absolute())

        selected = self._select_preferred_candidate(candidates, role=role, task_dir=task_dir)
        if selected is not None:
            return str(selected.absolute())

        # Multiple candidates - return None to indicate ambiguity
        return None

    def _select_preferred_candidate(self, candidates: List[Path], role: str | None, task_dir: Path | None) -> Optional[Path]:
        """Select the most likely MD input when a task has auxiliary step files."""
        if not role or task_dir is None:
            return None

        task_base = self._normalize_task_stem(task_dir.name)

        def score(path: Path) -> tuple[int, str]:
            rel_parts = path.relative_to(task_dir).parts if path.is_relative_to(task_dir) else path.parts
            stem = self._normalize_task_stem(path.stem)
            name = path.name.lower()
            parent = path.parent.name.lower()
            value = 100

            if role == 'structure':
                if parent == task_dir.name.lower():
                    value -= 20
                if stem == task_base:
                    value -= 40
                if stem.startswith(task_base) or task_base.startswith(stem):
                    value -= 10
                if name.startswith('step'):
                    value += 40
                if path.suffix.lower() == '.pdb':
                    value -= 5

            elif role == 'topology':
                if parent == 'prod':
                    value -= 30
                if name == 'md.tpr':
                    value -= 50
                elif path.suffix.lower() == '.tpr':
                    value -= 20
                if name.startswith('step'):
                    value += 30

            elif role == 'trajectory':
                if parent == 'prod':
                    value -= 30
                if name == 'md.xtc' and parent == 'prod':
                    value -= 50
                elif path.suffix.lower() == '.xtc':
                    value -= 20
                if name.startswith('step'):
                    value += 30

            return value, "/".join(rel_parts)

        ranked = sorted(candidates, key=score)
        if len(ranked) < 2:
            return ranked[0]
        if score(ranked[0])[0] == score(ranked[1])[0]:
            return None
        return ranked[0]

    def _normalize_task_stem(self, value: str) -> str:
        lowered = value.lower()
        for suffix in ('_sd_run2', '_run2', '_sd'):
            if lowered.endswith(suffix):
                lowered = lowered[: -len(suffix)]
        return lowered

    def _validate_task(self, task: DiscoveredTask, required_files: Optional[List[str]] = None) -> None:
        """
        Validate task completeness and uniqueness.

        Updates task.validation_status and task.validation_messages in-place.

        Args:
            task: DiscoveredTask object to validate
        """
        messages = []

        # Check for missing files
        required_types = self._normalize_required_files(required_files)
        missing = [file_type for file_type in task.input_files.missing_files() if file_type in required_types]
        if missing:
            for file_type in missing:
                messages.append(f"Missing {file_type} file")

        # Check for ambiguity by looking at search results
        # We need to re-scan to detect multiple candidates
        task_dir = Path(task.task_root)
        search_locations = []
        input_dir = task_dir / 'input'
        prod_dir = task_dir / 'prod'
        if input_dir.exists():
            search_locations.append(input_dir)
        search_locations.append(task_dir)
        if prod_dir.exists():
            search_locations.append(prod_dir)

        ambiguous_types = []

        # Check structure ambiguity
        if 'structure' in required_types:
            structure_candidates = self._find_all_files_by_extension(
                search_locations, STRUCTURE_EXTENSIONS
            )
            if len(structure_candidates) > 1 and not task.input_files.structure_path:
                ambiguous_types.append('structure')
                messages.append(f"Ambiguous structure: found {len(structure_candidates)} candidates")

        # Check topology ambiguity
        if 'topology' in required_types:
            topology_candidates = self._find_all_files_by_extension(
                search_locations, TOPOLOGY_EXTENSIONS
            )
            if len(topology_candidates) > 1 and not task.input_files.topology_path:
                ambiguous_types.append('topology')
                messages.append(f"Ambiguous topology: found {len(topology_candidates)} candidates")

        # Check trajectory ambiguity
        if 'trajectory' in required_types:
            trajectory_candidates = self._find_all_files_by_extension(
                search_locations, TRAJECTORY_EXTENSIONS
            )
            if len(trajectory_candidates) > 1 and not task.input_files.trajectory_path:
                ambiguous_types.append('trajectory')
                messages.append(f"Ambiguous trajectory: found {len(trajectory_candidates)} candidates")

        # Determine validation status
        if ambiguous_types:
            task.validation_status = "ambiguous"
        elif missing:
            task.validation_status = "invalid"
        else:
            task.validation_status = "valid"

        task.validation_messages = messages

    def _normalize_required_files(self, required_files: Optional[List[str]]) -> set[str]:
        """Normalize the required file types for validation."""
        if required_files is None:
            return set(REQUIRED_FILE_TYPES)

        normalized = set(required_files)
        invalid = normalized - set(REQUIRED_FILE_TYPES)
        if invalid:
            raise TaskDiscoveryError('<required_files>', f"unsupported required_files: {sorted(invalid)}")

        return normalized

    def _find_all_files_by_extension(
        self,
        search_locations: List[Path],
        extensions: List[str]
    ) -> List[Path]:
        """
        Find ALL files with given extensions (for ambiguity detection).

        Args:
            search_locations: List of directories to search
            extensions: List of file extensions

        Returns:
            List of all matching files
        """
        candidates = []
        for location in search_locations:
            for ext in extensions:
                matches = list(location.glob(f"*{ext}"))
                candidates.extend(matches)

        # Remove duplicates
        return list(set(candidates))


# ============================================================================
# Manifest Writer
# ============================================================================

class ManifestWriter:
    """
    Writes discovery reports to disk in various formats.
    """

    def write_jsonl(self, report: DiscoveryReport, output_path: Path) -> None:
        """
        Write manifest in JSONL format (one task per line).

        Args:
            report: DiscoveryReport to write
            output_path: Path to output JSONL file

        Raises:
            ManifestWriteError: If write fails
        """
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, 'w') as f:
                for task in report.all_tasks:
                    line = json.dumps(task.to_dict())
                    f.write(line + '\n')

            logger.info(f"Wrote JSONL manifest to: {output_path}")

        except Exception as e:
            raise ManifestWriteError(str(output_path), str(e))

    def write_csv(self, report: DiscoveryReport, output_path: Path) -> None:
        """
        Write manifest in CSV format.

        Metadata and tags are serialized as JSON strings.

        Args:
            report: DiscoveryReport to write
            output_path: Path to output CSV file

        Raises:
            ManifestWriteError: If write fails
        """
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)

            fieldnames = [
                'task_id',
                'task_root',
                'source_root',
                'structure_path',
                'topology_path',
                'trajectory_path',
                'task_file_path',
                'validation_status',
                'validation_messages',
                'tags',
                'metadata',
                'discovered_at'
            ]

            with open(output_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()

                for task in report.all_tasks:
                    row = {
                        'task_id': task.task_id,
                        'task_root': task.task_root,
                        'source_root': task.source_root,
                        'structure_path': task.input_files.structure_path or '',
                        'topology_path': task.input_files.topology_path or '',
                        'trajectory_path': task.input_files.trajectory_path or '',
                        'task_file_path': task.task_file_path or '',
                        'validation_status': task.validation_status,
                        'validation_messages': json.dumps(task.validation_messages),
                        'tags': json.dumps(task.tags),
                        'metadata': json.dumps(task.metadata),
                        'discovered_at': task.discovered_at
                    }
                    writer.writerow(row)

            logger.info(f"Wrote CSV manifest to: {output_path}")

        except Exception as e:
            raise ManifestWriteError(str(output_path), str(e))


# ============================================================================
# High-Level API
# ============================================================================

def discover_tasks(root: Path, task_depth: int = 1, required_files: Optional[List[str]] = None) -> DiscoveryReport:
    """
    Discover all tasks in a batch root directory.

    Args:
        root: Path to batch root directory
        task_depth: Directory depth below ``root`` where task directories live.
            Use ``1`` for flat layouts like ``root/task_a`` and ``2`` for nested
            layouts like ``root/pdb_id/variant``.
        required_files: File types that must be present for a task to be valid.
            Defaults to ``['structure', 'topology', 'trajectory']``.

    Returns:
        DiscoveryReport with all discovered tasks

    Raises:
        TaskDiscoveryError: If root directory does not exist

    Example:
        >>> report = discover_tasks(Path("/data/batch"))
        >>> nested = discover_tasks(Path("/data/nested_batch"), task_depth=2)
        >>> print(f"Found {report.num_valid} valid tasks")
        >>> for task in nested.valid_tasks:
        ...     print(f"  - {task.task_id}")
    """
    discoverer = TaskDiscoverer()
    return discoverer.discover_tasks(root, task_depth=task_depth, required_files=required_files)


def discover_tasks_from_list(
    task_list: List[Dict[str, Any]],
    source_root: str = '<task_list>',
    required_files: Optional[List[str]] = None
) -> DiscoveryReport:
    """
    Build a DiscoveryReport from a manual list of task dictionaries.

    Args:
        task_list: Explicit task dictionaries. Supported keys include ``task_id``
            or ``system_id``, ``structure`` or ``structure_path``, ``topology``
            or ``topology_path``, ``trajectory`` or ``trajectory_raw``, plus
            optional ``task_root``, ``tags``, and ``metadata``.
        source_root: Base directory used to resolve relative paths from the task
            list. Keep the default when the task list is already absolute or
            synthetic.
        required_files: File types that must be present for a task to be valid.
            Defaults to ``['structure', 'topology', 'trajectory']``.

    Returns:
        DiscoveryReport with all discovered tasks

    Example:
        >>> report = discover_tasks_from_list([
        ...     {
        ...         "task_id": "1ao7_standard",
        ...         "structure": "data/1ao7/standard/md.gro",
        ...         "topology": "data/1ao7/standard/md.tpr",
        ...         "trajectory_raw": "data/1ao7/standard/md.xtc",
        ...     }
        ... ], source_root="/workspace")
        >>> print(report.total_tasks)
    """
    discoverer = TaskDiscoverer()
    return discoverer.discover_tasks_from_list(task_list, source_root=source_root, required_files=required_files)


def write_manifest(
    report: DiscoveryReport,
    output_path: Path,
    format: str = "jsonl"
) -> None:
    """
    Write discovery report to manifest file.

    Args:
        report: DiscoveryReport to write
        output_path: Path to output file
        format: Format ("jsonl" or "csv")

    Raises:
        ValueError: If format is not supported
        ManifestWriteError: If write fails

    Example:
        >>> report = discover_tasks(Path("/data/batch"))
        >>> write_manifest(report, Path("manifest.jsonl"), format="jsonl")
        >>> write_manifest(report, Path("manifest.csv"), format="csv")
    """
    writer = ManifestWriter()

    if format == "jsonl":
        writer.write_jsonl(report, output_path)
    elif format == "csv":
        writer.write_csv(report, output_path)
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'jsonl' or 'csv'.")
