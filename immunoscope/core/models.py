"""
ImmunoScope Core Data Models

Defines data structures for task discovery and manifest generation.
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, Any
from datetime import datetime
import json


@dataclass
class TaskInputFiles:
    """
    Represents the input files for a single task.

    Attributes:
        structure_path: Path to structure file (.pdb, .gro)
        topology_path: Path to topology file (.tpr, .top)
        trajectory_path: Path to trajectory file (.xtc, .trr)
    """
    structure_path: Optional[str] = None
    topology_path: Optional[str] = None
    trajectory_path: Optional[str] = None

    def is_complete(self) -> bool:
        """Check if all core files are present."""
        return all([
            self.structure_path is not None,
            self.topology_path is not None,
            self.trajectory_path is not None
        ])

    def missing_files(self) -> list[str]:
        """Return list of missing file types."""
        missing = []
        if self.structure_path is None:
            missing.append("structure")
        if self.topology_path is None:
            missing.append("topology")
        if self.trajectory_path is None:
            missing.append("trajectory")
        return missing

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class DiscoveredTask:
    """
    Represents a discovered task with validation status.

    Attributes:
        task_id: Unique identifier for the task
        task_root: Absolute path to task directory
        source_root: Absolute path to batch root directory
        input_files: TaskInputFiles object
        task_file_path: Path to task.yaml if present
        validation_status: One of "valid", "invalid", "ambiguous"
        validation_messages: List of validation messages/warnings
        tags: List of tags from task metadata
        metadata: Additional metadata from task file
        discovered_at: ISO 8601 timestamp of discovery
    """
    task_id: str
    task_root: str
    source_root: str
    input_files: TaskInputFiles
    validation_status: str  # "valid", "invalid", "ambiguous"
    validation_messages: list[str] = field(default_factory=list)
    task_file_path: Optional[str] = None
    tags: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    discovered_at: str = field(default_factory=lambda: datetime.now().isoformat())

    def is_valid(self) -> bool:
        """Check if task is valid."""
        return self.validation_status == "valid"

    def is_invalid(self) -> bool:
        """Check if task is invalid."""
        return self.validation_status == "invalid"

    def is_ambiguous(self) -> bool:
        """Check if task is ambiguous."""
        return self.validation_status == "ambiguous"

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        data = asdict(self)
        # Convert input_files to dict
        data['input_files'] = self.input_files.to_dict()
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> 'DiscoveredTask':
        """Create DiscoveredTask from dictionary."""
        # Reconstruct TaskInputFiles
        input_files_data = data.pop('input_files', {})
        input_files = TaskInputFiles(**input_files_data)

        return cls(
            input_files=input_files,
            **data
        )


@dataclass
class DiscoveryReport:
    """
    Represents the complete discovery report for a batch.

    Attributes:
        source_root: Absolute path to batch root directory
        discovered_at: ISO 8601 timestamp of discovery
        total_tasks: Total number of tasks discovered
        valid_tasks: List of valid tasks
        invalid_tasks: List of invalid tasks
        ambiguous_tasks: List of ambiguous tasks
        all_tasks: List of all tasks
    """
    source_root: str
    discovered_at: str
    total_tasks: int
    valid_tasks: list[DiscoveredTask]
    invalid_tasks: list[DiscoveredTask]
    ambiguous_tasks: list[DiscoveredTask]
    all_tasks: list[DiscoveredTask]

    @property
    def num_valid(self) -> int:
        """Number of valid tasks."""
        return len(self.valid_tasks)

    @property
    def num_invalid(self) -> int:
        """Number of invalid tasks."""
        return len(self.invalid_tasks)

    @property
    def num_ambiguous(self) -> int:
        """Number of ambiguous tasks."""
        return len(self.ambiguous_tasks)

    @property
    def success_rate(self) -> float:
        """Success rate (valid / total)."""
        if self.total_tasks == 0:
            return 0.0
        return self.num_valid / self.total_tasks

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'source_root': self.source_root,
            'discovered_at': self.discovered_at,
            'total_tasks': self.total_tasks,
            'num_valid': self.num_valid,
            'num_invalid': self.num_invalid,
            'num_ambiguous': self.num_ambiguous,
            'success_rate': self.success_rate,
            'valid_tasks': [t.to_dict() for t in self.valid_tasks],
            'invalid_tasks': [t.to_dict() for t in self.invalid_tasks],
            'ambiguous_tasks': [t.to_dict() for t in self.ambiguous_tasks],
            'all_tasks': [t.to_dict() for t in self.all_tasks]
        }

    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> 'DiscoveryReport':
        """Create DiscoveryReport from dictionary."""
        # Reconstruct task lists
        valid_tasks = [DiscoveredTask.from_dict(t) for t in data.get('valid_tasks', [])]
        invalid_tasks = [DiscoveredTask.from_dict(t) for t in data.get('invalid_tasks', [])]
        ambiguous_tasks = [DiscoveredTask.from_dict(t) for t in data.get('ambiguous_tasks', [])]
        all_tasks = [DiscoveredTask.from_dict(t) for t in data.get('all_tasks', [])]

        return cls(
            source_root=data['source_root'],
            discovered_at=data['discovered_at'],
            total_tasks=data['total_tasks'],
            valid_tasks=valid_tasks,
            invalid_tasks=invalid_tasks,
            ambiguous_tasks=ambiguous_tasks,
            all_tasks=all_tasks
        )

    def print_summary(self) -> None:
        """Print a human-readable summary."""
        print("=" * 80)
        print("Task Discovery Summary")
        print("=" * 80)
        print(f"Source root: {self.source_root}")
        print(f"Discovered at: {self.discovered_at}")
        print(f"Total tasks: {self.total_tasks}")
        print(f"  Valid: {self.num_valid} ({self.success_rate*100:.1f}%)")
        print(f"  Invalid: {self.num_invalid}")
        print(f"  Ambiguous: {self.num_ambiguous}")
        print("=" * 80)
