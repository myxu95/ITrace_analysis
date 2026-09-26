"""
Pipeline Context - Unified data passing mechanism

This module provides the core data structure for passing data between
pipeline nodes, eliminating the need for file path string manipulation.
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, Any, Optional, List
from pathlib import Path
import json
import copy as copy_module
import logging

logger = logging.getLogger(__name__)


@dataclass
class PipelineContext:
    """
    Pipeline execution context - manages all data and state.

    This is the central data structure that flows through the pipeline,
    carrying input files, intermediate results, and metadata. It eliminates
    the need for modules to guess file paths or rely on naming conventions.

    Attributes:
        system_id: Unique identifier for this task/system
        topology: Path to topology file (e.g., md.tpr)
        trajectory_raw: Path to raw trajectory file (e.g., md.xtc)
        structure_pdb: Optional PDB structure file
        trajectory_processed: Path to PBC-processed trajectory
        trajectory_aligned: Path to aligned trajectory
        selections: Atom selection strings (e.g., {'protein': 'protein'})
        results: Analysis results (e.g., {'rmsd': {...}, 'rmsf': {...}})
        metadata: Additional metadata
        output_dir: Base output directory for this task
        errors: List of error messages
        warnings: List of warning messages
        should_stop: Flag to stop pipeline execution
        temporary_files: List of temporary files for cleanup

    Example:
        >>> context = PipelineContext(
        ...     system_id="1ao7",
        ...     topology="data/1ao7/md.tpr",
        ...     trajectory_raw="data/1ao7/md.xtc"
        ... )
        >>> context.add_result('rmsd', {'mean': 2.34, 'std': 0.56})
        >>> output_path = context.get_output_path("rmsd.xvg")
    """

    # System identification
    system_id: str

    # Input files
    topology: str
    trajectory_raw: str
    structure_pdb: Optional[str] = None

    # Processed files
    trajectory_processed: Optional[str] = None
    trajectory_aligned: Optional[str] = None
    index_file: Optional[str] = None

    # Atom selections
    selections: Dict[str, str] = field(default_factory=dict)

    # Analysis results
    results: Dict[str, Any] = field(default_factory=dict)

    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)

    # Output configuration
    output_dir: Optional[str] = None

    # Error tracking
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # Flow control
    should_stop: bool = False

    # Temporary file tracking (for cleanup)
    temporary_files: List[str] = field(default_factory=list)

    # Index manager (lazy initialization)
    _index_manager: Optional['IndexManager'] = field(default=None, init=False, repr=False)

    @property
    def index_manager(self) -> 'IndexManager':
        """
        Get IndexManager instance for this context.

        The IndexManager is lazily initialized on first access.

        Returns:
            IndexManager instance

        Example:
            >>> phla_group_id = context.index_manager.get_group_id('pHLA')
            >>> tcr_group_id = context.index_manager.get_group_id('TCR')
        """
        if self._index_manager is None:
            from md_analysis.analysis.topology import IndexManager
            self._index_manager = IndexManager(self)
        return self._index_manager

    def get_output_path(self, filename: str) -> str:
        """
        Get full path for an output file.

        Automatically creates the output directory if it doesn't exist.

        Args:
            filename: Name of the output file

        Returns:
            Full path to the output file

        Example:
            >>> context.get_output_path("rmsd.xvg")
            './output/1ao7/rmsd.xvg'
        """
        if self.output_dir is None:
            self.output_dir = f"./output/{self.system_id}"

        output_path = Path(self.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        return str(output_path / filename)

    def get_subdir_path(self, subdir: str, filename: str) -> str:
        """
        Get path for file in a subdirectory.

        This is the base method for all specialized path getters.

        Args:
            subdir: Subdirectory path (can be nested like "analysis/rmsd")
            filename: Name of the output file

        Returns:
            Full path to the output file

        Example:
            >>> context.get_subdir_path("analysis/rmsd", "rmsd_protein.xvg")
            './output/1ao7/analysis/rmsd/rmsd_protein.xvg'
        """
        if self.output_dir is None:
            self.output_dir = f"./output/{self.system_id}"

        subdir_path = Path(self.output_dir) / subdir
        subdir_path.mkdir(parents=True, exist_ok=True)

        return str(subdir_path / filename)

    def get_preprocessing_path(self, filename: str) -> str:
        """
        Get path for preprocessing output file.

        Args:
            filename: Name of the preprocessing file

        Returns:
            Full path in preprocessing subdirectory

        Example:
            >>> context.get_preprocessing_path("processed.xtc")
            './output/1ao7/preprocessing/processed.xtc'
        """
        return self.get_subdir_path("preprocessing", filename)

    def get_analysis_path(self, analysis_type: str, filename: str) -> str:
        """
        Get path for analysis result file.

        Args:
            analysis_type: Type of analysis (e.g., "rmsd", "rmsf", "angles")
            filename: Name of the result file

        Returns:
            Full path in analysis subdirectory

        Example:
            >>> context.get_analysis_path("rmsd", "rmsd_protein.xvg")
            './output/1ao7/analysis/rmsd/rmsd_protein.xvg'
            >>> context.get_analysis_path("angles", "docking_angles.csv")
            './output/1ao7/analysis/angles/docking_angles.csv'
        """
        return self.get_subdir_path(f"analysis/{analysis_type}", filename)

    def get_plot_path(self, filename: str) -> str:
        """
        Get path for plot file.

        Args:
            filename: Name of the plot file

        Returns:
            Full path in plots subdirectory

        Example:
            >>> context.get_plot_path("rmsd_overview.png")
            './output/1ao7/plots/rmsd_overview.png'
        """
        return self.get_subdir_path("plots", filename)

    def get_quality_path(self, filename: str) -> str:
        """
        Get path for quality control file.

        Args:
            filename: Name of the quality file

        Returns:
            Full path in quality subdirectory

        Example:
            >>> context.get_quality_path("energy_quality.json")
            './output/1ao7/quality/energy_quality.json'
        """
        return self.get_subdir_path("quality", filename)

    def get_index_path(self, filename: str) -> str:
        """
        Get path for GROMACS index file.

        Args:
            filename: Name of the index file

        Returns:
            Full path in indices subdirectory

        Example:
            >>> context.get_index_path("protein.ndx")
            './output/1ao7/indices/protein.ndx'
        """
        return self.get_subdir_path("indices", filename)

    def get_temp_path(self, filename: str) -> str:
        """
        Get path for temporary file.

        Automatically registers the file for cleanup.

        Args:
            filename: Name of the temporary file

        Returns:
            Full path in temp subdirectory

        Example:
            >>> temp_file = context.get_temp_path("intermediate.xtc")
            >>> # File is automatically added to context.temporary_files
        """
        temp_path = self.get_subdir_path("temp", filename)
        if temp_path not in self.temporary_files:
            self.temporary_files.append(temp_path)
        return temp_path

    def copy(self) -> 'PipelineContext':
        """
        Create a deep copy of this context.

        Useful for parallel node execution where each node needs
        its own independent context.

        Returns:
            Deep copy of this context

        Example:
            >>> context_copy = context.copy()
            >>> context_copy.system_id = "modified"
            >>> context.system_id  # Original unchanged
            '1ao7'
        """
        return copy_module.deepcopy(self)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert context to dictionary.

        Returns:
            Dictionary representation of the context

        Example:
            >>> data = context.to_dict()
            >>> data['system_id']
            '1ao7'
        """
        return asdict(self)

    def save(self, filepath: str):
        """
        Save context to JSON file.

        Args:
            filepath: Path to save the JSON file

        Example:
            >>> context.save("./output/1ao7/context.json")
        """
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

        logger.info(f"Context saved to {filepath}")

    @classmethod
    def load(cls, filepath: str) -> 'PipelineContext':
        """
        Load context from JSON file.

        Args:
            filepath: Path to the JSON file

        Returns:
            PipelineContext instance

        Example:
            >>> context = PipelineContext.load("./output/1ao7/context.json")
        """
        with open(filepath, 'r') as f:
            data = json.load(f)

        logger.info(f"Context loaded from {filepath}")
        return cls(**data)

    def add_result(self, key: str, value: Any):
        """
        Add analysis result.

        Args:
            key: Result identifier (e.g., 'rmsd', 'rmsf')
            value: Result data (any JSON-serializable object)

        Example:
            >>> context.add_result('rmsd', {
            ...     'mean': 2.34,
            ...     'std': 0.56,
            ...     'output_file': './output/1ao7/rmsd.xvg'
            ... })
        """
        self.results[key] = value
        logger.debug(f"Added result: {key}")

    def get_result(self, key: str, default: Any = None) -> Any:
        """
        Get analysis result.

        Args:
            key: Result identifier
            default: Default value if key not found

        Returns:
            Result data or default value

        Example:
            >>> rmsd_data = context.get_result('rmsd')
            >>> mean_rmsd = rmsd_data['mean']
        """
        return self.results.get(key, default)

    def has_errors(self) -> bool:
        """
        Check if context has any errors.

        Returns:
            True if errors exist, False otherwise
        """
        return len(self.errors) > 0

    def add_error(self, error: str):
        """
        Add error message.

        Args:
            error: Error message

        Example:
            >>> context.add_error("RMSD calculation failed: missing topology")
        """
        self.errors.append(error)
        logger.error(f"Error added to context: {error}")

    def add_warning(self, warning: str):
        """
        Add warning message.

        Args:
            warning: Warning message

        Example:
            >>> context.add_warning("Low number of frames detected: 500")
        """
        self.warnings.append(warning)
        logger.warning(f"Warning added to context: {warning}")

    def set_selection(self, name: str, selection_string: str):
        """
        Set atom selection string.

        Args:
            name: Selection identifier (e.g., 'protein', 'backbone')
            selection_string: MDAnalysis selection string

        Example:
            >>> context.set_selection('protein', 'protein')
            >>> context.set_selection('backbone', 'name CA C N O')
        """
        self.selections[name] = selection_string
        logger.debug(f"Selection set: {name} = {selection_string}")

    def get_selection(self, name: str, default: str = "protein") -> str:
        """
        Get atom selection string.

        Args:
            name: Selection identifier
            default: Default selection if not found

        Returns:
            Selection string
        """
        return self.selections.get(name, default)

    def __repr__(self) -> str:
        """String representation for debugging."""
        status = "ERROR" if self.has_errors() else "OK"
        n_results = len(self.results)
        return (f"PipelineContext(system_id='{self.system_id}', "
                f"status={status}, results={n_results})")


@dataclass
class ProcessingResult:
    """
    Standardized result from processing modules.

    This is the standard return type for all Core Modules (Layer 1),
    ensuring consistent output format.

    Attributes:
        success: Whether the processing succeeded
        output_file: Main output file path
        temporary_files: List of temporary files created
        processing_stats: Statistics (e.g., n_frames, processing_time)
        metadata: Additional metadata
        error_message: Error message if failed

    Example:
        >>> result = ProcessingResult(
        ...     success=True,
        ...     output_file="./output/1ao7/rmsd.xvg",
        ...     processing_stats={'n_frames': 1000, 'time_sec': 45.2},
        ...     metadata={'version': '1.0.0'}
        ... )
    """

    success: bool
    output_file: str
    temporary_files: Optional[List[str]] = None
    processing_stats: Optional[Dict[str, Any]] = None
    metadata: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    def __repr__(self) -> str:
        """String representation."""
        status = "SUCCESS" if self.success else "FAILED"
        return f"ProcessingResult(status={status}, output={self.output_file})"
