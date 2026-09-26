"""
Standardized Index Generation Module for GROMACS

This module provides a standardized interface for generating GROMACS index files
following the 6 design principles:
1. Clear Inputs (dataclass)
2. Clear Outputs (structured result)
3. Clear Side Effects (tracked)
4. Clear Errors (specific exceptions)
5. Testable (mockable dependencies)
6. Schedulable (progress callbacks, cancellable)

Author: ITrace Development Team
Date: 2026-03-18
"""

import subprocess
import logging
import threading
from pathlib import Path
from typing import Optional, Dict, List, Callable, Any
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

import MDAnalysis as mda

logger = logging.getLogger(__name__)


# ============================================================================
# Error Definitions (Principle 4: Clear Errors)
# ============================================================================

class IndexGenerationError(Exception):
    """Base exception for index generation errors."""
    pass


class InvalidInputError(IndexGenerationError):
    """Invalid input parameter error (user error)."""
    def __init__(self, param_name: str, value: Any, reason: str):
        super().__init__(f"Invalid {param_name}: {reason} (got: {value})")
        self.param_name = param_name
        self.value = value
        self.reason = reason


class TopologyFileError(IndexGenerationError):
    """Topology file access/format error."""
    def __init__(self, file_path: str, reason: str):
        super().__init__(f"Topology file error ({file_path}): {reason}")
        self.file_path = file_path


class GROMACSCommandError(IndexGenerationError):
    """GROMACS command execution error."""
    def __init__(self, command: str, stderr: str, suggestion: str = None):
        msg = f"GROMACS command failed: {command}\nError: {stderr}"
        if suggestion:
            msg += f"\nSuggestion: {suggestion}"
        super().__init__(msg)
        self.command = command
        self.stderr = stderr


class ComponentNotFoundError(IndexGenerationError):
    """Requested component not found in topology."""
    def __init__(self, component: str, available: List[str]):
        super().__init__(
            f"Component '{component}' not found. "
            f"Available: {', '.join(available)}"
        )
        self.component = component
        self.available = available


# ============================================================================
# Input Definitions (Principle 1: Clear Inputs)
# ============================================================================

class IndexGenerationMethod(Enum):
    """Index generation method."""
    PDB_BASED = "pdb_based"          # Use PDB chain IDs
    TOPOLOGY_BASED = "topology_based"  # Use topology chain splitting
    SEQUENCE_BASED = "sequence_based"  # Use sequence matching (CDR)
    CUSTOM = "custom"                 # User-defined selection


@dataclass
class ComponentDefinition:
    """Definition of a component to generate index for."""
    name: str                        # Component name (e.g., "pHLA", "TCR")
    chains: Optional[List[str]] = None  # Chain IDs (for PDB-based)
    group_ids: Optional[List[int]] = None  # Group IDs (for topology-based)
    selection: Optional[str] = None     # Custom MDAnalysis selection
    sequence: Optional[str] = None      # Sequence to match (for CDR)
    ca_only: bool = False               # Select only CA atoms


@dataclass
class IndexGenerationInput:
    """
    Standardized input for index generation.

    Attributes:
        topology: Path to topology file (.tpr, .gro, .pdb)
        method: Index generation method
        components: List of components to generate
        output_file: Output index file path
        standardized_pdb: Optional pre-standardized PDB file
        gmx_executable: GROMACS executable (default: "gmx")
        auto_standardize: Auto-standardize PDB if method is PDB_BASED
    """
    topology: str
    method: IndexGenerationMethod
    components: List[ComponentDefinition]
    output_file: str

    # Optional parameters
    standardized_pdb: Optional[str] = None
    gmx_executable: str = "gmx"
    auto_standardize: bool = True

    def validate(self) -> None:
        """Validate input parameters."""
        # Check topology file exists
        topology_path = Path(self.topology)
        if not topology_path.exists():
            raise InvalidInputError(
                'topology',
                self.topology,
                'File does not exist'
            )

        # Check topology file extension
        valid_extensions = ['.tpr', '.gro', '.pdb']
        if topology_path.suffix.lower() not in valid_extensions:
            raise InvalidInputError(
                'topology',
                self.topology,
                f'Invalid extension. Must be one of: {", ".join(valid_extensions)}'
            )

        # Check components not empty
        if not self.components:
            raise InvalidInputError(
                'components',
                self.components,
                'At least one component must be specified'
            )

        # Check output file path is writable
        output_path = Path(self.output_file)
        output_dir = output_path.parent
        if not output_dir.exists():
            try:
                output_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise InvalidInputError(
                    'output_file',
                    self.output_file,
                    f'Cannot create output directory: {e}'
                )

        # Validate components based on method
        for comp in self.components:
            if self.method == IndexGenerationMethod.PDB_BASED:
                if not comp.chains:
                    raise InvalidInputError(
                        'components',
                        comp.name,
                        f'Component "{comp.name}" requires chains for PDB-based method'
                    )
            elif self.method == IndexGenerationMethod.TOPOLOGY_BASED:
                if not comp.group_ids:
                    raise InvalidInputError(
                        'components',
                        comp.name,
                        f'Component "{comp.name}" requires group_ids for topology-based method'
                    )
            elif self.method == IndexGenerationMethod.SEQUENCE_BASED:
                if not comp.sequence:
                    raise InvalidInputError(
                        'components',
                        comp.name,
                        f'Component "{comp.name}" requires sequence for sequence-based method'
                    )
            elif self.method == IndexGenerationMethod.CUSTOM:
                if not comp.selection:
                    raise InvalidInputError(
                        'components',
                        comp.name,
                        f'Component "{comp.name}" requires selection for custom method'
                    )


# ============================================================================
# Output Definitions (Principle 2: Clear Outputs)
# ============================================================================

@dataclass
class ComponentIndexInfo:
    """Information about a generated component index."""
    name: str                        # Component name
    group_id: int                    # Group ID in index file
    atom_count: int                  # Number of atoms
    residue_count: Optional[int] = None  # Number of residues (if available)
    chains: Optional[List[str]] = None   # Chain IDs (if PDB-based)


@dataclass
class IndexGenerationResult:
    """
    Standardized output for index generation.

    Attributes:
        success: Whether generation succeeded
        output_file: Path to generated index file
        components: List of generated component information
        processing_stats: Statistics about processing
        metadata: Additional metadata
        error_message: Error message if failed
        temporary_files: List of temporary files created
    """
    success: bool
    output_file: Optional[str] = None
    components: List[ComponentIndexInfo] = field(default_factory=list)

    # Statistics
    processing_stats: Dict[str, Any] = field(default_factory=dict)
    # Example: {'n_components': 5, 'total_atoms': 5000, 'processing_time_sec': 1.2}

    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)
    # Example: {'timestamp': '2026-03-18T10:00:00', 'method': 'pdb_based'}

    # Error information
    error_message: Optional[str] = None

    # Side effects
    temporary_files: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'success': self.success,
            'output_file': self.output_file,
            'components': [
                {
                    'name': c.name,
                    'group_id': c.group_id,
                    'atom_count': c.atom_count,
                    'residue_count': c.residue_count,
                    'chains': c.chains
                }
                for c in self.components
            ],
            'processing_stats': self.processing_stats,
            'metadata': self.metadata,
            'error_message': self.error_message,
            'temporary_files': self.temporary_files
        }


# ============================================================================
# Side Effect Tracking (Principle 3: Clear Side Effects)
# ============================================================================

@dataclass
class SideEffectTracker:
    """Track all side effects during index generation."""
    files_created: List[str] = field(default_factory=list)
    files_deleted: List[str] = field(default_factory=list)
    commands_executed: List[Dict] = field(default_factory=list)

    def track_file_creation(self, filepath: str):
        """Record file creation."""
        self.files_created.append(filepath)
        logger.debug(f"File created: {filepath}")

    def track_file_deletion(self, filepath: str):
        """Record file deletion."""
        self.files_deleted.append(filepath)
        logger.debug(f"File deleted: {filepath}")

    def track_command(self, command: List[str], returncode: int, stderr: str = ""):
        """Record command execution."""
        self.commands_executed.append({
            'command': ' '.join(command),
            'returncode': returncode,
            'stderr': stderr,
            'timestamp': datetime.now().isoformat()
        })
        logger.debug(f"Command executed: {' '.join(command)}")

    def get_summary(self) -> Dict[str, Any]:
        """Get summary of side effects."""
        return {
            'files_created': len(self.files_created),
            'files_deleted': len(self.files_deleted),
            'commands_executed': len(self.commands_executed),
            'details': {
                'files_created': self.files_created,
                'files_deleted': self.files_deleted,
                'commands': [c['command'] for c in self.commands_executed]
            }
        }


# ============================================================================
# Index Generator (Principles 5 & 6: Testable & Schedulable)
# ============================================================================

class IndexGenerator:
    """
    Standardized GROMACS index file generator.

    This class follows the 6 design principles:
    1. Clear Inputs: Uses IndexGenerationInput dataclass
    2. Clear Outputs: Returns IndexGenerationResult
    3. Clear Side Effects: Tracks all file operations and commands
    4. Clear Errors: Specific exception types with context
    5. Testable: Supports dependency injection and mocking
    6. Schedulable: Progress callbacks and cancellation support
    """

    def __init__(self):
        """Initialize index generator."""
        self._cancel_flag = threading.Event()
        self.progress_callback: Optional[Callable[[float, str], None]] = None
        self.side_effects = SideEffectTracker()

    def set_progress_callback(self, callback: Callable[[float, str], None]):
        """
        Set progress callback function.

        Args:
            callback: Function(progress: float, message: str)
                      progress: 0.0-1.0
                      message: Current status message
        """
        self.progress_callback = callback

    def cancel(self):
        """Cancel ongoing operation."""
        self._cancel_flag.set()
        logger.info("Index generation cancelled by user")

    def is_cancelled(self) -> bool:
        """Check if operation is cancelled."""
        return self._cancel_flag.is_set()

    def _report_progress(self, progress: float, message: str):
        """Report progress to callback."""
        if self.progress_callback:
            self.progress_callback(progress, message)
        logger.info(f"[{progress*100:.0f}%] {message}")

    def generate(self, input_params: IndexGenerationInput) -> IndexGenerationResult:
        """
        Generate GROMACS index file.

        Args:
            input_params: Index generation input parameters

        Returns:
            IndexGenerationResult with generation results

        Example:
            >>> generator = IndexGenerator()
            >>> input_params = IndexGenerationInput(
            ...     topology="md.tpr",
            ...     method=IndexGenerationMethod.PDB_BASED,
            ...     components=[
            ...         ComponentDefinition(name="pHLA", chains=['A', 'B', 'C']),
            ...         ComponentDefinition(name="TCR", chains=['D', 'E'])
            ...     ],
            ...     output_file="components.ndx"
            ... )
            >>> result = generator.generate(input_params)
            >>> if result.success:
            ...     print(f"Generated: {result.output_file}")
        """
        start_time = datetime.now()
        self.side_effects = SideEffectTracker()  # Reset tracker

        try:
            # Step 1: Validate input
            self._report_progress(0.0, "Validating input parameters...")
            input_params.validate()

            if self.is_cancelled():
                return IndexGenerationResult(
                    success=False,
                    error_message="Operation cancelled by user"
                )

            # Step 2: Prepare topology
            self._report_progress(0.2, "Preparing topology file...")
            topology_file = self._prepare_topology(input_params)

            if self.is_cancelled():
                return IndexGenerationResult(
                    success=False,
                    error_message="Operation cancelled by user"
                )

            # Step 3: Generate index based on method
            self._report_progress(0.4, f"Generating index using {input_params.method.value}...")

            if input_params.method == IndexGenerationMethod.PDB_BASED:
                components_info = self._generate_pdb_based_index(
                    input_params, topology_file
                )
            elif input_params.method == IndexGenerationMethod.TOPOLOGY_BASED:
                components_info = self._generate_topology_based_index(
                    input_params, topology_file
                )
            elif input_params.method == IndexGenerationMethod.SEQUENCE_BASED:
                components_info = self._generate_sequence_based_index(
                    input_params, topology_file
                )
            elif input_params.method == IndexGenerationMethod.CUSTOM:
                components_info = self._generate_custom_index(
                    input_params, topology_file
                )
            else:
                raise InvalidInputError(
                    'method',
                    input_params.method,
                    f'Unsupported method'
                )

            if self.is_cancelled():
                return IndexGenerationResult(
                    success=False,
                    error_message="Operation cancelled by user"
                )

            # Step 4: Validate output
            self._report_progress(0.9, "Validating generated index...")
            self._validate_output(input_params.output_file)

            # Calculate statistics
            end_time = datetime.now()
            processing_time = (end_time - start_time).total_seconds()
            total_atoms = sum(c.atom_count for c in components_info)

            # Build result
            result = IndexGenerationResult(
                success=True,
                output_file=str(Path(input_params.output_file).resolve()),
                components=components_info,
                processing_stats={
                    'n_components': len(components_info),
                    'total_atoms': total_atoms,
                    'processing_time_sec': processing_time
                },
                metadata={
                    'timestamp': start_time.isoformat(),
                    'method': input_params.method.value,
                    'topology': str(Path(input_params.topology).resolve()),
                    'gmx_executable': input_params.gmx_executable
                },
                temporary_files=self.side_effects.files_created.copy()
            )

            self._report_progress(1.0, "Index generation complete!")
            logger.info(f"Successfully generated index with {len(components_info)} components")

            return result

        except InvalidInputError as e:
            logger.error(f"Input validation failed: {e}")
            return IndexGenerationResult(
                success=False,
                error_message=str(e)
            )

        except IndexGenerationError as e:
            logger.error(f"Index generation failed: {e}")
            return IndexGenerationResult(
                success=False,
                error_message=str(e)
            )

        except Exception as e:
            logger.exception("Unexpected error during index generation")
            return IndexGenerationResult(
                success=False,
                error_message=f"Unexpected error: {str(e)}"
            )

    def _prepare_topology(self, input_params: IndexGenerationInput) -> Path:
        """
        Prepare topology file for index generation.

        Returns:
            Path to topology file to use (standardized PDB if available)
        """
        topology_path = Path(input_params.topology)

        # If standardized PDB provided, use it
        if input_params.standardized_pdb:
            std_pdb_path = Path(input_params.standardized_pdb)
            if std_pdb_path.exists():
                logger.info(f"Using pre-standardized PDB: {std_pdb_path.name}")
                return std_pdb_path

        # If auto-standardize and PDB-based method
        if (input_params.auto_standardize and
            input_params.method == IndexGenerationMethod.PDB_BASED):

            # If already PDB, use it directly
            if topology_path.suffix.lower() == '.pdb':
                logger.info(f"Using PDB topology: {topology_path.name}")
                return topology_path

            # Otherwise, convert to PDB
            logger.info(f"Converting {topology_path.suffix} to PDB...")
            pdb_path = self._convert_to_pdb(
                topology_path,
                input_params.gmx_executable
            )
            return pdb_path

        # Use original topology
        return topology_path

    def _convert_to_pdb(self, topology_file: Path, gmx_executable: str) -> Path:
        """Convert TPR/GRO to PDB using gmx editconf."""
        output_pdb = topology_file.parent / f"{topology_file.stem}_converted.pdb"

        cmd = [
            gmx_executable, 'editconf',
            '-f', str(topology_file),
            '-o', str(output_pdb)
        ]

        try:
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60
            )

            self.side_effects.track_command(cmd, process.returncode, process.stderr)

            if process.returncode != 0:
                raise GROMACSCommandError(
                    ' '.join(cmd),
                    process.stderr,
                    "Check that GROMACS is properly installed and topology file is valid"
                )

            if not output_pdb.exists():
                raise TopologyFileError(
                    str(topology_file),
                    "PDB conversion succeeded but file was not created"
                )

            self.side_effects.track_file_creation(str(output_pdb))
            logger.info(f"Converted to PDB: {output_pdb.name}")

            return output_pdb

        except subprocess.TimeoutExpired:
            raise GROMACSCommandError(
                ' '.join(cmd),
                "Command timed out after 60 seconds",
                "Try with a smaller topology file or increase timeout"
            )

    def _generate_pdb_based_index(
        self,
        input_params: IndexGenerationInput,
        topology_file: Path
    ) -> List[ComponentIndexInfo]:
        """Generate index using PDB chain IDs."""
        output_path = Path(input_params.output_file)

        # Build gmx make_ndx commands
        commands = []
        group_id = 18  # Start after default groups

        for comp in input_params.components:
            chain_str = ' '.join(comp.chains)
            commands.append(f"chain {chain_str}")
            commands.append("")  # Empty line to save selection
            commands.append(f"name {group_id} {comp.name}")
            group_id += 1

        commands.append("q")
        command_string = '\n'.join(commands) + '\n'

        # Execute gmx make_ndx
        cmd = [
            input_params.gmx_executable, 'make_ndx',
            '-f', str(topology_file),
            '-o', str(output_path)
        ]

        try:
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate(input=command_string, timeout=30)

            self.side_effects.track_command(cmd, process.returncode, stderr)

            if process.returncode != 0:
                raise GROMACSCommandError(
                    ' '.join(cmd),
                    stderr,
                    "Check chain IDs match those in topology file"
                )

            self.side_effects.track_file_creation(str(output_path))

            # Parse output to get component info
            components_info = self._parse_index_file(output_path, topology_file)

            return components_info

        except subprocess.TimeoutExpired:
            raise GROMACSCommandError(
                ' '.join(cmd),
                "Command timed out",
                "Index generation taking too long"
            )

    def _generate_topology_based_index(
        self,
        input_params: IndexGenerationInput,
        topology_file: Path
    ) -> List[ComponentIndexInfo]:
        """Generate index using topology group IDs."""
        output_path = Path(input_params.output_file)

        # Build gmx make_ndx commands
        commands = []
        name_group_id = 18  # Group ID for naming

        for comp in input_params.components:
            # Select existing groups and combine
            group_str = ' | '.join(str(gid) for gid in comp.group_ids)
            commands.append(group_str)
            commands.append("")  # Empty line to save selection
            commands.append(f"name {name_group_id} {comp.name}")
            name_group_id += 1

        commands.append("q")
        command_string = '\n'.join(commands) + '\n'

        # Execute gmx make_ndx
        cmd = [
            input_params.gmx_executable, 'make_ndx',
            '-f', str(topology_file),
            '-o', str(output_path)
        ]

        try:
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate(input=command_string, timeout=30)

            self.side_effects.track_command(cmd, process.returncode, stderr)

            if process.returncode != 0:
                raise GROMACSCommandError(
                    ' '.join(cmd),
                    stderr,
                    "Check group IDs match those in topology file"
                )

            self.side_effects.track_file_creation(str(output_path))

            # Parse output
            components_info = self._parse_index_file(output_path, topology_file)

            return components_info

        except subprocess.TimeoutExpired:
            raise GROMACSCommandError(
                ' '.join(cmd),
                "Command timed out"
            )

    def _generate_sequence_based_index(
        self,
        input_params: IndexGenerationInput,
        topology_file: Path
    ) -> List[ComponentIndexInfo]:
        """Generate index using sequence matching (for CDR regions)."""
        output_path = Path(input_params.output_file)

        try:
            u = mda.Universe(str(topology_file))

            components_info = []

            with open(output_path, 'w') as f:
                for comp in input_params.components:
                    # Find sequence in topology
                    atoms = self._find_sequence_in_topology(
                        u, comp.sequence, comp.chains[0] if comp.chains else None
                    )

                    if atoms is None:
                        logger.warning(f"Sequence not found for {comp.name}")
                        continue

                    # Filter to CA atoms if requested
                    if comp.ca_only:
                        atoms = atoms.select_atoms("name CA")

                    # Write group to file
                    f.write(f"[ {comp.name} ]\n")
                    atom_indices = atoms.indices + 1  # 1-indexed
                    for i in range(0, len(atom_indices), 10):
                        line = ' '.join(str(idx) for idx in atom_indices[i:i+10])
                        f.write(f"{line}\n")
                    f.write("\n")

                    components_info.append(ComponentIndexInfo(
                        name=comp.name,
                        group_id=len(components_info),
                        atom_count=len(atoms),
                        residue_count=len(atoms.residues),
                        chains=comp.chains
                    ))

            self.side_effects.track_file_creation(str(output_path))

            return components_info

        except Exception as e:
            raise IndexGenerationError(f"Sequence-based index generation failed: {e}")

    def _generate_custom_index(
        self,
        input_params: IndexGenerationInput,
        topology_file: Path
    ) -> List[ComponentIndexInfo]:
        """Generate index using custom MDAnalysis selections."""
        output_path = Path(input_params.output_file)

        try:
            u = mda.Universe(str(topology_file))

            components_info = []

            with open(output_path, 'w') as f:
                for comp in input_params.components:
                    # Use custom selection
                    atoms = u.select_atoms(comp.selection)

                    if len(atoms) == 0:
                        logger.warning(f"No atoms selected for {comp.name}")
                        continue

                    # Write group to file
                    f.write(f"[ {comp.name} ]\n")
                    atom_indices = atoms.indices + 1  # 1-indexed
                    for i in range(0, len(atom_indices), 10):
                        line = ' '.join(str(idx) for idx in atom_indices[i:i+10])
                        f.write(f"{line}\n")
                    f.write("\n")

                    components_info.append(ComponentIndexInfo(
                        name=comp.name,
                        group_id=len(components_info),
                        atom_count=len(atoms),
                        residue_count=len(atoms.residues) if atoms.residues else None
                    ))

            self.side_effects.track_file_creation(str(output_path))

            return components_info

        except Exception as e:
            raise IndexGenerationError(f"Custom index generation failed: {e}")

    def _find_sequence_in_topology(
        self,
        universe: mda.Universe,
        sequence: str,
        chain_id: Optional[str]
    ) -> Optional[mda.AtomGroup]:
        """Find atoms corresponding to sequence in topology."""
        if chain_id:
            chain_atoms = universe.select_atoms(f"chainID {chain_id} and name CA")
        else:
            chain_atoms = universe.select_atoms("name CA")

        # Get chain sequence
        chain_sequence = ''.join([
            mda.lib.util.convert_aa_code(r.resname)
            for r in chain_atoms.residues
        ])

        # Find sequence position
        pos = chain_sequence.find(sequence)
        if pos == -1:
            return None

        # Get residues
        start_idx = pos
        end_idx = pos + len(sequence)
        matched_residues = chain_atoms.residues[start_idx:end_idx]

        # Get all atoms in these residues
        resid_start = matched_residues[0].resid
        resid_end = matched_residues[-1].resid

        if chain_id:
            selection = f"chainID {chain_id} and resid {resid_start}-{resid_end}"
        else:
            selection = f"resid {resid_start}-{resid_end}"

        return universe.select_atoms(selection)

    def _parse_index_file(
        self,
        index_file: Path,
        topology_file: Path
    ) -> List[ComponentIndexInfo]:
        """Parse generated index file to extract component information."""
        components_info = []

        try:
            with open(index_file, 'r') as f:
                lines = f.readlines()

            current_group = None
            group_id = 0

            for line in lines:
                line = line.strip()

                # Group header: [ group_name ]
                if line.startswith('[') and line.endswith(']'):
                    group_name = line[1:-1].strip()

                    # Count atoms in this group
                    atom_indices = []
                    idx = lines.index(line + '\n') + 1
                    while idx < len(lines):
                        next_line = lines[idx].strip()
                        if next_line.startswith('['):
                            break
                        if next_line:
                            atom_indices.extend([int(x) for x in next_line.split()])
                        idx += 1

                    components_info.append(ComponentIndexInfo(
                        name=group_name,
                        group_id=group_id,
                        atom_count=len(atom_indices)
                    ))

                    group_id += 1

            return components_info

        except Exception as e:
            logger.warning(f"Failed to parse index file: {e}")
            return []

    def _validate_output(self, output_file: str):
        """Validate generated index file."""
        output_path = Path(output_file)

        if not output_path.exists():
            raise IndexGenerationError("Index file was not created")

        if output_path.stat().st_size == 0:
            raise IndexGenerationError("Index file is empty")

        # Check file format
        with open(output_path, 'r') as f:
            content = f.read()
            if '[' not in content or ']' not in content:
                raise IndexGenerationError("Index file has invalid format")
