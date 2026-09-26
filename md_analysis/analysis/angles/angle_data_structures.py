"""
Angle Analysis Data Structures

Standard input/output structures for docking angle analysis,
following ITrace module design standards.

Author: ITrace Development Team
Date: 2026-03-18
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, List, Any
from pathlib import Path
import numpy as np


@dataclass
class DockingAngleInput:
    """
    Standardized input for docking angle analysis.

    Follows ITrace module design standard: Clear Inputs

    Parameters
    ----------
    topology : str
        Path to topology file (TPR, PDB, etc.)
    trajectory : str, optional
        Path to trajectory file (XTC, TRR, etc.)
        If None, analyze single frame from topology
    mhc_selection : str, optional
        MDAnalysis selection for MHC chain (e.g., 'chainID A')
        If None and auto_identify_chains=True, will auto-detect
    tcr_alpha_selection : str, optional
        MDAnalysis selection for TCR alpha chain
    tcr_beta_selection : str, optional
        MDAnalysis selection for TCR beta chain
    auto_identify_chains : bool, default=True
        Whether to automatically identify chains using ANARCI
    use_anarci : bool, default=True
        Whether to use ANARCI for chain identification
    stride : int, default=1
        Use every Nth frame for trajectory analysis
    print_each_frame : bool, default=False
        Whether to print per-frame angle results to the terminal
    output_dir : str, optional
        Output directory for results
    """
    # Required parameters
    topology: str
    trajectory: Optional[str] = None

    # Chain selections (optional - supports auto identification)
    mhc_selection: Optional[str] = None
    tcr_alpha_selection: Optional[str] = None
    tcr_beta_selection: Optional[str] = None

    # Auto chain identification
    auto_identify_chains: bool = True
    use_anarci: bool = True

    # Analysis parameters
    stride: int = 1
    print_each_frame: bool = False
    output_dir: Optional[str] = None

    def validate(self) -> None:
        """
        Validate input parameters.

        Raises
        ------
        FileNotFoundError
            If topology or trajectory files do not exist
        ValueError
            If parameters are invalid or inconsistent
        """
        # File existence checks
        if not Path(self.topology).exists():
            raise FileNotFoundError(f"Topology file not found: {self.topology}")

        if self.trajectory and not Path(self.trajectory).exists():
            raise FileNotFoundError(f"Trajectory file not found: {self.trajectory}")

        # Parameter range checks
        if self.stride < 1:
            raise ValueError(f"stride must be >= 1, got {self.stride}")

        # Mode validation
        if not self.auto_identify_chains:
            # Manual mode: require all three selections
            if not (self.mhc_selection and self.tcr_alpha_selection and self.tcr_beta_selection):
                raise ValueError(
                    "When auto_identify_chains=False, must provide all three selections: "
                    "mhc_selection, tcr_alpha_selection, tcr_beta_selection"
                )

    def to_dict(self) -> Dict[str, Any]:
        """Convert input parameters to dictionary."""
        return asdict(self)


@dataclass
class DockingAngleResult:
    """
    Standardized output for docking angle analysis.

    Follows ITrace module design standard: Clear Outputs

    Attributes
    ----------
    success : bool
        Whether the analysis completed successfully
    crossing_angle : float, optional
        Crossing angle for single frame (degrees)
    incident_angle : float, optional
        Incident angle for single frame (degrees)
    times : np.ndarray, optional
        Time points for trajectory analysis (ps)
    crossing_angles : np.ndarray, optional
        Crossing angles over trajectory (degrees)
    incident_angles : np.ndarray, optional
        Incident angles over trajectory (degrees)
    statistics : dict, optional
        Statistical summary with keys:
        - crossing_mean, crossing_std, crossing_min, crossing_max
        - incident_mean, incident_std, incident_min, incident_max
        - n_frames
    output_files : list
        List of generated output file paths
    metadata : dict
        Analysis metadata with keys:
        - module_version: version string
        - timestamp: ISO format timestamp
        - topology, trajectory: input file paths
        - stride: frame stride used
        - chain_identifications: chain identification results (if auto mode)
        - alignment_quality: MHC sequence alignment quality metrics
    error_message : str, optional
        Error message if analysis failed
    """
    # Execution status
    success: bool

    # Single frame results
    crossing_angle: Optional[float] = None
    incident_angle: Optional[float] = None

    # Trajectory results
    times: Optional[np.ndarray] = None
    crossing_angles: Optional[np.ndarray] = None
    incident_angles: Optional[np.ndarray] = None

    # Statistics
    statistics: Optional[Dict[str, Any]] = None

    # Output files (clear side effects)
    output_files: List[str] = field(default_factory=list)

    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)

    # Error information (if failed)
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert result to dictionary for serialization.

        Returns
        -------
        dict
            Dictionary representation with numpy arrays converted to lists
        """
        result = {
            'success': self.success,
            'crossing_angle': self.crossing_angle,
            'incident_angle': self.incident_angle,
            'statistics': self.statistics,
            'output_files': self.output_files,
            'metadata': self.metadata,
            'error_message': self.error_message
        }

        # Convert numpy arrays to lists (JSON compatible)
        if self.times is not None:
            result['times'] = self.times.tolist()
        if self.crossing_angles is not None:
            result['crossing_angles'] = self.crossing_angles.tolist()
        if self.incident_angles is not None:
            result['incident_angles'] = self.incident_angles.tolist()

        return result

    def get_summary(self) -> str:
        """
        Get human-readable summary of results.

        Returns
        -------
        str
            Formatted summary string
        """
        if not self.success:
            return f"Analysis failed: {self.error_message}"

        lines = ["Docking Angle Analysis Results", "=" * 40]

        if self.crossing_angle is not None:
            # Single frame
            lines.extend([
                f"Crossing angle: {self.crossing_angle:.2f}°",
                f"Incident angle: {self.incident_angle:.2f}°"
            ])
        elif self.statistics:
            # Trajectory
            stats = self.statistics

            # Add frames analyzed if available
            if 'n_frames' in stats:
                lines.append(f"Frames analyzed: {stats['n_frames']}")
                lines.append("")

            # Crossing angle statistics
            if 'crossing_mean' in stats:
                lines.append("Crossing angle:")
                lines.append(f"  Mean: {stats['crossing_mean']:.2f}°")
                if 'crossing_std' in stats:
                    lines[-1] += f" ± {stats['crossing_std']:.2f}°"
                if 'crossing_min' in stats and 'crossing_max' in stats:
                    lines.append(f"  Range: [{stats['crossing_min']:.2f}, {stats['crossing_max']:.2f}]°")
                lines.append("")

            # Incident angle statistics
            if 'incident_mean' in stats:
                lines.append("Incident angle:")
                lines.append(f"  Mean: {stats['incident_mean']:.2f}°")
                if 'incident_std' in stats:
                    lines[-1] += f" ± {stats['incident_std']:.2f}°"
                if 'incident_min' in stats and 'incident_max' in stats:
                    lines.append(f"  Range: [{stats['incident_min']:.2f}, {stats['incident_max']:.2f}]°")

        if self.output_files:
            lines.extend(["", "Output files:"])
            for f in self.output_files:
                lines.append(f"  - {f}")

        return "\n".join(lines)
