"""
RMSD Calculator - Refactored version with standardized input/output.

This module provides RMSD calculation with clear inputs, outputs, and error handling.
"""

import numpy as np
import pandas as pd
import subprocess
from dataclasses import dataclass
from typing import Optional, Tuple
from pathlib import Path
import logging

from ...core.context import ProcessingResult
from ...core.exceptions import InputValidationError, ProcessingError

# Optional MDAnalysis import
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    mda = None
    rms = None
    align = None

logger = logging.getLogger(__name__)


@dataclass
class RMSDInput:
    """
    Standardized input for RMSD calculation.

    Attributes:
        topology: Path to topology file (.tpr, .gro, .pdb)
        trajectory: Path to trajectory file (.xtc, .trr)
        selection: Atom selection string (MDAnalysis syntax)
        reference_frame: Reference frame index (default: 0)
        output_file: Optional output file path
        method: Calculation method ('mdanalysis' or 'gromacs')
        gmx_executable: GROMACS executable command
    """
    topology: str
    trajectory: str
    selection: str = "protein and name CA"
    reference_frame: int = 0
    output_file: Optional[str] = None
    method: str = "mdanalysis"
    gmx_executable: str = "gmx"

    def validate(self):
        """
        Validate input parameters.

        Raises:
            InputValidationError: If validation fails
        """
        # Check topology file
        if not Path(self.topology).exists():
            raise InputValidationError(
                param_name="topology",
                value=self.topology,
                reason="file does not exist"
            )

        # Check trajectory file
        if not Path(self.trajectory).exists():
            raise InputValidationError(
                param_name="trajectory",
                value=self.trajectory,
                reason="file does not exist"
            )

        # Check method
        if self.method not in ["mdanalysis", "gromacs"]:
            raise InputValidationError(
                param_name="method",
                value=self.method,
                reason="must be 'mdanalysis' or 'gromacs'"
            )

        # Check MDAnalysis availability if method is mdanalysis
        if self.method == "mdanalysis" and not HAS_MDANALYSIS:
            raise InputValidationError(
                param_name="method",
                value="mdanalysis",
                reason="MDAnalysis is not installed"
            )

        # Check reference frame
        if self.reference_frame < 0:
            raise InputValidationError(
                param_name="reference_frame",
                value=self.reference_frame,
                reason="must be >= 0"
            )


@dataclass
class RMSDResult:
    """
    Standardized output from RMSD calculation.

    Attributes:
        success: Whether calculation succeeded
        output_file: Path to output file (if saved)
        times: Time array (ps)
        rmsd_values: RMSD values array (nm or Angstrom)
        mean_rmsd: Mean RMSD value
        std_rmsd: Standard deviation of RMSD
        min_rmsd: Minimum RMSD value
        max_rmsd: Maximum RMSD value
        n_frames: Number of frames processed
        metadata: Additional metadata
        error_message: Error message if failed
    """
    success: bool
    output_file: Optional[str] = None
    times: Optional[np.ndarray] = None
    rmsd_values: Optional[np.ndarray] = None
    mean_rmsd: Optional[float] = None
    std_rmsd: Optional[float] = None
    min_rmsd: Optional[float] = None
    max_rmsd: Optional[float] = None
    n_frames: Optional[int] = None
    metadata: Optional[dict] = None
    error_message: Optional[str] = None

    def to_dict(self) -> dict:
        """Convert to dictionary (without numpy arrays)."""
        return {
            'success': self.success,
            'output_file': self.output_file,
            'mean_rmsd': self.mean_rmsd,
            'std_rmsd': self.std_rmsd,
            'min_rmsd': self.min_rmsd,
            'max_rmsd': self.max_rmsd,
            'n_frames': self.n_frames,
            'metadata': self.metadata,
            'error_message': self.error_message
        }


class RMSDCalculator:
    """
    RMSD calculation module (refactored version).

    This calculator provides RMSD calculation with standardized input/output
    following the six design principles.
    """

    def __init__(self):
        """Initialize RMSDCalculator."""
        pass

    def calculate(self, input_params: RMSDInput) -> RMSDResult:
        """
        Calculate RMSD with standardized input/output.

        This is the main entry point following the design principles:
        1. Clear inputs (RMSDInput)
        2. Clear outputs (RMSDResult)
        3. Input validation
        4. Error handling
        5. No side effects (only reads files, writes to specified output)

        Args:
            input_params: Standardized input parameters

        Returns:
            RMSDResult with calculation results

        Example:
            >>> calculator = RMSDCalculator()
            >>> input_params = RMSDInput(
            ...     topology="md.tpr",
            ...     trajectory="md_processed.xtc",
            ...     selection="protein and name CA"
            ... )
            >>> result = calculator.calculate(input_params)
            >>> print(f"Mean RMSD: {result.mean_rmsd:.3f} nm")
        """
        # Validate inputs
        try:
            input_params.validate()
        except InputValidationError as e:
            logger.error(f"Input validation failed: {e}")
            return RMSDResult(
                success=False,
                error_message=str(e)
            )

        # Execute calculation
        try:
            if input_params.method == "mdanalysis":
                return self._calculate_mdanalysis(input_params)
            else:
                return self._calculate_gromacs(input_params)

        except Exception as e:
            logger.exception(f"RMSD calculation failed: {e}")
            return RMSDResult(
                success=False,
                error_message=f"Calculation error: {str(e)}"
            )

    def _calculate_mdanalysis(self, params: RMSDInput) -> RMSDResult:
        """Calculate RMSD using MDAnalysis."""
        try:
            # Load universe
            universe = mda.Universe(params.topology, params.trajectory)
            logger.info(f"Loaded {len(universe.trajectory)} frames")

            # Select atoms
            atoms = universe.select_atoms(params.selection)
            if len(atoms) == 0:
                raise ProcessingError(
                    step="atom selection",
                    reason=f"selection '{params.selection}' matched 0 atoms"
                )

            logger.info(f"Selected {len(atoms)} atoms")

            # Set reference frame
            if params.reference_frame >= len(universe.trajectory):
                raise InputValidationError(
                    param_name="reference_frame",
                    value=params.reference_frame,
                    reason=f"exceeds trajectory length ({len(universe.trajectory)})"
                )

            universe.trajectory[params.reference_frame]
            reference_coords = atoms.positions.copy()

            # Calculate RMSD for each frame
            rmsd_values = []
            times = []

            for ts in universe.trajectory:
                # MDAnalysis returns RMSD in Angstrom; convert it to nm here.
                rmsd_val = rms.rmsd(atoms.positions, reference_coords, superposition=True) / 10.0
                rmsd_values.append(rmsd_val)
                times.append(ts.time)

            times = np.array(times)
            rmsd_values = np.array(rmsd_values)

            # Save to file if requested
            if params.output_file:
                self._save_data(times, rmsd_values, params.output_file)

            # Calculate statistics
            mean_rmsd = float(np.mean(rmsd_values))
            std_rmsd = float(np.std(rmsd_values))
            min_rmsd = float(np.min(rmsd_values))
            max_rmsd = float(np.max(rmsd_values))

            logger.info(f"RMSD calculation completed. Mean: {mean_rmsd:.3f} nm, Std: {std_rmsd:.3f} nm")

            return RMSDResult(
                success=True,
                output_file=params.output_file,
                times=times,
                rmsd_values=rmsd_values,
                mean_rmsd=mean_rmsd,
                std_rmsd=std_rmsd,
                min_rmsd=min_rmsd,
                max_rmsd=max_rmsd,
                n_frames=len(times),
                metadata={
                    'method': 'mdanalysis',
                    'selection': params.selection,
                    'reference_frame': params.reference_frame,
                    'n_atoms': len(atoms)
                }
            )

        except Exception as e:
            logger.exception(f"MDAnalysis RMSD calculation failed: {e}")
            return RMSDResult(
                success=False,
                error_message=str(e)
            )

    def _calculate_gromacs(self, params: RMSDInput) -> RMSDResult:
        """Calculate RMSD using GROMACS."""
        # TODO: Implement GROMACS method
        raise NotImplementedError("GROMACS method not yet implemented in refactored version")

    def _save_data(self, times: np.ndarray, rmsd_values: np.ndarray, output_file: str):
        """Save RMSD data to file."""
        df = pd.DataFrame({
            'Time (ps)': times,
            'RMSD (nm)': rmsd_values
        })

        # Determine format based on extension
        if output_file.endswith('.csv'):
            df.to_csv(output_file, index=False)
        elif output_file.endswith('.xvg'):
            # GROMACS XVG format
            with open(output_file, 'w') as f:
                f.write("# RMSD calculation\n")
                f.write("@    title \"RMSD\"\n")
                f.write("@    xaxis  label \"Time (ps)\"\n")
                f.write("@    yaxis  label \"RMSD (nm)\"\n")
                f.write("@TYPE xy\n")
                for t, r in zip(times, rmsd_values):
                    f.write(f"{t:12.3f} {r:12.6f}\n")
        else:
            # Default to CSV
            df.to_csv(output_file, index=False)

        logger.info(f"RMSD data saved to {output_file}")
