"""
Post-PBC Quality Validator

This module provides validation tools for assessing trajectory quality after PBC correction.
It checks trajectory integrity, PBC correction quality, and structural stability.

Classes:
    PostPBCValidator: Main validator class for post-PBC quality assessment

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import os
import logging
from typing import Dict, Optional, Tuple
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms


class PostPBCValidator:
    """
    Post-PBC Quality Validator

    Validates trajectory quality after PBC correction by checking:
    1. Trajectory integrity (frame count, file size, corrupted frames)
    2. PBC correction quality (center-of-mass drift, protein integrity)
    3. Structural stability (radius of gyration variation, frame-to-frame jumps)

    Parameters
    ----------
    expected_n_frames : int, optional
        Expected number of frames in the trajectory
    max_com_drift : float, default=1.0
        Maximum allowed center-of-mass drift in nm
    max_rg_std_ratio : float, default=0.15
        Maximum allowed ratio of Rg standard deviation to mean
    max_frame_jump_rmsd : float, default=0.5
        Maximum allowed frame-to-frame RMSD jump in nm
    min_file_size_mb : float, default=1.0
        Minimum expected file size in MB

    Examples
    --------
    >>> validator = PostPBCValidator(max_com_drift=1.0, max_rg_std_ratio=0.15)
    >>> results = validator.comprehensive_validation("md_pbc.xtc", "md.tpr")
    >>> print(f"Overall status: {results['overall_status']}")
    >>> print(f"Quality grade: {results['quality_grade']}")
    """

    def __init__(self,
                 expected_n_frames: Optional[int] = None,
                 max_com_drift: float = 1.0,
                 max_rg_std_ratio: float = 0.15,
                 max_frame_jump_rmsd: float = 0.5,
                 min_file_size_mb: float = 1.0):
        """Initialize PostPBCValidator with quality thresholds"""
        self.expected_n_frames = expected_n_frames
        self.max_com_drift = max_com_drift
        self.max_rg_std_ratio = max_rg_std_ratio
        self.max_frame_jump_rmsd = max_frame_jump_rmsd
        self.min_file_size_mb = min_file_size_mb

        self.logger = logging.getLogger(__name__)

    def validate_trajectory_integrity(self,
                                     trajectory: str,
                                     topology: str) -> Dict:
        """
        Validate trajectory file integrity

        Checks:
        - File existence and size
        - Frame count matches expectation
        - Ability to read all frames (detect corruption)

        Parameters
        ----------
        trajectory : str
            Path to trajectory file
        topology : str
            Path to topology file

        Returns
        -------
        dict
            Validation results with keys:
            - status: 'pass' or 'fail'
            - n_frames: number of frames read
            - file_size_mb: trajectory file size in MB
            - issues: list of detected issues
        """
        results = {
            'status': 'pass',
            'n_frames': 0,
            'file_size_mb': 0.0,
            'issues': []
        }

        # Check file existence
        if not os.path.exists(trajectory):
            results['status'] = 'fail'
            results['issues'].append(f"Trajectory file not found: {trajectory}")
            return results

        # Check file size
        file_size_mb = os.path.getsize(trajectory) / (1024 * 1024)
        results['file_size_mb'] = file_size_mb

        if file_size_mb < self.min_file_size_mb:
            results['status'] = 'fail'
            results['issues'].append(
                f"File size {file_size_mb:.2f} MB < minimum {self.min_file_size_mb} MB"
            )

        # Load trajectory and count frames
        try:
            u = mda.Universe(topology, trajectory)
            n_frames = len(u.trajectory)
            results['n_frames'] = n_frames

            # Check expected frame count
            if self.expected_n_frames is not None:
                if n_frames != self.expected_n_frames:
                    results['issues'].append(
                        f"Frame count {n_frames} != expected {self.expected_n_frames}"
                    )
                    results['status'] = 'fail'

            # Try to read all frames to detect corruption
            corrupted_frames = []
            for ts in u.trajectory:
                if ts.positions is None or len(ts.positions) == 0:
                    corrupted_frames.append(ts.frame)

            if corrupted_frames:
                results['status'] = 'fail'
                results['issues'].append(
                    f"Detected {len(corrupted_frames)} corrupted frames: {corrupted_frames[:10]}"
                )

        except Exception as e:
            results['status'] = 'fail'
            results['issues'].append(f"Failed to load trajectory: {str(e)}")

        return results

    def validate_pbc_correction_quality(self,
                                       trajectory: str,
                                       topology: str,
                                       selection: str = "protein") -> Dict:
        """
        Validate PBC correction quality

        Checks:
        - Center-of-mass drift (should be small after centering)
        - Protein integrity (atoms within reasonable coordinate range)
        - No molecular breakage across periodic boundaries

        Parameters
        ----------
        trajectory : str
            Path to PBC-corrected trajectory file
        topology : str
            Path to topology file
        selection : str, default="protein"
            MDAnalysis selection string for analysis

        Returns
        -------
        dict
            Validation results with keys:
            - status: 'pass' or 'fail'
            - com_drift_nm: center-of-mass drift in nm
            - coordinate_range_nm: coordinate range (min, max) in nm
            - issues: list of detected issues
        """
        results = {
            'status': 'pass',
            'com_drift_nm': 0.0,
            'coordinate_range_nm': (0.0, 0.0),
            'issues': []
        }

        try:
            u = mda.Universe(topology, trajectory)
            selected_atoms = u.select_atoms(selection)

            if len(selected_atoms) == 0:
                results['status'] = 'fail'
                results['issues'].append(f"No atoms selected with: {selection}")
                return results

            # Calculate center-of-mass trajectory
            com_positions = []
            for ts in u.trajectory:
                com = selected_atoms.center_of_mass()
                com_positions.append(com)

            com_positions = np.array(com_positions)

            # Calculate COM drift (max distance from mean position)
            com_mean = np.mean(com_positions, axis=0)
            com_distances = np.linalg.norm(com_positions - com_mean, axis=1)
            com_drift_angstrom = np.max(com_distances)
            com_drift_nm = com_drift_angstrom / 10.0

            results['com_drift_nm'] = com_drift_nm

            if com_drift_nm > self.max_com_drift:
                results['status'] = 'fail'
                results['issues'].append(
                    f"COM drift {com_drift_nm:.3f} nm > threshold {self.max_com_drift} nm"
                )

            # Check coordinate range
            all_positions = []
            for ts in u.trajectory:
                all_positions.append(selected_atoms.positions)

            all_positions = np.vstack(all_positions)
            coord_min = np.min(all_positions) / 10.0  # Convert to nm
            coord_max = np.max(all_positions) / 10.0

            results['coordinate_range_nm'] = (coord_min, coord_max)

            # Check if coordinates are within reasonable range
            # After PBC correction, should be within a box (roughly -5 to 5 nm from center)
            if coord_min < -50.0 or coord_max > 50.0:
                results['issues'].append(
                    f"Unusual coordinate range: [{coord_min:.1f}, {coord_max:.1f}] nm"
                )

        except Exception as e:
            results['status'] = 'fail'
            results['issues'].append(f"PBC quality check failed: {str(e)}")

        return results

    def validate_structural_stability(self,
                                     trajectory: str,
                                     topology: str,
                                     selection: str = "protein",
                                     stride: int = 1) -> Dict:
        """
        Validate structural stability throughout trajectory

        Checks:
        - Radius of gyration (Rg) variation (std/mean ratio)
        - Frame-to-frame RMSD jumps (detect sudden changes)

        Parameters
        ----------
        trajectory : str
            Path to trajectory file
        topology : str
            Path to topology file
        selection : str, default="protein"
            MDAnalysis selection string for analysis
        stride : int, default=1
            Stride for sampling frames (use >1 for faster validation)

        Returns
        -------
        dict
            Validation results with keys:
            - status: 'pass' or 'fail'
            - rg_mean_nm: mean radius of gyration in nm
            - rg_std_nm: Rg standard deviation in nm
            - rg_std_ratio: std/mean ratio
            - max_frame_jump_rmsd_nm: maximum frame-to-frame RMSD in nm
            - issues: list of detected issues
        """
        results = {
            'status': 'pass',
            'rg_mean_nm': 0.0,
            'rg_std_nm': 0.0,
            'rg_std_ratio': 0.0,
            'max_frame_jump_rmsd_nm': 0.0,
            'issues': []
        }

        try:
            u = mda.Universe(topology, trajectory)
            selected_atoms = u.select_atoms(selection)

            if len(selected_atoms) == 0:
                results['status'] = 'fail'
                results['issues'].append(f"No atoms selected with: {selection}")
                return results

            # Calculate radius of gyration for all frames
            rg_values = []
            for ts in u.trajectory[::stride]:
                rg = selected_atoms.radius_of_gyration()
                rg_values.append(rg)

            rg_values = np.array(rg_values) / 10.0  # Convert to nm

            rg_mean = np.mean(rg_values)
            rg_std = np.std(rg_values)
            rg_std_ratio = rg_std / rg_mean if rg_mean > 0 else 0.0

            results['rg_mean_nm'] = rg_mean
            results['rg_std_nm'] = rg_std
            results['rg_std_ratio'] = rg_std_ratio

            if rg_std_ratio > self.max_rg_std_ratio:
                results['status'] = 'fail'
                results['issues'].append(
                    f"Rg std/mean ratio {rg_std_ratio:.3f} > threshold {self.max_rg_std_ratio}"
                )

            # Calculate frame-to-frame RMSD to detect jumps
            # Use a subset of frames for efficiency
            frame_indices = list(range(0, len(u.trajectory), max(stride, 10)))
            if len(frame_indices) > 100:
                frame_indices = frame_indices[:100]  # Limit to 100 frames

            frame_rmsd_jumps = []
            for i in range(len(frame_indices) - 1):
                u.trajectory[frame_indices[i]]
                pos1 = selected_atoms.positions.copy()

                u.trajectory[frame_indices[i + 1]]
                pos2 = selected_atoms.positions.copy()

                # Calculate RMSD between consecutive frames
                rmsd_jump = rms.rmsd(pos1, pos2, superposition=True)
                frame_rmsd_jumps.append(rmsd_jump)

            if frame_rmsd_jumps:
                max_jump = np.max(frame_rmsd_jumps) / 10.0  # Convert to nm
                results['max_frame_jump_rmsd_nm'] = max_jump

                if max_jump > self.max_frame_jump_rmsd:
                    results['status'] = 'fail'
                    results['issues'].append(
                        f"Max frame-to-frame RMSD {max_jump:.3f} nm > threshold {self.max_frame_jump_rmsd} nm"
                    )

        except Exception as e:
            results['status'] = 'fail'
            results['issues'].append(f"Structural stability check failed: {str(e)}")

        return results

    def comprehensive_validation(self,
                                trajectory: str,
                                topology: str,
                                selection: str = "protein",
                                stride: int = 1) -> Dict:
        """
        Run comprehensive validation including all checks

        Parameters
        ----------
        trajectory : str
            Path to PBC-corrected trajectory file
        topology : str
            Path to topology file
        selection : str, default="protein"
            MDAnalysis selection string for analysis
        stride : int, default=1
            Stride for sampling frames in stability checks

        Returns
        -------
        dict
            Comprehensive validation results with keys:
            - overall_status: 'pass' or 'fail'
            - overall_grade: quality grade ('A', 'B', 'C', 'D')
            - integrity: results from validate_trajectory_integrity()
            - pbc_quality: results from validate_pbc_correction_quality()
            - stability: results from validate_structural_stability()
            - all_issues: combined list of all issues
        """
        self.logger.info(f"Running comprehensive post-PBC validation on {trajectory}")

        # Run all validation checks
        integrity = self.validate_trajectory_integrity(trajectory, topology)
        pbc_quality = self.validate_pbc_correction_quality(trajectory, topology, selection)
        stability = self.validate_structural_stability(trajectory, topology, selection, stride)

        # Combine results
        all_issues = (
            integrity.get('issues', []) +
            pbc_quality.get('issues', []) +
            stability.get('issues', [])
        )

        # Determine overall status
        overall_status = 'pass'
        if (integrity['status'] == 'fail' or
            pbc_quality['status'] == 'fail' or
            stability['status'] == 'fail'):
            overall_status = 'fail'

        # Assign quality grade based on metrics
        grade = self._assign_quality_grade(integrity, pbc_quality, stability)

        results = {
            'overall_status': overall_status,
            'overall_grade': grade,
            'integrity': integrity,
            'pbc_quality': pbc_quality,
            'stability': stability,
            'all_issues': all_issues
        }

        self.logger.info(f"Validation complete: Status={overall_status}, Grade={grade}")

        return results

    def _assign_quality_grade(self,
                             integrity: Dict,
                             pbc_quality: Dict,
                             stability: Dict) -> str:
        """
        Assign quality grade (A/B/C/D) based on validation results

        Grading criteria:
        - Grade A: All checks pass with excellent metrics
        - Grade B: All checks pass with good metrics
        - Grade C: All checks pass with acceptable metrics
        - Grade D: One or more checks fail
        """
        # Grade D if any check fails
        if (integrity['status'] == 'fail' or
            pbc_quality['status'] == 'fail' or
            stability['status'] == 'fail'):
            return 'D'

        # Score based on metrics (0-100)
        score = 100.0

        # COM drift penalty (0-20 points)
        com_drift = pbc_quality.get('com_drift_nm', 0.0)
        if com_drift > self.max_com_drift * 0.5:
            score -= 20 * (com_drift / self.max_com_drift)

        # Rg stability penalty (0-30 points)
        rg_std_ratio = stability.get('rg_std_ratio', 0.0)
        if rg_std_ratio > self.max_rg_std_ratio * 0.5:
            score -= 30 * (rg_std_ratio / self.max_rg_std_ratio)

        # Frame jump penalty (0-20 points)
        max_jump = stability.get('max_frame_jump_rmsd_nm', 0.0)
        if max_jump > self.max_frame_jump_rmsd * 0.5:
            score -= 20 * (max_jump / self.max_frame_jump_rmsd)

        # Assign grade
        if score >= 90:
            return 'A'
        elif score >= 75:
            return 'B'
        elif score >= 60:
            return 'C'
        else:
            return 'D'
