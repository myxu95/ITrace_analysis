"""
Angle Analysis Module for MD Simulations

This module provides docking angle calculation for TCR-pMHC complexes.

Public API (3 classes only):
-----------------------------
- DockingAngleAnalyzer: Main analyzer class
- DockingAngleInput: Standardized input parameters
- DockingAngleResult: Standardized output results

All internal implementation details are hidden within analyzer.py module.

Key Features:
-------------
- Sequence-based MHC region detection (robust to non-standard PDB numbering)
- ANARCI integration for TCR disulfide identification
- Geometric angle definitions (Crossing + Incident)
- Pure Python implementation (MDAnalysis + NumPy)
- Support for both static structures and MD trajectories
- Conforms to ITrace 6 design principles

Example Usage:
--------------
Basic trajectory analysis:
>>> from md_analysis.analysis.angles import DockingAngleAnalyzer, DockingAngleInput
>>> analyzer = DockingAngleAnalyzer()
>>> result = analyzer.analyze(DockingAngleInput(
...     topology='md.tpr',
...     trajectory='md_pbc.xtc',
...     stride=10,
...     output_dir='./results'
... ))
>>> print(f"Crossing: {result.statistics['crossing_mean']:.2f}°")

With progress tracking:
>>> def progress_callback(progress, message):
...     print(f"[{progress*100:.0f}%] {message}")
>>> analyzer.set_progress_callback(progress_callback)
>>> result = analyzer.analyze(input_params)

Notes:
------
- Refactored 2026-03-18: Consolidated from 9 files to 3 files
- All internal classes are hidden (prefixed with _)
- Only 3 classes are exported in public API
- Old modules archived for reference
"""

# Public API - Only 3 classes exported through __all__.
from .angle_data_structures import DockingAngleInput, DockingAngleResult
from .analyzer import DockingAngleAnalyzer, _GeometryUtils
from .cdr3_geometry import CDR3GeometryAnalyzer, CDR3GeometryResult


class PrincipalAxesCalculator:
    """
    Backward-compatible principal-axis helper.

    The refactored docking-angle module keeps this outside ``__all__`` so the
    compact public API remains three primary classes, while older tests and
    notebooks can still import the numerical utility directly.
    """

    def calculate_inertia_tensor(self, positions, masses):
        """Calculate a mass-weighted 3D inertia tensor."""
        import numpy as np

        coords = np.asarray(positions, dtype=float)
        weights = np.asarray(masses, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError("positions must have shape (N, 3)")
        if weights.shape[0] != coords.shape[0]:
            raise ValueError("masses length must match positions")

        center = np.average(coords, axis=0, weights=weights)
        shifted = coords - center
        x, y, z = shifted.T
        tensor = np.array([
            [np.sum(weights * (y * y + z * z)), -np.sum(weights * x * y), -np.sum(weights * x * z)],
            [-np.sum(weights * y * x), np.sum(weights * (x * x + z * z)), -np.sum(weights * y * z)],
            [-np.sum(weights * z * x), -np.sum(weights * z * y), np.sum(weights * (x * x + y * y))],
        ])
        return tensor

    def calculate_principal_axes(self, positions, masses):
        """Return principal axes and moments sorted in descending order."""
        import numpy as np

        tensor = self.calculate_inertia_tensor(positions, masses)
        moments, vectors = np.linalg.eigh(tensor)
        order = np.argsort(moments)[::-1]
        return vectors[:, order].T, moments[order]


def angle_between_vectors(v1, v2, degrees=True):
    """Backward-compatible wrapper around the refactored geometry utility."""
    return _GeometryUtils.angle_between_vectors(v1, v2, degrees=degrees)


def project_vector_onto_plane(vector, plane_normal):
    """Backward-compatible wrapper around the refactored geometry utility."""
    return _GeometryUtils.project_vector_onto_plane(vector, plane_normal)


def signed_angle_between_vectors(v1, v2, normal, degrees=True):
    """Calculate signed angle from v1 to v2 around a reference normal."""
    import numpy as np

    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    normal_norm = normal / np.linalg.norm(normal)
    unsigned = _GeometryUtils.angle_between_vectors(v1_norm, v2_norm, degrees=degrees)
    sign = np.sign(np.dot(np.cross(v1_norm, v2_norm), normal_norm))
    return float(unsigned * (sign if sign != 0 else 1.0))

__all__ = [
    'DockingAngleAnalyzer',
    'DockingAngleInput',
    'DockingAngleResult',
]

# Version info
__version__ = '4.3.0'
__refactoring_date__ = '2026-03-18'
__specification__ = 'Consolidated module with hidden internals'
