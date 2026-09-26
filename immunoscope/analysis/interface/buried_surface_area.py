"""
Buried surface area (BSA) analyzer.

Following the platform layering rules, this module only owns the calculation:
- compute per-frame SASA for component A, component B, and complex AB
- derive buried surface area and interface ratio
- return standard result objects without pipeline/CLI responsibilities
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Optional

import MDAnalysis as mda
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import freesasa

    FREESASA_AVAILABLE = True
except ImportError:  # pragma: no cover - runtime dependency
    FREESASA_AVAILABLE = False


@dataclass(frozen=True)
class BuriedSurfaceAreaResult:
    """Standard BSA time-series result."""

    times: np.ndarray
    sasa_a: np.ndarray
    sasa_b: np.ndarray
    sasa_ab: np.ndarray
    buried_surface_area: np.ndarray
    interface_ratio: np.ndarray
    selection_a: str
    selection_b: str
    probe_radius: float
    stride: int
    time_unit: str

    def to_timeseries_frame(self) -> pd.DataFrame:
        """Export the per-frame time-series table."""

        return pd.DataFrame(
            {
                f"time_{self.time_unit}": self.times,
                "buried_surface_area": self.buried_surface_area,
                "interface_ratio": self.interface_ratio,
                "sasa_a": self.sasa_a,
                "sasa_b": self.sasa_b,
                "sasa_ab": self.sasa_ab,
            }
        )

    def to_summary(self) -> dict[str, object]:
        """Export summary statistics."""

        def _stats(values: np.ndarray) -> dict[str, float]:
            if values.size == 0:
                return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0}
            return {
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }

        return {
            "selection_a": self.selection_a,
            "selection_b": self.selection_b,
            "probe_radius": float(self.probe_radius),
            "stride": int(self.stride),
            "time_unit": self.time_unit,
            "n_frames": int(self.times.size),
            "time_start": float(self.times[0]) if self.times.size else 0.0,
            "time_end": float(self.times[-1]) if self.times.size else 0.0,
            "buried_surface_area": _stats(self.buried_surface_area),
            "interface_ratio": _stats(self.interface_ratio),
            "sasa_a": _stats(self.sasa_a),
            "sasa_b": _stats(self.sasa_b),
            "sasa_ab": _stats(self.sasa_ab),
        }


class BuriedSurfaceAreaAnalyzer:
    """Calculate buried surface area time series between two components."""

    def __init__(self, topology: str, trajectory: str):
        if not FREESASA_AVAILABLE:
            raise ValueError("FreeSASA is not available. Install freesasa first.")

        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    def calculate(
        self,
        selection_a: str,
        selection_b: str,
        probe_radius: float = 1.4,
        stride: int = 1,
        time_unit: str = "ps",
    ) -> BuriedSurfaceAreaResult:
        """Calculate BSA and interface ratio."""

        atoms_a = self.universe.select_atoms(selection_a)
        atoms_b = self.universe.select_atoms(selection_b)
        atoms_ab = self.universe.select_atoms(f"({selection_a}) or ({selection_b})")
        if len(atoms_a) == 0:
            raise ValueError(f"selection_a selected no atoms: {selection_a}")
        if len(atoms_b) == 0:
            raise ValueError(f"selection_b selected no atoms: {selection_b}")

        indices = range(0, len(self.universe.trajectory), stride)
        n_frames = len(indices)
        times = np.zeros(n_frames, dtype=float)
        sasa_a = np.zeros(n_frames, dtype=float)
        sasa_b = np.zeros(n_frames, dtype=float)
        sasa_ab = np.zeros(n_frames, dtype=float)

        parameters = freesasa.Parameters()
        parameters.setProbeRadius(probe_radius)

        logger.info(
            "Starting BSA calculation: frames=%s stride=%s selection_a=%s selection_b=%s",
            n_frames,
            stride,
            selection_a,
            selection_b,
        )

        for i, frame_index in enumerate(indices):
            ts = self.universe.trajectory[frame_index]
            times[i] = ts.time / 1000.0 if time_unit == "ns" else ts.time
            sasa_a[i] = self._calculate_group_sasa(atoms_a, parameters)
            sasa_b[i] = self._calculate_group_sasa(atoms_b, parameters)
            sasa_ab[i] = self._calculate_group_sasa(atoms_ab, parameters)

        buried_surface_area = (sasa_a + sasa_b - sasa_ab) / 2.0
        denominator = sasa_a + sasa_b
        interface_ratio = np.divide(
            2.0 * buried_surface_area,
            denominator,
            out=np.zeros_like(buried_surface_area),
            where=denominator > 0,
        )

        return BuriedSurfaceAreaResult(
            times=times,
            sasa_a=sasa_a,
            sasa_b=sasa_b,
            sasa_ab=sasa_ab,
            buried_surface_area=buried_surface_area,
            interface_ratio=interface_ratio,
            selection_a=selection_a,
            selection_b=selection_b,
            probe_radius=probe_radius,
            stride=stride,
            time_unit=time_unit,
        )

    @staticmethod
    def _calculate_group_sasa(
        atoms: mda.AtomGroup,
        parameters: "freesasa.Parameters",
    ) -> float:
        """Build an in-memory FreeSASA structure and calculate total surface area."""

        structure = freesasa.Structure()
        for atom in atoms:
            chain_id = BuriedSurfaceAreaAnalyzer._get_chain_id(atom)
            residue_number = str(int(atom.resid))
            structure.addAtom(
                str(atom.name).strip(),
                str(atom.resname).strip(),
                residue_number,
                chain_id,
                float(atom.position[0]),
                float(atom.position[1]),
                float(atom.position[2]),
            )
        result = freesasa.calc(structure, parameters)
        return float(result.totalArea())

    @staticmethod
    def _get_chain_id(atom) -> str:
        """Recover a stable single-character chain ID when possible."""

        for attr in ("chainID", "segid"):
            value: Optional[str] = getattr(atom, attr, None)
            if value:
                value = str(value).strip()
                if value:
                    return value[0]
        return "A"
