"""TCR-MHC center-of-mass distance analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import MDAnalysis as mda

    HAS_MDANALYSIS = True
except ImportError:  # pragma: no cover - runtime dependency
    mda = None
    HAS_MDANALYSIS = False


@dataclass
class COMDistanceInput:
    """Center-of-mass distance analysis input."""

    tcr_selection: str = "protein and chainid A B"
    mhc_selection: str = "protein and chainid C D E"
    stride: int = 1
    time_unit: str = "ps"

    def validate(self) -> None:
        if not self.tcr_selection or not self.tcr_selection.strip():
            raise ValueError("tcr_selection cannot be empty")
        if not self.mhc_selection or not self.mhc_selection.strip():
            raise ValueError("mhc_selection cannot be empty")
        if self.stride < 1:
            raise ValueError("stride must be >= 1")
        if not self.time_unit:
            raise ValueError("time_unit cannot be empty")


@dataclass(frozen=True)
class COMDistanceResult:
    """Center-of-mass distance analysis result."""

    times: np.ndarray
    distances: np.ndarray
    mean: float
    std: float
    min: float
    max: float
    tcr_selection: str
    mhc_selection: str
    stride: int
    time_unit: str

    def to_timeseries_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                f"time_{self.time_unit}": self.times,
                "com_distance_angstrom": self.distances,
            }
        )

    def to_summary(self) -> dict[str, float | str | int]:
        return {
            "mean": float(self.mean),
            "std": float(self.std),
            "min": float(self.min),
            "max": float(self.max),
            "n_frames": int(self.times.size),
            "stride": int(self.stride),
            "time_unit": self.time_unit,
            "tcr_selection": self.tcr_selection,
            "mhc_selection": self.mhc_selection,
        }


class COMDistanceCalculator:
    """Calculate per-frame TCR-MHC center-of-mass distance."""

    def __init__(self, topology: str, trajectory: str):
        if not HAS_MDANALYSIS:
            raise ImportError("COMDistanceCalculator requires MDAnalysis. Install it first.")

        if not Path(topology).exists():
            raise FileNotFoundError(f"Topology file not found: {topology}")
        if not Path(trajectory).exists():
            raise FileNotFoundError(f"Trajectory file not found: {trajectory}")

        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    def calculate(self, input_params: COMDistanceInput) -> COMDistanceResult:
        input_params.validate()

        tcr_atoms = self.universe.select_atoms(input_params.tcr_selection)
        mhc_atoms = self.universe.select_atoms(input_params.mhc_selection)
        if len(tcr_atoms) == 0:
            raise ValueError(f"tcr_selection selected no atoms: {input_params.tcr_selection}")
        if len(mhc_atoms) == 0:
            raise ValueError(f"mhc_selection selected no atoms: {input_params.mhc_selection}")

        times: list[float] = []
        distances: list[float] = []

        for ts in self.universe.trajectory[:: input_params.stride]:
            tcr_com = tcr_atoms.center_of_mass()
            mhc_com = mhc_atoms.center_of_mass()
            distance = float(np.linalg.norm(tcr_com - mhc_com))
            frame_time = float(getattr(ts, "time", ts.frame))
            times.append(frame_time)
            distances.append(distance)

        if not distances:
            raise ValueError("No frames are available for calculation")

        distance_array = np.asarray(distances, dtype=float)
        time_array = np.asarray(times, dtype=float)
        return COMDistanceResult(
            times=time_array,
            distances=distance_array,
            mean=float(np.mean(distance_array)),
            std=float(np.std(distance_array)),
            min=float(np.min(distance_array)),
            max=float(np.max(distance_array)),
            tcr_selection=input_params.tcr_selection,
            mhc_selection=input_params.mhc_selection,
            stride=int(input_params.stride),
            time_unit=input_params.time_unit,
        )
