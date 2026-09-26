"""Per-residue sidechain dihedral entropy analyzer.

Uses MDAnalysis ``Janin`` analysis to extract χ1 and χ2 angles per residue
per frame, then summarises each residue's rotamer occupancy as:

  * `chi1_entropy` / `chi2_entropy` — Shannon entropy (bits) of the angle
    distribution binned into rotamer wells (60° wide, 6 bins per 360°).
    Higher values mean the sidechain visits multiple rotamers; lower values
    indicate a locked conformation.
  * `rotamer_diversity` — number of rotamer wells visited ≥5% of the time
    (0–6). A pragmatic complement to entropy that maps directly to
    "rotameric flexibility" in design discussions.

Residues lacking χ2 (ALA / CYS / GLY / PRO / SER / THR / VAL) are skipped
by ``Janin`` (warning is emitted by MDAnalysis).
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import math
from typing import Dict, Iterable, List, Optional, Tuple

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.dihedrals import Janin

logger = logging.getLogger(__name__)


_ResKey = Tuple[str, int, str]

# 6 bins of 60° each across [0, 360).
_N_BINS = 6
_BIN_WIDTH = 360.0 / _N_BINS
# Diversity threshold: a bin counts as "visited" only if it holds ≥5% of frames.
_DIVERSITY_THRESHOLD = 0.05


@dataclass(frozen=True)
class ChiDihedralResult:
    """Per-residue chi-entropy result averaged over the trajectory."""

    residue_frame: pd.DataFrame
    summary: dict
    stride: int


class ChiDihedralAnalyzer:
    """Calculate per-residue χ1 / χ2 entropy along an MD trajectory."""

    def __init__(self, topology: str, trajectory: str):
        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    def calculate(
        self,
        selection: str = "protein",
        stride: int = 1,
    ) -> ChiDihedralResult:
        atoms = self.universe.select_atoms(selection)
        if len(atoms) == 0:
            raise ValueError(f"selection selected no atoms: {selection}")

        indices = list(range(0, len(self.universe.trajectory), stride))
        n_frames = len(indices)
        logger.info(
            "Starting chi-dihedral analysis: stride=%s frames=%s",
            stride, n_frames,
        )

        try:
            janin = Janin(atoms).run(
                start=indices[0], stop=indices[-1] + 1, step=stride,
            )
        except ValueError as exc:
            raise ValueError(
                f"Janin analysis failed (possibly altloc / missing atoms): {exc}"
            ) from exc

        angles = janin.results.angles  # shape (n_sampled_frames, n_residues, 2)
        if angles.ndim != 3 or angles.shape[2] != 2:
            raise ValueError(
                f"Unexpected Janin angle shape: {angles.shape}; expected (frames, residues, 2)"
            )

        # Janin removes residues without χ2 from its residue list. Recover the
        # residues that actually contributed angles via its atomgroup.
        residue_group = janin.atomgroup.residues if hasattr(janin, "atomgroup") else None
        if residue_group is None or len(residue_group) != angles.shape[1]:
            # Fall back to the AtomGroup's residues filtered by select_remove default.
            skipped = {"ALA", "CYS", "GLY", "PRO", "SER", "THR", "VAL"}
            residue_group = [
                res for res in atoms.residues
                if str(res.resname).strip().upper() not in skipped
            ]

        records: List[Dict[str, object]] = []
        n_residues_used = min(len(residue_group), angles.shape[1])
        for res_idx in range(n_residues_used):
            residue = residue_group[res_idx]
            chain_id = self._chain_id(residue)
            resid = int(residue.resid)
            resname = str(residue.resname).strip()

            chi1 = angles[:, res_idx, 0]
            chi2 = angles[:, res_idx, 1]
            chi1_entropy = self._entropy(chi1)
            chi2_entropy = self._entropy(chi2)
            rotamer_diversity = self._rotamer_diversity(chi1, chi2)

            records.append(
                {
                    "chain_id": chain_id,
                    "resid": resid,
                    "resname": resname,
                    "chi1_entropy": float(chi1_entropy),
                    "chi2_entropy": float(chi2_entropy),
                    "rotamer_diversity": int(rotamer_diversity),
                    "n_frames_observed": int(len(chi1)),
                }
            )

        residue_frame = pd.DataFrame.from_records(records).sort_values(
            ["chain_id", "resid"]
        ).reset_index(drop=True)

        summary = self._build_summary(residue_frame, selection, stride, n_frames)

        return ChiDihedralResult(
            residue_frame=residue_frame,
            summary=summary,
            stride=stride,
        )

    @staticmethod
    def _entropy(angles: np.ndarray) -> float:
        """Shannon entropy (bits) of a rotamer-binned angle distribution."""
        if len(angles) == 0:
            return float("nan")
        wrapped = np.mod(angles, 360.0)
        bins = np.floor(wrapped / _BIN_WIDTH).astype(int)
        bins = np.clip(bins, 0, _N_BINS - 1)
        counts = np.bincount(bins, minlength=_N_BINS).astype(float)
        total = counts.sum()
        if total <= 0:
            return float("nan")
        probs = counts / total
        nonzero = probs[probs > 0.0]
        entropy = float(-(nonzero * np.log2(nonzero)).sum())
        return abs(entropy)  # avoid -0.0

    @staticmethod
    def _rotamer_diversity(chi1: np.ndarray, chi2: np.ndarray) -> int:
        """Number of (χ1, χ2) joint bins occupied ≥5% of the time."""
        if len(chi1) == 0:
            return 0
        chi1_bins = np.floor(np.mod(chi1, 360.0) / _BIN_WIDTH).astype(int)
        chi2_bins = np.floor(np.mod(chi2, 360.0) / _BIN_WIDTH).astype(int)
        chi1_bins = np.clip(chi1_bins, 0, _N_BINS - 1)
        chi2_bins = np.clip(chi2_bins, 0, _N_BINS - 1)
        joint = chi1_bins * _N_BINS + chi2_bins
        counts = np.bincount(joint, minlength=_N_BINS * _N_BINS).astype(float)
        total = counts.sum()
        if total <= 0:
            return 0
        return int((counts / total >= _DIVERSITY_THRESHOLD).sum())

    @staticmethod
    def _chain_id(residue) -> str:
        for attr in ("segid", "chainID"):
            try:
                value = getattr(residue.atoms[0], attr, None)
            except IndexError:
                value = None
            if value:
                value = str(value).strip()
                if value:
                    return value[0]
        return "A"

    @staticmethod
    def _build_summary(
        residue_frame: pd.DataFrame,
        selection: str,
        stride: int,
        n_frames: int,
    ) -> dict:
        if residue_frame.empty:
            return {
                "selection": selection,
                "stride": int(stride),
                "n_frames": int(n_frames),
                "n_residues": 0,
            }
        return {
            "selection": selection,
            "stride": int(stride),
            "n_frames": int(n_frames),
            "n_residues": int(len(residue_frame)),
            "mean_chi1_entropy": float(residue_frame["chi1_entropy"].mean(skipna=True)),
            "mean_chi2_entropy": float(residue_frame["chi2_entropy"].mean(skipna=True)),
            "mean_rotamer_diversity": float(
                residue_frame["rotamer_diversity"].mean(skipna=True)
            ),
        }
