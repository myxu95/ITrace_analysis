"""Per-residue secondary structure (DSSP) analyzer.

Runs MDAnalysis's DSSP implementation across a trajectory and reduces the
per-frame H/E/L codes to:

  * `secondary_structure` — the dominant code over the sampled frames
  * `ss_propensity_helix/sheet/loop` — fraction of frames in each bucket
  * `ss_stability` — fraction of frames matching the dominant code

For mutation design these answer "is the residue inside a rigid SS element
the substitution would disturb, or in a flexible loop where exchanges are
relatively safe?" — a question RMSF magnitude alone cannot answer.

The MDAnalysis DSSP reports three codes: ``H`` (helix), ``E`` (extended /
sheet), and ``-`` (coil / loop). We aggregate them per residue.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import logging
from typing import Dict, List, Optional, Tuple

import MDAnalysis as mda
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    from MDAnalysis.analysis.dssp import DSSP as _MDADSSP

    DSSP_AVAILABLE = True
except ImportError:  # pragma: no cover
    _MDADSSP = None
    DSSP_AVAILABLE = False


_ResKey = Tuple[str, int, str]


@dataclass(frozen=True)
class DSSPResult:
    """Per-residue secondary structure summary."""

    residue_frame: pd.DataFrame
    summary: dict
    stride: int


class DSSPAnalyzer:
    """Per-residue DSSP aggregation over an MD trajectory."""

    _HELIX_CODES = {"H", "G", "I"}    # MDAnalysis DSSP only emits H, but stay permissive
    _SHEET_CODES = {"E", "B"}

    def __init__(self, topology: str, trajectory: str):
        if not DSSP_AVAILABLE:
            raise ValueError(
                "MDAnalysis DSSP is not available. Upgrade MDAnalysis to a version "
                "that ships `MDAnalysis.analysis.dssp`."
            )
        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    def calculate(
        self,
        selection: str = "protein",
        stride: int = 1,
    ) -> DSSPResult:
        atoms = self.universe.select_atoms(selection)
        if len(atoms) == 0:
            raise ValueError(f"selection selected no atoms: {selection}")
        residues = atoms.residues
        n_residues = len(residues)

        indices = list(range(0, len(self.universe.trajectory), stride))
        n_frames = len(indices)
        logger.info(
            "Starting DSSP per-residue analysis: residues=%s frames=%s stride=%s",
            n_residues, n_frames, stride,
        )

        # Run MDAnalysis DSSP, passing only the relevant atom group.
        dssp = _MDADSSP(atoms).run(start=indices[0], stop=indices[-1] + 1, step=stride)
        codes = dssp.results.dssp  # shape (n_sampled_frames, n_residues)
        if codes.shape[1] != n_residues:
            # MDAnalysis maps to whatever residues the AtomGroup covers; we trust
            # `atoms.residues` ordering and align by index.
            logger.warning(
                "DSSP residue count (%s) differs from AtomGroup residues (%s); "
                "aligning by minimum length.",
                codes.shape[1], n_residues,
            )

        # Per-residue accumulation
        records: List[Dict[str, object]] = []
        residues_used = min(codes.shape[1], n_residues)
        for res_idx in range(residues_used):
            residue = residues[res_idx]
            chain_id = self._chain_id(residue)
            resid = int(residue.resid)
            resname = str(residue.resname).strip()

            row = codes[:, res_idx].tolist()
            counter = Counter(row)
            n_obs = sum(counter.values())
            if n_obs == 0:
                continue

            helix = sum(counter[c] for c in self._HELIX_CODES)
            sheet = sum(counter[c] for c in self._SHEET_CODES)
            loop = n_obs - helix - sheet

            ss_propensity_helix = helix / n_obs
            ss_propensity_sheet = sheet / n_obs
            ss_propensity_loop = loop / n_obs

            # Dominant code mapped to canonical labels
            buckets = {
                "helix": ss_propensity_helix,
                "sheet": ss_propensity_sheet,
                "loop": ss_propensity_loop,
            }
            dominant_label = max(buckets, key=buckets.get)
            ss_stability = buckets[dominant_label]
            secondary_structure = {"helix": "H", "sheet": "E", "loop": "L"}[dominant_label]

            records.append(
                {
                    "chain_id": chain_id,
                    "resid": resid,
                    "resname": resname,
                    "secondary_structure": secondary_structure,
                    "ss_propensity_helix": float(ss_propensity_helix),
                    "ss_propensity_sheet": float(ss_propensity_sheet),
                    "ss_propensity_loop": float(ss_propensity_loop),
                    "ss_stability": float(ss_stability),
                    "n_frames_observed": int(n_obs),
                }
            )

        residue_frame = pd.DataFrame.from_records(records).sort_values(
            ["chain_id", "resid"]
        ).reset_index(drop=True)

        summary = self._build_summary(residue_frame, selection, stride, n_frames)

        return DSSPResult(residue_frame=residue_frame, summary=summary, stride=stride)

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
                "ss_counts": {},
            }
        counts = residue_frame["secondary_structure"].value_counts().to_dict()
        return {
            "selection": selection,
            "stride": int(stride),
            "n_frames": int(n_frames),
            "n_residues": int(len(residue_frame)),
            "ss_counts": {str(k): int(v) for k, v in counts.items()},
            "mean_ss_stability": float(residue_frame["ss_stability"].mean(skipna=True)),
            "mean_helix_propensity": float(residue_frame["ss_propensity_helix"].mean(skipna=True)),
            "mean_sheet_propensity": float(residue_frame["ss_propensity_sheet"].mean(skipna=True)),
            "mean_loop_propensity": float(residue_frame["ss_propensity_loop"].mean(skipna=True)),
        }
