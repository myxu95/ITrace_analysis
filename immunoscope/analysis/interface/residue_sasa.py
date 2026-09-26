"""Per-residue SASA analyzer.

Computes solvent-accessible surface area for every residue in two forms:

  * `sasa_bound`   — SASA of the residue in the full complex
  * `sasa_unbound` — SASA of the residue when only its own group (TCR or
                     pMHC) is present

The pair (bound, unbound) reveals how much surface each residue loses to the
interface (`delta_sasa = sasa_unbound - sasa_bound`) and what fraction is
still exposed when complexed (`relative_exposure = sasa_bound / sasa_unbound`).

For mutation design these answer "is there room to put a bigger residue here?"
in a way BSA at the chain level cannot: BSA tells you how much total interface
exists, residue-level SASA tells you which positions are buried by the partner
versus merely framing the interface.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Dict, Iterable, List, Optional, Tuple

import MDAnalysis as mda
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import freesasa

    FREESASA_AVAILABLE = True
except ImportError:  # pragma: no cover
    FREESASA_AVAILABLE = False


# Residue key in the per-frame accumulator: (chain_id, resid, resname).
_ResKey = Tuple[str, int, str]


@dataclass(frozen=True)
class ResidueSasaResult:
    """Per-residue SASA result averaged over the trajectory."""

    residue_frame: pd.DataFrame
    summary: dict
    stride: int


class ResidueSasaAnalyzer:
    """Calculate per-residue bound / unbound SASA along an MD trajectory."""

    def __init__(self, topology: str, trajectory: str):
        if not FREESASA_AVAILABLE:
            raise ValueError("FreeSASA is not available. Install freesasa first.")
        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    def calculate(
        self,
        selection_a: str,
        selection_b: str,
        probe_radius: float = 1.4,
        stride: int = 1,
    ) -> ResidueSasaResult:
        atoms_a = self.universe.select_atoms(selection_a)
        atoms_b = self.universe.select_atoms(selection_b)
        if len(atoms_a) == 0:
            raise ValueError(f"selection_a selected no atoms: {selection_a}")
        if len(atoms_b) == 0:
            raise ValueError(f"selection_b selected no atoms: {selection_b}")
        atoms_complex = atoms_a | atoms_b

        parameters = freesasa.Parameters()
        parameters.setProbeRadius(probe_radius)

        # Accumulators: per-residue running sum + frame count.
        bound_total: Dict[_ResKey, float] = {}
        bound_side: Dict[_ResKey, float] = {}
        unbound_total: Dict[_ResKey, float] = {}
        seen_frames: Dict[_ResKey, int] = {}
        side: Dict[_ResKey, str] = {}   # "A" or "B"

        indices = list(range(0, len(self.universe.trajectory), stride))
        logger.info(
            "Starting residue SASA calculation: frames=%s stride=%s probe=%.2f",
            len(indices), stride, probe_radius,
        )

        for frame_index in indices:
            self.universe.trajectory[frame_index]

            complex_areas = self._residue_areas(atoms_complex, parameters)
            a_areas = self._residue_areas(atoms_a, parameters)
            b_areas = self._residue_areas(atoms_b, parameters)

            for key, area in complex_areas.items():
                bound_total[key] = bound_total.get(key, 0.0) + area.total
                bound_side[key] = bound_side.get(key, 0.0) + area.sideChain
                seen_frames[key] = seen_frames.get(key, 0) + 1

            for key, area in a_areas.items():
                unbound_total[key] = unbound_total.get(key, 0.0) + area.total
                side[key] = "A"
            for key, area in b_areas.items():
                unbound_total[key] = unbound_total.get(key, 0.0) + area.total
                side[key] = "B"

        n_frames = len(indices)
        residue_frame = self._build_frame(
            bound_total, bound_side, unbound_total, seen_frames, side, n_frames
        )

        summary = self._build_summary(
            residue_frame, selection_a, selection_b, probe_radius, stride, n_frames
        )

        return ResidueSasaResult(
            residue_frame=residue_frame,
            summary=summary,
            stride=stride,
        )

    # ------------------------------------------------------------------ #
    # Per-frame SASA
    # ------------------------------------------------------------------ #

    def _residue_areas(
        self,
        atoms: mda.AtomGroup,
        parameters: "freesasa.Parameters",
    ) -> Dict[_ResKey, "freesasa.ResidueArea"]:
        structure = freesasa.Structure()
        keys_in_order: List[_ResKey] = []
        seen: Dict[Tuple[str, int], _ResKey] = {}
        for atom in atoms:
            chain_id = self._get_chain_id(atom)
            resid = int(atom.resid)
            resname = str(atom.resname).strip()
            key = (chain_id, resid, resname)
            if (chain_id, resid) not in seen:
                seen[(chain_id, resid)] = key
                keys_in_order.append(key)
            structure.addAtom(
                str(atom.name).strip(),
                resname,
                str(resid),
                chain_id,
                float(atom.position[0]),
                float(atom.position[1]),
                float(atom.position[2]),
            )

        result = freesasa.calc(structure, parameters)
        residue_areas = result.residueAreas()  # {chain: {resid_str: ResidueArea}}

        flat: Dict[_ResKey, "freesasa.ResidueArea"] = {}
        for chain_id, resid, resname in keys_in_order:
            chain_table = residue_areas.get(chain_id)
            if not chain_table:
                continue
            entry = chain_table.get(str(resid))
            if entry is None:
                continue
            flat[(chain_id, resid, resname)] = entry
        return flat

    # ------------------------------------------------------------------ #
    # Aggregation
    # ------------------------------------------------------------------ #

    @staticmethod
    def _build_frame(
        bound_total: Dict[_ResKey, float],
        bound_side: Dict[_ResKey, float],
        unbound_total: Dict[_ResKey, float],
        seen_frames: Dict[_ResKey, int],
        side: Dict[_ResKey, str],
        n_frames: int,
    ) -> pd.DataFrame:
        records = []
        keys = set(bound_total.keys()) | set(unbound_total.keys())
        for key in keys:
            chain_id, resid, resname = key
            frames = seen_frames.get(key, 0)
            sasa_bound = bound_total.get(key, 0.0) / frames if frames else float("nan")
            sasa_side = bound_side.get(key, 0.0) / frames if frames else float("nan")
            unbound_frames = n_frames  # unbound computed every frame
            sasa_unbound = (
                unbound_total.get(key, 0.0) / unbound_frames if unbound_frames else float("nan")
            )
            if np.isnan(sasa_bound) or np.isnan(sasa_unbound):
                delta = float("nan")
                rel_exp = float("nan")
            else:
                delta = sasa_unbound - sasa_bound
                rel_exp = sasa_bound / sasa_unbound if sasa_unbound > 0 else float("nan")
            burial_state = ResidueSasaAnalyzer._classify_burial(delta, rel_exp, sasa_unbound)
            records.append(
                {
                    "chain_id": chain_id,
                    "resid": resid,
                    "resname": resname,
                    "side": side.get(key, ""),
                    "sasa_bound": float(sasa_bound),
                    "sasa_unbound": float(sasa_unbound),
                    "delta_sasa": float(delta),
                    "relative_exposure": float(rel_exp),
                    "sasa_sidechain_bound": float(sasa_side),
                    "burial_state": burial_state,
                    "n_frames_observed": int(frames),
                }
            )
        return pd.DataFrame.from_records(records).sort_values(
            ["chain_id", "resid"]
        ).reset_index(drop=True)

    @staticmethod
    def _classify_burial(delta: float, rel_exposure: float, sasa_unbound: float) -> str:
        if np.isnan(delta) or np.isnan(rel_exposure):
            return "unknown"
        if sasa_unbound < 10.0:
            return "core"
        if delta < 10.0:
            return "exposed"
        if rel_exposure < 0.3:
            return "interface_core"
        return "interface_rim"

    @staticmethod
    def _build_summary(
        residue_frame: pd.DataFrame,
        selection_a: str,
        selection_b: str,
        probe_radius: float,
        stride: int,
        n_frames: int,
    ) -> dict:
        counts = residue_frame["burial_state"].value_counts().to_dict()
        return {
            "selection_a": selection_a,
            "selection_b": selection_b,
            "probe_radius": float(probe_radius),
            "stride": int(stride),
            "n_frames": int(n_frames),
            "n_residues": int(len(residue_frame)),
            "burial_state_counts": {state: int(n) for state, n in counts.items()},
            "mean_sasa_bound": float(residue_frame["sasa_bound"].mean(skipna=True)),
            "mean_sasa_unbound": float(residue_frame["sasa_unbound"].mean(skipna=True)),
            "mean_delta_sasa": float(residue_frame["delta_sasa"].mean(skipna=True)),
            "mean_relative_exposure": float(
                residue_frame["relative_exposure"].mean(skipna=True)
            ),
        }

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

    @staticmethod
    def _get_chain_id(atom) -> str:
        for attr in ("chainID", "segid"):
            value: Optional[str] = getattr(atom, attr, None)
            if value:
                value = str(value).strip()
                if value:
                    return value[0]
        return "A"
