"""RRCS (Residue-Residue Contact Score) analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O"}


@dataclass(frozen=True)
class RRCSPairSpec:
    """RRCS calculation specification for one residue pair."""

    chain_id_1: str
    resid_1: int
    resname_1: str
    residue_label_1: str
    chain_id_2: str
    resid_2: int
    resname_2: str
    residue_label_2: str
    atom_indices_1: np.ndarray
    atom_indices_2: np.ndarray
    atom_indices_1_no_backbone: np.ndarray
    atom_indices_2_no_backbone: np.ndarray
    selection_scope: str = "interface"


class RRCSAnalyzer:
    """Calculate per-frame RRCS time series for residue pairs."""

    def __init__(self, topology: str, trajectory: str):
        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    def calculate_pair_timeseries(
        self,
        pair_specs: Iterable[RRCSPairSpec],
        radius_min: float = 3.23,
        radius_max: float = 4.63,
        stride: int = 1,
        exclude_neighbor_backbone: bool = True,
        use_prefilter: bool = True,
    ) -> pd.DataFrame:
        """
        Calculate per-frame residue-pair RRCS.

        RRCS follows the gmx_RRCS convention:
        - non-hydrogen atom-pair distance <= radius_min contributes 1
        - non-hydrogen atom-pair distance >= radius_max contributes 0
        - intermediate distances decay linearly
        - residue-pair score is the sum over all atom-pair scores
        """
        if radius_min <= 0 or radius_max <= 0 or radius_max <= radius_min:
            raise ValueError("Invalid radius_min/radius_max; expected 0 < radius_min < radius_max")
        if stride < 1:
            raise ValueError("stride must be >= 1")

        pair_specs = list(pair_specs)
        if not pair_specs:
            raise ValueError("pair_specs must not be empty")

        rows: list[dict[str, object]] = []
        sampled_frames = 0

        for ts in self.universe.trajectory[::stride]:
            sampled_frames += 1
            positions = self.universe.atoms.positions
            time_ps = float(getattr(ts, "time", ts.frame))
            box = getattr(ts, "dimensions", None)

            for spec in pair_specs:
                atom_indices_1, atom_indices_2 = self._resolve_active_atom_indices(
                    spec=spec,
                    exclude_neighbor_backbone=exclude_neighbor_backbone,
                )
                if len(atom_indices_1) == 0 or len(atom_indices_2) == 0:
                    rrcs_score = 0.0
                else:
                    coords_1 = positions[atom_indices_1]
                    coords_2 = positions[atom_indices_2]
                    if use_prefilter and not self._passes_prefilter(coords_1, coords_2, radius_max):
                        rrcs_score = 0.0
                    else:
                        distances = distance_array(coords_1, coords_2, box=box)
                        rrcs_score = float(self._compute_rrcs_score(distances, radius_min, radius_max))

                rows.append(
                    {
                        "frame": int(ts.frame),
                        "time_ps": time_ps,
                        "chain_id_1": spec.chain_id_1,
                        "resid_1": int(spec.resid_1),
                        "resname_1": spec.resname_1,
                        "residue_label_1": spec.residue_label_1,
                        "chain_id_2": spec.chain_id_2,
                        "resid_2": int(spec.resid_2),
                        "resname_2": spec.resname_2,
                        "residue_label_2": spec.residue_label_2,
                        "selection_scope": spec.selection_scope,
                        "rrcs": rrcs_score,
                    }
                )

        if sampled_frames == 0:
            raise ValueError("Trajectory does not contain any sampleable frames")

        return pd.DataFrame(rows)

    @staticmethod
    def summarize_pair_timeseries(rrcs_timeseries: pd.DataFrame) -> pd.DataFrame:
        """Summarize per-frame RRCS values at pair level."""
        if rrcs_timeseries.empty:
            return pd.DataFrame(
                columns=[
                    "chain_id_1",
                    "resid_1",
                    "resname_1",
                    "residue_label_1",
                    "chain_id_2",
                    "resid_2",
                    "resname_2",
                    "residue_label_2",
                    "selection_scope",
                    "mean_rrcs",
                    "median_rrcs",
                    "max_rrcs",
                    "rrcs_nonzero_frames",
                    "rrcs_nonzero_fraction",
                    "total_frames",
                ]
            )

        group_columns = [
            "chain_id_1",
            "resid_1",
            "resname_1",
            "residue_label_1",
            "chain_id_2",
            "resid_2",
            "resname_2",
            "residue_label_2",
            "selection_scope",
        ]

        summary = (
            rrcs_timeseries.groupby(group_columns, as_index=False)
            .agg(
                mean_rrcs=("rrcs", "mean"),
                median_rrcs=("rrcs", "median"),
                max_rrcs=("rrcs", "max"),
                rrcs_nonzero_frames=("rrcs", lambda values: int(np.count_nonzero(np.asarray(values) > 0.0))),
                total_frames=("rrcs", "size"),
            )
        )
        summary["rrcs_nonzero_fraction"] = (
            summary["rrcs_nonzero_frames"] / summary["total_frames"].replace(0, np.nan)
        ).fillna(0.0)

        return summary.sort_values(
            by=["mean_rrcs", "max_rrcs", "chain_id_1", "resid_1", "chain_id_2", "resid_2"],
            ascending=[False, False, True, True, True, True],
        ).reset_index(drop=True)

    @staticmethod
    def _resolve_active_atom_indices(
        spec: RRCSPairSpec,
        exclude_neighbor_backbone: bool,
    ) -> tuple[np.ndarray, np.ndarray]:
        if (
            exclude_neighbor_backbone
            and spec.chain_id_1 == spec.chain_id_2
            and abs(int(spec.resid_1) - int(spec.resid_2)) < 5
        ):
            atom_indices_1 = spec.atom_indices_1_no_backbone
            atom_indices_2 = spec.atom_indices_2_no_backbone
            if len(atom_indices_1) and len(atom_indices_2):
                return atom_indices_1, atom_indices_2
        return spec.atom_indices_1, spec.atom_indices_2

    @staticmethod
    def _passes_prefilter(coords_1: np.ndarray, coords_2: np.ndarray, axis_cutoff: float) -> bool:
        diff = np.abs(coords_1[:, None, :] - coords_2[None, :, :])
        return bool(np.any(np.all(diff < axis_cutoff, axis=2)))

    @staticmethod
    def _compute_rrcs_score(distances: np.ndarray, radius_min: float, radius_max: float) -> float:
        score_matrix = np.where(
            distances <= radius_min,
            1.0,
            np.where(
                distances >= radius_max,
                0.0,
                (radius_max - distances) / (radius_max - radius_min),
            ),
        )
        return float(np.sum(score_matrix))


_PAIR_TOKEN_RE = re.compile(r"^(?P<chain>[A-Za-z0-9]+):(?P<start>-?\d+)(?:-(?P<end>-?\d+))?$")


def parse_rrcs_pair_file(pair_file: str | Path) -> list[tuple[str, int, str, int]]:
    """
    Parse an RRCS residue-pair file.

    Format:
    - one pair group per line, left and right sides separated by `$`
    - each side accepts one or more tokens separated by spaces or commas
    - token format is `CHAIN:RESID` or `CHAIN:START-END`
    - comment lines start with `#`

    Examples:
    - `D:95 $ C:4`
    - `D:95-100 $ C:4-6`
    - `D:95,D:96 $ A:150 A:151`
    """
    path = Path(pair_file)
    lines = path.read_text(encoding="utf-8").splitlines()
    pairs: list[tuple[str, int, str, int]] = []

    for line_no, raw_line in enumerate(lines, start=1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if "$" not in line:
            raise ValueError(f"RRCS pair file line {line_no} is missing '$' separator: {raw_line}")
        left_part, right_part = [part.strip() for part in line.split("$", 1)]
        left_tokens = _split_rrcs_side_tokens(left_part, line_no)
        right_tokens = _split_rrcs_side_tokens(right_part, line_no)
        left_residues = [res for token in left_tokens for res in _expand_rrcs_token(token, line_no)]
        right_residues = [res for token in right_tokens for res in _expand_rrcs_token(token, line_no)]
        for chain_1, resid_1 in left_residues:
            for chain_2, resid_2 in right_residues:
                pairs.append((chain_1, resid_1, chain_2, resid_2))

    if not pairs:
        raise ValueError(f"RRCS pair file is empty or contains no valid residue pairs: {path}")
    return pairs


def _split_rrcs_side_tokens(side: str, line_no: int) -> list[str]:
    normalized = side.replace(",", " ")
    tokens = [token.strip() for token in normalized.split() if token.strip()]
    if not tokens:
        raise ValueError(f"RRCS pair file line {line_no} contains an empty residue side")
    return tokens


def _expand_rrcs_token(token: str, line_no: int) -> list[tuple[str, int]]:
    match = _PAIR_TOKEN_RE.match(token)
    if not match:
        raise ValueError(f"RRCS pair file line {line_no} has an invalid token: {token}")
    chain = match.group("chain")
    start = int(match.group("start"))
    end = int(match.group("end") or start)
    if end < start:
        raise ValueError(f"RRCS pair file line {line_no} has an invalid range: {token}")
    return [(chain, resid) for resid in range(start, end + 1)]


def build_residue_heavy_atom_indices(residue) -> tuple[np.ndarray, np.ndarray]:
    """Return all heavy-atom indices and non-backbone heavy-atom indices."""
    heavy_atoms = residue.atoms.select_atoms("not name H*")
    heavy_no_backbone = residue.atoms.select_atoms("not name H* and not name N CA C O")
    return heavy_atoms.indices.copy(), heavy_no_backbone.indices.copy()


def residue_chain_id(residue) -> str:
    """Return a stable chain identifier for a residue."""
    segid = str(getattr(residue.segment, "segid", "")).strip()
    chain_id = str(getattr(residue, "chainID", "")).strip()
    return chain_id or segid or "UNKNOWN"
