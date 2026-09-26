"""
Residue-residue contact frequency analysis for pHLA-TCR trajectories.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
import logging

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd
from ..interactions.occupancy_metrics import encode_segments, frame_set_to_segments, summarize_frames

logger = logging.getLogger(__name__)


class ResidueContactFrequencyAnalyzer:
    """
    Calculate residue-residue contact frequencies between two selections.

    A contact is defined as: any heavy-atom pair between two residues has a
    minimum distance below the specified cutoff in a frame.
    """

    def __init__(self, topology: str, trajectory: str):
        self.topology = topology
        self.trajectory = trajectory
        self.universe = mda.Universe(topology, trajectory)

    def calculate_residue_contact_frequencies(
        self,
        selection1: str,
        selection2: str,
        cutoff: float = 4.5,
        stride: int = 1,
        heavy_atoms_only: bool = True,
        min_frequency: float = 0.0,
    ) -> pd.DataFrame:
        """
        Calculate contact frequency for every residue pair across a trajectory.
        """
        atoms1 = self._select_atoms(selection1, heavy_atoms_only=heavy_atoms_only)
        atoms2 = self._select_atoms(selection2, heavy_atoms_only=heavy_atoms_only)

        if len(atoms1.residues) == 0 or len(atoms2.residues) == 0:
            raise ValueError("Selections must contain at least one residue each")

        residues1 = list(atoms1.residues)
        residues2 = list(atoms2.residues)

        res1_atom_indices = self._build_local_atom_index_map(atoms1, residues1)
        res2_atom_indices = self._build_local_atom_index_map(atoms2, residues2)

        total_frames = 0
        contact_counts: Dict[tuple, int] = {}
        min_distances: Dict[tuple, float] = {}
        contact_frame_sets: Dict[tuple, set[int]] = {}

        for ts in self.universe.trajectory[::stride]:
            total_frames += 1
            pos1 = atoms1.positions
            pos2 = atoms2.positions
            frame_distances = distance_array(pos1, pos2, box=ts.dimensions)

            for i, residue1 in enumerate(residues1):
                atom_idx1 = res1_atom_indices[i]
                if not atom_idx1:
                    continue

                for j, residue2 in enumerate(residues2):
                    atom_idx2 = res2_atom_indices[j]
                    if not atom_idx2:
                        continue

                    pair_key = (
                        self._residue_chain_id(residue1),
                        int(residue1.resid),
                        residue1.resname,
                        self._residue_chain_id(residue2),
                        int(residue2.resid),
                        residue2.resname,
                    )

                    sub_distances = frame_distances[np.ix_(atom_idx1, atom_idx2)]
                    min_distance = float(np.min(sub_distances))

                    previous_min = min_distances.get(pair_key)
                    if previous_min is None or min_distance < previous_min:
                        min_distances[pair_key] = min_distance

                    if min_distance < cutoff:
                        contact_counts[pair_key] = contact_counts.get(pair_key, 0) + 1
                        contact_frame_sets.setdefault(pair_key, set()).add(int(ts.frame))

        if total_frames == 0:
            raise ValueError("Trajectory does not contain any frames")

        rows: List[Dict[str, object]] = []
        for pair_key, contact_frames in contact_counts.items():
            frequency = contact_frames / total_frames
            if frequency < min_frequency:
                continue

            chain1, resid1, resname1, chain2, resid2, resname2 = pair_key
            frames_seen = contact_frame_sets.get(pair_key, set())
            segments = frame_set_to_segments(frames_seen)
            segment_summary = summarize_frames(frames_seen)
            rows.append(
                {
                    "chain_id_1": chain1,
                    "resid_1": resid1,
                    "resname_1": resname1,
                    "residue_label_1": f"{resname1}{resid1}",
                    "chain_id_2": chain2,
                    "resid_2": resid2,
                    "resname_2": resname2,
                    "residue_label_2": f"{resname2}{resid2}",
                    "contact_frames": int(contact_frames),
                    "total_frames": int(total_frames),
                    "contact_frequency": float(frequency),
                    "min_distance_observed": float(min_distances[pair_key]),
                    "frame_segments": encode_segments(segments),
                    "n_segments": int(segment_summary["n_segments"]),
                    "max_consecutive_frames": int(segment_summary["max_consecutive_frames"]),
                    "first_frame": segment_summary["first_frame"],
                    "last_frame": segment_summary["last_frame"],
                    "cutoff_angstrom": float(cutoff),
                    "selection_1": selection1,
                    "selection_2": selection2,
                }
            )

        if not rows:
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
                    "contact_frames",
                    "total_frames",
                    "contact_frequency",
                    "min_distance_observed",
                    "frame_segments",
                    "n_segments",
                    "max_consecutive_frames",
                    "first_frame",
                    "last_frame",
                    "cutoff_angstrom",
                    "selection_1",
                    "selection_2",
                ]
            )

        return pd.DataFrame(rows).sort_values(
            by=["contact_frequency", "chain_id_1", "resid_1", "chain_id_2", "resid_2"],
            ascending=[False, True, True, True, True],
        ).reset_index(drop=True)

    @staticmethod
    def summarize_by_residue(
        contacts: pd.DataFrame,
        residue_side: int,
        frequency_column: str = "contact_frequency",
    ) -> pd.DataFrame:
        """
        Aggregate partner frequencies by residue on one side of the pair table.
        """
        if residue_side not in (1, 2):
            raise ValueError("residue_side must be 1 or 2")

        prefix = f"_{residue_side}"
        required = [
            f"chain_id{prefix}",
            f"resid{prefix}",
            f"resname{prefix}",
            f"residue_label{prefix}",
            frequency_column,
        ]
        missing = [column for column in required if column not in contacts.columns]
        if missing:
            raise ValueError(f"Missing required columns for summary: {missing}")

        grouped = (
            contacts.groupby(
                [
                    f"chain_id{prefix}",
                    f"resid{prefix}",
                    f"resname{prefix}",
                    f"residue_label{prefix}",
                ],
                as_index=False,
            )[frequency_column]
            .agg(["sum", "max", "count"])
            .reset_index()
        )
        grouped.columns = [
            f"chain_id{prefix}",
            f"resid{prefix}",
            f"resname{prefix}",
            f"residue_label{prefix}",
            "frequency_sum",
            "frequency_max",
            "n_partner_residues",
        ]
        return grouped.sort_values(
            by=["frequency_sum", "frequency_max"],
            ascending=[False, False],
        ).reset_index(drop=True)

    def _select_atoms(self, selection: str, heavy_atoms_only: bool = True):
        final_selection = selection
        if heavy_atoms_only:
            final_selection = f"({selection}) and not name H*"

        atoms = self.universe.select_atoms(final_selection)
        if len(atoms) == 0:
            raise ValueError(f"Selection returned no atoms: {selection}")
        return atoms

    @staticmethod
    def _build_local_atom_index_map(atoms, residues) -> List[List[int]]:
        atom_to_local = {atom.index: i for i, atom in enumerate(atoms)}
        residue_indices: List[List[int]] = []
        for residue in residues:
            local_indices = [
                atom_to_local[atom.index]
                for atom in residue.atoms
                if atom.index in atom_to_local
            ]
            residue_indices.append(local_indices)
        return residue_indices

    @staticmethod
    def _residue_chain_id(residue) -> str:
        segid = str(getattr(residue, "segid", "")).strip()
        chain_id = str(getattr(residue, "chainID", "")).strip()
        return chain_id or segid or "UNKNOWN"
