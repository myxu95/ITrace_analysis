"""
pi-pi and cation-pi residue-pair analyzer for pHLA-TCR complexes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.exceptions import NoDataError
from .occupancy_metrics import encode_segments, frame_set_to_segments, summarize_frames


AROMATIC_RING_ATOMS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HID": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HIE": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HIP": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HSD": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HSE": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HSP": ("CG", "ND1", "CD2", "CE1", "NE2"),
}

CATION_GROUP_ATOMS = {
    "LYS": ("NZ",),
    "ARG": ("NE", "CZ", "NH1", "NH2"),
    "HIP": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HSP": ("CG", "ND1", "CD2", "CE1", "NE2"),
}


@dataclass
class _RingGroup:
    residue_index: int
    atom_indices: np.ndarray
    chain_id: str
    resid: int
    resname: str
    residue_label: str


@dataclass
class _CationGroup:
    residue_index: int
    atom_indices: np.ndarray
    chain_id: str
    resid: int
    resname: str
    residue_label: str


@dataclass
class _PairAccumulator:
    chain_id_1: str
    resid_1: int
    resname_1: str
    residue_label_1: str
    chain_id_2: str
    resid_2: int
    resname_2: str
    residue_label_2: str
    min_distance_observed: float
    min_geometry_angle_observed: float
    frames_seen: set[int]


class _BasePiAnalyzer:
    def __init__(self, topology: str, trajectory: str, reference_structure: Optional[str] = None):
        self.universe = mda.Universe(topology, trajectory)
        self.reference_universe = mda.Universe(reference_structure) if reference_structure else None

    def _get_residue_metadata(self, atom) -> dict:
        residue = atom.residue
        try:
            chain_id = atom.chainID
        except NoDataError:
            if self.reference_universe is None:
                raise
            ref_atom = self.reference_universe.atoms[atom.index]
            ref_residue = ref_atom.residue
            return {
                "chain_id": ref_atom.chainID,
                "resid": int(ref_residue.resid),
                "resname": str(ref_residue.resname),
            }
        return {
            "chain_id": chain_id,
            "resid": int(residue.resid),
            "resname": str(residue.resname),
        }

    def _build_ring_groups(self, chain_ids: Iterable[str]) -> list[_RingGroup]:
        chain_ids = set(chain_ids)
        groups: list[_RingGroup] = []
        for residue in self.universe.residues:
            try:
                chain_id = residue.atoms[0].chainID
            except NoDataError:
                if self.reference_universe is None:
                    continue
                chain_id = self.reference_universe.atoms[residue.atoms[0].index].chainID
            if chain_id not in chain_ids:
                continue
            ring_atoms = AROMATIC_RING_ATOMS.get(str(residue.resname))
            if not ring_atoms:
                continue
            atom_indices = [atom.index for atom in residue.atoms if atom.name in ring_atoms]
            if len(atom_indices) < 5:
                continue
            meta = self._get_residue_metadata(residue.atoms[0])
            groups.append(
                _RingGroup(
                    residue_index=int(residue.ix),
                    atom_indices=np.asarray(atom_indices, dtype=int),
                    chain_id=meta["chain_id"],
                    resid=meta["resid"],
                    resname=meta["resname"],
                    residue_label=f"{meta['resname']}{meta['resid']}",
                )
            )
        return groups

    def _build_cation_groups(self, chain_ids: Iterable[str]) -> list[_CationGroup]:
        chain_ids = set(chain_ids)
        groups: list[_CationGroup] = []
        for residue in self.universe.residues:
            try:
                chain_id = residue.atoms[0].chainID
            except NoDataError:
                if self.reference_universe is None:
                    continue
                chain_id = self.reference_universe.atoms[residue.atoms[0].index].chainID
            if chain_id not in chain_ids:
                continue
            cation_atoms = CATION_GROUP_ATOMS.get(str(residue.resname))
            if not cation_atoms:
                continue
            atom_indices = [atom.index for atom in residue.atoms if atom.name in cation_atoms]
            if len(atom_indices) == 0:
                continue
            meta = self._get_residue_metadata(residue.atoms[0])
            groups.append(
                _CationGroup(
                    residue_index=int(residue.ix),
                    atom_indices=np.asarray(atom_indices, dtype=int),
                    chain_id=meta["chain_id"],
                    resid=meta["resid"],
                    resname=meta["resname"],
                    residue_label=f"{meta['resname']}{meta['resid']}",
                )
            )
        return groups

    @staticmethod
    def _ring_centroid_and_normal(positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        centroid = positions.mean(axis=0)
        centered = positions - centroid
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        normal = vh[-1]
        normal = normal / np.linalg.norm(normal)
        return centroid, normal

    @staticmethod
    def _acute_angle(v1: np.ndarray, v2: np.ndarray) -> float:
        cos_theta = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)
        theta = float(np.degrees(np.arccos(cos_theta)))
        return min(theta, 180.0 - theta)

    @staticmethod
    def _empty_result() -> pd.DataFrame:
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
                "min_geometry_angle_observed",
                "frame_segments",
                "n_segments",
                "max_consecutive_frames",
                "first_frame",
                "last_frame",
                "interaction_family",
                "distance_cutoff_angstrom",
            ]
        )


class PiStackingPairAnalyzer(_BasePiAnalyzer):
    """Calculate residue-level pi-pi interaction occupancy across the pHLA-TCR interface."""

    def calculate_interface_pairs(
        self,
        chain_mapping: Dict[str, str],
        stride: int = 1,
        distance_cutoff: float = 6.5,
        parallel_angle_cutoff: float = 30.0,
        tshape_angle_min: float = 60.0,
        tshape_angle_max: float = 120.0,
    ) -> pd.DataFrame:
        phla_chains = {chain_mapping["mhc_alpha"], chain_mapping["b2m"], chain_mapping["peptide"]}
        tcr_chains = {chain_mapping["tcr_alpha"], chain_mapping["tcr_beta"]}
        phla_rings = self._build_ring_groups(phla_chains)
        tcr_rings = self._build_ring_groups(tcr_chains)
        total_frames = len(range(0, len(self.universe.trajectory), stride))
        if not phla_rings or not tcr_rings:
            return self._empty_result()

        pair_map: dict[tuple, _PairAccumulator] = {}
        for ts in self.universe.trajectory[::stride]:
            frame_index = int(ts.frame)
            ring_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

            def ring_geom(group: _RingGroup):
                if group.residue_index not in ring_cache:
                    ring_cache[group.residue_index] = self._ring_centroid_and_normal(self.universe.atoms[group.atom_indices].positions)
                return ring_cache[group.residue_index]

            for phla_group in phla_rings:
                phla_centroid, phla_normal = ring_geom(phla_group)
                for tcr_group in tcr_rings:
                    tcr_centroid, tcr_normal = ring_geom(tcr_group)
                    centroid_distance = float(np.linalg.norm(phla_centroid - tcr_centroid))
                    if centroid_distance > distance_cutoff:
                        continue
                    plane_angle = self._acute_angle(phla_normal, tcr_normal)
                    if not (plane_angle <= parallel_angle_cutoff or tshape_angle_min <= plane_angle <= tshape_angle_max):
                        continue
                    key = (phla_group.chain_id, phla_group.resid, tcr_group.chain_id, tcr_group.resid)
                    if key not in pair_map:
                        pair_map[key] = _PairAccumulator(
                            chain_id_1=phla_group.chain_id,
                            resid_1=phla_group.resid,
                            resname_1=phla_group.resname,
                            residue_label_1=phla_group.residue_label,
                            chain_id_2=tcr_group.chain_id,
                            resid_2=tcr_group.resid,
                            resname_2=tcr_group.resname,
                            residue_label_2=tcr_group.residue_label,
                            min_distance_observed=centroid_distance,
                            min_geometry_angle_observed=plane_angle,
                            frames_seen={frame_index},
                        )
                    else:
                        acc = pair_map[key]
                        acc.frames_seen.add(frame_index)
                        acc.min_distance_observed = min(acc.min_distance_observed, centroid_distance)
                        acc.min_geometry_angle_observed = min(acc.min_geometry_angle_observed, plane_angle)

        if not pair_map:
            return self._empty_result()

        rows = []
        for acc in pair_map.values():
            contact_frames = len(acc.frames_seen)
            segments = frame_set_to_segments(acc.frames_seen)
            segment_summary = summarize_frames(acc.frames_seen)
            rows.append(
                {
                    "chain_id_1": acc.chain_id_1,
                    "resid_1": acc.resid_1,
                    "resname_1": acc.resname_1,
                    "residue_label_1": acc.residue_label_1,
                    "chain_id_2": acc.chain_id_2,
                    "resid_2": acc.resid_2,
                    "resname_2": acc.resname_2,
                    "residue_label_2": acc.residue_label_2,
                    "contact_frames": contact_frames,
                    "total_frames": total_frames,
                    "contact_frequency": float(contact_frames / total_frames) if total_frames else 0.0,
                    "min_distance_observed": acc.min_distance_observed,
                    "min_geometry_angle_observed": acc.min_geometry_angle_observed,
                    "frame_segments": encode_segments(segments),
                    "n_segments": int(segment_summary["n_segments"]),
                    "max_consecutive_frames": int(segment_summary["max_consecutive_frames"]),
                    "first_frame": segment_summary["first_frame"],
                    "last_frame": segment_summary["last_frame"],
                    "interaction_family": "pi_pi",
                    "distance_cutoff_angstrom": distance_cutoff,
                }
            )
        return pd.DataFrame(rows).sort_values(
            by=["contact_frequency", "min_distance_observed", "min_geometry_angle_observed"],
            ascending=[False, True, True],
        ).reset_index(drop=True)


class CationPiPairAnalyzer(_BasePiAnalyzer):
    """Calculate residue-level cation-pi interaction occupancy across the pHLA-TCR interface."""

    def calculate_interface_pairs(
        self,
        chain_mapping: Dict[str, str],
        stride: int = 1,
        distance_cutoff: float = 6.0,
        normal_angle_cutoff: float = 60.0,
    ) -> pd.DataFrame:
        phla_chains = {chain_mapping["mhc_alpha"], chain_mapping["b2m"], chain_mapping["peptide"]}
        tcr_chains = {chain_mapping["tcr_alpha"], chain_mapping["tcr_beta"]}
        phla_rings = self._build_ring_groups(phla_chains)
        tcr_rings = self._build_ring_groups(tcr_chains)
        phla_cations = self._build_cation_groups(phla_chains)
        tcr_cations = self._build_cation_groups(tcr_chains)
        total_frames = len(range(0, len(self.universe.trajectory), stride))
        if (not phla_rings and not tcr_rings) or (not phla_cations and not tcr_cations):
            return self._empty_result()

        pair_map: dict[tuple, _PairAccumulator] = {}
        for ts in self.universe.trajectory[::stride]:
            frame_index = int(ts.frame)
            ring_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}
            cation_cache: dict[int, np.ndarray] = {}

            def ring_geom(group: _RingGroup):
                if group.residue_index not in ring_cache:
                    ring_cache[group.residue_index] = self._ring_centroid_and_normal(self.universe.atoms[group.atom_indices].positions)
                return ring_cache[group.residue_index]

            def cation_center(group: _CationGroup):
                if group.residue_index not in cation_cache:
                    cation_cache[group.residue_index] = self.universe.atoms[group.atom_indices].positions.mean(axis=0)
                return cation_cache[group.residue_index]

            for ring_group, cation_group, ring_on_phla in (
                *((ring, cat, True) for ring in phla_rings for cat in tcr_cations),
                *((ring, cat, False) for ring in tcr_rings for cat in phla_cations),
            ):
                ring_centroid, ring_normal = ring_geom(ring_group)
                cat_center = cation_center(cation_group)
                centroid_distance = float(np.linalg.norm(ring_centroid - cat_center))
                if centroid_distance > distance_cutoff:
                    continue
                vector = cat_center - ring_centroid
                if np.linalg.norm(vector) == 0:
                    continue
                normal_angle = self._acute_angle(vector, ring_normal)
                if normal_angle > normal_angle_cutoff:
                    continue
                if ring_on_phla:
                    left_chain, left_resid, left_resname, left_label = (
                        ring_group.chain_id,
                        ring_group.resid,
                        ring_group.resname,
                        ring_group.residue_label,
                    )
                    right_chain, right_resid, right_resname, right_label = (
                        cation_group.chain_id,
                        cation_group.resid,
                        cation_group.resname,
                        cation_group.residue_label,
                    )
                else:
                    left_chain, left_resid, left_resname, left_label = (
                        cation_group.chain_id,
                        cation_group.resid,
                        cation_group.resname,
                        cation_group.residue_label,
                    )
                    right_chain, right_resid, right_resname, right_label = (
                        ring_group.chain_id,
                        ring_group.resid,
                        ring_group.resname,
                        ring_group.residue_label,
                    )
                key = (left_chain, left_resid, right_chain, right_resid)
                if key not in pair_map:
                    pair_map[key] = _PairAccumulator(
                        chain_id_1=left_chain,
                        resid_1=left_resid,
                        resname_1=left_resname,
                        residue_label_1=left_label,
                        chain_id_2=right_chain,
                        resid_2=right_resid,
                        resname_2=right_resname,
                        residue_label_2=right_label,
                        min_distance_observed=centroid_distance,
                        min_geometry_angle_observed=normal_angle,
                        frames_seen={frame_index},
                    )
                else:
                    acc = pair_map[key]
                    acc.frames_seen.add(frame_index)
                    acc.min_distance_observed = min(acc.min_distance_observed, centroid_distance)
                    acc.min_geometry_angle_observed = min(acc.min_geometry_angle_observed, normal_angle)

        if not pair_map:
            return self._empty_result()

        rows = []
        for acc in pair_map.values():
            contact_frames = len(acc.frames_seen)
            segments = frame_set_to_segments(acc.frames_seen)
            segment_summary = summarize_frames(acc.frames_seen)
            rows.append(
                {
                    "chain_id_1": acc.chain_id_1,
                    "resid_1": acc.resid_1,
                    "resname_1": acc.resname_1,
                    "residue_label_1": acc.residue_label_1,
                    "chain_id_2": acc.chain_id_2,
                    "resid_2": acc.resid_2,
                    "resname_2": acc.resname_2,
                    "residue_label_2": acc.residue_label_2,
                    "contact_frames": contact_frames,
                    "total_frames": total_frames,
                    "contact_frequency": float(contact_frames / total_frames) if total_frames else 0.0,
                    "min_distance_observed": acc.min_distance_observed,
                    "min_geometry_angle_observed": acc.min_geometry_angle_observed,
                    "frame_segments": encode_segments(segments),
                    "n_segments": int(segment_summary["n_segments"]),
                    "max_consecutive_frames": int(segment_summary["max_consecutive_frames"]),
                    "first_frame": segment_summary["first_frame"],
                    "last_frame": segment_summary["last_frame"],
                    "interaction_family": "cation_pi",
                    "distance_cutoff_angstrom": distance_cutoff,
                }
            )
        return pd.DataFrame(rows).sort_values(
            by=["contact_frequency", "min_distance_observed", "min_geometry_angle_observed"],
            ascending=[False, True, True],
        ).reset_index(drop=True)
