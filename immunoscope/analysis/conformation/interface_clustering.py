"""pHLA-TCR key interface-state clustering analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import cdist, squareform


def _format_region_name(region: str) -> str:
    mapping = {
        "peptide": "peptide",
        "alpha1_helix": "MHC α1 helix",
        "alpha2_helix": "MHC α2 helix",
    }
    return mapping.get(region, region.replace("_", " "))


def _format_tcr_chain_name(tcr_name: str) -> str:
    mapping = {
        "TCR_alpha": "TCR α CDR3",
        "TCR_beta": "TCR β CDR3",
    }
    return mapping.get(tcr_name, tcr_name.replace("_", " "))


@dataclass
class InterfaceClusteringResult:
    """Interface clustering result container."""

    frame_assignments: pd.DataFrame
    cluster_summary: pd.DataFrame
    cluster_difference_matrix: pd.DataFrame
    cluster_feature_summary: pd.DataFrame
    cluster_signatures: pd.DataFrame
    cluster_feature_digest: pd.DataFrame
    summary: dict[str, Any]


class InterfaceClusteringAnalyzer:
    """Key interface-state clustering analyzer."""

    def __init__(self, topology_file: str, trajectory_file: str, structure_pdb: str):
        self.structure_universe = mda.Universe(structure_pdb)
        self.trajectory_universe = mda.Universe(topology_file, trajectory_file)
        _apply_protein_pbc_transformations(self.trajectory_universe)

        if self.structure_universe.atoms.n_atoms != self.trajectory_universe.atoms.n_atoms:
            raise ValueError(
                "structure_pdb and topology/trajectory atom counts differ; stable mapping cannot be built"
            )

    def calculate(
        self,
        chain_mapping: dict[str, str],
        cdr_detection: dict[str, Any],
        stride: int = 10,
        contact_cutoff_angstrom: float = 4.5,
        distance_cutoff: float = 0.35,
        linkage_method: str = "average",
        geometry_weight: float = 0.45,
        sidechain_weight: float = 0.25,
        interaction_weight: float = 0.30,
    ) -> InterfaceClusteringResult:
        contact_stride = _derive_contact_stride(stride)
        selections = self._prepare_selections(chain_mapping, cdr_detection)
        frame_records = self._extract_frame_records(
            selections=selections,
            geometry_stride=stride,
            contact_stride=contact_stride,
            contact_cutoff_angstrom=contact_cutoff_angstrom,
        )
        if len(frame_records) < 2:
            raise ValueError("At least 2 sampled frames are required for clustering")

        geometry_features = np.vstack([record["geometry"] for record in frame_records])
        sidechain_features = np.vstack([record["sidechain"] for record in frame_records])
        interaction_features = np.vstack([record["interaction"] for record in frame_records])

        d_geometry = _normalized_distance_matrix(geometry_features, metric="euclidean")
        d_sidechain = _normalized_distance_matrix(sidechain_features, metric="euclidean")
        d_interaction = _normalized_distance_matrix(interaction_features, metric="cityblock")

        d_joint = (
            geometry_weight * d_geometry
            + sidechain_weight * d_sidechain
            + interaction_weight * d_interaction
        )
        np.fill_diagonal(d_joint, 0.0)

        labels = self._cluster_frames(
            d_joint,
            distance_cutoff=distance_cutoff,
            linkage_method=linkage_method,
        )
        frame_assignments = self._build_frame_assignments(frame_records, labels)
        cluster_summary = self._build_cluster_summary(frame_records, frame_assignments, d_joint, labels)
        (
            cluster_difference_matrix,
            cluster_feature_summary,
            cluster_signatures,
            cluster_feature_digest,
        ) = self._build_cluster_explanations(
            frame_records=frame_records,
            labels=labels,
            joint_distance=d_joint,
            interaction_labels=selections["interaction_labels"],
            interaction_display_labels=selections["interaction_display_labels"],
            interaction_descriptions=selections["interaction_descriptions"],
        )
        cluster_summary = self._attach_structural_descriptors(cluster_summary, cluster_signatures)
        summary = self._build_summary(
            frame_assignments=frame_assignments,
            cluster_summary=cluster_summary,
            stride=stride,
            contact_stride=contact_stride,
            contact_cutoff_angstrom=contact_cutoff_angstrom,
            distance_cutoff=distance_cutoff,
            linkage_method=linkage_method,
            geometry_weight=geometry_weight,
            sidechain_weight=sidechain_weight,
            interaction_weight=interaction_weight,
        )

        return InterfaceClusteringResult(
            frame_assignments=frame_assignments,
            cluster_summary=cluster_summary,
            cluster_difference_matrix=cluster_difference_matrix,
            cluster_feature_summary=cluster_feature_summary,
            cluster_signatures=cluster_signatures,
            cluster_feature_digest=cluster_feature_digest,
            summary=summary,
        )

    def _prepare_selections(
        self,
        chain_mapping: dict[str, str],
        cdr_detection: dict[str, Any],
    ) -> dict[str, Any]:
        mhc_alpha_chain = chain_mapping.get("mhc_alpha")
        b2m_chain = chain_mapping.get("b2m")
        peptide_chain = chain_mapping.get("peptide")
        tcr_alpha_chain = chain_mapping.get("tcr_alpha")
        tcr_beta_chain = chain_mapping.get("tcr_beta")
        if not all([mhc_alpha_chain, peptide_chain, tcr_alpha_chain, tcr_beta_chain]):
            raise ValueError("chain_mapping is missing mhc_alpha/peptide/tcr_alpha/tcr_beta")

        cdr_chains = cdr_detection.get("chains", {})
        cdr3_sets: dict[str, set[int]] = {}
        for chain_name in ("TCR_alpha", "TCR_beta"):
            chain_info = cdr_chains.get(chain_name, {})
            cdr3_info = chain_info.get("cdrs", {}).get(3)
            if not cdr3_info:
                raise ValueError(f"{chain_name} is missing CDR3 detection results")
            start, end = cdr3_info["residue_range"]
            cdr3_sets[chain_name] = set(range(int(start), int(end) + 1))

        ref_positions = self.structure_universe.atoms.positions.copy()
        framework_ca_indices: list[int] = []
        geometry_indices: list[int] = []
        sidechain_groups: list[tuple[str, np.ndarray, np.ndarray]] = []

        peptide_residue_keys: list[tuple[str, int]] = []
        alpha1_residue_keys: list[tuple[str, int]] = []
        alpha2_residue_keys: list[tuple[str, int]] = []

        peptide_residue_atoms: dict[tuple[str, int], np.ndarray] = {}
        alpha1_residue_atoms: dict[tuple[str, int], np.ndarray] = {}
        alpha2_residue_atoms: dict[tuple[str, int], np.ndarray] = {}
        residue_details: dict[tuple[str, int], str] = {}
        cdr3_residue_atoms: list[tuple[str, str, int, str, np.ndarray]] = []
        complex_atom_indices: list[int] = []

        for residue in self.structure_universe.residues:
            chain_id = _residue_chain_id(residue)
            residue_key = (chain_id, int(residue.resid))
            if chain_id in {mhc_alpha_chain, b2m_chain, peptide_chain, tcr_alpha_chain, tcr_beta_chain}:
                complex_atom_indices.extend(residue.atoms.indices.tolist())

            if chain_id in (tcr_alpha_chain, tcr_beta_chain):
                tcr_name = "TCR_alpha" if chain_id == tcr_alpha_chain else "TCR_beta"
                ca_atoms = residue.atoms.select_atoms("name CA")
                if len(ca_atoms) and residue.resid not in cdr3_sets[tcr_name]:
                    framework_ca_indices.extend(ca_atoms.indices.tolist())

                if residue.resid in cdr3_sets[tcr_name]:
                    backbone = residue.atoms.select_atoms("backbone and not name H*")
                    if len(backbone):
                        geometry_indices.extend(backbone.indices.tolist())
                    sidechain = residue.atoms.select_atoms("not backbone and not name H*")
                    if len(ca_atoms) and len(sidechain):
                        label = f"{chain_id}:{int(residue.resid)}:{residue.resname}"
                        sidechain_groups.append((label, ca_atoms.indices.copy(), sidechain.indices.copy()))
                    heavy_atoms = residue.atoms.select_atoms("not name H*")
                    cdr3_residue_atoms.append(
                        (tcr_name, chain_id, int(residue.resid), residue.resname, heavy_atoms.indices.copy())
                    )
                    residue_details[(chain_id, int(residue.resid))] = residue.resname

            elif chain_id == peptide_chain:
                peptide_residue_keys.append(residue_key)
                heavy_atoms = residue.atoms.select_atoms("not name H*")
                peptide_residue_atoms[residue_key] = heavy_atoms.indices.copy()
                residue_details[residue_key] = residue.resname
                backbone = residue.atoms.select_atoms("backbone and not name H*")
                if len(backbone):
                    geometry_indices.extend(backbone.indices.tolist())
                ca_atoms = residue.atoms.select_atoms("name CA")
                sidechain = residue.atoms.select_atoms("not backbone and not name H*")
                if len(ca_atoms) and len(sidechain):
                    label = f"{chain_id}:{int(residue.resid)}:{residue.resname}"
                    sidechain_groups.append((label, ca_atoms.indices.copy(), sidechain.indices.copy()))

            elif chain_id == mhc_alpha_chain:
                heavy_atoms = residue.atoms.select_atoms("not name H*")
                if 50 <= residue.resid <= 86:
                    alpha1_residue_keys.append(residue_key)
                    alpha1_residue_atoms[residue_key] = heavy_atoms.indices.copy()
                    residue_details[residue_key] = residue.resname
                elif 140 <= residue.resid <= 176:
                    alpha2_residue_keys.append(residue_key)
                    alpha2_residue_atoms[residue_key] = heavy_atoms.indices.copy()
                    residue_details[residue_key] = residue.resname

        interaction_pairs: list[tuple[np.ndarray, np.ndarray]] = []
        interaction_labels: list[str] = []
        interaction_display_labels: list[str] = []
        interaction_descriptions: list[str] = []
        for src_tcr_name, src_chain, src_resid, src_resname, source_atoms in cdr3_residue_atoms:
            for target_key in peptide_residue_keys:
                interaction_pairs.append((source_atoms, peptide_residue_atoms[target_key]))
                interaction_labels.append(
                    f"{src_chain}:{src_resid}__{target_key[0]}:{target_key[1]}[peptide]"
                )
                interaction_display_labels.append(
                    f"{'TCRα' if src_tcr_name == 'TCR_alpha' else 'TCRβ'} {src_chain}:{src_resid} ↔ peptide {target_key[0]}:{target_key[1]}"
                )
                interaction_descriptions.append(
                    f"{_format_tcr_chain_name(src_tcr_name)} residue {src_chain}:{src_resid} {src_resname} contacts "
                    f"{_format_region_name('peptide')} residue {target_key[0]}:{target_key[1]} {residue_details[target_key]}"
                )
            for target_key in alpha1_residue_keys:
                interaction_pairs.append((source_atoms, alpha1_residue_atoms[target_key]))
                interaction_labels.append(
                    f"{src_chain}:{src_resid}__{target_key[0]}:{target_key[1]}[alpha1_helix]"
                )
                interaction_display_labels.append(
                    f"{'TCRα' if src_tcr_name == 'TCR_alpha' else 'TCRβ'} {src_chain}:{src_resid} ↔ α1 {target_key[0]}:{target_key[1]}"
                )
                interaction_descriptions.append(
                    f"{_format_tcr_chain_name(src_tcr_name)} residue {src_chain}:{src_resid} {src_resname} contacts "
                    f"{_format_region_name('alpha1_helix')} residue {target_key[0]}:{target_key[1]} {residue_details[target_key]}"
                )
            for target_key in alpha2_residue_keys:
                interaction_pairs.append((source_atoms, alpha2_residue_atoms[target_key]))
                interaction_labels.append(
                    f"{src_chain}:{src_resid}__{target_key[0]}:{target_key[1]}[alpha2_helix]"
                )
                interaction_display_labels.append(
                    f"{'TCRα' if src_tcr_name == 'TCR_alpha' else 'TCRβ'} {src_chain}:{src_resid} ↔ α2 {target_key[0]}:{target_key[1]}"
                )
                interaction_descriptions.append(
                    f"{_format_tcr_chain_name(src_tcr_name)} residue {src_chain}:{src_resid} {src_resname} contacts "
                    f"{_format_region_name('alpha2_helix')} residue {target_key[0]}:{target_key[1]} {residue_details[target_key]}"
                )

        if not framework_ca_indices:
            raise ValueError("No TCR framework CA atoms were found for alignment")
        if not geometry_indices:
            raise ValueError("No interface geometry feature atoms were found")
        if not sidechain_groups:
            raise ValueError("No key side-chain feature atoms were found")
        if not interaction_pairs:
            raise ValueError("No CDR3-to-peptide/alpha1/alpha2 interface relation features were found")

        return {
            "reference_positions": ref_positions,
            "framework_ca_indices": np.array(sorted(set(framework_ca_indices)), dtype=int),
            "geometry_indices": np.array(geometry_indices, dtype=int),
            "sidechain_groups": sidechain_groups,
            "interaction_pairs": interaction_pairs,
            "interaction_labels": interaction_labels,
            "interaction_display_labels": interaction_display_labels,
            "interaction_descriptions": interaction_descriptions,
            "complex_atom_indices": np.array(sorted(set(complex_atom_indices)), dtype=int),
        }

    def _extract_frame_records(
        self,
        selections: dict[str, Any],
        geometry_stride: int,
        contact_stride: int,
        contact_cutoff_angstrom: float,
    ) -> list[dict[str, Any]]:
        ref_positions = selections["reference_positions"]
        framework_indices = selections["framework_ca_indices"]
        ref_framework = ref_positions[framework_indices]
        geometry_indices = selections["geometry_indices"]
        contact_window_radius = max(0, geometry_stride // 2)

        contact_samples: list[dict[str, Any]] = []
        for ts in self.trajectory_universe.trajectory[::contact_stride]:
            mobile_positions = self.trajectory_universe.atoms.positions.copy()
            aligned_positions = _align_positions(
                all_positions=mobile_positions,
                mobile_reference=mobile_positions[framework_indices],
                target_reference=ref_framework,
            )
            interaction_binary, interaction_distance = self._compute_interaction_observables(
                aligned_positions,
                selections["interaction_pairs"],
                cutoff_angstrom=contact_cutoff_angstrom,
            )
            contact_samples.append(
                {
                    "frame": int(ts.frame),
                    "interaction": interaction_binary,
                    "interaction_distance": interaction_distance,
                }
            )

        records: list[dict[str, Any]] = []
        for ts in self.trajectory_universe.trajectory[::geometry_stride]:
            mobile_positions = self.trajectory_universe.atoms.positions.copy()
            aligned_positions = _align_positions(
                all_positions=mobile_positions,
                mobile_reference=mobile_positions[framework_indices],
                target_reference=ref_framework,
            )
            contact_feature, contact_distance = self._aggregate_contact_window(
                anchor_frame=int(ts.frame),
                contact_samples=contact_samples,
                contact_window_radius=contact_window_radius,
            )
            records.append(
                {
                    "frame": int(ts.frame),
                    "time_ps": float(ts.time),
                    "geometry": aligned_positions[geometry_indices].reshape(-1),
                    "aligned_positions": aligned_positions,
                    "sidechain": self._compute_sidechain_feature(
                        aligned_positions,
                        selections["sidechain_groups"],
                    ),
                    "interaction": contact_feature,
                    "interaction_distance": contact_distance,
                }
            )
        return records

    @staticmethod
    def _compute_sidechain_feature(
        positions: np.ndarray,
        sidechain_groups: list[tuple[str, np.ndarray, np.ndarray]],
    ) -> np.ndarray:
        feature: list[float] = []
        for _, ca_indices, sidechain_indices in sidechain_groups:
            ca_pos = positions[ca_indices].mean(axis=0)
            sidechain_com = positions[sidechain_indices].mean(axis=0)
            direction = sidechain_com - ca_pos
            norm = np.linalg.norm(direction)
            direction = direction / norm if norm > 1e-8 else np.zeros(3, dtype=float)
            feature.extend(direction.tolist())
        return np.array(feature, dtype=float)

    @staticmethod
    def _compute_interaction_observables(
        positions: np.ndarray,
        interaction_pairs: list[tuple[np.ndarray, np.ndarray]],
        cutoff_angstrom: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        binary = np.zeros(len(interaction_pairs), dtype=float)
        distances = np.zeros(len(interaction_pairs), dtype=float)
        for idx, (source_indices, target_indices) in enumerate(interaction_pairs):
            min_distance = cdist(positions[source_indices], positions[target_indices]).min()
            distances[idx] = float(min_distance)
            binary[idx] = 1.0 if min_distance <= cutoff_angstrom else 0.0
        return binary, distances

    @staticmethod
    def _aggregate_contact_window(
        anchor_frame: int,
        contact_samples: list[dict[str, Any]],
        contact_window_radius: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        window_binary = [
            sample["interaction"]
            for sample in contact_samples
            if abs(sample["frame"] - anchor_frame) <= contact_window_radius
        ]
        window_distances = [
            sample["interaction_distance"]
            for sample in contact_samples
            if abs(sample["frame"] - anchor_frame) <= contact_window_radius
        ]
        if not window_binary:
            nearest = min(contact_samples, key=lambda sample: abs(sample["frame"] - anchor_frame))
            return nearest["interaction"].copy(), nearest["interaction_distance"].copy()
        return np.mean(np.vstack(window_binary), axis=0), np.mean(np.vstack(window_distances), axis=0)

    @staticmethod
    def _cluster_frames(
        distance_matrix: np.ndarray,
        distance_cutoff: float,
        linkage_method: str,
    ) -> np.ndarray:
        linkage_matrix = linkage(squareform(distance_matrix, checks=False), method=linkage_method)
        raw_labels = fcluster(linkage_matrix, t=distance_cutoff, criterion="distance")
        mapping = {old: new for new, old in enumerate(sorted(set(raw_labels.tolist())))}
        labels = np.array([mapping[label] for label in raw_labels], dtype=int)
        return labels

    @staticmethod
    def _build_frame_assignments(frame_records: list[dict[str, Any]], labels: np.ndarray) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "frame": [record["frame"] for record in frame_records],
                "time_ps": [record["time_ps"] for record in frame_records],
                "cluster_id": labels,
            }
        )

    @staticmethod
    def _build_cluster_summary(
        frame_records: list[dict[str, Any]],
        frame_assignments: pd.DataFrame,
        joint_distance: np.ndarray,
        labels: np.ndarray,
    ) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        ordered = frame_assignments.sort_values("time_ps").reset_index(drop=True)
        total_frames = len(ordered)
        if total_frames > 1:
            sample_interval = float(np.median(np.diff(ordered["time_ps"].to_numpy(dtype=float))))
        else:
            sample_interval = 0.0
        segments = _build_cluster_segments(ordered)
        for cluster_id in sorted(frame_assignments["cluster_id"].unique().tolist()):
            member_indices = np.where(labels == cluster_id)[0]
            submatrix = joint_distance[np.ix_(member_indices, member_indices)]
            if len(member_indices) == 1:
                medoid_local_idx = 0
                mean_distance = 0.0
            else:
                mean_distances = submatrix.mean(axis=1)
                medoid_local_idx = int(np.argmin(mean_distances))
                mean_distance = float(mean_distances[medoid_local_idx])
            cluster_rows = frame_assignments.iloc[member_indices].sort_values("time_ps").reset_index(drop=True)
            cluster_segments = segments.get(int(cluster_id), [])
            representative_frame = int(cluster_rows.iloc[0]["frame"])
            mean_dwell_time = 0.0
            transitions_in = 0
            transitions_out = 0
            if cluster_segments:
                longest = max(cluster_segments, key=lambda item: item["length"])
                representative_frame = int(longest["representative_frame"])
                mean_dwell_time = float(np.mean([item["duration_ps"] for item in cluster_segments]))
                transitions_in = int(sum(1 for item in cluster_segments if item["start_index"] > 0))
                transitions_out = int(sum(1 for item in cluster_segments if item["end_index"] < total_frames - 1))
            rows.append(
                {
                    "cluster_id": int(cluster_id),
                    "population_percent": float(len(member_indices) / total_frames * 100.0),
                    "frames": int(len(member_indices)),
                    "representative_frame": int(representative_frame),
                    "medoid_frame": int(cluster_rows.iloc[medoid_local_idx]["frame"]),
                    "first_appearance_time_ps": float(cluster_rows["time_ps"].min()),
                    "mean_dwell_time_ps": mean_dwell_time if mean_dwell_time > 0 else sample_interval,
                    "transitions_in": transitions_in,
                    "transitions_out": transitions_out,
                    "main_structural_descriptor": "",
                }
            )
        columns = [
            "cluster_id",
            "population_percent",
            "frames",
            "representative_frame",
            "medoid_frame",
            "first_appearance_time_ps",
            "mean_dwell_time_ps",
            "transitions_in",
            "transitions_out",
            "main_structural_descriptor",
        ]
        return (
            pd.DataFrame(rows)[columns]
            .sort_values(["frames", "cluster_id"], ascending=[False, True])
            .reset_index(drop=True)
        )

    @staticmethod
    def _build_cluster_explanations(
        frame_records: list[dict[str, Any]],
        labels: np.ndarray,
        joint_distance: np.ndarray,
        interaction_labels: list[str],
        interaction_display_labels: list[str],
        interaction_descriptions: list[str],
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        interaction_matrix = np.vstack([record["interaction"] for record in frame_records])
        interaction_distance_matrix = np.vstack([record["interaction_distance"] for record in frame_records])
        cluster_ids = sorted(set(labels.tolist()))

        difference_rows: list[list[float]] = []
        for cluster_i in cluster_ids:
            row_values: list[float] = []
            members_i = np.where(labels == cluster_i)[0]
            for cluster_j in cluster_ids:
                members_j = np.where(labels == cluster_j)[0]
                if cluster_i == cluster_j:
                    row_values.append(0.0)
                else:
                    sub = joint_distance[np.ix_(members_i, members_j)]
                    row_values.append(float(sub.mean()) if sub.size else 0.0)
            difference_rows.append(row_values)
        cluster_difference_matrix = pd.DataFrame(
            difference_rows,
            index=[f"cluster_{cluster_id}" for cluster_id in cluster_ids],
            columns=[f"cluster_{cluster_id}" for cluster_id in cluster_ids],
        )

        feature_rows: list[dict[str, Any]] = []
        signature_rows: list[dict[str, Any]] = []
        digest_rows: list[dict[str, Any]] = []
        for cluster_id in cluster_ids:
            cluster_mask = labels == cluster_id
            cluster_matrix = interaction_matrix[cluster_mask]
            cluster_distance = interaction_distance_matrix[cluster_mask]
            other_matrix = interaction_matrix[~cluster_mask]
            occupancy = cluster_matrix.mean(axis=0)
            mean_distance = cluster_distance.mean(axis=0)
            occupancy_other = (
                other_matrix.mean(axis=0) if len(other_matrix) else np.zeros_like(occupancy)
            )
            delta = occupancy - occupancy_other
            feature_order = np.lexsort((mean_distance, -occupancy))
            for rank, idx in enumerate(feature_order[:15], start=1):
                feature_rows.append(
                    {
                        "cluster_id": int(cluster_id),
                        "pair_label": interaction_labels[int(idx)],
                        "pair_display": interaction_display_labels[int(idx)],
                        "pair_description": interaction_descriptions[int(idx)],
                        "contact_occupancy": float(occupancy[int(idx)]),
                        "mean_pair_distance_angstrom": float(mean_distance[int(idx)]),
                        "rank": int(rank),
                    }
                )
            signature_order = np.argsort(-(occupancy * np.clip(delta, a_min=0.0, a_max=None)))
            for idx in signature_order[:10]:
                signature_rows.append(
                    {
                        "cluster_id": int(cluster_id),
                        "pair_label": interaction_labels[int(idx)],
                        "pair_display": interaction_display_labels[int(idx)],
                        "pair_description": interaction_descriptions[int(idx)],
                        "contact_occupancy": float(occupancy[int(idx)]),
                        "mean_pair_distance_angstrom": float(mean_distance[int(idx)]),
                        "occupancy_delta_vs_others": float(delta[int(idx)]),
                        "signature_score": float(occupancy[int(idx)] * max(delta[int(idx)], 0.0)),
                    }
                )
            for rank, idx in enumerate(feature_order[:5], start=1):
                digest_rows.append(
                    {
                        "cluster_id": int(cluster_id),
                        "digest_type": "key_feature",
                        "pair_label": interaction_labels[int(idx)],
                        "pair_display": interaction_display_labels[int(idx)],
                        "pair_description": interaction_descriptions[int(idx)],
                        "contact_occupancy": float(occupancy[int(idx)]),
                        "mean_pair_distance_angstrom": float(mean_distance[int(idx)]),
                        "score": float(occupancy[int(idx)]),
                        "rank": int(rank),
                    }
                )
            positive_signature_indices = [
                int(idx)
                for idx in signature_order
                if float(occupancy[int(idx)] * max(delta[int(idx)], 0.0)) > 0.0
            ]
            for rank, idx in enumerate(positive_signature_indices[:5], start=1):
                digest_rows.append(
                    {
                        "cluster_id": int(cluster_id),
                        "digest_type": "state_signature",
                        "pair_label": interaction_labels[int(idx)],
                        "pair_display": interaction_display_labels[int(idx)],
                        "pair_description": interaction_descriptions[int(idx)],
                        "contact_occupancy": float(occupancy[int(idx)]),
                        "mean_pair_distance_angstrom": float(mean_distance[int(idx)]),
                        "score": float(occupancy[int(idx)] * max(delta[int(idx)], 0.0)),
                        "rank": int(rank),
                    }
                )
        return (
            cluster_difference_matrix,
            pd.DataFrame(feature_rows),
            pd.DataFrame(signature_rows),
            pd.DataFrame(digest_rows),
        )

    @staticmethod
    def _build_summary(
        frame_assignments: pd.DataFrame,
        cluster_summary: pd.DataFrame,
        stride: int,
        contact_stride: int,
        contact_cutoff_angstrom: float,
        distance_cutoff: float,
        linkage_method: str,
        geometry_weight: float,
        sidechain_weight: float,
        interaction_weight: float,
    ) -> dict[str, Any]:
        return {
            "n_frames": int(len(frame_assignments)),
            "n_clusters": int(frame_assignments["cluster_id"].nunique()),
            "stride": int(stride),
            "contact_stride": int(contact_stride),
            "contact_cutoff_angstrom": float(contact_cutoff_angstrom),
            "distance_cutoff": float(distance_cutoff),
            "linkage_method": linkage_method,
            "weights": {
                "geometry": float(geometry_weight),
                "sidechain": float(sidechain_weight),
                "interaction": float(interaction_weight),
            },
            "largest_cluster": cluster_summary.sort_values("frames", ascending=False).iloc[0].to_dict(),
        }

    @staticmethod
    def _attach_structural_descriptors(
        cluster_summary: pd.DataFrame,
        cluster_signatures: pd.DataFrame,
    ) -> pd.DataFrame:
        summary = cluster_summary.copy()
        descriptors: dict[int, str] = {}
        if not cluster_signatures.empty:
            for cluster_id, group in cluster_signatures.groupby("cluster_id"):
                top_pairs = group.sort_values(
                    ["signature_score", "contact_occupancy"],
                    ascending=[False, False],
                )["pair_display"].head(3).tolist()
                descriptors[int(cluster_id)] = " | ".join(top_pairs)
        summary["main_structural_descriptor"] = summary["cluster_id"].map(
            lambda cid: descriptors.get(int(cid), "No dominant signature")
        )
        return summary

    def export_cluster_structures(
        self,
        chain_mapping: dict[str, str],
        cdr_detection: dict[str, Any],
        frame_assignments: pd.DataFrame,
        cluster_summary: pd.DataFrame,
        stride: int,
        output_dir: str | Path,
    ) -> list[dict[str, Any]]:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        selections = self._prepare_selections(chain_mapping, cdr_detection)
        ref_positions = selections["reference_positions"]
        framework_indices = selections["framework_ca_indices"]
        ref_framework = ref_positions[framework_indices]
        complex_indices = selections["complex_atom_indices"]

        sampled_frames = {
            int(ts.frame): _align_positions(
                all_positions=self.trajectory_universe.atoms.positions.copy(),
                mobile_reference=self.trajectory_universe.atoms.positions.copy()[framework_indices],
                target_reference=ref_framework,
            )
            for ts in self.trajectory_universe.trajectory[::stride]
        }

        structure_rows: list[dict[str, Any]] = []
        for _, row in cluster_summary.iterrows():
            cluster_id = int(row["cluster_id"])
            cluster_rows = frame_assignments[frame_assignments["cluster_id"] == cluster_id].sort_values("time_ps")
            representative_frame = int(row["representative_frame"])
            medoid_frame = int(row["medoid_frame"])
            representative_file = output_path / f"cluster_{cluster_id}_representative.pdb"
            medoid_file = output_path / f"cluster_{cluster_id}_medoid.pdb"
            average_file = output_path / f"cluster_{cluster_id}_average.pdb"

            self._write_positions_to_pdb(sampled_frames[representative_frame], complex_indices, representative_file)
            self._write_positions_to_pdb(sampled_frames[medoid_frame], complex_indices, medoid_file)

            cluster_positions = [sampled_frames[int(frame)] for frame in cluster_rows["frame"].tolist()]
            average_positions = np.mean(np.stack(cluster_positions, axis=0), axis=0)
            self._write_positions_to_pdb(average_positions, complex_indices, average_file)

            structure_rows.append(
                {
                    "cluster_id": cluster_id,
                    "representative_pdb": str(representative_file),
                    "medoid_pdb": str(medoid_file),
                    "average_pdb": str(average_file),
                }
            )
        return structure_rows

    def _write_positions_to_pdb(
        self,
        positions: np.ndarray,
        atom_indices: np.ndarray,
        output_file: str | Path,
    ) -> None:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        original_positions = self.structure_universe.atoms.positions.copy()
        try:
            self.structure_universe.atoms.positions = positions
            selected_atoms = self.structure_universe.atoms[atom_indices]
            with mda.Writer(str(output_path), n_atoms=selected_atoms.n_atoms) as writer:
                writer.write(selected_atoms)
        finally:
            self.structure_universe.atoms.positions = original_positions


def _normalized_distance_matrix(features: np.ndarray, metric: str) -> np.ndarray:
    matrix = cdist(features, features, metric=metric)
    max_value = float(matrix.max())
    if max_value > 1e-12:
        matrix = matrix / max_value
    np.fill_diagonal(matrix, 0.0)
    return matrix


def _residue_chain_id(residue: Any) -> str:
    atoms = residue.atoms
    if hasattr(atoms, "chainIDs") and len(atoms.chainIDs):
        chain_id = str(atoms.chainIDs[0]).strip()
        if chain_id:
            return chain_id
    return str(getattr(residue.segment, "segid", "")).strip()


def _apply_protein_pbc_transformations(universe: Any) -> None:
    try:
        from MDAnalysis.transformations import center_in_box, unwrap

        protein = universe.select_atoms("protein")
        if protein.n_atoms == 0:
            return
        universe.trajectory.add_transformations(
            unwrap(protein),
            center_in_box(protein, center="geometry", wrap=False),
        )
    except Exception:
        return


def _derive_contact_stride(geometry_stride: int) -> int:
    if geometry_stride <= 1:
        return 1
    return int(max(1, math.ceil(geometry_stride / 2.0)))


def _build_cluster_segments(frame_assignments: pd.DataFrame) -> dict[int, list[dict[str, Any]]]:
    ordered = frame_assignments.sort_values("time_ps").reset_index(drop=True)
    if ordered.empty:
        return {}
    if len(ordered) > 1:
        sample_interval = float(np.median(np.diff(ordered["time_ps"].to_numpy(dtype=float))))
    else:
        sample_interval = 0.0

    segments: dict[int, list[dict[str, Any]]] = {}
    start_index = 0
    current_cluster = int(ordered.iloc[0]["cluster_id"])
    for idx in range(1, len(ordered) + 1):
        if idx == len(ordered) or int(ordered.iloc[idx]["cluster_id"]) != current_cluster:
            segment = ordered.iloc[start_index:idx]
            midpoint = int(len(segment) // 2)
            start_time = float(segment.iloc[0]["time_ps"])
            end_time = float(segment.iloc[-1]["time_ps"])
            duration = (end_time - start_time + sample_interval) if len(segment) > 0 else 0.0
            segments.setdefault(current_cluster, []).append(
                {
                    "start_index": int(start_index),
                    "end_index": int(idx - 1),
                    "length": int(len(segment)),
                    "start_time_ps": start_time,
                    "end_time_ps": end_time,
                    "duration_ps": float(duration),
                    "representative_frame": int(segment.iloc[midpoint]["frame"]),
                }
            )
            if idx < len(ordered):
                start_index = idx
                current_cluster = int(ordered.iloc[idx]["cluster_id"])
    return segments


def _align_positions(
    all_positions: np.ndarray,
    mobile_reference: np.ndarray,
    target_reference: np.ndarray,
) -> np.ndarray:
    mobile_center = mobile_reference.mean(axis=0)
    target_center = target_reference.mean(axis=0)
    mobile_centered = mobile_reference - mobile_center
    target_centered = target_reference - target_center
    covariance = mobile_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T
    return (all_positions - mobile_center) @ rotation + target_center


def write_cluster_id_vs_time(frame_assignments: pd.DataFrame, output_file: str | Path) -> None:
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    work = frame_assignments.sort_values("time_ps")
    cluster_ids = sorted(work["cluster_id"].unique().tolist())
    palette = plt.cm.Set2(np.linspace(0, 1, max(len(cluster_ids), 3)))
    fig, ax = plt.subplots(figsize=(12, 2.6))
    for idx, cluster_id in enumerate(cluster_ids):
        group = work[work["cluster_id"] == cluster_id]
        ax.scatter(group["time_ps"], np.full(len(group), cluster_id), s=18, color=palette[idx], label=f"Cluster {cluster_id}")
    ax.set_title("A. Cluster ID vs time", fontsize=12, fontweight="bold")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Cluster")
    ax.set_yticks(cluster_ids)
    ax.grid(alpha=0.15, axis="x")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_cluster_population_over_time(frame_assignments: pd.DataFrame, output_file: str | Path) -> None:
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    work = frame_assignments.sort_values("time_ps").reset_index(drop=True)
    if work.empty:
        return
    cluster_ids = sorted(work["cluster_id"].unique().tolist())
    palette = plt.cm.Set2(np.linspace(0, 1, max(len(cluster_ids), 3)))
    window = max(5, len(work) // 20)
    window = min(window, len(work))
    fig, axes = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    for idx, cluster_id in enumerate(cluster_ids):
        is_cluster = (work["cluster_id"] == cluster_id).astype(float)
        rolling = is_cluster.rolling(window=window, min_periods=1, center=True).mean()
        cumulative = is_cluster.expanding(min_periods=1).mean()
        axes[0].plot(work["time_ps"], rolling, color=palette[idx], linewidth=2, label=f"Cluster {cluster_id}")
        axes[1].plot(work["time_ps"], cumulative, color=palette[idx], linewidth=2, label=f"Cluster {cluster_id}")
    axes[0].set_title("B1. State population over time", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Rolling population")
    axes[0].grid(alpha=0.15)
    axes[1].set_title("B2. Cumulative population", fontsize=12, fontweight="bold")
    axes[1].set_xlabel("Time (ps)")
    axes[1].set_ylabel("Cumulative fraction")
    axes[1].grid(alpha=0.15)
    axes[0].legend(ncol=min(4, len(cluster_ids)), fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
