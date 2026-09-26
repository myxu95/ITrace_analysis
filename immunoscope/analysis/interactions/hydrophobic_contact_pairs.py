"""
Hydrophobic-contact residue-pair analyzer for pHLA-TCR complexes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import MDAnalysis as mda
import pandas as pd
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib.distances import capped_distance
from .occupancy_metrics import encode_segments, frame_set_to_segments, summarize_frames


HYDROPHOBIC_RESNAMES = ("ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO", "CYS")


@dataclass
class _HydrophobicAccumulator:
    chain_id_1: str
    resid_1: int
    resname_1: str
    residue_label_1: str
    chain_id_2: str
    resid_2: int
    resname_2: str
    residue_label_2: str
    min_distance_observed: float
    frames_seen: set[int]


class HydrophobicContactPairAnalyzer:
    """Calculate residue-level hydrophobic-contact occupancy across the pHLA-TCR interface."""

    def __init__(self, topology: str, trajectory: str, reference_structure: Optional[str] = None):
        self.universe = mda.Universe(topology, trajectory)
        self.reference_universe = mda.Universe(reference_structure) if reference_structure else None

    def calculate_interface_pairs(
        self,
        chain_mapping: Dict[str, str],
        stride: int = 1,
        distance_cutoff: float = 4.5,
    ) -> pd.DataFrame:
        phla_chains = {
            chain_mapping["mhc_alpha"],
            chain_mapping["b2m"],
            chain_mapping["peptide"],
        }
        tcr_chains = {
            chain_mapping["tcr_alpha"],
            chain_mapping["tcr_beta"],
        }
        residue_clause = " or ".join(f"resname {name}" for name in HYDROPHOBIC_RESNAMES)
        sidechain_clause = "not (name N or name CA or name C or name O or name OXT or name H*)"
        phla_atoms = self.universe.select_atoms(
            f"({' or '.join(f'chainID {cid}' for cid in sorted(phla_chains))}) and ({residue_clause}) and {sidechain_clause}"
        )
        tcr_atoms = self.universe.select_atoms(
            f"({' or '.join(f'chainID {cid}' for cid in sorted(tcr_chains))}) and ({residue_clause}) and {sidechain_clause}"
        )

        total_frames = len(range(0, len(self.universe.trajectory), stride))
        if len(phla_atoms) == 0 or len(tcr_atoms) == 0:
            return self._empty_result()

        pair_map: dict[tuple, _HydrophobicAccumulator] = {}
        for ts in self.universe.trajectory[::stride]:
            frame_index = int(ts.frame)
            pairs, distances = capped_distance(
                phla_atoms.positions,
                tcr_atoms.positions,
                max_cutoff=distance_cutoff,
                return_distances=True,
                box=ts.dimensions,
            )
            for (idx_a, idx_b), distance in zip(pairs, distances):
                atom_a = phla_atoms[int(idx_a)]
                atom_b = tcr_atoms[int(idx_b)]
                phla_meta = self._get_residue_metadata(atom_a)
                tcr_meta = self._get_residue_metadata(atom_b)
                self._accumulate_pair(pair_map, phla_meta, tcr_meta, frame_index, float(distance))

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
                    "frame_segments": encode_segments(segments),
                    "n_segments": int(segment_summary["n_segments"]),
                    "max_consecutive_frames": int(segment_summary["max_consecutive_frames"]),
                    "first_frame": segment_summary["first_frame"],
                    "last_frame": segment_summary["last_frame"],
                    "interaction_family": "hydrophobic_contact",
                    "distance_cutoff_angstrom": distance_cutoff,
                }
            )

        return pd.DataFrame(rows).sort_values(
            by=["contact_frequency", "min_distance_observed"],
            ascending=[False, True],
        ).reset_index(drop=True)

    def _accumulate_pair(self, pair_map, phla_meta, tcr_meta, frame_index: int, distance: float):
        key = (
            phla_meta["chain_id"],
            int(phla_meta["resid"]),
            tcr_meta["chain_id"],
            int(tcr_meta["resid"]),
        )
        if key not in pair_map:
            pair_map[key] = _HydrophobicAccumulator(
                chain_id_1=phla_meta["chain_id"],
                resid_1=int(phla_meta["resid"]),
                resname_1=phla_meta["resname"],
                residue_label_1=f"{phla_meta['resname']}{int(phla_meta['resid'])}",
                chain_id_2=tcr_meta["chain_id"],
                resid_2=int(tcr_meta["resid"]),
                resname_2=tcr_meta["resname"],
                residue_label_2=f"{tcr_meta['resname']}{int(tcr_meta['resid'])}",
                min_distance_observed=distance,
                frames_seen={frame_index},
            )
        else:
            acc = pair_map[key]
            acc.min_distance_observed = min(acc.min_distance_observed, distance)
            acc.frames_seen.add(frame_index)

    def _get_residue_metadata(self, atom) -> dict:
        residue = atom.residue
        try:
            chain_id = atom.chainID
        except NoDataError:
            if self.reference_universe is None:
                raise
            ref_atom = self.reference_universe.atoms[atom.index]
            ref_residue = ref_atom.residue
            chain_id = ref_atom.chainID
            resid = int(ref_residue.resid)
            resname = str(ref_residue.resname)
            return {"chain_id": chain_id, "resid": resid, "resname": resname}

        return {
            "chain_id": chain_id,
            "resid": int(residue.resid),
            "resname": str(residue.resname),
        }

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
                "frame_segments",
                "n_segments",
                "max_consecutive_frames",
                "first_frame",
                "last_frame",
                "interaction_family",
                "distance_cutoff_angstrom",
            ]
        )
