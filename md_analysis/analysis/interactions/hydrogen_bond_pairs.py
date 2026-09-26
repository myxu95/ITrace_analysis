"""
Hydrogen-bond residue-pair analyzer for pHLA-TCR complexes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.exceptions import NoDataError
from .occupancy_metrics import encode_segments, frame_set_to_segments, summarize_frames


@dataclass
class _ResiduePairAccumulator:
    chain_id_1: str
    resid_1: int
    resname_1: str
    residue_label_1: str
    chain_id_2: str
    resid_2: int
    resname_2: str
    residue_label_2: str
    min_distance_observed: float
    angle_sum: float
    event_count: int
    frames_seen: set[int]


class HydrogenBondPairAnalyzer:
    """Calculate residue-level hydrogen-bond occupancy across the pHLA-TCR interface."""

    def __init__(self, topology: str, trajectory: str, reference_structure: Optional[str] = None):
        self.universe = mda.Universe(topology, trajectory)
        self.reference_universe = mda.Universe(reference_structure) if reference_structure else None

    def calculate_interface_pairs(
        self,
        chain_mapping: Dict[str, str],
        stride: int = 1,
        distance_cutoff: float = 3.5,
        angle_cutoff: float = 150.0,
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
        all_selection = " or ".join(
            f"chainID {chain_id}" for chain_id in sorted(phla_chains | tcr_chains)
        )
        phla_selection = " or ".join(f"chainID {chain_id}" for chain_id in sorted(phla_chains))
        tcr_selection = " or ".join(f"chainID {chain_id}" for chain_id in sorted(tcr_chains))
        hydrogens_selection = f"({all_selection}) and name H*"

        analyzer = HBA(
            universe=self.universe,
            donors_sel=all_selection,
            hydrogens_sel=hydrogens_selection,
            acceptors_sel=all_selection,
            between=[phla_selection, tcr_selection],
            d_h_cutoff=1.2,
            d_a_cutoff=distance_cutoff,
            d_h_a_angle_cutoff=angle_cutoff,
        )
        analyzer.run(step=stride)

        hbonds = analyzer.results.hbonds if hasattr(analyzer, "results") else analyzer.hbonds
        total_frames = len(range(0, len(self.universe.trajectory), stride))
        if hbonds is None or len(hbonds) == 0:
            return self._empty_result()

        pair_map: dict[tuple, _ResiduePairAccumulator] = {}
        for row in hbonds:
            frame_index = int(row[0])
            donor_atom = self.universe.atoms[int(row[1])]
            acceptor_atom = self.universe.atoms[int(row[3])]
            distance = float(row[4])
            angle = float(row[5])

            donor_meta = self._get_residue_metadata(donor_atom)
            acceptor_meta = self._get_residue_metadata(acceptor_atom)

            donor_chain = donor_meta["chain_id"]
            acceptor_chain = acceptor_meta["chain_id"]
            donor_side = self._classify_side(donor_chain, phla_chains, tcr_chains)
            acceptor_side = self._classify_side(acceptor_chain, phla_chains, tcr_chains)
            if donor_side == acceptor_side or "unknown" in {donor_side, acceptor_side}:
                continue

            if donor_side == "phla":
                phla_meta = donor_meta
                tcr_meta = acceptor_meta
            else:
                phla_meta = acceptor_meta
                tcr_meta = donor_meta

            key = (
                phla_meta["chain_id"],
                int(phla_meta["resid"]),
                tcr_meta["chain_id"],
                int(tcr_meta["resid"]),
            )
            if key not in pair_map:
                pair_map[key] = _ResiduePairAccumulator(
                    chain_id_1=phla_meta["chain_id"],
                    resid_1=int(phla_meta["resid"]),
                    resname_1=phla_meta["resname"],
                    residue_label_1=f"{phla_meta['resname']}{int(phla_meta['resid'])}",
                    chain_id_2=tcr_meta["chain_id"],
                    resid_2=int(tcr_meta["resid"]),
                    resname_2=tcr_meta["resname"],
                    residue_label_2=f"{tcr_meta['resname']}{int(tcr_meta['resid'])}",
                    min_distance_observed=distance,
                    angle_sum=angle,
                    event_count=1,
                    frames_seen={frame_index},
                )
            else:
                acc = pair_map[key]
                acc.min_distance_observed = min(acc.min_distance_observed, distance)
                acc.angle_sum += angle
                acc.event_count += 1
                acc.frames_seen.add(frame_index)

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
                    "mean_angle_observed": float(acc.angle_sum / acc.event_count),
                    "frame_segments": encode_segments(segments),
                    "n_segments": int(segment_summary["n_segments"]),
                    "max_consecutive_frames": int(segment_summary["max_consecutive_frames"]),
                    "first_frame": segment_summary["first_frame"],
                    "last_frame": segment_summary["last_frame"],
                    "interaction_family": "hbond",
                    "distance_cutoff_angstrom": distance_cutoff,
                    "angle_cutoff_deg": angle_cutoff,
                }
            )

        if not rows:
            return self._empty_result()

        return pd.DataFrame(rows).sort_values(
            by=["contact_frequency", "min_distance_observed"],
            ascending=[False, True],
        ).reset_index(drop=True)

    @staticmethod
    def _classify_side(chain_id: str, phla_chains: set[str], tcr_chains: set[str]) -> str:
        if chain_id in phla_chains:
            return "phla"
        if chain_id in tcr_chains:
            return "tcr"
        return "unknown"

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
                "mean_angle_observed",
                "frame_segments",
                "n_segments",
                "max_consecutive_frames",
                "first_frame",
                "last_frame",
                "interaction_family",
                "distance_cutoff_angstrom",
                "angle_cutoff_deg",
            ]
        )

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
