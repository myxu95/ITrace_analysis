"""
Region-aware residue-level RMSF analysis module.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis.analysis.rms import RMSF

from ..topology.complex_residue_semantics import ComplexResidueSemanticAnnotator


TCR_REGION_ORDER = [
    "CDR1_alpha",
    "CDR2_alpha",
    "CDR3_alpha",
    "non_cdr_alpha",
    "CDR1_beta",
    "CDR2_beta",
    "CDR3_beta",
    "non_cdr_beta",
]

PHLA_REGION_ORDER = [
    "peptide",
    "alpha1_helix",
    "alpha2_helix",
    "non_groove",
    "beta2m",
]


@dataclass
class ResidueRMSFResult:
    """Region-aware RMSF result container."""

    residue_frame: pd.DataFrame
    region_summary: pd.DataFrame
    summary: dict
    stride: int
    time_unit: str


@dataclass
class ResidueRMSFAtomSetsResult:
    """Per-residue backbone / sidechain RMSF result.

    Sidechain RMSF is NaN for GLY (no sidechain heavy atoms). rmsf_ratio is
    sidechain / backbone (NaN where either is missing or backbone == 0).
    """

    residue_frame: pd.DataFrame
    summary: dict
    stride: int


class ResidueRMSFAnalyzer:
    """Calculate and annotate pHLA-TCR residue-level RMSF."""

    def __init__(
        self,
        topology_file: str,
        trajectory_file: str,
        structure_pdb: str,
    ):
        self.topology_file = str(topology_file)
        self.trajectory_file = str(trajectory_file)
        self.structure_pdb = str(structure_pdb)
        self.trajectory_universe = mda.Universe(self.topology_file, self.trajectory_file)
        self.structure_universe = mda.Universe(self.structure_pdb)

    def calculate(
        self,
        chain_mapping: dict,
        cdr_detection: Optional[dict] = None,
        stride: int = 1,
        selection: Optional[str] = None,
        time_unit: str = "ps",
    ) -> ResidueRMSFResult:
        """Calculate residue-level RMSF and add region annotations."""
        if stride < 1:
            raise ValueError("stride must be >= 1")

        selected_chain_ids = [chain_id for chain_id in chain_mapping.values() if chain_id]
        if selection is None:
            selected_chain_clause = " or ".join(f"chainID {chain_id}" for chain_id in selected_chain_ids)
            selection = f"name CA and ({selected_chain_clause})"

        structure_atoms = self.structure_universe.select_atoms(selection)
        if len(structure_atoms) == 0:
            raise ValueError(f"No RMSF selection atoms found: {selection}")

        trajectory_atoms = self.trajectory_universe.atoms[structure_atoms.indices]
        if len(trajectory_atoms) != len(structure_atoms):
            raise ValueError("Structure reference and trajectory topology atom indices do not match; RMSF cannot be mapped")

        rmsf_values = RMSF(trajectory_atoms).run(step=stride).results.rmsf
        if len(rmsf_values) != len(structure_atoms):
            raise RuntimeError("RMSF result length does not match atom count")

        annotator = ComplexResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

        records = []
        for atom, rmsf in zip(structure_atoms, rmsf_values):
            semantics = annotator.annotate_residue(atom.chainID, int(atom.resid))
            region_group = self._resolve_region_group(semantics)
            records.append(
                {
                    "chain_id": atom.chainID,
                    "resid": int(atom.resid),
                    "resname": str(atom.resname),
                    "atom_name": str(atom.name),
                    "component": semantics.component,
                    "complex_side": semantics.complex_side,
                    "phla_region": semantics.phla_region,
                    "mhc_region": semantics.mhc_region,
                    "mhc_subregion": semantics.mhc_subregion,
                    "tcr_chain": semantics.tcr_chain,
                    "tcr_region": semantics.tcr_region,
                    "tcr_region_detailed": semantics.tcr_region_detailed,
                    "region_group": region_group,
                    "rmsf_angstrom": float(rmsf),
                }
            )

        residue_frame = pd.DataFrame.from_records(records)
        region_summary = self._summarize_regions(residue_frame)

        n_frames = len(self.trajectory_universe.trajectory[::stride])
        start_time = float(self.trajectory_universe.trajectory[0].time)
        end_time = float(self.trajectory_universe.trajectory[-1].time)
        time_divisor = 1000.0 if time_unit == "ns" else 1.0
        summary = {
            "selection": selection,
            "stride": int(stride),
            "n_residues": int(len(residue_frame)),
            "n_frames": int(n_frames),
            "time_unit": time_unit,
            "time_start": start_time / time_divisor,
            "time_end": end_time / time_divisor,
            "mean_rmsf_angstrom": float(residue_frame["rmsf_angstrom"].mean()),
            "max_rmsf_angstrom": float(residue_frame["rmsf_angstrom"].max()),
            "tcr_mean_rmsf_angstrom": float(
                residue_frame.loc[residue_frame["component"].isin(["TCR_alpha", "TCR_beta"]), "rmsf_angstrom"].mean()
            ),
            "phla_mean_rmsf_angstrom": float(
                residue_frame.loc[residue_frame["component"].isin(["HLA_alpha", "peptide"]), "rmsf_angstrom"].mean()
            ),
        }

        return ResidueRMSFResult(
            residue_frame=residue_frame,
            region_summary=region_summary,
            summary=summary,
            stride=stride,
            time_unit=time_unit,
        )

    def calculate_atom_sets(
        self,
        chain_mapping: dict,
        cdr_detection: Optional[dict] = None,
        stride: int = 1,
    ) -> ResidueRMSFAtomSetsResult:
        """Compute per-residue RMSF split into backbone and sidechain.

        Backbone = {N, CA, C, O}. Sidechain = protein heavy atoms not in backbone
        (hydrogens excluded). Per-residue value is the mean over the residue's
        atoms in that set. GLY sidechain → NaN. rmsf_ratio = sidechain / backbone.
        """
        if stride < 1:
            raise ValueError("stride must be >= 1")

        selected_chain_ids = [chain_id for chain_id in chain_mapping.values() if chain_id]
        if not selected_chain_ids:
            raise ValueError("chain_mapping has no chain IDs to analyse")
        chain_clause = " or ".join(f"chainID {chain_id}" for chain_id in selected_chain_ids)

        backbone_sel = f"(name N or name CA or name C or name O) and ({chain_clause})"
        sidechain_sel = (
            f"protein and not (name N or name CA or name C or name O) "
            f"and not (name H* or type H) and ({chain_clause})"
        )

        backbone_per_residue = self._atom_set_rmsf(backbone_sel, stride=stride)
        sidechain_per_residue = self._atom_set_rmsf(sidechain_sel, stride=stride)

        annotator = ComplexResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

        records = []
        seen_keys = set(backbone_per_residue.keys()) | set(sidechain_per_residue.keys())
        for key in seen_keys:
            chain_id, resid, resname = key
            bb = backbone_per_residue.get(key)
            sc = sidechain_per_residue.get(key)
            bb_mean = bb["mean"] if bb else float("nan")
            bb_n = bb["n_atoms"] if bb else 0
            sc_mean = sc["mean"] if sc else float("nan")
            sc_n = sc["n_atoms"] if sc else 0
            if np.isnan(bb_mean) or np.isnan(sc_mean) or bb_mean == 0.0:
                ratio = float("nan")
            else:
                ratio = float(sc_mean / bb_mean)
            semantics = annotator.annotate_residue(chain_id, resid)
            region_group = self._resolve_region_group(semantics)
            records.append(
                {
                    "chain_id": chain_id,
                    "resid": resid,
                    "resname": resname,
                    "component": semantics.component,
                    "complex_side": semantics.complex_side,
                    "tcr_chain": semantics.tcr_chain,
                    "tcr_region": semantics.tcr_region,
                    "tcr_region_detailed": semantics.tcr_region_detailed,
                    "region_group": region_group,
                    "rmsf_backbone": float(bb_mean),
                    "rmsf_sidechain": float(sc_mean),
                    "rmsf_ratio": float(ratio),
                    "n_backbone_atoms": int(bb_n),
                    "n_sidechain_atoms": int(sc_n),
                }
            )

        residue_frame = pd.DataFrame.from_records(records).sort_values(
            ["chain_id", "resid"]
        ).reset_index(drop=True)

        n_frames = len(self.trajectory_universe.trajectory[::stride])
        bb_series = residue_frame["rmsf_backbone"]
        sc_series = residue_frame["rmsf_sidechain"]
        ratio_series = residue_frame["rmsf_ratio"]
        summary = {
            "backbone_selection": backbone_sel,
            "sidechain_selection": sidechain_sel,
            "stride": int(stride),
            "n_residues": int(len(residue_frame)),
            "n_frames": int(n_frames),
            "mean_rmsf_backbone": float(bb_series.mean(skipna=True)),
            "mean_rmsf_sidechain": float(sc_series.mean(skipna=True)),
            "median_rmsf_ratio": float(ratio_series.median(skipna=True)),
            "n_rotamer_tolerant": int((ratio_series >= 2.0).sum()),
            "n_anchored": int((ratio_series < 1.0).sum()),
        }
        return ResidueRMSFAtomSetsResult(
            residue_frame=residue_frame,
            summary=summary,
            stride=stride,
        )

    def _atom_set_rmsf(self, selection: str, stride: int) -> dict:
        """Run RMSF on a selection then aggregate per residue.

        Returns dict keyed by (chain_id, resid, resname) → {mean, n_atoms}.
        """
        structure_atoms = self.structure_universe.select_atoms(selection)
        if len(structure_atoms) == 0:
            return {}

        trajectory_atoms = self.trajectory_universe.atoms[structure_atoms.indices]
        if len(trajectory_atoms) != len(structure_atoms):
            raise ValueError("Structure and trajectory atom indices do not align")

        rmsf_values = RMSF(trajectory_atoms).run(step=stride).results.rmsf
        if len(rmsf_values) != len(structure_atoms):
            raise RuntimeError("RMSF result length does not match selection size")

        per_residue: dict = {}
        for atom, value in zip(structure_atoms, rmsf_values):
            key = (str(atom.chainID), int(atom.resid), str(atom.resname))
            bucket = per_residue.setdefault(key, {"sum": 0.0, "n": 0})
            bucket["sum"] += float(value)
            bucket["n"] += 1
        return {
            key: {"mean": b["sum"] / b["n"], "n_atoms": b["n"]}
            for key, b in per_residue.items()
        }

    @staticmethod
    def _resolve_region_group(semantics) -> str:
        if semantics.component in {"TCR_alpha", "TCR_beta"}:
            if semantics.tcr_region_detailed:
                return semantics.tcr_region_detailed
            if semantics.tcr_chain:
                return f"non_cdr_{semantics.tcr_chain}"
            return "tcr_other"
        if semantics.component == "peptide":
            return "peptide"
        if semantics.component == "HLA_alpha":
            return semantics.mhc_subregion or "non_groove"
        if semantics.component == "beta2m":
            return "beta2m"
        return "unknown"

    @staticmethod
    def _summarize_regions(residue_frame: pd.DataFrame) -> pd.DataFrame:
        summary = (
            residue_frame.groupby(["region_group"], as_index=False)
            .agg(
                n_residues=("resid", "count"),
                mean_rmsf_angstrom=("rmsf_angstrom", "mean"),
                median_rmsf_angstrom=("rmsf_angstrom", "median"),
                max_rmsf_angstrom=("rmsf_angstrom", "max"),
                min_rmsf_angstrom=("rmsf_angstrom", "min"),
            )
        )
        order = {name: idx for idx, name in enumerate(TCR_REGION_ORDER + PHLA_REGION_ORDER)}
        summary["sort_order"] = summary["region_group"].map(lambda value: order.get(value, 999))
        summary = summary.sort_values(["sort_order", "region_group"]).drop(columns=["sort_order"]).reset_index(drop=True)
        return summary


def write_tcr_rmsf_profile(residue_frame: pd.DataFrame, output_file: Path) -> None:
    """Plot TCR alpha/beta residue RMSF profiles."""
    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    plotted = False
    colors = {"alpha": "#2f5d50", "beta": "#c37a3b"}
    for chain in ["alpha", "beta"]:
        chain_df = residue_frame[residue_frame["tcr_chain"] == chain].copy()
        if chain_df.empty:
            continue
        chain_df = chain_df.sort_values("resid")
        ax.plot(
            chain_df["resid"],
            chain_df["rmsf_angstrom"],
            linewidth=1.8,
            color=colors[chain],
            label=f"TCR {chain}",
        )
        plotted = True
    if not plotted:
        ax.text(0.5, 0.5, "No TCR RMSF profile", ha="center", va="center", fontsize=13)
        ax.axis("off")
    else:
        ax.set_title("TCR RMSF profile", fontsize=14, weight="bold")
        ax.set_xlabel("Residue ID")
        ax.set_ylabel("RMSF (Å)")
        ax.grid(alpha=0.18, linestyle="--")
        ax.legend(frameon=False)
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_phla_rmsf_profile(residue_frame: pd.DataFrame, output_file: Path) -> None:
    """Plot peptide and HLA alpha RMSF profiles."""
    fig, axes = plt.subplots(2, 1, figsize=(9.8, 6.4), sharex=False)
    mhc_df = residue_frame[residue_frame["component"] == "HLA_alpha"].copy().sort_values("resid")
    peptide_df = residue_frame[residue_frame["component"] == "peptide"].copy().sort_values("resid")

    if mhc_df.empty:
        axes[0].text(0.5, 0.5, "No HLA alpha RMSF profile", ha="center", va="center", fontsize=12)
        axes[0].axis("off")
    else:
        axes[0].plot(mhc_df["resid"], mhc_df["rmsf_angstrom"], color="#658c7c", linewidth=1.8)
        axes[0].set_title("MHC alpha RMSF profile", fontsize=13, weight="bold")
        axes[0].set_ylabel("RMSF (Å)")
        axes[0].grid(alpha=0.18, linestyle="--")
        axes[0].set_axisbelow(True)

    if peptide_df.empty:
        axes[1].text(0.5, 0.5, "No peptide RMSF profile", ha="center", va="center", fontsize=12)
        axes[1].axis("off")
    else:
        axes[1].bar(peptide_df["resid"].astype(str), peptide_df["rmsf_angstrom"], color="#e0a34d")
        axes[1].set_title("Peptide RMSF profile", fontsize=13, weight="bold")
        axes[1].set_ylabel("RMSF (Å)")
        axes[1].set_xlabel("Residue ID")
        axes[1].grid(axis="y", alpha=0.18, linestyle="--")
        axes[1].set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_region_rmsf_summary(region_summary: pd.DataFrame, output_file: Path) -> None:
    """Plot a region-level RMSF summary bar chart."""
    fig, ax = plt.subplots(figsize=(11.2, 5.2))
    if region_summary.empty:
        ax.text(0.5, 0.5, "No region RMSF summary", ha="center", va="center", fontsize=13)
        ax.axis("off")
    else:
        top = region_summary.copy()
        ax.bar(
            top["region_group"],
            top["mean_rmsf_angstrom"],
            color=["#2f5d50" if value in TCR_REGION_ORDER else "#c37a3b" for value in top["region_group"]],
        )
        ax.set_title("Region-level RMSF summary", fontsize=14, weight="bold")
        ax.set_ylabel("Mean RMSF (Å)")
        ax.grid(axis="y", alpha=0.18, linestyle="--")
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)
