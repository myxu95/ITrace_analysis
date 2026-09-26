"""Backbone φ/ψ (Ramachandran) per-residue analysis.

Computes, per residue, the time-averaged φ/ψ, the populations of the four
Ramachandran regions (alpha / beta / left_alpha / other), and an angular
spread (flexibility proxy). Mirrors the logic in the `ims run` dihedrals
module (`cli/commands/run.py:_run_module_dihedrals`) but lives here as a
reusable function so it can be called both by the pipeline and to backfill
older analysis runs that predate the `dihedrals` module.

Outputs (written under `<case>/analysis/dihedrals/`):
  residue_dihedrals.csv     — per-residue rows
  dihedrals_summary.json    — region_summary + most_flexible/rigid + means
                              (the file the `dihedrals` view reads)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

_REGIONS = ("alpha", "beta", "left_alpha", "other")


def _classify(phi: float, psi: float) -> str:
    """Simplified Ramachandran region classifier (Lovell et al. 2003).

    `phi`/`psi` are in DEGREES (MDAnalysis Ramachandran returns degrees).
    """
    if -180 <= phi <= -30 and -70 <= psi <= 50:
        return "alpha"
    if -180 <= phi <= -40 and 90 <= psi <= 180:
        return "beta"
    if 30 <= phi <= 90 and -20 <= psi <= 100:
        return "left_alpha"
    return "other"


def _circ_mean_deg(deg: np.ndarray) -> float:
    """Circular mean of an angle array (degrees), NaN-aware."""
    d = deg[~np.isnan(deg)]
    if d.size == 0:
        return float("nan")
    r = np.radians(d)
    return float(np.degrees(np.arctan2(np.mean(np.sin(r)), np.mean(np.cos(r)))))


def _circ_std_deg(deg: np.ndarray) -> float:
    """Circular standard deviation (degrees), NaN-aware.

    Bounded and wrap-aware: a residue oscillating across the ±180° seam is
    NOT reported as hugely variable, unlike a naive linear std.
    """
    d = deg[~np.isnan(deg)]
    if d.size == 0:
        return 0.0
    r = np.radians(d)
    resultant = np.hypot(np.mean(np.cos(r)), np.mean(np.sin(r)))
    resultant = min(1.0, max(1e-12, float(resultant)))
    return float(np.degrees(np.sqrt(-2.0 * np.log(resultant))))


def compute_backbone_dihedrals(
    topology: str, trajectory: str, stride: int = 1
) -> Tuple[List[Dict], Dict]:
    """Run Ramachandran on a trajectory; return (per_residue_rows, summary)."""
    import MDAnalysis as mda
    from MDAnalysis.analysis.dihedrals import Ramachandran

    u = mda.Universe(topology, trajectory)
    protein = u.select_atoms("protein and name CA")
    rama = Ramachandran(protein).run(step=stride)
    angles = rama.results.angles  # (n_frames, n_residues, 2) -> (phi, psi)
    n_frames = int(angles.shape[0])
    # `rama.ag1.residues` are exactly the residues for which φ/ψ were computed
    # (chain termini are dropped), so it aligns 1:1 with the angle columns.
    # `rama.atomgroup.residues` includes the dropped termini and is off by the
    # number of chains — using it mis-maps residues and IndexErrors.
    ag1 = getattr(rama, "ag1", None)
    residues = list(ag1.residues) if ag1 is not None else list(rama.atomgroup.residues)
    n_cols = int(angles.shape[1])
    if len(residues) != n_cols:
        residues = residues[:n_cols]

    rows: List[Dict] = []
    for ri, residue in enumerate(residues):
        phi_arr = angles[:, ri, 0]
        psi_arr = angles[:, ri, 1]
        if np.all(np.isnan(phi_arr)) or np.all(np.isnan(psi_arr)):
            continue
        regions = [
            _classify(p, q)
            for p, q in zip(phi_arr, psi_arr)
            if not (np.isnan(p) or np.isnan(q))
        ]
        n = len(regions) or 1
        counts = {r: regions.count(r) / n for r in _REGIONS}
        chain = (
            getattr(residue, "chainID", None)
            or getattr(residue, "segid", None)
            or "-"
        )
        phi_std = _circ_std_deg(phi_arr)
        psi_std = _circ_std_deg(psi_arr)
        rows.append(
            {
                "chain_id": chain,
                "resid": int(residue.resid),
                "resname": residue.resname,
                "phi_mean_deg": _circ_mean_deg(phi_arr),
                "phi_std_deg": phi_std,
                "psi_mean_deg": _circ_mean_deg(psi_arr),
                "psi_std_deg": psi_std,
                "frac_alpha": float(counts["alpha"]),
                "frac_beta": float(counts["beta"]),
                "frac_left_alpha": float(counts["left_alpha"]),
                "frac_other": float(counts["other"]),
                "dominant_region": max(counts, key=counts.get),
                # Flexibility proxy: mean circular spread of φ and ψ (degrees).
                "angular_spread_deg": float(0.5 * (phi_std + psi_std)),
            }
        )

    summary: Dict = {
        "n_frames": n_frames,
        "n_residues_analyzed": len(rows),
        "stride": stride,
        "region_summary": {},
    }
    if rows:
        for region in _REGIONS:
            fracs = [r[f"frac_{region}"] for r in rows]
            summary["region_summary"][region] = {
                "mean_fraction": float(np.mean(fracs)),
                "n_residues_dominant": sum(
                    1 for r in rows if r["dominant_region"] == region
                ),
            }
        spreads = [r["angular_spread_deg"] for r in rows]
        summary["most_flexible_residues"] = sorted(
            rows, key=lambda r: -r["angular_spread_deg"]
        )[:10]
        summary["most_rigid_residues"] = sorted(
            rows, key=lambda r: r["angular_spread_deg"]
        )[:10]
        summary["mean_angular_spread_deg"] = float(np.mean(spreads))
    return rows, summary


def write_backbone_dihedrals(out_dir: Path, rows: List[Dict], summary: Dict) -> Path:
    import csv

    out_dir.mkdir(parents=True, exist_ok=True)
    if rows:
        with (out_dir / "residue_dihedrals.csv").open(
            "w", newline="", encoding="utf-8"
        ) as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    summary_path = out_dir / "dihedrals_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary_path


# MATCHED (topology, trajectory) pairs — these must have the same atom count,
# so they are resolved as pairs, never independently. Ordered by preference:
# the stripped/aligned analysis input first (what the analysis modules used).
_MATCHED_PAIRS = (
    ("preparation/analysis_input/analysis_structure.pdb",
     "preparation/analysis_input/analysis_trajectory.xtc"),
    ("preparation/processed_trajectory_converted.pdb",
     "preparation/processed_trajectory.xtc"),
    ("01_preprocess/md_processed_converted.pdb",
     "01_preprocess/md_processed.xtc"),
)


def _resolve_topology_trajectory(case: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Resolve a matched (topology, trajectory) pair for the case.

    Prefers the authoritative analysis pair recorded in run_manifest.json
    (`analysis_topology`/`analysis_trajectory`); falls back to the matched
    on-disk pairs above. Topology and trajectory MUST come from the same
    source — pairing mismatched files yields an atom-count error.
    """
    manifest = case / "run_manifest.json"
    if manifest.exists():
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            topo, traj = _scan_manifest_pair(data)
            if topo and traj and (case_root(topo)).exists() and (case_root(traj)).exists():
                return Path(topo), Path(traj)
        except Exception:  # noqa: BLE001 — fall through to disk pairs
            pass
    for topo_rel, traj_rel in _MATCHED_PAIRS:
        topo, traj = case / topo_rel, case / traj_rel
        if topo.exists() and traj.exists():
            return topo, traj
    return None, None


def case_root(p: str) -> Path:
    return Path(p)


def _scan_manifest_pair(data) -> Tuple[Optional[str], Optional[str]]:
    """Pull a matched analysis topology+trajectory out of a manifest dict."""
    found: Dict[str, str] = {}

    def walk(d):
        if isinstance(d, dict):
            for k, v in d.items():
                kl = k.lower()
                if isinstance(v, str):
                    if kl == "analysis_topology":
                        found.setdefault("topo", v)
                    elif kl == "analysis_trajectory":
                        found.setdefault("traj", v)
                walk(v)
        elif isinstance(d, list):
            for x in d:
                walk(x)

    walk(data)
    return found.get("topo"), found.get("traj")


def backfill_case(case_dir: str, stride: int = 1) -> Dict:
    """Compute + write backbone dihedrals for an existing analysis case dir."""
    case = Path(case_dir)
    topo, traj = _resolve_topology_trajectory(case)
    if not topo or not traj:
        return {"status": "failed", "error": "matched topology/trajectory not found"}
    rows, summary = compute_backbone_dihedrals(str(topo), str(traj), stride=stride)
    out_dir = case / "analysis" / "dihedrals"
    path = write_backbone_dihedrals(out_dir, rows, summary)
    return {
        "status": "completed",
        "summary_json": str(path),
        "n_residues": len(rows),
        "n_frames": summary.get("n_frames", 0),
        "topology": str(topo),
        "trajectory": str(traj),
    }


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Backfill backbone φ/ψ dihedrals for a case")
    ap.add_argument("case_dir")
    ap.add_argument("--stride", type=int, default=1)
    res = backfill_case(ap.parse_args().case_dir, stride=ap.parse_args().stride)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()


__all__ = [
    "compute_backbone_dihedrals",
    "write_backbone_dihedrals",
    "backfill_case",
]
