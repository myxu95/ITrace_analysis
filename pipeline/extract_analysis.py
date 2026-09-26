"""Curate the heavy per-system analysis output into compact web assets.

The batch run (``pipeline.run_analysis``) writes a large tree of CSV/JSON/PNG
files per trajectory under an analysis output root. For the website we only need
a small, flat summary plus a handful of figures, so this script:

  1. reads the per-stage JSON summaries for each trajectory,
  2. distils them into a single compact ``analysis.json``, and
  3. copies a curated set of figures into ``web_data/<id>/analysis/``.

Both products live under ``web_data/<id>/analysis/`` and are served as static
assets by the backend (``/data/<id>/analysis/...``).

Usage::

    python -m pipeline.extract_analysis --source /home/xmy/immuno_analysis
    python -m pipeline.extract_analysis --source /home/xmy/immuno_analysis --ids 1ao7_run2
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
from pathlib import Path

from pipeline import config
from pipeline import cdr_contacts
from pipeline import interface_metrics
from pipeline import struct_metrics

# Identifies the analysis toolkit that produced the per-stage data, so the
# served analysis.json records which code computed it. The toolkit is
# archived in this repository as `md_analysis/`; no independently verifiable
# commit hash is recorded for the code that produced the published dataset.
ANALYSIS_PROVENANCE = {
    "tool": "md_analysis",
    "version": "0.1.0",
}

# Per-stage JSON summary locations, relative to <source>/<stage>/<traj_id>/.
JSON_PATHS = {
    "identity": "analysis/identity/biological_identity.json",
    "rmsf": "analysis/rmsf/rmsf_summary.json",
    "bsa": "analysis/interface/interface_summary.json",
    "contact_occupancy": "analysis/contacts/occupancy/occupancy_summary.json",
    "contact_annotation": "analysis/contacts/contact_annotation_summary.json",
    "hbond": "analysis/interactions/hydrogen_bonds/hbond_annotation_summary.json",
    "hbond_occupancy": "analysis/interactions/hydrogen_bonds/occupancy/occupancy_summary.json",
    "saltbridge": "analysis/interactions/salt_bridges/salt_bridge_annotation_summary.json",
    "saltbridge_occupancy": "analysis/interactions/salt_bridges/occupancy/occupancy_summary.json",
    "hydrophobic": "analysis/interactions/hydrophobic_contacts/hydrophobic_annotation_summary.json",
    "hydrophobic_occupancy": "analysis/interactions/hydrophobic_contacts/occupancy/occupancy_summary.json",
    "pipi": "analysis/interactions/pi_interactions/pi_pi_annotation_summary.json",
    "pipi_occupancy": "analysis/interactions/pi_interactions/occupancy/occupancy_summary.json",
    "cationpi": "analysis/interactions/cation_pi_interactions/cation_pi_annotation_summary.json",
    "cationpi_occupancy": "analysis/interactions/cation_pi_interactions/occupancy/occupancy_summary.json",
}

# Curated figures: (stage, source-relative-path) -> destination file name.
FIGURES = {
    ("rmsf", "analysis/rmsf/region_rmsf_summary.png"): "rmsf_region.png",
    ("rmsf", "analysis/rmsf/phla_rmsf_profile.png"): "rmsf_phla.png",
    ("rmsf", "analysis/rmsf/tcr_rmsf_profile.png"): "rmsf_tcr.png",
    ("bsa", "analysis/interface/bsa_timeseries.png"): "bsa_timeseries.png",
    ("contact", "analysis/contacts/occupancy/contact_persistent_timeline.png"): "contact_timeline.png",
    ("contact", "analysis/contacts/pep_TCRa_heatmap.png"): "contact_pep_tcra.png",
    ("contact", "analysis/contacts/pep_TCRb_heatmap.png"): "contact_pep_tcrb.png",
}

# Interaction stages keyed by their compact name -> (annotation key, occupancy key).
INTERACTIONS = {
    "hbond": ("hbond", "hbond_occupancy"),
    "saltbridge": ("saltbridge", "saltbridge_occupancy"),
    "hydrophobic": ("hydrophobic", "hydrophobic_occupancy"),
    "pipi": ("pipi", "pipi_occupancy"),
    "cationpi": ("cationpi", "cationpi_occupancy"),
}


def load_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def stat_block(d: dict | None, keys=("mean", "std", "min", "max")) -> dict | None:
    if not isinstance(d, dict):
        return None
    return {k: d.get(k) for k in keys}


def angle_stats(csv_path: Path) -> dict | None:
    if not csv_path.exists():
        return None
    times, crossing, incident = [], [], []
    try:
        with csv_path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    c = float(row["Crossing(deg)"])
                    i = float(row["Incident(deg)"])
                except (KeyError, ValueError):
                    continue
                try:
                    t = float(row.get("Time(ps)", "nan"))
                except ValueError:
                    t = float("nan")
                crossing.append(c)
                incident.append(i)
                times.append(t)
    except OSError:
        return None
    if not crossing:
        return None

    def ms(values):
        return {
            "mean": round(statistics.fmean(values), 2),
            "std": round(statistics.pstdev(values), 2) if len(values) > 1 else 0.0,
        }

    # downsample the time series to ~150 points for an interactive plot
    n = len(crossing)
    step = max(1, n // 150)
    idx = list(range(0, n, step))
    have_time = times[0] == times[0]  # False if NaN
    series = {
        "time_ns": [round(times[k] / 1000.0, 2) for k in idx] if have_time else list(range(len(idx))),
        "crossing": [round(crossing[k], 2) for k in idx],
        "incident": [round(incident[k], 2) for k in idx],
    }
    # Crossing/incident angles are reported as-is; tilted and reverse-polarity
    # dockings with lower (or negative) angles are genuinely observed, so no value
    # is excluded or flagged against a reference band.
    return {"crossing_deg": ms(crossing), "incident_deg": ms(incident),
            "n_frames": n, "series": series}


def occupancy_classes(occ: dict | None) -> list[dict]:
    """Pull the per-interaction-class occupancy down to the fields the UI shows."""
    if not isinstance(occ, dict):
        return []
    out = []
    for c in occ.get("interaction_classes", []) or []:
        out.append(
            {
                "class": c.get("interaction_class"),
                "n_pairs": c.get("n_pairs"),
                "mean_occupancy": round(c.get("mean_occupancy", 0.0), 3),
                "stable_fraction": round(c.get("stable_pair_fraction", 0.0), 3),
            }
        )
    return out


def build_analysis(traj_id: str, source: Path) -> dict | None:
    """Read all per-stage summaries for one trajectory into a compact dict."""

    def sj(key: str):
        stage = {
            "identity": "identity",
            "rmsf": "rmsf",
            "bsa": "bsa",
            "contact_occupancy": "contact",
            "contact_annotation": "contact",
            "hbond": "hbond", "hbond_occupancy": "hbond",
            "saltbridge": "saltbridge", "saltbridge_occupancy": "saltbridge",
            "hydrophobic": "hydrophobic", "hydrophobic_occupancy": "hydrophobic",
            "pipi": "pipi", "pipi_occupancy": "pipi",
            "cationpi": "cationpi", "cationpi_occupancy": "cationpi",
        }[key]
        return load_json(source / stage / traj_id / JSON_PATHS[key])

    identity = sj("identity")
    rmsf = sj("rmsf")
    bsa = sj("bsa")
    contact_occ = sj("contact_occupancy")
    contact_ann = sj("contact_annotation")

    # If essentially nothing ran for this trajectory, skip it.
    if not any([identity, rmsf, bsa, contact_occ, contact_ann]):
        return None

    out: dict = {"traj_id": traj_id, "stages_present": []}

    if identity:
        out["identity"] = identity
        out["stages_present"].append("identity")
    if rmsf:
        out["rmsf"] = {
            "mean_rmsf_angstrom": rmsf.get("mean_rmsf_angstrom"),
            "max_rmsf_angstrom": rmsf.get("max_rmsf_angstrom"),
            "tcr_mean_rmsf_angstrom": rmsf.get("tcr_mean_rmsf_angstrom"),
            "phla_mean_rmsf_angstrom": rmsf.get("phla_mean_rmsf_angstrom"),
            "n_residues": rmsf.get("n_residues"),
            "n_frames": rmsf.get("n_frames"),
        }
        out["stages_present"].append("rmsf")
    if bsa:
        out["bsa"] = {
            "buried_surface_area": stat_block(bsa.get("buried_surface_area")),
            "interface_ratio": stat_block(bsa.get("interface_ratio")),
            "n_frames": bsa.get("n_frames"),
        }
        out["stages_present"].append("bsa")
    if contact_occ or contact_ann:
        out["contact"] = {
            "n_pairs": (contact_occ or {}).get("n_pairs"),
            "n_stable_pairs": (contact_occ or {}).get("n_stable_pairs"),
            "n_transient_pairs": (contact_occ or {}).get("n_transient_pairs"),
            "classes": occupancy_classes(contact_occ),
            "annotation": contact_ann,
        }
        out["stages_present"].append("contact")

    interactions = {}
    for name, (ann_key, occ_key) in INTERACTIONS.items():
        ann = sj(ann_key)
        occ = sj(occ_key)
        if not ann and not occ:
            continue
        interactions[name] = {
            "n_total_pairs": (ann or {}).get("n_total_pairs"),
            "n_cdr3_pairs": (ann or {}).get("n_cdr3_pairs"),
            "n_stable_pairs": (occ or {}).get("n_stable_pairs"),
            "classes": occupancy_classes(occ),
        }
        out["stages_present"].append(name)
    if interactions:
        out["interactions"] = interactions

    angle = angle_stats(source / "angle" / traj_id / "angles" / "docking_angles.csv")
    if angle:
        out["angle"] = angle
        out["stages_present"].append("angle")

    # NOTE: the interface-clustering stage is deliberately NOT ingested. Its
    # descriptors (n_clusters, largest_cluster population/dwell) are not reproducible
    # across replicas — the dominant-cluster population varies by a median of 29.6
    # percentage points between independent runs, 34.6% of the library range, because
    # the clustering cuts a dendrogram at a fixed absolute cutoff over a
    # per-trajectory max-normalised distance matrix. See pipeline.interface_descriptors.

    return out


def copy_figures(traj_id: str, source: Path, dest_dir: Path) -> list[str]:
    copied = []
    for (stage, rel), name in FIGURES.items():
        src = source / stage / traj_id / rel
        if src.exists():
            shutil.copy2(src, dest_dir / name)
            copied.append(name)
    return copied


def _merge_peptide_struct(analysis: dict, peptide_hla: dict) -> None:
    """Add HLA-contact / SASA / anchor flag to each peptide_table row.

    NOTE: ``struct_metrics`` now patches analysis.json directly with the identical
    rule (it is CSV-free and no longer needs the md_analysis source root); this
    legacy helper is kept for the full-``extract_analysis`` path and mirrors it.
    An anchor is buried (low absolute per-residue SASA < ``ANCHOR_SASA_MAX_NM2`` —
    the saturating ``hla_contact`` cannot resolve burial) AND ignored by the TCR
    (``tcr_contact < 0.20``): the class-I P2 and C-terminal (PΩ) pocket anchors.
    """
    table = (analysis.get("interface") or {}).get("peptide_table")
    struct_metrics.merge_peptide_table(table, peptide_hla)


def _sync_identity_tcr(analysis: dict, meta: dict) -> None:
    """Overwrite the served identity.tcr_identity with the ANARCI-correct meta.tcr.

    The upstream identity stage fills tcr_identity from a regex CDR detector that
    mislabels CDR3 (e.g. reports a CDR1 fragment) and leaves V/J genes empty. The
    UI reads meta.tcr, but analysis.json is a downloadable artifact, so we mirror
    the corrected values into it to keep the two consistent.
    """
    tcr = meta.get("tcr") or {}
    ident = (analysis.get("identity") or {}).get("tcr_identity")
    chains = tcr.get("chains") or {}
    if ident is None or not chains:
        return
    for pos in ("alpha", "beta"):
        c = chains.get(pos)
        if not c:
            continue
        if c.get("cdr3"):
            ident[f"cdr3_{pos}_sequence"] = c["cdr3"]
        if c.get("v_gene"):
            ident[f"{pos}_v_gene"] = c["v_gene"]
        if c.get("j_gene"):
            ident[f"{pos}_j_gene"] = c["j_gene"]
        ident[f"{pos}_genotype_confidence"] = "high"
    ident["source"] = "anarci_imgt"
    ident["genotype_source"] = "anarci"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Curate analysis output into web_data assets.")
    parser.add_argument("--source", type=Path, required=True, help="Analysis output root (from run_analysis).")
    parser.add_argument("--ids", default=None, help="Comma-separated trajectory ids (default: all in web_data).")
    args = parser.parse_args(argv)

    if args.ids:
        ids = [s.strip() for s in args.ids.split(",")]
    else:
        # List the trajectories that actually exist in WEB_DATA (the assets we are
        # curating), NOT SOURCE_ROOT — otherwise a wrong IMMUNO_SOURCE_ROOT silently
        # produces a mismatched id list and skips every live system.
        ids = sorted(d.name for d in config.WEB_DATA.iterdir()
                     if d.is_dir() and (d / config.OUT_META).exists())

    n_ok = 0
    for traj_id in ids:
        analysis = build_analysis(traj_id, args.source)
        if analysis is None:
            print(f"[skip] {traj_id}: no analysis output found")
            continue
        dest_dir = config.WEB_DATA / traj_id / "analysis"
        dest_dir.mkdir(parents=True, exist_ok=True)
        figures = copy_figures(traj_id, args.source, dest_dir)
        analysis["figures"] = figures
        analysis["provenance"] = ANALYSIS_PROVENANCE

        # Curate the residue-level contact matrix as a downloadable CSV (bravo §8.1.7).
        contact_src = (args.source / "contact" / traj_id
                       / "analysis/contacts/residue_contact_frequencies.csv")
        if contact_src.exists():
            shutil.copy2(contact_src, dest_dir / "contacts.csv")

        # CDR contact decomposition (correct ANARCI/IMGT boundaries) supersedes
        # the regex-based CDR labels in the contact stage. Requires annotate_tcr
        # to have written meta.tcr first; silently skipped otherwise.
        meta_path = config.WEB_DATA / traj_id / config.OUT_META
        if meta_path.exists():
            meta = json.loads(meta_path.read_text())
            _sync_identity_tcr(analysis, meta)
            tcr_cdr = cdr_contacts.decompose(traj_id, meta, args.source)
            if tcr_cdr is not None:
                analysis["tcr_cdr"] = tcr_cdr
            interface = interface_metrics.compute(traj_id, meta, args.source)
            if interface is not None:
                analysis["interface"] = interface
            # region-keyed RMSF flexibility profile (re-embed; reads the rmsf CSV)
            rmsf_prof = interface_metrics.rmsf_profile(traj_id, args.source)
            if rmsf_prof is not None:
                analysis["rmsf_profile"] = rmsf_prof
        # Per-region RMSD is trajectory-based and produced by the separate
        # split_rmsd step; re-embed its sidecar so re-running this curation keeps it.
        rmsd_sidecar = dest_dir / "rmsd_regions.json"
        if rmsd_sidecar.exists():
            rr = json.loads(rmsd_sidecar.read_text())
            analysis["rmsd_regions"] = rr
            if rr.get("fnat"):          # Fnat/Q rides in the split_rmsd sidecar; surface it top-level
                analysis["fnat"] = rr["fnat"]
        # Trajectory-derived structural metrics (anchor contacts / SASA / bulge / Rg).
        struct_sidecar = dest_dir / "struct_metrics.json"
        if struct_sidecar.exists():
            struct = json.loads(struct_sidecar.read_text())
            if struct.get("geometry"):
                analysis["geometry"] = struct["geometry"]
            if struct.get("bsa_decomposition"):
                analysis["bsa_decomposition"] = struct["bsa_decomposition"]
                # Cross-check against md_analysis's one-sided BSA (self-consistency gate).
                recon = struct_metrics.bsa_reconciliation(analysis.get("bsa"), struct["bsa_decomposition"])
                if recon:
                    analysis["bsa_decomposition"]["reconciliation"] = recon
            _merge_peptide_struct(analysis, struct.get("peptide_hla") or {})
        # Concerted-motion (DCCM + dynamic network) is trajectory-based; re-embed its
        # sidecar and re-register its PNGs so re-running this curation keeps them.
        cm_sidecar = dest_dir / "concerted_motion.json"
        if cm_sidecar.exists():
            cm = json.loads(cm_sidecar.read_text())
            analysis["concerted_motion"] = cm
            for png in (cm.get("figures") or []):
                if (dest_dir / png).exists() and png not in analysis["figures"]:
                    analysis["figures"].append(png)
        # Essential-dynamics (PC1 cosine + cross-replica RMSIP) and peptide dPCA are
        # trajectory-based; re-embed their sidecars so re-running curation keeps them.
        for _sc_name, _sc_key in (("essential_dynamics.json", "essential_dynamics"),
                                  ("peptide_dihedrals.json", "peptide_dpca"),
                                  ("tcr_cdr3_dpca.json", "tcr_cdr3_dpca"),
                                  ("coupled_states.json", "coupled_states")):
            _sc = dest_dir / _sc_name
            if _sc.exists():
                analysis[_sc_key] = json.loads(_sc.read_text())

        # NOTE: per-trajectory confidence scoring (dynamics_trust + a composite A/B/C
        # trust_tier) is intentionally NOT surfaced — a public trajectory database
        # should not grade its own entries. The underlying honesty signals remain
        # available as objective fields (essential_dynamics PC1 cosine content / RMSIP,
        # bsa_decomposition reconciliation, tcr_cdr reliability); we just never roll
        # them into a score.

        (dest_dir / "analysis.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")
        n_ok += 1
        print(f"[ok]   {traj_id}: {len(analysis['stages_present'])} stages, {len(figures)} figures")

    print(f"\nCurated analysis for {n_ok}/{len(ids)} trajectories into {config.WEB_DATA}/<id>/analysis/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
