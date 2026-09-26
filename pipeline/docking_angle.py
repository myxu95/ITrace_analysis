"""Correct TCR-pMHC crossing / incident angle, recomputed from the structure.

The md_analysis `angle` stage mis-resolves the groove / TCR axes and reports clear
diagonal dockings as near-parallel (e.g. the four A6/HLA-A2 complexes came out
2.7-8.5 deg instead of the canonical 22-69 deg range, and inconsistently). This recomputes both angles
per frame with chains identified by ROLE (not by chain id, which is not consistent
across the PDB) and overwrites the `angle` block of each analysis.json. See
docs/docking_angle.md for the algorithm and its validation.

The crossing angle is stored as a DIRECTED angle in [0, 180] deg: forward (canonical)
dockings keep their acute 0-90 value, while reverse-polarity dockings are reflected
past 90 deg (180 - acute) so they read as the >90 deg / "toward 180" crossings they
are, rather than being folded into the canonical band. The unsigned acute magnitude
is kept alongside as `crossing_magnitude_deg`.

    IMMUNO_WEB_DATA=.../immuno-dyn python -m pipeline.docking_angle [--ids a,b,c]
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from . import config

# Polarity is only meaningful when the inter-domain vector has a real in-plane
# component along the groove; near-perpendicular frames (crossing > ~80 deg) carry
# a noisy sign, so they are excluded from the polarity vote rather than counted.
POLARITY_MIN_COS = 0.174          # cos(80 deg)
MHC_PLATFORM = 180                 # first N residues of the MHC heavy chain (alpha1+alpha2)
TCR_VARIABLE = 110                 # first N residues of each TCR chain (variable domain)
SERIES_POINTS = 150                # downsampled time series for the detail plot

# Antigen chain for structures whose upstream peptide detector left peptide_chain
# null (non-standard antigens, e.g. lipopeptides), so the angle is computed here
# instead of skipped and left with the stale md_analysis value.
_MANUAL_PEPTIDE_CHAIN = {"7byd": "C"}   # GLY-GLY-ALA-ILE lipopeptide on Mamu-B*05104


def _chain_ca(top, chid):
    """Ordered Cα atom indices for a chain (by chain id)."""
    for ch in top.chains:
        if ch.chain_id == chid:
            return [next((a.index for a in r.atoms if a.name == "CA"), None) for r in ch.residues]
    return []


def _mhc_chain(meta, pep_ch, a_ch, b_ch):
    """MHC heavy chain = the largest chain that is not peptide / TCRα / TCRβ."""
    cand = [c for c in (meta.get("chains") or []) if c["id"] not in (pep_ch, a_ch, b_ch)]
    return max(cand, key=lambda c: c["n_residues"])["id"] if cand else None


def _frame_angles(X, mhc_ca, va, vb, pep):
    """(crossing, incident, polarity) for one frame. crossing/incident in degrees;
    polarity = +1 if the Vα→Vβ vector points N→C along the groove (canonical
    docking), −1 if reversed (reverse-polarity docking), 0 if undefined."""
    G = X[mhc_ca]
    c = G.mean(0)
    _, _, vt = np.linalg.svd(G - c)
    normal = vt[2]                                   # out-of-platform direction
    pv = X[pep[-1]] - X[pep[0]]                       # peptide N->C (groove long direction)
    long_axis = vt[0] if abs(vt[0] @ pv) >= abs(vt[1] @ pv) else vt[1]
    if long_axis @ pv < 0:
        long_axis = -long_axis
    va_c, vb_c = X[va].mean(0), X[vb].mean(0)
    tcr = vb_c - va_c                                # Vα domain -> Vβ domain (inter-domain vector)
    tnorm = np.linalg.norm(tcr)
    if tnorm < 1e-6:
        return 0.0, 0.0, 0.0
    # Crossing ("TCR twist"): in-plane angle of the inter-domain vector to the groove long axis.
    tcr_p = tcr - (tcr @ normal) * normal            # project onto groove plane
    pnorm = np.linalg.norm(tcr_p)
    if pnorm < 1e-6:
        crossing, polarity = 0.0, 0.0
    else:
        proj = tcr_p @ long_axis                     # signed component along groove N->C
        cosa = abs(proj) / (pnorm * np.linalg.norm(long_axis))
        crossing = float(np.degrees(np.arccos(np.clip(cosa, -1.0, 1.0))))
        # Polarity: which way the inter-domain vector points along the groove. The
        # acute crossing magnitude folds canonical (proj>0) and reversed (proj<0)
        # together; this sign keeps them distinguishable, and downstream it reflects
        # reversed frames past 90 deg into the directed 0-180 crossing convention.
        # (Dataset: ~96% forward; the reversed set is the known reverse-polarity
        # B17/H2-Db, 5sws/5swz/7jwi and 9gv7 dockings.) Undefined (0) for
        # near-perpendicular frames where the sign is noise.
        polarity = 0.0 if cosa < POLARITY_MIN_COS else (1.0 if proj >= 0 else -1.0)
    # Incident (tilt): how far the TCR tips off vertical over the groove — the angle
    # between the platform normal and the pMHC->TCR centre-of-mass "approach" vector.
    # (The previous version measured the inter-domain vector's own out-of-plane
    # elevation — a small "roll", not the conventional tilt of the TCR over the pMHC.)
    approach = 0.5 * (va_c + vb_c) - 0.5 * (c + X[pep].mean(0))
    anorm = np.linalg.norm(approach)
    incident = 0.0 if anorm < 1e-6 else float(
        np.degrees(np.arccos(np.clip(abs(approach @ normal) / anorm, 0.0, 1.0))))
    return crossing, incident, polarity


def compute(traj_id: str, meta: dict) -> dict | None:
    import mdtraj as md
    pep_ch = meta.get("peptide_chain") or _MANUAL_PEPTIDE_CHAIN.get(meta.get("pdb_id"))
    tc = (meta.get("tcr") or {}).get("chains") or {}
    a_ch = (tc.get("alpha") or {}).get("structural_chain")
    b_ch = (tc.get("beta") or {}).get("structural_chain")
    if not (pep_ch and a_ch and b_ch):
        return None
    mhc_ch = _mhc_chain(meta, pep_ch, a_ch, b_ch)
    if not mhc_ch:
        return None

    d = config.WEB_DATA / traj_id
    topo = d / config.OUT_TOPOLOGY
    traj = d / config.OUT_TRAJ
    if not topo.exists():
        return None
    t = md.load(str(traj), top=str(topo)) if traj.exists() else md.load(str(topo))
    top = t.topology

    def cas(chid):
        return [i for i in _chain_ca(top, chid) if i is not None]
    mhc_ca = cas(mhc_ch)[:MHC_PLATFORM]
    va, vb = cas(a_ch)[:TCR_VARIABLE], cas(b_ch)[:TCR_VARIABLE]
    pep = cas(pep_ch)
    if len(mhc_ca) < 20 or len(va) < 20 or len(vb) < 20 or len(pep) < 2:
        return None

    cross = np.empty(t.n_frames)
    inc = np.empty(t.n_frames)
    pol = np.empty(t.n_frames)
    for f in range(t.n_frames):
        cross[f], inc[f], pol[f] = _frame_angles(t.xyz[f], mhc_ca, va, vb, pep)
    # Polarity vote over frames where the sign is defined (exclude near-perpendicular).
    defined = pol != 0
    reversed_fraction = float((pol < 0).sum() / defined.sum()) if defined.any() else 0.0
    reversed_polarity = reversed_fraction > 0.5

    # Directed crossing angle in [0, 180]: forward frames keep their acute angle;
    # reverse-polarity frames are reflected past 90 deg (180 - acute), so a
    # reverse-polarity docking reads as a >90 deg ("toward 180") crossing instead of
    # being folded into the canonical 0-90 band. The acute magnitude is retained
    # separately as crossing_magnitude_deg. (pol: +1 forward, -1 reversed, 0 = sign
    # undefined near-perpendicular -> left at the acute value, which is ~90 there.)
    directed = np.where(pol < 0, 180.0 - cross, cross)

    def ms(a):
        return {"mean": round(float(a.mean()), 2), "std": round(float(a.std()), 2),
                "min": round(float(a.min()), 2), "max": round(float(a.max()), 2)}

    dur = meta.get("duration_ns")
    times = (np.linspace(0, dur, t.n_frames) if dur else np.arange(t.n_frames, dtype=float))
    idx = np.linspace(0, t.n_frames - 1, min(SERIES_POINTS, t.n_frames)).astype(int)
    return {
        "crossing_deg": ms(directed),            # directed 0-180 (>90 = reverse polarity)
        "crossing_magnitude_deg": ms(cross),     # unsigned acute 0-90 (provenance)
        "incident_deg": ms(inc),
        "n_frames": int(t.n_frames),
        "series": {
            "time_ns": [round(float(times[i]), 2) for i in idx],
            "crossing": [round(float(directed[i]), 2) for i in idx],
            "crossing_magnitude": [round(float(cross[i]), 2) for i in idx],
            "incident": [round(float(inc[i]), 2) for i in idx],
        },
        "reversed_polarity": bool(reversed_polarity),
        "reversed_polarity_fraction": round(reversed_fraction, 3),
        "method": "crossing:groove_pca_peptide_axis_tcr_interdomain_directed0to180;incident:pmhc_tcr_approach_tilt;polarity:interdomain_sign_along_groove",
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all).")
    args = ap.parse_args(argv)
    only = set(args.ids.split(",")) if args.ids else None

    web = config.WEB_DATA
    n_ok = n_skip = 0
    for meta_path in sorted(web.glob("*/" + config.OUT_META)):
        tid = meta_path.parent.name
        if only and tid not in only:
            continue
        an_path = meta_path.parent / "analysis" / "analysis.json"
        if not an_path.exists():
            n_skip += 1
            continue
        meta = json.loads(meta_path.read_text())
        try:
            angle = compute(tid, meta)
        except Exception as exc:  # never let one bad trajectory abort the batch
            print(f"[warn] {tid}: {exc}")
            angle = None
        if angle is None:
            n_skip += 1
            continue
        analysis = json.loads(an_path.read_text())
        analysis["angle"] = angle
        an_path.write_text(json.dumps(analysis, indent=2))
        n_ok += 1
    print(f"Docking angle recomputed for {n_ok} trajectories, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
