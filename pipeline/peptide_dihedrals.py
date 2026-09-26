"""Peptide dihedral-PCA (dPCA) + 2D free-energy landscape.

Alignment-free peptide conformational heterogeneity from THREE MFA-balanced feature
blocks — backbone (φ/ψ, {cos,sin}), side-chain (χ1, {cos,sin}), and global shape
(pairwise Cα–Cα distances) — each normalised by its first singular value so no block
dominates; PCA to 2D, and build a −kT ln P free-energy landscape over (PC1, PC2). This
is fuller than a pure-backbone dPCA: it sees backbone plasticity, side-chain rotamer
substates, AND overall bulge/extension shape together. Substates (bulge up/down,
terminus tilt, rotamer flips) and their populations come from k-means on the 2D scores.
Directly addresses peptide cross-reactivity / escape plasticity for an immunology audience.

Purely trajectory-derived (topology.pdb + traj.xtc); globs WEB_DATA, no --source.
Writes a sidecar ``analysis/peptide_dihedrals.json`` (re-embedded by
extract_analysis under ``peptide_dpca``) and patches analysis.json directly.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.peptide_dihedrals [--ids a,b]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md
from scipy import ndimage
from scipy.stats import gaussian_kde

from . import config
from .struct_metrics import HELIX1, HELIX2, HELIX1_WIDE, HELIX2_WIDE, _ca_in_spans

SIDECAR = "peptide_dihedrals.json"
KT_KCAL = 0.0019872041 * 310.0   # kcal/mol at 310 K
GRID = 80                         # KDE evaluation grid (per axis) for the free-energy map
                                  # (the frontend bicubic-upsamples this into smooth iso-lines;
                                  #  a finer native grid only shrinks the facets, not the steps)
FEL_CAP = 6.0                     # cap ΔG (kcal/mol) FAR above the plotted range (~2.4): the
                                  # far-field stays smooth in the visible band (a low, hard cap
                                  # made a blocky plateau cliff → jagged outer contours). The
                                  # frontend colour-clips the far-field to a uniform light fill.
BASIN_ZMAX = 2.5                  # only accept ΔG minima in the well-sampled band as basins
                                  # (was tied to the old low cap; now an explicit threshold)
KMAX = 3                          # k-means substates


def _rotamer(a: float | None) -> str | None:
    """3-well rotamer label from a side-chain χ angle (degrees): g+ (~+60), t (~180),
    g- (~-60). None passes through (residue lacks that χ)."""
    if a is None:
        return None
    a = ((a + 180.0) % 360.0) - 180.0
    if 0.0 <= a < 120.0:
        return "g+"
    if -120.0 <= a < 0.0:
        return "g-"
    return "t"


def _groove_ca(top, mhc_chain):
    """Cα atom indices of the MHC α1/α2 groove helices (canonical spans, widen if sparse).
    The superposition target so basin representatives overlay in a common groove frame."""
    if not mhc_chain:
        return []
    ca = _ca_in_spans(top, mhc_chain, (HELIX1, HELIX2))
    if len(ca) < 20:
        ca = _ca_in_spans(top, mhc_chain, (HELIX1_WIDE, HELIX2_WIDE))
    return ca


def _groove_atoms(top, mhc_chain):
    """ALL atoms of the α1/α2 groove-helix residues (for the 3D-context reference PDB)."""
    idx = []
    for ch in top.chains:
        if ch.chain_id != mhc_chain:
            continue
        for r in ch.residues:
            if any(lo <= r.resSeq <= hi for lo, hi in (HELIX1, HELIX2)):
                idx += [a.index for a in r.atoms]
    return idx


def _kde_grid(x, y, gx, gy, periodic=False, period=360.0):
    """ΔG = -kT ln ρ (kcal/mol, min-shifted + capped) from a Gaussian KDE of (x, y),
    evaluated on the gx×gy grid. periodic=True wraps by ±period (φ/ψ torus). Returns a
    Plotly z[y][x] list."""
    XX, YY = np.meshgrid(gx, gy)
    kde = gaussian_kde(np.vstack([x, y]))
    if periodic:
        dens = np.zeros(XX.size)
        for dx in (-period, 0.0, period):
            for dy in (-period, 0.0, period):
                dens += kde(np.vstack([XX.ravel() + dx, YY.ravel() + dy]))
        dens = dens.reshape(XX.shape)
    else:
        dens = kde(np.vstack([XX.ravel(), YY.ravel()])).reshape(XX.shape)
    Fe = -KT_KCAL * np.log(dens + 1e-12)
    Fe = np.minimum(Fe - Fe.min(), FEL_CAP)
    return [[round(float(Fe[iy, ix]), 2) for ix in range(len(gx))] for iy in range(len(gy))]


def _detect_basins(pc1, pc2, gx, gy, zlist, max_basins=4, min_pop=0.08):
    """Conformational basins = local minima of the ΔG grid. For each basin return its
    PC-space centre, ΔG floor, population (fraction of frames nearest that centre), and
    ``frame_frac`` — the representative frame's relative time position in [0,1], so the
    detail viewer can seek to the matching playback frame regardless of how many frames
    the playback file (traj_view.xtc / traj.xtc) actually carries."""
    Z = np.asarray(zlist, dtype=float)                 # (ny, nx): Z[iy, ix]
    ny, nx = Z.shape
    win = max(5, (nx // 8) | 1)                        # odd window ~1/8 axis: skip KDE ripple
    is_min = ndimage.minimum_filter(Z, size=win, mode="nearest") == Z
    is_min &= Z < BASIN_ZMAX                            # only well-sampled minima are basins
    if not is_min.any():
        return []
    lab, n = ndimage.label(is_min)
    frames = np.column_stack([pc1, pc2])
    centres = []
    for c in range(1, n + 1):
        ys, xs = np.where(lab == c)
        k = int(np.argmin(Z[ys, xs]))                  # deepest cell of the component
        iy, ix = int(ys[k]), int(xs[k])
        centres.append((float(Z[iy, ix]), float(gx[ix]), float(gy[iy])))
    if not centres:
        return []
    cen = np.array([[cx, cy] for _, cx, cy in centres])
    d2 = ((frames[:, None, :] - cen[None, :, :]) ** 2).sum(axis=2)   # (F, nb)
    assign = d2.argmin(axis=1)
    nF = len(pc1)
    out = []
    for i, (dG, cx, cy) in enumerate(centres):
        members = np.where(assign == i)[0]
        pop = members.size / nF
        rep = int(members[np.argmin(d2[members, i])]) if members.size else int(d2[:, i].argmin())
        out.append({"pc1": round(cx, 3), "pc2": round(cy, 3),
                    "dG": round(dG, 2), "pop": round(float(pop), 3),
                    "frame": rep, "frame_frac": round(rep / max(1, nF - 1), 4)})
    keep = [b for b in out if b["pop"] >= min_pop] or [min(out, key=lambda b: b["dG"])]
    keep.sort(key=lambda b: (b["dG"], -b["pop"]))       # deepest / most-populated first
    return keep[:max_basins]


def _peptide_chain(meta: dict, top) -> str | None:
    """meta.peptide_chain, else the shortest chain (peptides are 8–13 aa)."""
    pc = meta.get("peptide_chain")
    if pc:
        return pc
    best, best_n = None, 1e9
    for ch in top.chains:
        n = sum(1 for _ in ch.residues)
        if 4 <= n <= 25 and n < best_n:
            best, best_n = ch.chain_id, n
    return best


def _kmeans(X: np.ndarray, k: int, iters: int = 40, seed_stride: int = 0):
    """Tiny deterministic Lloyd k-means (no sklearn dep). Returns labels."""
    n = X.shape[0]
    if n < k:
        return np.zeros(n, dtype=int)
    # deterministic spread-out seeds: evenly spaced by PC1 rank
    order = np.argsort(X[:, 0])
    cent = X[order[np.linspace(0, n - 1, k).astype(int)]].copy()
    labels = np.zeros(n, dtype=int)
    for _ in range(iters):
        d = ((X[:, None, :] - cent[None, :, :]) ** 2).sum(axis=2)
        new = d.argmin(axis=1)
        if np.array_equal(new, labels):
            break
        labels = new
        for c in range(k):
            m = labels == c
            if m.any():
                cent[c] = X[m].mean(axis=0)
    return labels


def compute(traj_id: str, meta: dict) -> dict | None:
    wd = config.WEB_DATA / traj_id
    tp, xp = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
    if not (tp.exists() and xp.exists()):
        return None
    t = md.load(str(xp), top=str(tp))
    if t.n_frames < 10:
        return None
    pep = _peptide_chain(meta, t.topology)
    if not pep:
        return None
    pep_atoms = t.topology.select(f"chainid {next((c.index for c in t.topology.chains if c.chain_id == pep), -1)}")
    if pep_atoms.size == 0:
        return None
    # Groove-align the whole complex on the MHC α1/α2 helices (frame 0 = reference) so the
    # per-basin representative peptides overlay in one common frame and their differences read
    # clearly. dPCA features below (dihedrals + Cα distances) are rotation-invariant, so this
    # changes only the saved coordinates, not the landscape.
    mhc_chain = next((cid for cid, role in (meta.get("chain_roles") or {}).items() if role == "MHC"), None)
    groove_ca = _groove_ca(t.topology, mhc_chain)
    if len(groove_ca) >= 4:
        t.superpose(t, 0, atom_indices=np.asarray(groove_ca))
    else:
        mhc_chain = None
    sub = t.atom_slice(pep_atoms)
    phi_idx, phi = md.compute_phi(sub)   # (n_phi,4) atom quartets ; (F, n_phi) radians
    psi_idx, psi = md.compute_psi(sub)
    try:
        chi1_idx, chi1 = md.compute_chi1(sub)   # (F, n_chi1) side-chain χ1 (residues with one)
    except Exception:                    # noqa: BLE001
        chi1_idx, chi1 = None, None
    try:
        chi2_idx, chi2 = md.compute_chi2(sub)   # χ2 → side chains resolved past Cβ (TCR reads these)
    except Exception:                    # noqa: BLE001
        chi2_idx, chi2 = None, None
    # pairwise Cα–Cα distances (|i−j| ≥ 2; skip the near-constant i,i+1 bond) — a
    # superposition-free readout of the peptide's global backbone shape (bulge / extension
    # / kink), in nm; complements the local φ/ψ angles.
    dpair, ca = None, sub.topology.select("name CA")
    if ca.size >= 3:
        ii, jj = np.triu_indices(ca.size, k=2)
        if ii.size:
            pos = sub.xyz[:, ca, :]      # (F, nCA, 3)
            dpair = np.linalg.norm(pos[:, ii, :] - pos[:, jj, :], axis=2)   # (F, n_pair)

    # THREE MFA-balanced blocks — backbone (φ/ψ, {cos,sin}), side-chain (χ1+χ2, {cos,sin}),
    # global shape (Cα distances). Each is mean-centred and divided by its first singular
    # value so no block dominates the joint PCA. Gly/Ala carry no χ1 and drop out of that
    # block; a peptide too short for either simply uses the blocks it has. χ2 is included so
    # the basins resolve side-chain rotamers past Cβ (what the TCR actually reads), not just χ1.
    def _norm(B):
        if B is None or not B.shape[1]:
            return None
        B = B - B.mean(axis=0)
        s1 = float(np.linalg.svd(B, compute_uv=False)[0])
        return B / s1 if s1 > 1e-9 else B
    def _ang(*angs):
        cols = []
        for a in angs:
            if a is not None and a.shape[1]:
                cols += [np.cos(a), np.sin(a)]
        return _norm(np.concatenate(cols, axis=1)) if cols else None
    parts = [b for b in (_ang(phi, psi), _ang(chi1, chi2), _norm(dpair)) if b is not None]
    if not parts:
        return None
    Xc = np.concatenate(parts, axis=1)           # blocks already mean-centred + normalised
    if Xc.shape[1] < 2:
        return None
    n_chi1 = int(chi1.shape[1]) if chi1 is not None else 0
    n_chi2 = int(chi2.shape[1]) if chi2 is not None else 0
    n_dist = int(dpair.shape[1]) if dpair is not None else 0
    # PCA via SVD (robust, no covariance blow-up)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2
    var_frac = (var / var.sum()).tolist() if var.sum() > 0 else [0, 0]
    scores = U[:, :2] * S[:2]                     # (F, 2) PC scores
    if scores.shape[1] < 2:
        return None
    pc1, pc2 = scores[:, 0], scores[:, 1]

    # 2D free-energy landscape via Gaussian KDE (smooth + hole-free, unlike a sparse
    # histogram): F = -kT ln rho over a fine grid, shifted so min = 0 and capped so deep
    # basins are prominent and the under-sampled far-field is a flat high-energy plateau.
    pad1 = (pc1.max() - pc1.min()) * 0.08 or 0.1
    pad2 = (pc2.max() - pc2.min()) * 0.08 or 0.1
    gx = np.linspace(pc1.min() - pad1, pc1.max() + pad1, GRID)
    gy = np.linspace(pc2.min() - pad2, pc2.max() + pad2, GRID)
    try:
        z = _kde_grid(pc1, pc2, gx, gy)
    except (np.linalg.LinAlgError, ValueError):
        return None
    # 5 decimals, NOT 2: PC ranges can be ~0.2 wide, so 2-dp rounding collapses the 80 grid
    # coords to ~20 unique values (a degenerate, duplicated axis) → Plotly draws stair-stepped
    # contours. 5 dp keeps every grid line distinct so the iso-lines render smooth.
    xc = [round(float(v), 5) for v in gx]
    yc = [round(float(v), 5) for v in gy]

    # Basin centres (ΔG local minima) + the representative frame for each, so the detail
    # page can mark the wells and pop up each basin's conformation.
    try:
        basins = _detect_basins(pc1, pc2, gx, gy, z)
    except Exception:   # noqa: BLE001
        basins = []

    # Persist each basin's representative peptide conformation (groove-aligned, all-atom incl.
    # side chains) as a small PDB, plus a one-off MHC α1/α2 groove reference, so the detail
    # page can overlay all basins in 3D on a common groove frame and read the peptide
    # (backbone + side-chain rotamer) differences.
    an_dir = wd / "analysis"
    an_dir.mkdir(parents=True, exist_ok=True)
    for old in list(an_dir.glob("basin_*.pdb")) + [an_dir / "groove_ref.pdb"]:
        try:
            old.unlink()
        except OSError:
            pass
    for k, b in enumerate(basins):
        try:
            sub[b["frame"]].save_pdb(str(an_dir / f"basin_{k}.pdb"))
            b["pdb"] = f"basin_{k}.pdb"
        except Exception:   # noqa: BLE001
            pass
    groove_ref = None
    if mhc_chain:
        try:
            gidx = _groove_atoms(t.topology, mhc_chain)
            if gidx:
                t[0].atom_slice(np.asarray(gidx)).save_pdb(str(an_dir / "groove_ref.pdb"))
                groove_ref = "groove_ref.pdb"
        except Exception:   # noqa: BLE001
            groove_ref = None

    # Per-basin side-chain fingerprint: each peptide residue's χ1/χ2 (+ 3-well rotamer label)
    # at that basin's representative frame — so the 3D card can call out which residue's
    # rotamer flips between conformational states (the recognition-relevant, side-chain part).
    top = sub.topology
    chi1_res = ({top.atom(int(chi1_idx[j][0])).residue.index: j for j in range(chi1_idx.shape[0])}
                if chi1 is not None and chi1.shape[1] else {})
    chi2_res = ({top.atom(int(chi2_idx[j][0])).residue.index: j for j in range(chi2_idx.shape[0])}
                if chi2 is not None and chi2.shape[1] else {})
    for b in basins:
        fr = int(b["frame"])
        sc = []
        for rr in sorted(chi1_res):
            c1 = float(np.degrees(chi1[fr, chi1_res[rr]]))
            c2 = float(np.degrees(chi2[fr, chi2_res[rr]])) if rr in chi2_res else None
            sc.append({"position": rr + 1, "resname": top.residue(rr).name,
                       "chi1": round(c1, 1), "chi2": round(c2, 1) if c2 is not None else None,
                       "rot1": _rotamer(c1), "rot2": _rotamer(c2)})
        b["sidechain"] = sc

    # substates via k-means on the 2D scores
    labels = _kmeans(scores, KMAX)
    pops = sorted((float((labels == c).mean()) for c in range(KMAX)), reverse=True)
    pops = [round(p, 3) for p in pops if p > 0]
    n_states = sum(1 for p in pops if p >= 0.10)

    return {
        "peptide_chain": pep,
        "pc1_var_frac": round(float(var_frac[0]), 3),
        "pc2_var_frac": round(float(var_frac[1]), 3),
        "features": {"n_phi": int(phi.shape[1]), "n_psi": int(psi.shape[1]),
                     "n_chi1": n_chi1, "n_chi2": n_chi2, "n_ca_dist": n_dist},
        "n_substates": int(n_states),
        "substate_populations": pops,
        "fel": {"z": z, "x": xc, "y": yc,
                "x_label": "PC1", "y_label": "PC2", "unit": "kcal/mol"},
        "basins": basins,
        "groove_ref": groove_ref,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all in WEB_DATA).")
    args = ap.parse_args(argv)
    only = set(args.ids.split(",")) if args.ids else None

    n_ok = n_skip = 0
    for meta_path in sorted(config.WEB_DATA.glob("*/" + config.OUT_META)):
        tid = meta_path.parent.name
        if only and tid not in only:
            continue
        meta = json.loads(meta_path.read_text())
        try:
            res = compute(tid, meta)
        except Exception as exc:   # noqa: BLE001
            print(f"[warn] {tid}: {exc}")
            res = None
        if res is None:
            n_skip += 1
            continue
        an_dir = meta_path.parent / "analysis"
        an_dir.mkdir(parents=True, exist_ok=True)
        (an_dir / SIDECAR).write_text(json.dumps(res))
        an_json = an_dir / "analysis.json"
        if an_json.exists():
            data = json.loads(an_json.read_text())
            data["peptide_dpca"] = res
            an_json.write_text(json.dumps(data))
        n_ok += 1
    print(f"Peptide dPCA: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
