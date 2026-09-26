#!/usr/bin/env python
"""Freeze one pooled dihedral-PCA landscape per side, shared by an entry's three replicas.

WHY THIS IS RECOMPUTED AND NOT READ.  The release stores a dihedral-PCA landscape
for every trajectory separately, and each one is expressed in its own basis: the
components come from a singular value decomposition of that trajectory's own
feature matrix, so PC1 of run 1 is a different linear combination of torsions and
distances from PC1 of run 2, the eigenvectors are not served, and the sign of each
is arbitrary.  For 1nam the three peptide surfaces carry 39%, 17% and 22% of the
variance on their first component and resolve 4, 2 and 3 basins.  Stacking those
three pictures on one pair of axes would be meaningless -- the axes do not mean
the same thing -- which is why this script exists.

WHAT IT DOES INSTEAD.  It rebuilds the same three feature blocks the pipeline
builds -- backbone phi/psi, side-chain chi1/chi2, pairwise Calpha-Calpha
distances, as {cos, sin} pairs for the angles -- on all three replicas at once,
block-normalises the pooled matrix, and takes ONE decomposition of it.  Every
frame of all three runs is then a point in one common plane, and the three can be
compared.  The features are internal coordinates and so invariant to how each run
happens to be superposed, which is what makes pooling legitimate; the runs share a
topology, so the columns correspond one-to-one.

The recipe is not reimplemented.  ``pipeline.peptide_dihedrals`` is imported and
its ``_kde_grid`` and ``_detect_basins`` are called directly, so the free-energy
convention (-kT ln rho at 310 K, shifted to its own minimum, capped at 6 kcal/mol),
the 80x80 grid, the 2.5 kcal/mol basin ceiling, the 8% population floor and the
four-basin cap are the served ones and not new choices.  Running the per-run path
through this file reproduces the served variance fractions to the fourth decimal,
which is checked below and printed.

WHAT IS WRITTEN.  For each side -- peptide, and CDR3 alpha+beta -- the pooled
surface on the common grid, its basins, and for each replica: its own surface on
that same grid, and how its frames distribute over the pooled basins.  That last
number is the point of the exercise.  Three replicas whose frames land in the same
basins in similar proportions have found the same wells; a replica whose frames
sit in a basin the other two never visit would be the finding.

This quantity is NOT part of the release.  It is a validation figure computed from
the deposited frames, and the script that computes it is deposited with the paper.

    IMMUNO_WEB_DATA=... python build_pooled_fel.py [pdb_id]

Reads the served trajectories read-only; writes ``pooled_fel.json`` here.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("IMMUNO_WEB_DATA", "/home/xmy/work/data/immunotrace/web_data_1000")
sys.path.insert(0, "/tmp/mdtraj_pkg")          # mdtraj, kept out of the conda envs
sys.path.insert(0, "/home/xmy/work/Immuno-Dyn")

import numpy as np
import mdtraj as md

from pipeline.peptide_dihedrals import (
    GRID, _detect_basins, _groove_ca, _kde_grid, _peptide_chain)
from pipeline.tcr_cdr3 import FW_MARGIN, _res_atoms, _res_ca, _tcr_struct_chain

R = Path(os.environ["IMMUNO_WEB_DATA"])
HERE = Path(__file__).resolve().parent
OUT = HERE / "pooled_fel.json"

EXEMPLAR = "1nam"               # the entry Figure 4 draws; see build_landscape.py
RUNS = (1, 2, 3)
PAD = 0.08                      # grid margin, as in the pipeline
ISO = 1.0                       # kcal/mol; the level each replica's outline is drawn at


# ---------------------------------------------------------------- feature blocks

def _norm(B):
    """Mean-centre and divide by the first singular value -- the pipeline's block
    balancing, verbatim, so no block dominates the joint decomposition."""
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


def _blocks(sub):
    """Raw, UNNORMALISED blocks for one atom selection: the normalisation is left
    to the caller because it must be done once over the pooled frames."""
    _, phi = md.compute_phi(sub)
    _, psi = md.compute_psi(sub)
    try:
        _, chi1 = md.compute_chi1(sub)
    except Exception:                                    # noqa: BLE001
        chi1 = None
    try:
        _, chi2 = md.compute_chi2(sub)
    except Exception:                                    # noqa: BLE001
        chi2 = None
    ca = sub.topology.select("name CA")
    dpair = None
    if ca.size >= 3:
        ii, jj = np.triu_indices(ca.size, k=2)
        if ii.size:
            pos = sub.xyz[:, ca, :]
            dpair = np.linalg.norm(pos[:, ii, :] - pos[:, jj, :], axis=2)
    return dict(phi=phi, psi=psi, chi1=chi1, chi2=chi2, dpair=dpair)


def _stack(blocks):
    """Concatenate a list of per-run block dicts frame-wise, keeping the keys."""
    out = {}
    for k in ("phi", "psi", "chi1", "chi2", "dpair"):
        vals = [b[k] for b in blocks]
        out[k] = None if any(v is None for v in vals) else np.concatenate(vals, axis=0)
    return out


def _design(b):
    parts = [x for x in (_ang(b["phi"], b["psi"]), _ang(b["chi1"], b["chi2"]),
                         _norm(b["dpair"])) if x is not None]
    return np.concatenate(parts, axis=1)


def _scores(b):
    """First two components of the block-balanced design matrix, sign-fixed.

    The sign of a singular vector is arbitrary.  It does not matter within one
    decomposition, but it would make this file's output depend on the LAPACK build,
    so it is pinned: each component is flipped so that its largest-magnitude
    loading is positive."""
    X = _design(b)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    var = S ** 2
    vf = (var / var.sum()) if var.sum() > 0 else np.zeros_like(var)
    sc = U[:, :2] * S[:2]
    for k in range(2):
        if Vt[k][np.argmax(np.abs(Vt[k]))] < 0:
            sc[:, k] *= -1.0
    return sc[:, 0], sc[:, 1], float(vf[0]), float(vf[1]), X.shape[1]


# ---------------------------------------------------------------- selections

def _peptide_sub(tid):
    wd = R / tid
    meta = json.loads((wd / "meta.json").read_text())
    t = md.load(str(wd / "traj.xtc"), top=str(wd / "topology.pdb"))
    pep = _peptide_chain(meta, t.topology)
    ci = next((c.index for c in t.topology.chains if c.chain_id == pep), -1)
    sel = t.topology.select(f"chainid {ci}")
    mhc = next((c for c, r in (meta.get("chain_roles") or {}).items() if r == "MHC"), None)
    g = _groove_ca(t.topology, mhc)
    if len(g) >= 4:
        t.superpose(t, 0, atom_indices=np.asarray(g))
    return t.atom_slice(sel)


def _cdr3_sub(tid):
    wd = R / tid
    meta = json.loads((wd / "meta.json").read_text())
    tc = (json.loads((wd / "analysis" / "analysis.json").read_text()).get("tcr_cdr") or {})
    cdr = tc.get("cdr_resids") or {}
    a_cdr, b_cdr = cdr.get("alpha") or {}, cdr.get("beta") or {}
    ach, bch = _tcr_struct_chain(meta, "alpha"), _tcr_struct_chain(meta, "beta")
    t = md.load(str(wd / "traj.xtc"), top=str(wd / "topology.pdb"))
    top0 = t.topology
    atoms = sorted(_res_atoms(top0, ach, a_cdr.get("cdr3") or [])
                   + _res_atoms(top0, bch, b_cdr.get("cdr3") or []))

    def fw(chain, cd):
        allcdr = set(cd.get("cdr1") or []) | set(cd.get("cdr2") or []) | set(cd.get("cdr3") or [])
        cut = max(cd.get("cdr3") or [0]) + FW_MARGIN
        return _res_ca(top0, chain, lambda s: s <= cut and s not in allcdr)

    fw_ca = fw(ach, a_cdr) + fw(bch, b_cdr)
    if len(fw_ca) >= 8:
        t.superpose(t, 0, atom_indices=np.asarray(fw_ca))
    return t.atom_slice(np.asarray(atoms))


# ---------------------------------------------------------------- one side

def side(name, subs, stored):
    """Pool the replicas of one side and report where each of them lands."""
    per = [_blocks(s) for s in subs]
    widths = {tuple((0 if v is None else v.shape[1]) for v in
                    (b["phi"], b["psi"], b["chi1"], b["chi2"], b["dpair"])) for b in per}
    assert len(widths) == 1, f"{name}: replicas disagree on feature width {widths}"
    nfr = [s.n_frames for s in subs]

    # solo pass, for the check that this file reproduces what is served
    solo = []
    for b, st in zip(per, stored):
        _, _, v1, v2, _ = _scores(b)
        solo.append(dict(recomputed=[round(v1, 4), round(v2, 4)],
                         stored=[round(st["pc1_var_frac"], 4), round(st["pc2_var_frac"], 4)],
                         basins=len(st["basins"])))

    pc1, pc2, v1, v2, ncol = _scores(_stack(per))
    p1 = (pc1.max() - pc1.min()) * PAD or 0.1
    p2 = (pc2.max() - pc2.min()) * PAD or 0.1
    gx = np.linspace(pc1.min() - p1, pc1.max() + p1, GRID)
    gy = np.linspace(pc2.min() - p2, pc2.max() + p2, GRID)
    z = _kde_grid(pc1, pc2, gx, gy)
    basins = _detect_basins(pc1, pc2, gx, gy, z)
    cen = np.array([[b["pc1"], b["pc2"]] for b in basins])

    runs, off = [], 0
    for i, n in enumerate(nfr):
        a, b_ = off, off + n
        off = b_
        rx, ry = pc1[a:b_], pc2[a:b_]
        # where this replica's OWN served basins land in the pooled plane.  The
        # release stores each basin's representative frame index, and that frame
        # is one of the frames just projected, so the served decomposition can be
        # placed on the shared map without being recomputed.
        own = [dict(k=k + 1, pop=bb.get("pop", 0.0), dG=bb.get("dG", 0.0),
                    pc1=round(float(rx[bb["frame"]]), 4),
                    pc2=round(float(ry[bb["frame"]]), 4))
               for k, bb in enumerate(stored[i]["basins"])
               if 0 <= bb.get("frame", -1) < n]
        d2 = ((np.column_stack([rx, ry])[:, None, :] - cen[None, :, :]) ** 2).sum(axis=2)
        assign = d2.argmin(axis=1)
        runs.append(dict(
            run=RUNS[i], n_frames=int(n),
            # the replica's own surface on the POOLED grid, so its outline can be
            # drawn over the pooled one and the two are the same map
            z=_kde_grid(rx, ry, gx, gy),
            occ=[round(float((assign == k).mean()), 4) for k in range(len(basins))],
            own_basins=own,
            centroid=[round(float(rx.mean()), 4), round(float(ry.mean()), 4)],
        ))

    # how far apart the three replicas' own occupied regions are, in the only
    # currency that is comparable here: the fraction of the pooled ISO-level
    # region each pair shares
    def region(zz):
        return np.asarray(zz, dtype=float) <= ISO
    regs = [region(r["z"]) for r in runs]
    jac = []
    for i in range(len(regs)):
        for j in range(i + 1, len(regs)):
            u = (regs[i] | regs[j]).sum()
            jac.append(round(float((regs[i] & regs[j]).sum() / u), 3) if u else 0.0)

    return dict(
        name=name, n_features=int(ncol), n_frames=[int(v) for v in nfr],
        pc1_var_frac=round(v1, 4), pc2_var_frac=round(v2, 4),
        grid=dict(x=[round(float(v), 5) for v in gx], y=[round(float(v), 5) for v in gy]),
        z=z, basins=basins, runs=runs, iso=ISO, jaccard=jac, solo=solo,
    )


def main():
    pid = sys.argv[1] if len(sys.argv) > 1 else EXEMPLAR
    tids = [f"{pid}_run{r}" for r in RUNS]
    an = [json.loads((R / t / "analysis" / "analysis.json").read_text()) for t in tids]
    out = dict(pdb_id=pid, runs=list(RUNS), grid_n=GRID,
               peptide=side("peptide", [_peptide_sub(t) for t in tids],
                            [a["peptide_dpca"] for a in an]),
               cdr3=side("cdr3", [_cdr3_sub(t) for t in tids],
                         [a["tcr_cdr3_dpca"] for a in an]))
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    for k in ("peptide", "cdr3"):
        d = out[k]
        print(f"{k:9s} {d['n_features']} features, {sum(d['n_frames'])} pooled frames, "
              f"PC1 {d['pc1_var_frac']:.0%} PC2 {d['pc2_var_frac']:.0%}, "
              f"{len(d['basins'])} pooled basins")
        for s, r in zip(d["solo"], d["runs"]):
            print(f"    run {r['run']}  solo var {s['recomputed']} vs served {s['stored']}"
                  f"  ({s['basins']} solo basins)   pooled occupancy "
                  + " ".join(f"{v:.0%}" for v in r["occ"]))
        print(f"    pairwise overlap of the <={ISO:.0f} kcal/mol regions: "
              + ", ".join(f"{v:.0%}" for v in d["jaccard"]))
    print(f"{OUT.name}   {OUT.stat().st_size / 1024:.0f} kB")


if __name__ == "__main__":
    main()
