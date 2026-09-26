"""Compute a canonical viewing orientation for each pHLA-TCR complex.

Writes ``<web_data>/<id>/orient.json`` with two unit vectors used by the
thumbnail renderer to place a fixed camera in the standard immunology view
(pHLA at the bottom, peptide groove horizontal, TCR on top, facing the viewer):

    up   = pHLA-centre -> TCR-centre   (so the TCR sits at the top)
    view = camera direction, perpendicular to both `up` and the peptide axis

Chain roles are assigned geometrically, which is robust to the chain ordering
and to MHC/TCR size overlap that a residue-count heuristic gets wrong:
  - peptide  : the shortest chain (already recorded in meta.json)
  - beta-2 m : the ~100-residue chain (very conserved length)
  - MHC heavy: of the rest, the chain whose centroid is closest to beta-2m
               (the heavy chain and beta-2m together form the MHC molecule)
  - TCR a/b  : the remaining two chains

Usage:
    python -m pipeline.compute_orient                 # all trajectories
    python -m pipeline.compute_orient 1ao7_run2 ...   # only the given ids
"""
from __future__ import annotations

import json
import sys

import mdtraj as md
import numpy as np
from tqdm import tqdm

from . import config


def _ca_coords_by_chain(traj: md.Trajectory):
    """Return {chain_index: (n,3) CA coords} for the single frame."""
    top = traj.topology
    out = {}
    for ci, chain in enumerate(top.chains):
        idx = [a.index for a in chain.atoms if a.name == "CA"]
        if idx:
            out[ci] = traj.xyz[0, idx]
    return out


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def compute_one(traj_id: str) -> dict | None:
    out_dir = config.WEB_DATA / traj_id
    topo = out_dir / config.OUT_TOPOLOGY
    if not topo.exists():
        return None
    traj = md.load(str(topo))
    cas = _ca_coords_by_chain(traj)
    if len(cas) < 4:
        return None  # not a standard pHLA-TCR complex

    sizes = {ci: len(c) for ci, c in cas.items()}
    centroids = {ci: c.mean(axis=0) for ci, c in cas.items()}

    # peptide = shortest chain
    pep = min(sizes, key=lambda ci: sizes[ci])
    rest = [ci for ci in cas if ci != pep]
    # beta-2 microglobulin = chain whose length is closest to 100
    b2m = min(rest, key=lambda ci: abs(sizes[ci] - 100))
    rest2 = [ci for ci in rest if ci != b2m]
    # MHC heavy = of the remaining chains, the one closest to beta-2m
    mhc = min(rest2, key=lambda ci: np.linalg.norm(centroids[ci] - centroids[b2m]))
    tcr = [ci for ci in rest2 if ci != mhc]

    phla_atoms = np.concatenate([cas[mhc], cas[b2m], cas[pep]])
    tcr_atoms = np.concatenate([cas[ci] for ci in tcr])
    up = _unit(tcr_atoms.mean(axis=0) - phla_atoms.mean(axis=0))

    # peptide long axis (first principal component)
    p = cas[pep] - cas[pep].mean(axis=0)
    pep_axis = _unit(np.linalg.svd(p, full_matrices=False)[2][0])
    # horizontal axis orthogonal to `up`, then the viewing direction
    x = _unit(pep_axis - np.dot(pep_axis, up) * up)
    view = _unit(np.cross(x, up))

    # tight bounding box of the whole complex in the (x, up, view) view frame,
    # so the thumbnail renderer can frame the structure with minimal whitespace.
    # mdtraj coords are nm; report centre/extents in angstrom to match Mol*.
    allca = np.concatenate(list(cas.values()))
    ctr = allca.mean(axis=0)
    rel = allca - ctr
    hx, hu, hv = rel @ x, rel @ up, rel @ view

    def mid_half(a):
        return float((a.min() + a.max()) / 2), float((a.max() - a.min()) / 2)

    mx, half_w = mid_half(hx)
    mu, half_h = mid_half(hu)
    mv, _ = mid_half(hv)
    center = (ctr + mx * x + mu * up + mv * view) * 10.0

    return {
        "up": [float(v) for v in up],
        "view": [float(v) for v in view],
        "center": [float(v) for v in center],
        "half_w": half_w * 10.0,
        "half_h": half_h * 10.0,
    }


def main(argv: list[str]) -> int:
    ids = [a for a in argv if not a.startswith("--")]
    if not ids:
        ids = [d.name for d in sorted(config.WEB_DATA.glob("*")) if d.is_dir()]

    ok = skipped = 0
    for tid in tqdm(ids, unit="traj"):
        try:
            o = compute_one(tid)
        except Exception as exc:  # noqa: BLE001
            tqdm.write(f"FAILED {tid}: {exc}")
            o = None
        if o is None:
            skipped += 1
            continue
        with (config.WEB_DATA / tid / "orient.json").open("w") as fh:
            json.dump(o, fh)
        ok += 1
    print(f"\nOrientations: {ok} written, {skipped} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
