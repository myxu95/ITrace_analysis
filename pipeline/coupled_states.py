"""Coupled states — the joint occupancy matrix of peptide backbone conformational states
(``peptide_dpca`` basins) × TCR CDR3 loop states (``tcr_cdr3_dpca`` basins), for the same
trajectory frames. Quantifies whether the two sides of the recognition interface move in a
CORRELATED way (conformational selection / induced fit) or independently.

Per frame we already have — after the peptide + CDR3 dPCAs — a peptide state and a CDR3
state; the joint histogram M[peptide, cdr3] (co-occupancy), its marginals, the per-cell
enrichment over independence M/(p_i·c_j), and the normalised mutual information between the
two state variables (0 = independent, 1 = fully locked) summarise the coupling.

Depends on ``peptide_dpca`` + ``tcr_cdr3_dpca`` already being in analysis.json; recomputes the
per-frame scores self-contained and assigns each frame to the nearest stored basin centre.
Writes sidecar ``analysis/coupled_states.json`` (re-embedded under ``coupled_states``) and
patches analysis.json. A second pass adds cross-replica reproducibility.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.coupled_states [--ids a,b]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md

from . import config
from .peptide_dihedrals import _peptide_chain
from .tcr_cdr3 import _res_atoms, _tcr_struct_chain

SIDECAR = "coupled_states.json"


def _scores(sub) -> np.ndarray | None:
    """2D dPCA scores for a residue subset — identical feature construction to the peptide /
    CDR3 dPCA (backbone φ/ψ + side-chain χ1/χ2 + Cα distances, each MFA-normalised)."""
    _, phi = md.compute_phi(sub)
    _, psi = md.compute_psi(sub)
    try:
        _, chi1 = md.compute_chi1(sub)
    except Exception:   # noqa: BLE001
        chi1 = None
    try:
        _, chi2 = md.compute_chi2(sub)
    except Exception:   # noqa: BLE001
        chi2 = None
    ca = sub.topology.select("name CA")
    dpair = None
    if ca.size >= 3:
        ii, jj = np.triu_indices(ca.size, k=2)
        if ii.size:
            pos = sub.xyz[:, ca, :]
            dpair = np.linalg.norm(pos[:, ii, :] - pos[:, jj, :], axis=2)

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
    Xc = np.concatenate(parts, axis=1)
    if Xc.shape[1] < 2:
        return None
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    return U[:, :2] * S[:2]


def _assign(scores: np.ndarray, basins: list) -> np.ndarray | None:
    """Voronoi-assign each frame to the nearest stored basin centre (b.pc1, b.pc2)."""
    cen = np.array([[b["pc1"], b["pc2"]] for b in basins if b.get("pc1") is not None], dtype=float)
    if cen.size == 0 or scores is None:
        return None
    return ((scores[:, None, :] - cen[None, :, :]) ** 2).sum(axis=2).argmin(axis=1)


def _verdict(nmi: float) -> str:
    if nmi < 0.05:
        return "independent"
    if nmi < 0.20:
        return "weak"
    if nmi < 0.40:
        return "moderate"
    return "strong"


def compute(traj_id: str, meta: dict) -> dict | None:
    wd = config.WEB_DATA / traj_id
    tp, xp = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
    an_json = wd / "analysis" / "analysis.json"
    if not (tp.exists() and xp.exists() and an_json.exists()):
        return None
    a = json.loads(an_json.read_text())
    pd, cd = a.get("peptide_dpca") or {}, a.get("tcr_cdr3_dpca") or {}
    p_bas, c_bas = pd.get("basins") or [], cd.get("basins") or []
    if not p_bas or not c_bas:
        return None
    tc = (a.get("tcr_cdr") or {}).get("cdr_resids") or {}
    a_cdr3 = (tc.get("alpha") or {}).get("cdr3") or []
    b_cdr3 = (tc.get("beta") or {}).get("cdr3") or []
    if not a_cdr3 or not b_cdr3:
        return None

    t = md.load(str(xp), top=str(tp))
    if t.n_frames < 10:
        return None
    top0 = t.topology

    pep = _peptide_chain(meta, top0)
    pep_atoms = top0.select(f"chainid {next((c.index for c in top0.chains if c.chain_id == pep), -1)}") if pep else np.array([])
    ach, bch = _tcr_struct_chain(meta, "alpha"), _tcr_struct_chain(meta, "beta")
    if pep_atoms.size == 0 or not ach or not bch:
        return None
    cdr3_atoms = sorted(_res_atoms(top0, ach, a_cdr3) + _res_atoms(top0, bch, b_cdr3))
    if len(cdr3_atoms) < 20:
        return None

    p_asn = _assign(_scores(t.atom_slice(pep_atoms)), p_bas)
    c_asn = _assign(_scores(t.atom_slice(np.asarray(cdr3_atoms))), c_bas)
    if p_asn is None or c_asn is None:
        return None

    nP, nC = len(p_bas), len(c_bas)
    M = np.zeros((nP, nC))
    for p, c in zip(p_asn, c_asn):
        M[p, c] += 1
    n_frames = int(M.sum())
    if n_frames < 10:
        return None
    M /= M.sum()
    pm, cm = M.sum(axis=1), M.sum(axis=0)
    eps = 1e-12
    enr = M / (np.outer(pm, cm) + eps)
    # normalised mutual information (0 = independent, 1 = one state fully predicts the other)
    mi = float(sum(M[i, j] * np.log((M[i, j] + eps) / (pm[i] * cm[j] + eps))
                   for i in range(nP) for j in range(nC) if M[i, j] > 0))
    Hp = float(-sum(p * np.log(p) for p in pm if p > 0))
    Hc = float(-sum(c * np.log(c) for c in cm if c > 0))
    degenerate = nP < 2 or nC < 2 or min(Hp, Hc) <= eps
    nmi = 0.0 if degenerate else mi / min(Hp, Hc)
    # strongest coupled pair = the cell most enriched over independence (with real occupancy)
    ti, tj = int(np.argmax(np.where(M >= 0.05, enr, 0))) // nC, int(np.argmax(np.where(M >= 0.05, enr, 0))) % nC
    top_pair = {"peptide": ti, "cdr3": tj, "occ": round(float(M[ti, tj]), 4),
                "enrichment": round(float(enr[ti, tj]), 2)}

    return {
        "n_peptide": nP, "n_cdr3": nC, "n_frames": n_frames,
        "matrix": [[round(float(v), 4) for v in row] for row in M],
        "peptide_pop": [round(float(v), 4) for v in pm],
        "cdr3_pop": [round(float(v), 4) for v in cm],
        "enrichment": [[round(float(v), 2) for v in row] for row in enr],
        "nmi": round(nmi, 3),
        "coupling": "n/a" if degenerate else _verdict(nmi),
        "reason": ("peptide single-state" if nP < 2 else "CDR3 single-state") if degenerate else None,
        "top_pair": top_pair,
        "replica_agreement": None,   # filled by the cross-replica second pass in main()
    }


def _write(an_dir: Path, res: dict) -> None:
    an_dir.mkdir(parents=True, exist_ok=True)
    (an_dir / SIDECAR).write_text(json.dumps(res))
    an_json = an_dir / "analysis.json"
    if an_json.exists():
        data = json.loads(an_json.read_text())
        data["coupled_states"] = res
        an_json.write_text(json.dumps(data))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all in WEB_DATA).")
    args = ap.parse_args(argv)
    only = set(args.ids.split(",")) if args.ids else None

    written: dict[str, tuple[Path, dict]] = {}
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
        written[tid] = (meta_path.parent / "analysis", res)
        n_ok += 1

    # cross-replica agreement: group by complex (pdb prefix), flag coupling as reproducible if
    # a majority of the complex's replicas are at least weakly coupled.
    by_complex: dict[str, list[str]] = {}
    for tid in written:
        by_complex.setdefault(tid.split("_run")[0], []).append(tid)
    for tids in by_complex.values():
        vals = [written[t][1] for t in tids]
        n = len(vals)
        n_coupled = sum(1 for v in vals if v["coupling"] not in ("n/a", "independent"))
        agree = {"n": n, "n_coupled": n_coupled, "reproducible": n_coupled * 2 >= n and n_coupled > 0}
        for v in vals:
            v["replica_agreement"] = agree

    for _tid, (an_dir, res) in written.items():
        _write(an_dir, res)
    print(f"Coupled states: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
