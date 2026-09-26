"""Essential-dynamics honesty gate: PC1 cosine content + cross-replica RMSIP.

Two trust metrics that condition every dynamics conclusion in the database:

  * **PC1 cosine content** (Hess 2002) — if the largest principal component's
    projection looks like a half-cosine (content ≈ 1), the 200 ns is mostly random
    diffusion, NOT converged sampling. Per trajectory.
  * **Subspace RMSIP** — root-mean-square inner product of the top-10 Cartesian
    Cα PCA eigenvectors between a complex's independent replicas. High RMSIP ⇒ the
    replicas explored the *same* essential motions (reproducible dynamics). Per
    complex (needs ≥2 replicas with a matching Cα count).

Purely trajectory-derived (topology.pdb + traj.xtc); globs WEB_DATA and processes
one complex (all its replicas) at a time so RMSIP can be computed. Writes a
per-trajectory sidecar ``analysis/essential_dynamics.json`` (re-embedded by
extract_analysis) and patches analysis.json directly.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.essential_dynamics [--ids a,b]
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import mdtraj as md

from . import config

N_MODES = 10            # top eigenvectors compared for RMSIP
SIDECAR = "essential_dynamics.json"

# Trust thresholds, PINNED to the 735-trajectory web_data distribution (2026-07-01):
#   pc1_cosine_content  median 0.445, p90 0.765  -> >0.7 flags the diffusion-like tail (~17%)
#   subspace_rmsip      min 0.599, median 0.760  -> <0.7 flags partial overlap (~5%); <0.5 never
#                       occurs in-set (kept as a defensive floor).
COSINE_DIFFUSIVE = 0.7          # PC1 cosine content above this ⇒ PC1 resembles random diffusion (Hess 2002)
RMSIP_HI, RMSIP_LO = 0.7, 0.5   # cross-replica essential-subspace overlap bands


def _ca_pca(traj_path: Path, top_path: Path):
    """Top-N Cartesian Cα PCA modes + PC1 cosine content for one trajectory.

    Returns (pc1_cosine_content, pc1_variance_fraction, n_ca, eigvecs[3N, N_MODES])
    or None. Eigenvectors are mass-unweighted Cartesian (3N) modes, descending.
    """
    t = md.load(str(traj_path), top=str(top_path))
    if t.n_frames < 4:
        return None
    ca = t.topology.select("name CA")
    if ca.size < 8:
        return None
    sub = t.atom_slice(ca)
    sub.superpose(sub, 0)
    F = sub.n_frames
    disp = (sub.xyz - sub.xyz.mean(axis=0)).reshape(F, -1)   # (F, 3N)
    cov = disp.T @ disp / F                                   # (3N, 3N)
    evals, evecs = np.linalg.eigh(cov)                        # ascending
    order = evals[::-1]
    vecs = evecs[:, ::-1]                                     # columns = modes, descending
    var_frac = float(order[0] / order.sum()) if order.sum() > 0 else 0.0
    # PC1 cosine content (Hess 2002), discrete form, scale-invariant
    proj = disp @ vecs[:, 0]
    k = np.arange(F)
    denom = float(proj @ proj)
    cc = (2.0 / F) * (np.cos(np.pi * k / F) @ proj) ** 2 / denom if denom > 0 else None
    return (round(float(cc), 3) if cc is not None else None,
            round(var_frac, 3), int(ca.size), vecs[:, :N_MODES])


def _rmsip(a: np.ndarray, b: np.ndarray) -> float:
    """RMSIP of two mode sets (3N × k): sqrt(mean over i,j of (a_i·b_j)^2) ∈ [0,1]."""
    m = a.T @ b                       # (k, k) inner products
    return float(np.sqrt((m ** 2).sum() / m.shape[0]))


def process_complex(traj_ids: list[str]) -> dict:
    """Compute per-replica cosine content + cross-replica RMSIP for one complex."""
    per = {}
    for tid in traj_ids:
        wd = config.WEB_DATA / tid
        tp, xp = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
        if not (tp.exists() and xp.exists()):
            continue
        try:
            res = _ca_pca(xp, tp)
        except Exception as exc:   # noqa: BLE001
            print(f"[warn] {tid}: {exc}")
            res = None
        if res:
            cc, vf, n_ca, vecs = res
            per[tid] = {"cosine": cc, "var_frac": vf, "n_ca": n_ca, "vecs": vecs}

    # pairwise RMSIP among replicas that share a Cα count
    by_n = defaultdict(list)
    for tid, d in per.items():
        by_n[d["n_ca"]].append(tid)
    rmsip_of = {tid: [] for tid in per}
    for tids in by_n.values():
        for i in range(len(tids)):
            for j in range(i + 1, len(tids)):
                r = _rmsip(per[tids[i]]["vecs"], per[tids[j]]["vecs"])
                rmsip_of[tids[i]].append(r)
                rmsip_of[tids[j]].append(r)

    out = {}
    for tid, d in per.items():
        rs = rmsip_of[tid]
        out[tid] = {
            "pc1_cosine_content": d["cosine"],
            "pc1_variance_fraction": d["var_frac"],
            "subspace_rmsip": round(float(np.mean(rs)), 3) if rs else None,
            "n_replicas_compared": len(rs) + 1 if rs else 1,
            "n_modes": N_MODES,
        }
    return out


def dynamics_trust(pc1_cosine, rmsip, cdr_reliable) -> dict:
    """Aggregate dynamics-trust verdict for one trajectory.

    Rolls the three honesty signals — PC1 cosine content, cross-replica RMSIP,
    and CDR-decomposition reliability — into a single
    {"level": "high"|"medium"|"low", "reasons": [...]}. A missing input (stage not
    run / single replica) adds an 'unknown' reason but is NOT counted as a demerit,
    so absence of evidence never masquerades as low trust.
    """
    reasons, demerits = [], 0
    if pc1_cosine is not None and pc1_cosine > COSINE_DIFFUSIVE:
        reasons.append(f"PC1 cosine content {pc1_cosine:.2f} — dominant motion resembles random "
                       "diffusion, not converged sampling"); demerits += 1
    if rmsip is None:
        reasons.append("no cross-replica RMSIP (single replica or mismatched Cα)")
    elif rmsip < RMSIP_LO:
        reasons.append(f"cross-replica RMSIP {rmsip:.2f} — replicas explore different essential "
                       "motions"); demerits += 2
    elif rmsip < RMSIP_HI:
        reasons.append(f"cross-replica RMSIP {rmsip:.2f} — only partial overlap of replica motions")
        demerits += 1
    if cdr_reliable is False:
        reasons.append("CDR contact decomposition unreliable — chain/recognition ratios excluded")
        demerits += 1
    level = "high" if demerits == 0 else "medium" if demerits <= 2 else "low"
    return {"level": level, "reasons": reasons}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None,
                    help="Comma-separated traj ids; expanded to all replicas of their complexes.")
    args = ap.parse_args(argv)

    # group every trajectory by complex (pdb_id)
    complexes = defaultdict(list)
    for meta_path in sorted(config.WEB_DATA.glob("*/" + config.OUT_META)):
        tid = meta_path.parent.name
        pdb = json.loads(meta_path.read_text()).get("pdb_id") or tid.rsplit("_", 1)[0]
        complexes[pdb].append(tid)

    if args.ids:
        want = {i.strip() for i in args.ids.split(",")}
        pdbs = {tid.rsplit("_", 1)[0] for tid in want}
        complexes = {p: v for p, v in complexes.items() if p in pdbs}

    n_ok = n_skip = 0
    for pdb, tids in complexes.items():
        result = process_complex(sorted(tids))
        for tid, block in result.items():
            an_dir = config.WEB_DATA / tid / "analysis"
            an_dir.mkdir(parents=True, exist_ok=True)
            (an_dir / SIDECAR).write_text(json.dumps(block))
            an_json = an_dir / "analysis.json"
            if an_json.exists():
                data = json.loads(an_json.read_text())
                data["essential_dynamics"] = block
                an_json.write_text(json.dumps(data))
            n_ok += 1
        n_skip += len(tids) - len(result)
    print(f"Essential dynamics: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
