"""TCR CDR3 conformational landscape — the recognition-critical, flexible loops of the TCR,
mirroring the peptide dPCA (``pipeline/peptide_dihedrals.py``). Combines CDR3α + CDR3β into
ONE alignment-free dPCA (backbone φ/ψ + side-chain χ1/χ2 + intra/inter-loop Cα distances),
builds a −kT·ln P free-energy landscape, detects conformational basins, and saves
framework-superposed representative structures so the detail page can overlay them and read
which loop residue / rotamer flips between states.

Superposing on the Vα/Vβ FRAMEWORK (the stable Ig β-sandwich, CDRs excluded) isolates the
loops' INTERNAL conformational change from the rigid-body TCR docking (already covered by the
crossing / incident angles). A rigid, preorganised CDR3 → high specificity; a plastic one →
adaptability / cross-reactivity — the TCR-side complement to the peptide plasticity card.

Reads the served ``analysis.json``'s ``tcr_cdr.cdr_resids`` (ANARCI/IMGT → author resSeq,
already computed); only reliable annotations carrying both α/β CDR3 loops are processed.
Writes a sidecar ``analysis/tcr_cdr3_dpca.json`` (re-embedded under ``tcr_cdr3_dpca``) and
patches analysis.json directly.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.tcr_cdr3 [--ids a,b]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md

from . import config
from .peptide_dihedrals import _kde_grid, _detect_basins, _rotamer, _kmeans, GRID, KMAX

SIDECAR = "tcr_cdr3_dpca.json"
FW_MARGIN = 12   # V-domain framework cutoff = cdr3_end + this (resSeq); excludes the C domain


def _tcr_struct_chain(meta: dict, which: str) -> str | None:
    """Structural chain letter for the alpha/beta TCR chain. Prefer meta.tcr.chains
    (canonical, same source as cdr_resids), fall back to chain_roles."""
    chains = (meta.get("tcr") or {}).get("chains") or {}
    c = (chains.get(which) or {}).get("structural_chain")
    if c:
        return c
    greek = {"alpha": "α", "beta": "β"}[which]
    for cid, role in (meta.get("chain_roles") or {}).items():
        rl = str(role)
        if "TCR" in rl and (greek in rl or which in rl.lower()):
            return cid
    return None


def _res_atoms(top, chain_letter, resseqs) -> list[int]:
    rs = set(resseqs)
    return [a.index for ch in top.chains if ch.chain_id == chain_letter
            for r in ch.residues if r.resSeq in rs for a in r.atoms]


def _res_ca(top, chain_letter, keep_pred) -> list[int]:
    out = []
    for ch in top.chains:
        if ch.chain_id != chain_letter:
            continue
        for r in ch.residues:
            if keep_pred(r.resSeq):
                ca = next((a.index for a in r.atoms if a.name == "CA"), None)
                if ca is not None:
                    out.append(ca)
    return out


def _ref_atoms(top, chain_letter, keep_pred) -> list[int]:
    return [a.index for ch in top.chains if ch.chain_id == chain_letter
            for r in ch.residues if keep_pred(r.resSeq) for a in r.atoms]


def compute(traj_id: str, meta: dict) -> dict | None:
    wd = config.WEB_DATA / traj_id
    tp, xp = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
    an_json = wd / "analysis" / "analysis.json"
    if not (tp.exists() and xp.exists() and an_json.exists()):
        return None
    tc = (json.loads(an_json.read_text()).get("tcr_cdr") or {})
    # NB: tc["reliable"] gates the CONTACT decomposition (which CDR reads peptide vs MHC — the
    # recognition-summary ratios), NOT the CDR3 residue positions. The cdr_resids come straight
    # from ANARCI/IMGT and are correct even when the contact reliability is off (verified: for
    # the flagged run3 replicas they select exactly the ANARCI CDR3 sequence). So we do NOT gate
    # the CDR3 conformational landscape on it — only require a resolvable CDR3 on both chains.
    contact_reliable = bool(tc.get("reliable"))
    cdr = tc.get("cdr_resids") or {}
    a_cdr = cdr.get("alpha") or {}
    b_cdr = cdr.get("beta") or {}
    a_cdr3, b_cdr3 = a_cdr.get("cdr3") or [], b_cdr.get("cdr3") or []
    if not a_cdr3 or not b_cdr3:
        return None
    ach, bch = _tcr_struct_chain(meta, "alpha"), _tcr_struct_chain(meta, "beta")
    if not ach or not bch:
        return None

    t = md.load(str(xp), top=str(tp))
    if t.n_frames < 10:
        return None
    top0 = t.topology

    # CDR3α + CDR3β atoms (author resSeq); sorted → topology order (α chain precedes β).
    cdr3_atoms = sorted(_res_atoms(top0, ach, a_cdr3) + _res_atoms(top0, bch, b_cdr3))
    if len(cdr3_atoms) < 20:
        return None
    n_alpha_res = sum(1 for ch in top0.chains if ch.chain_id == ach
                      for r in ch.residues if r.resSeq in set(a_cdr3))

    # Framework superpose: Vα/Vβ β-sandwich Cα (resSeq ≤ cdr3_end + margin, ALL CDRs excluded)
    # so the loops' internal motion isn't aligned out; frame 0 is the reference.
    def _fw_ca(chain, cd):
        allcdr = set(cd.get("cdr1") or []) | set(cd.get("cdr2") or []) | set(cd.get("cdr3") or [])
        cut = max(cd.get("cdr3") or [0]) + FW_MARGIN
        return _res_ca(top0, chain, lambda s: s <= cut and s not in allcdr)
    fw_ca = _fw_ca(ach, a_cdr) + _fw_ca(bch, b_cdr)
    framework_ref = None
    if len(fw_ca) >= 8:
        t.superpose(t, 0, atom_indices=np.asarray(fw_ca))

    sub = t.atom_slice(np.asarray(cdr3_atoms))
    phi_idx, phi = md.compute_phi(sub)
    psi_idx, psi = md.compute_psi(sub)
    try:
        chi1_idx, chi1 = md.compute_chi1(sub)
    except Exception:   # noqa: BLE001
        chi1_idx, chi1 = None, None
    try:
        chi2_idx, chi2 = md.compute_chi2(sub)
    except Exception:   # noqa: BLE001
        chi2_idx, chi2 = None, None
    # Cα–Cα distances across BOTH loops (|i−j| ≥ 2): intra-loop shape + inter-loop arrangement
    dpair, ca = None, sub.topology.select("name CA")
    if ca.size >= 3:
        ii, jj = np.triu_indices(ca.size, k=2)
        if ii.size:
            pos = sub.xyz[:, ca, :]
            dpair = np.linalg.norm(pos[:, ii, :] - pos[:, jj, :], axis=2)

    # three MFA-balanced blocks (identical treatment to the peptide dPCA)
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
    n_chi1 = int(chi1.shape[1]) if chi1 is not None else 0
    n_chi2 = int(chi2.shape[1]) if chi2 is not None else 0
    n_dist = int(dpair.shape[1]) if dpair is not None else 0
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2
    var_frac = (var / var.sum()).tolist() if var.sum() > 0 else [0, 0]
    scores = U[:, :2] * S[:2]
    if scores.shape[1] < 2:
        return None
    pc1, pc2 = scores[:, 0], scores[:, 1]

    pad1 = (pc1.max() - pc1.min()) * 0.08 or 0.1
    pad2 = (pc2.max() - pc2.min()) * 0.08 or 0.1
    gx = np.linspace(pc1.min() - pad1, pc1.max() + pad1, GRID)
    gy = np.linspace(pc2.min() - pad2, pc2.max() + pad2, GRID)
    try:
        z = _kde_grid(pc1, pc2, gx, gy)
    except (np.linalg.LinAlgError, ValueError):
        return None
    xc = [round(float(v), 5) for v in gx]   # 5 dp: coarse rounding degenerates the grid → jagged
    yc = [round(float(v), 5) for v in gy]

    try:
        basins = _detect_basins(pc1, pc2, gx, gy, z)
    except Exception:   # noqa: BLE001
        basins = []

    an_dir = wd / "analysis"
    an_dir.mkdir(parents=True, exist_ok=True)
    for old in list(an_dir.glob("cdr3_basin_*.pdb")) + [an_dir / "cdr3_framework_ref.pdb"]:
        try:
            old.unlink()
        except OSError:
            pass
    for k, b in enumerate(basins):
        try:
            sub[b["frame"]].save_pdb(str(an_dir / f"cdr3_basin_{k}.pdb"))
            b["pdb"] = f"cdr3_basin_{k}.pdb"
        except Exception:   # noqa: BLE001
            pass
    # framework reference for the 3D context: V-domain of both chains EXCEPT CDR3 (keeps the
    # β-sandwich + CDR1/2 scaffold; CDR3 is drawn separately, coloured by basin).
    try:
        ref_idx = []
        for chain, cd in ((ach, a_cdr), (bch, b_cdr)):
            cut = max(cd.get("cdr3") or [0]) + FW_MARGIN
            cdr3set = set(cd.get("cdr3") or [])
            ref_idx += _ref_atoms(top0, chain, lambda s, c=cut, x=cdr3set: s <= c and s not in x)
        if ref_idx:
            t[0].atom_slice(np.asarray(sorted(ref_idx))).save_pdb(str(an_dir / "cdr3_framework_ref.pdb"))
            framework_ref = "cdr3_framework_ref.pdb"
    except Exception:   # noqa: BLE001
        framework_ref = None

    # per-basin side-chain fingerprint (χ1/χ2 + rotamer), labelled by loop + resSeq
    stop = sub.topology
    chi1_res = ({stop.atom(int(chi1_idx[j][0])).residue.index: j for j in range(chi1_idx.shape[0])}
                if chi1 is not None and chi1.shape[1] else {})
    chi2_res = ({stop.atom(int(chi2_idx[j][0])).residue.index: j for j in range(chi2_idx.shape[0])}
                if chi2 is not None and chi2.shape[1] else {})
    for b in basins:
        fr = int(b["frame"])
        sc = []
        for rr in sorted(chi1_res):
            ro = stop.residue(rr)
            c1 = float(np.degrees(chi1[fr, chi1_res[rr]]))
            c2 = float(np.degrees(chi2[fr, chi2_res[rr]])) if rr in chi2_res else None
            sc.append({"loop": "α" if rr < n_alpha_res else "β", "resid": int(ro.resSeq),
                       "resname": ro.name, "chi1": round(c1, 1),
                       "chi2": round(c2, 1) if c2 is not None else None,
                       "rot1": _rotamer(c1), "rot2": _rotamer(c2)})
        b["sidechain"] = sc

    labels = _kmeans(scores, KMAX)
    pops = sorted((float((labels == c).mean()) for c in range(KMAX)), reverse=True)
    pops = [round(p, 3) for p in pops if p > 0]
    n_states = sum(1 for p in pops if p >= 0.10)

    return {
        "chains": {"alpha": ach, "beta": bch},
        "cdr3_len": {"alpha": len(a_cdr3), "beta": len(b_cdr3)},
        "contact_reliable": contact_reliable,   # whether the CDR CONTACT decomposition was reliable
        "pc1_var_frac": round(float(var_frac[0]), 3),
        "pc2_var_frac": round(float(var_frac[1]), 3),
        "features": {"n_phi": int(phi.shape[1]), "n_psi": int(psi.shape[1]),
                     "n_chi1": n_chi1, "n_chi2": n_chi2, "n_ca_dist": n_dist},
        "n_substates": int(n_states),
        "substate_populations": pops,
        "fel": {"z": z, "x": xc, "y": yc,
                "x_label": "PC1", "y_label": "PC2", "unit": "kcal/mol"},
        "basins": basins,
        "framework_ref": framework_ref,
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
            data["tcr_cdr3_dpca"] = res
            an_json.write_text(json.dumps(data))
        n_ok += 1
    print(f"TCR CDR3 dPCA: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
