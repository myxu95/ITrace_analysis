"""Interface interaction fingerprint — the pMHC↔TCR interface split into the five
specific interaction TYPES, per residue pair, over the trajectory.

For every inter-partner residue pair (pMHC side = MHC/peptide/β2m vs TCR side) this
detects, per frame, whether the pair forms a **hydrogen bond, salt bridge,
hydrophobic contact, π–π stack or cation–π interaction**, using standard PLIP-style
geometric criteria on the served HEAVY-ATOM trajectory (topology.pdb + traj.xtc).
It reports, per pair per type, the occupancy (fraction of frames present) and the
on/off frame segments — so the Explore/Detail interface viewer can filter contacts
by type and show each type's dynamics. Self-contained: mdtraj only, no source CSVs.

The served trajectory is protein heavy atoms only (no explicit H), so H-bonds use a
heavy-atom donor–acceptor distance criterion (no D–H···A angle); the other four
types do not require hydrogens.

    IMMUNO_WEB_DATA=.../web_data python -m pipeline.interface_interactions [--ids a,b]
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import mdtraj as md

from . import config
from . import chain_roles as _chain_roles

# ---- geometric cutoffs (Å) — PLIP-style, adapted for a heavy-atom trajectory ----
CUT = {"hbond": 3.5, "saltbridge": 4.0, "hydrophobic": 4.0, "pipi": 5.5, "cationpi": 6.0}
PRUNE_CA_A = 16.0     # ignore residue pairs whose Cα never come within this
MIN_OCC = 0.02        # drop pair/type present in <2% of frames (single-frame noise)

_PMHC_ROLES = {"MHC", "peptide", "β2m"}
_TCR_ROLES = {"TCRα", "TCRβ"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
AROMATIC_RING = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}
CATION_ATOMS = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",), "HIS": ("ND1", "NE2")}
ANION_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
_ROLE_SIDE = None


def _segments(mask: np.ndarray) -> list[list[int]]:
    """Boolean per-frame presence -> list of [start, end] inclusive frame segments."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    segs, s = [], idx[0]
    for a, b in zip(idx[:-1], idx[1:]):
        if b != a + 1:
            segs.append([int(s), int(a)]); s = b
    segs.append([int(s), int(idx[-1])])
    return segs


def _residue_features(traj, roles_by_chain, chain_letter):
    """Per interface residue -> dict of feature-atom index lists + side + label."""
    feats = {}
    for res in traj.topology.residues:
        letter = chain_letter(res.chain.index)
        role = roles_by_chain.get(letter)
        side = "pmhc" if role in _PMHC_ROLES else "tcr" if role in _TCR_ROLES else None
        if side is None:
            continue
        rn = res.name
        by = {a.name: a.index for a in res.atoms}
        ca = by.get("CA")
        if ca is None:
            continue
        no = [a.index for a in res.atoms if a.element is not None and a.element.symbol in ("N", "O")]
        cat = [by[n] for n in CATION_ATOMS.get(rn, ()) if n in by]
        ani = [by[n] for n in ANION_ATOMS.get(rn, ()) if n in by]
        phob = [a.index for a in res.atoms
                if rn in HYDROPHOBIC_RES and a.element is not None and a.element.symbol == "C"
                and a.name not in ("C", "CA")]
        ring = [by[n] for n in AROMATIC_RING.get(rn, ()) if n in by]
        feats[res.index] = {
            "side": side, "ca": ca, "no": no, "cat": cat, "ani": ani,
            "phob": phob, "ring": ring if len(ring) >= 4 else [],
            "label": f"{rn}{res.resSeq}", "chain": letter, "resid": res.resSeq, "role": role,
        }
    return feats


def _min_dist_pairs(traj, pair_atom_index):
    """Given list of (respair_key, [(ai,aj),...]) build one flat atom-pair array,
    compute distances (nm) for all frames, and return {respair_key: min-dist(frame)}."""
    flat, owner = [], []
    for key, aps in pair_atom_index:
        for ai, aj in aps:
            flat.append((ai, aj)); owner.append(key)
    if not flat:
        return {}
    d = md.compute_distances(traj, np.array(flat, dtype=np.int32)) * 10.0   # nm -> Å
    cols_by_key = defaultdict(list)
    for c, key in enumerate(owner):
        cols_by_key[key].append(c)
    return {key: d[:, cols].min(axis=1) for key, cols in cols_by_key.items()}


def _centroid(xyz, idx_list):
    return xyz[:, idx_list, :].mean(axis=1) * 10.0   # (F,3) in Å


def compute(traj_id: str, meta: dict) -> dict | None:
    wd = config.WEB_DATA / traj_id
    top_path, traj_path = wd / config.OUT_TOPOLOGY, wd / config.OUT_TRAJ
    if not (top_path.exists() and traj_path.exists()):
        return None
    traj = md.load(str(traj_path), top=str(top_path))
    F = traj.n_frames
    if F < 2:
        return None
    roles = _chain_roles.role_map(meta) or {}
    if not roles:
        return None
    # map mdtraj chain.index -> chain letter via meta.chains order (file order)
    order = [c.get("id") for c in (meta.get("chains") or [])]
    chain_letter = (lambda i: order[i] if i < len(order) else None)
    feats = _residue_features(traj, roles, chain_letter)
    pmhc = [ri for ri, f in feats.items() if f["side"] == "pmhc"]
    tcr = [ri for ri, f in feats.items() if f["side"] == "tcr"]
    if not pmhc or not tcr:
        return None

    # prune inter-partner residue pairs by Cα distance (any frame within PRUNE_CA_A)
    ca_p = np.array([feats[r]["ca"] for r in pmhc]); ca_t = np.array([feats[r]["ca"] for r in tcr])
    xyz = traj.xyz
    cand = []
    for i, rp in enumerate(pmhc):
        dca = np.linalg.norm(xyz[:, ca_p[i], :][:, None, :] - xyz[:, ca_t, :], axis=2).min(axis=0) * 10.0
        for j, rt in enumerate(tcr):
            if dca[j] <= PRUNE_CA_A:
                cand.append((rp, rt))
    if not cand:
        return None

    # per-type atom/centroid distance -> per-pair present-mask
    present = {t: {} for t in CUT}   # type -> {pairkey: bool-mask(F)}

    def _atom_type(t, get_atoms_a, get_atoms_b):
        pai = []
        for rp, rt in cand:
            aa, bb = get_atoms_a(feats[rp]), get_atoms_b(feats[rt])
            if aa and bb:
                pai.append(((rp, rt), [(x, y) for x in aa for y in bb]))
        for key, dmin in _min_dist_pairs(traj, pai).items():
            m = dmin <= CUT[t]
            if m.any():
                present[t][key] = m

    # symmetric atom sets (a pair can donate/accept from either side)
    _atom_type("hbond", lambda f: f["no"], lambda f: f["no"])
    _atom_type("hydrophobic", lambda f: f["phob"], lambda f: f["phob"])
    # salt bridge: cation(a)-anion(b) OR anion(a)-cation(b)
    for ga, gb in ((lambda f: f["cat"], lambda f: f["ani"]), (lambda f: f["ani"], lambda f: f["cat"])):
        pai = []
        for rp, rt in cand:
            aa, bb = ga(feats[rp]), gb(feats[rt])
            if aa and bb:
                pai.append(((rp, rt), [(x, y) for x in aa for y in bb]))
        for key, dmin in _min_dist_pairs(traj, pai).items():
            m = dmin <= CUT["saltbridge"]
            if m.any():
                present["saltbridge"][key] = present["saltbridge"].get(key, np.zeros(F, bool)) | m

    # π–π: aromatic ring centroids within cutoff
    cent = {ri: _centroid(xyz, feats[ri]["ring"]) for ri in feats if feats[ri]["ring"]}
    for rp, rt in cand:
        if rp in cent and rt in cent:
            d = np.linalg.norm(cent[rp] - cent[rt], axis=1)
            m = d <= CUT["pipi"]
            if m.any():
                present["pipi"][(rp, rt)] = m
    # cation–π: cation atoms vs aromatic ring centroid, both directions
    for rp, rt in cand:
        for cat_res, ring_res in ((rp, rt), (rt, rp)):
            cats, ring = feats[cat_res]["cat"], (cent.get(ring_res))
            if cats and ring is not None:
                catxyz = xyz[:, cats, :] * 10.0
                d = np.linalg.norm(catxyz - ring[:, None, :], axis=2).min(axis=1)
                m = d <= CUT["cationpi"]
                if m.any():
                    key = (rp, rt)
                    present["cationpi"][key] = present["cationpi"].get(key, np.zeros(F, bool)) | m

    # assemble per-pair records
    pairs = {}
    for t, d in present.items():
        for (rp, rt), mask in d.items():
            occ = float(mask.mean())
            if occ < MIN_OCC:
                continue
            key = (rp, rt)
            rec = pairs.get(key)
            if rec is None:
                fa, fb = feats[rp], feats[rt]
                rec = pairs[key] = {
                    "a_chain": fa["chain"], "a_resid": fa["resid"], "a_label": fa["label"], "a_role": fa["role"],
                    "b_chain": fb["chain"], "b_resid": fb["resid"], "b_label": fb["label"], "b_role": fb["role"],
                    "types": {},
                }
            rec["types"][t] = {"occ": round(occ, 3), "seg": _segments(mask)}
    out = sorted(pairs.values(),
                 key=lambda r: -max(v["occ"] for v in r["types"].values()))
    counts = {t: sum(1 for r in out if t in r["types"]) for t in CUT}
    return {"n_frames": F, "cutoffs_angstrom": CUT, "counts": counts, "pairs": out}


def _write(traj_id: str, result: dict) -> None:
    an = config.WEB_DATA / traj_id / "analysis"
    an.mkdir(parents=True, exist_ok=True)
    (an / "interactions_pairs.json").write_text(json.dumps(result))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all).")
    args = ap.parse_args(argv)
    web = config.WEB_DATA
    metas = sorted(web.glob("*/" + config.OUT_META))
    if args.ids:
        want = {i.strip() for i in args.ids.split(",")}
        metas = [m for m in metas if m.parent.name in want]
    n_ok = n_skip = 0
    for mp in metas:
        tid = mp.parent.name
        try:
            meta = json.loads(mp.read_text())
            res = compute(tid, meta)
        except Exception as exc:   # noqa: BLE001
            print(f"[warn] {tid}: {exc}"); res = None
        if res is None:
            n_skip += 1; continue
        _write(tid, res); n_ok += 1
    print(f"Interface interactions: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
