"""Per-region RMSD relative to the MHC peptide-binding groove (bravo 1.2).

A single whole-complex RMSD hides the parts that matter for recognition. This
superposes every frame on the conserved MHC alpha1/alpha2 groove helices and
reports how much each functional region moves *relative to that platform*:

  * peptide backbone        -- does the peptide drift in the groove?
  * CDR3 loops (alpha+beta)  -- mobility of the specificity loops
  * all CDR loops            -- CDR1/2/3 together
  * TCR chains               -- TCR rigid-body + internal motion
  * interface (heavy atoms)  -- is the contact interface maintained?

Computed on the 200-frame web trajectory (protein-only, fast). The conventional
self-aligned whole-complex RMSD already lives in rmsd.json (GROMACS), so it is
not recomputed here.

Writes a sidecar ``analysis/rmsd_regions.json`` (re-embedded by extract_analysis)
and patches analysis.json directly.

    IMMUNO_WEB_DATA=.../immuno-dyn IMMUNO_SOURCE_ROOT=.../parallel3_preprocess \
        python -m pipeline.split_rmsd --source /home/xmy/immuno_analysis_run3
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import mdtraj as md

from . import config
from . import cdr_contacts

GROOVE_SUBREGIONS = {"alpha1_helix", "alpha2_helix"}
SERIES_POINTS = 100


def _atom_maps(top):
    """(chain_id, resSeq) -> CA index / backbone indices / heavy-atom indices."""
    ca: dict = {}
    bb: dict = {}
    heavy: dict = {}
    for atom in top.atoms:
        key = (atom.residue.chain.chain_id, atom.residue.resSeq)
        if atom.element.symbol != "H":
            heavy.setdefault(key, []).append(atom.index)
        if atom.name == "CA":
            ca[key] = atom.index
        if atom.name in ("N", "CA", "C", "O"):
            bb.setdefault(key, []).append(atom.index)
    return ca, bb, heavy


def _groove_residues(source_root: Path, traj_id: str) -> set:
    out = set()
    p = source_root / "rmsf" / traj_id / "analysis/rmsf/residue_rmsf.csv"
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("mhc_subregion") in GROOVE_SUBREGIONS:
                try:
                    out.add((r["chain_id"], int(r["resid"])))
                except (ValueError, KeyError):
                    pass
    return out


def _interface_residues(source_root: Path, traj_id: str) -> set:
    out = set()
    p = source_root / "contact" / traj_id / cdr_contacts.CONTACT_REL
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        for r in csv.DictReader(fh):
            for s in ("1", "2"):
                try:
                    out.add((r[f"chain_id_{s}"], int(r[f"resid_{s}"])))
                except (ValueError, KeyError):
                    pass
    return out


def _cdr_residues(source_root: Path, traj_id: str, meta: dict) -> dict:
    """{cdr_key: set((chain_id, resid))} unioned over alpha/beta chains."""
    out = {"cdr1": set(), "cdr2": set(), "cdr3": set()}
    tcr = meta.get("tcr") or {}
    summ = source_root / "contact" / traj_id / cdr_contacts.SUMMARY_REL
    if not summ.exists():
        return out
    topo = Path(json.loads(summ.read_text()).get("topology", ""))
    if not topo.exists():
        return out
    chains = cdr_contacts.parse_topology_chains(topo)
    for c in (tcr.get("chains") or {}).values():
        letter = c.get("structural_chain")
        if letter not in chains:
            continue
        cmap = cdr_contacts.cdr_resids_for_chain(chains[letter], c)
        for key in ("cdr1", "cdr2", "cdr3"):
            out[key] |= {(letter, rid) for rid in cmap.get(key, ())}
    return out


def _stats(series: np.ndarray) -> dict:
    return {"mean": round(float(series.mean()), 3),
            "std": round(float(series.std()), 3),
            "max": round(float(series.max()), 3)}


FNAT_CUTOFF_NM = 0.45   # 4.5 Å — matches the project's heavy-atom contact definition


def native_contacts(t, meta: dict, interface: set) -> dict | None:
    """Fnat / Q — does the crystallographic binding mode survive? (CAPRI/DockQ readout)

    Native contacts = inter-component (pMHC ↔ TCR) interface residue pairs in
    heavy-atom contact (< 4.5 Å) at frame 0. Fnat = fraction retained per frame;
    Q = soft Best–Hummer fraction (smooth over the cutoff). Pairs are drawn from
    the contact-derived interface residue set, split into pMHC vs TCR side.
    """
    tcr_chains = {c.get("structural_chain")
                  for c in (meta.get("tcr", {}).get("chains") or {}).values()}
    if not tcr_chains:
        return None
    res_index = {(r.chain.chain_id, r.resSeq): r.index for r in t.topology.residues}
    tcr_side = sorted({res_index[k] for k in interface if k[0] in tcr_chains and k in res_index})
    pmhc_side = sorted({res_index[k] for k in interface if k[0] not in tcr_chains and k in res_index})
    if not tcr_side or not pmhc_side:
        return None
    pairs = np.array([(p, q) for p in pmhc_side for q in tcr_side], dtype=np.int32)
    if pairs.shape[0] == 0:
        return None
    d, _ = md.compute_contacts(t, pairs, scheme="closest-heavy")   # (n_frames, n_pairs), nm
    native = d[0] < FNAT_CUTOFF_NM
    n_native = int(native.sum())
    if n_native < 3:
        return None
    dn = d[:, native]
    fnat = (dn < FNAT_CUTOFF_NM).mean(axis=1)
    r0 = dn[0]
    beta, lam = 50.0, 1.8                                          # nm⁻¹, unitless (Best–Hummer)
    q = (1.0 / (1.0 + np.exp(beta * (dn - lam * r0)))).mean(axis=1)

    duration = meta.get("duration_ns")
    n = t.n_frames
    time_ns = np.linspace(0.0, duration, n) if duration else np.arange(n, dtype=float)
    step = max(1, n // SERIES_POINTS)
    return {
        "n_native_contacts": n_native,
        "fnat_mean": round(float(fnat.mean()), 3),
        "fnat_min": round(float(fnat.min()), 3),
        "q_mean": round(float(q.mean()), 3),
        "series": {
            "time_ns": [round(float(x), 2) for x in time_ns[::step]],
            "fnat": [round(float(v), 3) for v in fnat[::step]],
            "q": [round(float(v), 3) for v in q[::step]],
        },
    }


def compute(traj_id: str, meta: dict, source_root: Path) -> dict | None:
    wd = config.WEB_DATA / traj_id
    top_path, traj_path = wd / "topology.pdb", wd / "traj.xtc"
    if not top_path.exists() or not traj_path.exists():
        return None
    pep_chain = meta.get("peptide_chain")
    tcr_chains = {c.get("structural_chain")
                  for c in (meta.get("tcr", {}).get("chains") or {}).values()}

    t = md.load(str(traj_path), top=str(top_path))
    if t.n_frames < 2:
        return None
    ca, bb, heavy = _atom_maps(t.topology)

    groove = _groove_residues(source_root, traj_id)
    groove_idx = [ca[k] for k in groove if k in ca]
    if len(groove_idx) < 4:
        return None  # no stable reference frame
    t.superpose(t, 0, atom_indices=np.array(groove_idx, dtype=np.int32))
    ref = t.xyz[0]

    def region_rmsd(indices):
        idx = np.array(sorted(set(indices)), dtype=np.int32)
        if idx.size == 0:
            return None
        diff = t.xyz[:, idx, :] - ref[idx]
        return np.sqrt((diff ** 2).sum(axis=2).mean(axis=1)) * 10.0  # nm -> Angstrom

    cdr = _cdr_residues(source_root, traj_id, meta)
    interface = _interface_residues(source_root, traj_id)

    # Approximate each TCR chain's variable domain as residues up to ~15 past its
    # CDR3 (FR4 end), so the constant domain's large lever-arm swing is excluded.
    v_cutoff: dict = {}
    for letter, rid in cdr["cdr3"]:
        v_cutoff[letter] = max(v_cutoff.get(letter, 0), rid + 15)

    regions = {
        "peptide_backbone": [i for k, idxs in bb.items() if k[0] == pep_chain for i in idxs],
        "cdr3_loops": [ca[k] for k in cdr["cdr3"] if k in ca],
        "cdr_loops": [ca[k] for key in ("cdr1", "cdr2", "cdr3") for k in cdr[key] if k in ca],
        "tcr_variable": [idx for k, idx in ca.items()
                         if k[0] in tcr_chains and k[1] <= v_cutoff.get(k[0], 0)],
        "tcr_chains": [idx for k, idx in ca.items() if k[0] in tcr_chains],
        "interface_heavy": [i for k, idxs in heavy.items() if k in interface for i in idxs],
    }

    duration = meta.get("duration_ns")
    n = t.n_frames
    time_ns = (np.linspace(0.0, duration, n) if duration else np.arange(n)).tolist()
    step = max(1, n // SERIES_POINTS)

    out = {"reference": "mhc_groove_helices", "n_frames": n, "n_groove_residues": len(groove_idx),
           "regions": {}}
    for name, idxs in regions.items():
        series = region_rmsd(idxs)
        if series is None:
            continue
        out["regions"][name] = {
            **_stats(series),
            "series": {
                "time_ns": [round(x, 2) for x in time_ns[::step]],
                "rmsd": [round(float(v), 3) for v in series[::step]],
            },
        }
    fnat = native_contacts(t, meta, interface)
    if fnat:
        out["fnat"] = fnat
    return out if (out["regions"] or out.get("fnat")) else None


def _iter_traj_dirs(source_root: Path):
    contact_root = source_root / "contact"
    if not contact_root.is_dir():
        return
    for d in sorted(contact_root.iterdir()):
        if d.is_dir():
            yield d.name


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--ids", default=None, help="Comma-separated traj ids (default: all in source).")
    args = ap.parse_args(argv)
    only = set(args.ids.split(",")) if args.ids else None

    n_ok = n_skip = 0
    for traj_id in _iter_traj_dirs(args.source):
        if only and traj_id not in only:
            continue
        meta_path = config.WEB_DATA / traj_id / config.OUT_META
        if not meta_path.exists():
            n_skip += 1
            continue
        meta = json.loads(meta_path.read_text())
        try:
            result = compute(traj_id, meta, args.source)
        except Exception as exc:  # never let one bad trajectory abort the batch
            print(f"[warn] {traj_id}: {exc}")
            result = None
        if result is None:
            n_skip += 1
            continue
        an_dir = config.WEB_DATA / traj_id / "analysis"
        an_dir.mkdir(parents=True, exist_ok=True)
        (an_dir / "rmsd_regions.json").write_text(json.dumps(result))
        an_json = an_dir / "analysis.json"
        if an_json.exists():
            data = json.loads(an_json.read_text())
            data["rmsd_regions"] = result
            if result.get("fnat"):
                data["fnat"] = result["fnat"]
            an_json.write_text(json.dumps(data))
        n_ok += 1

    print(f"Split RMSD: {n_ok} written, {n_skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
