#!/usr/bin/env python
"""Freeze the basis-free replica-agreement table: how far apart are an entry's
three replicas, measured in Angstrom rather than in anyone's choice of basis?

Figure 4's landscape panels draw one entry's three replicas on one free-energy
surface.  That surface is a projection and it cannot settle the question it
raises.  Each served landscape is expressed in its OWN dPCA basis -- the
eigenvectors are not stored and the component signs are arbitrary -- so the three
cannot simply be laid on top of one another; and even after they are re-derived
in a shared basis (build_pooled_fel.py) the two components carry about a quarter
of the variance, and distance in that plane tracks true RMSD only weakly
(Spearman +0.2 to +0.5 over the entries tested).  This table is the measurement
the picture cannot make, and it is made in Angstrom.

Every trajectory's basins are served as all-atom PDBs -- ``basin_k.pdb`` for the
peptide, ``cdr3_basin_k.pdb`` for the CDR3 loops -- one representative frame
each.  Pool an entry's basins over its three replicas, take every pairwise
heavy-atom RMSD after optimal superposition, and ask two questions of each basin:

  d_self   how far is the closest basin of the SAME replica?  This is the scale
           on which one trajectory is itself willing to call two conformations
           different -- the resolution the run works at, set by the data.
  d_cross  how far is the closest basin of a DIFFERENT replica?  How far you
           have to go to find the same conformation in another run.

MINIMUM AGAINST MINIMUM.  An earlier version of this table compared d_cross with
the median over all within-replica PAIRS, and that comparison is broken however
the data come out: d_cross is a minimum over roughly five candidates and a
minimum is smaller than a median of the same distribution by construction.  It
made the replicas look closer to each other than to themselves.  They are not.
Both columns here are minima, so the ratio means something.

The third question is a shape rather than a scale:

  hit      the basin's overall nearest neighbour is in another replica
  exp      the number of such hits expected if the replica labels were shuffled,
           summed as len(other) / (n - 1) over basins -- about 0.74, because with
           three replicas most of the other basins belong to another run

hit well below exp says each replica's basins cluster with their own siblings:
the runs carry a small systematic offset from one another rather than
interleaving.  Both facts belong in the table.  The scale of the offset is what
matters for reuse, and it is what the ratio column reports.

near_w is d_cross weighted by basin population, so the wells that actually hold
the frames dominate; near_max is the worst basin in the entry, kept because a
median can hide one replica that went somewhere the others never did.
within_pair is the old all-pairs within-replica median, kept because it is the
spread of one run's basins and is worth having -- but it is NOT the comparator
for d_cross, and nothing should pair the two.

Superposition is mdtraj's, over heavy atoms only: the served representatives are
written from trajectories aligned on different things (the peptide on the MHC
groove, the CDR3 loops on the Valpha/Vbeta framework), so the comparison must do
its own optimal fit rather than trust the stored frame.

ok=0 rows are entries whose basin sets do not share an atom count -- the known
run3 chain-mislabelling defect, not a conformational statement.  They are written
out rather than dropped so the count stays visible.

Read-only against the release.  Writes replica_basins.tsv.
"""
import csv
import itertools
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, "/tmp/mdtraj_pkg")           # mdtraj, kept out of the conda envs

import numpy as np
import mdtraj as md

R = Path(os.environ.get("IMMUNO_WEB_DATA",
                        "/home/xmy/work/data/immunotrace/web_data_1000"))
HERE = Path(__file__).resolve().parent
OUT = HERE / "replica_basins.tsv"

SIDES = {"peptide": ("peptide_dpca", "basin_"),
         "cdr3": ("tcr_cdr3_dpca", "cdr3_basin_")}
RUNS = (1, 2, 3)


def _load(pid, side):
    """-> (run labels, populations, basin structures) pooled over the replicas."""
    key, pre = SIDES[side]
    runs, pops, trs = [], [], []
    for run in RUNS:
        d = R / f"{pid}_run{run}" / "analysis"
        if not (d / "analysis.json").exists():
            continue
        try:
            an = json.loads((d / "analysis.json").read_text())
        except Exception:
            continue
        for k, b in enumerate((an.get(key) or {}).get("basins") or []):
            p = d / f"{pre}{k}.pdb"
            if p.exists():
                runs.append(run)
                pops.append(float(b.get("pop", 0.0)))
                trs.append(md.load(str(p)))
    return np.array(runs), np.array(pops), trs


def entry(pid, side):
    runs, pops, trs = _load(pid, side)
    if len(set(runs.tolist())) < 2:
        return None
    blank = dict(pdb_id=pid, side=side, ok=0, n_rep=len(set(runs.tolist())),
                 n_basins="", d_self="", d_cross="", near_w="", near_max="",
                 within_pair="", hit="", exp="")
    if len({t.n_atoms for t in trs}) != 1:
        return blank

    n = len(trs)
    heavy = trs[0].topology.select("not element H")
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        D[i, j] = D[j, i] = float(md.rmsd(trs[j], trs[i], 0,
                                          atom_indices=heavy)[0]) * 10.0

    d_self, d_cross, wt, hit, exp = [], [], [], 0, 0.0
    for i in range(n):
        own = [j for j in range(n) if j != i and runs[j] == runs[i]]
        oth = [j for j in range(n) if runs[j] != runs[i]]
        if own:
            d_self.append(D[i, own].min())
        d_cross.append(D[i, oth].min())
        wt.append(pops[i])
        nb = min([j for j in range(n) if j != i], key=lambda j: D[i, j])
        hit += int(runs[nb] != runs[i])
        exp += len(oth) / (n - 1.0)

    dc = np.array(d_cross)
    w = np.array(wt)
    within = [D[i, j] for i, j in itertools.combinations(range(n), 2)
              if runs[i] == runs[j]]
    return dict(pdb_id=pid, side=side, ok=1, n_rep=len(set(runs.tolist())), n_basins=n,
                d_self=f"{np.median(d_self):.4f}" if d_self else "",
                d_cross=f"{np.median(dc):.4f}",
                near_w=f"{float((dc * w).sum() / (w.sum() or 1.0)):.4f}",
                near_max=f"{dc.max():.4f}",
                within_pair=f"{np.median(within):.4f}" if within else "",
                hit=hit, exp=f"{exp:.4f}")


def main():
    ids = sorted({p.name.rsplit("_run", 1)[0] for p in R.iterdir() if "_run" in p.name})
    rows = [r for pid in ids for side in SIDES
            for r in [entry(pid, side)] if r is not None]

    with OUT.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)

    print(f"{OUT.name}   {len(rows)} rows, {len({r['pdb_id'] for r in rows})} complexes")
    q = lambda a, p: np.percentile(a, p)
    for side in SIDES:
        s = [r for r in rows if r["side"] == side and r["ok"] and r["d_self"]]
        bad = [r for r in rows if r["side"] == side and not r["ok"]]
        ds = np.array([float(r["d_self"]) for r in s])
        dc = np.array([float(r["d_cross"]) for r in s])
        rat = dc / ds
        H = sum(int(r["hit"]) for r in rows if r["side"] == side and r["ok"])
        E = sum(float(r["exp"]) for r in rows if r["side"] == side and r["ok"])
        N = sum(int(r["n_basins"]) for r in rows if r["side"] == side and r["ok"])
        print(f"  {side:8s} {len(s)} usable, {len(bad)} atom-count mismatch")
        print(f"    closest own-run basin     median {np.median(ds):.2f} A  "
              f"IQR {q(ds, 25):.2f}-{q(ds, 75):.2f}")
        print(f"    closest other-run basin   median {np.median(dc):.2f} A  "
              f"IQR {q(dc, 25):.2f}-{q(dc, 75):.2f}")
        print(f"    ratio cross/self          median {np.median(rat):.2f}  "
              f"IQR {q(rat, 25):.2f}-{q(rat, 75):.2f};  "
              f"{int((rat <= 1).sum())} of {rat.size} at or below 1")
        print(f"    nearest neighbour in another run: {H}/{N} = {H / N:.1%}  "
              f"(shuffled: {E / N:.1%})")


if __name__ == "__main__":
    main()
