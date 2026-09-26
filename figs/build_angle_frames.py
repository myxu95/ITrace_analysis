#!/usr/bin/env python
"""Pool the docking angles of a complex over every frame, not the plot series.

``analysis.json`` stores the angle as a 150-point downsample for the web plot.
Figure 3 pools the three replicas of a complex and quotes p5/p50/p95 of that
pool, and until now it pooled the downsample -- 3 x 150 points described as a
3 x 1001-frame sample in three places.  The angle is slow enough that the two
agree to the printed precision, but the figure states a sampling basis, so it
should use the one it states.

This recomputes the angle at full frame resolution straight from the published
trajectory, with the same per-frame geometry the pipeline uses (imported, not
copied, so the two cannot drift apart), and freezes the pooled percentiles in
``angle_percentiles.tsv``.  ``build_replicas.py`` reads that file; it is a
separate step because the recompute is ~20 min over 735 trajectories while the
rest of the table build is seconds.

Run it AFTER any trajectory repair lands in web_data_1000 -- it reads
``traj.xtc``, so a stale download file produces a stale percentile.

Columns: pdb_id, n_replicas, n_frames (pooled), inc_p5/p50/p95,
cross_p5/p50/p95.  Angles in degrees; crossing is the directed 0-180 value
(>90 = reverse polarity), matching ``angle.crossing_deg``.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/data/work/Immuno-Dyn")
from pipeline import config, docking_angle as da  # noqa: E402

HERE = Path(__file__).resolve().parent
COMP = HERE / "composition.tsv"
OUT = HERE / "angle_percentiles.tsv"


def quant(values, p):
    """Same definition as repdata.quant, so the two tables cannot disagree."""
    v = sorted(values)
    k = (len(v) - 1) * p
    f = int(k)
    c = min(f + 1, len(v) - 1)
    return v[f] + (v[c] - v[f]) * (k - f)


def series(traj_id):
    """-> (incident, directed_crossing) at full frame resolution, or raise."""
    import mdtraj as md

    d = config.WEB_DATA / traj_id
    meta = json.loads((d / "meta.json").read_text())
    pep_ch = meta.get("peptide_chain") or da._MANUAL_PEPTIDE_CHAIN.get(meta.get("pdb_id"))
    tc = (meta.get("tcr") or {}).get("chains") or {}
    a_ch = (tc.get("alpha") or {}).get("structural_chain")
    b_ch = (tc.get("beta") or {}).get("structural_chain")
    mhc_ch = da._mhc_chain(meta, pep_ch, a_ch, b_ch)
    if not (pep_ch and a_ch and b_ch and mhc_ch):
        raise RuntimeError(f"{traj_id}: chain roles incomplete")

    t = md.load(str(d / config.OUT_TRAJ), top=str(d / config.OUT_TOPOLOGY))
    top = t.topology

    def cas(chid):
        return [i for i in da._chain_ca(top, chid) if i is not None]

    mhc_ca = cas(mhc_ch)[:da.MHC_PLATFORM]
    va, vb = cas(a_ch)[:da.TCR_VARIABLE], cas(b_ch)[:da.TCR_VARIABLE]
    pep = cas(pep_ch)
    if len(mhc_ca) < 20 or len(va) < 20 or len(vb) < 20 or len(pep) < 2:
        raise RuntimeError(f"{traj_id}: too few CA atoms to define the geometry")

    cross = np.empty(t.n_frames)
    inc = np.empty(t.n_frames)
    pol = np.empty(t.n_frames)
    for f in range(t.n_frames):
        cross[f], inc[f], pol[f] = da._frame_angles(t.xyz[f], mhc_ca, va, vb, pep)
    # The fold is per FRAME, not per trajectory: docking_angle.compute reports a
    # trajectory-level reversed_polarity vote but transforms each frame by its own
    # sign, so a trajectory that changes polarity mid-run carries both branches.
    # (pol: +1 forward, -1 reversed, 0 = undefined near-perpendicular, left acute.)
    directed = np.where(pol < 0, 180.0 - cross, cross)
    return inc, directed


def main():
    rows = []
    ids = [r["pdb_id"] for r in csv.DictReader(COMP.open(), delimiter="\t")]
    for n, pid in enumerate(ids, 1):
        runs = sorted(p.parent.name for p in config.WEB_DATA.glob(f"{pid}_run*/meta.json"))
        if not runs:
            raise RuntimeError(f"{pid}: no trajectories in {config.WEB_DATA}")
        inc, cross = [], []
        for tid in runs:
            i, c = series(tid)
            inc.append(i)
            cross.append(c)
        inc = np.concatenate(inc)
        cross = np.concatenate(cross)
        row = {"pdb_id": pid, "n_replicas": len(runs), "n_frames": inc.size}
        for tag, v in (("inc", inc), ("cross", cross)):
            for p in (5, 50, 95):
                row[f"{tag}_p{p}"] = f"{quant(v.tolist(), p / 100):.3f}"
        rows.append(row)
        print(f"[{n}/{len(ids)}] {pid} {len(runs)} runs {inc.size} frames "
              f"inc {row['inc_p5']}/{row['inc_p50']}/{row['inc_p95']}", flush=True)

    with OUT.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    print(f"{OUT.name} {len(rows)} complexes, "
          f"{sum(r['n_replicas'] for r in rows)} trajectories")


if __name__ == "__main__":
    main()
