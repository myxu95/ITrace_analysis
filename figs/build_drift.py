#!/usr/bin/env python
"""Freeze the per-trajectory equilibration table that Figure 4 draws from.

The unit here is the TRAJECTORY, not the complex.  Every other frozen table in
this directory aggregates the three replicas of an entry, because the questions
they answer are about entries; this one asks whether each of the 735 deposited
trajectories had settled by the time it was stored, which is a question about
the trajectory itself and cannot survive being averaged with its siblings.

Three numbers per trajectory, all read from the served ``rmsd.json`` -- the
backbone RMSD against the energy-minimised starting structure, 1001 points over
200 ns:

  tail90_rmsd_nm   the mean over the last 90% of the run (20-200 ns).  The
                   first 10% is where a minimised structure relaxes into the
                   thermostat, so including it would measure the setup rather
                   than the sampling.  A mean, not a final value: one frame is
                   one frame.

  tail90_band_nm   the 5th-to-95th percentile band of the RMSD over that same
                   window: the range the value actually moves through once the
                   run has settled.  This is the quantity the figure reads on,
                   not the mean.  How far a structure sits from its own starting
                   coordinates is set as much by how big it is and by what the
                   fit was done on as by how it behaves; how wide a band it
                   wanders through, in the same nanometres, is a property of the
                   sampling.  A 5-95 band rather than max-min because over 901
                   points max-min is set by whichever single frame was most
                   extreme, and rather than a standard deviation because a band
                   is a range and reads directly off the same axis.

  slope_nm_per_ns  the least-squares slope over the final 20% (160-200 ns).
                   The mean above says how far the structure has moved; this
                   says whether it was still moving when we stopped.  They are
                   different questions and, in this library, uncorrelated ones,
                   which is why the figure plots them against each other rather
                   than reporting either alone.

The final 20% is a deliberately short window.  A slope fitted over the whole run
would be dominated by the early rise no matter how flat the end was, and would
report drift for every trajectory that equilibrated normally.

One provenance flag travels with the rows:

  dup_series   6g9q -- run1 and run2 carry bit-identical RMSD series.  This is
               the same defect the replica table records as ``dup_analysis``,
               seen here in the raw curve rather than in the derived analysis,
               and it means the release holds two independent runs of this
               complex and not three.  The pair is kept in the table and in the
               figure, because the figure is a census of what was deposited,
               but it must be disclosed: two of the 735 points coincide exactly.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
HERE = Path(__file__).resolve().parent
OUT = HERE / "drift.tsv"

TAIL_FRAC = 0.10        # leading fraction dropped before averaging
SLOPE_FRAC = 0.20       # trailing fraction the slope is fitted over
DUP_SERIES = {"6g9q"}   # complexes known to carry a duplicated RMSD series


def quant(s, p):
    """p-th quantile of the already-sorted ``s``, linearly interpolated."""
    i = p * (len(s) - 1)
    lo = int(i)
    return s[lo] if lo + 1 >= len(s) else s[lo] + (s[lo + 1] - s[lo]) * (i - lo)


def slope(t, r):
    """Least-squares slope of r on t, written out so the table has no deps."""
    n = len(t)
    mt = sum(t) / n
    mr = sum(r) / n
    sxx = sum((x - mt) ** 2 for x in t)
    return sum((x - mt) * (y - mr) for x, y in zip(t, r)) / sxx


def main():
    rows, seen, short = [], defaultdict(list), []
    for p in sorted(DATA.glob("*_run*/rmsd.json")):
        tid = p.parent.name
        pid, run = tid.rsplit("_run", 1)
        d = json.loads(p.read_text())
        t = [x / 1000.0 for x in d["time"]]          # ps as served -> ns
        r = d["rmsd"]                                 # nm
        if len(r) != 1001:
            short.append((tid, len(r)))
        i_tail = round(TAIL_FRAC * (len(r) - 1))
        i_slope = round((1.0 - SLOPE_FRAC) * (len(r) - 1))
        tail = r[i_tail:]
        srt = sorted(tail)
        seen[tuple(r)].append(tid)
        rows.append(dict(pdb_id=pid, run=int(run),
                         tail90_rmsd_nm=f"{sum(tail) / len(tail):.6f}",
                         tail90_band_nm=f"{quant(srt, 0.95) - quant(srt, 0.05):.6f}",
                         slope_nm_per_ns=f"{slope(t[i_slope:], r[i_slope:]):.8f}"))

    dup = {v[0].rsplit("_run", 1)[0] for v in seen.values() if len(v) > 1}
    for row in rows:
        row["dup_series"] = "dup_series" if row["pdb_id"] in dup else ""

    rows.sort(key=lambda r: (r["pdb_id"], r["run"]))
    with OUT.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)

    print(f"{OUT.name}   {len(rows)} trajectories, "
          f"{len({r['pdb_id'] for r in rows})} complexes")
    # a hand-written flag outlives the defect it names; check both directions
    if dup - DUP_SERIES:
        print(f"   !! dup_series NOT DECLARED: {sorted(dup - DUP_SERIES)}"
              f"  -- new occurrence, add it before trusting this table")
    if DUP_SERIES - dup:
        print(f"   !! dup_series NO LONGER PRESENT: {sorted(DUP_SERIES - dup)}"
              f"  -- repaired upstream? drop it from the set and redraw")
    for v in (v for v in seen.values() if len(v) > 1):
        print(f"   duplicated RMSD series: {v}")
    if short:
        print("   NOT 1001 POINTS:", short)


if __name__ == "__main__":
    main()
