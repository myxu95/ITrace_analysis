#!/usr/bin/env python
"""Freeze the two RMSD series Figure 4 draws: the library envelope and one exemplar.

The other frozen tables in this directory reduce each trajectory to a scalar.
This one keeps the curve, because the question is about shape over time -- does
the RMSD rise and then stop rising -- and no scalar carries that.

Two tables come out of one pass over the served ``rmsd.json`` files (the
backbone RMSD against the energy-minimised starting structure, 1001 points over
200 ns, nm as stored):

  rmsd_envelope.tsv   the pointwise 5th, 25th, 50th, 75th and 95th percentile
                      of all 735 trajectories at each of the 1001 time points.
                      Not a mean and a standard deviation: the distribution is
                      right-skewed at every time point -- a handful of large,
                      mobile complexes sit far above a tight majority -- so a
                      symmetric interval would put its lower edge below curves
                      that exist and its upper edge below curves that also
                      exist.  Percentiles are read off the same axis as the
                      data and make no claim about the shape of the spread.

  rmsd_exemplar.tsv   the three replicas of one complex, long format, one row
                      per point.  The exemplar's identity travels in the table
                      rather than in a constant in the figure, so the two can
                      never disagree.

Why an exemplar at all.  A per-trajectory summary (see build_drift.py) answers
"did they settle" across the whole library but cannot show what one entry looks
like, and a user downloads one entry.  Drawing a single complex alone would
invite the obvious objection that it was chosen for looking good; drawing it on
top of the library's own percentile band answers that objection inside the
figure, since the reader can see where this entry sits among the 735.

EXEMPLAR is a deliberate, disclosed choice, made on published numbers and not on
appearance: 1nam's three replicas all sit near the library median on both
equilibration axes, and all three carry a well-resolved landscape on both sides
of the interface, so the same entry can serve this figure and the companion
landscape figure.  Changing it is a one-line edit here plus a rebuild.

Provenance note: 6g9q run1 and run2 carry bit-identical RMSD series, a known
defect of that entry's stored analyses (see build_drift.py).  Both stay in the
envelope -- it is a census of what was released, and one duplicated series
among 735 moves no percentile -- but the exemplar must never be 6g9q, and the
build asserts that.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

DATA = Path("/home/xmy/work/data/immunotrace/web_data_1000")
HERE = Path(__file__).resolve().parent
OUT_ENV = HERE / "rmsd_envelope.tsv"
OUT_EX = HERE / "rmsd_exemplar.tsv"

EXEMPLAR = "1nam"
PCTS = [5, 25, 50, 75, 95]
NPT = 1001
DUP_SERIES = {"6g9q"}


def quant(s, p):
    """p-th quantile of the already-sorted ``s``, linearly interpolated."""
    i = p * (len(s) - 1)
    lo = int(i)
    return s[lo] if lo + 1 >= len(s) else s[lo] + (s[lo + 1] - s[lo]) * (i - lo)


def main():
    assert EXEMPLAR not in DUP_SERIES, EXEMPLAR
    series, grid, short = {}, None, []
    for p in sorted(DATA.glob("*_run*/rmsd.json")):
        tid = p.parent.name
        d = json.loads(p.read_text())
        if len(d["rmsd"]) != NPT:
            short.append((tid, len(d["rmsd"])))
            continue
        t = [x / 1000.0 for x in d["time"]]          # ps as served -> ns
        if grid is None:
            grid = t
        elif t != grid:
            raise SystemExit(f"{tid}: time grid differs from the rest of the library")
        series[tid] = d["rmsd"]                       # nm

    ids = sorted(series)
    n = len(ids)
    with OUT_ENV.open("w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["time_ns"] + [f"p{q:02d}_nm" for q in PCTS])
        for i, tv in enumerate(grid):
            col = sorted(series[k][i] for k in ids)
            wr.writerow([f"{tv:.1f}"] + [f"{quant(col, q / 100.0):.6f}" for q in PCTS])

    runs = sorted(int(k.rsplit("_run", 1)[1]) for k in ids
                  if k.rsplit("_run", 1)[0] == EXEMPLAR)
    assert runs == [1, 2, 3], (EXEMPLAR, runs)
    with OUT_EX.open("w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["pdb_id", "run", "time_ns", "rmsd_nm"])
        for r in runs:
            for tv, v in zip(grid, series[f"{EXEMPLAR}_run{r}"]):
                wr.writerow([EXEMPLAR, r, f"{tv:.1f}", f"{v:.6f}"])

    print(f"{OUT_ENV.name}   {len(grid)} time points, {n} trajectories, "
          f"{len({k.rsplit('_run', 1)[0] for k in ids})} complexes")
    print(f"{OUT_EX.name}   exemplar {EXEMPLAR}, runs {runs}")
    dup = sorted({a.rsplit("_run", 1)[0] for i, a in enumerate(ids)
                  for b in ids[i + 1:] if series[a] == series[b]
                  and a.rsplit("_run", 1)[0] == b.rsplit("_run", 1)[0]})
    if set(dup) - DUP_SERIES:
        print(f"   !! duplicated series NOT DECLARED: {sorted(set(dup) - DUP_SERIES)}")
    if DUP_SERIES - set(dup):
        print(f"   !! declared duplicate NO LONGER PRESENT: {sorted(DUP_SERIES - set(dup))}")
    if short:
        print("   NOT 1001 POINTS, excluded:", short)


if __name__ == "__main__":
    main()
