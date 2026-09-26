"""Shared loader for rmsd_envelope.tsv and rmsd_exemplar.tsv, built by build_rmsdtrace.py.

Units are the ones the pipeline stores -- nm for the RMSD, ns for the time --
so the tables read the same as the served ``rmsd.json`` they were built from.
Conversion to Angstrom belongs to the figure.

The exemplar's identity is read out of its table rather than declared here, so
rebuilding with a different EXEMPLAR in build_rmsdtrace.py needs no edit on this
side and the figure cannot end up naming one complex while drawing another.
"""
import csv

from figstyle import HERE

TSV_ENV = HERE / "rmsd_envelope.tsv"
TSV_EX = HERE / "rmsd_exemplar.tsv"

N_TRAJ = 735
N_COMPLEX = 245
N_POINT = 1001
PCTS = [5, 25, 50, 75, 95]


def envelope():
    """-> (time_ns, {5: [...], 25: [...], 50: [...], 75: [...], 95: [...]}).

    Pointwise percentiles of all N_TRAJ trajectories at each of the N_POINT time
    points, nm.  Nothing is excluded: the band is a census of the release.
    """
    t, q = [], {p: [] for p in PCTS}
    for r in csv.DictReader(TSV_ENV.open(), delimiter="\t"):
        t.append(float(r["time_ns"]))
        for p in PCTS:
            q[p].append(float(r[f"p{p:02d}_nm"]))
    assert len(t) == N_POINT, len(t)
    return t, q


def exemplar():
    """-> (pdb_id, time_ns, {run: [rmsd_nm, ...]}) for the one drawn complex."""
    pid, t, runs = None, [], {}
    for r in csv.DictReader(TSV_EX.open(), delimiter="\t"):
        pid = r["pdb_id"]
        k = int(r["run"])
        if k not in runs:
            runs[k] = []
        runs[k].append(float(r["rmsd_nm"]))
        if k == 1:
            t.append(float(r["time_ns"]))
    assert sorted(runs) == [1, 2, 3], sorted(runs)
    assert all(len(v) == N_POINT for v in runs.values())
    assert len(t) == N_POINT
    return pid, t, runs


def smooth(v, k):
    """Centred running mean over +-k points, edge-padded so length is preserved.

    The stored series is one frame per 200 ps and the eye reads its high
    frequency as thickness rather than as signal; the smoothed line is drawn
    over the raw one, never instead of it.
    """
    n = len(v)
    pad = [v[0]] * k + list(v) + [v[-1]] * k
    return [sum(pad[i:i + 2 * k + 1]) / (2 * k + 1) for i in range(n)]


def quant(s, p):
    """p-th quantile of the already-sorted ``s``, linearly interpolated."""
    i = p * (len(s) - 1)
    lo = int(i)
    return s[lo] if lo + 1 >= len(s) else s[lo] + (s[lo + 1] - s[lo]) * (i - lo)
