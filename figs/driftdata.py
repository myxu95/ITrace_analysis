"""Shared loader for drift.tsv, built by build_drift.py.

One row per trajectory, in the units the pipeline stores: nm for the mean and
for the band, nm/ns for the slope.  Conversion to Angstrom belongs to the figure, not to the
table -- the table should read the same as the served data it was built from.
"""
import csv

from figstyle import HERE

TSV = HERE / "drift.tsv"

N_TRAJ = 735
N_COMPLEX = 245


def trajectories():
    """-> [{pdb_id, run, tail90, band, slope, dup}], one per deposited trajectory.

    Nothing is excluded.  The figure this feeds is a census of what was
    released, so the one complex whose two runs carry an identical RMSD series
    stays in and is disclosed by ``dup`` rather than being quietly dropped.
    """
    rows = [dict(pdb_id=r["pdb_id"], run=int(r["run"]),
                 tail90=float(r["tail90_rmsd_nm"]),
                 band=float(r["tail90_band_nm"]),
                 slope=float(r["slope_nm_per_ns"]),
                 dup=bool(r["dup_series"]))
            for r in csv.DictReader(TSV.open(), delimiter="\t")]
    assert len(rows) == N_TRAJ, len(rows)
    assert len({r["pdb_id"] for r in rows}) == N_COMPLEX
    return rows


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = sum((x - ma) ** 2 for x in a) ** 0.5
    db = sum((y - mb) ** 2 for y in b) ** 0.5
    return num / (da * db)
