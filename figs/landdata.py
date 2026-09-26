"""Shared loader for landscape.tsv and landscape_exemplar.json, built by build_landscape.py.

One row per trajectory, in the units the pipeline stores: variance fractions and
normalised mutual information are dimensionless, the free-energy grid is kcal/mol
above each surface's own minimum.

The vocabulary is fixed here and must not drift.  A BASIN is a local minimum of
the two-dimensional free-energy surface; it is what the coupling matrix is
indexed by.  The pipeline also stores a k-means SUBSTATE count over the same
scores, which is a different number, and neither this module nor anything drawn
from it carries that one.
"""
import csv
import json

from figstyle import HERE

TSV = HERE / "landscape.tsv"
JSON = HERE / "landscape_exemplar.json"

N_TRAJ = 735
N_COMPLEX = 245

# coupled_states.py, strongest first; "n/a" is a degenerate matrix, not a failure
VERDICTS = ["strong", "moderate", "weak", "independent", "n/a"]


def trajectories():
    """-> [{pdb_id, run, pep_basins, cdr3_basins, pep_pc1, ..., nmi, coupling, ...}].

    Where the matrix is degenerate -- one side resolved a single basin, so
    there is nothing for the mutual information to be computed between -- the
    pipeline writes ``nmi`` 0.0 and ``coupling`` "n/a", and ``reason`` names the
    side.  A zero there means "not defined", not "measured and found
    independent", which is why the verdict and not the number is what should be
    counted.  Those rows stay in: how many there are is a property of the
    release.
    """
    rows = [dict(pdb_id=r["pdb_id"], run=int(r["run"]),
                 pep_basins=int(r["pep_basins"]), cdr3_basins=int(r["cdr3_basins"]),
                 pep_pc1=float(r["pep_pc1"]), pep_pc2=float(r["pep_pc2"]),
                 cdr3_pc1=float(r["cdr3_pc1"]), cdr3_pc2=float(r["cdr3_pc2"]),
                 nmi=float(r["nmi"]) if r["nmi"] else None,
                 coupling=r["coupling"], reason=r["coupling_reason"],
                 reproducible=bool(r["replica_reproducible"]))
            for r in csv.DictReader(TSV.open(), delimiter="\t")]
    assert len(rows) == N_TRAJ, len(rows)
    assert len({r["pdb_id"] for r in rows}) == N_COMPLEX
    assert set(r["coupling"] for r in rows) <= set(VERDICTS)
    return rows


def exemplar():
    """-> the frozen plottable landscape of the one drawn trajectory."""
    d = json.loads(JSON.read_text())
    for side in ("peptide", "cdr3"):
        assert len(d[side]["fel"]["z"]) == len(d[side]["fel"]["y"])
        assert len(d[side]["fel"]["z"][0]) == len(d[side]["fel"]["x"])
    m = d["coupled"]["matrix"]
    assert len(m) == len(d["peptide"]["basins"])
    assert len(m[0]) == len(d["cdr3"]["basins"])
    return d


def quant(s, p):
    """p-th quantile of the already-sorted ``s``, linearly interpolated."""
    i = p * (len(s) - 1)
    lo = int(i)
    return s[lo] if lo + 1 >= len(s) else s[lo] + (s[lo + 1] - s[lo]) * (i - lo)
