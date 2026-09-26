"""Shared loader for the two dynamics tables Figure 3 draws from.

Built by build_dynamics.py.  Counting unit is the COMPLEX: a complex is
simulated three times from one starting structure, so each row here is already
the mean over that complex's replicas and every spread the figures show is
spread across complexes, never across replicas.
"""
import csv

from figstyle import HERE

BULGE_TSV = HERE / "dynamics_bulge.tsv"
PEP9_TSV = HERE / "dynamics_pep9.tsv"

N_TOTAL = 245
FOCUS_LENGTH = 9          # the only length numerous enough to read position by position

# absolute thresholds behind peptide_table.anchor, from pipeline/struct_metrics.py.
# They are fixed cut-offs on measured quantities, not quantiles of this library.
ANCHOR_TCR_MAX = 0.20     # occupancy: the TCR essentially never touches the residue
ANCHOR_SASA_MAX = 0.20    # nm^2: the residue is buried in the groove


def bulge():
    """-> [{pdb_id, peptide_length, bulge}] , one row per complex."""
    rows = []
    for r in csv.DictReader(BULGE_TSV.open(), delimiter="\t"):
        rows.append({"pdb_id": r["pdb_id"],
                     "peptide_length": int(r["peptide_length"]),
                     "bulge": float(r["bulge_angstrom"])})
    assert len(rows) == N_TOTAL, len(rows)
    return rows


def pep9():
    """-> {position: {field: [value per complex]}} for the 9-mer complexes."""
    out, seen = {}, set()
    for r in csv.DictReader(PEP9_TSV.open(), delimiter="\t"):
        p = int(r["position"])
        d = out.setdefault(p, {"tcr": [], "hla": [], "sasa": [], "pdb_id": []})
        d["tcr"].append(float(r["tcr_contact"]))
        d["hla"].append(float(r["hla_contact"]))
        d["sasa"].append(float(r["sasa_nm2"]))
        d["pdb_id"].append(r["pdb_id"])
        seen.add(r["pdb_id"])
    assert sorted(out) == list(range(1, FOCUS_LENGTH + 1)), sorted(out)
    assert all(len(v["tcr"]) == len(seen) for v in out.values())
    return out


def quant(values, p):
    """Linear-interpolated quantile; one definition for every Figure 3 box."""
    v = sorted(values)
    k = (len(v) - 1) * p
    f = int(k)
    c = min(f + 1, len(v) - 1)
    return v[f] + (v[c] - v[f]) * (k - f)
