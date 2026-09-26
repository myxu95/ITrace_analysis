"""Shared loader for composition.tsv (one row per complex, 245 rows).

Built by build_composition.py.  Counting unit is the COMPLEX, not the
trajectory: the three replicas of a complex are the same molecule.
"""
import csv

from figstyle import HERE

TSV = HERE / "composition.tsv"
N_TOTAL = 245

# display names; the archive taxonomy itself is unchanged
CAT_NAME = {
    "viral": "Viral",
    "tumor/self": "Tumour / self",
    "synthetic": "Synthetic",
    "murine": "Murine self",
    "bacterial": "Bacterial",
    "other": "Parasitic",
}
CAT_ORDER = ["viral", "tumor/self", "synthetic", "murine", "bacterial", "other"]


def load():
    rows = list(csv.DictReader(TSV.open(), delimiter="\t"))
    for r in rows:
        r["peptide_length"] = int(r["peptide_length"])
        r["release_year"] = int(r["release_year"])
    assert len(rows) == N_TOTAL, len(rows)
    return rows
