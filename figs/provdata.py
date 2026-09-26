"""Shared vocabulary and derived counts for the assay-provenance figure and table.

The four evidence classes are named in ONE place because they are printed in
two: ``fig2/fig_provenance_evidence.py`` draws them as a donut key and
``make_table_provenance.py`` prints them as four table rows.  Renaming a class
in only one of those puts a figure and a table in the same paper into open
disagreement, which is the kind of error no reviewer has to look for -- it
just reads as carelessness.  So the names live here and both scripts import
them.

The stratified share in ``era_tcr_share`` is computed rather than written down
for the same reason: it was previously hard-coded prose, it was labelled as
IEDB coverage when it is nothing of the kind, and it was derived from the
citation year, which 28 of the 245 complexes do not have -- they were dropped
in silence.  PDB release year is complete for all 245, so the split below is
over the whole library.
"""
import csv
from collections import Counter

import compdata
from figstyle import HERE

PROV = HERE.parents[1] / "Scientific_data" / "provenance"
N_TOTAL = 245

# (key in provenance_master.assay_evidence_level, display name), strongest first.
#
# One class is assigned per complex, strongest match wins (merge_provenance.py:23),
# so the four are mutually exclusive and exhaustive and their shares add to 100%.
#
# All four names have the same grammatical shape -- a noun phrase whose only
# varying part is WHAT the retrieved assay matched: the receptor itself, the
# epitope, or just the paper.  That is the whole meaning of the axis, so it
# belongs in the labels rather than in a legend the reader has to go find.  The
# two "only" qualifiers carry the exclusivity for the same reason.
CLASSES = [
    ("TCR-level",       "Assay for this TCR"),
    ("epitope-level",   "Assay for the epitope only"),
    ("literature-only", "Assay in the publication only"),
    ("none",            "No assay retrieved"),
]
KEYS = [k for k, _ in CLASSES]
NAME = dict(CLASSES)

EARLY, LATE = 2010, 2019       # released up to EARLY / from LATE onward


def master():
    rows = list(csv.DictReader(open(PROV / "provenance_master.csv", encoding="utf-8-sig")))
    assert len(rows) == N_TOTAL, len(rows)
    return rows


def levels(rows=None):
    rows = rows if rows is not None else master()
    lvl = Counter(r["assay_evidence_level"] for r in rows)
    assert sum(lvl[k] for k in KEYS) == N_TOTAL, lvl     # the four really do partition
    return lvl


def anchored(rows=None):
    """Complexes with a retrievable assay of any kind -> (n, percent)."""
    rows = rows if rows is not None else master()
    n = sum(1 for r in rows if r["assay_evidence_level"] != "none")
    return n, 100.0 * n / N_TOTAL


def n_kd(rows=None):
    """Complexes carrying at least one quantitative KD."""
    rows = rows if rows is not None else master()
    return sum(1 for r in rows if int(r["n_KD_values"] or 0) > 0)


def era_tcr_share():
    """Share with an assay for the TCR itself, early era vs late -> (pct, pct).

    Stratified on PDB release year, which every complex has, so nothing is
    dropped.  This is what shows the weakest class to be a curation artefact:
    the receptor-level share collapses with recency while the structures
    themselves did not get any less studied.
    """
    year = {r["pdb_id"]: r["release_year"] for r in compdata.load()}
    lvl = {r["pdb_id"]: r["assay_evidence_level"] for r in master()}
    assert set(year) == set(lvl), "provenance and composition disagree on the id set"
    out = []
    for lo, hi in ((0, EARLY), (LATE, 9999)):
        sub = [p for p, y in year.items() if lo <= y <= hi]
        out.append(100.0 * sum(1 for p in sub if lvl[p] == "TCR-level") / len(sub))
    return tuple(out)


if __name__ == "__main__":
    rows = master()
    lvl = levels(rows)
    for k in KEYS:
        print(f"  {NAME[k]:32s} {lvl[k]:4d}  {100 * lvl[k] / N_TOTAL:5.1f}%")
    n, pct = anchored(rows)
    print(f"\n  anchored {n}/{N_TOTAL} ({pct:.1f}%), quantitative KD {n_kd(rows)}")
    print("  TCR-level share  <=%d: %.1f%%   >=%d: %.1f%%" % (EARLY, era_tcr_share()[0],
                                                              LATE, era_tcr_share()[1]))
