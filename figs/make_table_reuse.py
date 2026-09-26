#!/usr/bin/env python
"""Receptor reuse: the entries that share one T-cell receptor.

Grouped by what the set holds fixed and what it varies, because that is the
question a reader brings to it - a within-receptor peptide series is a different
object from a receptor seen twice on the same complex.  Everything is read from
Scientific_data/provenance/receptor_reuse.tsv, which build_reuse_table.py
generates from the clonotype assignment, so the table cannot drift from it.
"""
from pathlib import Path
import csv, collections

from figstyle import *          # canvas, type scale, ink, save

# The tables keep the pre-2026-09 type scale.  The figures were enlarged because
# they are reduced on placement; these tables are placed at their drawn width,
# and every column here is measured in mm against text set at 7 pt.
FS_BODY, FS_SMALL, FS_HEAD = 7.0, 6.5, 7.0
F_BODY, F_SMALL, F_HEAD = font(FS_BODY), font(FS_SMALL), font(FS_HEAD, bold=True)

PROV = Path(__file__).resolve().parents[2] / "Scientific_data" / "provenance"
SRC = PROV / "receptor_reuse.tsv"

# achromatic, like Table 1 and Table_lineage: grouping is carried by rules and
# weight.  Hue in this manuscript encodes data, and nothing here is data.
F_NAME, F_ID = F_BODY, font(FS_SMALL, mono=True)

# column geometry, in mm from the left text edge
# the two numeric columns are set by their headings, not their digits: "entries"
# and "peptides" are 8.0 and 9.9 mm wide against 2.4 mm of number.
X_MHC, X_N, X_NPEP, X_PDB = 26.5, 58.5, 70.5, 73.0
LH, GAP_BAND = 3.55, 3.1

# the order the groups are printed in: the long peptide series first, because
# that is what most of the reuse is, and the repeat determinations last
GROUPS = [
    ("peptide", "One receptor, a series of peptides"),
    ("peptide and MHC substitution",
     "One receptor, a series of peptides, and the MHC altered as well"),
    ("peptide and MHC allele", "One receptor, peptide and allele both varied"),
    ("MHC allele", "One receptor, one peptide, more than one allele"),
    ("MHC substitution",
     "One receptor, one peptide, one allele, MHC altered in the groove"),
    ("neither", "One receptor and one complex, determined more than once"),
]


def alleles(field):
    """'HLA-B*35:01 HLA-B*35:08' -> 'HLA-B*35:01, :08'.

    Only the trailing field is elided, and only when everything before it is
    shared; HLA-B*08:01 and HLA-B*44:05 stay written out, because there the
    locus and the first field are the thing that differs.
    """
    a = field.split()
    if len(a) == 1:
        return a[0]
    stem = a[0].rsplit(":", 1)[0]
    if all(x.startswith(stem + ":") for x in a[1:]):
        return a[0] + ", " + ", ".join(":" + x.rsplit(":", 1)[1] for x in a[1:])
    # the locus prefix is written once; B*44:05 after HLA-B*08:01 is unambiguous
    return a[0] + ", " + ", ".join(x.replace("HLA-", "", 1) for x in a[1:])


def load():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8"), delimiter="\t"))
    grp = collections.OrderedDict((k, []) for k, _ in GROUPS)
    for r in rows:
        grp[r["varies"]].append(r)
    for v in grp.values():
        v.sort(key=lambda r: (-int(r["n_structures"]), r["receptor_name"]))
    return rows, grp


LEGEND = {   # prose harvested by build_legends.py; see figstyle.LEAN
    "title": "Receptors that appear in more than one entry",
    "subtitle": "Forty-one of the 156 receptors in the library are present in two or more "
                "complexes. Each is a series in which the receptor is held fixed and the "
                "peptide, the allele, the MHC construct, or nothing, is varied.",
    "notes": [],
}


def draw(ax, top, w):
    rows, grp = load()
    ns = sum(int(r["n_structures"]) for r in rows)
    big = [r for r in rows if int(r["n_structures"]) >= 5]
    var = [r["receptor_name"] for r in rows if r["lineage_role"] == "derivative"]
    alt = [r for r in rows if r["mhc_substitutions"]]

    LEGEND["notes"] = [
        f"{len(rows)} receptors covering {ns} of the 245 complexes and {3 * ns} of the 735 "
        "trajectories; the remaining 115 receptors are each present once. Two entries carry "
        "the same receptor when their alpha and beta variable domains are identical over "
        "the IMGT numbering, which is what ANARCI returns for the deposited sequences; "
        "constant domains, expression tags and cloning remnants are therefore not part of "
        "the comparison, and neither are the CDR3 strings on their own, which collide "
        "between receptors that differ elsewhere.",
        "The V and J calls, peptide sequences, CDR3 loops and per-row notes are in "
        "receptor_reuse.tsv, where peptides are counted as distinct sequences.",
        f"{len(var)} of these are engineered variants rather than natural receptors "
        f"({', '.join(sorted(var))}); they are listed as their own receptors, and the wild "
        "types they derive from are in the lineage table.",
        "The allele column reports the allele the construct was built from. Within an allele "
        f"the MHC still varies across {len(alt)} of these receptors, and mhc_substitutions in "
        "the source table gives the residues, recovered by aligning each deposited mature "
        "heavy chain against the unmodified form of its own allele; the comparison reproduces "
        "every substitution the depositions annotate and contradicts none. Most are "
        "engineered, two are not: 1mwa and 2ol3 carry the natural mouse alleles H-2Kbm3 and "
        "H-2Kbm8, which the resource records under H-2Kb. Substitutions above residue 182 sit "
        "outside the peptide-binding platform and are marked, and a set is not counted as an "
        "MHC comparison on their account.",
        "Two constructs, 2uwe and 6q3s, were labelled by sequence match with the rare "
        "natural alleles HLA-A*02:354 and HLA-A*02:624; 2uwe reproduces A*02:354 exactly, "
        "and 6q3s carries one substitution more than A*02:624. Both were built on "
        "HLA-A*02:01 and are reported as such, since either label would read as an allelic "
        "comparison neither study performed. The matches are kept in mhc_alleles_as_matched.",
    ]

    y = header(ax, 0, top, w, LEGEND["title"], LEGEND["subtitle"])
    y -= 3.6

    for x, t, ha in ((0.0, "Receptor", "left"), (X_MHC, "MHC", "left"),
                     (X_N, "entries", "right"), (X_NPEP, "peptides", "right"),
                     (X_PDB, "PDB identifiers", "left")):
        ax.text(x, y, t, fontproperties=F_HEAD, color=INK, ha=ha, va="baseline")
    y -= 1.6
    rule(ax, 0, w, y, lw=0.9)
    y -= LH

    for key, label in GROUPS:
        block = grp[key]
        if not block:
            continue
        y -= GAP_BAND
        ax.text(0, y, label, fontproperties=font(FS_BODY, bold=True), color=INK,
                ha="left", va="baseline")
        n = sum(int(r["n_structures"]) for r in block)
        ax.text(w, y, f"{len(block)} receptor{'s' * (len(block) > 1)}, {n} entries",
                fontproperties=F_BODY,
                color=SECONDARY, ha="right", va="baseline")
        y -= 1.3
        rule(ax, 0, w, y, lw=0.35, color=INK)
        y -= LH * 0.92

        for r in block:
            ax.text(0, y, r["receptor_name"], fontproperties=F_NAME, color=INK,
                    ha="left", va="baseline")
            ax.text(X_MHC, y, alleles(r["mhc_alleles"]), fontproperties=F_BODY,
                    color=SECONDARY, ha="left", va="baseline")
            ax.text(X_N, y, r["n_structures"], fontproperties=F_BODY, color=INK,
                    ha="right", va="baseline")
            ax.text(X_NPEP, y, r["n_peptides"], fontproperties=F_BODY, color=INK,
                    ha="right", va="baseline")
            ax.text(X_PDB, y, r["pdb_ids"], fontproperties=F_ID, color=INK,
                    ha="left", va="baseline")
            y -= LH

    y -= 1.0
    rule(ax, 0, w, y, lw=0.9)
    return notes(ax, 0, w, y - 1.2, LEGEND["notes"], gap=1.6, rule_lw=0.0)


if __name__ == "__main__":
    render(draw, "Table_reuse")
    rows, _ = load()
    # each cell must clear the next column: the name clears X_MHC, the allele
    # clears the right-aligned entry count, and the identifiers clear the page.
    lim = {"receptor_name": (X_MHC - 1.5, F_NAME),
           "mhc": (X_N - text_w("10", F_BODY) - X_MHC - 1.0, F_BODY),
           "pdb_ids": (W - X_PDB, F_ID)}
    bad = []
    for r in rows:
        for k, (lo, fp) in lim.items():
            t = alleles(r["mhc_alleles"]) if k == "mhc" else r[k]
            if text_w(t, fp) > lo:
                bad.append((r["clonotype_id"], k, round(text_w(t, fp), 1), lo))
    print("   column overflow:", bad or "none")
    print("   widest PDB list:", max((round(text_w(r["pdb_ids"], F_ID), 1),
                                      r["receptor_name"]) for r in rows),
          f"(available {W - X_PDB:.1f} mm)")
    heads = [(0.0, "Receptor", "left"), (X_MHC, "MHC", "left"),
             (X_N, "entries", "right"), (X_NPEP, "peptides", "right"),
             (X_PDB, "PDB identifiers", "left")]
    span = [(x if h == "left" else x - text_w(t, F_HEAD),
             x + text_w(t, F_HEAD) if h == "left" else x) for x, t, h in heads]
    print("   heading gaps:", [round(b[0] - a[1], 1)
                               for a, b in zip(span, span[1:])])
