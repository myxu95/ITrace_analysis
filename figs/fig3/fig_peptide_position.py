#!/usr/bin/env python
"""Peptide position profile across the 154 nine-mer complexes -- drawn on its own.

Two rows over one shared P1-P9 axis: how often the TCR is in contact with each
peptide residue, and how much of that residue is solvent accessible.  They are
the two fields that define an anchor, so putting them on the same ruler lets
the reader see the definition instead of being told the conclusion.

Boxes, not bars.  Both quantities are strongly bimodal across the library -- a
residue is either continuously touched or essentially never touched -- so a
mean lands in a gap where almost no complex sits.  The box gives the quartiles
and the whiskers the 10th-90th percentile, which is what actually describes a
library of 154 independent complexes.

The MHC side is not drawn.  Contact occupancy to the MHC saturates (median
0.95-1.00 at every position), so it separates nothing; solvent accessibility is
the axis that resolves burial, and it is the one plotted.

Only 9-mers.  They are 154 of the 245 complexes and the only length numerous
enough for a per-position reading; a shared position index across lengths would
compare P5 of an 8-mer with P5 of a 13-mer, which are not the same site.
"""
import sys
from pathlib import Path

# figstyle / dyndata live one level up, shared by every figure folder
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from figstyle import *
from dyndata import (ANCHOR_SASA_MAX, ANCHOR_TCR_MAX, FOCUS_LENGTH, N_TOTAL,
                     pep9, quant)

LEGEND = {}   # prose harvested by build_legends.py; see figstyle.LEAN

# A rotated row title is as long as the string, so the row has to be at least
# that tall or the title overruns it and collides with the next row's.
# "TCR contact occupancy" measures 37 mm at the current body size.
PH = 39.0                       # height of each of the two plot rows
AXW = 13.0                      # left tick gutter: rotated title + widest tick label
ROWGAP = 7.0                    # between the lower edge of one row and the top of the next
BOXW = 7.6                      # box width cap, mm; narrowed to fit a tight slot
WHISK = (0.10, 0.90)            # whisker percentiles; n = 154 per position

ROWS = [
    ("tcr",  "TCR contact occupancy", PRIMARY,  PRIMARY_D,
     1.0, [0.0, 0.25, 0.50, 0.75, 1.0], ["0", "0.25", "0.50", "0.75", "1.00"]),
    ("sasa", "Residue SASA (nm²)",    THIRD,    THIRD_D,
     0.66, [0.0, 0.2, 0.4, 0.6],       ["0", "0.2", "0.4", "0.6"]),
]


def draw(ax, top, w):
    d = pep9()
    pos = sorted(d)
    n = len(d[1]["tcr"])
    anchor = {p: 100.0 * sum(1 for t, s in zip(d[p]["tcr"], d[p]["sasa"])
                             if t < ANCHOR_TCR_MAX and s < ANCHOR_SASA_MAX) / n
              for p in pos}
    m = {p: quant(d[p]["tcr"], 0.5) for p in pos}
    hot = [p for p in pos if m[p] >= 0.80]
    cold = [p for p in pos if p not in hot]
    top_anchor = sorted(pos, key=lambda p: -anchor[p])[:2]

    LEGEND.update(
        title="Peptide position profile",
        subtitle=f"The {n} nine-residue complexes of the {N_TOTAL}, position by position. "
                 f"Box, interquartile range across complexes; rule, median; whiskers, "
                 f"{WHISK[0]:.0%}–{WHISK[1]:.0%} percentiles. Upper row, the maximum contact "
                 f"occupancy between that peptide residue and any TCR residue; lower row, its "
                 f"time-averaged solvent-accessible surface area. Both are per-complex means "
                 f"over the three replicas.",
        notes=[f"Median TCR occupancy is {min(m[p] for p in hot):.2f}–"
               f"{max(m[p] for p in hot):.2f} at P" + ", P".join(str(p) for p in hot)
               + f" and {min(m[p] for p in cold):.2f}–{max(m[p] for p in cold):.2f} at P"
               + ", P".join(str(p) for p in cold)
               + ", so the central residues are in contact in nearly every complex while the "
                 "flanks are in contact in some and not others; the interquartile ranges at P1 "
                 "and P3 span most of the 0–1 scale, which is the variation a mean would hide.",
               f"peptide_table flags a position anchor when its solvent accessibility is "
               f"below {ANCHOR_SASA_MAX:.2f} nm² and its TCR occupancy below "
               f"{ANCHOR_TCR_MAX:.2f}; both are fixed cut-offs on measured quantities, not "
               f"quantiles of this library. The flag is not drawn -- the two rows are the "
               f"quantities it is computed from. The two commonest anchors are "
               + " and ".join(f"P{p} ({anchor[p]:.0f}%)"
                              for p in sorted(top_anchor, key=lambda p: -anchor[p]))
               + ", so the canonical anchor positions are a majority of the library rather "
                 "than a property of all of it.",
               "Contact occupancy to the MHC is not plotted because it saturates: its median "
               f"is {min(quant(d[p]['hla'], 0.5) for p in pos):.2f}–"
               f"{max(quant(d[p]['hla'], 0.5) for p in pos):.2f} across the nine positions, so "
               "solvent accessibility, not MHC contact, is the quantity that resolves burial."],
    )
    y = header(ax, 0, top, w, LEGEND["title"], LEGEND["subtitle"])

    px0, pw = AXW, w - AXW
    sw = pw / len(pos)
    bw = min(BOXW, sw * 0.72)   # nine boxes have to clear each other at any width
    cx = {p: px0 + (i + 0.5) * sw for i, p in enumerate(pos)}

    bases = []
    ytop = y - 5.5
    for key, title, fc, ec, vmax, ticks, labels in ROWS:
        base = ytop - PH
        bases.append(base)
        sc = PH / vmax

        for t, lab in zip(ticks, labels):
            yy = base + t * sc
            if t:                           # zero is the baseline itself
                ax.plot([px0, px0 + pw], [yy, yy], lw=0.25, color=HAIR, zorder=1,
                        solid_capstyle="butt")
            ax.text(px0 - 1.6, yy, lab, fontproperties=F_SMALL, color=SECONDARY,
                    ha="right", va="center", zorder=5)
        rule(ax, px0, px0 + pw, base, lw=0.7)
        axis_y(ax, base + PH / 2, title)

        for p in pos:
            v = d[p][key]
            box(ax, cx[p], bw,
                *[quant(v, q) * sc + base
                  for q in (WHISK[0], 0.25, 0.50, 0.75, WHISK[1])],
                tint(fc, 0.55), ec, med_c=ec)
        ytop = base - ROWGAP

    # one set of position labels, under the lower row and shared by both
    yl = bases[1] - 2.2
    for p in pos:
        ax.text(cx[p], yl, f"P{p}", fontproperties=F_BODY, color=INK,
                ha="center", va="top", zorder=5)

    yx = yl - line_h(F_BODY) * 2.0
    axis_x(ax, px0 + pw / 2, yx, "Peptide position")
    return notes(ax, 0, w, yx - line_h(F_BODY) * 0.60, LEGEND["notes"])


if __name__ == "__main__":
    render(draw, "fig_peptide_position")
