#!/usr/bin/env python
"""Peptide bulge height against peptide length, all 245 complexes -- on its own.

Bulge is the time-averaged maximum height of a peptide Calpha above the plane
of the alpha1/alpha2 groove-helix rim, so it is the one scalar that says how
far out of the groove a peptide sits.  Length is the obvious thing to plot it
against: a groove of fixed length has to accommodate a peptide of variable
length somewhere, and the height records where.

Every complex is in the plot.  Boxes carry the quartiles and the whiskers the
full range, not a percentile, because these groups are small enough that the
range is the honest summary; the two lengths represented by a single complex
are drawn as that one point and named, rather than being dropped for having no
distribution -- which is what an n >= 5 rule would have done to them.

No complex carries a 5-, 6- or 7-residue peptide, so the axis breaks between 4
and 8 exactly as it does in the peptide-length figure, and for the same reason.
"""
import sys
from collections import defaultdict
from pathlib import Path

# figstyle / dyndata live one level up, shared by every figure folder
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from figstyle import *
from dyndata import N_TOTAL, bulge, quant

LEGEND = {}   # prose harvested by build_legends.py; see figstyle.LEAN

PH = 40.0                       # plot height
AXW = 11.0                      # left tick gutter
BOXW_MAX = 9.0                  # cap on box width so wide panels stay elegant
GAP_U = 0.50                    # an elided run gets half the width of a real slot
BREAK_PAD = 1.5                 # clear space each side of the ellipsis, mm
VMAX = 15.0
TICKS = [0, 3, 6, 9, 12, 15]
DOT_R = 0.85                    # radius of a single-complex marker, mm
# THE COUNT STRIP.  How many complexes carry each length used to be a row of
# digits under a stub head, which read as a table glued to the plot rather than
# as part of it -- and it hid the one fact about this axis worth seeing, that
# nearly two thirds of the library is a 9-mer.  It is a second plot on the same
# slots instead, sharing the length ticks, which now sit under both.
#
# PH2 is not a free choice.  This panel and the peptide-position panel beside it
# hang from one top edge and are read as one row, so their x titles have to land
# on one line -- two axis titles a centimetre apart under two plots of the same
# width is the first thing the eye finds and it reads as a slip.  Both panels put
# the title the same distance under their own lower baseline, so aligning them is
# just making the stacks equal:
#
#   here     6.0 + PH(40.0) + CGAP(6.0) + PH2
#   position 5.5 + PH(39.0) + ROWGAP(7.0) + PH(39.0)  = 90.5
#
# which fixes PH2 at 38.5 and spends the slack this panel had in its row on the
# strip rather than on white.  Changing PH, CGAP, or either constant in
# fig_peptide_position.py breaks the alignment; the arithmetic above is what to
# redo, not a number to nudge until it looks right.
PH2 = 38.5                      # plot height of the count strip -- see above
# Enough that the tallest bar's count clears the bulge baseline above it: the
# label sits 1.0 mm over a bar that fills PH2, and is 2.8 mm tall itself.
CGAP = 6.0


def draw(ax, top, w):
    by = defaultdict(list)
    ids = defaultdict(list)
    for r in bulge():
        by[r["peptide_length"]].append(r["bulge"])
        ids[r["peptide_length"]].append(r["pdb_id"])
    lo, hi = min(by), max(by)
    missing = [L for L in range(lo, hi + 1) if L not in by]
    solo = {L: (ids[L][0], by[L][0]) for L in by if len(by[L]) == 1}
    med = {L: quant(by[L], 0.5) for L in by}
    iqr = [(quant(by[L], .25), quant(by[L], .75)) for L in (8, 9, 10, 11)]
    assert all(a[1] < b[0] for a, b in zip(iqr, iqr[1:])), iqr   # the note claims this
    lows = sorted(zip(by[13], ids[13]))

    LEGEND.update(
        title="Peptide bulge by length",
        subtitle=f"All {N_TOTAL} complexes, grouped by the length of the simulated peptide. "
                 f"Bulge is the time-averaged maximum height of a peptide Cα above the plane "
                 f"of the α1/α2 groove-helix rim, averaged over the three replicas. Box, "
                 f"interquartile range; rule, median; whiskers, full range. Lengths carrying "
                 f"one complex are drawn as that complex and named. The lower strip counts the "
                 f"complexes at each length on the same slots. The axis is broken between "
                 f"{lo} and {missing[-1] + 1}, where no complex sits.",
        notes=[f"The median rises from {med[8]:.1f} Å at 8 residues through {med[9]:.1f} Å at "
               f"9 and {med[10]:.1f} Å at 10 to {med[11]:.1f} Å at 11, and the interquartile "
               f"ranges of those four groups do not overlap at all, so length accounts for most "
               f"of the variation in this descriptor across the library.",
               f"The 13-residue group is the exception, spanning "
               f"{min(by[13]):.1f}–{max(by[13]):.1f} Å over its {len(by[13])} complexes: "
               f"{lows[0][1].upper()} ({lows[0][0]:.1f} Å) and {lows[1][1].upper()} ({lows[1][0]:.1f} Å) keep a "
               f"long peptide close to the groove rather than arching it, while the other "
               f"{len(by[13]) - 2} fall in {lows[2][0]:.1f}–{lows[-1][0]:.1f} Å.",
               f"{solo[4][0].upper()} is the {lo}-residue lipopeptide and {solo[12][0].upper()} the one "
               f"12-residue complex; neither supports a distribution, and both are shown so "
               f"that the figure accounts for every complex in the library."],
    )
    y = header(ax, 0, top, w, LEGEND["title"], LEGEND["subtitle"])

    base = y - 6.0 - PH
    sc = PH / VMAX
    px0, pw = AXW, w - AXW

    # slot axis: every run of unoccupied lengths collapses to one narrow elision
    slots = []
    for L in range(lo, hi + 1):
        if L in by:
            slots.append(L)
        elif slots[-1] is not None:
            slots.append(None)
    sw = pw / sum(1.0 if s is not None else GAP_U for s in slots)
    bw = min(sw * 0.62, BOXW_MAX)

    x, cxs, breaks = px0, {}, []
    for s in slots:
        u = sw * (1.0 if s is not None else GAP_U)
        if s is None:
            breaks.append(x + u / 2)
        else:
            cxs[s] = x + u / 2
        x += u

    for t in TICKS:
        yy = base + t * sc
        if t:                               # zero is the baseline itself
            ax.plot([px0, px0 + pw], [yy, yy], lw=0.25, color=HAIR, zorder=1,
                    solid_capstyle="butt")
        ax.text(px0 - 1.6, yy, str(t), fontproperties=F_SMALL, color=SECONDARY,
                ha="right", va="center", zorder=5)
    axis_y(ax, base + PH / 2, "Peptide bulge (Å)")

    # baseline, interrupted at each break
    half = text_w("…", F_BODY) / 2 + BREAK_PAD
    edges = [px0] + [e for b in breaks for e in (b - half, b + half)] + [px0 + pw]
    for x0, x1 in zip(edges[::2], edges[1::2]):
        rule(ax, x0, x1, base, lw=0.7)
    for b in breaks:
        ax.text(b, base, "…", fontproperties=F_BODY, color=INK,
                ha="center", va="center", zorder=6)

    for L, cx in cxs.items():
        v = by[L]
        if L in solo:
            pid, val = solo[L]
            ax.add_patch(plt.Circle((cx, base + val * sc), DOT_R, facecolor=SECOND,
                                    edgecolor=SECOND_D, lw=0.5, zorder=5))
            ax.text(cx, base + val * sc + DOT_R + 1.0, pid.upper(),
                    fontproperties=font(FS_BODY, bold=True), color=SECOND_D,
                    ha="center", va="baseline", zorder=6)
            continue
        box(ax, cx, bw, *[quant(v, q) * sc + base for q in (0, .25, .5, .75, 1)],
            tint(PRIMARY, 0.55), PRIMARY_D, med_c=PRIMARY_D)

    # count strip -----------------------------------------------------
    # Same slots, same bar width, same break: the strip is the box plot's own
    # x axis drawn a second time as a quantity.  A bar takes the colour of the
    # mark above it, so the two lengths carrying one complex read as one thing
    # in both rows rather than as a box that lost its box.
    base2 = base - CGAP - PH2
    nmax = max(len(v) for v in by.values())
    sc2 = PH2 / nmax
    for L, cx in cxs.items():
        n = len(by[L])
        h = max(n * sc2, 0.4)
        # Solid, not the boxes' tint: a count is a central quantity, and the
        # tinted version of it came out lighter than the box outlines above it,
        # so the strip read as a shadow of the plot rather than as a plot.  This
        # is also the fill every other count bar in the paper carries.
        bar(ax, cx - bw / 2, base2, bw, h,
            SECOND if L in solo else PRIMARY, ends="top", r=0.45)
        # Every bar carries its count.  With 154 against 1 on one linear axis
        # the small groups are slivers by construction -- which is the point of
        # drawing it -- so the number is what makes them readable, and it rides
        # its own bar rather than sitting in a row with the others.
        ax.text(cx, base2 + h + 1.0, str(n), fontproperties=F_BODY, color=INK,
                ha="center", va="baseline", zorder=5)
    axis_y(ax, base2 + PH2 / 2, "Complexes")
    for x0, x1 in zip(edges[::2], edges[1::2]):
        rule(ax, x0, x1, base2, lw=0.7)
    for b in breaks:
        ax.text(b, base2, "\u2026", fontproperties=F_BODY, color=INK,
                ha="center", va="center", zorder=6)

    # the length ticks belong to both plots, so they sit under the lower one
    yl = base2 - 2.2
    for L, cx in cxs.items():
        ax.text(cx, yl, str(L), fontproperties=F_BODY, color=INK,
                ha="center", va="top", zorder=5)
    yx = yl - line_h(F_BODY) * 1.05
    axis_x(ax, px0 + pw / 2, yx - line_h(F_BODY) * 0.95, "Peptide length")

    return notes(ax, 0, w, yx - line_h(F_BODY) * 1.45, LEGEND["notes"])


if __name__ == "__main__":
    render(draw, "fig_bulge_by_length")
