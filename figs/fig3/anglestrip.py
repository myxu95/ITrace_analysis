"""The ranked-band drawing shared by the two angle-range figures.

Incident and crossing angle want exactly the same picture -- the whole library
as one band -- but they cannot share an x-axis, because ranking by one angle
scrambles the other.  So they are two figures, and the drawing lives here once
rather than being copied into both scripts.

The line is each complex's median angle over its pooled 3 x 1001 frames and the
shading is that complex's 5th to 95th percentile, so the height of the band at
any point along the axis is how much that one complex's angle fluctuates over
its trajectories, and the rise of the line across the axis is the range the
library covers.  Percentiles rather than the full range: over three thousand
frames the extremes describe one frame, and the question the figure answers is
where the angle normally sits.

Ranking by the median is what makes the line a line -- it is monotonic by
construction, and says nothing on its own.  The band edges are not, and their
unevenness against that smooth line is the actual reading: two complexes with
the same median need not be equally steady.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from figstyle import (F_BODY, F_SMALL, FS_BODY, HAIR, INK, PRIMARY, PRIMARY_D,
                      SECOND, SECOND_D, SECONDARY, LW_TREND, axis_x, axis_y,
                      header, line_h, rect, rule, tint)

PH = 46.0                   # plot height, mm
AXW = 12.0                  # left tick gutter
RANK_TICKS = (1, 50, 100, 150, 200)
LW_MED = LW_TREND           # the median line (shared with Fig 4a,b,d)
# How far the band is tinted toward white.  Raised from 0.70: the band is the
# fluctuation of one complex, drawn behind the median line, and at 0.70 it read
# as a second filled quantity competing with the line rather than as the spread
# around it.  Not raised further -- HAIR, the gridline tone, is #e1e7ef, and a
# band tinted past ~0.80 stops being separable from the grid it sits on.
BAND_F = 0.78
KEY_W, KEY_H, KEY_GAP = 6.4, 2.8, 1.5   # class key: swatch size and leading


def spans(flags):
    """Consecutive equal flags as ``(first, last, flag)`` index ranges."""
    out, i = [], 0
    while i < len(flags):
        j = i
        while j + 1 < len(flags) and flags[j + 1] == flags[i]:
            j += 1
        out.append((i, j, flags[i]))
        i = j + 1
    return out


def strip(ax, top, w, rows, key, vmax, ticks, ylab, legend, exception=(),
          classes=None):
    """Median line over a percentile band, ranked; returns the lowest ink.

    ``exception`` names the complexes drawn in the contrast colour -- the house
    rule is one hue plus one for the exception class, never a hue per complex.
    Each run of them gets its own band and its own line rather than a recoloured
    stretch of one continuous curve: the two classes are separate populations,
    and joining them would draw a slope where there is only a change of kind.
    ``classes`` names the two, ordinary first, and asks for the key that says so
    -- a second colour with nothing to explain it is just a stray mark.
    """
    y = header(ax, 0, top, w, legend["title"], legend["subtitle"])
    base = y - 6.0 - PH
    sc = PH / vmax
    px0, pw = AXW, w - AXW

    for t in ticks:
        yy = base + t * sc
        if t:
            ax.plot([px0, px0 + pw], [yy, yy], lw=0.25, color=HAIR, zorder=1,
                    solid_capstyle="butt")
        ax.text(px0 - 1.6, yy, str(t), fontproperties=F_SMALL, color=SECONDARY,
                ha="right", va="center", zorder=5)
    axis_y(ax, base + PH / 2, ylab)

    order = sorted(rows, key=lambda r: r[f"{key}_p50"])
    sw = pw / len(order)
    xs = [px0 + (i + 0.5) * sw for i in range(len(order))]
    lo, mid, hi = ([base + r[f"{key}_p{p}"] * sc for r in order]
                   for p in (5, 50, 95))

    flags = [r["pdb_id"] in exception for r in order]
    for i0, i1, odd in spans(flags):
        s = slice(i0, i1 + 1)
        # half a slot of overhang at each end, so that consecutive runs abut
        # exactly and the field is filled edge to edge: a complex owns its slot
        # on the rank axis, not just the point at the centre of it
        band = ([xs[i0] - sw / 2] + xs[s] + [xs[i1] + sw / 2],
                [lo[i0]] + lo[s] + [lo[i1]],
                [hi[i0]] + hi[s] + [hi[i1]])
        ax.fill_between(*band, linewidth=0, zorder=2,
                        facecolor=tint(SECOND if odd else PRIMARY, BAND_F))
        ax.plot(xs[s], mid[s], lw=LW_MED, zorder=4, solid_capstyle="round",
                solid_joinstyle="round", color=SECOND_D if odd else PRIMARY_D)

    if classes:
        # top left: the only corner of either angle axis that no band reaches
        n_odd = sum(flags)
        ky = base + PH - 1.6 - KEY_H
        for lab, n, hue, ln in ((classes[0], len(order) - n_odd, PRIMARY, PRIMARY_D),
                                (classes[1], n_odd, SECOND, SECOND_D)):
            rect(ax, px0 + 3.0, ky, KEY_W, KEY_H, tint(hue, BAND_F), z=5)
            ax.plot([px0 + 3.0, px0 + 3.0 + KEY_W], [ky + KEY_H / 2] * 2,
                    lw=LW_MED, color=ln, zorder=6, solid_capstyle="butt")
            ax.text(px0 + 3.0 + KEY_W + 1.6, ky + KEY_H / 2, f"{lab} ({n})",
                    fontproperties=F_BODY, color=INK, ha="left", va="center",
                    zorder=6)
            ky -= KEY_H + KEY_GAP

    rule(ax, px0, px0 + pw, base, lw=0.7)
    last = len(order)
    for t in RANK_TICKS + (last,):
        cx = px0 + (t - 0.5) * sw
        ax.plot([cx, cx], [base, base - 0.9], lw=0.5, color=INK, zorder=5,
                solid_capstyle="butt")
        # the final label would hang off the page if it were centred on its tick
        ax.text(cx, base - 1.7, str(t), fontproperties=F_SMALL, color=SECONDARY,
                ha="right" if t == last else "center", va="top", zorder=5)
    yl = base - 1.7 - line_h(F_SMALL) * 0.95
    axis_x(ax, px0 + pw / 2, yl - line_h(F_BODY) * 0.95, "Complex (ranked)")
    return yl - line_h(F_BODY) * 0.95
