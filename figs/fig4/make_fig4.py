#!/usr/bin/env python
"""Figure 4 -- technical validation of the deposited trajectories.

Four panels answering one question in two registers.  a, b and c are library
statements: every complex in the archive, one mark each, so a reader can see
the whole distribution rather than a sentence about it.  d is an entry
statement: one released complex opened up, so a reader can see what a single
download actually contains.  The two registers are the two things a data
descriptor has to establish -- that the collection holds together, and that an
individual record is worth downloading.

Each library panel answers its quantity's own version of "is this spread
small".  a can answer it by inspection: a reader knows what losing a tenth of
the native contacts means.  b and c cannot, so each is drawn against something.

b was a second ranked strip of the incident angle and said nothing Figure 3d
had not already said, so it now uses 3d's band instead of repeating it.  Every
complex is one point: the angular interval its frames occupy against the
disagreement between its three replica means, with equality drawn in.  The band
is a total and contains the disagreement plotted against it, which makes this a
decomposition and not an independent test; the manuscript carries a second
yardstick with the between-run part removed, and the docstring of
repdata.angle_variation says why both are quoted.  The sectors above the
scatter are those two medians drawn as angles at true scale, because nobody
reads 11 degrees of incident-angle travel off a number.

c is drawn against arithmetic rather than against another measurement.  Ten
modes drawn at
random out of 3N would overlap by sqrt(10/3N) -- 0.06 to 0.07 across this
library -- and the drawn band sits an order of magnitude above it.  That line is
why c is here and not in the supplement.  Most of c is empty on purpose: the
axis has to reach down to the chance band for the comparison to be visible, and
the gap between the band and the data is the statement.  The same is true of
b's upper left, which is the region where replicas would disagree by more than
the trajectories move.

d shows what the entry is for.  a, b and c say the trajectory is sound; d says
it resolves structure, which is the reason to download it.

The comparison of the crystallographic B-factor against our RMSF, which stood
here as a fifth and sixth panel, is a single-entry illustration that cannot
carry a library-wide claim; it is Figure 5 instead, where it can be shown next
to the same statistic computed over the whole release.

Colour carries panel identity here.  Blue is Q, ochre is the incident angle,
teal is RMSIP, and d is drawn in a fourth ramp, so a reader flicking back to a
number knows which panel it came from without reading the letter.  Inside b the
two tones of ochre carry the weight rule the strips carry between a segment and
its tick: the interval the frames cover is the light context, the disagreement
between runs is the dark quantity under test.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from figstyle import *
from landdata import exemplar as land_exemplar
from matplotlib.colors import LinearSegmentedColormap
from math import cos, floor, log10, radians, sin

from repdata import angle_variation, icc, quant, replicas

LEGEND = {}
DROP_FLAGGED = False
N_SET = 245

# panel letters: 8 pt bold lower case, which is what the journal sets
F_PANEL = font(FS_TAG, bold=True)
# 10 modes out of 3N atoms overlap by sqrt(10/3N) by chance alone
RAND_LO, RAND_HI = 0.063, 0.074


def sf2(v):
    """Two significant figures, which is the manuscript's ceiling everywhere.

    Round first, then take the exponent from the rounded value, so 9.96 prints
    as 10 and not as 10.0.
    """
    if v == 0:
        return "0"
    e = floor(log10(abs(v)))
    d = max(0, 1 - e)
    r = round(v, d)
    e2 = floor(log10(abs(r))) if r else e
    if e2 != e:
        d = max(0, 1 - e2)
        r = round(v, d)
    return f"{r:.{d}f}"

ROW1 = 56.0             # a, b
ROW2 = 44.0             # c, d
VGAP = 5.2
COLGAP = 7.0

AXW = 11.0              # tick gutter of a strip panel
MGAP = 3.2              # field to marginal
MW = 11.0               # marginal width
# Right gutter: the marginal's tick numbers, then its rotated title.  Its parts
# are RTICK + the widest number + RTITLE + half the rotated title's own height,
# and it was set at 10.6 against an arithmetic that took text_w at face value
# and left the title's box 2.4 mm from a number measured ~4% short -- under a
# millimetre of real air, which is what a reader sees as the two touching.
RAXW = 12.0             # right gutter: the marginal's tick numbers + its title
RTICK = 2.2             # marginal edge to its tick numbers
RTITLE = 3.6            # those numbers to the centre of the rotated title
SAFE = 1.10             # text_w is a metric estimate and runs ~4% short
NBINS = 30
PGUT = 7.2              # left gutter of a surface panel, for its rotated title

# WEIGHT.  Two hundred and forty-five vertical segments have more ink in them
# than anything else on the page, so they are drawn light and the marks a
# reader is meant to read first -- the complex mean, the median rule -- are
# drawn dark.  Light here is a tint toward white, not an alpha: the figure is
# exported without an alpha channel, and a line at alpha a over white is the
# same colour as the same line tinted 1 - a, so SEG_T = 0.48 is what the eye
# receives from alpha 0.52.  Tinting keeps that appearance in a file that has
# no transparency to lose in conversion.
SEG_T = 0.48            # replica-range segments   (~ alpha 0.52)
# The marginal's fill.  Matched to BAND_F in fig3/anglestrip.py, because the
# two are the same kind of mark -- a filled area standing for a spread, with a
# dark line read off it -- and a reader crossing between the figures should not
# meet the same idea at two weights.  The segments keep their own value: they
# are hairlines, not a fill, and the same tint on a 0.4 mm line puts it under
# the paper.
BAR_T = 0.78            # marginal silhouette      (= fig3 BAND_F)
LW_SEG = 0.40           # segment
LW_MEAN = LW_TREND      # complex mean: the darkest, heaviest mark of the three
LW_MED = 0.45           # median rule: same dark tone, deliberately thinner

FELN = 4                # basins drawn
BASIN = ["#5d4585", "#9c4a52", "#96762e", "#3f7a52"]
ZMAX = 6.0
FILL = [i * 0.2 for i in range(int(ZMAX / 0.2) + 1)]
ISO = [1.0, 2.0, 3.0, 4.0, 5.0]
FELMAP = LinearSegmentedColormap.from_list(
    "fel", [tint(PRIMARY_D, f) for f in (0.02, 0.30, 0.62, 0.86, 0.97)])


# ------------------------------------------------------------------ helpers
def letter(ax, s, x, y):
    ax.text(x, y, s, fontproperties=F_PANEL, color=INK, ha="left", va="bottom",
            zorder=99)


def strip(ax, x0, base, w, h, stem, vmax, ticks, fmt, hue, hue_d, ylab,
          medlab, rangelab=None, band=None):
    """One quantity, every complex, ranked by its own mean.

    A segment is the lowest to the highest of a complex's three replica means
    and a tick is the complex mean, so segment height IS the disagreement
    between runs and the height of the field is the range the library covers.
    The two are drawn on one axis because that is the comparison the panel
    exists to make: a reader should not have to hold a number from one plot in
    their head while looking at another.

    The marginal bins those same segment heights on that same axis, so "how far
    replicas disagree" and "how far the library spreads" are read off one ruler.
    One ruler, but not one quantity: down the field the axis reads as the value
    itself and down the marginal it reads as the disagreement between replicas,
    which is why the marginal carries its own rotated title on the right rather
    than borrowing the field's.
    """
    rows = replicas(stem, drop_flagged=DROP_FLAGGED)
    # vmax is hand-set per panel and the data under it is recomputed, so check
    # rather than trust: a segment past the top is not clipped, it is drawn
    # over whatever panel sits above this one.
    hi = max(v for r in rows for v in r["values"])
    assert hi <= vmax, (f"{stem}: {hi:.3f} exceeds the field ceiling {vmax} -- "
                        f"raise vmax and its ticks")
    order = sorted(rows, key=lambda r: r["mean"])
    fw = w - AXW - MGAP - MW - RAXW
    mx0 = x0 + AXW + MGAP + fw
    sc = h / vmax
    Y = lambda v: base + v * sc

    if band:                      # the chance-level reference, drawn under all.
        # Neutral, not the panel hue: it is not a measurement of this library,
        # it is the number this library has to beat.
        # Across the field only.  It used to run on through the marginal, but
        # the marginal's axis is not the value -- it is the disagreement
        # between replicas -- so a chance level for the value has no meaning
        # over it, and in d, where every spread is small, the grey bar landed
        # on top of the silhouette and the two read as one mark.
        rect(ax, x0 + AXW, Y(band[0]), fw, (band[1] - band[0]) * sc,
             tint(SECONDARY, 0.58), z=1)
    for t in ticks:
        if t:
            ax.plot([x0 + AXW, mx0 + MW], [Y(t)] * 2, lw=0.25, color=HAIR,
                    zorder=2, solid_capstyle="butt")
        ax.text(x0 + AXW - 1.4, Y(t), fmt(t), fontproperties=F_SMALL,
                color=SECONDARY, ha="right", va="center", zorder=6)
    axis_y(ax, base + h / 2, ylab, x=x0 + 2.0)

    sw = fw / len(order)
    mxs, mys = [], []
    for i, r in enumerate(order):
        cx = x0 + AXW + (i + 0.5) * sw
        lo, hi = min(r["values"]), max(r["values"])
        c = tint(hue, SEG_T)
        ax.plot([cx, cx], [Y(lo), Y(hi)], lw=LW_SEG, color=c, zorder=3,
                solid_capstyle="butt")
        mxs.append(cx)
        mys.append(Y(r["mean"]))
    # The complex means are monotonic by construction -- the rank axis IS their
    # order -- so they join into one continuous curve instead of 245 detached
    # ticks.  Same mark, same weight and same reading as the ranked median of
    # Figure 3c,d, which is the other place this library is drawn end to end.
    ax.plot(mxs, mys, lw=LW_MEAN, color=hue_d, zorder=4,
            solid_capstyle="round", solid_joinstyle="round")
    rule(ax, x0 + AXW, x0 + AXW + fw, base, lw=0.55, color=SECONDARY)
    for t, lab in ((0.5, "1"), (len(order) - 0.5, str(len(order)))):
        cx = x0 + AXW + t * sw
        ax.text(cx, base - 1.4, lab, fontproperties=F_SMALL, color=SECONDARY,
                ha="left" if lab == "1" else "right", va="top", zorder=6)

    spreads = [r["spread"] for r in rows]
    hist = [0] * NBINS
    for s in spreads:
        hist[min(int(s / vmax * NBINS), NBINS - 1)] += 1
    hmax = max(hist)
    bh = h / NBINS
    # One silhouette rather than thirty rectangles.  The bars carried a reading
    # they do not have: a rectangle with its own edges says each bin is a thing
    # to be looked at, when the only thing here worth looking at is the shape
    # of the whole distribution.  Drawn as a frequency polygon through the bin
    # centres and closed onto the axis -- the same histogram, no smoothing
    # model laid over it, and the same mark Fig 3c,d uses for a spread.
    ys = [base + (i + 0.5) * bh for i in range(NBINS)]
    xs = [mx0 + n / hmax * MW for n in hist]
    ax.fill_betweenx([base] + ys + [base + h], mx0, [mx0] + xs + [mx0],
                     linewidth=0, zorder=3, facecolor=tint(hue, BAR_T))
    med = quant(spreads, 0.5)
    ax.plot([mx0, mx0 + MW], [Y(med)] * 2, lw=LW_MED, color=hue_d, zorder=6,
            dashes=(1.6, 1.2), solid_capstyle="butt")
    ytxt = base + h - line_h(F_SMALL) * 0.25
    ax.text(mx0 + MW, ytxt, medlab(med), fontproperties=F_SMALL, color=hue_d,
            ha="right", va="top", zorder=7)
    ax.plot([mx0, mx0 + 2.2], [ytxt - line_h(F_SMALL) * 0.38] * 2, lw=LW_MED,
            color=hue_d, zorder=7, dashes=(1.6, 1.2), solid_capstyle="butt")
    rule(ax, mx0, mx0 + MW, base, lw=0.55, color=SECONDARY)
    # The marginal shares the field's gridlines but not its meaning -- down the
    # field the axis is the value, down here it is the disagreement between
    # replicas -- so it carries its own numbers on its own side.  A rotated
    # title over an unnumbered axis names a quantity without letting anyone
    # read one off it.
    tw = max(text_w(fmt(t), F_SMALL) for t in ticks) * SAFE
    for t in ticks:
        ax.text(mx0 + MW + RTICK, Y(t), fmt(t), fontproperties=F_SMALL,
                color=SECONDARY, ha="left", va="center", zorder=6)
    if rangelab:
        axis_y(ax, base + h / 2, rangelab, x=mx0 + MW + RTICK + tw + RTITLE)
    return dict(rows=rows, spreads=spreads, med=med, fw=fw, mx0=mx0,
                x_mid=x0 + AXW + fw / 2, m_mid=mx0 + MW / 2)


def surface(ax, x0, base, w, h, d, title, ylab=True):
    """One dihedral-PCA free-energy surface with its deepest basins numbered."""
    fel = d["fel"]
    X, Y_, Z = fel["x"], fel["y"], fel["z"]
    sx, sy = w / (X[-1] - X[0]), h / (Y_[-1] - Y_[0])
    xs = [x0 + (v - X[0]) * sx for v in X]
    ys = [base + (v - Y_[0]) * sy for v in Y_]
    zc = [[min(max(v, 0.0), ZMAX) for v in row] for row in Z]
    ax.contourf(xs, ys, zc, levels=FILL, cmap=FELMAP, zorder=2)
    ax.contour(xs, ys, zc, levels=ISO, colors="white", linewidths=0.22, zorder=3)
    for i, b in enumerate(d["basins"][:FELN], 1):
        bx = x0 + (b["pc1"] - X[0]) * sx
        by = base + (b["pc2"] - Y_[0]) * sy
        ax.add_patch(plt.Circle((bx, by), 1.25, facecolor=BASIN[i - 1],
                                edgecolor="white", lw=0.6, zorder=7))
        ax.text(bx, by, str(i), fontproperties=font(FS_TINY, bold=True),
                color="white", ha="center", va="center", zorder=8)
    ax.add_patch(plt.Rectangle((x0, base), w, h, facecolor="none",
                               edgecolor=SECONDARY, lw=0.45, zorder=9))
    ax.text(x0 + w / 2, base + h + 1.1, title, fontproperties=F_HEAD, color=INK,
            ha="center", va="baseline", zorder=6)
    ax.text(x0 + w / 2, base - 1.3, f"PC1 ({d['pc1_var_frac']:.0%})",
            fontproperties=F_SMALL, color=SECONDARY, ha="center", va="top",
            zorder=6)
    if ylab:
        axis_y(ax, base + h / 2, f"PC2 ({d['pc2_var_frac']:.0%})", x=x0 - 4.6)


def cbar(ax, x0, y, w, h, cmap, lo, hi, ticks, unit, fmt="{:g}", clip=False):
    n = 90
    for i in range(n):
        rect(ax, x0 + i * w / n, y, w / n + 0.04, h, cmap(i / (n - 1.0)), z=3)
    ax.add_patch(plt.Rectangle((x0, y), w, h, facecolor="none",
                               edgecolor=SECONDARY, lw=0.4, zorder=6))
    for v in ticks:
        cx = x0 + (v - lo) / (hi - lo) * w
        ax.text(cx, y - 1.0, ("\u2265 " if clip and v == ticks[-1] else "") + fmt.format(v),
                fontproperties=F_SMALL,
                color=SECONDARY, va="top", zorder=6,
                ha="left" if v == ticks[0] else "right" if v == ticks[-1]
                else "center")
    ax.text(x0 + w + 1.6, y + h / 2, unit, fontproperties=F_SMALL, color=SECONDARY,
            ha="left", va="center", zorder=6)


# ------------------------------------------------------------------ panels
# Each panel draws itself from its own origin and hands back what it measured
# together with the y of its lowest ink, which is the contract render() wants.
# The assembled plate and the four separate panel files therefore come out of
# one drawing of each panel rather than two, and cannot drift apart.
def _foot(ax, P, base):
    """The two x titles under a strip panel; returns the y of its lowest ink."""
    # A full line, not most of one.  The rank axis carries "1" and "245" set
    # va="top" from the baseline, so a title dropped 0.95 of a line lands
    # inside them -- which it did not at 7 pt, when the title was narrow
    # enough to sit between the two.
    # "Complex, ranked", not "Complex, ranked (245)".  The two titles are
    # centred on their own fields, whose midpoints are 24.3 mm apart in a and b
    # and 22.8 mm in the narrower c.  With the count in it the field title
    # needs 20.4 mm of half-width against Count's 5.5, so the pair wanted 25.9
    # and ran into each other in all three panels -- worst in c, where they
    # overlapped by three millimetres.  Without it the pair wants 20.8 and
    # clears in c by 2.0 mm.  Nothing is lost: the rank axis already carries 1
    # and 245 as its end ticks.
    y = base - 1.4 - line_h(F_SMALL) * 1.55
    axis_x(ax, P["x_mid"], y, "Complex, ranked")
    axis_x(ax, P["m_mid"], y, "Count")
    return y - line_h(F_HEAD) * 0.9


def panel_a(ax, x0, top, w):
    """The fraction of native contacts, Q, for every complex."""
    letter(ax, "a", x0, top - 2.6)
    h = ROW1 - 9.0
    b = top - 5.4 - h
    P = strip(ax, x0, b, w, h, "q", 1.0, (0, .2, .4, .6, .8, 1.0),
              lambda t: f"{t:.1f}", PRIMARY, PRIMARY_D,
              "Fraction of native contacts, *Q*",
              sf2, rangelab="Replica range (*Q*)")
    return P, _foot(ax, P, b)


# Panel b's two angular quantities, both in degrees.  BAND is the whole
# interval a complex's frames occupy, REP is how far the three replica means
# disagree about where that interval sits.  The weight rule of this figure
# applies to the pair as it does to a segment and its tick: the context is
# light, the quantity being validated is dark.
BAND_C = tint(SECOND, 0.66)
BAND_E = tint(SECOND_D, 0.42)
REP_C = SECOND_D
WEDH = 11.0             # the schematic block at the top of panel b
WGAP = 3.4              # schematic to scatter
WDROP = 1.5             # air under the schematic's baseline
DOT_R = 0.42            # one complex in the scatter
XLO, XHI = 5.0, 70.0            # degrees covered   (log)
YLO, YHI = 0.25, 35.0           # replica range (°) (log)
XT = (5, 10, 20, 50)
YT = (0.3, 1, 3, 10, 30)


def _degfmt(t):
    return f"{t:g}"


def wedge(ax, vx, vy, r, deg, fc, ec, z):
    """A true-scale angular sector opening to the right from (vx, vy)."""
    n = 24
    pts = [(vx, vy)] + [(vx + r * cos(radians(deg * i / n)),
                         vy + r * sin(radians(deg * i / n))) for i in range(n + 1)]
    ax.add_patch(plt.Polygon(pts, closed=True, facecolor=fc, edgecolor=ec,
                             lw=0.4, zorder=z, joinstyle="round"))


def panel_b(ax, x0, top, w):
    """How much the incident angle changes, and how much of that is run to run.

    Figure 3d already draws where each complex's incident angle sits and how
    wide a band its frames occupy.  This panel does not redraw that band -- it
    uses it as the denominator.  The scatter puts, for every complex, the
    angular interval its frames cover against the disagreement between its
    three replica means, on one pair of degree axes with the line of equality
    drawn in.  Below that line the run-to-run offset is smaller than the motion
    the trajectories already show, which is the reproducibility statement the
    panel exists to make; the dashed line is where the library's median sits.

    The band is a total and therefore contains the disagreement plotted against
    it, so the comparison is a decomposition rather than an independent test.
    That is stated rather than hidden, and the manuscript quotes a second
    yardstick with the between-run part removed -- the replica range against the
    SD inside a single trajectory -- next to this one.

    The schematic above the scatter is the same two medians drawn as angles at
    true scale, because a reader has no intuition for what 11 degrees of
    incident-angle travel looks like and a number on an axis does not give them
    one.  It is drawn at the largest radius the column allows for exactly that
    reason: at a smaller radius both sectors collapse to slivers and the
    schematic stops being a measurement.
    """
    letter(ax, "b", x0, top - 2.6)
    h = ROW1 - 9.0
    b = top - 5.4 - h
    D = angle_variation(drop_flagged=DROP_FLAGGED)
    m_band = quant([d["band"] for d in D], 0.5)
    m_rep = quant([d["spread"] for d in D], 0.5)
    m_ratio = quant([d["ratio"] for d in D], 0.5)

    gut = 9.6
    fw = w - gut - 1.0
    sh = h - WEDH - WGAP
    fx = x0 + gut

    # ---- the schematic: the two medians as true-scale sectors
    wy = b + sh + WGAP + WDROP
    lab_hi = f"{sf2(m_band)}\u00b0 covered"
    lab_lo = f"{sf2(m_rep)}\u00b0 run to run"
    lw_ = max(text_w(lab_hi, F_SMALL), text_w(lab_lo, F_SMALL)) * SAFE
    r = fw - lw_ - 2.2
    y_hi = wy + r * sin(radians(m_band))
    y_lo = wy + r * sin(radians(m_rep))
    # Each label sits at its own sector's edge.  If the two medians are close
    # enough that one line of type will not fit between them, the lower label
    # drops just far enough to clear -- a couple of tenths of a millimetre reads
    # as still attached to its ray, a couple of millimetres does not.
    y_lab = min(y_lo, y_hi - line_h(F_SMALL))
    assert y_lo - y_lab <= 1.5, (
        f"schematic radius {r:.1f} mm would detach the lower label by "
        f"{y_lo - y_lab:.1f} mm -- shorten the labels or widen the column")
    wedge(ax, fx, wy, r, m_band, BAND_C, BAND_E, 3)
    wedge(ax, fx, wy, r, m_rep, tint(REP_C, 0.10), REP_C, 4)
    ax.plot([fx, fx + r + 1.0], [wy] * 2, lw=0.45, color=SECONDARY, zorder=5,
            solid_capstyle="butt")
    ax.text(fx + r + 1.6, y_hi, lab_hi, fontproperties=F_SMALL, color=BAND_E,
            ha="left", va="center", zorder=6)
    ax.text(fx + r + 1.6, y_lab, lab_lo, fontproperties=F_SMALL, color=REP_C,
            ha="left", va="center", zorder=6)

    # ---- the scatter
    lx = log10(XLO)
    sx = fw / (log10(XHI) - lx)
    ly = log10(YLO)
    sy = sh / (log10(YHI) - ly)
    X = lambda v: fx + (log10(v) - lx) * sx
    Y = lambda v: b + (log10(v) - ly) * sy
    for t in YT:
        ax.plot([fx, fx + fw], [Y(t)] * 2, lw=0.25, color=HAIR, zorder=2,
                solid_capstyle="butt")
        ax.text(fx - 1.4, Y(t), _degfmt(t), fontproperties=F_SMALL,
                color=SECONDARY, ha="right", va="center", zorder=6)
    axis_y(ax, b + sh / 2, "Replica range (\u00b0)", x=x0 + 1.8)

    # equality first, under the data: it is not a measurement of this library,
    # it is the line the library has to stay under, so it is drawn neutral for
    # the same reason c's chance band is.
    xe = min(XHI, YHI)
    ax.plot([X(XLO), X(xe)], [Y(XLO), Y(xe)], lw=0.55,
            color=tint(SECONDARY, 0.40), zorder=3, solid_capstyle="butt")
    # Labelled at the low end, where the field is empty: at the high end the
    # equality line runs through the densest part of the cloud.
    ax.text(X(XLO * 1.10), Y(XLO * 1.10) + 0.5, "equal", fontproperties=F_SMALL,
            color=SECONDARY, ha="left", va="bottom", zorder=7)
    ax.plot([X(XLO), X(XHI)], [Y(XLO * m_ratio), Y(XHI * m_ratio)], lw=LW_MED,
            color=REP_C, zorder=6, dashes=(1.6, 1.2), solid_capstyle="butt")
    # The label hangs under the dashed line at the right edge, but the line
    # descends leftwards across the width of the text, so the clearance has to
    # be measured at the text's LEFT end or the line cuts through it.  The
    # panel is rendered at two widths, so this cannot be a fixed offset.
    lab = f"median {sf2(m_ratio)}\u00d7"
    lw_r = text_w(lab, F_SMALL) * SAFE
    slope = (Y(XHI * m_ratio) - Y(XLO * m_ratio)) / (X(XHI) - X(XLO))
    ax.text(X(XHI), Y(XHI * m_ratio) - slope * lw_r - 0.8, lab,
            fontproperties=F_SMALL, color=REP_C, ha="right", va="top", zorder=7)

    for d in D:
        assert XLO <= d["band"] <= XHI and YLO <= d["spread"] <= YHI, (
            f"{d['pdb_id']}: {d['band']:.1f}\u00b0 / {d['spread']:.2f}\u00b0 "
            f"is outside panel b's field -- widen it")
        c = tint(SECOND, 0.26)
        ax.add_patch(plt.Circle((X(d["band"]), Y(d["spread"])), DOT_R,
                                facecolor=c, edgecolor="none", zorder=5))

    rule(ax, fx, fx + fw, b, lw=0.55, color=SECONDARY)
    for t in XT:
        ax.text(X(t), b - 1.4, _degfmt(t), fontproperties=F_SMALL,
                color=SECONDARY, ha="center", va="top", zorder=6)
    y = b - 1.4 - line_h(F_SMALL) * 1.55
    axis_x(ax, fx + fw / 2, y, "Angular range covered (\u00b0)")
    P = dict(rows=D, band=m_band, rep=m_rep, ratio=m_ratio,
             under=sum(1 for d in D if d["ratio"] < 1.0))
    return P, y - line_h(F_HEAD) * 0.9


def panel_c(ax, x0, top, w):
    """One released entry opened up: its two dihedral-PCA landscapes."""
    letter(ax, "d", x0, top - 2.6)
    d = land_exemplar()
    h = ROW2 - 11.0
    b = top - 5.4 - h
    fw = (w - 4.0 - PGUT * 2) / 2
    surface(ax, x0 + PGUT, b, fw, h, d["peptide"], "Peptide")
    surface(ax, x0 + PGUT * 2 + fw + 4.0, b, fw, h, d["cdr3"], "CDR3\u03b1+\u03b2")
    ycb = b - 1.3 - line_h(F_SMALL) * 0.95 - 2.6
    # 14.0, not 12.0: the unit label is set 1.6 mm past the bar's right end and
    # runs 9.8 mm, and this panel is now the right-hand column, so the bar has
    # to end far enough inside the plate for the label to stay on it.
    cbar(ax, x0 + PGUT, ycb, w - 14.0 - PGUT, 2.0, FELMAP, 0, ZMAX, (0, 3, 6),
         "kcal/mol", clip=True)
    return d, ycb - 1.0 - line_h(F_SMALL) * 1.5


def panel_d(ax, x0, top, w):
    """Subspace overlap between replicas, against what chance would give."""
    letter(ax, "c", x0, top - 2.6)
    h = ROW2 - 11.0
    b = top - 5.4 - h
    P = strip(ax, x0, b, w, h, "rmsip", 0.9, (0, .3, .6, .9),
              lambda t: f"{t:.1f}", THIRD, THIRD_D, "Subspace RMSIP",
              sf2, rangelab="Replica range (RMSIP)",
              band=(RAND_LO, RAND_HI))
    ax.text(x0 + AXW + 0.8, b + (RAND_HI + 0.012) * h / 0.9, "chance",
            fontproperties=F_SMALL, color=SECONDARY, ha="left", va="bottom",
            zorder=7)
    return P, _foot(ax, P, b)


# ------------------------------------------------------------------ figure
def draw(ax, top, w):
    col = (w - COLGAP) / 2
    dcol = col + 3.0

    # row 1: a Q, b the incident angle's variation.  row 2: c RMSIP against
    # chance, d the entry's two landscapes.  The letters run a-b-c-d in the
    # order the text first cites them.
    A, f1 = panel_a(ax, 0, top, col)
    B, _ = panel_b(ax, col + COLGAP, top, col)
    y = f1 - VGAP
    ecol = w - dcol - COLGAP
    E, fd = panel_d(ax, 0, y, ecol)
    d, fc = panel_c(ax, ecol + COLGAP, y, dcol)
    floor = min(fc, fd)

    pid, run = d["traj_id"].rsplit("_run", 1)
    # numbered in order of depth, capped at FELN; the caption counts what is
    # actually drawn rather than promising the cap
    nb = tuple(min(len(d[k]["basins"]), FELN) for k in ("peptide", "cdr3"))
    wd = {1: "one", 2: "two", 3: "three", 4: "four"}
    nbtxt = (f"the {wd[nb[0]]} basins each surface resolves are numbered"
             if nb[0] == nb[1] else
             f"the {wd[nb[0]]} basins of the peptide surface and the "
             f"{wd[nb[1]]} of the CDR3 surface are numbered")

    fn = [r["spread"] for r in A["rows"]]
    rat = [d["ratio"] for d in B["rows"]]
    rm = [r["mean"] for r in E["rows"]]
    rs = [r["spread"] for r in E["rows"]]
    LEGEND.update(
        title="Technical validation of the deposited trajectories",
        subtitle=(
            f"a and c, all {N_SET} complexes, ranked left to right by their "
            f"own means. A segment spans the lowest to the highest of the three "
            f"replica means and the tick is that mean, so segment height is the "
            f"disagreement between independent runs; the narrow panel bins those heights "
            f"on the same axis. a, *Q*, the fraction of native contacts retained "
            f"(median range {sf2(quant(fn,.5))}, interquartile range {sf2(quant(fn,.25))}\u2013"
            f"{sf2(quant(fn,.75))}; ICC {icc([r['values'] for r in A['rows']]):.2f}). "
            f"b, incident-angle variation, one point per complex: the angular "
            f"interval the complex's pooled frames occupy, p5 to p95, against the "
            f"range of its three replica means, both in degrees on logarithmic axes. "
            f"The replica range is the smaller of the two in "
            # "245 of 245" reads as a coincidence; when the comparison holds
            # everywhere the sentence should say so.
            f"{'all ' + str(N_SET) if B['under'] == N_SET else str(B['under']) + f' of {N_SET}'} "
            f"complexes, by a median factor of {sf2(quant(rat,.5))} (interquartile "
            f"range {sf2(quant(rat,.25))}\u2013{sf2(quant(rat,.75))}); the solid line is "
            f"equality and the dashed line that median factor. The sectors above the "
            f"scatter are the two medians drawn as angles at true scale, "
            f"{sf2(B['band'])}\u00b0 covered against {sf2(B['rep'])}\u00b0 run to run. "
            f"c, overlap of the top ten "
            f"C\u03b1 principal modes between replicas (median {sf2(quant(rm,.5))}, range "
            f"{sf2(min(rm))}\u2013{sf2(max(rm))}; median between-replica range {sf2(quant(rs,.5))}); "
            f"the shaded band is what ten modes drawn at random would give, "
            f"\u221a(10/3N) = {sf2(RAND_LO)}\u2013{sf2(RAND_HI)} here. "
            f"d, dihedral-PCA free-energy surfaces of the peptide and the "
            f"CDR3\u03b1+\u03b2 loops of one released entry, {pid.upper()} run {run}, on the "
            f"first two dihedral principal components with the variance each carries; "
            f"\u0394G is measured above each surface's own minimum, iso-lines "
            f"every 1 kcal/mol, and {nbtxt} by depth."),
        notes=[
            f"A basin in d is a local minimum within {d['basin_zmax']:.1f} kcal/mol "
            f"of the global minimum holding at least {d['basin_min_pop']:.0%} of the "
            f"frames; the surfaces are estimated on an "
            f"{len(d['peptide']['fel']['x'])}\u00d7{len(d['peptide']['fel']['y'])} grid and "
            f"truncated at {ZMAX:.0f} kcal/mol for display. The release stores the same two "
            f"surfaces and basin set for all 735 trajectories. d shows "
            f"what a single download contains and is not a validation statistic.",
        ],
    )
    return notes(ax, 0, w, floor, LEGEND["notes"])


# The width each panel occupies inside the plate.  A panel file is rendered at
# that same width, so the four dropped side by side at their own size rebuild
# the plate; anyone wanting a different arrangement scales vector art.
_COL = (W - COLGAP) / 2
PANELS = {"a": (panel_a, _COL),
          "b": (panel_b, _COL),
          "c": (panel_d, W - (_COL + 3.0) - COLGAP),
          "d": (panel_c, _COL + 3.0)}

# A panel file gets the same 2.4 mm of white on all four sides that render()
# already leaves above and below.  Inside the plate every panel is flush with
# its own column -- the letter starts at x = 0 and the marginal ends at x = w --
# which is right there and wrong in a standalone file, where ink on the trim
# edge reads as clipped and is clipped by any placement frame that rounds down.
# The drawn width is unchanged, so a panel still sets 1:1 against the plate.
PAD = 2.4

if __name__ == "__main__":
    render(draw, "fig4")
    out = Path(__file__).resolve().parent / "panels"
    for k, (fn, wid) in PANELS.items():
        render(lambda ax, top, w, _f=fn: _f(ax, PAD, top, w - 2 * PAD)[1],
               f"fig4{k}", width=wid + 2 * PAD, pad=PAD, out=out)
