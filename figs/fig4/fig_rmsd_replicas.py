#!/usr/bin/env python
"""What one entry's three replicas look like, and where that entry sits in the library.

The companion figure reduces every trajectory to two scalars and plots all 735
at once.  That answers "did the release settle" but not "what am I downloading",
and a user downloads one entry: three trajectories of one complex, started from
the same minimised structure and differing only in their randomly drawn initial
velocities.  This figure draws those three curves in full.

Drawing one complex on its own would invite the obvious objection that it was
chosen for looking good, so the library is drawn behind it: at each of the 1001
time points, the 5th, 25th, 50th, 75th and 95th percentile of all 735 deposited
trajectories.  The reader can then see for themselves where this entry falls,
and the choice of exemplar stops being something the text has to be trusted on.

Percentiles rather than a mean with a standard-deviation ribbon.  At every time
point the distribution is right-skewed -- a minority of large, mobile complexes
sit far above a tight majority -- so a symmetric ribbon would draw an upper edge
below curves that exist and a lower edge below the smallest RMSD in the set.

What replica agreement is and is not.  Three trajectories of the same molecule
are not expected to superpose frame by frame; molecular dynamics is chaotic and
identical curves would mean the runs were not independent.  What is expected is
that they occupy the same range once settled, which is why the right-hand
marginal draws each replica's own 5th-to-95th percentile band over the settled
window on the same axis as the traces.  Three bands that overlap are the result;
one band displaced from the other two would be the finding.

The vertical dashes mark 20 ns, where that window opens.  The first 10% of a run
is a minimised structure relaxing into the thermostat, which belongs to the
setup rather than to the sampling, and it is dropped here for the same reason
and at the same point as in the per-trajectory table.

Three hues, not three shades of one.  The three curves cross each other dozens
of times over 200 ns, and at every crossing a reader has to be able to say which
curve came out the other side; graded shades of one blue cannot do that, because
the shade a line appears to have depends on what it is lying on top of.  The
hues carry no order -- the replicas are exchangeable, differing only in their
initial velocities -- they are only labels.  The marginal bars take the same
three, left to right, so the marginal doubles as the key for the traces.
"""
import sys
from pathlib import Path

# figstyle / rmsddata / driftdata live one level up, shared by every figure folder
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from compdata import load as composition
from driftdata import trajectories
from figstyle import *
from rmsddata import N_COMPLEX, N_TRAJ, envelope, exemplar, quant, smooth

LEGEND = {}   # prose harvested by build_legends.py; see figstyle.LEAN

PH = 46.0                       # main field height
AXW = 10.0                      # left tick gutter
GAP = 5.5                       # main field to marginal
MW = 15.0                       # right marginal width
KEYH = 4.6                      # key row above the field
YMAX = 8.0                      # Angstrom; the library's 95th percentile peaks at 7.7
YTICKS = [0, 2, 4, 6, 8]
XTICKS = [0, 50, 100, 150, 200]
TAIL = 20.0                     # ns; start of the settled window, as in build_drift.py
SMOOTH = 5                      # +-5 stored points = a 2 ns running mean
# replicas 1-3: three hues at the same weight, so no replica reads as dominant
REP = [(PRIMARY, PRIMARY_D), (SECOND, SECOND_D), (THIRD, THIRD_D)]
RAW_F = 0.62                    # tint of the stored series, under its own running mean
BAR_F = 0.55                    # tint of the marginal bars
SW_W, SW_H = 4.6, 1.8           # key swatch


def draw(ax, top, w):
    t, q = envelope()
    pid, _, runs = exemplar()
    A = 10.0                                            # nm -> Angstrom, everywhere
    qs = {p: [v * A for v in smooth(q[p], 8)] for p in q}
    raw = {k: [v * A for v in runs[k]] for k in runs}
    sm = {k: smooth(raw[k], SMOOTH) for k in raw}
    keep = [i for i, v in enumerate(t) if v >= TAIL]

    ent = next(r for r in composition() if r["pdb_id"] == pid)
    lib = trajectories()
    ex = {r["run"]: r for r in lib if r["pdb_id"] == pid}
    tails = sorted(r["tail90"] * A for r in lib)
    bands = sorted(r["band"] * A for r in lib)
    pct = lambda s, v: 100.0 * sum(1 for x in s if x < v) / len(s)
    mean = lambda v: sum(v) / len(v)
    tail_ex = {k: mean([raw[k][i] for i in keep]) for k in raw}
    band_ex = {k: (lambda s: quant(s, 0.95) - quant(s, 0.05))(
        sorted(raw[k][i] for i in keep)) for k in raw}
    i20 = keep[0]
    spread = max(tail_ex.values()) - min(tail_ex.values())
    # the three replicas' settled ranges, and how much of that is shared
    los = [quant(sorted(raw[k][i] for i in keep), 0.05) for k in sorted(raw)]
    his = [quant(sorted(raw[k][i] for i in keep), 0.95) for k in sorted(raw)]
    overlap = (min(his) - max(los)) / (max(his) - min(los))
    assert max(qs[95]) < YMAX and max(max(v) for v in raw.values()) < YMAX

    LEGEND.update(
        title="Replica agreement in backbone RMSD",
        subtitle=f"The three replicas of {pid} — {ent['peptide_seq']} on {ent['mhc_allele']} — "
                 f"drawn against the whole library. Backbone RMSD to the energy-minimised "
                 f"starting structure, 1001 stored frames over 200 ns. Grey, the pointwise 5th–"
                 f"95th and 25th–75th percentile of all {N_TRAJ} deposited trajectories at each "
                 f"time point, with their median; the three replicas in three colours, the thin "
                 f"line the stored series and the heavy line a "
                 f"{2 * SMOOTH * 0.2:.0f} ns running mean. "
                 f"Dashes mark {TAIL:.0f} ns, where the settled window opens. Right, on the same "
                 f"axis, each replica's own 5th–95th percentile band and median over "
                 f"{TAIL:.0f}–200 ns; each bar carries its trace's colour.",
        notes=[f"The replicas share a starting structure and differ only in their randomly "
               f"drawn initial velocities, so they are not expected to superpose frame by frame "
               f"— identical curves would mean the runs were not independent. What is expected "
               f"is a shared range, and the three settled bands overlap over "
               f"{overlap:.0%} of their combined extent. Their means over {TAIL:.0f}–200 ns are "
               + ", ".join(f"{tail_ex[k]:.2f}" for k in sorted(tail_ex))
               + f" Å, a spread of {spread:.2f} Å, and their bands are "
               + ", ".join(f"{band_ex[k]:.2f}" for k in sorted(band_ex)) + " Å.",
               f"The exemplar is a released entry, not a curated one. Across the library the "
               f"mean RMSD over the settled window has median {quant(tails, .5):.2f} Å "
               f"(interquartile range {quant(tails, .25):.2f}–{quant(tails, .75):.2f} Å) and the "
               f"band median {quant(bands, .5):.2f} Å ({quant(bands, .25):.2f}–"
               f"{quant(bands, .75):.2f} Å); these three trajectories fall at the "
               + ", ".join(f"{pct(tails, ex[k]['tail90'] * A):.0f}" for k in sorted(ex))
               + " percentile on the first and the "
               + ", ".join(f"{pct(bands, ex[k]['band'] * A):.0f}" for k in sorted(ex))
               + " percentile on the second.",
               f"The library median rises from {qs[50][0]:.1f} Å at the start to "
               f"{qs[50][i20]:.1f} Å by {TAIL:.0f} ns and then to {qs[50][-1]:.1f} Å over the "
               f"remaining {200 - TAIL:.0f} ns, so most of the movement away from the starting "
               f"coordinates happens in the first tenth of a run and the rest of it is "
               f"fluctuation. The envelope widens upward rather than shifting: its 5th "
               f"percentile ends at {qs[5][-1]:.1f} Å against {qs[95][-1]:.1f} Å for the 95th, "
               f"which is the spread of complex sizes and interface mobilities in the archive, "
               f"not a spread in how well the runs equilibrated.",
               f"Every trajectory in the release is in the envelope and none is excluded; the "
               f"{N_TRAJ} series are the same ones served as rmsd.json, so both the band and "
               f"the exemplar can be recomputed from the deposited data."],
    )
    y = header(ax, 0, top, w, LEGEND["title"], LEGEND["subtitle"])

    px0 = AXW
    pw = w - AXW - GAP - MW
    mx0 = px0 + pw + GAP
    sx = pw / 200.0
    X = [px0 + v * sx for v in t]

    # key above the field: the library's three grey elements, then the replicas
    # as one entry whose three segments run left to right in replica order, the
    # same order the marginal puts them in and the only place the mapping is set
    yk = y - KEYH
    kx = px0
    def swatch(fn, lab, wide=SW_W):
        nonlocal kx
        fn(kx)
        ax.text(kx + wide + 1.5, yk + SW_H / 2, lab, fontproperties=F_SMALL, color=INK,
                ha="left", va="center", zorder=6)
        kx += wide + 1.5 + text_w(lab, F_SMALL) + 5.2

    for f, lab in ((0.90, f"{N_TRAJ} trajectories, 5–95%"), (0.78, "25–75%")):
        swatch(lambda x, f=f: rect(ax, x, yk, SW_W, SW_H, tint(SECONDARY, f), z=5), lab)
    swatch(lambda x: ax.plot([x, x + SW_W], [yk + SW_H / 2] * 2, lw=0.6,
                             color=tint(SECONDARY, 0.30), zorder=6,
                             solid_capstyle="butt"), "median")
    seg = (SW_W * 1.5 - 1.0) / 3.0
    swatch(lambda x: [ax.plot([x + k * (seg + 0.5), x + k * (seg + 0.5) + seg],
                              [yk + SW_H / 2] * 2, lw=1.1, color=REP[k][1], zorder=6,
                              solid_capstyle="butt") for k in range(3)],
           f"{pid} replicas 1, 2, 3", wide=SW_W * 1.5)
    assert kx < w, kx

    base = yk - 3.4 - PH
    sy = PH / YMAX
    Y = lambda v: base + v * sy

    for v in YTICKS:
        if v:
            ax.plot([px0, mx0 + MW], [Y(v)] * 2, lw=0.25, color=HAIR, zorder=1,
                    solid_capstyle="butt")
        ax.text(px0 - 1.6, Y(v), str(v), fontproperties=F_SMALL, color=SECONDARY,
                ha="right", va="center", zorder=7)
    axis_y(ax, base + PH / 2, "Backbone RMSD (Å)")

    ax.fill_between(X, [Y(v) for v in qs[5]], [Y(v) for v in qs[95]],
                    color=tint(SECONDARY, 0.90), lw=0, zorder=2)
    ax.fill_between(X, [Y(v) for v in qs[25]], [Y(v) for v in qs[75]],
                    color=tint(SECONDARY, 0.78), lw=0, zorder=3)
    ax.plot(X, [Y(v) for v in qs[50]], lw=0.6, color=tint(SECONDARY, 0.30), zorder=4,
            solid_capstyle="round")

    for k in sorted(raw):                       # stored series, under its own mean
        ax.plot(X, [Y(v) for v in raw[k]], lw=0.18,
                color=tint(REP[k - 1][0], RAW_F), zorder=5, solid_capstyle="butt")
    for k in sorted(sm):
        ax.plot(X, [Y(v) for v in sm[k]], lw=0.8, color=REP[k - 1][1],
                zorder=6, solid_capstyle="round")

    ax.plot([px0 + TAIL * sx] * 2, [base, base + PH], lw=0.5, color=INK, zorder=7,
            dashes=(1.8, 1.4), solid_capstyle="butt")

    rule(ax, px0, px0 + pw, base, lw=0.7)
    for v in XTICKS:
        cx = px0 + v * sx
        ax.plot([cx, cx], [base, base - 0.9], lw=0.5, color=INK, zorder=7,
                solid_capstyle="butt")
        ax.text(cx, base - 1.7, str(v), fontproperties=F_SMALL, color=SECONDARY,
                ha="right" if v == XTICKS[-1] else "center", va="top", zorder=7)

    # the marginal shares the RMSD axis with the field, which is why it sits here:
    # each replica's settled range read against the same ruler as its trace
    bw = MW / 3.8
    for k in sorted(raw):
        s = sorted(raw[k][i] for i in keep)
        lo, hi, med = quant(s, 0.05), quant(s, 0.95), quant(s, 0.5)
        cx = mx0 + (k - 0.5) * (MW / 3)
        rect(ax, cx - bw / 2, Y(lo), bw, (hi - lo) * sy,
             tint(REP[k - 1][0], BAR_F), z=4)
        ax.plot([cx - bw / 2, cx + bw / 2], [Y(med)] * 2, lw=0.9,
                color=REP[k - 1][1], zorder=6, solid_capstyle="butt")
        ax.text(cx, base - 1.7, str(k), fontproperties=F_SMALL, color=SECONDARY,
                ha="center", va="top", zorder=7)
    rule(ax, mx0, mx0 + MW, base, lw=0.7)
    ax.text(mx0 + MW / 2, base + PH + 1.4, f"{TAIL:.0f}–200 ns", fontproperties=F_HEAD,
            color=INK, ha="center", va="baseline", zorder=6)

    yl = base - 1.7 - line_h(F_SMALL) * 0.95
    yt = yl - line_h(F_BODY) * 0.95
    axis_x(ax, px0 + pw / 2, yt, "Simulation time (ns)")
    axis_x(ax, mx0 + MW / 2, yt, "Replica")
    return notes(ax, 0, w, yt - line_h(F_BODY) * 0.60, LEGEND["notes"])


if __name__ == "__main__":
    render(draw, "fig_rmsd_replicas")
