#!/usr/bin/env python
"""The conformational landscapes every entry carries, and how the two sides couple.

A trajectory is only useful if something can be read off it, and what this
archive reads off each one is where the two mobile halves of the interface sit
and whether they sit there together.  Three blocks, all recomputed by the
pipeline from the deposited frames and all served with the entry.

Left and centre, a dihedral principal-component analysis of each half -- the
peptide, then the CDR3 alpha and beta loops taken together -- rendered as a free
energy surface, kcal/mol above that surface's own minimum.  Dihedral and not
Cartesian: torsions are internal coordinates, so the surface cannot pick up the
rigid-body drift that a Cartesian PCA of a solvated complex spends its first
components on.  The CDR3 loops are superposed on the Valpha/Vbeta framework, so
what is measured is loop motion relative to the domain carrying it.

Right, the joint occupancy of the two basin sets over the 1001 frames.  This is
the block that makes the pair of landscapes more than two separate pictures: two
independent halves would fill the matrix in proportion to the row and column
totals, and this one does not -- half its cells are empty.

BASIN, NOT SUBSTATE.  Each dPCA block stores two decompositions of the same
scores: basins, which are local minima of the surface drawn here, and a k-means
substate count over the raw scores, which is a different number.  The matrix is
built by assigning each frame to the nearest basin centre, so basins are what
index it, and basins are what this figure numbers, draws and names.

Numbers, not areas, in the matrix.  Sixteen cells is small enough to read
exactly, and the reader needs the exact zeros: an empty cell is the evidence
that the two halves are not moving independently, and a pale square is a weaker
claim than 0.00.

Four colours for four basins, and the same four in the matrix margin.  The
basins are a nominal index -- 1 is only the deepest, and nothing follows from
1 being next to 2 -- so they take four separable hues rather than four steps of
one, and the matrix carries its row and column indices as the same numbered
discs.  That is the whole reason the colours are here: a cell of the matrix is a
pair of basins, and a reader should be able to put a finger on a well in one
landscape and find the row it indexes without counting.  The four deliberately
avoid the blue of the surfaces they sit on and the blue/ochre/teal that mean
replica 1, 2 and 3 in the companion figure, so no marker here can be read as a
replica.

One trajectory is drawn, not one complex.  The matrix is a property of a single
run -- it counts that run's frames -- and averaging three of them would blur
four basins into a smear.  The legend carries the cross-replica agreement, which
is the quantity that says whether the drawn run is representative of its two
siblings, and the library distribution behind it.
"""
import sys
from pathlib import Path

# figstyle / landdata live one level up, shared by every figure folder
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from collections import Counter

from figstyle import *
from landdata import N_COMPLEX, N_TRAJ, VERDICTS, exemplar, quant, trajectories
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe

LEGEND = {}   # prose harvested by build_legends.py; see figstyle.LEAN

PH = 42.0                       # landscape height
AXW = 8.0                       # left gutter of each landscape: the rotated PC2 title
GAP = 10.0                      # between blocks
MATW = 23.0                     # joint-occupancy matrix width
ZMAX = 6.0                      # the pipeline clips the surface here; do not re-clip
FILL = [i * 0.2 for i in range(int(ZMAX / 0.2) + 1)]    # fine, so the wash is smooth
ISO = [i * 1.0 for i in range(1, int(ZMAX) + 1)]       # coarse, so the shape is readable
CBH = 2.6                       # colour-bar height
CBTICKS = [0, 2, 4, 6]
DOT_R = 1.45                    # numbered basin marker; holds a 6.5 pt digit
MDOT_R = 1.15                   # the same marker as a matrix row/column index
# basin 1-4.  Four hues and not four shades: the index is nominal.  All four are
# dark enough to hold white type over the whole range of the surface beneath
# them, and none is the blue/ochre/teal that mean replica 1/2/3 elsewhere in
# Figure 4.
BASIN = ["#5d4585", "#9c4a52", "#96762e", "#3f7a52"]
STOPS = (0.02, 0.30, 0.62, 0.86, 0.97)      # dark at the minimum, white at the rim
CMAP = LinearSegmentedColormap.from_list("fel", [tint(PRIMARY_D, f) for f in STOPS])
# coupled_states.py::_verdict -- fixed cut-offs on the mutual information, not
# quantiles of this library, so the words mean the same thing outside it
CUTS = [0.05, 0.20, 0.40]


def landscape(ax, x0, base, pw, d, title):
    """One free-energy surface with its basins numbered in stored (ΔG) order."""
    fel = d["fel"]
    X, Y, Z = fel["x"], fel["y"], fel["z"]
    sx = pw / (X[-1] - X[0])
    sy = PH / (Y[-1] - Y[0])
    xs = [x0 + (v - X[0]) * sx for v in X]
    ys = [base + (v - Y[0]) * sy for v in Y]
    zc = [[min(max(v, 0.0), ZMAX) for v in row] for row in Z]
    ax.contourf(xs, ys, zc, levels=FILL, cmap=CMAP, zorder=2)
    ax.contour(xs, ys, zc, levels=ISO, colors="white", linewidths=0.25, zorder=3)
    for i, b in enumerate(d["basins"], 1):
        bx, by = x0 + (b["pc1"] - X[0]) * sx, base + (b["pc2"] - Y[0]) * sy
        ax.add_patch(plt.Circle((bx, by), DOT_R, facecolor=BASIN[(i - 1) % len(BASIN)],
                                edgecolor="white", lw=0.7, zorder=7))
        ax.text(bx, by, str(i), fontproperties=font(FS_SMALL, bold=True), color="white",
                ha="center", va="center", zorder=8)
    ax.add_patch(plt.Rectangle((x0, base), pw, PH, facecolor="none", edgecolor=SECONDARY,
                               lw=0.5, zorder=9))
    ax.text(x0 + pw / 2, base + PH + 1.4, title, fontproperties=F_HEAD, color=INK,
            ha="center", va="baseline", zorder=6)
    ax.text(x0 + pw / 2, base - 1.8, f"PC1 ({d['pc1_var_frac']:.0%})",
            fontproperties=F_HEAD, color=INK, ha="center", va="top", zorder=6)
    axis_y(ax, base + PH / 2, f"PC2 ({d['pc2_var_frac']:.0%})", x=x0 - AXW + 3.2)


def draw(ax, top, w):
    d = exemplar()
    pid = d["traj_id"].rsplit("_run", 1)[0]
    run = d["traj_id"].rsplit("_run", 1)[1]
    cs = d["coupled"]
    M = cs["matrix"]
    nr, nc = len(M), len(M[0])
    rows = trajectories()

    verd = Counter(r["coupling"] for r in rows)
    defined = [r for r in rows if r["coupling"] != "n/a"]
    why = Counter(r["reason"] for r in rows if r["coupling"] == "n/a")
    nmi = sorted(r["nmi"] for r in defined)
    pep_b = Counter(r["pep_basins"] for r in rows)
    cdr_b = Counter(r["cdr3_basins"] for r in rows)
    v_pep = sorted(r["pep_pc1"] + r["pep_pc2"] for r in rows)
    v_cdr = sorted(r["cdr3_pc1"] + r["cdr3_pc2"] for r in rows)
    repro = sum(1 for r in rows if r["reproducible"])
    tp = cs["top_pair"]
    empty = sum(1 for row in M for v in row if v == 0)
    multi = sum(n for k, n in pep_b.items() if k > 1)
    multi_c = sum(n for k, n in cdr_b.items() if k > 1)

    LEGEND.update(
        title="Interface conformational landscapes",
        subtitle=f"{pid} run {run}, one of the {N_TRAJ} deposited trajectories. Left and "
                 f"centre, the dihedral-PCA free-energy surface of the peptide and of the "
                 f"CDR3α+β loops: backbone φ/ψ, side-chain χ1/χ2 and Cα–Cα distances, block-"
                 f"normalised, projected on their first two components and binned over the "
                 f"1001 frames; shading is ΔG above each surface's own minimum, iso-lines every "
                 f"1 kcal/mol, and the surface is clipped at {ZMAX:.0f}. Markers number the "
                 f"basins in order of increasing ΔG, one colour each. Right, the fraction of "
                 f"frames in each peptide-basin × CDR3-basin pair; the margin repeats the "
                 f"markers, so a cell is the pair of basins named on its row and column.",
        notes=[f"A basin is a local minimum of the surface drawn here, and it is what the "
               f"matrix is indexed by: every frame is assigned to its nearest basin centre. "
               f"The pipeline also stores a k-means substate count over the same component "
               f"scores; that is a different quantity and is not what is shown.",
               f"{empty} of the {nr * nc} cells are empty. Peptide basin {tp['peptide'] + 1} and "
               f"CDR3 basin {tp['cdr3'] + 1} are the most enriched pair, holding "
               f"{tp['occ']:.1%} of frames against the {tp['occ'] / tp['enrichment']:.1%} that "
               f"independent halves would give — an enrichment of {tp['enrichment']:.1f}×. "
               f"Normalised mutual information between the two assignments is "
               f"{cs['nmi']:.2f}; the pipeline calls ≥{CUTS[2]:.2f} strong, {CUTS[1]:.2f}–"
               f"{CUTS[2]:.2f} moderate, {CUTS[0]:.2f}–{CUTS[1]:.2f} weak and below "
               f"{CUTS[0]:.2f} independent, on fixed cut-offs rather than on quantiles of this "
               f"library. All three replicas of this entry are coupled "
               f"({cs['replica_agreement']['n_coupled']} of "
               f"{cs['replica_agreement']['n']}).",
               f"Both surfaces and the matrix are stored for all {N_TRAJ} trajectories of the "
               f"{N_COMPLEX} complexes. {multi} resolve more than one peptide basin and "
               f"{multi_c} more than one CDR3 basin, which leaves {verd['n/a']} trajectories "
               f"whose matrix is degenerate — {why['peptide single-state']} with a single "
               f"peptide basin, {why['CDR3 single-state']} with a single CDR3 basin — and no "
               f"mutual information to compute. Of the remaining {len(defined)}, "
               + ", ".join(f"{verd[v]} are {v}" for v in VERDICTS[:-1])
               + f", with a median of {quant(nmi, .5):.2f} (interquartile range "
                 f"{quant(nmi, .25):.2f}–{quant(nmi, .75):.2f}).",
               f"The first two components carry a median {quant(v_pep, .5):.0%} of the peptide's "
               f"dihedral variance and {quant(v_cdr, .5):.0%} of the loops', so each surface is "
               f"a projection and not the whole space; basins that are distinct on it are "
               f"distinct, while frames that coincide on it need not be. Coupling is drawn for "
               f"one run because the matrix counts that run's frames, and it reproduces: for "
               f"{repro} of the {N_TRAJ} trajectories the entry's three replicas agree on "
               f"whether the two halves are coupled."],
    )
    y = header(ax, 0, top, w, LEGEND["title"], LEGEND["subtitle"])

    base = y - 5.5 - PH
    pw = (w - AXW * 2 - GAP * 2 - MATW) / 2
    landscape(ax, AXW, base, pw, d["peptide"], "Peptide")
    x2 = AXW * 2 + pw + GAP
    landscape(ax, x2, base, pw, d["cdr3"], "CDR3 loops (α + β)")

    # the matrix: numbers, because the exact zeros are the evidence
    mx0 = x2 + pw + GAP
    cw = MATW / nc
    ch = min(PH / nr, cw * 1.15)
    mtop = base + PH
    vmax = max(v for row in M for v in row)
    for i, row in enumerate(M):
        for j, v in enumerate(row):
            f = 0.96 - min(v / vmax, 1.0) * 0.88
            rect(ax, mx0 + j * cw, mtop - (i + 1) * ch, cw, ch, tint(PRIMARY_D, f),
                 ec="white", lw=0.6, z=3)
            ax.text(mx0 + (j + .5) * cw, mtop - (i + .5) * ch, f"{v:.2f}",
                    fontproperties=F_SMALL, color="white" if f < 0.45 else INK,
                    ha="center", va="center", zorder=5)
    # the margin indices are the landscape markers, shrunk: the matrix is indexed
    # by the basins drawn to its left, and saying so with the same object costs
    # nothing and removes a step of counting
    def mdot(cx, cy, k):
        ax.add_patch(plt.Circle((cx, cy), MDOT_R, facecolor=BASIN[k % len(BASIN)],
                                edgecolor="white", lw=0.6, zorder=5))
        ax.text(cx, cy, str(k + 1), fontproperties=font(FS_SMALL, bold=True),
                color="white", ha="center", va="center", zorder=6)

    for j in range(nc):
        mdot(mx0 + (j + .5) * cw, mtop + 1.0 + MDOT_R, j)
    for i in range(nr):
        mdot(mx0 - 1.4 - MDOT_R, mtop - (i + .5) * ch, i)
    ymat = mtop - nr * ch
    ax.text(mx0 + MATW / 2, mtop + 1.6 + 2 * MDOT_R, "Joint occupancy",
            fontproperties=F_HEAD, color=INK, ha="center", va="baseline", zorder=6)
    ax.text(mx0 + MATW / 2, ymat - 1.8, "CDR3 basin", fontproperties=F_HEAD, color=INK,
            ha="center", va="top", zorder=6)
    axis_y(ax, mtop - nr * ch / 2, "Peptide basin", x=mx0 - 5.6)
    ax.text(mx0 + MATW / 2, ymat - 1.8 - line_h(F_HEAD) * 1.45,
            f"NMI {cs['nmi']:.2f} · {cs['coupling']}", fontproperties=F_HEAD,
            color=INK, ha="center", va="top", zorder=6)

    # ΔG scale, spanning the two surfaces it applies to
    ycb = base - 1.8 - line_h(F_HEAD) * 1.7 - CBH
    cbw = x2 + pw - AXW
    n = 96
    for i in range(n):
        rect(ax, AXW + i * cbw / n, ycb, cbw / n + 0.03, CBH,
             CMAP(i / (n - 1.0)), z=3)
    ax.add_patch(plt.Rectangle((AXW, ycb), cbw, CBH, facecolor="none", edgecolor=SECONDARY,
                               lw=0.4, zorder=6))
    for v in CBTICKS:
        cx = AXW + v / ZMAX * cbw
        ax.text(cx, ycb - 1.2, f"≥{v}" if v == ZMAX else str(v), fontproperties=F_SMALL,
                color=SECONDARY, ha="left" if v == 0 else "right" if v == ZMAX else "center",
                va="top", zorder=6)
    yc = ycb - 1.2 - line_h(F_SMALL) * 0.95 - line_h(F_HEAD) * 0.95
    axis_x(ax, AXW + cbw / 2, yc, "ΔG above the free-energy minimum (kcal/mol)")
    return notes(ax, 0, w, yc - line_h(F_BODY) * 0.60, LEGEND["notes"])


if __name__ == "__main__":
    render(draw, "fig_interface_pca")
