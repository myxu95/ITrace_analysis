"""Table 1 figure: analysis category -> where it lives in a trajectory record.

Scientific Data single-column text width (141 mm, measured from a published PDF).
Three columns: category | summary key inside analysis.json | detailed sidecar files.
"""
from figstyle import *          # canvas, type scale, ink, save

# The tables keep the pre-2026-09 type scale.  The figures were enlarged because
# they are reduced on placement; these tables are placed at their drawn width,
# and every column here is measured in mm against text set at 7 pt.
FS_BODY, FS_SMALL, FS_HEAD = 7.0, 6.5, 7.0
F_BODY, F_SMALL, F_HEAD = font(FS_BODY), font(FS_SMALL), font(FS_HEAD, bold=True)

F_CAT = F_BODY                        # category text
F_KEY = font(FS_SMALL, mono=True)     # analysis.json keys
F_FILE = font(FS_SMALL, mono=True)    # sidecar file names
F_NOTE = F_BODY

# fully achromatic: the three columns are separated by typeface and weight,
# not by hue.  Both text colours are dark -- INK for the identifiers, dark
# slate for the secondary column; nothing here is set in grey.
KEY, FILE, HEADC, NOTEC, DASH = INK, SECONDARY, INK, SECONDARY, SECONDARY

# (category, [summary keys in analysis.json], [detailed sidecar files])
#
# Row order follows the Methods 2.6 bullet list, not the order the keys are
# physically written into analysis.json.  The two disagree on exactly two
# entries: the pipeline writes tcr_cdr and geometry early, but the narrative
# puts the CDR-loop decomposition after the collective-motion analyses and the
# geometric metrics near the end, because that pass runs over the trajectory
# on its own.  A reader mapping the bullets onto this table should be able to
# do it line by line, so the table follows the text.
ROWS = [
    ("Biological identity",              ["identity"],                    ["meta.json"]),
    ("Flexibility",                      ["rmsf", "rmsf_profile"],        []),
    ("Interface size",                   ["bsa"],                         []),
    ("Interface partitioning",           ["bsa_decomposition"],           ["struct_metrics.json"]),
    ("Contacts",                         ["contact"],                     ["contacts.csv"]),
    ("Interaction inventory",            ["interactions"],                ["interactions_pairs.json"]),
    ("TCR docking geometry",             ["angle"],                       []),
    ("Regional RMSD",                    ["rmsd_regions"],                ["rmsd_regions.json"]),
    ("Fraction of native contacts (*Q*)", ["fnat"],                       ["rmsd_regions.json"]),
    ("Concerted motion",                 ["concerted_motion"],            ["concerted_motion.json"]),
    ("Essential dynamics",               ["essential_dynamics"],          ["essential_dynamics.json"]),
    ("CDR-loop decomposition",           ["tcr_cdr"],                     []),
    ("Peptide conformational landscape", ["peptide_dpca"],                ["peptide_dihedrals.json",
                                                                          "basin_*.pdb", "groove_ref.pdb"]),
    ("CDR3 conformational landscape",    ["tcr_cdr3_dpca"],               ["tcr_cdr3_dpca.json",
                                                                          "cdr3_basin_*.pdb",
                                                                          "cdr3_framework_ref.pdb"]),
    ("Peptide–CDR3 state coupling",      ["coupled_states"],              ["coupled_states.json"]),
    ("Geometric metrics",                ["geometry"],                    ["struct_metrics.json"]),
    ("Interface hotspots",               ["interface"],                   []),
]

LEGEND = {   # prose harvested by build_legends.py; see figstyle.LEAN
    "title": "Where each analysis category is recorded",
    "subtitle": "analysis.json carries the summary keys listed in the middle column. The "
                "detailed time series, matrices and per-residue tables are in the files at "
                "right; an em dash means the category has no sidecar beyond the summary.",
    "notes": [],
}

C1, C2, GAP = 48.0, 32.0, 4.0
C3 = W - C1 - C2 - 2 * GAP
X1, X2, X3 = 0.0, C1 + GAP, C1 + C2 + 2 * GAP
SEP = "  "


tw = text_w


def pack(items, avail, fp):
    """One identifier per line - unambiguous, and every item still fits its column."""
    return list(items)


packed = [(cat, pack(k, C2, F_KEY), pack(f, C3, F_FILE)) for cat, k, f in ROWS]

LH, PADR = 3.25, 2.0
heights = [max(1, len(k), len(f)) * LH + 2 * PADR for _, k, f in packed]
H_HEAD, H_NOTE = 6.0, 0.0   # the note is legend prose now; see LEGEND
H = 1.2 + H_HEAD + sum(heights) + H_NOTE + 1.2

fig, ax = canvas(H, flip=True)

y = 1.2
ax.plot([0, W], [y, y], lw=0.9, color=INK, solid_capstyle="butt")
for x, t in ((X1, "Analysis category"), (X2, "Key in analysis.json"), (X3, "Detailed record files")):
    ax.text(x, y + H_HEAD / 2 + 0.15, t, fontproperties=F_HEAD, color=HEADC, va="center")
y += H_HEAD
ax.plot([0, W], [y, y], lw=0.45, color=INK, solid_capstyle="butt")

for i, ((cat, klines, flines), h) in enumerate(zip(packed, heights)):
    if i:
        ax.plot([0, W], [y, y], lw=0.22, color=HAIR, solid_capstyle="butt", zorder=1)
    cx = X1
    for t, fp in parts(cat, F_CAT):   # *Q* is a symbol, so it is set italic
        ax.text(cx, y + PADR + LH / 2 + 0.1, t, fontproperties=fp, color=INK,
                va="center", zorder=2)
        cx += tw(t + "|", fp) - tw("|", fp)
    for col, lines, fp, c in ((X2, klines, F_KEY, KEY), (X3, flines, F_FILE, FILE)):
        if not lines:
            ax.text(col, y + PADR + LH / 2 + 0.1, "—", fontproperties=fp, color=DASH, va="center", zorder=2)
            continue
        y0 = y + PADR + LH / 2
        for j, ln in enumerate(lines):
            ax.text(col, y0 + j * LH + 0.1, ln, fontproperties=fp, color=c, va="center", zorder=2)
    y += h

ax.plot([0, W], [y, y], lw=0.9, color=INK, solid_capstyle="butt")

save(fig, "Table1")
print(f"   cols {C1}/{C2}/{C3:.1f}")
bad = [(c, l, round(tw(l, fp), 1), lim) for c, k, f in packed
       for lines, fp, lim in ((k, F_KEY, C2), (f, F_FILE, C3)) for l in lines if tw(l, fp) > lim]
print("overflow:", bad or "none")
print("widest category:", max((round(parts_w(c, F_CAT), 1), c) for c, _, _ in packed), f"/ {C1}")
