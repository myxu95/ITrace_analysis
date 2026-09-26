"""Shared drawing style for every ITrace Scientific Data figure.

One import so that all figures agree on canvas width, typeface, type scale,
ink colours, palette and export format.  A figure script starts with

    from figstyle import *

and ends with ``save(fig, "fig_name")`` rather than calling savefig itself.

Journal constraints encoded here
--------------------------------
Scientific Data is single-column: the text block measures 141 mm across
(measured from a published PDF), so 141 mm is the FULL width, not a half
width -- Nature's generic 183 mm does not apply.  Text height is ~247 mm.
Figure lettering must be <= 7 pt.  Files must be RGB with no alpha channel
on a white opaque background.

Two house rules
---------------
LEGIBILITY RULE.  Nothing a reader has to read is set below 6.5 pt and
nothing is set in grey.  Grey is for rules and hairlines only; every glyph
that carries information is INK or, at the very lightest, SECONDARY -- which
is a dark slate, not a grey.  Figures are drawn one per file at their own
natural size; if a panel needs 5 pt to fit, the panel is too small and gets
its own figure instead of being shrunk.

HUE RULE.  Colour encodes data, never typography.  Tables, archive trees and
flow charts stay achromatic.  A hue appears only when it carries a variable,
and the default is one colour plus one contrast colour for the exception
class -- not a different hue per bar.

HOUSE PALETTE.  Figure 1 is the schematic that opens the paper, so it, not
the data figures, sets the colour of the article.  Every hue below is sampled
from it: the periwinkle of its headings and badge, the ochre of its rendered
complex, the teal of its contour panel, and the three light steps it uses for
box strokes and fills.  Its geometry is borrowed too -- ``bar()`` caps a bar
with the same rounded end as a Fig 1 pill.  ``palette()`` is the one exception
and says why in its own docstring.
"""
from pathlib import Path

import colorsys
import inspect
import math

import matplotlib
matplotlib.use("Agg")
# SVG is the editable deliverable, so its text stays text rather than being
# converted to outlines; PDF and PNG remain the fixed-appearance pair.
matplotlib.rcParams["svg.fonttype"] = "none"
# Type 42, not matplotlib's default Type 3.  A Type 3 font ships its glyphs as
# uninterpreted drawing operators, so the text is vector but is no longer text:
# it cannot be searched, re-encoded or re-set by a typesetter, and production
# desks reject it for exactly that reason.  Type 42 embeds the TrueType outlines
# with their character codes intact.  Nothing about the rendered page changes.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MPath
from matplotlib.textpath import TextPath

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------- canvas
MM = 1 / 25.4
W = 141.0                 # single-column text width, the full width here
H_MAX = 247.0             # text height; a figure taller than this spills
DPI = 600

# ---------------------------------------------------------------- type
SANS = "Nimbus Sans"      # Helvetica clone; Sci Data figures are Helvetica/Arial
MONO = "Noto Sans Mono"   # file names, JSON keys

# THE FIGURE SCALE.  Raised 2026-09-03.  These figures are assembled several to
# a sheet and then placed into the manuscript at less than their drawn width, so
# a panel set at the 7 pt journal ceiling printed at roughly 4 pt -- legible in
# the source file and not on the page.  The scale is therefore set for what
# survives that reduction rather than for what the standalone panel looks like,
# and every figure in the set is now drawn at one width (W) so the reduction is
# the same for all of them and no two panels can print at different sizes.
# FS_FLOOR is unchanged: it guards the smallest mark, not the nominal scale.
# The four table scripts pin themselves back to the old 7 / 6.5 pt -- their
# columns are measured in mm against text set at 7 pt and do not hold at 10.
FS_TAG = 12.0             # panel letter (bold lower-case a, b, c)
FS_TITLE = 11.0           # figure / panel title
FS_HEAD = 10.0            # column and axis titles
FS_BODY = 10.0            # the default for everything
FS_SMALL = 8.5            # tick labels, file names
FS_TINY = 8.5             # kept as an alias so nothing renders below the floor
FS_FLOOR = 6.5

def font(size=FS_BODY, bold=False, mono=False):
    if size < FS_FLOOR:
        raise ValueError(f"{size} pt is below the {FS_FLOOR} pt legibility floor")
    return FontProperties(family=MONO if mono else SANS, size=size,
                          weight="bold" if bold else "normal")

F_TAG = font(FS_TAG, bold=True)
F_TITLE = font(FS_TITLE, bold=True)
F_HEAD = font(FS_HEAD, bold=True)
F_BODY = font(FS_BODY)
F_SMALL = font(FS_SMALL)
F_TINY = F_SMALL

def text_w(s, fp):
    """Rendered width of a string in mm, for fitting text into columns."""
    if not s:
        return 0.0
    return TextPath((0, 0), s, prop=fp,
                    size=fp.get_size_in_points()).get_extents().width / 72.0 / MM

def line_h(fp, lead=1.32):
    return fp.get_size_in_points() * lead / 72.0 / MM


# A symbol is a variable and is set in italic wherever it appears -- Q on an
# axis is the same Q as in the running text, and the two have to look alike.
# The markup keeps a label one readable literal at the call site instead of a
# list of runs: "Native-contact fraction, *Q*".
def parts(s, fp):
    """Split "roman *italic* roman" into [(text, fontproperties), ...]."""
    it = FontProperties(family=fp.get_family(), size=fp.get_size(),
                        weight=fp.get_weight(), style="italic")
    return [(t, it if i % 2 else fp)
            for i, t in enumerate(s.split("*")) if t]


def parts_w(s, fp):
    """Width in mm of a marked-up string, measured run by run.

    Measured with a sentinel so a run ending in a space keeps its advance,
    which ``text_w`` alone drops.
    """
    return sum(text_w(t + "|", f) - text_w("|", f) for t, f in parts(s, fp))

# ---------------------------------------------------------------- ink
# text colours -- both dark enough to read at 6.5 pt in print
INK = "#12212c"           # primary text, heavy rules
SECONDARY = "#3b4956"     # subtitles, units, secondary columns (dark slate)
SUBTLE = SECONDARY        # legacy name

# line colours -- NEVER used for text.  Tinted toward Fig 1's periwinkle so
# that a hairline in a data figure and a divider in the schematic agree.
MUTE = "#9ba6b5"          # tick marks, empty-cell dashes
HAIR = "#e1e7ef"          # row separators, light rules
# NEUTRAL / NEUTRAL_EDGE are the chord's J sectors and their leader lines, and
# stay a true neutral: against a full-wheel V palette a periwinkle band stops
# reading as "no hue here" and starts reading as one more category.
NEUTRAL = ("#b9c3cb", "#d5dce1")   # alternating neutral bands (fills)
NEUTRAL_EDGE = "#8f99a3"
FILEC = SECONDARY         # legacy name

# ---------------------------------------------------------------- hue
# Sampled from figs/fig1.png -- see HOUSE PALETTE above.
PRIMARY = "#3a68ae"       # default data colour: Fig 1's heading / badge blue
PRIMARY_D = "#254779"     # its darker tone, for emphasis
# Desaturated one notch from the Figure 1 sample (S 0.48 -> 0.37).  At the
# sampled saturation a panel drawn entirely in ochre outweighs the same panel
# drawn in PRIMARY, because warm hues advance; taking the chroma out levels
# the two without changing which hue the reader is looking at.
SECOND = "#b17d53"        # the contrast colour: exceptions, a second class
SECOND_D = "#7c5332"
# A third data hue, for the case SECOND must not be spent on: two measurements
# of the same population shown side by side, where neither is an exception and
# reusing SECOND would make one of them read as one.  It is ACCENT's hue -- Fig
# 1's contour teal -- lifted to the lightness and saturation PRIMARY and SECOND
# share, so the three sit at the same weight on a page; ACCENT itself stays the
# darker single-highlight tone it was sampled as.
# Greyed for the same reason and one more: teal sits 35 degrees from PRIMARY on
# the wheel, and at full chroma two panels side by side read as two blues in
# competition rather than as two quantities.  S 0.34 -> 0.26 settles it back.
THIRD = "#4c7f80"
THIRD_D = "#315858"
ACCENT = "#417071"        # rare single-highlight case: Fig 1's contour teal

# The light steps are Fig 1's own rounded boxes, used unchanged, so a bar
# track here and a panel there are literally the same colour.
STROKE = "#95aad3"        # Fig 1 box stroke -- edges, leader lines
BAND = "#cad4e9"          # its palest stroke -- alternating bands
GROUND = "#eef4fb"        # Fig 1 box fill -- bar tracks, plot grounds

def tint(hexc, f):
    """Blend a colour f of the way toward white (0 = unchanged, 1 = white)."""
    r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
    m = lambda v: int(round(v + (255 - v) * f))
    return f"#{m(r):02x}{m(g):02x}{m(b):02x}"

def palette(n, start=12.0):
    """n maximally separated hues, for many-category panels (chord sectors).

    This is the ONE place the house palette is deliberately not applied.  A
    chord ring is read by tracing a ribbon across it, so its sectors need to
    be told apart by hue alone; confining them to Figure 1's blue-to-ochre arc
    made a thirty-gene ring beautiful and unusable.  Full wheel, full
    saturation.  Lightness and saturation alternate so that neighbouring
    categories stay distinguishable even when their sectors touch.
    """
    out = []
    for i in range(n):
        h = (i * 360.0 / n + start) % 360 / 360.0
        light = 0.44 if i % 2 == 0 else 0.58
        sat = 0.70 if i % 2 == 0 else 0.62
        cap = 0.30 if i % 2 == 0 else 0.50
        r, g, b = colorsys.hls_to_rgb(h, light, sat)
        # equalise PERCEIVED lightness: HLS lightness is not luminance, so a
        # yellow at L=0.58 prints far paler than a blue at the same L
        y = 0.2126 * r + 0.7152 * g + 0.0722 * b
        if y > cap:
            r, g, b = (v * cap / y for v in (r, g, b))
        out.append(f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}")
    return out

# a sixth-order categorical ramp, all of it inside Fig 1's family
CAT = [PRIMARY, SECOND, "#407268", "#755d98", "#a58a40", "#6b7d94"]

# ---------------------------------------------------------------- drawing
def canvas(h, w=W, flip=False):
    """A figure whose data coordinates are millimetres.

    Origin is bottom-left; pass ``flip=True`` for top-down layout (tables,
    trees, flow charts), where y grows downwards from the top edge.
    """
    fig = plt.figure(figsize=(w * MM, h * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0) if flip else ax.set_ylim(0, h)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig, ax

def rect(ax, x, y, w, h, fc, ec="none", lw=0.0, z=3):
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                               lw=lw, zorder=z))

_K = 0.5523                       # cubic control offset for a quarter circle
_ENDS = {"top": (2, 3), "right": (1, 2), "left": (0, 3),
         "bottom": (0, 1), "all": (0, 1, 2, 3), "none": ()}

def bar(ax, x, y, w, h, fc, ends="all", r=0.55, ec="none", lw=0.0, z=3):
    """A bar whose free end is capped like a Figure 1 pill.

    ``ends`` names the side that is free -- "top" for a column standing on a
    baseline, "right" for a row growing out of an axis, "all" for a floating
    swatch.  The seated end stays square, so the rounding reads as a cap and
    not as a pill that has come loose from its axis.  The radius is clamped to
    the bar, so a one-count bar shrinks to a lozenge instead of overshooting.
    """
    r = min(r, abs(w) / 2, abs(h) / 2)
    if r <= 0.03:
        return rect(ax, x, y, w, h, fc, ec, lw, z)
    corners = _ENDS[ends]
    pts = [(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    verts, codes = [], []
    for i, b in enumerate(pts):
        if i not in corners:
            verts.append(b)
            codes.append(MPath.LINETO)
            continue
        out = []
        for other in (pts[i - 1], pts[(i + 1) % 4]):       # in-edge, out-edge
            dx, dy = other[0] - b[0], other[1] - b[1]
            d = math.hypot(dx, dy) or 1.0
            out.append((dx / d, dy / d))
        (ux, uy), (vx, vy) = out
        verts += [(b[0] + ux * r, b[1] + uy * r),
                  (b[0] + ux * r * (1 - _K), b[1] + uy * r * (1 - _K)),
                  (b[0] + vx * r * (1 - _K), b[1] + vy * r * (1 - _K)),
                  (b[0] + vx * r, b[1] + vy * r)]
        codes += [MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    codes[0] = MPath.MOVETO
    verts.append(verts[0])
    codes.append(MPath.CLOSEPOLY)
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=fc, edgecolor=ec,
                           lw=lw, zorder=z))

def box(ax, cx, w, lo, q1, med, q3, hi, fc, ec, lw=0.5, z=3, med_c=INK, med_lw=0.9):
    """A vertical box-and-whisker, in millimetre coordinates.

    ``lo q1 med q3 hi`` are already y positions in mm -- the caller owns the
    scale, exactly as it owns a bar's height.  Used where a distribution is
    bimodal enough that a mean would misdescribe it: a box says how the library
    is spread, a bar only says where its centre landed.  The whisker caps are
    half the box width and the median is a full-width rule, so the box stays
    readable when the fill is saturated and when the box is only a millimetre
    tall.
    """
    ax.plot([cx, cx], [lo, hi], lw=lw, color=ec, zorder=z, solid_capstyle="butt")
    for yy in (lo, hi):
        ax.plot([cx - w / 4, cx + w / 4], [yy, yy], lw=lw, color=ec, zorder=z,
                solid_capstyle="butt")
    rect(ax, cx - w / 2, q1, w, max(q3 - q1, 0.05), fc, ec=ec, lw=lw, z=z + 1)
    ax.plot([cx - w / 2, cx + w / 2], [med, med], lw=med_lw, color=med_c,
            zorder=z + 2, solid_capstyle="butt")


def rule(ax, x0, x1, y, lw=0.6, color=INK, z=4):
    ax.plot([x0, x1], [y, y], lw=lw, color=color, zorder=z, solid_capstyle="butt")

# The ranked central-tendency line -- Fig 3c,d's median angle and Fig 4a,b,d's
# complex mean.  One constant, because the two figures draw the same kind of
# mark and a reader crossing between them should not have to decide whether a
# difference in weight means something.
#
# Raised from 1.10.  Both figures set this line over a light spread mark -- a
# percentile band in Fig 3, 245 replica segments in Fig 4 -- and at 1.10 the
# line was carrying less weight than the thing it is drawn on top of, which
# inverts the reading the house rule asks for: light for the spread, dark and
# heavy for the central tendency.  The line is the quantity a reader takes a
# number off; the spread is context around it.
LW_TREND = 1.50


def axis_y(ax, y, s, x=2.2):
    """Rotated axis title in the left gutter, centred on the plot height.

    Every axis title in the set goes through this function and :func:`axis_x`,
    so the treatment cannot drift apart again the way it had: they are the one
    class of text a reader must read before the data means anything, so they
    are the one class set bold, and they are always ink, never secondary.
    """
    p = parts(s, F_HEAD)
    if len(p) == 1:
        ax.text(x, y, s, fontproperties=F_HEAD, color=INK,
                ha="center", va="center", rotation=90, zorder=5)
        return
    # Rotated, the baseline runs bottom to top, so the runs are laid out along
    # y and centred as a block on the plot height; rotation_mode="anchor" makes
    # va the across-the-stroke direction, which is what centres them on x.
    yy = y - parts_w(s, F_HEAD) / 2
    for t, f in p:
        ax.text(x, yy, t, fontproperties=f, color=INK, ha="left", va="center",
                rotation=90, rotation_mode="anchor", zorder=5)
        yy += text_w(t + "|", f) - text_w("|", f)


def axis_x(ax, x, y, s, ha="center"):
    """Axis title under the plot, centred on ``x``.  See :func:`axis_y`.

    ``x`` is the midpoint of the field the title names -- ``px0 + pw / 2`` for a
    main plot, the marginal's own midpoint for a marginal -- never the midpoint
    of the canvas, which sits half a tick gutter to the left of the plot it
    would be claiming to label.  Centred rather than set at the axis end so the
    pair with :func:`axis_y` reads the same way on both edges.
    """
    ax.text(x, y, s, fontproperties=F_HEAD, color=INK,
            ha=ha, va="baseline", zorder=5)


def tag(ax, letter, x, y):
    """Panel label: bold capital, baseline-aligned across a row."""
    ax.text(x, y, letter, fontproperties=F_TAG, color=INK,
            ha="left", va="bottom", zorder=99)

def italic(size=FS_BODY, bold=False):
    return FontProperties(family=SANS, size=size, style="italic",
                          weight="bold" if bold else "normal")

def runs(ax, x, y, parts, z=5):
    """Set several styled runs on one baseline; returns the end x.

    ``parts`` is a list of (text, fontproperties, colour) -- used so that a
    species name can be italic inside an otherwise roman sentence.
    """
    for txt, fp, col in parts:
        ax.text(x, y, txt, fontproperties=fp, color=col, ha="left",
                va="baseline", zorder=z)
        # measure with a sentinel so a trailing space keeps its advance
        x += text_w(txt + "|", fp) - text_w("|", fp)
    return x

# A figure that goes into a journal carries only what is needed to READ it --
# axes, categories, counts, keys.  Every descriptive sentence belongs to the
# figure legend, which is manuscript text, not pixels.  Scripts still build
# their prose (it is full of computed numbers); LEAN keeps it out of the image
# and ``build_legends.py`` harvests the same strings into FIGURE_LEGENDS.md, so
# the legend file can never drift from the figure it describes.
LEAN = True


def header(ax, x, top, width, title, subtitle=None, letter=None):
    """Panel letter only.  The title and the subtitle are legend prose.

    Returns the y of the first free line below the block, so the callers'
    layout arithmetic is unchanged.

    A journal figure is an image; what it shows is stated in the figure legend,
    which is manuscript text and is typeset by the publisher.  A title drawn
    into the top-left corner therefore prints twice -- once as pixels the
    copy-editor cannot touch, once as the real legend -- and in a multi-panel
    assembly it competes with the panel letter for the same corner.  So under
    ``LEAN`` nothing is drawn but the letter, and ``build_legends.py`` harvests
    the same ``title`` / ``subtitle`` strings into FIGURE_LEGENDS.md.

    Set ``LEAN = False`` to render an annotated copy for internal review: then
    the title and subtitle are drawn, because a reviewer looking at fourteen
    loose PNGs needs to know which one is which.
    """
    if not LEAN:
        yb = top - line_h(F_TITLE) * 0.80
        x_t = x
        if letter is not None:
            tag(ax, letter, x, yb)
            x_t = x + text_w(letter, F_TAG) + 1.8
        ax.text(x_t, yb, title, fontproperties=F_TITLE, color=INK,
                ha="left", va="baseline", zorder=9)
        y = yb - line_h(F_TITLE) * 0.42
        if subtitle:
            for i, ln in enumerate(wrap(subtitle.split("\n"), width - 0.4)):
                y -= line_h(F_BODY) * (0.92 if i == 0 else 1.0)
                ax.text(x, y, ln, fontproperties=F_BODY, color=SECONDARY,
                        ha="left", va="baseline", zorder=9)
            y -= line_h(F_BODY) * 0.34
        return y

    if letter is None:
        return top
    yb = top - line_h(F_TAG) * 0.80
    tag(ax, letter, x, yb)
    return yb - line_h(F_TAG) * 0.42


def wrap(paras, width, fp=None):
    """Greedy word wrap in millimetres -- scripts author whole sentences."""
    fp = fp or F_BODY
    out = []
    for para in paras:
        cur = ""
        for word in para.split():
            trial = f"{cur} {word}".strip()
            if cur and text_w(trial, fp) > width:
                out.append(cur)
                cur = word
            else:
                cur = trial
        if cur:
            out.append(cur)
    return out


def notes(ax, x, width, y, lines, gap=1.6, rule_lw=0.6):
    """Prose under the plot, above a rule.  Silent under ``LEAN``.

    Returns the new floor, so a caller can hand it straight back to render().
    """
    if LEAN or not lines:
        return y
    y -= gap
    rule(ax, x, x + width, y, lw=rule_lw)
    y -= line_h(F_BODY) * 1.55
    for ln in wrap(lines, width - 0.4):
        ax.text(x, y, ln, fontproperties=F_BODY, color=SECONDARY,
                ha="left", va="baseline", zorder=5)
        y -= line_h(F_BODY)
    return y + line_h(F_BODY) * 0.35


def render(draw, name, width=W, pad=2.4, out=None):
    """Draw once on a scratch axes to learn the height, then draw for real.

    ``draw(ax, top, width)`` lays the figure out downwards from ``top`` and
    returns the y of its lowest ink.  The canvas height therefore follows the
    content instead of being guessed, which is what lets every figure be
    rendered on its own at its own natural size.
    """
    scratch = plt.figure(figsize=(1, 1)).add_subplot()
    floor = draw(scratch, 400.0, width)
    plt.close(scratch.figure)
    h = 400.0 - floor + 2 * pad
    fig, ax = canvas(h, width)
    draw(ax, h - pad, width)
    save(fig, name, out)
    return fig


# ---------------------------------------------------------------- export
def _caller_dir():
    """Directory of the outermost script that is not this module.

    A figure writes next to its own script, so figures grouped into fig2/,
    fig2/ ... land in their own folder without every script having to plumb
    an output path through render().
    """
    for fr in inspect.stack():
        d = Path(fr.filename).resolve()
        if d != Path(__file__).resolve():
            return d.parent
    return HERE

SPLIT_MARK = ".split-by-type"


def _dest(out, ext):
    """Where a ``.<ext>`` file goes: ``out``, or ``out/<ext>/`` if it opts in.

    A folder opts in by holding a ``.split-by-type`` marker file.  Opting in is
    a property of the directory and not of the script, so a folder that has been
    split stays split for every figure rendered into it -- including scripts
    written later, which is the only version of this that does not rot the first
    time someone adds a panel.  Flat output stays the default: most figure
    folders hold few enough files that three subdirectories would be worse than
    none, and nothing outside the opted-in folder changes behaviour.
    """
    d = out / ext if (out / SPLIT_MARK).exists() else out
    d.mkdir(parents=True, exist_ok=True)
    return d


def save(fig, name, out=None):
    """Write <name>.pdf, <name>.svg and a 600 dpi RGB <name>.png without alpha.

    Output goes beside the calling script unless ``out`` says otherwise, and is
    filed into per-format subfolders if that directory carries the marker that
    ``_dest`` looks for.
    """
    out = Path(out) if out else _caller_dir()
    # dpi=DPI on the vector saves too.  Geometry and type are unaffected -- both
    # backends normalise the page to 72 dpi before drawing -- but an embedded
    # raster is resampled to whatever dpi the save is given, and without it that
    # is rcParams["figure.dpi"], i.e. 100.  Figure 5's two molecular renders were
    # going into the PDF and the SVG at 100 dpi, 555 px across a 141 mm plate,
    # from 1825 px sources; anything made from those files -- including a PNG
    # rasterised back out of the SVG -- inherited the loss.
    fig.savefig(_dest(out, "pdf") / f"{name}.pdf", dpi=DPI, facecolor="white")
    svg = _dest(out, "svg") / f"{name}.svg"
    fig.savefig(svg, dpi=DPI, facecolor="white")
    # The SVG names the font it was set in, which is a Helvetica clone that a
    # machine opening the file may not have.  Naming the substitutes explicitly
    # keeps the file rendering as drawn instead of falling back to a serif.
    svg.write_text(svg.read_text().replace(
        f"font-family: '{SANS}'", f"font-family: '{SANS}', Helvetica, Arial, sans-serif"))
    png = _dest(out, "png") / f"{name}.png"
    fig.savefig(png, dpi=DPI, facecolor="white")
    try:
        from PIL import Image
        im = Image.open(png)
        if im.mode != "RGB":                       # Sci Data: no alpha channel
            Image.alpha_composite(
                Image.new("RGBA", im.size, "white"), im.convert("RGBA")
            ).convert("RGB").save(png, dpi=(DPI, DPI))
    except ImportError:
        pass
    w_mm, h_mm = (v / MM for v in fig.get_size_inches())
    flag = "  OVER TEXT HEIGHT" if h_mm > H_MAX else ""
    where = "  ->  pdf/ svg/ png/" if (out / SPLIT_MARK).exists() else ""
    print(f"{name}.pdf/.png   {w_mm:.0f} x {h_mm:.0f} mm{flag}{where}")
    over = _overflow(fig, w_mm, h_mm)
    if over:
        print("  TEXT OUTSIDE THE CANVAS:")
        for t, x0, x1, y0, y1 in over:
            print(f"    {t[:58]!r}  x {x0:6.1f}..{x1:6.1f}  y {y0:6.1f}..{y1:6.1f}")

def _overflow(fig, w_mm, h_mm, tol=0.15):
    """Every text artist whose rendered box leaves the page, in mm.

    Measured against the page rather than against each axes' own coordinates,
    which is the same thing for a single-panel figure and the only correct
    thing for a composed one, where a sub-axes' millimetres start at its own
    corner and would report a panel-local position as if it were a page one.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bad = []
    for ax in fig.axes:
        for t in ax.texts:
            if not t.get_text().strip():
                continue
            b = t.get_window_extent(r)
            x0, x1 = (v / fig.dpi * 25.4 for v in (b.x0, b.x1))
            y0, y1 = (v / fig.dpi * 25.4 for v in (b.y0, b.y1))
            if x0 < -tol or x1 > w_mm + tol or y0 < -tol or y1 > h_mm + tol:
                bad.append((t.get_text(), x0, x1, y0, y1))
    return bad
