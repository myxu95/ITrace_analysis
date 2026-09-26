"""Set several finished panels side by side on one vector canvas.

Each panel keeps the drawing code it uses when it is rendered alone: it is
given its own sub-axes whose data coordinates are still millimetres, so the
composed figure is the same vector art as the standalone panels rather than a
sheet of pasted images, and a panel edited later re-composes without being
touched here.  What this module owns is only what a single panel cannot know --
the column split, the gutters, and the panel letters.
"""
from pathlib import Path

import matplotlib.pyplot as plt

from figstyle import F_TAG, H_MAX, INK, MM, W, line_h, save, tag


def measure(draw, width):
    """Natural height of a panel at ``width``, in mm, ink to ink."""
    scratch = plt.figure(figsize=(1, 1)).add_subplot()
    floor = draw(scratch, 400.0, width)
    plt.close(scratch.figure)
    return 400.0 - floor


def compose(name, rows, width=W, gx=5.0, gy=7.0, pad=2.4, lead=None, out=None):
    """``rows`` is a list of rows, each a list of ``(letter, draw)`` panels.

    Panels in a row split the width equally; a third element sets a relative
    weight when they should not.  Every panel in a row hangs from the same top
    edge, which is what makes a row read as a row.
    """
    lead = line_h(F_TAG) if lead is None else lead

    plan, total = [], 2 * pad + gy * (len(rows) - 1)
    for row in rows:
        wts = [c[2] if len(c) > 2 else 1.0 for c in row]
        free = width - gx * (len(row) - 1)
        ws = [free * v / sum(wts) for v in wts]
        hs = [measure(c[1], v) for c, v in zip(row, ws)]
        plan.append((ws, hs, row))
        total += lead + max(hs)

    fig = plt.figure(figsize=(width * MM, total * MM))
    fig.patch.set_facecolor("white")
    over = fig.add_axes([0, 0, 1, 1], zorder=99)
    over.set_xlim(0, width), over.set_ylim(0, total)
    over.set_aspect("equal"), over.axis("off"), over.patch.set_visible(False)

    top = total - pad
    for ws, hs, row in plan:
        x = 0.0
        for (letter, draw, *_), pw, ph in zip(row, ws, hs):
            ax = fig.add_axes([x / width, (top - lead - ph) / total,
                               pw / width, ph / total])
            ax.set_xlim(0, pw), ax.set_ylim(0, ph)
            ax.set_aspect("equal"), ax.axis("off"), ax.patch.set_visible(False)
            draw(ax, ph, pw)
            if letter:
                tag(over, letter, x, top - lead + 0.85)
            x += pw + gx
        top -= lead + max(hs) + gy

    save(fig, name, out)
    print(f"{name}: {width:.0f} x {total:.0f} mm"
          + (f"  OVER PAGE ({H_MAX:.0f} mm)" if total > H_MAX else ""))
    return fig
